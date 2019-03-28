#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw/basename/;
use File::Temp qw/tempdir/;
use File::Copy qw/cp mv/;

# Use bioperl goodness to make the gbk files
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

use threads;

$0=basename $0;

sub logmsg{print STDERR "$0: @_\n"}

exit main();

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s mlstdir=s numcpus=i annotationType=s)) or die $!;
  $$settings{tempdir}||=tempdir("$0.XXXXXX",TMPDIR=>1,CLEANUP=>1);
  $$settings{mlstdir} || die "ERROR: need --mlstdir\n".usage();
  $$settings{numcpus}||=0;
  $$settings{annotationType}||="quick";
    $$settings{annotationType}=lc($$settings{annotationType});

  my($asm) =@ARGV;
  die usage() if(!$asm);

  my @mlstLocus = glob("$$settings{mlstdir}/*.fna");
  my $numLoci = scalar(@mlstLocus);
  if(!$numLoci){
    die "ERROR: no loci were found in $$settings{mlstdir}";
  } elsif($numLoci < 7){
    logmsg "WARNING: there are only $numLoci loci in this scheme.";
  }
  
  # Make the blast database
  my $blastdb="$$settings{tempdir}/assembly.fna";
  cp($asm,$blastdb) or die "ERROR: could not copy $asm to $blastdb: $!";
  system("makeblastdb -dbtype nucl -in $blastdb 1>&2"); die if $?;

  my $numLociPerThread = int($numLoci / $$settings{numcpus}) + 1;
  my @thr;
  for(0..$$settings{numcpus}-1){
    my @threadLocus = splice(@mlstLocus, 0, $numLociPerThread);
    $thr[$_] = threads->new(\&annotationWorker, $blastdb, \@threadLocus, $settings);
    logmsg scalar(@threadLocus)." loci being compared in thread ".$thr[$_]->tid;
  }

  # While we are annotating, go ahead and read in the assembly
  my %seq;
  my $seqin=Bio::SeqIO->new(-file=>$asm);
  while(my $seq=$seqin->next_seq){
    # Add a source tag
    my $feature=Bio::SeqFeature::Generic->new(
      -start       =>1,
      -end         =>$seq->length,
      -primary     =>"source",
    );
    $seq->add_SeqFeature($feature);

    $seq{$seq->id}=$seq;
  }

  # Join threads and add the sequence features as we go.
  my $featureCounter=0;
  for(my $t=0;$t<@thr;$t++){
    my $tmp = $thr[$t]->join;

    my $percentDone = sprintf("%0.2f",($t+1)/$$settings{numcpus} * 100);
    logmsg "Adding features for ".scalar(@$tmp)." loci ($percentDone% done)";
    for my $feature(@$tmp){
      my $seqid = $feature->seq_id;
      $seq{$seqid}->add_SeqFeature($feature);
      $featureCounter++;
    }
  }
  logmsg "Done adding all $featureCounter features.";

  # Print the genbank file
  logmsg "Printing to stdout";
  my $outseq=Bio::SeqIO->new(-format=>"genbank");
  for my $seq( sort{ $b->length <=> $a->length } values(%seq)){
    $outseq->write_seq($seq);
  }

  return 0;
}

sub annotationWorker{
  my($blastdb, $loci, $settings)=@_;
  my $numLoci = scalar(@$loci);
  
  my @feature;

  for(my $i=0;$i<$numLoci;$i++){
    my $locusname=basename($$loci[$i], qw(.fna));
    my $blastResult="$$settings{tempdir}/$locusname.tsv";

    # Run the blast and sort by highest score
    system("blastn -query $$loci[$i] -db $blastdb -num_threads 1 -outfmt 6 > $blastResult");
    die if $?;

    # Create a seq feature
    open(my $blsFh, "sort -k12,12nr $blastResult | ") or die "ERROR: could not read and sort $blastResult: $!";
    while(<$blsFh>){
      chomp;
      my($allele,$contig,$identity,$alnLength,$mismatch, $gap, $qstart, $qend, $sstart, $send, $evalue, $score)=split /\t/;

      # check strand
      my $strand=1;
      my($feat_start,$feat_end)=($sstart,$send);
      if($sstart > $send){
        $strand=-1;
        ($feat_start,$feat_end)=($send,$sstart)
      }
      
      my $feature=Bio::SeqFeature::Generic->new(
        -start         =>$feat_start,
        -end           =>$feat_end,
        -strand        =>$strand,
        -primary       =>"gene",
        -source_tag    => $0,
        -display_name  =>$locusname,
        -score         =>$score,
        -seq_id        =>$contig,
        -tag           =>{ 
                           locus=>$locusname,
                           gene =>$allele,
                           note =>"Evidence:$0",
                         },
      );

      # If we're doing the quick method and not the smart method,
      # then just add the first feature with a good score
      if($$settings{annotationType} eq "quick"){
        push(@feature, $feature);
        last;
      } else {
        die "ERROR: I do not understand annotationType other than quick right now";
      }
    }
    close $blsFh;
  }

  return \@feature;
}

sub usage{
  "$0: annotate a genome quickly using blast and an MLST scheme.
  Usage: $0 assembly.fasta --mlstdir MLST > annotation.gbk
  --numcpus         1      number of cpus to use
  --mlstdir         ''     Location of MLST files. All 
                           files must have .fna extension.
  --tempdir         ''     Location of temporary files (optional)
  --annotationType  quick  The type of annotation.
                           Quick indicates that the first
                           hit will be one annotation per
                           MLST locus.
  "
}

