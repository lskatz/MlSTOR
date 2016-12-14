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
  
  my $richseqHash=annotateAsm($asm,$$settings{mlstdir},$settings);

  # Print the richseq to stdout
  my $outseq=Bio::SeqIO->new(-format=>"genbank");
  for my $seq(values(%$richseqHash)){
    $outseq->write_seq($seq);
  }

  return 0;
}

sub annotateAsm{
  my($asm,$mlstDir,$settings)=@_;

  # Make the blast database
  my $blastdb="$$settings{tempdir}/assembly.fna";
  cp($asm,$blastdb);
  system("makeblastdb -dbtype nucl -in $blastdb 1>&2"); die if $?;

  # Prepare the richseq
  my %seq;
  my $seqin=Bio::SeqIO->new(-file=>$asm);
  while(my $seq=$seqin->next_seq){
    $seq{$seq->id}=$seq;

    # Add a source tag
    my $feature=Bio::SeqFeature::Generic->new(
      -start         =>1,
      -end           =>$seq->length,
      -primary       =>"source",
    );
    $seq{$seq->id}->add_SeqFeature($feature);
  }

  my @locus=glob("$mlstDir/*.fna");
  my $numloci=@locus;
  for(my $i=0;$i<$numloci;$i++){
    my $locusname=basename($locus[$i], qw(.fna));
    my $blastResult="$$settings{tempdir}/$locusname.tsv";

    # Run the blast and sort by highest score
    # TODO do one blast per thread instead of
    # multithreading blast.
    system("blastn -query $locus[$i] -db $blastdb -num_threads $$settings{numcpus} -outfmt 6 | sort -k12,12nr > $blastResult");
    die if $?;

    open(my $blsFh, "<", $blastResult) or die "ERROR: could not read $blastResult: $!";
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
        -tag           =>{ 
                           locus=>$locusname,
                           gene =>$allele,
                           note =>"Evidence:$0",
                         },
      );

      # If we're doing the quick method and not the smart method,
      # then just add the first feature with a good score
      if($$settings{annotationType} eq "quick"){
        $seq{$contig}->add_SeqFeature($feature);
        last;
      }
    }
  }

  return \%seq;

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

