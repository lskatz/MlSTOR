# MlSTOR

MlSTOR: Multilocus Sequence Typing to find Orthologous Regions

## Synopsis

Annotate a genome using only an MLST scheme.   

## Usage

    annotateMLST.pl assembly.fasta --mlstdir MLST > annotation.gbk
      --numcpus         1      number of cpus to use
      --mlstdir         ''     Location of MLST files. All
                               files must have .fna extension.
      --tempdir         ''     Location of temporary files (optional)

## Setup

In order to annotate with an MLST scheme, you must have an MLST scheme.  The easiest way to get an MLST scheme is by using `mlst` at https://github.com/tseemann/mlst.

Example(s) of other schemes are below.  Contributions are welcome in the form of pull requests for these schemes.  Either in documentation format in this readme or via a tar.gz file similar to the section below for Escherichia.

### EnteroBase Eschericia

#### Prepackaged

    cd db
    tar zxvf Escherichia.ref.tar.gz

#### Straight from the source

    mkdir Eschericia.cgMLSTv1
    wget http://enterobase.warwick.ac.uk/schemes/Escherichia.cgMLSTv1/cgMLST_ref.fasta -O - | \
      perl -MBio::Perl -e '$in=Bio::SeqIO->new(-fh=>\*STDIN); while($seq=$in->next_seq){$locus = $seq->id; $out=Bio::SeqIO->new(-file=>">>Eschericia.cgMLSTv1/$locus.fna"); $out->write_seq($seq);}'

