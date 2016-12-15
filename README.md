#MlSTOR

MlSTOR: Multilocus Sequence Typing to find Orthologous Regions

##Synopsis

Annotate a genome using only an MLST scheme.  Finally a bioinformatics tool whose name you can be embarrassed by.  

##Usage

    annotateMLST.pl assembly.fasta --mlstdir MLST > annotation.gbk
      --numcpus         1      number of cpus to use
      --mlstdir         ''     Location of MLST files. All
                               files must have .fna extension.
      --tempdir         ''     Location of temporary files (optional)
      --annotationType  quick  The type of annotation.
                               Quick indicates that the first
                               hit will be one annotation per
                               MLST locus.

##Setup

In order to annotate with an MLST scheme, you must have an MLST scheme.  The easiest way to get an MLST scheme is by using `mlst` at https://github.com/tseemann/mlst.

On the other hand if you are in a hurry, you can get the core-genome MLST scheme from Pasteur Institut with this command:

    mkdir cgMLST
    cd cgMLST
    wget http://bigsdb.pasteur.fr/tmp/BIGSdb_7057_1481740268_42735.txt | cut -f 2 | tail -n +2 | xargs -P 12 -n 1 sh -c 'wget "http://bigsdb.pasteur.fr/perl/bigsdb/bigsdb.pl?db=pubmlst_listeria_seqdef_public&page=downloadAlleles&locus=$0" -O $0.fna'
