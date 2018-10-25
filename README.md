## **Mobile Genetic Elements Retrieving Tool** - **MGERT**

*MGERT* is a computational pipeline for easy retrieving of MGE's coding sequences of a particular family from genome assemblies.
*MGERT* utilizes several established bioinformatic tools combined into single pipeline which hides different technical quirks from the inexperienced user.

### Requirements


- [RepeatModeler 1.0.11 ](http://www.repeatmasker.org/RepeatModeler/)
- [RepeatMasker  open-4.0.7](http://www.repeatmasker.org/RMDownload.html)
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [RPS-blast v2.7.1+](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#RPSBFtp)
- [ORFinder v0.4.1](ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/)
- awk

- [Python](https://www.python.org/) 3.5 or higher
- Python libraries:
    - pandas v0.21.0
    - matplotlib v2.1.0
    - Biopython v1.70

### Short description

The pipeline includes five steps:
1. *de novo* search for all MGEs in the genome assembly with *RepeatModeler*. 
This step results in a set of consensus sequences for every MGE class/family found (in *fasta* format).
Note, that the classification of the consensuses is made by the REPET package, and you can retrieve only those MGEs that were classified.
2. collecting particular consensuses and search for their matches in the genome assembly (using *RepeatMasker*).
3. excising of found matches from the genome assembly according to coordinates in annotation table from previous step.
4. search only for those sequences that contain Conserved Domain (CD), ORF and CD in this ORF (via successive runs of RPS-blast and ORFinder)
5. adding flanking regions to each sequence with ORF that contains CD.

You may run the pipeline from any of these steps, for instance if you have your own MGEs library to search in a genome.
During the steps 3 & 4 the pipeline creates several diagnostic plots and calculates descriptive statistics on the found sequences.

### Usage examples


#### Preparation steps

   - First, run configuration script with the following command:

```
./MGERT.py --configure
```
This command will create a configuration file *config.json* with all the necessary paths (see "requirements" section) and filenames MGERT uses. MGERT will try to find all the paths automatically. Unless it couldn't find them, it will prompt a user to enter a path or a filename.


   - Then make your local Conserved Domain Database (CDD): put [PSSM](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#CD_PSSM) files (with *\*smp* extension) to your working directory and run  MGERT with the following flag:

```
./MGERT.py --make_cdd
```
This command will create a directory *LocalCDD* with all the necessary files inside it and the path to this CDD will be added to the *config.json*.

Now you can run different parts of the pipeline.

#### Full pipeline run starting from *de novo* MGE search

The shortest way is to run MGERT with  *all-default* parameters:

```
./MGERT.py -T L1
```
This command runs search and retrieving of [L1](https://en.wikipedia.org/wiki/LINE1) retrotransposons' ORFs and flanking regions in a genome assembly (note: this file is specified in the config file!). All search parameters set to default (see "Parameters" section).

#### Pipeline run from an arbitrary step

There are four possible stages to run the pipeline from:

- `cons` -

Consider the situation when you already have a repeat library called `L1_consensi.fasta` and therefore there is no need to run *de novo* part of the pipeline. In this case simply type the following command:

```
./MGERT.py -T L1  --from_stage seq --lib L1_consensi.fasta
```

### List of parameters

`-h, --help` -  show help message and exit

`--configure` - run the configuration script

`--make_cdd` - make local Conservwd Domain database (CDD)

`-f,  --from_stage [stage]` - choose one of the following steps from which the pipeline have to be run (default `rmod` - start from the *de novo* stage):

 - `cons` - collect consensi sequences;
 - `seq` - collect repeat sequences from the assembly;
 - `orf` - collect repeat sequences containing specified conserved domain and ORFs the same conserved domain;
 - `flank` - add flankig regions of specified length to ORFs.

`-t, --threads [integer]` - number of threads;

`-C, --censor [url]` - pass URL of Censor classification results to MGERT;

`-o, --ori` - if specified RepeatMasker will additionally output the \*ori file

`-m, --merge [integer]` - merge all repeats within X bp into a single entry

`-e, --e_value [real]` - set expectation value (E) for RPS-Blast, default 0.01;

`-c, --start_codon [integer]` - start codon to use. 0 = 'ATG' only; 1 = 'ATG' and alternative initiation codons; 2 = any sense codon. Default = 0;

`-l, --min_length [integer]` - set minimum length of ORF, default 1000 bp;

`-s, --strand [plus/minus/both]` - output ORFs on specified strand only. Default 'plus';

`-g, --genetic_code [integer]` - genetic code to use (1-31, Default 1). See [ncbi](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for details;

`-T, --mge_type [Penelope/BovB/RTE/CR1/L1/LINE etc]` - specify the type of MGE to search. This is a mandatory option;

`-le, --left_end [integer]` - length of ORFs' left flanking region. Default 500 bp;

`-re, --right_end [integer]` - length of ORFs' right flanking region. Default 500 bp;

`-L, --lib [fasta file]` - library for RepeatMasker (in fasta format). Default none;

`-r, --rm_tab [RepeatMasker table]` - specify repeat masker table to use. Default none;

`-v, --version` - show program's version number and exit.
