## **Mobile Genetic Elements Retrieving Tool** - **MGERT**

*MGERT* is a computational pipeline for easy retrieving of MGE's coding sequences of a particular family from genome assemblies.
*MGERT* utilizes several established bioinformatic tools combined into single pipeline which hides different technical quirks from an inexperienced user.

Table of contents:
  * [Requirements](#requirements)
  * [Short description](#short-description)
  * [Flowchart of MGERT pipeline](#flowchart-of-mgert-pipeline)
  * [Installation](#installation)
  * [Usage examples](#usage-examples)
    + [Preparation steps](#preparation-steps)
    + [Full pipeline run starting from *de novo* MGE search using RepeatModeler](#full-pipeline-run-starting-from--de-novo--mge-search-using-repeatmodeler)
    + [Pipeline run from an arbitrary step](#pipeline-run-from-an-arbitrary-step)
    + [Define the stage after which the pipeline should stop.](#define-the-stage-after-which-the-pipeline-should-stop)
  * [Output description](#output-description)
  * [List of arguments](#list-of-arguments)
    + [Required arguments](#required-arguments)
    + [Configuration arguments](#configuration-arguments)
    + [Optional arguments](#optional-arguments)


<a name="requirements"></a>
### Requirements


- [RepeatModeler 1.0.11 ](http://www.repeatmasker.org/RepeatModeler/)
- [RepeatMasker  open-4.0.7](http://www.repeatmasker.org/RMDownload.html)
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [RPS-blast v2.7.1+](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#RPSBFtp)
- [ORFfinder v0.4.1](https://www.ncbi.nlm.nih.gov/orffinder/) - direct FTP-link: 
 
  `ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz`

- [awk](https://en.wikipedia.org/wiki/AWK) (is a standard feature of most Unix-like operating systems)

- [Python](https://www.python.org/) 3.5 or higher
- Python libraries (should be easy to install via [pip](https://pypi.org/project/pip/) or [anaconda](https://www.anaconda.com/)):
    - [pandas v0.21.0](https://pandas.pydata.org/)
    - [matplotlib v2.1.0](https://matplotlib.org/)
    - [Biopython v1.70](https://biopython.org/)

<a name="short-description"></a>
### Short description

The pipeline includes five steps:
1. *de novo* search for all MGEs in the genome assembly with *RepeatModeler*. 
This step results in a set of consensus sequences for every MGE class/family found (in *fasta* format).
Note, that the classification of the consensuses is made by the RECON package (as the part of RepeatModeler pipeline), and you can retrieve only those MGEs that were classified.
2. collecting particular consensuses and search for their matches in the genome assembly (using *RepeatMasker*).
3. excising of found matches from the genome assembly according to coordinates in the annotation table from previous step.
4. search only for those sequences that contain Conserved Domains (CD), Open Reading Frame (ORF) and CD within this ORF (via successive runs of RPS-blast and ORFinder)
5. adding flanking regions to each sequence with the CD-encoding ORF.

You may run the pipeline from any of these steps, for instance, if you have your own MGEs library to search in the genome.
During the steps 3 & 4 the pipeline creates several diagnostic plots and calculates descriptive statistics on the found sequences (number, mean length, sd, median length, 25th and 75th length percentiles).


<a name="flowchart-of-mgert-pipeline"></a>
### Flowchart of MGERT pipeline

![flowchart](figures/flowchart.png)


<a name="installation"></a>
### Installation

Clone the repository and run installation script (requires administrator permissions):

```bash
git clone https://github.com/andrewgull/MGERT && cd MGERT
sudo ./install.sh

```

If you can't run installation script with sudo, you could place `MGERT.py` wherever you want, but you will be prompt 
to enter a path to `test_dataset.tgz` when running `MGERT.py --test` (see below for details).


<a name="usage-examples"></a>
### Usage examples


#### Preparation steps

   - First, run configuration script with the following command:

```bash
MGERT.py --configure
```
This command will create a configuration file *config.json* with all the necessary paths (see "Requirements" section) 
and filenames MGERT uses. MGERT will try to find all the paths automatically. Unless it couldn't find them, it will 
prompt a user to enter a path or a filename.

After the configuration step you may run MGERT with the option `--test` to check out whether everything works as 
it supposed to on a toy data set (despite the size of the dataset, it can take a while).

   - To validate ORFs of found TEs fast, you should create a local version of Conserved Domain Database (CDD). 
   To do this, download full Conserved Domain collection from the NCBI website: follow the [link to the Conserved Domain Database](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml), click on the **Conserved Domains** menu and choose **FTP** in the drop-down list. You will be redirected to the FTP site where you will find **cdd.tar.gz** archive. 
   You can download it using either browser or command line utility like `wget` (in the latter case use the direct link `ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz`).    
   To figure out what filename you need to extract corresponding [PSSM](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#CD_PSSM) file from the archive, go to [NCBI CDD](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd.shtml), type in the name of domain of interest (e.g. "RT")
   click **Search** button and see a list of related domains (e.g. "RT_like", "RT_pepA17", "RT_nLTR_like" etc).
   Clicking on any entry, you will see a short description, a hierarchy of related domains and their PSSM codes ("cd00304" for "RT_like") - and this code is exactly what you need.
   So, to extract PSSM file for RT-domain, run the following command:

```bash
tar -zxvf cdd.tar.gz cd00304.smp
```
      
   - Create a simple [CSV](https://en.wikipedia.org/wiki/Comma-separated_values) file (either comma or TAB delimited) that specifies PSSM-file - domain correspondence.
   Below you can see the correspondence file used for [Penelope retroelements](https://www.ncbi.nlm.nih.gov/pubmed/16093704) analysis:
   
   ```
   cd00304.smp  RT
   cd01648.smp  RT
   pfam00078.smp    RT
   pfam07727.smp    RT
   pfam13966.smp    RT
   cd01644.smp  RT
   cd01709.smp  RT
   cd10442.smp  EN
   cd00719.smp  EN 
   ```
    
 
   These files are used by MGERT to report only those ORFs that encodes for **both** domains (RT and EN) 
   regardless of what actual PSSM file produced a hit.
   
   - Finally, make local CDD: put [PSSM](https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#CD_PSSM) files (with *\*smp* extension) 
   to your working directory along with a CSV file specifying file - domain correspondence, and run MGERT with the following flag:

```bash
MGERT.py --make-cdd
```
This command will create a directory *LocalCDD* with all the necessary files inside it and the path to this CDD will be added to the *config.json*.


Now you can run the pipeline.


#### Full pipeline run starting from *de novo* MGE search using RepeatModeler

The shortest way is to run MGERT with  *all-default* parameters (see "Parameters" section):

```bash
MGERT.py --mge-type Penelope --assembly genome.fna.gz
```
This command runs search and retrieving of [Penelope](https://www.pnas.org/content/94/1/196) retrotransposons' ORFs and flanking regions in the genome assembly.

Note, that for subsequent runs of MGERT, you don't need to move and gzip the genome assembly file again. The only thing you 
should care about is that the name of the directory where the assembly file is, were the same as the name of assembly itself
(except extension) e.g. `./genome/genome.fna` 

#### Pipeline run from an arbitrary step

There are three possible steps to run the pipeline from (except the default one):

- consensus step

Let's consider the situation when you already have a repeat library called, say, `Penelope_consensi.fasta` and you want to find instances of the repeats from the library in your assembly, and therefore there is no need to run *de novo* part of the pipeline. In this case simply type in the following command:

```bash
MGERT.py --assembly genome.fna --mge-type Penelope  --from-stage cons --lib Penelope_consensi.fasta
```

If consensus library is not specified, it will be automatically generated from the RepeatModeler output.

Furthermore, after this step a table with descriptive statistics and a histogram of repeats' lengths will be generated (shown below)

```
count   81848.0
mean    1380.6
std     1987.7
min     12.0
25%     152.0
50%     446.0
75%     1945.0
max     27686.0
```

where:
 - *count* - number of found repeats/hits;
 - *mean* - mean length of found repeat/hit;
 - *std* - standard deviation;
 - *min* - minimal length of found repeat/hit;
 - *25%* - the 25th percentile (the 1st quartile) of length of found repeats/hits;
 - *50%* - the 50th percentile (the 2nd quartile or median) of length of found repeats/hits;
 - *75%* - the 75th percentile (the 3rd quartile) of length of found repeats/hits;
 - *max* - maximum length of found repeats/hits.

![histogram](figures/hist.png)

- coordinates step

In case when you have coordinates of repeats' matches (e.g. from previous step) - either *.out* or *.bed* file - you can run MGERT as follows:

```bash
MGERT.py --assembly genome.fna --mge-type MGE --from-stage coords --rm-table genome.fna.out
```

After this step a table with descriptive statistics and a histogram of repeats' lengths will be generated.

- ORFs step

In case when you have TEs sequences (in FASTA format, normally it's the output of the previous step) and want to find ORFs with conserved domains, run the following command:

```bash
MGERT.py --assembly genome.fna --mge-type MGE --from-stage orfs --sequence MGE_sequences.fasta
```

After this step a table with descriptive statistics and a histogram of repeats' lengths will be generated as well.

 - flanks step
 
Is useful if you want to add flanking regions of certain length to ORFs.

```bash
MGERT.py --assembly genome.fna --mge-type MGE --from-stage flanks
```

Note, that at this step input is taken automatically from the config file what is OK if the previous step was done,
otherwise MGERT will prompt you to enter the path to the fasta with ORFs.
Also, input ORFs must be the ones that ORFfinder produces, cause their headers contain all the information MGERT 
requires to excise them from a genome.

#### Define the stage after which the pipeline should stop.

Say, you want to run RepeatModeler only in order to check what types of TEs it will find. In this case run the command:

```bash
MGERT.py --assembly genome.fna.gz --to-stage rmod

MGERT.py --check-types ./genome/consensi.fa.classified

```

To stop MGERT after RepeatMasker run, use:

```bash
MGERT.py --assembly genome.fna.gz -T Penelope --to-stage coordinates

``` 

### Output description

MGERT outputs the following files:

 1. **rmod** step:
     - *RepModOut* directory with all the files produced by the RepeatModeler
     - *consensi.fa.classified* - a fasta file with all the TEs' consensus sequences found by the RepeatModeler
     - *All_TE_consensi.fa* - a fasta file with consensus sequences of specified TE (**input for the following step**)
     - *Unknown_consensi.fa.classified* - fasta file with unclassified consensus sequences
     - *Unknown_classified.fa* - a fasta file with unknown consensus sequences classified by the CENSOR software
 2. **consensus** step:
     - *TE_genome.fna.RMout.txt* - stdout of RepeatMasker
     - *TE_genome.fna.masked* - genome.fna with masked TE hits
     - *TE_genome.fna.out* - a table with TE hits coordinates (**input for the following step**)
     - *TE_genome.fna.ori.out* - a table with TE hits coordinates without RepeatMasker merging (**alternative input for the following step**)
     - *TE_genome.fna.out.gff* - an annotation of the TE hits in GFF format
     - *TE_genome.fna.tbl* - a table with number of bases masked by RepeatMasker (TE abundance percentage)
     - *genome.fna.fai* - genome index
     - *TE_genome.fna.out.bed* -  a table with TE hits coordinates in BED format (**alternative input for the following step**)
 3. **coordinates** step:
     - *TE_excised_matchesX.fa* - a fasta with TE hits retrieved from the genome. *X* stands for *--merge* value (**input for the following step**)
     - *TE_excised_matchesX.png* - a histogram of length distribution of corresponding sequences
     - *TE_excised_matchesX.stats* - a table with descriptive statistics
 4. **orfs** step:
     - *TE_matches_with_hits_eX.fa* - fasta file with those matches (from previous step) that contain CD hits. *X* stands for e-value cut off
     - *TE_matches_with_hits_eX.png* - a histogram of length distribution of corresponding sequences
     - *TE_matches_with_hits_eX.stats* - a table with descriptive statistics
     - *TE_cdsX_with_CD_eZ.fa* - a fasta file with TE's ORFs of min length *X* and e-value cut off *Z* (**input for the following step**)
     - *TE_cdsX_with_CD_eZ.faa* - a fasta file with TE's translated ORFs of min length *X* and e-value cut off *Z*
     - *TE_cdsX_with_CD_eZ.png* - a histogram of length distribution of corresponding nucleotide sequences
     - *TE_cdsX_with_CD_eZ.stats* - a table with descriptive statistics
 5. **flanks** step (optional):
     - *TE_cdsX_with_CD_eZ_extended_LaRb.fa* - a fasta file with TE's ORFs of min length *X* and e-value cut off *Z*, left flanking sequence of length *a* 
     and right flanking sequence of length *b*
 



<a name="list-of-arguments"></a>
### List of arguments

#### Required arguments

`-a, --assembly` - specify a genome assembly file (e.g. genome.fa.gz); this argument is mandatory on all stages since it
indicates where the working directory is.

`-T, --mge-type` - specify the type of MGE to search (e.g. L1/BovB/RTE/CR1/LINE/Penelope/DIRS etc.)


#### Configuration arguments

`-configure` - run the configuration script

`--make-cdd` - make local CDD

#### Optional arguments

`--test` - run self-test after configuration on a toy data set

`-cd, --cd-table` - specify a path to a comma or tab delimited table of SMP files and their grouping (e.g. domains.csv). CSV extension is mandatory.

`-f, --from-stage` - specify the step from which the pipeline should start: 'consensus' - get consensus sequences; 
'coords' - get sequences; 'orfs' - get ORFs; 'flanks' - add flanking sequences to CDS (default - *rmod*).

`-S, --to-stage` - specify the step (*rmod, consensus, coords, orfs* or *flanks*) at which the pipeline should finish (default - *flanks*)

`-k, --check-types` - print out all the types of MGE found in the RepeatModeler output (e.g. `./consensi.fa.classified`).

`-t, --threads` - set number of threads (default - all available CPUs).

`-C, --censor` - provide a path to CENSOR classification results (HTML file or URL).

`-o, --ori` - if specified MGERT will use the `*.ori` file to fetch the coordinates instead of `*_rm.out` file.

`-m, --merge` - merge all hits within that number of bp into a single entry. Default 0 bp (i.e. no merge).

`-e, --e-value` - set expectation value (E) for RPS-BLAST. Default 0.01.

`-c, --start-codon` - ORF start codon to use. 0 = 'ATG' only; 1 = 'ATG' and alternative initiation codons; 2 = any sense codon; Default 0.

`-l, --min-length` - set minimum length of ORF to be reported, default 1000 bp.

`-s, --strand` - output ORFs on specified strand only (e.g. plus/minus/both). Default 'plus'.

`-le, --left-end` - set length of left (5') flanking region. Default 0 bp.

`-re, -right-end` - set length of right (3') flanking region. Default 0 bp (if both *le* and *re* are set to 0, flanks mode will be omitted).

`-L, --rm-library` - specify a path to a library for RepeatMasker (in FASTA format). Use with `-f consensus` only. 
If consensus library is not specified, it will be compiled from RepeatModeler output automatically.

`-rm, --rm-table` - specify repeat masker table to use (with `*_rm.out` or `*.bed` extension). Use with `-f coords` option only.

`-sq, --sequence` - provide a path to a file of sequences where to look for conservative domains. Use with `-f orf` option only.

`-v, --version` - show program's version number and exit.
