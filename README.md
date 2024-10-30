# Corer

### Sequence-based Pangenomic Core Detection

* *k*-mer-based detection of a pangenome's core
* Assembled genomes or reads as input
* Customized core definition
* Low demands in run time and memory

## Publication

Schulz, T., Wittler, R., Stoye, J.: [Sequence-based pangenomic core detection](https://doi.org/10.1016/j.isci.2022.104413). iScience. (2022)

## Table of Contents

* [Requirements](https://github.com/gi-bielefeld/corer#requirements)
* [Compilation](https://github.com/gi-bielefeld/corer#compilation)
* [Usage](https://github.com/gi-bielefeld/corer#usage)
* [Test data](https://github.com/gi-bielefeld/corer#test-data)
* [Evaluation workflow](https://github.com/gi-bielefeld/corer#Evaluation-workflow)
* [FAQ](https://github.com/gi-bielefeld/corer#faq)
* [Contact](https://github.com/gi-bielefeld/corer#contact)
* [Licenses](https://github.com/gi-bielefeld/corer#Licenses)

## Requirements

Corer identifies the core of a given pangenome represented as a **compacted, colored de Bruijn graph** using the API of [Bifrost](https://github.com/pmelsted/bifrost) (version 1.2.1 or higher). Apart from the requirements of Bifrost (c++ and cmake), there are no further strict dependencies.

A provided evaluation workflow requires [snakemake](https://snakemake.readthedocs.io/en/stable/), and the packages
[Biopython](https://biopython.org) and [matplotlib](https://matplotlib.org) to be installed on your system.

## Compilation

```
cd <corer_directory>/src
make
```

By default, the installation creates:
* a binary (*Corer*)

You may want to make the binary accessible via your *PATH* variable.

Please note the installation instructions regarding the default maximum *k*-mer size of Bifrost from its README.
E.g., if your Bifrost libraries have been compiled to support a *k*-mer length of up to 63, change Corer's 
makefile accordingly (add `-DMAX_KMER_SIZE=64` to CFLAGS).

If during the compilation, the Bifrost library files are not found, make sure that the corresponding folder is found as include path by the C++ compiler. You may have to add
`-I/usr/local/include` (with the corresponding folder) to CFLAGS in the makefile.

## Usage:

```
Corer
```

displays the command line interface:
```
Corer [-hs] [-q QUORUM] [-d DELTA] [-i Graph_File_Prefix] [-o Output_File_Prefix] [-t Nb_Threads]

Extracting a pangenome's core.

Required parameters:
   -i   --igraph  Input graph file prefix
   -o   --ograph  Output graph file prefix

Optional parameters with required argument:
   -q   --quorum   Absolute quorum defining the core (default is 90% of pangenome size)
   -d   --delta    Maximum distance between core k-mers (default is 50)
   -t   --threads  Number of threads (default is 1)

Optional parameters without argument:
   -s   --snippets  Output unitig core snippets to stdout
   -h   --help      Display this help message
```

### Examples

1. **Get some testing data**

   You may want to get some toy data first to try out the program.

   ```
   wget 'https://hgdownload.cse.ucsc.edu/goldenPath/eboVir3/bigZips/160sequences.tar.gz'
   tar xvzf 160sequences.tar.gz
   ```

   These commands let you download and unpack a small dataset of 160 individual samples of different Ebola virus subspecies.

2. **Create a pangenome graph**

   Corer requires a Bifrost graph to perform a core prediction. First, create an input file list for Bifrost, e.g., by running
   
   ```
   ls -l *.fa | tr -s ' ' | cut -d' ' -f9 > list.txt
   ```
   
   In order to build the actual graph use Bifrost.

   ```
   path/to/bifrost/bin/Bifrost build -r list.txt -o ebolaPangenome -c
   ```

2. **Predict the core**

   Once the graph has been built, Corer can predict its core. The command

   ```
   Corer -i ebolaPangenome -o ebolaCore -q 128 -d 60
   ```

   predicts a core for the ebola pangenome, where a *core *k*-mer* needs to occur in at least 128 genomes and the maximum distance between two core
   *k*-mers is 61.

## Test data

Test data can be downloaded from public databases.

1. **Prokaryotic data sets**

   Four prokaryotic pangenomes may be downloaded from [NCBI](https://www.ncbi.nlm.nih.gov). Accession numbers can be found in the directory
   *experiments*.
   
3. **Arabidopsis pangenome**

   Assemblies of 18 accessions of Arabidopsis thaliana can be downloaded [here](http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/fasta/MASKED/).

   Corresponding read data sets for 17 of the above assemblies are available [here](https://www.ebi.ac.uk/ena/browser/view/PRJEB2457?show=reads).

## Evaluation workflow

The directory *experiments* contains an evaluation workflow that may be executed using
[snakemake](https://snakemake.readthedocs.io/en/stable/).

For execution of the workflow, proceed as follows:

* Install all programs to be tested on your system.

* Install gene annotation pipelines ([prokka](https://github.com/tseemann/prokka) for prokaryotic and [Augustus](http://bioinf.uni-greifswald.de/augustus/) for eukaryotic data sets).

* Install all tools needed for subsequent analyses.

* Download the testing data.

* Place your testing data into the provided subdirectories and add the locations of program 
  binaries and other dependencies into the configuration file *experiments/config.yaml*:

  ```
  # PLEASE ADJUST THE FOLLOWING PARAMETERS --------------------------------------

  #Program binaries that shall be used
   corer_bin: "../src/Corer"
   panaroo_bin: "/path/to/panaroo/bin"
   sibeliaz_bin: "/path/to/sibeliaz/bin"
   #Program binaries of gene annotation softwares
   prokka_bin: "/path/to/prokka/bin"
   augustus_bin: "/path/to/augustus/bin"
   #Bifrost binary that shall be used
   bifrost_bin: "/path/to/bifrost/bin"
   #Blastn binary that shall be used
   blastn_bin: "/path/to/blastn/bin"
   #Makeblastdb binary that shall be used
   makeblastdb_bin: "/path/to/makeblastdb/bin"
   #Path to BUSCO analysis script
   BUSCOscript: "/path/to/BUSCO/script.py"
   #Brassicales BUSCO database directory
   BUSCObrasDB: "/path/to/db/dir"


  #------------------------------------------------------------------------------
  ...
  ```
  
* Change into directory *experiments* and run the workflow.

  ```
  cd experiments
  snakemake
  ```

## FAQ

We recommend to have a look at the [FAQs of Bifrost](https://github.com/pmelsted/bifrost#faq).


## Contact

For any question, feedback or problem, please feel free to file an issue on [GitHub](https://github.com/gi-bielefeld/corer) or [contact](mailto:pangenomics-service@cebitec.uni-bielefeld.de) the developers and we will get back to you as soon as possible.

Corer is provided as a service of the [German Network for Bioinformatics Infrastructure (de.NBI)](https://www.denbi.de/). We would appriciate if you would participate in the evaluation of Corer by completing this [very short survey](https://www.surveymonkey.de/r/denbi-service?sc=bigi&tool=corer).

## Licenses

* The hash function library xxHash is BSD licensed (https://github.com/Cyan4973/xxHash)

* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)

* The libdivide library is zlib licensed (https://github.com/ridiculousfish/libdivide)

* The kseq library is copyrighted by Heng Li and released
  under the MIT license (http://lh3lh3.users.sourceforge.net/kseq.shtml)

* The CRoaring library is Apache 2.0 licensed (https://github.com/RoaringBitmap/CRoaring)

* Bifrost is BSD-2 licensed (https://github.com/pmelsted/bifrost)

* Corer is GNU GPLv3 licensed [LICENSE](https://gitlab.ub.uni-bielefeld.de/gi/plast/blob/master/LICENSE)
