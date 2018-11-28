# protTrace - A simulation based framework to estimate the evolutionary traceability of protein.
[![language: Python](https://img.shields.io/badge/language-Python-blue.svg?style=flat)](https://www.python.org/)
[![presented at: GCB2018](https://img.shields.io/badge/presented%20at-GCB2018-green.svg?style=flat)](http://gcb2018.de/)
[![published in: BioRxiv](https://img.shields.io/badge/published%20in-BioRxiv-ff69b4.svg?style=flat)](https://doi.org/10.1101/302109)
[![license: GPL-3.0](https://img.shields.io/badge/license-GNU--GPL3.0-lightgrey.svg)](https://opensource.org/licenses/GPL-3.0)

# Table of Contents
* [Scientific Context](#scientific-context)
* [Workflow](#Workflow)
* [Installation &amp; Usage](#installation--usage)
   * [protTrace &amp; Accessory Software](#prottrace--accessory-software)
   * [Configuring protTrace](#configuring-prottrace)
   * [Calling protTrace](#calling-prottrace)
* [Input Data](#input-data)
* [Test Run](#test-run)
* [WIKI](#wiki)
* [Bugs](#bugs)
* [Acknowledgements](#acknowledgements)
* [Code of Conduct &amp; License](#code-of-conduct--license)
* [Contact](#contact)
# Scientific context
*ProtTrace* is a simulation based approach to assess for a protein, the seed, over what evolutionary distances its orthologs can 
be found by means of sharing a significant sequence similarity. By doing so, it helps to differentiate between the true absence 
of an ortholog in a given species, and its non-detection due to a limited search sensitivity. *ProtTrace* was presented 2018 at the German Conference on Bioinformatics (GCB). The high resolution PDF of the corresponding poster is available from [HERE](https://github.com/BIONF/protTrace/wiki/images/Poster-ProtTrace.v2.pdf).
![Add Text](https://github.com/BIONF/protTrace/wiki/images/Poster-ProtTrace.v2.png "The evolutionary traceability of a protein")


# Workflow
The workflow of protTrace to infer the evolutionary traceability of a seed protein is shown in the figure below (mouse over to see details). It consists of three main steps
  1. **Parameterization:** The compilation of an orthologous group for this protein. In the standard setting, OMA orthologous groups are used. The sequences in the ortholog group are then used to infer the parameters of substitution and the insertion- and deletion process. 
  1. **Traceability calculation:** The in-silico evolution of the seed protein using the simulation software [REvolver](https://academic.oup.com/mbe/article/29/9/2133/1074669), and the determination of the traceability curve.
  1. **Visualization:** The inference of the traceability index for the protein in 233 species from all domains of life, and the generation of a colored tree.
A high resolution PDF of the image is available [HERE](https://github.com/BIONF/protTrace/wiki/images/Workflow-ProtTrace.v1.cap.pdf).

![Alt Text](https://github.com/BIONF/protTrace/wiki/images/Workflow-ProtTrace.v1.png "Workflow of the protTrace analysis** The workflow is distinguished into the categories Parameterization, Traceability calculation, and visualization. Boxes in green denote input files, boxes in orange represent meta-information, which is generated in the course of the analysis, and yellow boxes indicate output files that are generated as a result of an analysis. Arrows represent individual analysis steps, where the arrow style indicates whether the analysis step is obligatory (solid), or optional (dashed). Analysis steps that require the calling of external programs are indicated by the program name next to the corresponding arrow. Obligatory dependencies on 3d party software are represented by bold face black program names, those that are optional are indicated by grey font color.") 

# Installation & Usage
*protTrace* is written in Python 2.7, some helper scripts in Perl and R. Find below a the 3rd party software that is required by protTrace:
  * The ProtTrace package contains scripts written in different languages. In order to run ProtTrace you need the following resources:
  * Python v2.7.13 or higher. **Note, ProtTrace will not run under Python 3**
     * Install also the [DendroPy module ](https://www.dendropy.org/) (can be done via [Conda](https://github.com/BIONF/protTrace/wiki/EnvironmentSetUp).
  *. Perl v5 or higher including the following modules
     * Getopt::Long
     * List::Util
     * LWP::Simple
  * Java v1.7 or higher
  * R v3 or higher
  * [wget](https://www.gnu.org/software/wget/)
  
## protTrace & Accessory Software

| Program name | Version | Description | Mandatory | BioConda |
|------------- | ------- | ----------- | --------- | -------- |
|[MAFFT](https://mafft.cbrc.jp/alignment/software/) |v6 or higher|Multiple Sequence alignment|yes|[yes](https://bioconda.github.io/recipes/mafft/README.html)|
|[NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)|v2.7 or higher|Sequence similarity based search|yes|[yes](https://bioconda.github.io/recipes/blast/README.html)|
|[HMMER](http://hmmer.org/)|3.2 or higher|Sequence similarity based search using Hidden Markov Mode|yes|[yes](https://bioconda.github.io/recipes/hmmer/README.html)|
|[IQTREE](http://www.iqtree.org/)|1.6.7.1 or higher|Phylogenetic tree reconstruction|yes|[yes](https://anaconda.org/bioconda/iqtree)|
|[HaMStR OneSeq](https://github.com/BIONF/hamstr)|v1 or higher|targeted ortholog search|no|no|

For the start, we suggest to omit the optional use of HaMStR, since the use of this software comes along with some strict naming conventions.

Once that is out of the way (we suggest to use the [conda package management system](https://github.com/BIONF/protTrace/wiki/EnvironmentSetUp) for this) you can just clone this repository to get a copy of *PhyloProfile*.

```
git clone https://github.com/BIONF/protTrace
```

## Configuring protTrace
To configure protTrace simply move into the protTrace directory and run the [configure script](https://github.com/BIONF/protTrace/wiki#the-configuration-script)

```
perl bin/create_conf.pl -name=prog.conf -getOMA -getPfam
```

This will check if all dependencies are existing, it will allow you to set all parameters required for the protTrace run, and eventually
will download the required data from the [OMA database](https://omabrowser.org) and from the [Pfam database](http://pfam.xfam.org/). 
     * If you are confident that you have this data already available, you can omit either or both of the options **-getOMA** and **-getPfam**. You will then have to tell protTrace via the *create_conf.pl* script
where this data is located. 
     * **Make sure** to adhere to the [formatting requirements for the OMA data](https://github.com/BIONF/protTrace/wiki/PreparingOMA),
     and that you ran *hmmpress* on the Pfam database.

Once everything is set, you are ready to run protTest

## Calling protTest
Enter the protTest directory and type
```
python bin/protTrace.py -h
```
this should obtain
```
USAGE:  protTrace.py -i <omaIdsFile> | -f <fastaSeqsFile> -c <configFile> [-h]
        -i              Text file containing protein OMA ids (1 id per line)
        -f              List of input protein sequences in fasta format
        -c              Configuration file for setting program's dependencies
```
# Input Data
*protTest* can use either OMA protein ids, or a protein sequence in fasta format as input

In `toy_example/` you can find two files, test.ids and test.fasta for performing a test run with protTrace.

We describe the input in the section [Test Run](https://github.com/BIONF/protTrace/wiki#test-run) of our [WIKI](https://github.com/BIONF/protTrace/wiki/home).

# Test Run
We provide in the directory *toy_example* two files for testing protTrace
  1. *test.ids*: This file contains the OMA protein id of a yeast protein [DIM1](https://omabrowser.org/oma/info/YEAST05874/). To run this test:
     1. create a config file **prot.conf** using the *create_conf.pl* script. We recommend to leave all values as default for the start
     1. place the config file into the directory *toy_example*
     1. enter the directory *toy_example* and run protTrace by typing
     ```
     python ../bin/protTrace.py -i test.id -c prot.conf
     ```
     The output that will be generated by this run is described in the [WIKI](https://github.com/BIONF/protTrace/wiki#oma-id-as-input)
  1. *test.fasta*: This file contains the protein sequence of human ZNT3. 
     1. create or modify the config file **prog.conf** using the *create_conf.pl* script. Make sure to set in the section 
     [General Options](https://github.com/BIONF/protTrace/wiki/Config-File#general-options) the entry **species** to **HUMAN**
     1. place the config file into the directory *toy_example*
     1. enter the directory *toy_example* and run protTrace by typing
     ```
     python ../bin/protTrace.py -f test.fasta -c prot.conf
     ```
     The output that will be generated by this run is described in the [WIKI](https://github.com/BIONF/protTrace/wiki#protein-sequence-as-input)

# WIKI
Read the [WIKI](https://github.com/BIONF/protTrace/wiki/home) to explore the functionality of protTrace.

# Bugs
Any bug reports or comments, suggestions are highly appreciated. Please open an issue on GitHub or be in touch via email.

# Acknowledgements
We would like to thank the members of [Ebersberger group](http://www.bio.uni-frankfurt.de/43045195/ak-ebersberger) for many valuable suggestions and ...bug reports :)

# Contributors
* [Arpit Jain](https://github.com/aj87)
* [Ingo Ebersberger](https://github.com/BIONF)

# License
This tool is released under [GNU-GPL3.0 license](https://github.com/BIONF/protTrace/blob/master/LICENSE).

# How-To Cite
Arpit Jain, Arndt von Haeseler, Ingo Ebersberger The evolutionary Traceability of protein (2018) [BioRxiv](https://doi.org/10.1101/302109) 

# Contact
Ingo Ebersberger
ebersberger@bio.uni-frankfurt.de

