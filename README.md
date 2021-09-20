<img src="https://github.com/fethalen/orthoquoll/blob/master/orthoquoll_logo_700px.png" alt="ppp_logotype" width="200"/>

## About

This program evaluates the quality of one or more sets of orthologs that were
inferred using a tree-based orthology inference program such as
[PhyloPyPruner](https://github.com/fethalen/phylopypruner). In addition to the
usual metrics provided by the forementioned orthology inference program,
OrthoQuoll also calculates tree diameter statistics for all input orthologs.

## Installation

OrthoQuoll has two dependencies:
[FastTree](http://www.microbesonline.org/fasttree/) and
[MAFFT](https://mafft.cbrc.jp/alignment/software/). The easiest way to install
these is to install these via Anaconda. First, _if you don't have Anaconda
installed already_, follow
[these instructions](https://www.anaconda.com/products/individual). After
`conda` has been installed, simply put the following into your terminal:

```
$ conda install -c bioconda fasttree mafft
```

To run OrthoQuoll itself, clone this repository into your desired location,
`cd` into the base directory (`orthoquoll`) and run the program like so:

```
$ ./orthoquoll.py
```

In some cases, you may have to add the permission to execute the program first:

```
$ chmod +x orthoquoll.py
```

## Tutorial

The first thing you might want to do after you have installed the program is
to print the help menu by issuing the following command:

```
$ ./orthoquoll.py --help
usage: orthoquoll.py [-h] [--id STRING] [--realign] [--subdirs] [--no-trees] [--no-header] [--threads COUNT] [--overwrite] [--output PATH] PATH [PATH ...]

Extract statistics from a directory of alignments.

positional arguments:
  PATH             path to a directory that contains multiple FASTA files or one or more FASTA files

optional arguments:
  -h, --help       show this help message and exit
  --id STRING      give this supermatrix a custom name (default: unknown)
  --realign        realign all alignments using MAFFT's LINSI algorithm
  --subdirs        search for files in subdirectories
  --no-trees       do not infer phylogenetic trees and do not report tree diameter statistics (much faster)
  --no-header      if provided, do not include a header in the CSV output
  --threads COUNT  number of threads used for running MAFFT and FastTree (default: all available
  --overwrite      if provided, overwrite any existing file with the same output path
  --output PATH    path to the output file, pre-existing files are overwritten! (default: supermatrix_stats.csv)
```

The input to OrthoQuoll is a directory of FASTA files. The program expects that
each species is separated using a delimiter (either `|` or `@`). For example,
`Drosophila@16S` is a valid entry. There is an example directory included in
the base directory called `test_data`. The test data is made up of a small set
of orthologs from the animal group Lophotrochozoa that were inferred using
PhyloPyPruner.

In our example run, we will just utilize the base functionality of OrthoQuoll.
Move into the project's base directory and type the following:

```
$ ./orthoquoll.py test_data
OrthoQuoll 0.1.0
author: Felix Thalen <felix.thalen@uni-goettingen.de>

Generating trees using FastTree

Ortholog statistics:
  No. of alignments:                   10
  No. of sequences:                   819
  No. of OTUs:                         74
  Avg no. of sequences / alignment:    81
  Avg no. of OTUs / alignment:         56
  Avg sequence length (ungapped):     174
  Shortest sequence (ungapped):        52
  Longest sequence (ungapped):        232
  % missing data:                   26.60
  Concatenated alignment length:     1913
  Min tree diameter:                 1.21
  Max tree diameter:                 2.46
  Avg tree diameter:                 1.97
  Median tree diameter:              2.14

Wrote ortholog statistics to supermatrix_stats.csv

completed in 7.14 seconds
```

The output will also be saved to the file `supermatrix_stats.csv`. You can
change the name of the output file by using the flag `--output`. The results
written to this file is appended unless the flag `--overwrite` has been set.
You can use the flag `--id` to give your run a different name and you can use
`--no-header` to skip the header line in the output file. Here is an example
of what the `supermatrix_stats.csv` output file looks like:

```
id;alignments;sequences;otus;meanSequences;meanOtus;meanSeqLen;shortestSeq;longestSeq;pctMissingData;catAlignmentLen;minTreeDiameter;maxTreeDiameter;meanTreeDiameter;medianTreeDiameter
lophotrochozoa;10;819;74;81;56;174;52;232;26.6;1913;1.21;2.46;1.97;2.14
```

This file is most easily viewed using `column -t -s, supermatrix_stats.csv` or
by opening it in a spreadsheet program such as Excel or LibreOffice Calc.

© Animal Evolution and Biodiversity 2019
