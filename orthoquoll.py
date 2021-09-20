#!/usr/bin/env python3

# Author: Felix Thalen

"Extract statistics from a directory of alignments."

from multiprocessing import Pool
from multiprocessing import cpu_count
import argparse
import sys
import os
import statistics
import fasta
import fasttree
import mafft
import newick
import supermatrix
import time
import pathlib

FASTA_EXTENSIONS = {".fa", ".fas", ".fasta", ".fna", ".faa", ".fsa", ".ffn",
                    ".frn"}
NEWICK_EXTENSIONS = {".nw", ".newick", ".tre", ".tree", ".treefile"}
ALL_EXTENSIONS = FASTA_EXTENSIONS | NEWICK_EXTENSIONS
HEADER = "id;alignments;sequences;otus;meanSequences;meanOtus;meanSeqLen;\
shortestSeq;longestSeq;pctMissingData;catAlignmentLen;minTreeDiameter;\
maxTreeDiameter;meanTreeDiameter;medianTreeDiameter\n"
VERSION = "0.1.0"
RESET = "\033[0m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"

def parse_args():
    'Parse the user-provided arguments.'
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('alignments',
                       metavar='PATH',
                       nargs='+',
                       type=str,
                       help='path to a directory that contains multiple FASTA \
                             files or one or more FASTA files')
    parser.add_argument('--id',
                        metavar='STRING',
                        type=str,
                        default='unknown',
                        help='give this supermatrix a custom name (default: \
                              unknown)')
    parser.add_argument('--realign',
                        action='store_true',
                        default=False,
                        help='realign all alignments using MAFFT\'s LINSI \
                              algorithm')
    parser.add_argument('--subdirs',
                        action='store_true',
                        default=False,
                        help='search for files in subdirectories')
    parser.add_argument('--no-trees',
                        action='store_true',
                        default=False,
                        help='do not infer phylogenetic trees and do not \
                              report tree diameter statistics (much faster)')
    parser.add_argument('--no-header',
                        action='store_true',
                        default=False,
                        help='if provided, do not include a header in the \
                              CSV output')
    parser.add_argument('--threads',
                        metavar='COUNT',
                        default=cpu_count(),
                        help='number of threads used for running MAFFT and \
                              FastTree (default: all available')
    parser.add_argument('--overwrite',
                        action='store_true',
                        default=False,
                        help='if provided, overwrite any existing file with \
                              the same output path')
    parser.add_argument('--output',
                       metavar='PATH',
                       help='path to the output file, pre-existing files are \
                             overwritten! (default: supermatrix_stats.csv)',
                       default='supermatrix_stats.csv')
    return parser.parse_args(args=None if sys.argv[1:] else ['--help'])


def is_fasta(path):
    "Returns True if the provided path is the path to a FASTA file"
    filename, extension = os.path.splitext(path)
    filename = filename.split(".")[0]
    extension = extension.lower()

    if extension in FASTA_EXTENSIONS:
        return True
    return False


def is_newick(path):
    "Returns True if the provided path is the path to a Newick file."
    filename, extension = os.path.splitext(path)
    filename = filename.split(".")[0]
    extension = extension.lower()

    if extension in NEWICK_EXTENSIONS:
        return True
    return False


def fasta_files(directory, extensions=FASTA_EXTENSIONS, subdirs=False):
    """
    Returns a list of all of the FASTA files within the provided directory.
    Searches for files within all of the directories' subdirectorie(s) if the
    subdirs argument is provided. Produces an error message and exits with
    an exit code of 1, in case no files were found.
    """
    files = []
    file_pairs = dict()
    rootdir = pathlib.Path(directory)
    directories = [rootdir]

    if subdirs:
        directories += [path for path in rootdir.iterdir() if path.is_dir()]

    for subdirectory in directories:
        files += [path for path in subdirectory.iterdir() if path.suffix in extensions]

    if len(files) < 1:
       print('no FASTA file pairs found in the provided path')
       sys.exit()

    return files


def tree_stats(tree_files):
    """
    Takes a list of tree files as an input and returns a dictionary containing
    the following tree statistics: max, min, mean, and median.
    """
    diameters = []

    for tree_file in tree_files:
        node = newick.read(tree_file)
        diameters.append(node.tree_diameter())

    stats = {
        "min": min(diameters),
        "max": max(diameters),
        "mean": statistics.mean(diameters),
        "median": statistics.median(diameters)
    }
    return stats

def main():
    """
    Generate statistics from multiple sequence alignments (MSAs) in the
    provided directory.
    """
    START_TIME = time.time()
    args = parse_args()

    print("{}OrthoQuoll{} {}".format(BOLD, RESET, VERSION))
    print("author: Felix Thalen <felix.thalen@uni-goettingen.de>\n")

    alignments = []
    for alignment in args.alignments:
        alignment_path = pathlib.Path(alignment)
        if alignment_path.is_file() and alignment_path.suffix in FASTA_EXTENSIONS:
            alignments.append(alignment)
        elif alignment_path.is_dir():
            alignments += fasta_files(alignment, FASTA_EXTENSIONS, args.subdirs)

    msas = []
    for file in alignments:
        # cat_msas.msas.append(fasta.read(file))
        msas.append(fasta.read(str(file)))

    # Realign everything if '--realign' is set
    if args.realign:
        print("Realigning provided alignments using MAFFT (L-INS-i algorithm)")
        realigned_files = []

        with Pool(processes=args.threads) as pool:
            for alignment in pool.imap_unordered(mafft.run_mafft, alignments):
                realigned_files.append(alignment)

        alignments = realigned_files

    if args.overwrite and os.path.exists(args.output):
        os.remove(args.output)

    # Generate tree diameter statistics (if not unset)
    if not args.no_trees:
        # Generate trees from all of the alignments
        print("Generating trees using FastTree")
        trees = []

        with Pool(processes=args.threads) as pool:
            for tree in pool.imap_unordered(fasttree.run_fasttree, alignments):
                trees.append(tree)

        if not args.no_header:
            with open(args.output, 'a') as msas_stats_file:
                msas_stats_file.write(HEADER)

        tree_diameter_stats = tree_stats(trees)
    else:
        print("Skipping tree inference and tree diameter statistics")
        tree_diameter_stats = {
            "min": 0.0,
            "max": 0.0,
            "mean": 0.0,
            "median": 0.0
        }

    cat_msas = supermatrix.Supermatrix()

    for file in alignments:
        cat_msas.msas.append(fasta.read(file))

    cat_msas.report(args.id, args.output, tree_diameter_stats)

    # cleanup temporary files
    if args.realign and not args.no_trees:
        for path in realigned_files + trees:
            os.remove(path)
    elif not args.no_trees:
        for path in trees:
            os.remove(path)

    print("\ncompleted in {} seconds".format(
        round(time.time() - START_TIME, 2)))


if __name__ == '__main__':
    main()
