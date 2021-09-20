"Module for working with collections of MultipleSequenceAlignment objects."

from __future__ import absolute_import

UNDERLINE = "\033[4m"
RESET = "\033[0m"

class Supermatrix(object):
    "Represent a collection of MultipleSequenceAlignment objects."
    def __init__(self):
        self._msas = []

    def __len__(self):
        """Returns the length of all MultipleSequenceAlignment objects
        concatenated.
        """
        supermatrix_len = 0
        for msa in self.msas:
            supermatrix_len += msa.alignment_len()
        return supermatrix_len

    def __nonzero__(self):
        return True

    def __bool__(self):
        return True

    @property
    def msas(self):
        "A list of MultipleSequenceAlignment objects."
        return self._msas

    @msas.setter
    def msas(self, value):
        self._msas = value

    def missing_data(self):
        """Returns the percent missing data within the output orthologous
        MultipleSequenceAlignment object within this Summary object. Missing
        data is calculated as if these alignments were to be concatenated into
        a supermatrix, meaning that OTUs missing from an alignment are also
        considered.
        """
        pct_missing = 0.0
        no_of_alignments = 0
        no_of_otus = len(list(self.otus()))

        for msa in self.msas:
            no_of_alignments += 1
            otus_missing = no_of_otus - len(list(msa.otus()))
            pct_missing += msa.missing_data(otus_missing)

        if no_of_alignments > 0:
            return round((pct_missing / no_of_alignments) * 100, 1)
        else:
            return 0

    def sequences(self):
        """Returns the total number of sequences within all
        MultipleSequenceAlignment in this Supermatrix object.
        """
        no_of_seqs = 0
        for msa in self.msas:
            no_of_seqs += len(msa)
        return no_of_seqs

    def shortest_sequence(self):
        "Returns the shortest sequence within all of the Log objects."
        shortest = None
        for msa in self.msas:
            for sequence in msa.sequences:
                if not shortest or shortest > len(sequence.ungapped()):
                    shortest = len(sequence.ungapped())
        return shortest

    def longest_sequence(self):
        "Returns the longest sequence within all of the Log objects."
        longest = None
        for msa in self.msas:
            for sequence in msa.sequences:
                if not longest or longest < len(sequence.ungapped()):
                    longest = len(sequence.ungapped())
        return longest

    def avg_sequences(self):
        """Returns the average number of sequences per MultipleSequenceAlignment
        object.
        """
        no_of_files = 0
        no_of_seqs = 0
        for msa in self.msas:
            no_of_files += 1
            no_of_seqs += len(msa.sequences)

        if no_of_files > 0:
            return int(no_of_seqs / float(no_of_files))
        else:
            return 0

    def avg_seq_len(self):
        "Returns the average sequence length of all MSAs combined."
        seq_lens = 0
        sequences = 0
        for msa in self.msas:
            for sequence in msa.sequences:
                sequences += 1
                seq_lens += len(sequence.ungapped())
        if sequences > 0:
            return int(seq_lens / sequences)
        else:
            return 0

    def otus(self):
        "Returns a set of all OTUs within this SuperMatrix object."
        otus_in_supermatrix = set()
        for msa in self.msas:
            otus_in_supermatrix.update(msa.otus())
        return otus_in_supermatrix

    def avg_otus(self):
        "Returns the average number of OTUs in this summary."
        otus_total = 0
        ortholog_count = 0

        for msa in self.msas:
            otus = set(msa.iter_otus())
            otus_total += len(otus)
            ortholog_count += 1

        if ortholog_count > 0:
            return int(otus_total / ortholog_count)
        else:
            return 0

    def report(self, title, output, tree_stats):
        "Print statistics for this SuperMatrix object."
        report = """{}Ortholog statistics:{}
  No. of alignments:                {:5d}
  No. of sequences:                 {:5d}
  No. of OTUs:                      {:5d}
  Avg no. of sequences / alignment: {:5d}
  Avg no. of OTUs / alignment:      {:5d}
  Avg sequence length (ungapped):   {:5d}
  Shortest sequence (ungapped):     {:5d}
  Longest sequence (ungapped):      {:5d}
  % missing data:                   {:5.2f}
  Concatenated alignment length:    {:5d}
  Min tree diameter:                {:5.2f}
  Max tree diameter:                {:5.2f}
  Avg tree diameter:                {:5.2f}
  Median tree diameter:             {:5.2f}
  """.format(
    UNDERLINE,
    RESET,
    len(self.msas),
    self.sequences(),
    len(list(self.otus())),
    self.avg_sequences(),
    self.avg_otus(),
    self.avg_seq_len(),
    self.shortest_sequence(),
    self.longest_sequence(),
    self.missing_data(),
    len(self),
    tree_stats["min"],
    tree_stats["max"],
    tree_stats["mean"],
    tree_stats["median"],
  )

        row = "{};{};{};{};{};{};{};{};{};{};{};{};{};{};{}\n".format(
            title,
            len(self.msas),
            self.sequences(),
            len(self.otus()),
            self.avg_sequences(),
            self.avg_otus(),
            self.avg_seq_len(),
            self.shortest_sequence(),
            self.longest_sequence(),
            self.missing_data(),
            len(self),
            round(tree_stats["min"], 2),
            round(tree_stats["max"], 2),
            round(tree_stats["mean"], 2),
            round(tree_stats["median"], 2)
        )

        with open(output, "a") as msas_stats_file:
            msas_stats_file.write(row)

        print()
        print(report)
        print('Wrote ortholog statistics to {}'.format(output))
        return report
