"Store a tree and alignment pair"

import os

class FilePair(object):
    """
    Stores a corresponding tree and alignment file pair.
    """
    def __init__(self, tree_file="", alignment_file=""):
        self._tree = str(tree_file)
        self._alignment = str(alignment_file)

    def __str__(self):
        return self.filename

    def __nonzero__(self):
        return True

    def __bool__(self):
        return True

    @property
    def tree(self):
        "The path to a Newick-formatted file."
        return self._tree

    @tree.setter
    def tree(self, value):
        self._tree = value

    @property
    def alignment(self):
        "The path to an alignment file in FASTA format."
        return self._alignment

    @alignment.setter
    def alignment(self, value):
        self._alignment = value

    def has_alignment(self):
        return self.alignment

    def has_tree(self):
        return self.tree

    def is_filled(self):
        return self.tree and self.alignment

    def dirname(self):
        if self.alignment:
            return os.path.dirname(self.alignment)
        elif self.tree:
            return os.path.dirname(self.tree)
        else:
            return ""

    def filename(self):
        if self.alignment:
            return os.path.basename(os.path.splitext(self.alignment)[0])
        elif self.tree:
            return os.path.basename(os.path.splitext(self.tree)[0])
        else:
            return ""

    def alignment_ext(self):
        if self.alignment:
            return os.path.basename(os.path.splitext(self.alignment)[0])
        else:
            return ""
