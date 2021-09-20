# Module for running the tree inference program FastTree 2

import pathlib
import subprocess
import tempfile

def run_fasttree(alignment_in, tree_out=None):
    '''
    Runs FastTree on the provided alignment_in and writes the resulting Newick
    tree file to tree_out. If no output file has been specified, write the
    results to a temporary file. Returns the path to the output Newick file.
    '''
    alignment_path = pathlib.Path(alignment_in)

    if not alignment_path.is_file():
        raise FileNotFoundError('provided path ', alignment_in, ' \
is not a file')

    if tree_out:
        tree_path = pathlib.Path(tree_out)
        if tree_path.is_file():
            raise FileExistsError('provided output path ', tree_out, ' \
already exists')
        tree_out_path = tree_out
    else:
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tree_out_path = tmp.name

    fasttree_cmd = ["fasttree", "-out", tree_out_path, "-slow", "-gamma",
                    alignment_path]

    try:
        subprocess.run(fasttree_cmd, capture_output=True)
    except:
        raise ProcessLookupError('Couldn\'t run FastTree, please ensure \
that \'fasttree\' is installed and in your path')

    return tree_out_path
