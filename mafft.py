# Module for running the tree inference program FastTree 2

import pathlib
import subprocess
import tempfile

def run_mafft(path_in, path_out=None):
    '''
    Runs MAFFT, using the LINSI algorithm, on the provided alignment_path and
    outputs a file with the same filename as the alignment path but with the
    provided suffix added to the output. Returns the path to the realigned file.
    If no output file has been provided, write the results to a temporary file.
    '''
    alignment_in = pathlib.Path(path_in)

    if not alignment_in.is_file():
        raise FileNotFoundError('provided path ', alignment_in, ' \
is not a file')

    if path_out:
        alignment_out = pathlib.Path(path_out)
        if alignment_out.is_file():
            raise FileExistsError('provided output path ', alignment_out, ' \
already exists')
        alignment_out_path = alignment_out
    else:
        tmp = tempfile.NamedTemporaryFile(delete=False)
        alignment_out_path = tmp.name

    mafft_cmd = ['mafft', '--maxiterate', '1000', '--localpair', '--quiet',
                 alignment_in]

    try:
        with open(alignment_out_path, 'w') as output_file:
            subprocess.run(mafft_cmd, stdout=output_file)
    except:
        raise ProcessLookupError('Can\'t run MAFFT, please ensure that \
\'mafft\' is installed and is in your path')

    return alignment_out_path
