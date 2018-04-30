import os
import argparse
from pathlib import Path
import parmed as pmd

def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.

    Parameters
    ----------
    x : str
        Candidate file path

    Returns
    -------
    Path
        Validated path

    """

    if not os.path.isfile(x):
        # ArgumentTypeError gives a rejection message of the form:
        # error: argument input: <passed error message>
        if os.path.exists(x):
            raise argparse.ArgumentTypeError("{0} is not a file".format(x))
        else:
            raise argparse.ArgumentTypeError("{0} does not exist".format(x))

    return Path(x)


def str2bool(x):
    """
    'Type' for argparse - converts variants of (t)rue/(false), (y)es/(n)o and
    1/0 to True/False.

    Parameters
    ----------
    x : str
        Commandline choice

    Returns
    -------
    Bool
        Converted from string

    """

    if x.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif x.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def check_prmtop(top_filename):
    """
    'Type' for argparse - checks that file exists and is a valid Amber prmtop

    Parameters
    ----------
    top_filename : str
        Candidate file path

    Returns
    -------
    Path
        Validated path

    """

    # Check file exists (and convert string to Path object)
    top_path = extant_file(top_filename)

    try:
        # Use parmed to check file is a valid Amber topology
        pmd.amber.AmberParm(top_filename)
    except:
        raise argparse.ArgumentTypeError('{:s} is not a valid prmtop file.'.format(top_filename))

    return top_path
