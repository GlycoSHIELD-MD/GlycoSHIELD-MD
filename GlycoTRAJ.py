__author__ = "Matt Sikora"
__copyright__ = "Copyright 2021, GlycoSHIELD project"
__credits__ = ["Cyril Hanus"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Matt Sikora"
__email__ = "masikora@biophys.mpg.de"
__status__ = "Development"


import warnings
from argparse import ArgumentParser
from glycoshield import glycotraj


def run_glycotraj(maxframe, outname, pdblist, xtclist, chainlist, reslist):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        glycotraj(maxframe, outname, pdblist, xtclist, chainlist, reslist)
        print("OK")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--maxframe', dest='maxframe', help='number of frames to include for each glycan',required=True)
    parser.add_argument('--outname', dest='outname', help='root name of the output file',required=True)
    parser.add_argument('--pdblist', dest='pdblist', help='(input) list of coma-separated reference pdb files',required=True)
    parser.add_argument('--xtclist', dest='xtclist', help='(input) list of coma-separated trajectory files',required=True)
    parser.add_argument('--chainlist', dest='chainlist', help='(input) list of coma-separated chain descriptors where each glycan was attached',required=True)
    parser.add_argument('--reslist', dest='reslist', help='(input) list of coma-separated residue numbers where each glycan was attached',required=True)

    args = parser.parse_args()

    maxframe = int(args.maxframe)
    outname = args.outname
    pdblist = [i for i in args.pdblist.split(",")]
    xtclist = [i for i in args.xtclist.split(",")]
    chainlist = [i for i in args.chainlist.split(",")]
    reslist = [i for i in args.reslist.split(",")]

    run_glycotraj(maxframe, outname, pdblist, xtclist, chainlist, reslist)
