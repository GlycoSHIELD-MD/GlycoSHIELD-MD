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
from glycoshield import glycoshield


def run_glycoshield(protpdb, protxtc, inputfile, zmax, zmin, threshold, mode, dryrun, shuffle_sugar, ignorewarn,skip, path):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        gs = glycoshield(protpdb=protpdb, protxtc=protxtc, inputfile=inputfile,skip=skip, path=path)
        occupancy_single = gs.run()
        print(occupancy_single)
        print("OK")

def check_positive(value):
    """FROM https://stackoverflow.com/questions/64980270/how-to-allow-only-positive-integer-using-argparse"""
    try:
        value = int(value)
        if value <= 0:
            raise argparse.ArgumentTypeError("{} is not a positive integer".format(value))
    except ValueError:
        raise Exception("{} is not an integer".format(value))
    return value

if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument('--protpdb', dest='protpdb', help='pdb file for protein',required=True)
    parser.add_argument('--protxtc', dest='protxtc', help='xtc file for protein, optional if one needs multiple states of the protein',required=False)
    parser.add_argument('--inputfile', dest='inputfile', help='Input file with each line representing a sequon',required=True)
    parser.add_argument('--threshold', dest='threshold', help='threshold in A for checking for the overlap between sugar and protein.')
    parser.add_argument('--mode', dest='mode', help='do we take all atoms "All" or only Calpha "CG" when calculating distances. CG useful for large proteins but remember to change threshold! A good value here would be 0.7 for All atom and 3.5 for CG')
    parser.add_argument('--zmin', dest='zmin', help='position of the membrane below which we do not accept the glycan', default=-110000)
    parser.add_argument('--zmax', dest='zmax', help='position of the membrane above which we do not accept the glycan', default=110000)
    parser.add_argument('--dryrun', dest='dryrun', action='store_true', help="Dry run, do not write pdb/xtc files")
    parser.add_argument('--no-dryrun', dest='dryrun', action='store_false')
    parser.add_argument('--shuffle-sugar', dest='shuffle_sugar', action='store_true', help='randomize sugar frames when grafting onto protein')
    parser.add_argument('--no-shuffle-sugar', dest='shuffle_sugar', action='store_false')
    parser.add_argument('--ignorewarn', dest='ignorewarn', action='store_true', help="Ignore cases when no sugar is implanted")
    parser.add_argument('--no-ignorewarn', dest='ignorewarn', action='store_false')
    parser.add_argument('--skip', dest='skip', help='skip frames from the glycan trajectory, default = 1', default=1, type=check_positive)
    parser.add_argument('--path', dest='path', help='path for output files, default is current directory', default="./")
    parser.set_defaults(dryrun=False)
    parser.set_defaults(ignorewarn=False)
    parser.set_defaults(protxtc=None)
    parser.set_defaults(shuffle_sugar=False)
    parser.set_defaults(mode="All")

    args = parser.parse_args()

    protpdb = args.protpdb
    protxtc = args.protxtc
    inputfile = args.inputfile
    threshold = float(args.threshold)
    zmax = float(args.zmax)
    zmin = float(args.zmin)
    mode = args.mode
    dryrun = args.dryrun
    shuffle_sugar = args.shuffle_sugar
    ignorewarn = args.ignorewarn
    skip = int(args.skip)
    path= args.path

    run_glycoshield(protpdb, protxtc, inputfile, zmax, zmin, threshold, mode, dryrun, shuffle_sugar, ignorewarn,skip, path)
