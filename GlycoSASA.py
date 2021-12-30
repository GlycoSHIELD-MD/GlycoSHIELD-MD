from glycoshield import glycosasa, GMXTEST
from argparse import ArgumentParser
__author__ = "Matt Sikora"
__copyright__ = "Copyright 2021, GlycoSHIELD project"
__credits__ = ["Cyril Hanus"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Matt Sikora"
__email__ = "masikora@biophys.mpg.de"
__status__ = "Development"


import matplotlib
matplotlib.use('Agg')
import warnings
from argparse import ArgumentParser
from glycoshield import glycosasa, GMXTEST


def run_glycosasa(pdbfiles, xtcfiles, plottrace, probes, ndots, mode, keepoutput, endframe):
    with warnings.catch_warnings():
        # warnings.simplefilter('ignore')
        glycosasa(pdbfiles, xtcfiles, plottrace, probes, ndots, mode, keepoutput, endframe)
        print("OK")


if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument('--pdblist', dest='pdblist', help='coma-separated list of pdb references [no spaces! test1.pdb,test2.pdb ]')
    parser.add_argument('--xtclist', dest='xtclist', help='coma-separated list of xtc trajs [no spaces! test1.xtc,test2.xtc ]')
    parser.add_argument('--probelist', dest='probelist', help='coma-separated list of probe sizes for SASA calculation, in nm. can be a single one')
    parser.add_argument('--plottrace', dest='plottrace', action='store_true', help='plot per-residue shielding? default is no')
    parser.add_argument('--no-plottrace', dest='plottrace', action='store_false', help='do not plot per-residue shielding')
    parser.add_argument('--keepoutput', dest='keepoutput', action='store_true', help="do not delete tmp files, default is no")
    parser.add_argument('--no-keepoutput', dest='keepoutput', action='store_false', help="delete tmp files")
    parser.add_argument('--ndots', dest='ndots', help='number of dots used per atom by the Shrake-Rupley algorithm (see Gromacs documentation for details)', default=15)
    parser.add_argument('--mode', dest='mode', help='mode = avg or max, report max or avg Delta_sasa', default="max")
    parser.add_argument('--endframe', dest='endframe', help='last frame to read from the glycan trajectories. Assumes 1 ps time step', default=0)

    args = parser.parse_args()

    # Test if GROMACS can be run locally.
    GMXTEST()

    pdbfiles = [i for i in args.pdblist.split(",")]
    xtcfiles = [i for i in args.xtclist.split(",")]
    probes = [float(i) for i in args.probelist.split(",")]
    plottrace = args.plottrace
    if plottrace == "None":
        plottrace = None
    ndots = args.ndots
    if ndots == "None":
        ndots = None
    mode = args.mode
    parser.set_defaults(keepoutput=False)
    parser.set_defaults(plottrace=False)
    keepoutput = args.keepoutput
    endframe = args.endframe

    run_glycosasa(pdbfiles, xtcfiles, plottrace, probes, ndots, mode, keepoutput, endframe)
