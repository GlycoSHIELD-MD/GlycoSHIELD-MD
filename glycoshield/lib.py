__author__ = "Matt Sikora"
__copyright__ = "Copyright 2021, GlycoSHIELD project"
__credits__ = ["Cyril Hanus"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Matt Sikora"
__email__ = "masikora@biophys.mpg.de"
__status__ = "Development"

# here is where you would implement (move) the functionality of the 'glycoshield' Python modules
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import align
from argparse import ArgumentParser

      
def say_ok():
    print("OK!")
