#!/bin/bash

#python ../GlycoSHIELD-0.1.py --protpdb EC5.pdb --inputfile EC5_input --threshold 3.5 --mode CG --shuffle-sugar
python ../GlycoSHIELD.py --protpdb EC5.pdb --inputfile EC5_input --threshold 3.5 --mode CG --shuffle-sugar

#python ../GlycoTRAJ-0.1.py --maxframe 1172 --outname out_merged --pdblist A_463.pdb,A_492.pdb,A_533.pdb --xtclist A_463.xtc,A_492.xtc,A_533.xtc --chainlist A,A,A --reslist 463,492,533
python ../GlycoTRAJ.py --maxframe 1172 --outname out_merged --pdblist A_463.pdb,A_492.pdb,A_533.pdb --xtclist A_463.xtc,A_492.xtc,A_533.xtc --chainlist A,A,A --reslist 463,492,533

#python ../GlycoSASA-0.1.py --pdblist A_463.pdb,A_492.pdb,A_533.pdb --xtclist A_463.xtc,A_492.xtc,A_533.xtc --probelist 0.14,0.25 --endframe 1172 --plottrace
python ../GlycoSASA.py --pdblist A_463.pdb,A_492.pdb,A_533.pdb --xtclist A_463.xtc,A_492.xtc,A_533.xtc --probelist 0.14,0.25 --endframe 1172 --plottrace

