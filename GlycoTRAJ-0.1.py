__author__ = "Matt Sikora"
__copyright__ = "Copyright 2021, GlycoSHIELD project"
__credits__ = ["Cyril Hanus"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Matt Sikora"
__email__ = "masikora@biophys.mpg.de"
__status__ = "Development"

import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.base import AnalysisFromFunction
import sys,glob
import numpy as np

from argparse import ArgumentParser


def main(maxframe,outname,pdblist,xtclist,chainlist,reslist):
   """
   Read in multiple trajectories each containing a static protein structure and a single mobile glycan.
   Script will merge these trajectories into one (i.e. with static protein and all glycans moving)
   and trim it to the --maxframe frames.
   Output (pdb+xtc files) will be saved as --outname
   input is defined by lists of comma separated structure and trajectory files, and chain and residues:
   --pdblist and corresponding --xtclist, --chainlist, --residuelist (same order)
   e.g.
   python GlycoTRAJ-0.1.py --maxframe 10 --outname ttt --pdblist S_463.pdb,S_492.pdb,S_533.pdb --xtclist S_463.xtc,S_492.xtc,S_533.xtc --chainlist S,S,S --reslist 463,492,533
   """

   chains = []
   resids = []

   # Iterate over structures and trajs for single glycans to extract subtrajs
   for pdb,xtc,chain,resid in zip(pdblist,xtclist,chainlist,reslist):
            
      newxtc = pdb.replace(".pdb","_glycan.xtc")
      newpdb = pdb.replace(".pdb","_glycan.pdb")

      u = mda.Universe(pdb,xtc)
      
      # Select glycan (assume this is the only non-protein object, can be made more specific, here we also exclude palmitoylated residues)
      sugarsel = u.select_atoms("not protein and not resname CYSP")
      
      # save sub-trajs of same length
      with mda.Writer(newxtc, sugarsel.n_atoms) as W:
         for ts in u.trajectory[:maxframe]:
           W.write(sugarsel)
      
      sugarsel.write(newpdb)
      chains.append(chain)
      resids.append(resid)
      
      
   chains=sorted(np.unique(chains))
   resids=sorted(np.unique(resids))



   sels=[]
   coords=[]
   
   # open the protein frame and add to the trajectories
   pdb = pdblist[0]
   xtc = xtclist[0]
   u = mda.Universe(pdb,xtc,in_memory=True)
   
   # again - this selection can be made more specific
   protein = u.select_atoms('protein')
   
   # trim
   coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(),protein).run().results
   coordinates = coordinates[:maxframe,:,:]
   sels.append(protein)
   coords.append(coordinates)



   # Now go second time over the files and merge, renumbering the glycans
   for chain in chains:
      # last residue is? we do it per chain
      lastres=protein.select_atoms("segid {}".format(chain)).resids[-1]

      for resid in resids:
         pdb = chain+"_"+str(resid)+"_glycan.pdb"

         xtc = pdb.replace(".pdb",".xtc")
         u1 = mda.Universe(pdb,xtc,in_memory=True)
         sugar = u1.atoms
         sugar.segments.segids = chain
         # update residues
         sugar.residues.resids += lastres
         
         lastres += len(sugar.residues)
         
         coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(),sugar).run().results
         coords.append(coordinates)

         sels.append(u1)
         

   # Concatenate and save
   coords = np.column_stack(coords)
   u2=mda.Merge(*[i.atoms for i in sels])
   u2.load_new(coords, format=MemoryReader)
   with mda.Writer(outname+'.xtc', u2.atoms.n_atoms) as W:
      for ts in u2.trajectory:
         W.write(u2.atoms)
         
   #Save a reference
   u2.atoms.write(outname+'.pdb')

if __name__=="__main__":
   parser=ArgumentParser()
   parser.add_argument('--maxframe', dest='maxframe', help='number of frames to include for each glycan')
   parser.add_argument('--outname', dest='outname', help='root name of the output file')
   parser.add_argument('--pdblist', dest='pdblist', help='(input) list of coma-separated reference pdb files')
   parser.add_argument('--xtclist', dest='xtclist', help='(input) list of coma-separated trajectory files')
   parser.add_argument('--chainlist', dest='chainlist', help='(input) list of coma-separated chain descriptors where each glycan was attached')
   parser.add_argument('--reslist', dest='reslist', help='(input) list of coma-separated residue numbers where each glycan was attached')
   

   args = parser.parse_args()
   
   maxframe=int(args.maxframe)
   outname=args.outname
   
   pdblist=[i for i in args.pdblist.split(",")]
   xtclist=[i for i in args.xtclist.split(",")]
   chainlist=[i for i in args.chainlist.split(",")]
   reslist=[i for i in args.reslist.split(",")]
   
   main(maxframe,outname,pdblist,xtclist,chainlist,reslist)