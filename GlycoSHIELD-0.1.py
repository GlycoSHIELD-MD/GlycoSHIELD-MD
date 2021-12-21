__author__ = "Matt Sikora"
__copyright__ = "Copyright 2021, GlycoSHIELD project"
__credits__ = ["Cyril Hanus"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Matt Sikora"
__email__ = "masikora@biophys.mpg.de"
__status__ = "Development"

# demonstration of the import of the newly created glycoshield module
import glycoshield
# call a function from the glycoshield module
glycoshield.say_ok()

import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis import align
from argparse import ArgumentParser

         
def main2(protpdb,protxtc,inputfile,zmax,zmin,threshold,mode,dryrun,shuffle_sugar,ignorewarn):
   # parse input on every sugar to be added
   inputlines = inputparser(inputfile)
   
   #~ threshold=1.0 # threshold for the distance search 
   
   # when limiting no of atoms - expand max sugar span by this much. Only atoms withing such radius from the sequon are considered for clash checking.
   sugar_expand_factor = 1.5
   # Load protein
   # single frame:
   if protxtc is None:
      uprot=mda.Universe(protpdb,trajectory=True)
   else:
      uprot=mda.Universe(protpdb,protxtc)


   if mode=="All":
      selprot=uprot.select_atoms("all")
      sugar_sel_suffix = "all"
   elif mode=="CG":
      selprot=uprot.select_atoms("name CA")
      sugar_sel_suffix = 'name O5 or ( resname ANE5 and name O6 )'
   protframe=0
   # Iterate over protein states
   for prottp in uprot.trajectory:
      
      # counter for "entropy"
      occupancies = []
      
      # Iterate over sugars:
      for protchain,resids_on_protein,resids_on_sugar,sugarpdb,sugarxtc,pdbout,xtcout in inputlines:
         # load sugar
         usugar=mda.Universe(sugarpdb,sugarxtc,in_memory=True)
         # select (part of) the protein that is connected to the sugar
         tripep = usugar.select_atoms('protein and resid {} and name N CA CO'.format(" ".join([str(i) for i in resids_on_sugar])))
         # select sequon on protein
         sequon = uprot.select_atoms('segid {} and resid {} and name N CA CO'.format(protchain, " ".join([str(i) for i in resids_on_protein])))

         # the selection has to have same length for protein and peptide with the sugar!    
         # Now select part of the peptide that should be transplanted onto protein. By default we do not transplant the part which gets selected.
         # This way if we say that sugarresid is 1,2,3 and we are on a tripeptide - then nothing will get transplanted. 
         # And if we have a long loop 1,...,40 and we select only 1,40, then 2,...,39 will still be transplanted.
         
         peptide_to_transplant = [i for i in usugar.select_atoms('protein').residues.resids if i not in resids_on_sugar]      
         #~ print resids_on_sugar,uprot.atoms.resids
         
         if len(peptide_to_transplant)==0:
            sugar_without_protein = usugar.select_atoms('not protein')
         else:
            sugar_without_protein = usugar.select_atoms('not protein or ( protein and resid {} )'.format(" ".join([str(i) for i in peptide_to_transplant])))
         
         # We should probably be more lenient with the distance check between the first atoms of the peptide and the protein, i.e. dont check residues 2 and 39 in the example above
         if len(resids_on_sugar) != len(usugar.select_atoms('protein').residues.resids):
            buffer=1 # skip 1 AA after the first fitting AA and one before the last
            flanked_resids_on_sugar = sorted(list(set(usugar.select_atoms('protein').residues.resids)-set(resids_on_sugar)))[buffer:-buffer]
            sugar_without_protein_flanked = usugar.select_atoms('not protein or ( protein and resid {} )'.format(" ".join([str(i) for i in flanked_resids_on_sugar])))
         else:
            sugar_without_protein_flanked = sugar_without_protein
         
         
         if mode == "CG":
            sugar_without_protein_flanked = sugar_without_protein_flanked.select_atoms(sugar_sel_suffix)
            
         # Narrow the selection (especially important for large systems)
         
         sugar_root = sequon.center_of_mass()
         # Get size of the sugar:
         sizex=np.max(sugar_without_protein_flanked.positions[:,0])-np.min(sugar_without_protein_flanked.positions[:,0])
         sizey=np.max(sugar_without_protein_flanked.positions[:,1])-np.min(sugar_without_protein_flanked.positions[:,1])
         sizez=np.max(sugar_without_protein_flanked.positions[:,2])-np.min(sugar_without_protein_flanked.positions[:,2])
         maxdiam = np.max([sizex,sizey,sizez])*sugar_expand_factor
         
         # Generate list of coordinates.
         reduced_prot_positions = selprot.atoms.positions[(selprot.atoms.positions[:,0]<(sugar_root[0]+maxdiam)) & (selprot.atoms.positions[:,0]>(sugar_root[0]-maxdiam)) & 
                                                        (selprot.atoms.positions[:,1]<(sugar_root[1]+maxdiam)) & (selprot.atoms.positions[:,1]>(sugar_root[1]-maxdiam)) &
                                                        (selprot.atoms.positions[:,2]<(sugar_root[2]+maxdiam)) & (selprot.atoms.positions[:,2]>(sugar_root[2]-maxdiam))]
            
         # Placeholder for a reconstructed glycoprotein
         glycoprotein = mda.Merge(uprot.atoms,sugar_without_protein)
              
         coordinates = []
         # iterate over the traj to remove clashing frames. Can be sequential or random
         sugarframes = range(usugar.trajectory.n_frames)
         
         # randomize sugar conformers
         if shuffle_sugar:
            np.random.shuffle(list(sugarframes))
            
         for n_tp in sugarframes:
            _=usugar.trajectory[n_tp]
            #~ print tp.time
            # Align on the protein
            align.alignto(tripep,sequon)

            # check for overlap with a membrane: remove all conformers that have atoms below Zmin
            min=np.min(sugar_without_protein_flanked.positions[:,2])
            max=np.max(sugar_without_protein_flanked.positions[:,2])
            
            if ( not min<zmin ) and ( not max>zmax ) :
               
               # check for the overlap
               if threshold !=0:
                  distarr=distance_array(reduced_prot_positions,sugar_without_protein_flanked.positions)
                  if np.min(distarr)<threshold:
                     # overlap - do nothing
                     pass
                  else:
                     # OK - merge two systems copying protein and adding sugar.
                     coordinates.append(np.vstack((uprot.atoms.positions,sugar_without_protein.positions)))
               else:
                  # Threshold 0 - go ahead
                  coordinates.append(np.vstack((uprot.atoms.positions,sugar_without_protein.positions)))
         coordinates=np.array(coordinates)
         
         # save counters
         #~ occupancies.append(len(coordinates)/float(usugar.trajectory.n_frames))
         occupancies.append(len(coordinates))
         
         
         if len(coordinates)==0 and protxtc is None:
            # no frames selected. If this is a protein trajectory just ignore this.
            if not ignorewarn:
               
               raise BaseException("No frames were accepted for residues {} in chain {}!\nChange threshold or check if the sequon location is accessible!".format(resids_on_protein,protchain))
            else:
               print("WARNING!!!\nNo frames were accepted for residues {} in chain {}!\nThis means the sequon location is not accessible given the threshold\nYou decided to ignore this warning so we carry on!".format(resids_on_protein,protchain))
               continue
         # If iterating over protein traj - make sure to modify coordinates.
         if protxtc is not None:
            xtcout = xtcout.replace('.xtc','_protein_frame_{}.xtc'.format(protframe))
         
         # skip writing xtc,pdb
         if not dryrun:
            glycoprotein.load_new(coordinates, format=MemoryReader)
            # Write new reference
            glycoprotein.atoms.write(pdbout)
            
            # Write new trajectory
            with mda.coordinates.XTC.XTCWriter(xtcout,n_atoms=glycoprotein.atoms.n_atoms) as w:
               for tp in glycoprotein.trajectory:
                  w.write(glycoprotein.atoms)
                  
         del usugar,glycoprotein
         
      # Main output printout
      print(protframe," ".join([str(i) for i in occupancies]))
      
      # Increment frame counter
      protframe+=1
      
def inputparser(fname):
   inputs=[]
   with open(fname) as f:
      for line in f:
         if "#" not in line:
            l=line.split()
            protchain=l[0]
            protsequonresid=[int(i) for i in l[1].split(",")]
            sugarsequonresid=[int(i) for i in l[2].split(",")]
            sugarpdb=l[3]
            sugarxtc=l[4]
            pdbout=l[5]
            xtcout=l[6]
            inputs.append([protchain,protsequonresid,sugarsequonresid,sugarpdb,sugarxtc,pdbout,xtcout])
   return inputs
   
   
if __name__=="__main__":
   parser=ArgumentParser()
   parser.add_argument('--protpdb', dest='protpdb', help='pdb file for protein')
   parser.add_argument('--protxtc', dest='protxtc', help='xtc file for protein, optional if one needs multiple states of the protein')
   parser.add_argument('--inputfile', dest='inputfile', help='Input file with each line representing a sequon')
   parser.add_argument('--threshold',dest='threshold', help='threshold in A for checking for the overlap between sugar and protein.')
   parser.add_argument('--mode',dest='mode', help='do we take all atoms "All" or only Calpha "CG" when calculating distances. CG useful for large proteins but remember to change threshold! A good value here would be 0.7 for All atom and 3.5 for CG')
   parser.add_argument('--zmin',dest='zmin', help='position of the membrane below which we do not accept the glycan',default=-110000)
   parser.add_argument('--zmax',dest='zmax', help='position of the membrane above which we do not accept the glycan',default=110000)
   parser.add_argument('--dryrun', dest='dryrun', action='store_true',help="Dry run, do not write pdb/xtc files")
   parser.add_argument('--no-dryrun', dest='dryrun', action='store_false')
   parser.add_argument('--shuffle-sugar', dest='shuffle_sugar', action='store_true',help='randomize sugar frames when grafting onto protein')
   parser.add_argument('--no-shuffle-sugar', dest='shuffle_sugar', action='store_false')
   parser.add_argument('--ignorewarn', dest='ignorewarn', action='store_true',help="Ignore cases when no sugar is implanted")
   parser.add_argument('--no-ignorewarn', dest='ignorewarn', action='store_false')   
   parser.set_defaults(dryrun=False)
   parser.set_defaults(ignorewarn=False)
   parser.set_defaults(protxtc=None)
   parser.set_defaults(shuffle_sugar=False)
   
   parser.set_defaults(mode="All")

   args = parser.parse_args()
   
   protpdb = args.protpdb
   protxtc = args.protxtc
   inputfile=args.inputfile
   threshold=float(args.threshold)
   zmax=float(args.zmax)
   zmin=float(args.zmin)
   mode=args.mode
   dryrun=args.dryrun
   shuffle_sugar=args.shuffle_sugar
   ignorewarn=args.ignorewarn
   
   main2(protpdb,protxtc,inputfile,zmax,zmin,threshold,mode,dryrun,shuffle_sugar,ignorewarn)
