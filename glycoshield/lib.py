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
from tqdm import tqdm
      
def say_ok():
    print("OK!")


class glycoshield:
   def __init__(self,protpdb,protxtc,inputfile,threshold=3.5,mode="CG",zmin=None,zmax=None,dryrun=False,shuffle_sugar=True,ignorewarn=False):
      
      # set vars
      
      self.protpdb = protpdb
      self.protxtc = protxtc      
      self.inputfile = inputfile
      self.threshold = threshold
      self.mode = mode
      if zmin is None:
         self.zmin=-100000
      else:
         self.zmin=zmin
      if zmax is None:
         self.zmax=100000
      else:
         self.zmax=zmax
      
      self.dryrun = dryrun
      self.shuffle_sugar = shuffle_sugar
      self.ignorewarn = ignorewarn
      
      # parse input on every sugar to be added
      self.inputlines = self.inputparser()
      
      # when limiting no of atoms - expand max sugar span by this much. Only atoms withing such radius from the sequon are considered for clash checking.
      self.sugar_expand_factor = 1.5
      
      
      # Load structure
      # single frame:
      if self.protxtc is None:
         self.uprot=mda.Universe(self.protpdb,trajectory=True)
      else:
         self.uprot=mda.Universe(self.protpdb,self.protxtc)


      if self.mode=="All":
         self.selprot=self.uprot.select_atoms("all")
         self.sugar_sel_suffix = "all"
      elif self.mode=="CG":
         self.selprot=self.uprot.select_atoms("name CA")
         self.sugar_sel_suffix = 'name O5 or ( resname ANE5 and name O6 )'
      
   def run(self):
      # Iterate over protein states
      protframe = 0
      for prottp in self.uprot.trajectory:
         # counter for "entropy"
         occupancies = []
         
         # Iterate over sugars:
         for protchain,resids_on_protein,resids_on_sugar,sugarpdb,sugarxtc,pdbout,xtcout in self.inputlines:
            # load sugar
            self.usugar=mda.Universe(sugarpdb,sugarxtc,in_memory=True)
            # select (part of) the protein that is connected to the sugar
            tripep = self.usugar.select_atoms('protein and resid {} and name N CA CO'.format(" ".join([str(i) for i in resids_on_sugar])))
            # select sequon on protein
            sequon = self.uprot.select_atoms('segid {} and resid {} and name N CA CO'.format(protchain, " ".join([str(i) for i in resids_on_protein])))


            # the selection has to have same length for protein and peptide with the sugar!    
            # Now select part of the peptide that should be transplanted onto protein. By default we do not transplant the part which gets selected.
            # This way if we say that sugarresid is 1,2,3 and we are on a tripeptide - then nothing will get transplanted. 
            # And if we have a long loop 1,...,40 and we select only 1,40, then 2,...,39 will still be transplanted.
            
            # Loop grafting is experimental and not documented.
            
            peptide_to_transplant = [i for i in self.usugar.select_atoms('protein').residues.resids if i not in resids_on_sugar]      
            
            if len(peptide_to_transplant)==0:
               self.sugar_without_protein = self.usugar.select_atoms('not protein')
            else:
               self.sugar_without_protein = self.usugar.select_atoms('not protein or ( protein and resid {} )'.format(" ".join([str(i) for i in peptide_to_transplant])))
            
            # We should probably be more lenient with the distance check between the first atoms of the peptide and the protein, i.e. dont check residues 2 and 39 in the example above
            if len(resids_on_sugar) != len(self.usugar.select_atoms('protein').residues.resids):
               buffer=1 # skip 1 AA after the first fitting AA and one before the last
               flanked_resids_on_sugar = sorted(list(set(self.usugar.select_atoms('protein').residues.resids)-set(resids_on_sugar)))[buffer:-buffer]
               self.sugar_without_protein_flanked = self.usugar.select_atoms('not protein or ( protein and resid {} )'.format(" ".join([str(i) for i in flanked_resids_on_sugar])))
            else:
               self.sugar_without_protein_flanked = self.sugar_without_protein
            
            if self.mode == "CG":
               self.sugar_without_protein_flanked = self.sugar_without_protein_flanked.select_atoms(self.sugar_sel_suffix)
               
            
            # Initialise a list of sugar frames
            sugarframes=self.get_sugar_frames()

            # Narrow the selection (especially important for large systems)
            self.grid_search_around_sugar(sequon)

            coordinates = []
            # scan all conformations for clashes and append non-clashing ones to the coordinates:
            for n_tp in tqdm(sugarframes,desc="Structure {} residue {}".format(protframe,resids_on_protein[1])):
               
               self.process_sugar_conformer(n_tp,tripep,sequon,coordinates)
            coordinates=np.array(coordinates)
            
            
            
            
            # save counters
            occupancies.append(len(coordinates))
            
            
            if len(coordinates)==0 and self.protxtc is None:
               # no frames selected. If this is a protein trajectory just ignore this.
               if not self.ignorewarn:
                  
                  raise BaseException("No frames were accepted for residues {} in chain {}!\nChange threshold or check if the sequon location is accessible!".format(resids_on_protein,protchain))
               else:
                  print("WARNING!!!\nNo frames were accepted for residues {} in chain {}!\nThis means the sequon location is not accessible given the threshold\nYou decided to ignore this warning so we carry on!".format(resids_on_protein,protchain))
                  continue
            
            
            # Write output
            self.write_output(xtcout,pdbout,coordinates)
            
            # Cleanup
            del self.usugar
            
         # Main output printout [COULD BE REMOVED COMPLETELY AND INSTEAD SENT TO A VARIABLE]
         print(protframe," ".join([str(i) for i in occupancies]))
         
         # Increment frame counter
         protframe+=1
         
         
   def write_output(self,xtcout,pdbout,coordinates):
      """
      Save if not a dry run. Fixed output format if protein is a multiple frame trajectory
      
      [COULD ALREADY HERE TAKE OUTPUT OPTIONS .pdb OR .xtc AND A MAXIMUM NUMBER OF FRAMES TO WRITE?
      THE PDB TRAJECTORY WILL BE NEEDED FOR OUTPUT DOWNLOADABLE BY USERS, BUT CANNOT BE TOO BIG]
      """
      
      # If iterating over protein traj - make sure to modify coordinates.
      if self.protxtc is not None:
         xtcout = xtcout.replace('.xtc','_protein_frame_{}.xtc'.format(protframe))

      # if dry run - skip writing xtc,pdb
      if not self.dryrun:
         # Placeholder for a reconstructed glycoprotein
         glycoprotein = mda.Merge(self.uprot.atoms,self.sugar_without_protein)
         
         # fill it with the coordinates
         glycoprotein.load_new(coordinates, format=MemoryReader)
         # Write new reference pdb
         glycoprotein.atoms.write(pdbout)
         
         # Write new trajectory
         with mda.coordinates.XTC.XTCWriter(xtcout,n_atoms=glycoprotein.atoms.n_atoms) as w:
            for tp in glycoprotein.trajectory:
               w.write(glycoprotein.atoms)
               
      del glycoprotein
         
   def get_sugar_frames(self):
      # iterate over the traj to remove clashing frames. Can be sequential or random
      sugarframes = range(self.usugar.trajectory.n_frames)
      
      # randomize sugar conformers
      if self.shuffle_sugar:
         np.random.shuffle(list(sugarframes))
      return sugarframes         
         
   def process_sugar_conformer(self,n_tp,tripep,sequon,coordinates):
      """
      Process a single glycan frame. Check for clashes and whether it crosses zmin/zmax. If so, reject. Otherwise - add to 
      accepted
      """
      _=self.usugar.trajectory[n_tp]

      # Align on the protein
      align.alignto(tripep,sequon)

      # check for overlap with a membrane: remove all conformers that have atoms below Zmin
      min=np.min(self.sugar_without_protein_flanked.positions[:,2])
      max=np.max(self.sugar_without_protein_flanked.positions[:,2])
      
      if ( not min<self.zmin ) and ( not max>self.zmax ) :
         
         # check for the overlap
         if self.threshold !=0:
            distarr=distance_array(self.reduced_prot_positions,self.sugar_without_protein_flanked.positions)
            if np.min(distarr)<self.threshold:
               # overlap - do nothing
               pass
            else:
               # OK - merge two systems copying protein and adding sugar.
               coordinates.append(np.vstack((self.uprot.atoms.positions,self.sugar_without_protein.positions)))
         else:
            # Threshold 0 - go ahead
            coordinates.append(np.vstack((self.uprot.atoms.positions,self.sugar_without_protein.positions)))
            
            
            
   def grid_search_around_sugar(self,sequon):
      sugar_root = sequon.center_of_mass()
      # Get size of the sugar:
      sizex=np.max(self.sugar_without_protein_flanked.positions[:,0])-np.min(self.sugar_without_protein_flanked.positions[:,0])
      sizey=np.max(self.sugar_without_protein_flanked.positions[:,1])-np.min(self.sugar_without_protein_flanked.positions[:,1])
      sizez=np.max(self.sugar_without_protein_flanked.positions[:,2])-np.min(self.sugar_without_protein_flanked.positions[:,2])
      maxdiam = np.max([sizex,sizey,sizez])*self.sugar_expand_factor
      
      # Generate list of coordinates.
      self.reduced_prot_positions = self.selprot.atoms.positions[(self.selprot.atoms.positions[:,0]<(sugar_root[0]+maxdiam)) & (self.selprot.atoms.positions[:,0]>(sugar_root[0]-maxdiam)) & 
                                                     (self.selprot.atoms.positions[:,1]<(sugar_root[1]+maxdiam)) & (self.selprot.atoms.positions[:,1]>(sugar_root[1]-maxdiam)) &
                                                     (self.selprot.atoms.positions[:,2]<(sugar_root[2]+maxdiam)) & (self.selprot.atoms.positions[:,2]>(sugar_root[2]-maxdiam))]
            
   def inputparser(self):
      inputs=[]
      with open(self.inputfile) as f:
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
         


      

   