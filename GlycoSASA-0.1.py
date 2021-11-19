__author__ = "Matt Sikora"
__copyright__ = "Copyright 2021, GlycoSHIELD project"
__credits__ = ["Cyril Hanus"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Matt Sikora"
__email__ = "masikora@biophys.mpg.de"
__status__ = "Development"


import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.analysis.base import AnalysisFromFunction
import MDAnalysis.analysis.distances as distances
#from __future__ import print_function
import shlex
import os,subprocess,glob
from shutil import copyfile,move
import collections
from argparse import ArgumentParser

# Define modified amino acids (approximation: we turn them back to the original residues)
AMINO_ACID_VARIANTS = [
    ("ALA", ["AIB", "ALM", "AYA", "BNN", "CHG", "CSD", "DAL","DHA", "DNP", "FLA", "HAC", "MAA", "PRR", "TIH", "TPQ",]),
    ("ARG", ["AGM", "DAR", "HAR", "MMO", "ARM", "ARN", "HMR", "ACL"]),
    ("ASN", ["MEN", "DSG",]),
    ("ASP", ["DSP", "BHD", "2AS", "ASQ", "ASB", "ASA", "ASK", "ASH","ASL", "DAS"]),
    ("CYS", ["BCS", "BUC", "C5C", "C6C", "CCS", "CEA", "CME", "CSO", "CSP", "CSS", "CSW", "CSX", "CY1", "CY3", "CYG", "CYM", "CYP", "CYQ", "CYX", "DCY", "EFC", "OCS", "PEC", "PR3", "SCH", "SCS", "SCY", "SHC", "SMC", "SOC"]),
    ("GLU", ["GLH", "GGL", 'PCA', '5HP', 'DGL', 'CGU', 'GMA']),
    ("GLN", [("DGN", "X")]),
    ("GLY", ["GLZ", "SAR", 'NMC', 'GL3', 'GSC', 'MPQ', 'MSA']),
    ("HIS", ["DHI", "HID", "HIC", "HIE", "HIP", "HSD", "HSE","HSP", "MHS", "NEM", "NEP", "3AH"]),
    ("ILE", ['DIL', 'IIL']),
    ("LEU", ["BUG", "NLE", 'NLP', 'NLN', 'DLE', 'CLE', 'MLE']),
    ("LYS", ['LYM', 'ALY', 'LYZ', 'LYN', 'LLY', 'LLP', 'SHR', 'TRG', 'DLY', 'KCX']),
    ("MET", ["FME", "CXM", "OMT", "MSE"]),
    ("PHE", ["DAH", "DPN", "HPQ", "PHI", "PHL"]),
    ("PRO", ['DPR', 'HYP']),
    ("SER", ['OAS', 'MIS', 'SAC', 'SVA', 'SET', 'SEP', 'SEL', "DSN"]),
    ("THR", ["ALO", "BMT", "DTH", "THO", "TPO"]),
    ("TRP", ["DTR", "HTR", "LTR", "TPL", "TRO"]),
    ("TYR", ["DTY", "IYR", "PAQ", "PTR", "STY", "TYB", "TYM", "TYO", "TYQ", "TYS", "TYY"]),
    ("VAL", ["DIV", "DVA", "MVA"])
]

AMINO_ACID_VARIANTS_SUBSTITUTION = {}
for resline in AMINO_ACID_VARIANTS:
   res = resline[0]
   for modres in resline[1]:
      AMINO_ACID_VARIANTS_SUBSTITUTION[modres]=res



def get_SASA(pdbname,xtcname,index,groupname,selectgroups,outfile,outfileres,outfileatom,probe,ndots,endtime):
    """
    Use gromacs to calculate SASA.
    Tested with GMX 2019.x
    requires sourced GMXRC (gmx command in the PATH).
    
    selectgroups separated by ;
    
    The output consists of a text file that assigns a shielding value to each residue. These values multiplied by 100 are also 
    written into the b-factor column of the PDB file. Residues that are not visible to the probe at all are further marked by 0
    in the occupancy column of the PDB file. The output files are named: maxResidueSASA_probe_RPROBE.txt maxResidueSASA_probe_RPROBE.pdb, 
    where RPROBE is a probe radius.


   # there is also ndots argument, default to 15, regulating the number of dots in Shrake-Rupley algorithm. https://www.sciencedirect.com/science/article/abs/pii/0022283673900119?via%3Dihub
    """
    # Build a gromacs command
    cmd="gmx sasa -f {} -s {} -n {} -o {} -or {}".format(xtcname,pdbname,index,outfile,outfileres)
    cmd=cmd+" -oa {}".format(outfileatom)
    cmd=cmd+" -surface \'{}\' -xvg none -probe {} ".format(groupname,probe)
    cmd=cmd+" -output \'{}\' -ndots {}".format(selectgroups,ndots)
    cmd=cmd+" -e {}".format(endtime)
    print(cmd)
    
    # Execute
    command = shlex.split(cmd)
    proc = subprocess.Popen(command, stdout = subprocess.PIPE,stderr=subprocess.PIPE)
    out=[i for i in proc.stdout]
    err=[i for i in proc.stderr]
    proc.communicate()

    if proc.returncode!=0:
        for line in err:
            print(line)
            raise SystemExit("Error in Gromacs, stopping...")
            
    # Read in the gmx output files
    sasas=np.loadtxt(outfileres)
    sasaa=np.loadtxt(outfileatom)
    sasatime=np.loadtxt(outfile)
    return sasatime,sasas,sasaa,proc.returncode   

def GMXTEST():
    cmd='gmx --version --quiet | grep "GROMACS version"'
    command = shlex.split(cmd)
    try:
       proc = subprocess.Popen(command, stdout = subprocess.PIPE,stderr=subprocess.PIPE)
    except:
       raise BaseException("Gromacs not found or not working, stopping...")
    out=[i for i in proc.stdout]
    err=[i for i in proc.stderr]
    proc.communicate()

    if proc.returncode!=0:
        for line in err:
            print(line)
            raise SystemExit("Gromacs not found, stopping...")   
   
def main(pdbfiles,xtcfiles,plottrace,probes,ndots,mode,keepoutput,endtime):
   # Assumption is there is only a protein in the pdb file.
   # Everything else assumed to be a glycan.


   # Iterate over probe sizes
   for probe in probes:
      outrelativesasa=[]
      outrelativesasaaa=[]      
      
      
      # Iterate over trajectories for glycans
      # Load each glycan+protein trajectory
      for iglycan in range(len(xtcfiles)):
          
          pdb=pdbfiles[iglycan]
          xtc=xtcfiles[iglycan]
          u=mda.Universe(pdb,xtc)
          
          # Define some temporary file names.
          tmpindex='index0.ndx'
          tmpsel='prot; gly'
          tmpsys='system'
          baresel='prot'
          tmpsasa  = 'sasa1.xvg'
          tmpsasar = 'sasar.xvg'
          tmpsasaa = 'sasaa.xvg'
          tmppdb = 'test1.pdb'
          tmpbarepdb = 'test2.pdb'
          tmpbarextc = 'test2.xtc'
          
          # remove all lines but ATOM from the PDB (otherwise gmx sasa freezes)
          
          with open(tmppdb,'w') as f:
             with open(pdb,'r') as g:
                for line in g: 
                   if "ATOM" in line:
                      
                      # Replace modified AA with natural counterparts
                      myres = line[17:20]
                      if myres in AMINO_ACID_VARIANTS_SUBSTITUTION.keys():
                        line.replace(myres,AMINO_ACID_VARIANTS_SUBSTITUTION[myres])
                        
                      f.write(line)
                      
          # Re-open the pdb
          u=mda.Universe(tmppdb,xtc)
          
          sel_P=u.select_atoms('protein')
          sel_P.write(tmpbarepdb)
          
          
          # This is a dirty fix, to be corrected! Issue is that the CRYST1 is added by MDAnalysis and gmx sasa silently freezes on it
          with open(tmpbarepdb) as fff:
            with open('ggg','w') as ggg:
               for line in fff:
                  if "ATOM" in line:
                     ggg.write(line)
                     
          move('ggg',tmpbarepdb)
          # If the protein is static, bare SASA does not matter. But if we just calc relative sasa then we average over trajectory
          with mda.coordinates.XTC.XTCWriter(tmpbarextc,n_atoms=sel_P.atoms.n_atoms) as w:
               for tp in u.trajectory:
                  w.write(sel_P.atoms)
          
          sel_G=u.select_atoms('not protein')
          
          
          with mda.selections.gromacs.SelectionWriter(tmpindex) as w:
              w.write(sel_P,name='prot')
              w.write(sel_G,name='gly')
              w.write(u.atoms,name='system')
          
          sasatimes,sasar,sasaa,_ = get_SASA(tmppdb,xtc,tmpindex,tmpsys,tmpsel,tmpsasa,tmpsasar,tmpsasaa,probe,ndots,endtime)
          
          if keepoutput:
              out1name=pdb.replace(".pdb","")
              copyfile(tmpsasa,"raw_{}_probe_{}.xvg".format(out1name,probe))
              copyfile(tmpsasar,"raw_r_{}_probe_{}.xvg".format(out1name,probe))
              copyfile(tmpsasaa,"raw_a_{}_probe_{}.xvg".format(out1name,probe))
              
          # Now protein only, this will not change so we calc only once
          if iglycan==0:
              baresasatimes,baresasar,baresasaa,_ = get_SASA(tmpbarepdb,tmpbarextc,tmpindex,baresel,baresel,tmpsasa,tmpsasar,tmpsasaa,probe,ndots,endtime)
              
              
              if keepoutput:
                 copyfile(tmpsasa,"bare_probe_{}.xvg".format(probe))
                 copyfile(tmpsasar,"bare_r_probe_{}.xvg".format(probe))
                 copyfile(tmpsasaa,"bare_a_probe_{}.xvg".format(probe))
                 
              
          # Calculate SASA relative to bare protein
          
          # to eliminate warnings:
          baresasar_idx=np.where(baresasar[:,1]!=0)
          
          # write down which atoms do not see any probe. Occupancy=0 means that.
          occupancy_r=np.zeros(len(baresasar[:,1]))
          occupancy_r[baresasar_idx]=1.0
          
          relativesasa=np.zeros(len(baresasar[:,1]))
          relativesasa[baresasar_idx]=(baresasar[:,1][baresasar_idx]-sasar[:baresasar.shape[0],1][baresasar_idx])/baresasar[:,1][baresasar_idx]

          # all atom version
          
          # to eliminate warnings:
          baresasaa_idx=np.where(baresasaa[:,1]!=0)          
          
          # write down which atoms do not see any probe. Occupancy=0 means that.
          occupancy_a=np.zeros(len(baresasaa[:,1]))
          occupancy_a[baresasaa_idx]=1.0

          relativesasaaa=np.zeros(len(baresasaa[:,1]))
          relativesasaaa[baresasaa_idx]=(baresasaa[:,1][baresasaa_idx]-sasaa[:baresasaa.shape[0],1][baresasaa_idx])/baresasaa[:,1][baresasaa_idx]

          # Aappend relative SASA for each residue
          outrelativesasa.append(relativesasa)
          # and for each atom
          outrelativesasaaa.append(relativesasaaa)
          
          
          
          # cleanup
          os.remove(tmpindex)
          os.remove(tmpsasa)
          os.remove(tmpsasar)
          os.remove(tmpsasaa)
          for afile in glob.glob("#sasa*"):
             os.remove(afile)
          
  
         
      outrelativesasa=np.array(outrelativesasa)
      outrelativesasaaa=np.array(outrelativesasaaa)


      # mean or max?
      meanSASA=np.mean(outrelativesasa,axis=0)
      meanSASAaa=np.mean(outrelativesasaaa,axis=0)
      maxSASA=np.max(outrelativesasa,axis=0)
      maxSASAaa=np.max(outrelativesasaaa,axis=0)
      
      
      # adjust output to the mode # CHANGE - OUTPUT BOTH
      if mode=="avg":
         writeSASA=meanSASA
         writeSASAaa=meanSASAaa
      elif mode == "max":
         writeSASA=maxSASA
         writeSASAaa=maxSASAaa
         
         
         
      # make a per-residue plot of the sasa
      if plottrace:
         residues = np.unique(sel_P.resids)
         
         plt.clf()
         i=0
         for sasa in outrelativesasa:
             i+=1
             plt.plot(residues, sasa,lw=0.36,label='glycan {}'.format(i))
             
             plt.xlabel('residue No')
             plt.ylabel('Shielding')
         plt.plot(residues,maxSASA,c='k',ls='--',label='max')
         plt.plot(residues,meanSASA,c='gray',ls='-',label='avg')

         plt.legend()
         plt.ylim(0,1)
         plt.savefig('{}ResidueSASA_probe_{}.pdf'.format(mode,probe))


      # Move back to the reference frame
      # Re-open the pdb
      unw=mda.Universe(tmppdb)
      sel_Pnw=unw.select_atoms('protein')

      if mode == "avg":
         # write to tempfactors
         for iresidue in range(len(sel_P.residues)):
             sel_Pnw.residues[iresidue].atoms.tempfactors=meanSASA[iresidue]*100
             sel_Pnw.residues[iresidue].atoms.occupancies=occupancy_r[iresidue]
             
         sel_Pnw.write('avgResidueSASA_probe_{}.pdb'.format(probe))
         np.savetxt('avgResidueSASA_probe_{}.txt'.format(probe),meanSASA)
      
#      # now write the all-atom values
#      sel_Pnw.atoms.tempfactors=meanSASAaa*100
#      sel_Pnw.atoms.occupancies=occupancy_a
#      sel_Pnw.write('avgAtomSASA_probe_{}.pdb'.format(probe))
#      np.savetxt('avgAtomSASA_probe_{}.txt'.format(probe),meanSASAaa)
      if mode == "max":
         # now maxima:
         for iresidue in range(len(sel_P.residues)):
             sel_Pnw.residues[iresidue].atoms.tempfactors=maxSASA[iresidue]*100
             sel_Pnw.residues[iresidue].atoms.occupancies=occupancy_r[iresidue]
             
         sel_Pnw.write('maxResidueSASA_probe_{}.pdb'.format(probe))
         np.savetxt('maxResidueSASA_probe_{}.txt'.format(probe),maxSASA)
      
#      # now write the all-atom values
#      sel_Pnw.atoms.tempfactors=maxSASAaa*100
#      sel_Pnw.atoms.occupancies=occupancy_a
#      sel_Pnw.write('maxAtomSASA_probe_{}.pdb'.format(probe))
#      np.savetxt('maxAtomSASA_probe_{}.txt'.format(probe),maxSASAaa)

   # cleanup
   os.remove(tmpbarepdb)
   os.remove(tmppdb)
   os.remove(tmpbarextc)
      
if __name__=="__main__":
   parser=ArgumentParser()
   parser.add_argument('--pdblist', dest='pdblist', help='coma-separated list of pdb references [no spaces! test1.pdb,test2.pdb ]')
   parser.add_argument('--xtclist', dest='xtclist', help='coma-separated list of xtc trajs [no spaces! test1.xtc,test2.xtc ]')
   parser.add_argument('--probelist', dest='probelist', help='coma-separated list of probe sizes for SASA calculation, in nm. can be a single one')
   
   parser.add_argument('--plottrace', dest='plottrace', action='store_true', help='plot per-residue shielding? default is no')
   parser.add_argument('--no-plottrace', dest='plottrace', action='store_false', help='do not plot per-residue shielding')
      
   parser.add_argument('--keepoutput', dest='keepoutput', action='store_true',help="do not delete tmp files, default is no")
   parser.add_argument('--no-keepoutput', dest='keepoutput', action='store_false',help="delete tmp files")  
   
   parser.add_argument('--ndots', dest='ndots', help='number of dots used per atom by the Shrake-Rupley algorithm (see Gromacs documentation for details)',default=15)
   parser.add_argument('--mode', dest='mode', help='mode = avg or max, report max or avg Delta_sasa',default="max")
   parser.add_argument('--endframe', dest='endframe', help='last frame to read from the glycan trajectories. Assumes 1 ps time step',default=0)
   
   # GMX required, test if present
   GMXTEST()
   
   args = parser.parse_args()
   pdbfiles=[i for i in args.pdblist.split(",")]
   xtcfiles=[i for i in args.xtclist.split(",")]
   probes=[float(i) for i in args.probelist.split(",")]
   
   plottrace=args.plottrace
   if plottrace=="None":
      plottrace=None
   ndots=args.ndots
   if ndots=="None":
      ndots=None   
   
   mode=args.mode
   
   parser.set_defaults(keepoutput=False)
   parser.set_defaults(plottrace=False)
   
   
   keepoutput=args.keepoutput
   
   
   endframe=args.endframe
   main(pdbfiles,xtcfiles,plottrace,probes,ndots,mode,keepoutput,endframe)
