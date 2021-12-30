__author__ = "Matt Sikora"
__copyright__ = "Copyright 2021, GlycoSHIELD project"
__credits__ = ["Cyril Hanus"]
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Matt Sikora"
__email__ = "masikora@biophys.mpg.de"
__status__ = "Development"


# ~ import matplotlib
#~ matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.analysis import align
from argparse import ArgumentParser
from tqdm import tqdm
from packaging import version
import shlex
import os
import subprocess
import glob
from shutil import copyfile, move
from . import tables


class glycoshield:
    def __init__(self, protpdb, protxtc, inputfile, threshold=3.5, mode="CG", zmin=None, zmax=None, dryrun=False, shuffle_sugar=True, ignorewarn=False, pdbtraj=None, pdbtrajframes=0,verbose=False,path="./"):

        # set vars
        # Print output?
        self.verbose = verbose
        # where to store outputs
        self.path = path

        self.protpdb = protpdb
        self.protxtc = protxtc
        self.inputfile = inputfile
        self.threshold = threshold
        self.mode = mode
        if zmin is None:
            self.zmin = -100000
        else:
            self.zmin = zmin
        if zmax is None:
            self.zmax = 100000
        else:
            self.zmax = zmax

        self.dryrun = dryrun
        self.shuffle_sugar = shuffle_sugar
        self.ignorewarn = ignorewarn

        # parse input on every sugar to be added
        self.inputlines = self.inputparser()

        # Write essentials to a class variable
        self.chainlist = self.inputlines[:, 0]
        self.pdblist = self.inputlines[:, -2]
        self.xtclist = self.inputlines[:, -1]
        self.reslist = [i[1] for i in self.inputlines[:, 1]]

        # when limiting no of atoms - expand max sugar span by this much. Only atoms withing such radius from the sequon are considered for clash checking.
        self.sugar_expand_factor = 1.5

        # pdb output
        self.pdbtraj = pdbtraj
        self.pdbtrajframes = int(pdbtrajframes)

        # Load structure
        # single frame:
        if self.protxtc is None:
            self.uprot = mda.Universe(self.protpdb, trajectory=True)
        else:
            self.uprot = mda.Universe(self.protpdb, self.protxtc)

        if self.mode == "All":
            self.selprot = self.uprot.select_atoms("all")
            self.sugar_sel_suffix = "all"
        elif self.mode == "CG":
            self.selprot = self.uprot.select_atoms("name CA")
            self.sugar_sel_suffix = 'name O5 or ( resname ANE5 and name O6 )'

    def run(self):

        # will hold the number of grafted conformers per frame
        self.graftedconformers = []
        # Iterate over protein states

        protframe = 0
        for prottp in self.uprot.trajectory:
            # counter for "entropy"
            occupancies = []

            # Iterate over sugars:
            for protchain, resids_on_protein, resids_on_sugar, sugarpdb, sugarxtc, pdbout, xtcout in self.inputlines:
                # load sugar
                self.usugar = mda.Universe(sugarpdb, sugarxtc, in_memory=True)
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

                if len(peptide_to_transplant) == 0:
                    self.sugar_without_protein = self.usugar.select_atoms('not protein')
                else:
                    self.sugar_without_protein = self.usugar.select_atoms('not protein or ( protein and resid {} )'.format(" ".join([str(i) for i in peptide_to_transplant])))

                # We should probably be more lenient with the distance check between the first atoms of the peptide and the protein, i.e. dont check residues 2 and 39 in the example above
                if len(resids_on_sugar) != len(self.usugar.select_atoms('protein').residues.resids):
                    buffer = 1  # skip 1 AA after the first fitting AA and one before the last
                    flanked_resids_on_sugar = sorted(list(set(self.usugar.select_atoms('protein').residues.resids) - set(resids_on_sugar)))[buffer:-buffer]
                    self.sugar_without_protein_flanked = self.usugar.select_atoms('not protein or ( protein and resid {} )'.format(" ".join([str(i) for i in flanked_resids_on_sugar])))
                else:
                    self.sugar_without_protein_flanked = self.sugar_without_protein

                if self.mode == "CG":
                    self.sugar_without_protein_flanked = self.sugar_without_protein_flanked.select_atoms(self.sugar_sel_suffix)

                # Initialise a list of sugar frames
                sugarframes = self.get_sugar_frames()

                # Narrow the selection (especially important for large systems)
                self.grid_search_around_sugar(sequon)

                coordinates = []
                # scan all conformations for clashes and append non-clashing ones to the coordinates:
                for n_tp in tqdm(sugarframes, desc="Structure {} residue {}".format(protframe, resids_on_protein[1])):

                    self.process_sugar_conformer(n_tp, tripep, sequon, coordinates)
                coordinates = np.array(coordinates)

                # save counters
                occupancies.append(len(coordinates))

                if len(coordinates) == 0 and self.protxtc is None:
                    # no frames selected. If this is a protein trajectory just ignore this.
                    if not self.ignorewarn:

                        raise BaseException("No frames were accepted for residues {} in chain {}!\nChange threshold or check if the sequon location is accessible!".format(resids_on_protein, protchain))
                    else:
                        print("WARNING!!!\nNo frames were accepted for residues {} in chain {}!\nThis means the sequon location is not accessible given the threshold\nYou decided to ignore this warning so we carry on!".format(resids_on_protein, protchain))
                        continue

                # Write output
                self.write_output(xtcout, pdbout, coordinates, protframe, pdbtraj=self.pdbtraj, pdbtrajframes=self.pdbtrajframes)

                # Cleanup
                del self.usugar


            if self.verbose:
               # Main output printout
               print(protframe, " ".join([str(i) for i in occupancies]))

            # Increment frame counter
            protframe += 1
            self.graftedconformers.append(occupancies)

        return np.array(self.graftedconformers)

    def write_output(self, xtcout, pdbout, coordinates, protframe, pdbtraj=None, pdbtrajframes=0):
        """
        Save if not a dry run. Fixed output format if protein is a multiple frame trajectory

        pdbtraj - option to write trajectory as pdb file.
        pdbtrajframes - only write that many frames to the trajectory
        """

        # If iterating over protein traj - make sure to modify coordinates.
        if self.protxtc is not None:
            xtcout = xtcout.replace('.xtc', '_protein_frame_{}.xtc'.format(protframe))

        if self.protxtc is not None and pdbtraj is not None:
            pdbtraj = pdbtraj.replace('.pdb', '_protein_frame_{}.pdb'.format(protframe))


        # if dry run - skip writing xtc,pdb
        if not self.dryrun:
            # Placeholder for a reconstructed glycoprotein
            glycoprotein = mda.Merge(self.uprot.atoms, self.sugar_without_protein)

            # fill it with the coordinates
            glycoprotein.load_new(coordinates, format=MemoryReader)
            # Write new reference pdb
            glycoprotein.atoms.write(pdbout)

            # Write new trajectory
            with mda.coordinates.XTC.XTCWriter(xtcout, n_atoms=glycoprotein.atoms.n_atoms) as w:
                for tp in glycoprotein.trajectory:
                    w.write(glycoprotein.atoms)

            # pdb trajectory (LARGE). Here the name of the pdb file should be inherited after the xtc output files!
            if pdbtraj is not None:
                write_pdb_trajectory(glycoprotein, pdbtraj, pdbtrajframes)


            del glycoprotein

    def get_sugar_frames(self):
        # iterate over the traj to remove clashing frames. Can be sequential or random
        sugarframes = range(self.usugar.trajectory.n_frames)

        # randomize sugar conformers
        if self.shuffle_sugar:
            np.random.shuffle(list(sugarframes))
        return sugarframes

    def process_sugar_conformer(self, n_tp, tripep, sequon, coordinates):
        """
        Process a single glycan frame. Check for clashes and whether it crosses zmin/zmax. If so, reject. Otherwise - add to
        accepted
        """
        _ = self.usugar.trajectory[n_tp]

        # Align on the protein
        align.alignto(tripep, sequon)

        # check for overlap with a membrane: remove all conformers that have atoms below Zmin
        min = np.min(self.sugar_without_protein_flanked.positions[:, 2])
        max = np.max(self.sugar_without_protein_flanked.positions[:, 2])

        if (not min < self.zmin) and (not max > self.zmax):

            # check for the overlap
            if self.threshold != 0:
                distarr = distance_array(self.reduced_prot_positions, self.sugar_without_protein_flanked.positions)
                if np.min(distarr) < self.threshold:
                    # overlap - do nothing
                    pass
                else:
                    # OK - merge two systems copying protein and adding sugar.
                    coordinates.append(np.vstack((self.uprot.atoms.positions, self.sugar_without_protein.positions)))
            else:
                # Threshold 0 - go ahead
                coordinates.append(np.vstack((self.uprot.atoms.positions, self.sugar_without_protein.positions)))

    def grid_search_around_sugar(self, sequon):
        sugar_root = sequon.center_of_mass()
        # Get size of the sugar:
        sizex = np.max(self.sugar_without_protein_flanked.positions[:, 0]) - np.min(self.sugar_without_protein_flanked.positions[:, 0])
        sizey = np.max(self.sugar_without_protein_flanked.positions[:, 1]) - np.min(self.sugar_without_protein_flanked.positions[:, 1])
        sizez = np.max(self.sugar_without_protein_flanked.positions[:, 2]) - np.min(self.sugar_without_protein_flanked.positions[:, 2])
        maxdiam = np.max([sizex, sizey, sizez]) * self.sugar_expand_factor

        # Generate list of coordinates.
        self.reduced_prot_positions = self.selprot.atoms.positions[(self.selprot.atoms.positions[:, 0] < (sugar_root[0] + maxdiam)) & (self.selprot.atoms.positions[:, 0] > (sugar_root[0] - maxdiam)) &
                                                                   (self.selprot.atoms.positions[:, 1] < (sugar_root[1] + maxdiam)) & (self.selprot.atoms.positions[:, 1] > (sugar_root[1] - maxdiam)) &
                                                                   (self.selprot.atoms.positions[:, 2] < (sugar_root[2] + maxdiam)) & (self.selprot.atoms.positions[:, 2] > (sugar_root[2] - maxdiam))]

    def inputparser(self):
        inputs = []
        with open(self.inputfile) as f:
            for line in f:
                if "#" not in line:
                    l = line.split()
                    protchain = l[0]
                    protsequonresid = [int(i) for i in l[1].split(",")]
                    sugarsequonresid = [int(i) for i in l[2].split(",")]
                    sugarpdb = l[3]
                    sugarxtc = l[4]
                    pdbout = l[5]
                    xtcout = l[6]
                    inputs.append([protchain, protsequonresid, sugarsequonresid, sugarpdb, sugarxtc, pdbout, xtcout])
        return np.array(inputs, dtype=object)


def write_pdb_trajectory(universe, pdbtraj, pdbtrajframes):
    pdbtrajframes = np.min([universe.trajectory.n_frames, pdbtrajframes])
    with mda.Writer(pdbtraj, n_atoms=universe.atoms.n_atoms) as w:
        for tp in universe.trajectory[:pdbtrajframes]:
            w.write(universe.atoms)


def glycotraj(maxframe, outname, pdblist, xtclist, chainlist, reslist, pdbtraj=None, pdbtrajframes=0, path="./"):
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

    resid_per_chain = {}

    for achain, aresid in zip(chainlist, reslist):
        #print(achain, aresid)
        if achain not in resid_per_chain.keys():
            resid_per_chain[achain] = []
        if aresid not in resid_per_chain[achain]:
            resid_per_chain[achain].append(aresid)

    # Iterate over structures and trajs for single glycans to extract subtrajs
    for pdb, xtc, chain, resid in zip(pdblist, xtclist, chainlist, reslist):
        
        newxtc = pdb.replace(".pdb", "_glycan.xtc")
        newpdb = pdb.replace(".pdb", "_glycan.pdb")
      
        u = mda.Universe(pdb, xtc)
        
        # Select glycan (assume this is the only non-protein object, can be made more specific, here we also exclude palmitoylated residues)
        sugarsel = u.select_atoms("not protein and not resname CYSP")

        # save sub-trajs of same length
        with mda.Writer(newxtc, sugarsel.n_atoms) as W:
            for ts in u.trajectory[:maxframe]:
                W.write(sugarsel)

        sugarsel.write(newpdb)
        chains.append(chain)
        resids.append(resid)

    chains = sorted(np.unique(chains))
    # resids=sorted(np.unique(resids))

    sels = []
    coords = []
    
    # open the protein frame and add to the trajectories
    u = mda.Universe(pdblist[0],xtclist[0], in_memory=True)

    # again - this selection can be made more specific
    protein = u.select_atoms('protein')

    # trim
    coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(), protein).run().results
    if version.parse(mda.__version__) >= version.parse("2.0.0"):
        coordinates = coordinates["timeseries"]  # This is only needed in certain MDAnalysis versions

    coordinates = coordinates[:maxframe, :, :]
    sels.append(protein)
    coords.append(coordinates)

    # Now go second time over the files and merge, renumbering the glycans
    for chain in chains:
        # last residue is? we do it per chain
        lastres = protein.select_atoms("segid {}".format(chain)).resids[-1]

        for resid in resid_per_chain[chain]:
            pdb = path + chain + "_" + str(resid) + "_glycan.pdb"

            xtc = pdb.replace(".pdb", ".xtc")
            u1 = mda.Universe(pdb, xtc, in_memory=True)
            sugar = u1.atoms
            sugar.segments.segids = chain
            # update residues
            sugar.residues.resids += lastres

            lastres += len(sugar.residues)

            coordinates = AnalysisFromFunction(lambda ag: ag.positions.copy(), sugar).run().results
            if version.parse(mda.__version__) >= version.parse("2.0.0"):
                coordinates = coordinates["timeseries"]  # This is only needed in certain MDAnalysis versions
            coords.append(coordinates)

            sels.append(u1)

    # Concatenate and save

    coords = np.column_stack(coords)
    u2 = mda.Merge(*[i.atoms for i in sels])
    u2.load_new(coords, format=MemoryReader)
    with mda.Writer(outname + '.xtc', u2.atoms.n_atoms) as W:
        for ts in u2.trajectory:
            W.write(u2.atoms)

    # pdb trajectory (LARGE).
    if pdbtraj is not None:
        write_pdb_trajectory(u2, pdbtraj, pdbtrajframes)

    # Save a reference
    u2.atoms.write(outname + '.pdb')


def get_SASA(pdbname, xtcname, index, groupname, selectgroups, outfile, outfileres, outfileatom, probe, ndots, maxframe):
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
    cmd = "gmx sasa -f {} -s {} -n {} -o {} -or {}".format(xtcname, pdbname, index, outfile, outfileres)
    cmd = cmd + " -oa {}".format(outfileatom)
    cmd = cmd + " -surface \'{}\' -xvg none -probe {} ".format(groupname, probe)
    cmd = cmd + " -output \'{}\' -ndots {}".format(selectgroups, ndots)
    cmd = cmd + " -e {}".format(maxframe)
    #~ print(cmd)

    # Execute
    command = shlex.split(cmd)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = [i for i in proc.stdout]
    err = [i for i in proc.stderr]
    proc.communicate()

    if proc.returncode != 0:
        for line in err:
            print(line)
            raise SystemExit("Error in Gromacs, stopping...")

    # Read in the gmx output files
    sasas = np.loadtxt(outfileres)
    sasaa = np.loadtxt(outfileatom)
    sasatime = np.loadtxt(outfile)
    return sasatime, sasas, sasaa, proc.returncode


def GMXTEST():
    cmd = 'gmx --version --quiet | grep "GROMACS version"'
    command = shlex.split(cmd)
    try:
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except:
        raise BaseException("Gromacs not found or not working, stopping...")
    out = [i for i in proc.stdout]
    err = [i for i in proc.stderr]
    proc.communicate()

    if proc.returncode != 0:
        for line in err:
            print(line)
            raise SystemExit("Gromacs not found, stopping...")


def glycosasa(pdblist, xtclist, plottrace, probelist, ndots, mode, keepoutput, maxframe, path="./"):
    # GMX required, test if present
    GMXTEST()
    # Assumption is there is only a protein in the pdb file.
    # Everything else assumed to be a glycan.

    # list of outputs for each probe size
    outputs = []

    # Iterate over probe sizes
    for probe in probelist:
        outrelativesasa = []
        outrelativesasaaa = []

        # Iterate over trajectories for glycans
        # Load each glycan+protein trajectory
        for iglycan in range(len(xtclist)):

            pdb = pdblist[iglycan]
            xtc = xtclist[iglycan]
            u = mda.Universe(pdb, xtc)

            # Define some temporary file names.
            tmpindex = path + 'index0.ndx'
            tmpsel = 'prot; gly'
            tmpsys = 'system'
            baresel = 'prot'
            tmpsasa = path + 'sasa1.xvg'
            tmpsasar = path + 'sasar.xvg'
            tmpsasaa = path + 'sasaa.xvg'
            tmppdb = path + 'test1.pdb'
            tmpbarepdb = path + 'test2.pdb'
            tmpbarextc = path + 'test2.xtc'

            # remove all lines but ATOM from the PDB (otherwise gmx sasa freezes)

            with open(tmppdb, 'w') as f:
                with open(pdb, 'r') as g:
                    for line in g:
                        if "ATOM" in line:

                            # Replace modified AA with natural counterparts
                            myres = line[17:20]
                            if myres in tables.AMINO_ACID_VARIANTS_SUBSTITUTION.keys():
                                line.replace(myres, tables.AMINO_ACID_VARIANTS_SUBSTITUTION[myres])

                            f.write(line)

            # Re-open the pdb
            u = mda.Universe(tmppdb, xtc)

            sel_P = u.select_atoms('protein')
            sel_P.write(tmpbarepdb)

            # This is a dirty fix, to be corrected! Issue is that the CRYST1 is added by MDAnalysis and gmx sasa silently freezes on it
            with open(tmpbarepdb) as fff:
                with open('ggg', 'w') as ggg:
                    for line in fff:
                        if "ATOM" in line:
                            ggg.write(line)

            move('ggg', tmpbarepdb)
            # If the protein is static, bare SASA does not matter. But if we just calc relative sasa then we average over trajectory
            with mda.coordinates.XTC.XTCWriter(tmpbarextc, n_atoms=sel_P.atoms.n_atoms) as w:
                for tp in u.trajectory:
                    w.write(sel_P.atoms)

            sel_G = u.select_atoms('not protein')

            with mda.selections.gromacs.SelectionWriter(tmpindex) as w:
                w.write(sel_P, name='prot')
                w.write(sel_G, name='gly')
                w.write(u.atoms, name='system')

            sasatimes, sasar, sasaa, _ = get_SASA(tmppdb, xtc, tmpindex, tmpsys, tmpsel, tmpsasa, tmpsasar, tmpsasaa, probe, ndots, maxframe)

            if keepoutput:
                out1name = pdb.replace(".pdb", "")
                copyfile(tmpsasa, path + "raw_{}_probe_{}.xvg".format(out1name, probe))
                copyfile(tmpsasar, path + "raw_r_{}_probe_{}.xvg".format(out1name, probe))
                copyfile(tmpsasaa, path + "raw_a_{}_probe_{}.xvg".format(out1name, probe))

            # Now protein only, this will not change so we calc only once
            if iglycan == 0:
                baresasatimes, baresasar, baresasaa, _ = get_SASA(tmpbarepdb, tmpbarextc, tmpindex, baresel, baresel, tmpsasa, tmpsasar, tmpsasaa, probe, ndots, maxframe)

                if keepoutput:
                    copyfile(tmpsasa, path + "bare_probe_{}.xvg".format(probe))
                    copyfile(tmpsasar, path + "bare_r_probe_{}.xvg".format(probe))
                    copyfile(tmpsasaa, path + "bare_a_probe_{}.xvg".format(probe))

            # Calculate SASA relative to bare protein

            # to eliminate warnings:
            baresasar_idx = np.where(baresasar[:, 1] != 0)

            # write down which atoms do not see any probe. Occupancy=0 means that.
            occupancy_r = np.zeros(len(baresasar[:, 1]))
            occupancy_r[baresasar_idx] = 1.0

            relativesasa = np.zeros(len(baresasar[:, 1]))
            relativesasa[baresasar_idx] = (baresasar[:, 1][baresasar_idx] - sasar[:baresasar.shape[0], 1][baresasar_idx]) / baresasar[:, 1][baresasar_idx]

            # all atom version

            # to eliminate warnings:
            baresasaa_idx = np.where(baresasaa[:, 1] != 0)

            # write down which atoms do not see any probe. Occupancy=0 means that.
            occupancy_a = np.zeros(len(baresasaa[:, 1]))
            occupancy_a[baresasaa_idx] = 1.0

            relativesasaaa = np.zeros(len(baresasaa[:, 1]))
            relativesasaaa[baresasaa_idx] = (baresasaa[:, 1][baresasaa_idx] - sasaa[:baresasaa.shape[0], 1][baresasaa_idx]) / baresasaa[:, 1][baresasaa_idx]

            # Aappend relative SASA for each residue
            outrelativesasa.append(relativesasa)
            # and for each atom
            outrelativesasaaa.append(relativesasaaa)

            # cleanup
            os.remove(tmpindex)
            os.remove(tmpsasa)
            os.remove(tmpsasar)
            os.remove(tmpsasaa)
            for afile in glob.glob(path+"\#sasa*"):
                os.remove(afile)

        outrelativesasa = np.array(outrelativesasa)
        outrelativesasaaa = np.array(outrelativesasaaa)

        # mean or max?
        meanSASA = np.mean(outrelativesasa, axis=0)
        meanSASAaa = np.mean(outrelativesasaaa, axis=0)
        maxSASA = np.max(outrelativesasa, axis=0)
        maxSASAaa = np.max(outrelativesasaaa, axis=0)

        # adjust output to the mode # CHANGE - OUTPUT BOTH
        if mode == "avg":
            writeSASA = meanSASA
            writeSASAaa = meanSASAaa
        elif mode == "max":
            writeSASA = maxSASA
            writeSASAaa = maxSASAaa

        # make a per-residue plot of the sasa
        if plottrace:
            residues = np.unique(sel_P.resids)
            fig = plot_SASA(residues, outrelativesasa, maxSASA, meanSASA, probe, path=path)

        outputs.append([residues, outrelativesasa, maxSASA, meanSASA, probe, occupancy_r])

        # Move back to the reference frame
        # Re-open the pdb
        unw = mda.Universe(tmppdb)
        sel_Pnw = unw.select_atoms('protein')

        if mode == "avg":
            # write to tempfactors
            for iresidue in range(len(sel_P.residues)):
                sel_Pnw.residues[iresidue].atoms.tempfactors = meanSASA[iresidue] * 100
                sel_Pnw.residues[iresidue].atoms.occupancies = occupancy_r[iresidue]

            sel_Pnw.write(path + 'avgResidueSASA_probe_{}.pdb'.format(probe))
            np.savetxt(path + 'avgResidueSASA_probe_{}.txt'.format(probe), meanSASA)


        if mode == "max":
            # now maxima:
            for iresidue in range(len(sel_P.residues)):
                sel_Pnw.residues[iresidue].atoms.tempfactors = maxSASA[iresidue] * 100
                sel_Pnw.residues[iresidue].atoms.occupancies = occupancy_r[iresidue]

            sel_Pnw.write(path + 'maxResidueSASA_probe_{}.pdb'.format(probe))
            np.savetxt(path + 'maxResidueSASA_probe_{}.txt'.format(probe), maxSASA)



    # cleanup
    os.remove(tmpbarepdb)
    os.remove(tmppdb)
    os.remove(tmpbarextc)
    return outputs


def plot_SASA(residues, outrelativesasa, maxSASA, meanSASA, probe, path):

    plt.clf()
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    i = 0
    for sasa in outrelativesasa:
        i += 1
        ax.plot(residues, sasa, lw=0.8, label='glycan {}'.format(i))

        ax.set_xlabel('residue No', size=17)
        ax.set_ylabel('Shielding', size=17)
    ax.plot(residues, maxSASA, c='k', ls='--', label='Max')
    ax.plot(residues, meanSASA, c='gray', ls='-', label='Mean')

    ax.legend()
    ax.set_ylim(0, 1)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=14)

    plt.savefig(path + 'ResidueSASA_probe_{}.pdf'.format(probe))
    plt.show()
