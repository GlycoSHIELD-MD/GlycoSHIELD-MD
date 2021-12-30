import nglview
import MDAnalysis as mda
import numpy as np
import os


# color interpolation from here https://bsouthga.dev/posts/color-gradients-with-python
def hex_to_RGB(hex):
    ''' "#FFFFFF" -> [255,255,255] '''
    # Pass 16 to the integer function for change of base
    return [int(hex[i:i + 2], 16) for i in range(1, 6, 2)]


def RGB_to_hex(RGB):
    ''' [255,255,255] -> "#FFFFFF" '''
    # Components need to be integers for hex to make sense
    RGB = [int(x) for x in RGB]
    return "#" + "".join(["0{0:x}".format(v) if v < 16 else
                          "{0:x}".format(v) for v in RGB])


def color_dict(gradient):
    ''' Takes in a list of RGB sub-lists and returns dictionary of
      colors in RGB and hex form for use in a graphing function
      defined later on '''
    return {"hex": [RGB_to_hex(RGB) for RGB in gradient],
            "r": [RGB[0] for RGB in gradient],
            "g": [RGB[1] for RGB in gradient],
            "b": [RGB[2] for RGB in gradient]}


def linear_gradient(start_hex, finish_hex="#FFFFFF", n=10):
    ''' returns a gradient list of (n) colors between
      two hex colors. start_hex and finish_hex
      should be the full six-digit color string,
      inlcuding the number sign ("#FFFFFF") '''
    # Starting and ending colors in RGB form
    s = hex_to_RGB(start_hex)
    f = hex_to_RGB(finish_hex)
    # Initilize a list of the output colors with the starting color
    RGB_list = [s]
    # Calcuate a color at each evenly spaced value of t from 1 to n
    for t in range(1, n):
        # Interpolate RGB vector for color at the current value of t
        curr_vector = [
            int(s[j] + (float(t) / (n - 1)) * (f[j] - s[j]))
            for j in range(3)
        ]
        # Add it to our list of output colors
        RGB_list.append(curr_vector)

    return color_dict(RGB_list)


class NGL:
    def __init__(self, path="./tmp_files/"):
        self.pdbfiles = []
        self.xtcfiles = []
        self.n_frames = []
        self.sugarcolors = []
        self.path = path
        # Def protein and sugar colors
        self.proteincolor = "#aec1b0"

        self.startsugarcolor = "#009392"
        self.endsugarcolor = "#d0587e"

    def add_sugar(self, pdbfile, xtcfile, n_frames, color=None):
        # Can accept sugar color (hex)
        # assumes user will not change her/his mind and only def color of one sugar...
        self.pdbfiles.append(pdbfile)
        self.xtcfiles.append(xtcfile)
        self.n_frames.append(n_frames)
        if color is not None:
            self.sugarcolors.append(color)

    def subsample(self):
        nfile = 0

        for (pdbfile, xtcfile, framecount) in zip(self.pdbfiles, self.xtcfiles, self.n_frames):
            u = mda.Universe(pdbfile, xtcfile)
            # Write only sugar traj and prot reference
            sugar = u.select_atoms('not protein')
            protein = u.select_atoms('protein')
            with mda.Writer(self.path + "tmp_{}.pdb".format(nfile), u.atoms.n_atoms) as W:
                for ts in u.trajectory[:framecount]:
                    W.write(sugar)
            if nfile == 0:
                protein.write(self.path + 'tmp_prot.pdb')
            nfile += 1

    def build_representation(self):

        components = []
        if len(self.sugarcolors) != len(self.n_frames):
            # By now we know how many sugar types there are, so we can generate sugar colors,unless user provided colors
            lg = linear_gradient(start_hex=self.startsugarcolor, finish_hex=self.endsugarcolor, n=len(self.n_frames))
            sugarcolor = lg['hex']

        v1 = nglview.NGLWidget()
        component = v1.add_component(nglview.FileStructure(os.path.join(self.path, 'tmp_prot.pdb')), default_representation=False)
        components.append(component)

        v1.add_cartoon(selection='protein', color=self.proteincolor, component=0)
        for isugar in range(len(self.n_frames)):

            component = v1.add_component(nglview.FileStructure(os.path.join(self.path, 'tmp_{}.pdb'.format(isugar))), default_representation=False)
            components.append(component)

#             v1.add_component(self.path+'tmp_{}.pdb'.format(isugar),default_representation=False)
            v1.add_representation('licorice', selection='not protein', color=sugarcolor[isugar], component=isugar + 1)
            print(sugarcolor[isugar], isugar + 1, "A")
        v1.center()
        v1.control.zoom(-0.6)

        return v1


class NGLSASA:
    def __init__(self, path, pdbfile, occupancy, residues):
        self.path = path
        self.pdbfile = pdbfile
        self.proteincolor = 'gray'

        # Filter out residues to color gray
        idx_inv = np.where(occupancy < 1)
        idx_vis = np.where(occupancy > 0)
        self.invisible = " ".join([str(i) for i in residues[idx_inv]])
        self.visible = " ".join([str(i) for i in residues[idx_vis]])

    def build_representation(self):
        components = []
        v1 = nglview.NGLWidget()
        component = v1.add_component(nglview.FileStructure(os.path.join(self.path, self.pdbfile)), default_representation=False)
        components.append(component)
        # Colors: https://nglviewer.org/ngl/api/classes/colormakerregistry.html
        # types of the surface: https://nglviewer.org/ngl/api/manual/molecular-representations.html#surface
        # av looks good
        v1.add_representation('surface', selection=self.visible, color='bfactor', component=0, opacity=1, colorScale="OrRd", surfaceType='ms')
        # residues not accessible to a given probe
        v1.add_representation('surface', selection=self.invisible, color='gray', component=0, opacity=1, surfaceType='ms')

        #~ v1.add_cartoon(selection='protein and occupancy>0',color=self.proteincolor,component=0)
        v1.center()
        #~ scheme = nglview.color._ColorScheme([['red', '400-500']], label="aaa")

        return v1
