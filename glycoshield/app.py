import os
import glob
import base64
import getpass
import zipfile
import pathlib
import pandas as pd
import numpy as np
import streamlit as st
import MDAnalysis as mda
#import pymol2   # import done within function, only if needed
from . import lib
from .lib import glycoshield, glycotraj, glycosasa, clean_pdb


glycoshield_logo_still = "webapp/glycoshield_still.png"
glycoshield_logo_still = "webapp/GSlogo.png"
glycoshield_logo_anim = "webapp/glycoshield_anim.gif"
mpibp_logo = "webapp/mpibp-logo.png"
mpcdf_logo = "webapp/mpcdf-logo.png"
inserm_logo = "webapp/inserm-logo.png"
combined_logo = "webapp/glycoshield-logos-cropped-720.png"
glyco_logo_image_style = "width:256px;min_width:128px;vertical-align:middle;margin:24px 24px"


# --- functions for configuration management ---
def init_config():
    # we use Streamlit's session state to store variables and state between user interaction events
    cfg = st.session_state
    # set up directory and file names
    cfg["tutorial_dir"] = "TUTORIAL"
    cfg["glycan_library_dir"] = "GLYCAN_LIBRARY"
    cfg["work_dir"] = "webapp_work"
    cfg["output_dir"] = "webapp_output"
    cfg["pdb_input"] = ""

    cfg["pdbtrajfile"] = "merged_traj_pdb.pdb"
    cfg["pdbtrajfile_zip"] = cfg["pdbtrajfile"] + ".zip"

    pathlib.Path(cfg["work_dir"]).mkdir(exist_ok=True)
    pathlib.Path(cfg["output_dir"]).mkdir(exist_ok=True)
    cfg["output_zip"] = cfg["output_dir"] + ".zip"
    # flags to implement a finite state machine for the various steps
    cfg["glycoshield_done"] = False
    cfg["glycotraj_done"] = False
    cfg["glycosasa_done"] = False
    cfg["have_input"] = False
    cfg["have_sasa"] = False
    cfg["input_lines"] = ['#']   # raw input lines representing the present input
    # Table for the GUI showing the present input, somewhat redundant with 'input_lines',
    # first line is the labels of the table.
    cfg["input_table"] = _get_empty_input_table()
    cfg["init"] = True


def get_config():
    if "init" not in st.session_state:
        init_config()
    return st.session_state


def reset_webapp():
    cfg = get_config()
    remove_files = glob.glob(cfg["work_dir"] + "/*") + glob.glob(cfg["output_dir"] + "/*")
    for file in remove_files:
        os.unlink(file)
    init_config()


# --- functions defining the steps of the pipeline ---


def store_uploaded_file(uploaded_file):
    cfg = get_config()
    file_name = os.path.join(cfg["work_dir"], uploaded_file.name)
    with open(file_name, "wb") as f:
        f.write(uploaded_file.getbuffer())
    cfg["pdb_input"] = file_name
    cfg["have_input"] = True


def get_default_input():
    cfg = get_config()
    return os.path.join(cfg["tutorial_dir"], "EC5.pdb")


def use_default_input():
    cfg = get_config()
    cfg["pdb_input"] = get_default_input()
    cfg["have_input"] = True


def print_input_pdb():
    cfg = get_config()
    file_name = cfg["pdb_input"]
    st.write("Input PDB file: {}".format(file_name))

def clean_input_pdb():
    """Clean up the PDB to avoid problems later on."""
    cfg = get_config()
    file_name = cfg["pdb_input"]
    tmp_file_name = os.path.join(cfg["work_dir"], "tmp_clean.pdb")
    clean_pdb(file_name, tmp_file_name)
    # replace input with cleaned file
    os.rename(tmp_file_name,file_name)

def webapp_output_ready():
    cfg = get_config()
    if cfg["glycoshield_done"] and cfg["glycotraj_done"] and cfg["glycosasa_done"]:
        return True
    else:
        return False


def zip_webapp_output():
    # zip-directory function inspired from <https://stackoverflow.com/a/1855118>
    def zipdir(path, zip_fh):
        for root, dirs, files in os.walk(path):
            for file in files:
                zip_fh.write(
                    os.path.join(root, file),
                    os.path.relpath(os.path.join(root, file),
                                    os.path.join(path, '..')
                                    )
                )
    cfg = get_config()
    if webapp_output_ready():
        with zipfile.ZipFile(os.path.join(cfg["work_dir"], cfg["output_zip"]), 'w',
                             compression=zipfile.ZIP_DEFLATED, compresslevel=1) as zip_fh:
            zipdir(cfg["output_dir"], zip_fh)


def zip_pdb_trajectory():
    cfg = get_config()
    with zipfile.ZipFile(os.path.join(cfg["output_dir"], cfg["pdbtrajfile_zip"]), 'w',
                         compression=zipfile.ZIP_DEFLATED, compresslevel=1) as zip_fh:
        zip_fh.write(os.path.join(cfg["output_dir"], cfg["pdbtrajfile"]))


def get_webapp_output_pdbtraj():
    # Here implement a function to download just the multiframe pdb file
    cfg = get_config()
    if webapp_output_ready():
        zipfile = os.path.join(cfg["output_dir"], cfg["pdbtrajfile_zip"])
        # with open(zipfile, "rb") as f:
        #     data = f.read()
        # do not read data into buffer, rather pass file handle to streamlit, ultimately rely on garbage collector to close the file
        data = open(zipfile, "rb")
        size = os.path.getsize(zipfile) / 1024. / 1024.
    else:
        data = None
        size = 0
    return data, size


def get_webapp_output():
    cfg = get_config()
    if webapp_output_ready():
        zipfile = os.path.join(cfg["work_dir"], cfg["output_zip"])
        # with open(zipfile, "rb") as f:
        #     data = f.read()
        # do not read data into buffer, rather pass file handle to streamlit, ultimately rely on garbage collector to close the file
        data = open(zipfile, "rb")
        size = os.path.getsize(zipfile) / 1024. / 1024.
    else:
        data = None
        size = 0
    return data, size


def store_inputs(inputs):
    cfg = get_config()
    with open(os.path.join(cfg["work_dir"], "input_sugaring"), 'w') as f:
        f.write(inputs)


def run_glycoshield(bar, mode="CG", threshold=3.5, skip=1):
    "Exposed 'skip' is here to allow control from webapp for large files"
    cfg = get_config()
    pdbtraj = os.path.join(cfg["output_dir"], "test_pdb.pdb")
    pdbtrajframes = 30

    gs = glycoshield(
        protpdb=cfg["pdb_input"],
        protxtc=None,
        inputfile=os.path.join(cfg["work_dir"], "input_sugaring"),
        pdbtraj=pdbtraj,
        pdbtrajframes=pdbtrajframes,
        mode=mode,
        threshold=threshold,
        skip=skip  # DEBUG
    )
    occ = gs.run(streamlit_progressbar=bar)
    occn = []
    # normalise to initial no of frames used for grafting
    for aframe in occ:
        occn.append([])
        for isugar in range(len(aframe)):
            occn[-1].append(aframe[isugar] / float(gs.initialsugarframes[isugar]))
            #~ print(aframe[isugar],float(gs.initialsugarframes[isugar]))
        #~ print(aframe)
    # Capture and write down the occupancy for each site
    df = pd.DataFrame(occn, columns=('{}:{}'.format(gs.chainlist[i], gs.reslist[i]) for i in range(len(occ[0]))))

    st.write('Glycan occupancy')
    st.table(df)
    cfg["gs"] = gs
    cfg["occ"] = occ
    cfg["glycoshield_done"] = True


def check_glycoshield(bar=None):
    cfg = get_config()
    if cfg["glycoshield_done"]:
        # st.write("Done!")
        if bar is not None:
            bar.progress(1.0)
    return cfg["glycoshield_done"]


def run_glycotraj(bar_1, bar_2, pdbtrajframes=30):
    cfg = get_config()
    gs = cfg["gs"]
    occ = cfg["occ"]
    path = cfg["output_dir"]
    maxframe = np.min(occ[0])
    pdblist = gs.pdblist
    xtclist = gs.xtclist
    chainlist = gs.chainlist
    reslist = gs.reslist
    outname = os.path.join(path, "merged_traj")
    pdbtraj = os.path.join(path, cfg["pdbtrajfile"])

    actual_pdbtrajframes = glycotraj(
        maxframe,
        outname,
        pdblist,
        xtclist,
        chainlist,
        reslist,
        pdbtraj,
        pdbtrajframes,
        path,
        streamlit_progressbar_1=bar_1,
        streamlit_progressbar_2=bar_2,
    )
    cfg["glycotraj_actualpdbtrajframes"] = actual_pdbtrajframes
    cfg["glycotraj_done"] = True


def check_glycotraj(bar_1=None, bar_2=None):
    cfg = get_config()
    if cfg["glycotraj_done"]:
        # st.write("Done!")
        if bar_1 is not None:
            bar_1.progress(1.0)
        if bar_2 is not None:
            bar_2.progress(1.0)


def get_n_procs():
    return lib.get_n_procs()


def run_glycosasa(streamlit_progressbar=None, probelist=[0.14, 0.70],
        run_parallel=True, n_procs="auto"):
    cfg = get_config()
    gs = cfg["gs"]
    occ = cfg["occ"]
    path = cfg["output_dir"]
    maxframe = np.min(occ[0])
    maxframe = 10  # temporary
    pdblist = gs.pdblist
    xtclist = gs.xtclist
    # probelist = [0.14, 0.70]
    plottrace = True
    ndots = 15
    mode = "max"
    keepoutput = False
    sasas = glycosasa(
        pdblist=pdblist,
        xtclist=xtclist,
        plottrace=plottrace,
        probelist=probelist,
        ndots=ndots,
        mode=mode,
        keepoutput=keepoutput,
        maxframe=maxframe,
        path=path,
        run_parallel=run_parallel,
        n_procs=n_procs,
        streamlit_progressbar=streamlit_progressbar
    )
    cfg["sasas"] = sasas
    cfg["glycosasa_done"] = True


def check_glycosasa(bar):
    cfg = get_config()
    if cfg["glycosasa_done"]:
        # st.write("Done!")
        if bar is not None:
            bar.progress(1.0)
    return cfg["glycosasa_done"]


def visualize_sasa(pdb, probe, width=1200, height=800):
    from stmol import showmol
    import py3Dmol

    # get output for visualisation
    cfg = get_config()
    sasas_dict = {}
    for blob in cfg["sasas"]:
        probe = blob[4]
        sasas_dict[probe] = blob

    with open(pdb, 'r') as fp:
        data = fp.read()
    view = py3Dmol.view(
        data=data,
        style={'cartoon': {'color': 'gray'}},
        width=width,
        height=height
    )

    probe = float(probe)
    occupancy = np.array(sasas_dict[probe][5])
    residues = np.array(sasas_dict[probe][0])
    # maxSASA = np.array(sasas_dict[probe][2]) * 100. # This reflects the B-factor actual values

    cut = 100. # max for colormap displ. maxSASA is 0-1, so cut at 100 gives the same scale for all proteins
    mid = 50. # Midpoint

    sel_notocc = {'resi': residues[occupancy < 1].tolist()}  # Not accessible
    sel_occ = {'resi': residues[occupancy > 0].tolist()}
    view.addSurface(py3Dmol.SAS, {'opacity': 0.99, 'color': 'gray'}, sel_notocc)
    view.addSurface(py3Dmol.SAS, {'opacity': 0.99, 'colorscheme': {'prop': 'b', 'gradient': 'rwb', 'min': 0, 'max': cut, 'mid': mid}}, sel_occ)

    view.zoomTo()
    view.setBackgroundColor('white')
    showmol(
        view,
        width=width,
        height=height
    )


class visPy3Dmol:
    def __init__(self, path="./tmp_files/"):

        self.pdbfiles = []
        self.xtcfiles = []
        self.n_frames = []
        self.resampledfiles = []
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
        self.proteinfile = self.path + 'tmp_prot.pdb'

        for (pdbfile, xtcfile, framecount) in zip(self.pdbfiles, self.xtcfiles, self.n_frames):
            self.resampledfiles.append([])
            u = mda.Universe(pdbfile, xtcfile)
            # Write only sugar traj and prot reference
            sugar = u.select_atoms('not protein')
            protein = u.select_atoms('protein')

            # ~ with mda.Writer(resampledname, u.atoms.n_atoms) as W:
            iframe = 0
            for ts in u.trajectory[:framecount]:
                resampledname = self.path + "tmp_{}_{}.pdb".format(nfile, iframe)
                self.resampledfiles[nfile].append(resampledname)
                sugar.write(resampledname)
                iframe += 1
            if nfile == 0:
                protein.write(self.proteinfile)

            nfile += 1

    def visualize_brushes(self, height=800, width=1200):
        from .NGL import hex_to_RGB, RGB_to_hex, color_dict, linear_gradient
        from stmol import showmol
        import py3Dmol
        """
        Ideas:
        1. split pdb to protein and sugars (check how it was done with nglview). We could even directly read the A_XXX, B_XXx files?
        2. load sugars independently, xyz file format supports multimodel, maybe the way to go.
        Alternatiuve : addModel for each sugar frame and each sugar?
        https://mobile.twitter.com/david_koes/status/994891637217267712?lang=ca
        https://github.com/3dmol/3Dmol.js/issues/239
        Here they actually do it
        https://chem-workflows.com/articles/2021/10/15/2-virtual-screening/

        3. Steal coloring options from the NGL version.
        """
        components = []
        # Define colors
        if len(self.sugarcolors) != len(self.n_frames):
            # By now we know how many sugar types there are, so we can generate sugar colors,unless user provided colors
            lg = linear_gradient(start_hex=self.startsugarcolor, finish_hex=self.endsugarcolor, n=len(self.n_frames))
            sugarcolor = lg['hex']

        # Get output for visualisation
        cfg = get_config()
        view = py3Dmol.view(width=width, height=height)
        view.removeAllModels()
        #~ view.setViewStyle({'style':'outline','color':'black','width':0.0})
        view.addModel(open(self.proteinfile, 'r').read(), format='pdb')
        Prot = view.getModel()
        Prot.setStyle({'cartoon': {'arrows': True, 'tubes': True, 'style': 'oval', 'color': 'white'}})

        for isugar in range(len(self.n_frames)):
            for iframe in range(self.n_frames[isugar]):
                view.addModel(open(self.resampledfiles[isugar][iframe]).read(), format="pdb")
                zzz = view.getModel()
                zzz.setStyle({}, {'stick': {'color': sugarcolor[isugar], 'radius': 0.1, 'opacity': 1.0}})

        view.zoomTo()
        view.setBackgroundColor('white')
        showmol(
            view,
            width=width,
            height=height
        )


def visualize_brushes(pdb, height=100, width=100):
    # ~ from .NGL import hex_to_RGB, RGB_to_hex, color_dict, linear_gradient
    from stmol import showmol
    import py3Dmol
    """
    Ideas:
    1. split pdb to protein and sugars (check how it was done with nglview). We could even directly read the A_XXX, B_XXx files?
    2. load sugars independently, xyz file format supports multimodel, maybe the way to go.
    Alternatiuve : addModel for each sugar frame and each sugar?
    https://mobile.twitter.com/david_koes/status/994891637217267712?lang=ca
    https://github.com/3dmol/3Dmol.js/issues/239
    Here they actually do it
    https://chem-workflows.com/articles/2021/10/15/2-virtual-screening/

    3. Steal coloring options from the NGL version.
    """
    # Get output for visualisation
    cfg = get_config()

    with open(pdb, 'r') as fp:
        data = fp.read()
    view = py3Dmol.view(
        data=data,
        style={'cartoon': {'color': 'gray'}},
        width=width,
        height=height
    )

    view.zoomTo()
    view.setBackgroundColor('white')
    showmol(
        view,
        width=width,
        height=height
    )


def get_glycan_library_old():
    cfg = get_config()
    # lib = {}
    lib = []
    glycan_listdir = os.listdir(cfg["glycan_library_dir"])
    glycan_listdir.sort()
    for dir_raw in glycan_listdir:
        dir_path = os.path.join(cfg["glycan_library_dir"], dir_raw)
        if os.path.isdir(dir_path):
            # files = os.listdir(dir_path)
            # lib[dir_raw] = list(filter(lambda x: x.endswith(('.xtc', '.pdb')), files))
            lib.append(dir_raw)
    return lib


@st.cache_data
def get_glycan_library():
    lib = {}
    cfg = get_config()
    libdir = cfg["glycan_library_dir"]
    dirs = glob.glob(os.path.join(libdir, 'gs.*.*.*'))
    dirs.sort()
    for d in dirs:
        label = os.path.basename(d)
        gs, num, gtype, name = label.split('.')
        if gtype not in lib:
            lib[gtype] = {}
        thumb = os.path.join(d, "thumbnail.png")
        lib[gtype][name] = (d, label, thumb)
    return lib


@st.cache_data
def get_glycan_clickable_image_html(glycan_lib, glycan_type, n_cols=4):
    # create table with clickable images
    html_figs = []
    for image_label, (d, raw_label, image_file) in glycan_lib[glycan_type].items():
        image_data, _image_type = load_image(image_file)
        html_figs.append(
            clickable_image_html(image_label, image_data)
        )
    n_elem = len(html_figs)
    n_per_col = n_elem // n_cols + 1
    html = []
    html.append('<div class="container">')
    html.append('<div class="row">')
    for i in range(4):
        html.append('<div class="col">')
        for j in range(n_per_col):
            if len(html_figs) > 0:
                html.append(html_figs.pop(0))
        html.append('</div>')
    html.append('</div>')
    html.append('</div>')
    return "\n".join(html)


def get_chain_resids():
    cfg = get_config()
    output = {}
    if cfg["have_input"]:
        u = mda.Universe(cfg["pdb_input"])
        prot = u.select_atoms('protein')
        chains = np.unique(sorted(prot.atoms.segids))
        for chain in chains:
            sel = prot.select_atoms('segid ' + chain)
            output[chain] = np.unique(sel.resids)
    return output


def quit_binder_webapp():
    """Shut down a session running within a Docker container on Binder."""
    os.system("skill -u jovyan")


def create_input_line(chain, resid, glycan):
    cfg = get_config()
    resid_m = int(resid) - 1
    resid_p = int(resid) + 1
    glycan_pdb = os.path.join(
        os.path.join(cfg["glycan_library_dir"], glycan),
        "reference.pdb"
    )
    glycan_xtc = os.path.join(
        os.path.join(cfg["glycan_library_dir"], glycan),
        "conformations.xtc"
    )
    output_pdb = os.path.join(
        cfg["output_dir"],
        f"{chain}_{resid}.pdb"
    )
    output_xtc = os.path.join(
        cfg["output_dir"],
        f"{chain}_{resid}.xtc"
    )
    return f"{chain} {resid_m},{resid},{resid_p} 1,2,3 {glycan_pdb} {glycan_xtc} {output_pdb} {output_xtc}"


def add_input_line(line):
    cfg = get_config()
    if line not in cfg["input_lines"]:
        cfg["input_lines"].append(line)

def rem_input_line(line):
    cfg = get_config()
    try:
        cfg["input_lines"].remove(line)
    except:
        pass

def get_input_lines():
    cfg = get_config()
    return cfg["input_lines"]

def clear_input_lines():
    cfg = get_config()
    cfg["input_lines"] = ['#']


def get_input_table():
    cfg = get_config()
    return cfg["input_table"]

def clear_input_table():
    cfg = get_config()
    cfg["input_table"] = _get_empty_input_table()

def _get_empty_input_table():
    return [('Chain','Residue','Glycan Type','Glycan','Icon'),]

def add_input_row(row):
    cfg = get_config()
    if row not in cfg["input_table"]:
        cfg["input_table"].append(row)

def rem_input_row(row):
    cfg = get_config()
    try:
        cfg["input_table"].remove(row)
    except:
        pass


def get_input_table_html():
    cfg = get_config()
    table = cfg["input_table"]
    html = []
    # html.append('<br>\n')
    html.append('<table>\n')
    for i, row in enumerate(table):
        html.append('<tr>')
        if i == 0:
            td = 'th'
        else:
            td = 'td'
        for j, elem in enumerate(row):
            html.append(f'<{td}>')
            if (i > 0) and (j == 4):
                item = embed_image_into_html(elem)
            else:
                item = str(elem)
            html.append(item)
            html.append(f'</{td}>')
        html.append('</tr>\n')
    html.append('</table>\n')
    # html.append('<br>\n')
    return ''.join(html)


def display_image(image_file, streamlit_handle=st, image_style="", href=None):
    html_string = embed_image_into_html(image_file, image_style, href)
    streamlit_handle.markdown(
        html_string,
        unsafe_allow_html=True
    )

def embed_image_into_html(image_file, image_style="", href=None, image_height="96px"):
    html = []
    image_data, image_type = load_image(image_file)
    if href is not None:
        html.append(f'<a href="{href}" target="_blank">')
    if image_data is not None:
        html.append(f'<img src="data:image/{image_type};base64,{image_data}" height="{image_height}" style="{image_style}" alt="{image_file}"/>')
    if href is not None:
        html.append(f'</a>')
    return "".join(html)


def load_image(image_file):
    image_data = None
    image_type = None
    try:
        with open(image_file, "rb") as fp:
            _filename, extension = os.path.splitext(image_file)
            image_type = extension[1:].lower()
            assert(image_type in ('gif', 'png'))
            image_data = base64.b64encode(fp.read()).decode("utf-8")
    except:
        pass
    finally:
        return image_data, image_type


def clickable_image_html(image_label, image_data, image_type="png", image_height="96px"):
    return ("<figure>"
            f"<a href='#' id='{image_label}'><img height='{image_height}' src='data:image/{image_type};base64,{image_data}'></a>"
            f"<figcaption>{image_label}</figcaption>"
            "</figure>"
            )


def show_header(title="GlycoSHIELD Web Application", show_institute_logo=True,
                show_glycoshield_logo=True, enable_institute_links=False):
    if show_institute_logo:
        # if enable_institute_links:
        #     href = {
        #         "mpibp": "https://www.biophys.mpg.de/",
        #         "inserm": "https://www.inserm.fr/",
        #         "mpcdf": "https://www.mpcdf.mpg.de/"
        #     }
        # else:
        #     href = {
        #         "mpibp": None, "inserm": None, "mpcdf": None
        #     }
        # header_col1, header_col2, header_col3 = st.columns(3)
        # logo_image_style = "height:48px;vertical-align:middle;display:block;"
        # display_image(mpibp_logo, streamlit_handle=header_col1,
        #               image_style=logo_image_style + "margin-left:24px;margin-right:auto;",
        #               href=href["mpibp"])
        # display_image(inserm_logo, streamlit_handle=header_col2,
        #               image_style=logo_image_style + "margin-left:auto;margin-right:auto;",
        #               href=href["inserm"])
        # display_image(mpcdf_logo, streamlit_handle=header_col3,
        #               image_style=logo_image_style + "margin-left:auto;margin-right:24px;",
        #               href=href["mpcdf"])
        header_col = st.columns(1)
        logo_image_style = "vertical-align:middle;display:block;"
        display_image(combined_logo, streamlit_handle=header_col[0],
                      image_style=logo_image_style + "margin-left:auto;margin-right:auto;")
    st.title(title)
    if show_glycoshield_logo:
        display_image(glycoshield_logo_still, image_style=glyco_logo_image_style)


def show_sidebar():
    with st.sidebar:
        if st.button("Reset Web Application", help="Pushing this button restores the initial state of the application."):
            reset_webapp()
        if on_binder():
            label = "Shut Down Web Application"
            if st.button(label, help="Pushing this button shuts down the webapp, and you may close the browser tab."):
                quit_binder_webapp()
            st.write("")
            notebook_url = "../lab/tree/TutorialGlycoSHIELD.ipynb"
            display_image(image_file="webapp/tutorial-button.png", href=notebook_url, image_style="height:32px;")


def on_binder():
    return getpass.getuser() == "jovyan"


def cif_to_pdb(cif_file_name):
    # convert a file protein.cif to protein.pdb, using pymol2
    import pymol2
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(cif_file_name, 'myprotein')
        pdb_file_name = cif_file_name.lower().replace('.cif', '.pdb')
        pymol.cmd.save(pdb_file_name, selection='myprotein')
        return pdb_file_name
