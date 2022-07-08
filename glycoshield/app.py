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
from .lib import glycoshield, glycotraj, glycosasa


glycoshield_logo_still = "webapp/glycoshield_still.png"
glycoshield_logo_still = "webapp/GSlogo.png.png"
glycoshield_logo_anim = "webapp/glycoshield_anim.gif"
mpibp_logo = "webapp/mpibp-logo.png"
mpcdf_logo = "webapp/mpcdf-logo.png"
inserm_logo = "webapp/inserm-logo.png"
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
    cfg["input_lines"] = ['#']
    cfg["init"] = True


def get_config():
    if "init" not in st.session_state:
        init_config()
    return st.session_state


def reset_webapp():
    cfg = get_config()
    remove_files = glob.glob(cfg["work_dir"]+"/*") + glob.glob(cfg["output_dir"]+"/*")
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
       zip_fh.write(os.path.join(cfg["output_dir"],cfg["pdbtrajfile"]))

def get_webapp_output_pdbtraj():
    # Here implement a function to download just the multiframe pdb file
    cfg = get_config()
    if webapp_output_ready():
        zipfile = os.path.join(cfg["output_dir"], cfg["pdbtrajfile_zip"])
        with open(zipfile, "rb") as f:
            data = f.read()
        size = os.path.getsize(zipfile) / 1024. / 1024.
    else:
        data = ""
        size = 0
    return data, size
   

def get_webapp_output():
    cfg = get_config()
    if webapp_output_ready():
        zipfile = os.path.join(cfg["work_dir"], cfg["output_zip"])
        with open(zipfile, "rb") as f:
            data = f.read()
        size = os.path.getsize(zipfile) / 1024. / 1024.
    else:
        data = ""
        size = 0
    return data, size


def store_inputs(inputs):
    cfg = get_config()
    with open(os.path.join(cfg["work_dir"], "input_sugaring"), 'w') as f:
        f.write(inputs)


def run_glycoshield(bar,mode="CG",threshold=3.5):
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
    )
    occ = gs.run(streamlit_progressbar=bar)
    occn = []
    # normalise to initial no of frames used for grafting
    for aframe in occ:
        occn.append([])
        for isugar in range(len(aframe)):
            occn[-1].append(aframe[isugar]/float(gs.initialsugarframes[isugar]))
            #~ print(aframe[isugar],float(gs.initialsugarframes[isugar]))
        #~ print(aframe)
    # Capture and write down the occupancy for each site
    df=pd.DataFrame(occn,columns=('{}:{}'.format(gs.chainlist[i],gs.reslist[i]) for i in range(len(occ[0]))))
    
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


def run_glycotraj(bar_1, bar_2,pdbtrajframes = 30):
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
    
    glycotraj(
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
    cfg["glycotraj_done"] = True


def check_glycotraj(bar_1=None, bar_2=None):
    cfg = get_config()
    if cfg["glycotraj_done"]:
        # st.write("Done!")
        if bar_1 is not None:
            bar_1.progress(1.0)
        if bar_2 is not None:
            bar_2.progress(1.0)


def run_glycosasa(streamlit_progressbar=None, probelist=[0.14,0.70]):
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
        run_parallel=True,
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


def visualize(pdb_list):
    from stmol import showmol
    import py3Dmol
    view = py3Dmol.view()
    # view = py3Dmol.view(query="mmtf:1ycr")
    for pdb in pdb_list:
        view.addModel(open(pdb, 'r').read(), 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()
    view.setBackgroundColor('white')
    showmol(view, height=400, width=600)


def visualize_test(pdb):
    from stmol import showmol
    import py3Dmol
    # view = py3Dmol.view()
    # view = py3Dmol.view(query="mmtf:1ycr")
    # for pdb in pdb_list:
        # view.addModel(open(pdb, 'r').read(), 'pdb')
    # view.setStyle(chA, {'cartoon': {'color': 'spectrum'}})

    with open(pdb, 'r') as fp:
        data = fp.read()
    view = py3Dmol.view(
        data=data,
        style={'stick': {'colorscheme': 'greenCarbon'}},
        # style={'cartoon': {'color': 'spectrum'}},
        # query="chain:B"
    )
    chA = {'chain': 'A', 'opacity':0.7, 'color':'white'}

    view.addSurface(py3Dmol.VDW, chA)
    view.zoomTo()
    view.setBackgroundColor('white')
    showmol(view, height=800, width=800)


def visualize_sasa(pdb, height=800, width=1200):
    from stmol import showmol
    import py3Dmol
    with open(pdb, 'r') as fp:
        data = fp.read()
    view = py3Dmol.view(
        data=data,
        style={'stick': {'colorscheme': 'greenCarbon'}},
        width=width,
        height=height
    )
    chA = {'chain': 'A', 'opacity':0.7} #, 'color':'white'}
    view.addSurface(py3Dmol.VDW, chA)
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


@st.cache
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


@st.cache
def get_glycan_clickable_image_html(glycan_lib, glycan_type, n_cols=4):
    # create table with clickable images
    html_figs = []
    for image_label, (d, raw_label, image_file) in glycan_lib[glycan_type].items():
        image_data = load_image(image_file)
        html_figs.append(
            clickable_image_html(image_label, image_data)
        )
    n_elem = len(html_figs)
    n_per_col = n_elem//n_cols + 1
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


def display_image(image_file, streamlit_handle=st, image_style="", href=None):
    with open(image_file, "rb") as fp:
        _filename, extension = os.path.splitext(image_file)
        image_type = extension[1:]
        assert(image_type in ('gif','png'))
        image_data = base64.b64encode(fp.read()).decode("utf-8")
        html = []
        if href is not None:
            html.append(f'<a href="{href}" target="_blank">')
        html.append(f'<img src="data:image/{image_type};base64,{image_data}" style="{image_style}" alt="{image_file}"/>')
        if href is not None:
            html.append(f'</a>')
        streamlit_handle.markdown(
            "".join(html),
            unsafe_allow_html=True
        )


def load_image(image_file):
    with open(image_file, "rb") as fp:
        _filename, extension = os.path.splitext(image_file)
        image_type = extension[1:]
        assert(image_type in ('gif','png'))
        image_data = base64.b64encode(fp.read()).decode("utf-8")
    return image_data


def clickable_image_html(image_label, image_data, image_type="png"):
    return ("<figure>"
           f"<a href='#' id='{image_label}'><img height='96px' src='data:image/{image_type};base64,{image_data}'></a>"
           f"<figcaption>{image_label}</figcaption>"
            "</figure>"
    )

def show_header(title="GlycoSHIELD Web Application", show_institute_logo=True, show_glycoshield_logo=True):
    # logos on top, MPI-BP twice as a placeholder
    if show_institute_logo:
        header_col1, header_col2, header_col3 = st.columns(3)
        logo_image_style="width:192px;min_width:96px;vertical-align:middle;margin:24px 24px"
        display_image(mpibp_logo, streamlit_handle=header_col1, image_style=logo_image_style)
        display_image(inserm_logo, streamlit_handle=header_col2, image_style=logo_image_style)
        display_image(mpcdf_logo, streamlit_handle=header_col3, image_style=logo_image_style)
    st.title(title)
    if show_glycoshield_logo:
        display_image(glycoshield_logo_still, image_style=glyco_logo_image_style)


def show_sidebar():
    with st.sidebar:
        if st.button("Reset Web Application", help="Pushing this button restores the initial state of the application."):
            reset_webapp()
        if getpass.getuser() == "jovyan":
            label = "Shut Down Web Application"
            if st.button(label, help="Pushing this button shuts down the webapp, and you may close the browser tab."):
                quit_binder_webapp()
            st.write("")
            notebook_url = "../lab/tree/TutorialGlycoSHIELD.ipynb"
            display_image(image_file="webapp/glycoshield-tutorial.png", href=notebook_url)
