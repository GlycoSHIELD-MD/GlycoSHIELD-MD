import os
import sys
import shutil
import pathlib
import numpy as np
import streamlit as st
from glycoshield.lib import glycoshield, glycotraj, glycosasa


# --- functions for configuration management ---
def cfg_init():
    cfg = st.session_state
    # set up directory and file names
    cfg["tutorial_dir"] = "TUTORIAL"
    cfg["work_dir"] = "webapp_work"
    cfg["output_dir"] = "webapp_output"
    cfg["pdb_input"] = ""
    pathlib.Path(cfg["work_dir"]).mkdir(exist_ok=True)
    pathlib.Path(cfg["output_dir"]).mkdir(exist_ok=True)
    cfg["output_zip"] = cfg["output_dir"] + ".zip"
    # flags to implement a finite state machine for the various steps
    cfg["glycoshield_done"] = False
    cfg["glycotraj_done"] = False
    cfg["glycosasa_done"] = False
    cfg["have_input"] = False
    # 
    cfg["init"] = True

def cfg_get():
    if "init" not in st.session_state:
        cfg_init()
    return st.session_state


# --- functions defining the steps of the pipeline ---
def store_uploaded_files(uploaded_files):
    cfg = cfg_get()
    for file in uploaded_files:
        #bytes_data = file.read()
        # st.write("filename:", file.name)
        file_name = os.path.join(cfg["work_dir"], file.name)
        with open(file_name, "wb") as f:
            f.write(file.getbuffer())
        cfg["pdb_input"].append(file_name)
    if len(uploaded_files) > 0:
        cfg["have_input"] = True
    uploaded_files.clear()

def use_default_input():
    cfg = cfg_get()
    default_pdb = os.path.join(cfg["tutorial_dir"], "EC5.pdb")
    cfg["pdb_input"] = default_pdb
    cfg["have_input"] = True
    st.write("Using {}".format(default_pdb))

def webapp_output_ready():
    cfg = cfg_get()
    if cfg["glycoshield_done"] and cfg["glycotraj_done"] and cfg["glycosasa_done"]:
        return True
    else:
        return False

def zip_webapp_output():
    cfg = cfg_get()
    if webapp_output_ready():
        shutil.make_archive(
            os.path.join(cfg["work_dir"], cfg["output_zip"]).rstrip(".zip"),
            "zip",
            cfg["output_dir"]
        )

def get_webapp_output():
    cfg = cfg_get()
    if webapp_output_ready():
        zipfile = os.path.join(cfg["work_dir"], cfg["output_zip"])
        with open(zipfile, "rb") as f:
            data = f.read()
        size = os.path.getsize(zipfile)/1024./1024.
    else:
        data = ""
        size = 0
    return data, size

def store_inputs(inputs):
    cfg = cfg_get()
    with open(os.path.join(cfg["work_dir"], "input_sugaring"),'w') as f:
        f.write(inputs)

def run_glycoshield(bar):
    cfg = cfg_get()
    pdbtraj = os.path.join(cfg["output_dir"], "test_pdb.pdb")
    pdbtrajframes = 30
    gs = glycoshield(
            protpdb=cfg["pdb_input"],  # os.path.join(cfg["tutorial_dir"], "EC5.pdb"),
            protxtc=None,
            inputfile=os.path.join(cfg["work_dir"], "input_sugaring"),
            pdbtraj=pdbtraj,
            pdbtrajframes=pdbtrajframes,
    )
    occ = gs.run(streamlit_progressbar=bar)
    st.write(occ)
    cfg["gs"] = gs
    cfg["occ"] = occ
    cfg["glycoshield_done"] = True


def check_glycoshield():
    cfg = cfg_get()
    if cfg["glycoshield_done"]:
        st.write("Done!")
    return cfg["glycoshield_done"]


def run_glycotraj():
    cfg = cfg_get()
    gs = cfg["gs"]
    occ = cfg["occ"]
    path = cfg["output_dir"]
    maxframe=np.min(occ[0])
    pdblist=gs.pdblist
    xtclist=gs.xtclist
    chainlist=gs.chainlist
    reslist=gs.reslist
    outname=os.path.join(path, "merged_traj")
    pdbtraj=os.path.join(path, "test_merged_pdb.pdb")
    pdbtrajframes=30
    glycotraj(
        maxframe,
        outname,
        pdblist,
        xtclist,
        chainlist,
        reslist,
        pdbtraj,
        pdbtrajframes,
        path
    )
    cfg["glycotraj_done"] = True


def check_glycotraj():
    cfg = cfg_get()
    if cfg["glycotraj_done"]:
        st.write("Done!")


def run_glycosasa():
    cfg = cfg_get()
    gs = cfg["gs"]
    occ = cfg["occ"]
    path = cfg["output_dir"]
    maxframe=np.min(occ[0])
    maxframe=10 # temporary
    pdblist=gs.pdblist
    xtclist=gs.xtclist
    probelist=[0.14,0.70] # Possibly user will have an option to choose one value, makes it faster and easier to manage visualisation
    plottrace=True
    ndots=15
    mode="max"
    keepoutput=False
    sasas = glycosasa(
        pdblist=pdblist,
        xtclist=xtclist,
        plottrace=plottrace,
        probelist=probelist,
        ndots=ndots,
        mode=mode,
        keepoutput=keepoutput,
        maxframe=maxframe,
        path=path
    )
    cfg["sasas"] = sasas
    cfg["glycosasa_done"] = True


def check_glycosasa():
    cfg = cfg_get()
    if cfg["glycosasa_done"]:
        st.write("Done!")


def visualize(pdb_list):
    from stmol import showmol
    import py3Dmol
    view = py3Dmol.view()
    for pdb in pdb_list:
        view.addModel(open(pdb, 'r').read(),'pdb')
    view.setStyle({'cartoon':{'color':'spectrum'}})
    view.zoomTo()
    view.setBackgroundColor('white')
    showmol(view, height = 400, width=600)


# --- actual web application below ---
st.set_page_config(layout="wide")
st.title('GlycoSHIELD Web Application')


st.header("Upload")
uploaded_files = st.file_uploader(
    label="Upload input files for GlycoSHIELD",
    accept_multiple_files=True,
)
store_uploaded_files(uploaded_files)
if st.button("Use default PDB"):
    use_default_input()


# TODO make input config interactive and easy to access!!
st.header("Input")
#st.button("Add glycan ...")
inputs = st.text_area('Define inputs',
    '#\n'
    f'A 462,463,464 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {cfg_get()["output_dir"]}/A_463.pdb {cfg_get()["output_dir"]}/A_463.xtc\n'
    f'A 491,492,493 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {cfg_get()["output_dir"]}/A_492.pdb {cfg_get()["output_dir"]}/A_492.xtc\n'
)
store_inputs(inputs)


st.header("Run glycoSHIELD ...")
bar = st.progress(0)
if st.button("Run glycoSHIELD ..."):
    run_glycoshield(bar)

if check_glycoshield():
    pdb = [
        os.path.join(cfg_get()["output_dir"], "A_492.pdb"),
        os.path.join(cfg_get()["output_dir"], "A_463.pdb"),
    ]
    visualize(pdb_list=pdb)


st.header("Run glycoTRAJ ...")
if st.button("Run glycoTRAJ ..."):
    run_glycotraj()
check_glycotraj()


st.header("Run glycoSASA ...")
if st.button("Run glycoSASA ..."):
    run_glycosasa()
check_glycosasa()


st.header("Download")
zip_webapp_output()
data, size = get_webapp_output()
st.download_button(
   label=f"Download ZIP ({size:.1f} MB)",
   data=data,
   file_name=cfg_get()["output_zip"],
   mime="application/zip"
)
