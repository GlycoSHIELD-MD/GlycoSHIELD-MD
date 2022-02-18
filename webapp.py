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
    cfg["tutorial_dir"] = "TUTORIAL"
    cfg["work_dir"] = "webapp_work"
    cfg["output_dir"] = "webapp_output"
    pathlib.Path(cfg["work_dir"]).mkdir(exist_ok=True)
    pathlib.Path(cfg["output_dir"]).mkdir(exist_ok=True)
    cfg["output_zip"] = cfg["output_dir"] + ".zip"
    cfg["have_input"] = False
    cfg["have_output"] = False
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
    if len(uploaded_files) > 0:
        cfg["have_input"] = True
    uploaded_files.clear()

def zip_webapp_output():
    cfg = cfg_get()
    if cfg["have_output"]:
        shutil.make_archive(
            os.path.join(cfg["work_dir"], cfg["output_zip"]).rstrip(".zip"),
            "zip",
            cfg["output_dir"]
        )

def get_webapp_output():
    cfg = cfg_get()
    if cfg["have_output"]:
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

def run_processing(bar):
    cfg = cfg_get()
    pdbtraj = os.path.join(cfg["output_dir"], "test_pdb.pdb")
    pdbtrajframes = 30
    gs = glycoshield(
            protpdb=os.path.join(cfg["tutorial_dir"], "EC5.pdb"),
            protxtc=None,
            inputfile=os.path.join(cfg["work_dir"], "input_sugaring"),
            pdbtraj=pdbtraj,
            pdbtrajframes=pdbtrajframes,
    )
    occ = gs.run(streamlit_progressbar=bar)
    st.write(occ)
    cfg["have_output"] = True


# --- actual web application below ---
st.set_page_config(layout="wide")
st.title('GlycoSHIELD Web App')


st.header("Upload")
uploaded_files = st.file_uploader(
    label="Upload input files for GlycoSHIELD",
    accept_multiple_files=True,
)
store_uploaded_files(uploaded_files)


st.header("Input")
st.button("Add glycan ...")
inputs = st.text_area('Define inputs',
    '#\n'
    f'A 462,463,464 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {cfg_get()["output_dir"]}/A_463.pdb {cfg_get()["output_dir"]}/A_463.xtc\n'
    f'A 491,492,493 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {cfg_get()["output_dir"]}/A_492.pdb {cfg_get()["output_dir"]}/A_492.xtc\n'
)
store_inputs(inputs)


st.header("Run")
bar = st.progress(0)
if st.button("Run glycoSHIELD ..."):
    run_processing(bar)
    zip_webapp_output()


st.header("Download")
data, size = get_webapp_output()
st.download_button(
    label=f"Download ZIP ({size:.1f} MB)",
    data=data,
    file_name=cfg_get()["output_zip"],
    mime="application/zip"
)
