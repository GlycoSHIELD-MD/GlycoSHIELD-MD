import os
import sys
import shutil
import pathlib
import numpy as np
import streamlit as st


# --- functions for configuration management ---
def cfg_init():
    cfg = st.session_state
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
            cfg["output_zip"].rstrip(".zip"),
            "zip",
            cfg["output_dir"]
        )
        with open(cfg["output_zip"], "rb") as f:
            data = f.read()
    else:
        data = ""
    return data

def run_processing():
    cfg = cfg_get()
    if cfg["have_input"]:
        cfg["have_output"] = True


# --- actual web application below ---

st.title('GlycoSHIELD')


uploaded_files = st.file_uploader(
    label="Upload input files for GlycoSHIELD",
    accept_multiple_files=True,
)
store_uploaded_files(uploaded_files)

run_processing()

data = zip_webapp_output()
st.download_button(
    label="Download ZIP file of output data",
    data=data,
    file_name=cfg_get()["output_zip"],
    mime="application/zip"
)
