import os
import re
import sys
import getpass
import shutil
import pathlib
import numpy as np
import streamlit as st
import MDAnalysis as mda
from glycoshield.lib import glycoshield, glycotraj, glycosasa
import glycoshield.app as app


if __name__ == "X__main__":
    st.set_page_config(layout="wide")
    st.title('GlycoSHIELD Interactive Web Application')
    pdbs = [
        # 'GLYCAN_LIBRARY/Man5.pdb'
        'webapp_output/maxResidueSASA_probe_0.7.pdb'
    ]
    app.visualize_test(pdb=pdbs[0])



if __name__ == "__main__":
    st.set_page_config(layout="wide")
    st.title('GlycoSHIELD Interactive Web Application')

    st.header("Reset web application")
    if st.button("Reset"):
        app.reset_webapp()

    st.header("Define PDB file for input")
    if not app.get_config()["have_input"]:
        app.use_default_input()
    if st.button("Use default <EC5.pdb>"):
        app.use_default_input()
    uploaded_file = st.file_uploader(
        label="Upload PDB file",
        accept_multiple_files=False,
    )
    if uploaded_file is not None:
        app.store_uploaded_file(uploaded_file)
    app.print_input_pdb()

    st.header("Input")

    chain_resids = app.get_chain_resids()
    # st.write(chain_resids)
    glycan_lib = app.get_glycan_library()
    # st.write(glycan_lib)

    chain = st.selectbox("Chain", chain_resids.keys())

    if chain in chain_resids:
        resids = chain_resids[chain]
    else:
        resids = []
    resid = st.selectbox("Residue", resids)

    glycan = st.selectbox("Glycan", glycan_lib)

    new_line = app.create_input_line(chain, resid, glycan)

    st.text_area('New input line', new_line)

    col1, col2, col3, col4 = st.columns(4)

    if col1.button("Add"):
        app.add_input_line(new_line)

    if col2.button("Remove"):
        app.rem_input_line(new_line)

    inputs = st.text_area('All input lines',
                          "\n".join(app.get_input_lines())
                          # '#\n'
                          # f'A 462,463,464 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {get_config()["output_dir"]}/A_463.pdb {get_config()["output_dir"]}/A_463.xtc\n'
                          # f'A 491,492,493 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {get_config()["output_dir"]}/A_492.pdb {get_config()["output_dir"]}/A_492.xtc\n'
                          )

    if st.button("Clear inputs"):
        app.clear_input_lines()

    app.store_inputs(inputs)

    st.header("Run glycoSHIELD ...")
    bar = st.progress(0)
    if st.button("Run glycoSHIELD ..."):
        app.run_glycoshield(bar)

    # if app.check_glycoshield():
    #     pdb = set()
    #     input_lines = inputs.split('\n')
    #     for line in input_lines:
    #         items = line.split()
    #         if len(items) == 7:
    #             pdb.add(items[3])
    #             pdb.add(items[5])
    #     # pdb = [
    #     #     os.path.join(get_config()["output_dir"], "A_492.pdb"),
    #     #     os.path.join(get_config()["output_dir"], "A_463.pdb"),
    #     # ]
    #     app.visualize(pdb_list=list(pdb))

    st.header("Run glycoTRAJ ...")
    if st.button("Run glycoTRAJ ..."):
        app.run_glycotraj()
    app.check_glycotraj()

    st.header("Run glycoSASA ...")
    if st.button("Run glycoSASA ..."):
        app.run_glycosasa()

    if app.check_glycosasa():
        app.visualize_sasa(
            os.path.join(app.get_config()["output_dir"], "maxResidueSASA_probe_0.7.pdb")
        )

    st.header("Download")
    app.zip_webapp_output()
    data, size = app.get_webapp_output()
    st.download_button(
        label=f"Download ZIP ({size:.1f} MB)",
        data=data,
        file_name=app.get_config()["output_zip"],
        mime="application/zip"
    )

    # When running on Binder, offer a shutdown button
    if getpass.getuser() == "jovyan":
        label = "Quit Web Application"
        st.header(label)
        st.write("By pushing \"" + label + "\" the webapp will shut down, and you may close the browser tab.")
        if st.button(label):
            app.quit_binder_webapp()
