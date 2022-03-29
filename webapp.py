import os
import base64
# import re
# import sys
import getpass
# import shutil
# import pathlib
# import numpy as np
import streamlit as st
# import MDAnalysis as mda
import glycoshield.app as app


glycoshield_logo_still = "webapp/glycoshield_still.png"
glycoshield_logo_anim = "webapp/glycoshield_anim.gif"


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

    app.display_image_file(glycoshield_logo_still,
        width="10vw", min_width="92px")

    st.title('GlycoSHIELD Web Application')

    with st.sidebar:
        # st.header("Reset Web Application")
        if st.button("Reset Web Application"):
            app.reset_webapp()


    st.header("1. Define input PDB file")
    if not app.get_config()["have_input"]:
        app.use_default_input()
    st.write("By default, the application uses the file {}, alternatively you may upload a custom PDB file.".format(app.get_default_input()))
    uploaded_file = st.file_uploader(
        label="Upload PDB file",
        accept_multiple_files=False)
    if uploaded_file is not None:
        app.store_uploaded_file(uploaded_file)
    # if st.button("Use tutorial input file (EC5.pdb)"):
        # app.use_default_input()
    # app.print_input_pdb()


    st.header("2. Define GlycoSHIELD Input Lines")

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

    # st.text_area('Preview of new input line', new_line)

    col1, col2, col3, col4 = st.columns(4)

    if col1.button("Add"):
        app.add_input_line(new_line)

    if col2.button("Remove"):
        app.rem_input_line(new_line)

    if col3.button("Use default set of input lines"):
        app.clear_input_lines()
        default_input=[
            '#',
            f'A 462,463,464 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {app.get_config()["output_dir"]}/A_463.pdb {app.get_config()["output_dir"]}/A_463.xtc',
            f'A 491,492,493 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {app.get_config()["output_dir"]}/A_492.pdb {app.get_config()["output_dir"]}/A_492.xtc',
        ]
        for line in default_input:
            app.add_input_line(line)

    if col4.button("Clear all input lines"):
        app.clear_input_lines()

    inputs = st.text_area('Current input lines',
                          "\n".join(app.get_input_lines())
                          )

    app.store_inputs(inputs)


    st.header("3. Run glycoSHIELD and glycoTRAJ ...")

    progress_image_obj = st.empty()
    app.display_image_file(glycoshield_logo_still, progress_image_obj)
    glycoshield_progressbar = st.progress(0)
    glycostraj_progressbar_1 = st.progress(0)
    glycostraj_progressbar_2 = st.progress(0)

    if st.button("Run glycoSHIELD and glycoTRAJ ..."):
        app.display_image_file(glycoshield_logo_anim,
            progress_image_obj)
        app.run_glycoshield(glycoshield_progressbar)
        app.run_glycotraj(glycostraj_progressbar_1, glycostraj_progressbar_2)

    app.check_glycoshield(glycoshield_progressbar)
    app.check_glycotraj(glycostraj_progressbar_1, glycostraj_progressbar_2)
    app.display_image_file(glycoshield_logo_still, progress_image_obj)


    st.header("4. Run glycoSASA ...")
    glycosasa_progressbar = st.progress(0)
    if st.button("Run glycoSASA ..."):
        app.run_glycosasa(glycosasa_progressbar)

    if app.check_glycosasa(glycosasa_progressbar):
        st.image(os.path.join(app.get_config()["output_dir"], "ResidueSASA_probe_0.7.png"))
        app.visualize_sasa(
            os.path.join(app.get_config()["output_dir"], "maxResidueSASA_probe_0.7.pdb")
        )
        st.write("Hint: Use the Download button on the sidebar to download the output data as a Zip file.")


    # When running on Binder, offer a shutdown button
    with st.sidebar:
        if getpass.getuser() == "jovyan":
            label = "Quit Web Application"
            # st.header(label)
            # st.write("By pushing \"" + label + "\" the webapp will shut down, and you may close the browser tab.")
            if st.button(label, help="Pushing this button will shut down the webapp, and you may close the browser tab."):
                app.quit_binder_webapp()

    with st.sidebar:
        # st.header("Download Output")
        app.zip_webapp_output()
        data, size = app.get_webapp_output()
        st.download_button(
            label=f"Download Output ({size:.1f} MB)",
            data=data,
            file_name=app.get_config()["output_zip"],
            mime="application/zip"
        )
