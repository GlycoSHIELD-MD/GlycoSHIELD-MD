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

import streamlit_modal as modal
import streamlit.components.v1 as components
from st_click_detector import click_detector


glycoshield_logo_still = "webapp/glycoshield_still.png"
glycoshield_logo_anim = "webapp/glycoshield_anim.gif"
mpibp_logo = "webapp/mpibp-logo.png"
mpcdf_logo = "webapp/mpcdf-logo.png"


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

    with st.sidebar:
        # st.header("Reset Web Application")
        if st.button("Reset Web Application", help="Pushing this button restores the initial state of the application."):
            app.reset_webapp()

    # logos on top, MPI-BP twice as a placeholder
    header_col1, header_col2, header_col3 = st.columns(3)
    logo_image_style="width:256px;min_width:128px;"
    app.display_image(mpibp_logo, streamlit_handle=header_col1, image_style=logo_image_style)
    app.display_image(mpibp_logo, streamlit_handle=header_col2, image_style=logo_image_style)
    app.display_image(mpcdf_logo, streamlit_handle=header_col3, image_style=logo_image_style)

    st.title('GlycoSHIELD Web Application')
    app.display_image(glycoshield_logo_still, image_style=logo_image_style)

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
        # exclude the first and last residues
        resids = chain_resids[chain][1:-1]
    else:
        resids = []
    resid = st.selectbox("Residue", resids)

    st.write("Select Glycan")

    glycan_type = st.selectbox("Glycan Type", glycan_lib.keys())
    #png_tableau = {name:os.path.join(d, "thumbnail.png") for name, d in glycan_lib[glycan_type]}

    # create table with clickable images
    n_cols = 4

    html_figs = []
    for image_label, (d, raw_label, image_file) in glycan_lib[glycan_type].items():
        image_data = app.load_image(image_file)
        html_figs.append(
            app.clickable_image_html(image_label, image_data)
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
    html = "\n".join(html)

    # present the table to the user, allow to select by click
    clicked = click_detector(html)

    # catch the initial non-clicked case
    if clicked == "":
        clicked = list(glycan_lib[glycan_type].keys())[0]

    d, raw_label, image_file = glycan_lib[glycan_type][clicked]
    new_line = app.create_input_line(chain, resid, raw_label)

    st.write('Preview of input line')
    st.text(new_line)

    button_col1, button_col2, button_col3, button_col4 = st.columns(4)

    if button_col1.button("Add"):
        app.add_input_line(new_line)

    if button_col2.button("Remove"):
        app.rem_input_line(new_line)

    if button_col3.button("Use default set of input lines"):
        app.clear_input_lines()
        default_input=[
            '#',
            f'A 462,463,464 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {app.get_config()["output_dir"]}/A_463.pdb {app.get_config()["output_dir"]}/A_463.xtc',
            f'A 491,492,493 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {app.get_config()["output_dir"]}/A_492.pdb {app.get_config()["output_dir"]}/A_492.xtc',
        ]
        for line in default_input:
            app.add_input_line(line)

    if button_col4.button("Clear all input lines"):
        app.clear_input_lines()

    inputs = "\n".join(app.get_input_lines())
    st.write('Current input lines')
    st.text(inputs)

    app.store_inputs(inputs)


    st.header("3. Run glycoSHIELD and glycoTRAJ ...")

    progress_image_obj = st.empty()
    app.display_image(glycoshield_logo_still, progress_image_obj, image_style=logo_image_style)
    glycoshield_progressbar = st.progress(0)
    glycostraj_progressbar_1 = st.progress(0)
    glycostraj_progressbar_2 = st.progress(0)

    if st.button("Run glycoSHIELD and glycoTRAJ ..."):
        app.display_image(glycoshield_logo_anim, progress_image_obj, image_style=logo_image_style)
        app.run_glycoshield(glycoshield_progressbar)
        app.run_glycotraj(glycostraj_progressbar_1, glycostraj_progressbar_2)

    app.check_glycoshield(glycoshield_progressbar)
    app.check_glycotraj(glycostraj_progressbar_1, glycostraj_progressbar_2)
    app.display_image(glycoshield_logo_still, progress_image_obj, image_style=logo_image_style)


    st.header("4. Run glycoSASA ...")

    probe_values_str=st.text_input(label="Enter probe values (comma separated)", value="0.14, 0.70")
    probe_values = [float(x) for x in probe_values_str.split(',')]

    glycosasa_progressbar = st.progress(0)
    if st.button("Run glycoSASA ..."):
        app.run_glycosasa(glycosasa_progressbar, probelist=probe_values)
        st.write("Hint: Use the Download button on the sidebar to download the output data as a Zip file.")


    st.header("5. Analysis and Visualization ...")

    if app.check_glycosasa(glycosasa_progressbar):
        probe = st.radio(
            label="Select the probe value to display the plot and visualization for.",
            options=[str(x) for x in probe_values],
            index=0
        )
        st.image(os.path.join(app.get_config()["output_dir"], f"ResidueSASA_probe_{probe}.png"))
        app.visualize_sasa(
            os.path.join(app.get_config()["output_dir"], f"maxResidueSASA_probe_{probe}.pdb")
        )


    with st.sidebar:
        app.zip_webapp_output()
        data, size = app.get_webapp_output()
        st.download_button(
            label=f"Download Output ({size:.1f} MB)",
            help="Download the output data as a Zip file.",
            data=data,
            file_name=app.get_config()["output_zip"],
            mime="application/zip"
        )

    # When running on Binder, offer a shutdown button
    with st.sidebar:
        if getpass.getuser() == "jovyan":
            label = "Shut Down Web Application"
            if st.button(label, help="Pushing this button shuts down the webapp, and you may close the browser tab."):
                app.quit_binder_webapp()

            st.write("")
            notebook_url = "../lab/tree/TutorialGlycoSHIELD.ipynb"
            app.display_image(image_file="webapp/glycoshield-tutorial.png", href=notebook_url)
