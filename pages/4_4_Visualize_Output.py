import os
import base64
import getpass
import streamlit as st
import glycoshield.app as app

import streamlit_modal as modal
import streamlit.components.v1 as components
from st_click_detector import click_detector


st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
app.show_header(title="Analysis and Visualization", show_glycoshield_logo=False)

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

app.show_sidebar()
