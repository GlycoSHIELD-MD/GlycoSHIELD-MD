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
app.show_header(title="Download Output", show_glycoshield_logo=False)

app.zip_webapp_output()
data, size = app.get_webapp_output()
st.download_button(
    label=f"Download Output ({size:.1f} MB)",
    help="Download the output data as a Zip file.",
    data=data,
    file_name=app.get_config()["output_zip"],
    mime="application/zip"
)

app.show_sidebar()
