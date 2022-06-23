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
app.show_header()

app.show_sidebar()