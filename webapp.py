import streamlit as st
import glycoshield.app as app

st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
# use the global state object to transfer state between the pages
cfg = app.get_config()

app.show_header(title="Welcome to the GlycoSHIELD web application!")

st.markdown(
"""
This application enables users to interactively try the GlycoSHIELD method and implementation.
Users may upload their own PDB data file, alternatively a PDB dataset is provided for test purposes.

To use the GlycoSHIELD web application, navigate step-by-step through the pages accessible via the sidebar, starting from 1 and going to 5.
Each page offers additional help texts that will guide you through the process.

This web application exposes the GlycoSHIELD method and implementation as published in the following reference:

GlycoSHIELD, Matt Sikora, doi:xyz

In case you're using the web application for your research please cite our publication.
"""
)

app.show_sidebar()
