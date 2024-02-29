# Entry point for the Streamlit-based GlycoSHIELD web application

import streamlit as st
import glycoshield.app as app

st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
# We use a global state object to transfer state between the pages:
cfg = app.get_config()

app.show_header(title="Welcome to the GlycoSHIELD web application!", enable_institute_links=True)

st.markdown(
"""
This application enables users to graft glycan conformer arrays on specific
amino-acids of a selected protein structure.  50 glycan types are currently
available in the library.

To use the application, please navigate from step 1 to step 5 using the sidebar
on the left and follow the instructions provided at each step. On mobile devices,
use the arrow to show the side bar.

This application was developed as an open-access community resource and is
maintained by Klaus Reuter, Matt Sikora and Cyril Hanus.
Please support us by citing our paper
[Tsai et al., Cell, 187 (5), 1296-1311.e26, 2024](https://doi.org/10.1016/j.cell.2024.01.034).

The GlycoSHIELD source code as well as scripts to use GlycoSHIELD are available on
[GitLab](https://gitlab.mpcdf.mpg.de/dioscuri-biophysics/glycoshield-md).
"""
)

app.show_sidebar()
