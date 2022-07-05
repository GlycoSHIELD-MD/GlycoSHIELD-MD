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
This application enables users to graft glycan conformer arrays on specific amino-acids of a selected protein structure. 50 glycan types are currently available in the library.


To use the application, please navigate from step 1 to step 5 using the sidebar and follow instructions provided at each step.

This application was developed as an open-access community resource and is maintained by Klaus Reuter, Matt Sikora and Cyril Hanus. Please support us by citing our work: 
[https://doi.org/10.1101/2021.08.04.455134](https://doi.org/10.1101/2021.08.04.455134)

Full GlycoSHIELD source code as well as convenience scripts are available on [gitlab](https://gitlab.mpcdf.mpg.de/MPIBP-Hummer/glycoshield-md/-/tree/refactoring)

"""
)

app.show_sidebar()
