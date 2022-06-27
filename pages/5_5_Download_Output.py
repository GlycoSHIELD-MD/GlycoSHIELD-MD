import streamlit as st
import glycoshield.app as app

st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
app.show_header(title="Download Output", show_glycoshield_logo=False)

# Treat pdb trajectory separately
app.zip_pdb_trajectory()
data_pdb, size_pdb = app.get_webapp_output_pdbtraj()
st.download_button(
    label=f"Download glycosylated protein ({size_pdb:.1f} MB)",
    help="Download the multi-model pdb as a Zip file.",
    data=data_pdb,
    file_name=app.get_config()["pdbtrajfile_zip"],
    mime="application/zip"
)



app.zip_webapp_output()
data, size = app.get_webapp_output()
st.download_button(
    label=f"Download zip file ({size:.1f} MB)",
    help="Download the output data as a Zip file.",
    data=data,
    file_name=app.get_config()["output_zip"],
    mime="application/zip"
)

app.show_sidebar()
