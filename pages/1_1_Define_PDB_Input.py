import streamlit as st
import glycoshield.app as app


st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
app.show_header(title="Define PDB input file", show_glycoshield_logo=False)

if not app.get_config()["have_input"]:
    app.use_default_input()

st.write(
    "You may upload a PDB file using the uploader below. "
    f"In case no file is uploaded, the application uses the file {app.get_default_input()}."
)

uploaded_file = st.file_uploader(
    label="Upload PDB file",
    accept_multiple_files=False)
if uploaded_file is not None:
    app.store_uploaded_file(uploaded_file)
# if st.button("Use tutorial input file (EC5.pdb)"):
    # app.use_default_input()
# app.print_input_pdb()

app.show_sidebar()
