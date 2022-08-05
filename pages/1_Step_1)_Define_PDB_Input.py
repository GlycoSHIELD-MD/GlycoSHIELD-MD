import streamlit as st
import streamlit.components.v1 as components
import glycoshield.app as app


st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
app.show_header(title="Upload protein structure", show_glycoshield_logo=False)

cfg = app.get_config()
cfg["clicked_file"] = False


if not cfg["have_input"]:
    app.use_default_input()

st.write(
    "Upload a protein structure in Protein Data Bank (PDB) format using the uploader below. "
    "As a default, the IG-domain of Mouse N-cadherin is used (EC5, PDBid 3Q2W)."
)




uploaded_file = st.file_uploader(
    label="Upload PDB file",
    accept_multiple_files=False)
if uploaded_file is not None:
    app.store_uploaded_file(uploaded_file)
    cfg["clicked_file"] = True
if st.button("Use default protein", help="Use the 5th IG-domain of Mouse N-cadherin (EC5, PDBid 3Q2W)"):
    app.use_default_input()
    cfg["clicked_file"] = True
# app.print_input_pdb()

# NAVIGATION BY BUTTONS MIGHT WORK IN THE FUTURE https://github.com/streamlit/streamlit/issues/4832
#~ st.write("aaa [link](Step_2_Select_Glycans)",unsafe_allow_html=True)

if cfg["clicked_file"] == True:
    # Clean up PDB file...
    app.clean_input_pdb()
    st.write("PDB file ready, please go to Step 2 on the left")

else:
    st.write("Please upload a file or select the default...")

app.show_sidebar()
