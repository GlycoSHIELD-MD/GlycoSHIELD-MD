import streamlit as st
import glycoshield.app as app

st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
# use the global state object to transfer state between the pages
cfg = app.get_config()

app.show_header(title="Run GlycoSHIELD", show_glycoshield_logo=False)

st.markdown(
"""
Graft glycans on residues at step 2. Select mode and distance cutoff below..."""
)


progress_image_obj = st.empty()
app.display_image(app.glycoshield_logo_still, progress_image_obj, image_style=app.glyco_logo_image_style)
glycoshield_progressbar = st.progress(0)
glycostraj_progressbar_1 = st.progress(0)
glycostraj_progressbar_2 = st.progress(0)

# Set the mode and threshold values, no of frames to download (could be put in a single line)
st.write("""Set grafting mode
CG (coarse-grained): cut-off value defines maximal allowed distance between alpha carbons of amino acids and ring oxygens of glycans. This mode requires less calculation and may be used to mimic naturally occurring flexibility of the lateral chains of surface amino-acids. 
AA (all atoms): cut-off value defines maximal allowed distance between all protein and sugar atoms.""")
glycoshield_mode_str=st.selectbox(label="Set the mode for grafting",options=("CG","All"),index=0)
glycoshield_threshold_str=st.text_input(label="Set the grafting cutoff (Angstrom)",value="3.5")
glycoshield_threshold=float(glycoshield_threshold_str)
glycotraj_numpdbframes_str=st.text_input(label="Number of conformers for PDB download",value="30")
glycotraj_numpdbframes=int(glycotraj_numpdbframes_str)

if st.button("Run glycoSHIELD ..."):
    app.display_image(app.glycoshield_logo_anim, progress_image_obj, image_style=app.glyco_logo_image_style)
    app.run_glycoshield(glycoshield_progressbar,mode=glycoshield_mode_str, threshold=glycoshield_threshold)
    app.run_glycotraj(glycostraj_progressbar_1, glycostraj_progressbar_2,pdbtrajframes=glycotraj_numpdbframes)

app.check_glycoshield(glycoshield_progressbar)
app.check_glycotraj(glycostraj_progressbar_1, glycostraj_progressbar_2)
app.display_image(app.glycoshield_logo_still, progress_image_obj, image_style=app.glyco_logo_image_style)


probe_values_str=st.text_input(label="Enter probe values (comma separated)", value="0.14, 0.70")
probe_values = [float(x) for x in probe_values_str.split(',')]
cfg["probe_values"] = probe_values

glycosasa_progressbar = st.progress(0)
if st.button("Run glycoSASA ..."):
    app.run_glycosasa(glycosasa_progressbar, probelist=probe_values)
    cfg["have_sasa"] = True
    st.write("You may now proceed to page 4 and 5 to visualize and download the output data.")

app.show_sidebar()
