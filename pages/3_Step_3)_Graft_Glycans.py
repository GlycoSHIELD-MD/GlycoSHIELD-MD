import streamlit as st
import glycoshield.app as app
import os


st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
# use the global state object to transfer state between the pages
cfg = app.get_config()

app.show_header(title="Run GlycoSHIELD", show_glycoshield_logo=False)

st.markdown(
    "Graft glycans on the residues selected at Step 2. Select mode and distance cutoff below."
)


progress_image_obj = st.empty()
app.display_image(app.glycoshield_logo_still, progress_image_obj, image_style=app.glyco_logo_image_style)
glycoshield_progressbar = st.progress(0)
glycostraj_progressbar_1 = st.progress(0)
glycostraj_progressbar_2 = st.progress(0)

# Set the mode and threshold values, no of frames to download (could be put in a single line)
st.markdown("""Set parameters for grafting:""")
# Multi-column input
but1, but2, but3, but4 = st.columns(4)
with but1:
    #~ st.header("Set grafting mode")
    glycoshield_mode_str = st.selectbox(label="Set grafting mode", options=("CG", "All"), index=0, help="All (all atoms): graft a conformer if a distance between any of the protein and glycan atoms is larger than the cut-off; CG (coarse-grained): graft a conformer if a distance between protein alpha carbons and glycan ring oxygens is larger than the cut-off. Recommended for large protein structures!")
with but2:
    #~ st.header("Set grafting cut-off (Angstrom)")
    glycoshield_threshold_str = st.text_input(label="Set grafting cutoff (Angstrom)", value="3.5", help="All glycan conformers that come closer to the protein than this value will be removed")
with but3:
    #~ st.header("Set grafting mode")
    glycotraj_numpdbframes_str = st.text_input(label="Number of conformers for PDB download", value="30", help="Warning: large number of conformers can generate huge files!")
with but4:
    glycotraj_skip_str = st.text_input(label="Use only every n-th glycan conformer", value="10", help="Note: Smaller value will generate more complete shields, at a risk of exceeding the application memory for large proteins")

glycoshield_threshold = float(glycoshield_threshold_str)
glycotraj_numpdbframes = int(glycotraj_numpdbframes_str)
glycotraj_skip = int(glycotraj_skip_str)
if st.button("Run glycoSHIELD ..."):
    app.display_image(app.glycoshield_logo_anim, progress_image_obj, image_style=app.glyco_logo_image_style)
    app.run_glycoshield(glycoshield_progressbar, mode=glycoshield_mode_str, threshold=glycoshield_threshold, skip=glycotraj_skip)
    app.run_glycotraj(glycostraj_progressbar_1, glycostraj_progressbar_2, pdbtrajframes=glycotraj_numpdbframes)

app.check_glycoshield(glycoshield_progressbar)
app.check_glycotraj(glycostraj_progressbar_1, glycostraj_progressbar_2)
app.display_image(app.glycoshield_logo_still, progress_image_obj, image_style=app.glyco_logo_image_style)

# enable_viz = st.checkbox("Enable 3d visualization", value=False)
viz_options = ("off", "360x640", "640x480", "720x520", "1024x768", "1600x1200")
viz_select = st.selectbox("3d visualization", viz_options, index=0)
if viz_select == "off":
    enable_viz = False
    resolution = tuple(viz_options[1].split('x'))
else:
    enable_viz = True
    resolution = tuple(viz_select.split('x'))

if enable_viz:
    if cfg["glycotraj_done"]:
        VIS = app.visPy3Dmol(app.get_config()["output_dir"] + "/")
        chainlist = cfg["gs"].chainlist
        reslist = cfg["gs"].reslist
        actual_pdbtrajframes = cfg["glycotraj_actualpdbtrajframes"]

        # Parameters for vis (can we make it dynamic?)
        width = int(resolution[0])
        height = int(resolution[1])

        for (chain, resid) in zip(chainlist, reslist):
            VIS.add_sugar(app.get_config()["output_dir"] + '/{}_{}.pdb'.format(chain, resid),
                        app.get_config()["output_dir"] + '/{}_{}.xtc'.format(chain, resid),
                        actual_pdbtrajframes)
        VIS.subsample()
        VIS.visualize_brushes(height, width)
    else:
        st.write("You need to run GlycoSHIELD before you can see the visualization.")

app.show_sidebar()
