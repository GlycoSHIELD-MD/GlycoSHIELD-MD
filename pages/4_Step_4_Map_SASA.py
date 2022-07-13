import streamlit as st
import glycoshield.app as app
import os

st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
# use the global state object to transfer state between the pages
cfg = app.get_config()

app.show_header(title="Run SASA", show_glycoshield_logo=False)

st.markdown(
    "Map the Solvent Accessible Surface Area Reduction on the protein. Select SASA probe radii below."
)

progress_image_obj = st.empty()
app.display_image(app.glycoshield_logo_still, progress_image_obj, image_style=app.glyco_logo_image_style)

probe_values_str = st.text_input(label="Enter probe values (comma separated)", value="0.14, 0.70")
probe_values = [float(x) for x in probe_values_str.split(',')]
cfg["probe_values"] = probe_values

glycosasa_progressbar = st.progress(0)
if st.button("Run glycoSASA ..."):
    app.run_glycosasa(glycosasa_progressbar, probelist=probe_values)
    cfg["have_sasa"] = True
    st.write("You may now proceed to page 4 and 5 to visualize and download the output data.")

if cfg["have_sasa"]:
    probe_values = cfg["probe_values"]
    probe = st.radio(
        label="Select the probe value to display the plot and visualization for.",
        options=[str(x) for x in probe_values],
        index=0
    )
    st.image(os.path.join(app.get_config()["output_dir"], f"ResidueSASA_probe_{probe}.png"))
    app.visualize_sasa2(
        os.path.join(app.get_config()["output_dir"], f"maxResidueSASA_probe_{probe}.pdb"),
        probe
    )
else:
    st.write("You need to run the SASA calculation before you can see the visualization.")

app.show_sidebar()
