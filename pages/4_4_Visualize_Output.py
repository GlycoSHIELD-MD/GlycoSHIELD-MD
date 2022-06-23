import os
import streamlit as st
import glycoshield.app as app

st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
# use the global state object to transfer state between the pages
cfg = app.get_config()

app.show_header(title="Analysis and Visualization", show_glycoshield_logo=False)

if cfg["have_sasa"]:
    probe_values = cfg["probe_values"]
    probe = st.radio(
        label="Select the probe value to display the plot and visualization for.",
        options=[str(x) for x in probe_values],
        index=0
    )
    st.image(os.path.join(app.get_config()["output_dir"], f"ResidueSASA_probe_{probe}.png"))
    app.visualize_sasa(
        os.path.join(app.get_config()["output_dir"], f"maxResidueSASA_probe_{probe}.pdb")
    )
else:
    st.write("Please run some computation before using the visualization page.")

app.show_sidebar()
