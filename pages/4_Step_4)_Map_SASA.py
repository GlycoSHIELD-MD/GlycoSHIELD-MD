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

n_procs = st.number_input(label="Processes", min_value=1, max_value=app.get_n_procs(),
    help="Set the number of processes to be used for the SASA computation. Note that more processes may speed up the computation but will at the same time consume more memory."
)

# enable_viz = st.checkbox("Enable 3d visualization", value=False)

glycosasa_progressbar = st.progress(0)
if st.button("Run glycoSASA ..."):
    if n_procs > 1:
        run_parallel = True
    else:
        run_parallel = False
    app.run_glycosasa(glycosasa_progressbar, probelist=probe_values,
                      run_parallel=run_parallel, n_procs=n_procs)
    cfg["have_sasa"] = True
    st.write("You may now proceed to step 5 to download the output data.")

if cfg["have_sasa"]:
    probe_values = cfg["probe_values"]
    probe = st.radio(
        label="Select the probe value to display the plot and visualization for.",
        options=[str(x) for x in probe_values],
        index=0
    )
    st.image(os.path.join(app.get_config()["output_dir"], f"ResidueSASA_probe_{probe}.png"))


# enable_viz = st.checkbox("Enable 3d visualization", value=False)
viz_options = ("off", "360x640", "640x480", "720x520", "1024x768", "1600x1200")
viz_select = st.selectbox("3d visualization", viz_options, index=0)
if viz_select == "off":
    enable_viz = False
    resolution = tuple(viz_options[1].split('x'))
else:
    enable_viz = True
    resolution = tuple(viz_select.split('x'))

if cfg["have_sasa"] and enable_viz:
    # Parameters for vis (can we make it dynamic?)
    width = int(resolution[0])
    height = int(resolution[1])
    app.visualize_sasa(
        os.path.join(app.get_config()["output_dir"], f"maxResidueSASA_probe_{probe}.pdb"),
        probe,
        width,
        height
    )
    st.markdown(f"""<p style="background-color:#ffffff;color:#000000;font-size:24px;border-radius:2%;display:inline;text-align:center">Shielding:&nbsp&nbsp</p>
                <p style="background-color:#BB0103;color:#ffffff;font-size:24px;border-radius:2%;display:inline;text-align:center">&nbsp0&nbsp&nbsp&nbsp%&nbsp</p>
                <p style="background-color:#DB7A7B;color:#ffffff;font-size:24px;border-radius:2%;display:inline;text-align:center">&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp</p>
                <p style="background-color:#FFFEFE;color:#ffffff;font-size:24px;border-radius:2%;display:inline;text-align:center">&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp</p>
                <p style="background-color:#8F88D2;color:#ffffff;font-size:24px;border-radius:2%;display:inline;text-align:center">&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp</p>
                <p style="background-color:#1207A3;color:#ffffff;font-size:24px;border-radius:2%;display:inline;text-align:center">&nbsp100&nbsp%&nbsp</p>
                <p style="background-color:#999999;color:#ffffff;font-size:24px;border-radius:2%;display:inline;text-align:center">&nbspNot Accessible&nbsp</p>""",
                unsafe_allow_html=True)
else:
    st.write("You need to run the SASA calculation before you can see the visualization.")

app.show_sidebar()
