import streamlit as st
import glycoshield.app as app

st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
# use the global state object to transfer state between the pages
cfg = app.get_config()

app.show_header(title="Run glycoSHIELD, glycoTRAJ, and glycoSASA", show_glycoshield_logo=False)

progress_image_obj = st.empty()
app.display_image(app.glycoshield_logo_still, progress_image_obj, image_style=app.glyco_logo_image_style)
glycoshield_progressbar = st.progress(0)
glycostraj_progressbar_1 = st.progress(0)
glycostraj_progressbar_2 = st.progress(0)

if st.button("Run glycoSHIELD and glycoTRAJ ..."):
    app.display_image(app.glycoshield_logo_anim, progress_image_obj, image_style=app.glyco_logo_image_style)
    app.run_glycoshield(glycoshield_progressbar)
    app.run_glycotraj(glycostraj_progressbar_1, glycostraj_progressbar_2)

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
