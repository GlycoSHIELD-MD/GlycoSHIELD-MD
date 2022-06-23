import os
import base64
import getpass
import streamlit as st
import glycoshield.app as app

import streamlit_modal as modal
import streamlit.components.v1 as components
from st_click_detector import click_detector


st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
app.show_header(title="Run glycoSHIELD, glycoTRAJ, and glycoSASA", show_glycoshield_logo=False)
logo_image_style="width:256px;min_width:128px;"


progress_image_obj = st.empty()
app.display_image(app.glycoshield_logo_still, progress_image_obj, image_style=logo_image_style)
glycoshield_progressbar = st.progress(0)
glycostraj_progressbar_1 = st.progress(0)
glycostraj_progressbar_2 = st.progress(0)

if st.button("Run glycoSHIELD and glycoTRAJ ..."):
    app.display_image(app.glycoshield_logo_anim, progress_image_obj, image_style=logo_image_style)
    app.run_glycoshield(glycoshield_progressbar)
    app.run_glycotraj(glycostraj_progressbar_1, glycostraj_progressbar_2)

app.check_glycoshield(glycoshield_progressbar)
app.check_glycotraj(glycostraj_progressbar_1, glycostraj_progressbar_2)
app.display_image(app.glycoshield_logo_still, progress_image_obj, image_style=logo_image_style)


probe_values_str=st.text_input(label="Enter probe values (comma separated)", value="0.14, 0.70")
probe_values = [float(x) for x in probe_values_str.split(',')]

glycosasa_progressbar = st.progress(0)
if st.button("Run glycoSASA ..."):
    app.run_glycosasa(glycosasa_progressbar, probelist=probe_values)
    st.write("Hint: Use the Download button on the sidebar to download the output data as a Zip file.")

app.show_sidebar()
