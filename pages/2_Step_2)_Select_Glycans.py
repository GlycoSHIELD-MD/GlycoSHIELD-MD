import streamlit as st
import glycoshield.app as app
from st_click_detector import click_detector


st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
cfg = app.get_config()
cfg['have_inputs'] = False

app.show_header(title="Select glycans", show_glycoshield_logo=False)

st.markdown(
    """
To create GlycoSHIELD input, please define the protein chain, residue number and glycan type for each glycosylation site and add it to a PDB file using the "Add" button at the bottom of the page.
N and O-glycans available in the GlycoSHIELD library are listed as "complex" (C), "high-mannose" (M), "hybrid" (H), or "O-glycans" (O).

If you selected the default protein in the previous step, you can use a default glycosylation pattern (Man5) on the residues 463 and 492.
"""
)

st.write("Select Protein Chain and Residue")

chain_resids = app.get_chain_resids()
glycan_lib = app.get_glycan_library()

chain = st.selectbox("Chain", chain_resids.keys())

if chain in chain_resids:
    # exclude the first and last residues
    resids = chain_resids[chain][1:-1]
else:
    resids = []
resid = st.selectbox("Residue", resids)

st.write("Select Glycan")

glycan_type = st.selectbox("Glycan Type", glycan_lib.keys())

# present a table of clickable preview images to the user
html = app.get_glycan_clickable_image_html(glycan_lib, glycan_type)
clicked = click_detector(html)

# catch the initial non-clicked case
if clicked == "":
    clicked = list(glycan_lib[glycan_type].keys())[0]

d, raw_label, image_file = glycan_lib[glycan_type][clicked]
new_line = app.create_input_line(chain, resid, raw_label)

st.write('Preview of input line')
st.text(new_line)

button_col1, button_col2, button_col3, button_col4 = st.columns(4)

if button_col1.button("Add"):
    app.add_input_line(new_line)
    cfg['have_inputs'] = True

if button_col2.button("Remove"):
    app.rem_input_line(new_line)
    if len(app.get_input_lines()) == 0:
        cfg['have_inputs'] = False
if button_col3.button("Add default glycosylation", help="By default, we will apply Man5 glycan onto residues 463 and 492 of the N-cadherin domain"):
    app.clear_input_lines()
    # Update to reflect the real names of glycans (?)
    default_input = [
        '#',
        f'A 462,463,464 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {app.get_config()["output_dir"]}/A_463.pdb {app.get_config()["output_dir"]}/A_463.xtc',
        f'A 491,492,493 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {app.get_config()["output_dir"]}/A_492.pdb {app.get_config()["output_dir"]}/A_492.xtc',
    ]
    for line in default_input:
        app.add_input_line(line)
        cfg['have_inputs'] = True

if button_col4.button("Clear all input lines"):
    app.clear_input_lines()
    cfg['have_inputs'] = False

inputs = "\n".join(app.get_input_lines())
st.write('Current input lines')
st.text(inputs)

app.store_inputs(inputs)
if cfg['have_inputs']:

    st.write("If happy with the input, please go to Step 3 on the left")
else:
    st.write("Select at least one residue and glycan type for grafting...")

app.show_sidebar()
