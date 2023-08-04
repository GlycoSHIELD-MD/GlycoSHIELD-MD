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

st.markdown("#### Select Protein Chain and Residue")

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

# print("§§§", chain, resid, glycan_type, d, raw_label, clicked, image_file)
# TODO clean up redundancy that was introduced with the table
#                Chain  Residue  Glycan_Type   Glycan   Icon
new_table_row = (chain, resid,   glycan_type,  clicked, image_file)
new_line = app.create_input_line(chain, resid, raw_label)

st.markdown(f"**{chain}, {resid}, {glycan_type}, {clicked}**",
    help="Current selection of Chain, Residue, Glycan_Type, and Glycan.", unsafe_allow_html=True)

button_col1, button_col2, button_col3, button_col4 = st.columns(4)

if button_col1.button("Add",
        help="Add the current selection to the glycosylation."):
    app.add_input_line(new_line)
    app.add_input_row(new_table_row)
    cfg['have_inputs'] = True

if button_col2.button("Remove",
        help="Remove the current selection from the glycosylation."):
    app.rem_input_line(new_line)
    app.rem_input_row(new_table_row)
    if len(app.get_input_lines()) == 0:
        cfg['have_inputs'] = False

if button_col3.button("Add default glycosylation",
        help="Using this default, we will apply Man5 glycan onto residues 463 and 492 of the N-cadherin domain."):
    app.clear_input_lines()
    app.clear_input_table()
    # Update to reflect the real names of glycans (?)
    default_input = [
        '#',
        f'A 462,463,464 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {app.get_config()["output_dir"]}/A_463.pdb {app.get_config()["output_dir"]}/A_463.xtc',
        f'A 491,492,493 1,2,3 GLYCAN_LIBRARY/Man5.pdb GLYCAN_LIBRARY/Man5_dt1000.xtc {app.get_config()["output_dir"]}/A_492.pdb {app.get_config()["output_dir"]}/A_492.xtc',
    ]
    for line in default_input:
        app.add_input_line(line)
        cfg['have_inputs'] = True

if button_col4.button("Clear all input"):
    app.clear_input_lines()
    app.clear_input_table()
    cfg['have_inputs'] = False



table_html = app.get_input_table_html()
# print(table_html)
st.markdown("##")
st.markdown(table_html, unsafe_allow_html=True, help="Current manually selected inputs")
st.markdown("##")

# handling of string-based inputs for GlycoSHIELD below
inputs = "\n".join(app.get_input_lines())

if st.checkbox("Show advanced input options"):
    st.text_area("Preview of current input line",
                 value=new_line,
                 help="Preview of the input line as it would be generated currently by the above input generator.")
    inputs = st.text_area("All current GlycoSHIELD input lines (press Ctrl+Enter to apply changes after manual editing)",
                          value=inputs,
                          help="This text area shows all current input lines for GlycoSHIELD. You can edit them manually or copy-and-paste your input. Please press Ctrl+Enter to apply changes after manual editing!")

app.store_inputs(inputs)
if cfg['have_inputs']:
    st.write("If happy with the input, please go to Step 3 on the left!")
else:
    st.write("Select at least one residue and glycan type for grafting...")

app.show_sidebar()
