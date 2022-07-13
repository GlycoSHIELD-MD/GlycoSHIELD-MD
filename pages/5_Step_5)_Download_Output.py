import streamlit as st
import glycoshield.app as app

st.set_page_config(
    page_title="GlycoSHIELD",
    layout="wide"
)
# use the global state object to transfer state between the pages
cfg = app.get_config()

app.show_header(title="Download Output", show_glycoshield_logo=False)
st.markdown("""
You can download the multimodel file with the selected number of glycan conformers or the zipped folder containing a full GlycoSHIELD run output.

The glycan conformers can be visualised using the standard visualisation tools, e.g.
* [Pymol](https://pymol.org/2/):
  Open the output pdb file na type ```set all_states, on``` in the console to see all glycan conformers.
*  [VMD](https://www.ks.uiuc.edu/Research/vmd/):
   After opening the output file navigate to ```Graphics``` > ```Representations``` and look for a ```trajectory``` tab.
   There, in the "Draw Multiple Frames" field put ```1:1:X``` where X is the total number of conformers you would like to see.
""")


if app.webapp_output_ready():
    # always offer the comparably small PDB trajectory
    app.zip_pdb_trajectory()
    data_pdb, size_pdb = app.get_webapp_output_pdbtraj()
    st.download_button(
        label=f"Download zipped pdb file ({size_pdb:.1f} MB)",
        help="Download the zipped multi-model pdb file containing the protein and multiple glycan conformers.",
        data=data_pdb,
        file_name=app.get_config()["pdbtrajfile_zip"],
        mime="application/zip"
    )
    # optionally offer the potentially very large full output
    if st.button("Create zip file with full GlycoSHIELD output ...",
                 help="Create a zip file with the full output and offer it for download. Zipping may take a while ..."):
        app.zip_webapp_output()
        data, size = app.get_webapp_output()
        st.download_button(
            label=f"Download zip file with full GlycoSHIELD output ({size:.1f} MB)",
            help="Download the complete output data as a zip file.",
            data=data,
            file_name=app.get_config()["output_zip"],
            mime="application/zip"
        )
else:
    st.markdown("**Note**: Please perform steps 1 to 4 first in order to produce results for download.")

app.show_sidebar()
