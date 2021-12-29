# Data tables for glycoshield
# This file intentionally uses long lines.

# Define modified amino acids (approximation: we turn them back to the original residues)
AMINO_ACID_VARIANTS = [
    ("ALA", ["AIB", "ALM", "AYA", "BNN", "CHG", "CSD", "DAL", "DHA", "DNP", "FLA", "HAC", "MAA", "PRR", "TIH", "TPQ", ]),
    ("ARG", ["AGM", "DAR", "HAR", "MMO", "ARM", "ARN", "HMR", "ACL"]),
    ("ASN", ["MEN", "DSG", ]),
    ("ASP", ["DSP", "BHD", "2AS", "ASQ", "ASB", "ASA", "ASK", "ASH", "ASL", "DAS"]),
    ("CYS", ["BCS", "BUC", "C5C", "C6C", "CCS", "CEA", "CME", "CSO", "CSP", "CSS", "CSW", "CSX", "CY1", "CY3", "CYG", "CYM", "CYP", "CYQ", "CYX", "DCY", "EFC", "OCS", "PEC", "PR3", "SCH", "SCS", "SCY", "SHC", "SMC", "SOC"]),
    ("GLU", ["GLH", "GGL", 'PCA', '5HP', 'DGL', 'CGU', 'GMA']),
    ("GLN", [("DGN", "X")]),
    ("GLY", ["GLZ", "SAR", 'NMC', 'GL3', 'GSC', 'MPQ', 'MSA']),
    ("HIS", ["DHI", "HID", "HIC", "HIE", "HIP", "HSD", "HSE", "HSP", "MHS", "NEM", "NEP", "3AH"]),
    ("ILE", ['DIL', 'IIL']),
    ("LEU", ["BUG", "NLE", 'NLP', 'NLN', 'DLE', 'CLE', 'MLE']),
    ("LYS", ['LYM', 'ALY', 'LYZ', 'LYN', 'LLY', 'LLP', 'SHR', 'TRG', 'DLY', 'KCX']),
    ("MET", ["FME", "CXM", "OMT", "MSE"]),
    ("PHE", ["DAH", "DPN", "HPQ", "PHI", "PHL"]),
    ("PRO", ['DPR', 'HYP']),
    ("SER", ['OAS', 'MIS', 'SAC', 'SVA', 'SET', 'SEP', 'SEL', "DSN"]),
    ("THR", ["ALO", "BMT", "DTH", "THO", "TPO"]),
    ("TRP", ["DTR", "HTR", "LTR", "TPL", "TRO"]),
    ("TYR", ["DTY", "IYR", "PAQ", "PTR", "STY", "TYB", "TYM", "TYO", "TYQ", "TYS", "TYY"]),
    ("VAL", ["DIV", "DVA", "MVA"])
]

AMINO_ACID_VARIANTS_SUBSTITUTION = {}

def create_substitution_dict():
    global AMINO_ACID_VARIANTS_SUBSTITUTION
    for resline in AMINO_ACID_VARIANTS:
        res = resline[0]
        for modres in resline[1]:
            AMINO_ACID_VARIANTS_SUBSTITUTION[modres] = res

# create substitution table once at import time
create_substitution_dict()
