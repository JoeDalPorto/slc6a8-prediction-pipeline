# Canonical Human SLC6A8 (Creatine Transporter 1, UniProt P48029) - Full 635 Amino Acids
REF_SEQ = (
    "MAKKSAENGIYSVSGDEKKGPLIAPGPDGAPAKGDGPVGLGTPGGRLAVPPRETWTRQMD"
    "FIMSCVGFAVGLGNVWRFPYLCYKNGGGVFLIPYVLIALVGGIPIFFLEISLGQFMKAGS"
    "INVWNICPLFKGLGYASMVIVFYCNTYYIMVLAWGFYYLVKSFTTTLPWATCGHTWNTPD"
    "CVEIFRHEDCANASLANLTCDQLADRRSPVIEFWENKVLRLSGGLEVPGALNWEVTLCLL"
    "ACWVLVYFCVWKGVKSTGKIVYFTATFPYVVLVVLLVRGVLLPGALDGIIYYLKPDWSKL"
    "GSPQVWIDAGTQIFFSYAIGLGALTALGSYNRFNNNCYKDAIILALINSGTSFFAGFVVF"
    "SILGFMAAEQGVHISKVAESGPGLAFIAYPRAVTLMPVAPLWAALFFFMLLLLGLDSQFV"
    "GVEGFITGLLDLLPASYYFRFQREISVALCCALCFVIDLSMVTDGGMYVFQLFDYYSASG"
    "TTLLWQAFWECVVVAWVYGADRFMDDIACMIGYRPCPWMKWCWSFFTPLVCMGIFIFNVV"
    "YYEPLVYNNTYVYPWWGEAMGWAFALSSMLCVPLHLLGCLLRAKGTMAERWQHLTQPIWG"
    "LHHLEYRAQDADVRGLTTLTPVSESSKVVVVESVM"

)

# Fetch exact 635-AA sequence from NCBI/UniProt if available, otherwise keep exact 635-AA length
try:
    import requests
    r = requests.get("https://uniprot.org", timeout=8)
    if r.status_code == 200 and len(r.text.splitlines()) > 1:
        seq = "".join(r.text.strip().splitlines()[1:])
        if len(seq) == 635:
            REF_SEQ = seq
except Exception:
    pass

AA_MAP = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic Acid", "C": "Cysteine",
    "E": "Glutamic Acid", "Q": "Glutamine", "G": "Glycine", "H": "Histidine", "I": "Isoleucine",
    "L": "Leucine", "K": "Lysine", "M": "Methionine", "F": "Phenylalanine", "P": "Proline",
    "S": "Serine", "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine"
}

TM_RANGES = {
    "TM1": (61, 83),   "TM2": (95, 117),  "TM3": (153, 175), "TM4": (219, 241),
    "TM5": (253, 275), "TM6": (305, 327), "TM7": (351, 373), "TM8": (405, 427),
    "TM9": (455, 477), "TM10": (491, 513),"TM11": (545, 567),"TM12": (581, 603)
}

HELIX_MAP = {}
for tm_name, (start, end) in TM_RANGES.items():
    for pos in range(start, end + 1):
        HELIX_MAP[pos] = tm_name

LOOP_MAP = {}
for pos in range(1, len(REF_SEQ) + 1):
    if pos not in HELIX_MAP:
        if pos < 61:
            LOOP_MAP[pos] = "N-Terminal Cytoplasmic Tail"
        elif pos > 603:
            LOOP_MAP[pos] = "C-Terminal Cytoplasmic Tail"
        else:
            LOOP_MAP[pos] = "Interhelical Loop"
