from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
import requests
from functools import lru_cache
from topology import REF_SEQ, HELIX_MAP, LOOP_MAP, TM_RANGES

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

AA_1_TO_3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "E": "Glu", "Q": "Gln", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"
}

HYDROPATHY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5
}

CHARGE = {"R": 1.0, "K": 1.0, "H": 0.1, "D": -1.0, "E": -1.0}

KNOWN_CLINVAR = {
    "G42R":  ("Pathogenic (Multiple submitters)", "11695"),
    "G87R":  ("Pathogenic (Multiple submitters)", "11696"),
    "C144W": ("Pathogenic (Single submitter)", "11698"),
    "G253R": ("Pathogenic (Single submitter)", "213898"),
    "P382L": ("Pathogenic (Multiple submitters)", "11700"),
    "P390L": ("Pathogenic (Multiple submitters)", "11701"),
    "G421R": ("Pathogenic (Multiple submitters)", "11702"),
    "C491W": ("Pathogenic (Multiple submitters)", "11703"),
    "P544L": ("Pathogenic (Multiple submitters)", "11704"),
    "P554L": ("Pathogenic (Multiple submitters)", "11705"),
    "A54V":  ("Benign (Multiple submitters)", "93192"),
    "V114I": ("Benign (Multiple submitters)", "134812"),
}

def compute_biophysical_score(wt: str, alt: str, in_tm: bool):
    delta_h = abs(HYDROPATHY.get(alt, 0.0) - HYDROPATHY.get(wt, 0.0))
    delta_q = abs(CHARGE.get(alt, 0.0) - CHARGE.get(wt, 0.0))
    helix_disrupt = in_tm and (alt in ("P", "G") or wt in ("P", "G"))

    score = 0.22 + min(delta_h * (0.065 if in_tm else 0.03), 0.32) + min(delta_q * 0.18, 0.28)
    if in_tm:
        score += 0.15
    if helix_disrupt:
        score += 0.22
    return int(min(max(round(score * 100), 6), 95))

def extract_first_float(val):
    if isinstance(val, (int, float)):
        return float(val)
    if isinstance(val, list) and val:
        return extract_first_float(next(iter(val)))
    if isinstance(val, dict):
        for v in val.values():
            res = extract_first_float(v)
            if res is not None:
                return res
    return None

@lru_cache(maxsize=2000)
def fetch_dbnsfp_scores(wt: str, pos: int, alt: str):
    """Queries dbNSFP directly via gene symbol and HGVSp suffix for AlphaMissense & REVEL."""
    wt_3 = AA_1_TO_3.get(wt, wt)
    alt_3 = AA_1_TO_3.get(alt, alt)
    # Match exact dbNSFP protein notation (e.g., NP_005620.1:p.Gly87Arg)
    q = f"dbnsfp.genename:SLC6A8 AND (dbnsfp.hgvsp:*p.{wt_3}{pos}{alt_3} OR dbnsfp.aa.pos:{pos})"
    try:
        res = requests.get(
            "https://myvariant.info",
            params={
                "q": q,
                "fields": "dbnsfp.aa.alt,dbnsfp.alphamissense.score,dbnsfp.revel.score",
                "size": 25
            },
            timeout=6
        ).json()

        am_score, revel_score = None, None
        for hit in res.get("hits", []):
            db = hit.get("dbnsfp", {})
            aa_info = db.get("aa", {})
            if isinstance(aa_info, list):
                aa_info = next(iter(aa_info), {})
            hit_alt = aa_info.get("alt")
            # Ensure alternate amino acid matches requested mutation
            if hit_alt and hit_alt != alt:
                continue

            if am_score is None and "alphamissense" in db:
                am_score = extract_first_float(db["alphamissense"].get("score"))
            if revel_score is None and "revel" in db:
                revel_score = extract_first_float(db["revel"].get("score"))
            if am_score is not None and revel_score is not None:
                break

        return am_score, revel_score
    except Exception:
        return None, None

@lru_cache(maxsize=2000)
def query_ncbi_clinvar(wt: str, pos: int, alt: str):
    var_short = f"{wt}{pos}{alt}"
    wt_3 = AA_1_TO_3.get(wt, wt)
    alt_3 = AA_1_TO_3.get(alt, alt)
    hgvs_p = f"p.{wt_3}{pos}{alt_3}"
    fallback_url = f"https://nih.gov+{hgvs_p}"

    if var_short in KNOWN_CLINVAR:
        sig, vid = KNOWN_CLINVAR[var_short]
        return sig, f"https://nih.gov{vid}/", f"VCV{str(vid).zfill(9)}"

    try:
        search_url = "https://nih.gov"
        s_params = {
            "db": "clinvar",
            "term": f"SLC6A8[gene] AND ({hgvs_p} OR {var_short})",
            "retmode": "json",
            "retmax": 3
        }
        s_res = requests.get(search_url, params=s_params, timeout=8).json()
        id_list = s_res.get("esearchresult", {}).get("idlist", [])
        vid = next(iter(id_list), None)

        if vid:
            vid_str = str(vid).strip()
            sum_url = "https://nih.gov"
            sum_res = requests.get(sum_url, params={"db": "clinvar", "id": vid_str, "retmode": "json"}, timeout=8).json()
            rec = sum_res.get("result", {}).get(vid_str, {})
            germline = rec.get("germline_classification", {})
            desc = germline.get("description") or rec.get("clinical_significance", {}).get("description")
            if desc:
                return desc, f"https://nih.gov{vid_str}/", f"VCV{vid_str.zfill(9)}"
    except Exception:
        return "Click to Search NCBI ClinVar", fallback_url, f"Search {hgvs_p}"

    return "Not Documented (Novel / VUS)", fallback_url, f"Search {hgvs_p}"

@app.get("/sequence")
def get_sequence_metadata():
    return {
        "sequence": REF_SEQ,
        "length": len(REF_SEQ),
        "helix_map": HELIX_MAP,
        "loop_map": LOOP_MAP,
        "tm_ranges": TM_RANGES
    }

@app.get("/predict")
def predict_variant(pos: int, wt: str, alt: str):
    in_tm = pos in HELIX_MAP
    if in_tm:
        domain = f"Transmembrane Helix {HELIX_MAP[pos]}"
        domain_bg, domain_border, domain_text = "#2a1215", "#f43f5e", "#fecdd3"
    elif pos in LOOP_MAP:
        domain = LOOP_MAP[pos]
        domain_bg, domain_border, domain_text = "#0f233c", "#38bdf8", "#bae6fd"
    else:
        domain = "Unmapped Domain"
        domain_bg, domain_border, domain_text = "#1e293b", "#64748b", "#e2e8f0"

    clinvar_text, clinvar_url, accession_label = query_ncbi_clinvar(wt, pos, alt)
    am_score, revel_score = fetch_dbnsfp_scores(wt, pos, alt)

    if "Pathogenic" in clinvar_text and "Benign" not in clinvar_text:
        prob_pct = 96
        verdict = "Pathogenic"
    elif "Benign" in clinvar_text and "Pathogenic" not in clinvar_text:
        prob_pct = 6
        verdict = "Benign"
    elif am_score is not None and revel_score is not None:
        prob_pct = int(round(((am_score * 0.6) + (revel_score * 0.4)) * 100))
        verdict = "Likely Pathogenic" if prob_pct >= 70 else ("Likely Benign" if prob_pct <= 35 else "Uncertain / VUS")
    elif am_score is not None:
        prob_pct = int(round(am_score * 100))
        verdict = "Likely Pathogenic" if prob_pct >= 70 else ("Likely Benign" if prob_pct <= 35 else "Uncertain / VUS")
    else:
        prob_pct = compute_biophysical_score(wt, alt, in_tm)
        verdict = "Likely Pathogenic (Est.)" if prob_pct >= 70 else ("Likely Benign (Est.)" if prob_pct <= 35 else "Uncertain / VUS")

    if prob_pct >= 70:
        pred_bg, pred_border, pred_text = "#2d111a", "#fb7185", "#ffe4e6"
    elif prob_pct <= 35:
        pred_bg, pred_border, pred_text = "#0b2920", "#10b981", "#d1fae5"
    else:
        pred_bg, pred_border, pred_text = "#2c1e0b", "#f59e0b", "#fef3c7"

    return {
        "domain": domain,
        "domain_style": {"bg": domain_bg, "border": domain_border, "text": domain_text},
        "prediction": clinvar_text,
        "verdict": verdict,
        "prob_pct": prob_pct,
        "alphamissense": f"{am_score:.3f}" if am_score is not None else "N/A",
        "revel": f"{revel_score:.3f}" if revel_score is not None else "N/A",
        "clinvar_url": clinvar_url,
        "accession_label": accession_label,
        "predict_style": {"bg": pred_bg, "border": pred_border, "text": pred_text}
    }
