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

# Kyte-Doolittle Hydropathy Index
HYDROPATHY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5
}

# Formal Physiological Charge (pH 7.4)
CHARGE = {
    "R": 1.0, "K": 1.0, "H": 0.1,
    "D": -1.0, "E": -1.0
}

# Side-chain Molecular Volume (Angstroms^3)
VOLUME = {
    "G": 60, "A": 88, "S": 89, "C": 108, "D": 111, "P": 112, "N": 114,
    "T": 116, "E": 138, "V": 140, "Q": 143, "H": 153, "M": 162, "I": 166,
    "L": 166, "K": 168, "R": 173, "F": 189, "Y": 193, "W": 227
}

def compute_physicochemical_metrics(pos: int, wt: str, alt: str, in_tm: bool):
    """Calculates biophysical disruption metrics and pathogenicity probability %."""
    h_wt, h_alt = HYDROPATHY.get(wt, 0.0), HYDROPATHY.get(alt, 0.0)
    delta_h = round(h_alt - h_wt, 2)

    q_wt, q_alt = CHARGE.get(wt, 0.0), CHARGE.get(alt, 0.0)
    delta_q = round(q_alt - q_wt, 1)

    v_wt, v_alt = VOLUME.get(wt, 130), VOLUME.get(alt, 130)
    delta_v = v_alt - v_wt

    # Helix disruption flag (Proline/Glycine changes inside TM helices)
    helix_disrupt = in_tm and (alt in ("P", "G") or wt in ("P", "G"))

    # Weighted pathogenicity score (0.0 to 0.99)
    raw_score = 0.18
    raw_score += min(abs(delta_h) * (0.065 if in_tm else 0.03), 0.35)
    raw_score += min(abs(delta_q) * (0.18 if in_tm else 0.10), 0.30)
    raw_score += min(abs(delta_v) / 450.0, 0.15)
    if in_tm:
        raw_score += 0.14
    if helix_disrupt:
        raw_score += 0.22

    prob_pct = int(min(max(round(raw_score * 100), 4), 98))

    # Human-readable biophysical summary tags
    flags = []
    if helix_disrupt:
        flags.append("TM Helix Breaker")
    if in_tm and q_alt != 0.0 and q_wt == 0.0:
        flags.append("Membrane Charge Burial")
    elif abs(delta_q) >= 1.0:
        flags.append("Electrostatic Shift")
    if abs(delta_h) >= 3.5:
        flags.append("Severe Hydropathy Shift")
    if abs(delta_v) >= 65:
        flags.append("Steric Volume Clash")
    if not flags:
        flags.append("Conservative Substitution")

    return {
        "delta_h": f"{delta_h:+.2f}",
        "delta_q": f"{delta_q:+.1f}",
        "delta_v": f"{delta_v:+d} Å³",
        "biophysical_summary": " • ".join(flags),
        "prob_pct": prob_pct
    }

@app.get("/sequence")
def get_sequence_metadata():
    return {
        "sequence": REF_SEQ,
        "length": len(REF_SEQ),
        "helix_map": HELIX_MAP,
        "loop_map": LOOP_MAP,
        "tm_ranges": TM_RANGES
    }

@lru_cache(maxsize=2000)
def fetch_live_clinvar(wt: str, pos: int, alt: str):
    wt_3 = AA_1_TO_3.get(wt, wt)
    alt_3 = AA_1_TO_3.get(alt, alt)
    query = f"clinvar.gene.symbol:SLC6A8 AND clinvar.hgvs.protein:p.{wt_3}{pos}{alt_3}"
    url = f"https://myvariant.info{query}&fields=clinvar"
    
    try:
        resp = requests.get(url, timeout=5)
        data = resp.json()
        hits = data.get("hits", [])
        first_hit = next(iter(hits), None)
        
        if first_hit and "clinvar" in first_hit:
            clinvar = first_hit.get("clinvar", {})
            variant_id = clinvar.get("variant_id")
            rcv = clinvar.get("rcv", {})
            if isinstance(rcv, list):
                rcv = next(iter(rcv), {})
                
            sig = rcv.get("clinical_significance", "Classification Reported")
            review = rcv.get("review_status", "No assertion criteria")
            clinvar_url = f"https://nih.gov{variant_id}/" if variant_id else None
            accession_label = f"VCV{str(variant_id).zfill(9)}" if variant_id else "ClinVar Entry"
            return f"{sig} ({review})", clinvar_url, accession_label
            
    except Exception:
        return "ClinVar API Offline", None, None
        
    return "Not Documented in ClinVar (Novel / VUS)", None, None

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

    metrics = compute_physicochemical_metrics(pos, wt, alt, in_tm)
    clinvar_text, clinvar_url, accession_label = fetch_live_clinvar(wt, pos, alt)

    # Calibrate final pathogenicity % when ClinVar has an established classification
    prob_pct = metrics["prob_pct"]
    if "Pathogenic" in clinvar_text and "Benign" not in clinvar_text:
        prob_pct = max(prob_pct, 94)
    elif "Benign" in clinvar_text and "Pathogenic" not in clinvar_text:
        prob_pct = min(prob_pct, 12)

    if prob_pct >= 70 or ("Pathogenic" in clinvar_text and "Benign" not in clinvar_text):
        pred_bg, pred_border, pred_text = "#2d111a", "#fb7185", "#ffe4e6"
        verdict = "Likely Deleterious / Pathogenic"
    elif prob_pct <= 35 or ("Benign" in clinvar_text and "Pathogenic" not in clinvar_text):
        pred_bg, pred_border, pred_text = "#0b2920", "#10b981", "#d1fae5"
        verdict = "Likely Tolerated / Benign"
    else:
        pred_bg, pred_border, pred_text = "#2c1e0b", "#f59e0b", "#fef3c7"
        verdict = "Moderate Disruption / VUS"

    return {
        "domain": domain,
        "domain_style": {"bg": domain_bg, "border": domain_border, "text": domain_text},
        "prediction": clinvar_text,
        "verdict": verdict,
        "prob_pct": prob_pct,
        "metrics": metrics,
        "clinvar_url": clinvar_url,
        "accession_label": accession_label,
        "predict_style": {"bg": pred_bg, "border": pred_border, "text": pred_text}
    }
