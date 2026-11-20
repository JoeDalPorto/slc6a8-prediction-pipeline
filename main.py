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

def compute_fallback_heuristic(wt: str, alt: str, in_tm: bool):
    """Biophysical fallback when neither ClinVar nor AlphaMissense/REVEL has an entry."""
    delta_h = abs(HYDROPATHY.get(alt, 0.0) - HYDROPATHY.get(wt, 0.0))
    delta_q = abs(CHARGE.get(alt, 0.0) - CHARGE.get(wt, 0.0))
    helix_disrupt = in_tm and (alt in ("P", "G") or wt in ("P", "G"))

    score = 0.20 + min(delta_h * (0.06 if in_tm else 0.03), 0.30) + min(delta_q * 0.18, 0.28)
    if in_tm:
        score += 0.14
    if helix_disrupt:
        score += 0.22
    return int(min(max(round(score * 100), 5), 96))

def extract_first_float(val):
    """Safely pulls a single float from nested MyVariant list/dict fields."""
    if isinstance(val, (int, float)):
        return float(val)
    if isinstance(val, list) and val:
        return extract_first_float(val[0])
    if isinstance(val, dict):
        for v in val.values():
            res = extract_first_float(v)
            if res is not None:
                return res
    return None

@lru_cache(maxsize=2000)
def query_variant_databases(wt: str, pos: int, alt: str):
    """Queries MyVariant.info across all SLC6A8 transcripts for ClinVar + AlphaMissense + REVEL."""
    wt_3 = AA_1_TO_3.get(wt, wt)
    alt_3 = AA_1_TO_3.get(alt, alt)

    # Broadened wildcard query catches NP_005620.1:p.Gly42Arg, p.G42R, and dbNSFP records
    q = f"SLC6A8 AND (p.{wt_3}{pos}{alt_3} OR p.{wt}{pos}{alt})"
    url = "https://myvariant.info"
    params = {
        "q": q,
        "fields": "clinvar.variant_id,clinvar.rcv,dbnsfp.alphamissense,dbnsfp.revel",
        "size": 5
    }

    clinvar_text = "Not Documented in ClinVar (Novel / VUS)"
    clinvar_url = None
    accession_label = None
    am_score = None
    revel_score = None

    try:
        resp = requests.get(url, params=params, timeout=6)
        hits = resp.json().get("hits", [])

        for hit in hits:
            # 1. Extract ClinVar assertion if present
            if "clinvar" in hit and not clinvar_url:
                cv = hit["clinvar"]
                vid = cv.get("variant_id")
                rcv = cv.get("rcv", {})
                if isinstance(rcv, list):
                    rcv = next(iter(rcv), {})
                sig = rcv.get("clinical_significance")
                rev = rcv.get("review_status", "ClinVar record")
                if sig:
                    clinvar_text = f"{sig} ({rev})"
                    if vid:
                        clinvar_url = f"https://nih.gov{vid}/"
                        accession_label = f"VCV{str(vid).zfill(9)}"

            # 2. Extract AlphaMissense & REVEL scores from dbNSFP
            if "dbnsfp" in hit:
                db = hit["dbnsfp"]
                if am_score is None and "alphamissense" in db:
                    am_score = extract_first_float(db["alphamissense"].get("score"))
                if revel_score is None and "revel" in db:
                    revel_score = extract_first_float(db["revel"].get("score"))

    except Exception:
        clinvar_text = "Variant API Timeout / Offline"

    return clinvar_text, clinvar_url, accession_label, am_score, revel_score

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

    clinvar_text, clinvar_url, accession_label, am_score, revel_score = query_variant_databases(wt, pos, alt)

    # Prioritized Pathogenicity % Calibration:
    # 1. Established ClinVar Pathogenic/Benign assertion
    # 2. AlphaMissense / REVEL ensemble probability
    # 3. Biophysical heuristic fallback
    source_tag = "Biophysical Estimate"
    if "Pathogenic" in clinvar_text and "Benign" not in clinvar_text:
        prob_pct = 96
        source_tag = "ClinVar Confirmed"
    elif "Benign" in clinvar_text and "Pathogenic" not in clinvar_text:
        prob_pct = 6
        source_tag = "ClinVar Confirmed"
    elif am_score is not None and revel_score is not None:
        prob_pct = int(round(((am_score * 0.6) + (revel_score * 0.4)) * 100))
        source_tag = f"AlphaMissense ({am_score:.2f}) + REVEL ({revel_score:.2f})"
    elif am_score is not None:
        prob_pct = int(round(am_score * 100))
        source_tag = f"AlphaMissense ({am_score:.2f})"
    elif revel_score is not None:
        prob_pct = int(round(revel_score * 100))
        source_tag = f"REVEL ({revel_score:.2f})"
    else:
        prob_pct = compute_fallback_heuristic(wt, alt, in_tm)

    if prob_pct >= 70:
        pred_bg, pred_border, pred_text = "#2d111a", "#fb7185", "#ffe4e6"
        verdict = f"Pathogenic ({source_tag})"
    elif prob_pct <= 35:
        pred_bg, pred_border, pred_text = "#0b2920", "#10b981", "#d1fae5"
        verdict = f"Benign ({source_tag})"
    else:
        pred_bg, pred_border, pred_text = "#2c1e0b", "#f59e0b", "#fef3c7"
        verdict = f"Uncertain / VUS ({source_tag})"

    return {
        "domain": domain,
        "domain_style": {"bg": domain_bg, "border": domain_border, "text": domain_text},
        "prediction": clinvar_text,
        "verdict": verdict,
        "prob_pct": prob_pct,
        "clinvar_url": clinvar_url,
        "accession_label": accession_label,
        "predict_style": {"bg": pred_bg, "border": pred_border, "text": pred_text}
    }
