
SLC6A8-LeuT Variant Explorer (v1.0-beta)

An interactive, zero-dependency web application for mapping, visualizing, and evaluating missense variants across the canonical 635-amino-acid human creatine transporter (SLC6A8 / CRT1, UniProt P48029).
Live Interactive Tool: [github.io](https://joedalporto.github.io/slc6a8-prediction-pipeline/)

Key Features

• Interactive 12-TM LeuT-Fold Topology Model: Visualizes all 635 residues across transmembrane helices (TM1–TM12), extracellular loops (EL1–EL6), intracellular loops (IL1–IL5), and cytoplasmic tails. Click any residue node directly on the diagram or navigate via the sequence slider and -1 / +1 fine-tune steppers.

• Non-Contiguous Solute & Ion Binding Site Overlay: Highlights 18 critical functional residues across TM1, TM3, TM6, TM7, and TM8 that form the creatine substrate pocket and \(\text{Na}^{+}\)/\(\text{Cl}^{-}\) coordination network, with a quick-toggle switch in the header.

• ClinVar Reference & Direct NCBI Links: Cross-references curated SLC6A8 missense variants (Pathogenic, Benign, or VUS) and provides direct one-click links to matching NCBI ClinVar variation records and HGVS queries.

• Multi-Parameter Variant Scoring:

	• Predictive Score (0.000 – 1.000): Integrates curated clinical assertions with biophysical disruption modeling (Kyte-Doolittle hydropathy shift, formal charge change, helix-breaker penalties, and active-site involvement).
	
	• Folding Stability (\(\Delta\Delta G\), kcal/mol): Estimates thermodynamic folding disruption from steric volume clashes, membrane burial penalties, and helix-kink changes to help distinguish ER-retention/misfolding risks from fold-neutral catalytic variants.
	
	• Evolutionary Conservation (%): Position-specific conservation profile modeled across the SLC6 neurotransmitter transporter family.
	
• Session Variant Logger & CSV Export: Record evaluated mutations with custom tester notes during a session and export a 12-column .csv spreadsheet (Variant_1Letter, Variant_HGVS, Position, Structural_Domain, Transport_Site_Role, ClinVar_Status, Predictive_Score, Predictive_Verdict, Folding_ddG_kcal_mol, Stability_Impact, Conservation_Pct, Tester_Notes) for downstream analysis.

Architecture

The entire application is self-contained inside index.html using vanilla JavaScript, semantic CSS, and inline SVG rendering. No backend server, external CSS framework, or build step is required—enabling instant zero-latency execution directly in desktop and mobile browsers via GitHub Pages.

Repository Structure

• index.html — Complete standalone frontend, canonical UniProt P48029 sequence, 2D LeuT-fold SVG engine, biophysical scoring models, and CSV exporter.
• README.md — Project documentation and feature overview.

Disclaimer

This tool is intended for research and educational exploration only (v1.0-beta). Computational pathogenicity scores and folding stability estimates should not be used as standalone clinical diagnostic criteria.
