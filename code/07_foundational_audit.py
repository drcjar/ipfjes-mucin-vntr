# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

"""07_foundational_audit.py — verify foundational ASSIGNMENTS / ORIENTATIONS.

Complements verify_numbers.py. That script checks the data->code->manuscript NUMBER
chain but TRUSTS the haplogroup labels and allele orientations: a mislabel or allele
flip would still reconcile to 158/158. This audit instead re-derives each labelling,
orientation, and phenotype-coding decision and checks it against an INDEPENDENT
reference (gnomAD allele frequencies, Plender population frequencies, the OEM /
Reynolds 2023 paper, and biological expectation). A FLAG here means a foundational
assignment may be wrong even if every downstream statistic looks fine.
"""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings("ignore")
R = "."

rows = []
def audit(item, derived, expected, ok):
    rows.append((item, str(derived), expected, "OK" if ok else "FLAG"))

hap = pd.read_csv(f"{R}/muc5_haplogroups.csv")          # full genotyped panel (N=2912, UKB-European-dominated)
clin = pd.read_csv(f"{R}/flat_data_genotype_subset.csv"); clin["pid"] = clin.participant_id.astype(str)
hap["pid"] = hap.IID.str.split("_").str[0]
df = clin.merge(hap, on="pid").dropna(subset=["age", "ever_smoked", "peto_exposed"])
cases, ctrl = df[df.case == 1], df[df.case == 0]

# --- A. MUC5AC haplogroup tag assignment: frequency vs gnomAD/Plender (catches an H1/H3 swap) ---
h1, h2, h3 = (hap.H1_count.mean()/2, hap.H2_count.mean()/2, hap.H3_count.mean()/2)
audit("MUC5AC H1 freq (rs769768817 tag)", f"{h1:.2f}", "European H1 most common (~0.55)", 0.45 < h1 < 0.65)
audit("MUC5AC H3 freq (rs28542750 tag)", f"{h3:.2f}", "European H3 ~0.21 (gnomAD NFE A=0.21)", abs(h3-0.21) < 0.06)
audit("SWAP GUARD: H1 > H3 (H3 is minor/long clade)", f"H1={h1:.2f} vs H3={h3:.2f}", "H1>H3; H3<0.40 else H1/H3 swapped", h1 > h3 and h3 < 0.40)

# --- B. MUC5B discriminator orientation ---
b2 = hap.MUC5B_H2_count.mean()/2
audit("MUC5B H2 freq (common allele = H2)", f"{b2:.2f}", "Plender European H2 ~0.82", abs(b2-0.82) < 0.08)
rB = np.corrcoef(hap.MUC5B_H1_count, hap.rs35705950_T_count)[0, 1]
audit("rs35705950 T on H2 background: r(MUC5B_H1,T)<0", f"{rB:+.3f}", "negative (T on common H2)", rB < 0)

# --- C. rs35705950 coding AND case coding (joint orientation check) ---
tcc, tcn = 100*(cases.rs35705950_T_count > 0).mean(), 100*(ctrl.rs35705950_T_count > 0).mean()
audit("rs35705950 T-carrier: cases vs controls", f"{tcc:.0f}% vs {tcn:.0f}%",
      "cases>>controls (=> T=risk AND case=1=IPF)", tcc > tcn + 20)

# --- D. covariate definitions vs OEM / Reynolds 2023 ---
audit("Asbestos (Peto) cases vs controls", f"{100*cases.peto_exposed.mean():.0f}% vs {100*ctrl.peto_exposed.mean():.0f}%",
      "~67% vs 63% (OEM 2023; Table 1)", abs(100*cases.peto_exposed.mean()-66.8) < 4)
audit("Ever-smoked cases vs controls", f"{100*cases.ever_smoked.mean():.0f}% vs {100*ctrl.ever_smoked.mean():.0f}%",
      "~76% vs 69% (Table 1)", abs(100*cases.ever_smoked.mean()-75.7) < 4)
audit("Age: cases median", f"{cases.age.median():.0f}", "~77 (Table 1)", abs(cases.age.median()-77) <= 2)
audit("Cohort: cases / controls (genotyped IPFJES)", f"{len(cases)} / {len(ctrl)}", "428 / 401", len(cases) == 428 and len(ctrl) == 401)

# --- E. internal consistency ---
s = hap.H1_count + hap.H2_count + hap.H3_count
audit("MUC5AC counts sum to 2 per person", f"mean={s.mean():.3f}", "exactly 2", abs(s.mean()-2) < 0.01)

w = max(len(r[0]) for r in rows)
print(f"{'FOUNDATIONAL ASSIGNMENT':{w}}  {'DERIVED':22}  {'EXPECTED (independent ref)':38}  VERDICT")
for it, d, e, v in rows:
    print(f"{it:{w}}  {d:22}  {e:38}  {v}")
nflag = sum(v == "FLAG" for *_, v in rows)
print(f"\n{len(rows)} foundational checks · {len(rows)-nflag} OK · {nflag} FLAG")
print("NOTE: this audits LABELS/ORIENTATION (what verify_numbers trusts); the Plender Table S4"
      " r2 values remain the one external primary-source seal still to confirm.")
