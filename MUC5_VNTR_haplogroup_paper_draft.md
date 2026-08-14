# A *MUC5AC* H2 VNTR haplogroup association with idiopathic pulmonary fibrosis survives adjustment for the promoter variant rs35705950, whereas the *MUC5AC* H3 and *MUC5B* H1 associations do not: a tagging-SNP analysis of the IPF Job Exposures Study

**Authors:** Carl Reynolds¹²

**Affiliations:**
¹ Homerton University Hospital NHS Foundation Trust, London, England, UK
² Wolfson Institute of Population Health, Queen Mary University of London, London, England, UK

**Corresponding author:** Carl Reynolds, Homerton University Hospital NHS Foundation Trust, Homerton Row, London E9 6SR, UK. academic@carlreynolds.net

**Keywords:** idiopathic pulmonary fibrosis, *MUC5AC*, *MUC5B*, VNTR, haplogroup, gene-environment interaction, asbestos

**Take-home message (optional, for ERJ Open Research):** In a dedicated IPF case-control study, the *MUC5AC* H3 and *MUC5B* H1 associations were not separable from the *MUC5B* promoter variant rs35705950; a *MUC5AC* H2 haplogroup association did survive adjustment, including for every other independent signal at the locus, and warrants independent replication.

---

## Abstract

**Background:** A *MUC5B* promoter variant (rs35705950) is the strongest common genetic risk factor for idiopathic pulmonary fibrosis (IPF). The adjacent mucin genes *MUC5B* and *MUC5AC* contain variable number tandem repeat (VNTR) regions, resolved by long-read sequencing into SNP-tagged haplogroups that could influence airway mucus independently of rs35705950. Because the tagging SNPs lie in the same linkage disequilibrium block as rs35705950, whether an independent effect is resolvable with array data is unclear.

**Methods:** Using haplogroup tagging SNPs, we assigned *MUC5AC* (H1/H2/H3) and *MUC5B* (H1/H2) VNTR haplogroups from imputed genotype array chromosome 11 data in 798 IPF Job Exposures Study (IPFJES) participants (397 cases, 401 controls). Logistic regression adjusted for age, smoking, five principal components and centre, with and without rs35705950.

**Results:** *MUC5AC* H3 and *MUC5B* H1 were each associated with IPF (*MUC5AC* H3 OR 1.69, 95% CI 1.32–2.17, p=4.0×10⁻⁵; *MUC5B* H1 OR 0.72, 95% CI 0.52–0.99, p=0.041). Both tracked the rs35705950 T allele (r=0.36 and −0.10) and attenuated to non-significance on adjustment for it (H3 OR 1.17, 95% CI 0.89–1.54, p=0.27; *MUC5B* H1 OR 0.81, 95% CI 0.58–1.13, p=0.22). Adjusting rs35705950 for H3, by contrast, barely changed it (OR 3.98 to 3.77). A *MUC5AC* H2 haplogroup association persisted after adjustment (joint-model OR 0.66, 95% CI 0.49–0.89, p=0.006), and more weakly on conditioning for every other independent signal at the locus (OR 0.69, p=0.029).

**Conclusions:** In an IPF case-control study, a *MUC5AC* H2 haplogroup association persisted after adjustment for rs35705950 (OR 0.66, p=0.006), whereas the *MUC5AC* H3 and *MUC5B* H1 associations were not independent of it. The H2 signal is a single-cohort, unreplicated finding that needs individual-level replication; localising it to any co-inherited feature of the clade would require phased sequence.

---

## Introduction

Idiopathic pulmonary fibrosis (IPF) is a progressive and frequently fatal fibrotic lung disease (Raghu et al., 2022). The strongest common genetic risk factor is a gain-of-function promoter variant in the mucin gene *MUC5B* (rs35705950; minor allele frequency ~9% in Europeans). The minor (T) allele confers an odds ratio of approximately 4 per allele for IPF (Allen et al., 2017) and explains an estimated 2–11% of variation in IPF risk, against an estimated overall heritability of 28% (Moll et al., 2023). It increases *MUC5B* expression in the bronchiolar epithelium and is thought to predispose to fibrosis through mucociliary dysfunction (Seibold et al., 2011; Hancock et al., 2018) or epithelial cell stress (Dobrinskikh et al., 2023; Bridges et al., 2025).

Both *MUC5B* and *MUC5AC* lie within the 11p15.5 mucin gene cluster and encode gel-forming mucins with large central VNTR regions whose sequence and, at *MUC5AC*, repeat length varies between individuals. Plender et al. recently used long-read sequencing to characterise variation at both loci, defining three *MUC5AC* haplogroups (H1: short, H2: medium, H3: long) and two *MUC5B* haplogroups (H1: rare, H2: common), as well as identifying tagging SNPs (tSNPs) for haplogroup assignment in short-read datasets (Plender et al., 2024). The two classifications differ in kind. Both are phylogenetic clades defined from putatively neutral non-coding sequence: for *MUC5AC*, 5.6 kbp of intron plus 21 kbp of 3′-flanking intergenic sequence; for *MUC5B*, 26.5 kbp of intron. While the *MUC5AC* clades additionally correspond to distinct VNTR repeat lengths (H3 is ~670 amino acids longer than H1), the MUC5B protein is largely invariant in length (5,761–5,762 aa across common haplotypes). No length difference between the *MUC5B* haplogroups has been reported, so H1/H2 appear to index haplotype ancestry rather than repeat architecture.

There are biological grounds to ask whether VNTR haplogroup matters for IPF beyond rs35705950. The VNTR encodes the serine- and threonine-rich backbone carrying the O-glycans that set mucus viscoelasticity, and MUC5AC and MUC5B are the principal secreted determinants of mucociliary clearance (Plender et al., 2024; Hill et al., 2022; Abrami et al., 2024). Variation in VNTR structure could therefore act on that clearance mechanism (the route proposed to link rs35705950 to IPF; Seibold et al., 2011; Hancock et al., 2018) through protein architecture rather than expression. MUC5B is, moreover, the dominant secretory mucin in the distal airway epithelium implicated in IPF (Okuda et al., 2019), though it is *MUC5AC* whose haplogroups differ in repeat length.

Whether these long-read-defined VNTR haplogroups influence IPF risk has not been tested directly. In the *All of Us* phenome-wide scan, no *MUC5AC* haplogroup tag survived correction for any phenotype, and the *MUC5B* locus was represented there by the promoter variant rs35705950 rather than by a haplogroup tag (Plender et al., 2024). A genome-wide tandem-repeat study in IPF (Oketch et al., 2025) did report an associated repeat at 11p15.5 that lost significance once rs35705950 was included, but its catalogue was restricted to units of ≤20 bp: it indexed a shorter intergenic repeat, not the *MUC5AC* (24 bp) or *MUC5B* (87 bp) coding VNTRs, and did not test the haplogroup classification at all. The specific question left open is therefore whether the Plender haplogroups, or the coding-VNTR length they partly capture, carry any IPF risk separable from rs35705950. We address it in IPFJES (the IPF Job Exposures Study), a dedicated case-control study of occupational asbestos exposure in IPF (Reynolds et al., 2023), quantifying how much, if any, haplogroup association is independent of rs35705950.

---

## Methods

### Study population

IPFJES is a multicentre UK case-control study of men, recruiting from 21 hospitals in England, Scotland, and Wales (Reynolds et al., 2023). IPF cases were diagnosed according to the 2011 ATS/ERS/JRS/ALAT guidelines; recruitment ran from February 2017 to October 2019. Controls were participants without IPF, recruited from the same centres and frequency-matched for age. All participants provided written informed consent. Ethical approval was granted by the East Midlands – Nottingham 1 Research Ethics Committee (REC reference 17/EM/0021; IRAS project ID 203355). Of 902 IPFJES participants with clinical data, 829 also had typed-array and imputed chromosome 11 genotypes, with complete covariate data. One recruitment centre (n=31) enrolled cases only and caused complete separation when entered as a fixed effect, so it was excluded, giving an analytical sample of **798 participants (397 cases, 401 controls)** used for every analysis reported here.

### Genotyping and imputation

Participants were genotyped on the Biobank Affymetrix Axiom array (GRCh37). As previously described (Chin et al., 2026), sample-level quality control removed individuals with low genotyping call rate (<95%), high heterozygosity, sex mismatches, duplicates, and high relatedness. Imputation used the TOPMed reference panel (GRCh38); imputed variants with r²<0.5 or minor allele count ≤3 were excluded.

### Haplogroup assignment

*MUC5AC* VNTR haplogroup tagging SNPs (tSNPs) and their r² values were taken from Plender et al. Supplementary Table S4; allele orientation was verified against Table S6 (allele designations from *MUC5AC* upper-airway epithelial eQTL data), used for orientation only, as IPF risk direction was determined empirically. All 17 *MUC5AC* tSNPs were extracted from GRCh38 imputed data using PLINK v2. *MUC5AC* haplogroups were tagged by their concordance SNPs, each entered as an allele-dosage count: H3 by rs28542750 (r²=1.0) and H2 by rs879201833 (r²=0.888, the best available H2 tag), with H1 (the most common haplogroup, r²=0.899 for rs769768817) as the reference. An earlier parameterisation assigned H2 as the diplotype complement (2 − H1 − H3); the direct tag is used here as primary because it tags H2 directly rather than by subtraction, and the two agree closely (r=0.97) with comparable tagging quality (Supplementary Table S3). The joint model with the direct tag is therefore an approximate rather than exact H1-reference parameterisation (H2 and H3 are entered as dosages that do not sum to two). Tags were cross-checked against five H1-tagging and six H3-tagging SNPs (Supplementary Table S2). The haplogroup labels for these two tSNPs in Plender et al.'s main text appear transposed relative to their Supplementary Table S4. Their Table S4 resolves it, giving per-haplogroup r² of 0.111 (H1), 0.055 (H2) and 1.000 (H3) for rs28542750, and 0.899 (H1), 0.012 (H2) and 0.001 (H3) for rs769768817. We therefore assigned rs28542750 to H3 and rs769768817 to H1. Our own control frequencies corroborate this: rs28542750 (0.197) falls within the H3 tag set and matches the ~21% European H3 frequency, and rs769768817 (0.566) within the H1 tag set; this also fixes the orientation of both palindromes. The primary inference is in any case invariant to the labelling; only which clade is described as the long one would change.

The *MUC5B* H1/H2 discriminator (chr11:1,244,757 GRCh38; H2 r²=1.0, H1 r²=0.003) was taken from the Plender main text. Being a strand-ambiguous C/G variant, its orientation was resolved by frequency: the major allele (0.854 in controls) tags the common H2 haplogroup, and a minor-allele frequency of 0.146, far from 0.5, makes the assignment unambiguous despite the palindrome (consistent with Plender's overall H2 estimate of ~82%). This was corroborated by the inverse correlation of *MUC5B* H1 with the rs35705950 (chr11:1,219,991 GRCh38) risk allele (r=−0.10; see Results), consistent with the T allele being carried predominantly on the common *MUC5B* H2 background. All tSNPs passed Hardy-Weinberg testing in the IPFJES controls (minimum p=0.41); genotypes were hard-called from imputed data with no missingness. Their alleles, frequencies, and tagging r² are given in Supplementary Table S2.

### Covariates and exposure

Five genetic principal components, computed with PLINK v2 from LD-pruned typed-array data (314,996 variants) within the 798 analytical participants themselves, were used in all models. Occupational asbestos exposure was ascertained by Peto job-exposure matrix applied to lifetime occupational history (binary classification) (Reynolds et al., 2023); 519 of 798 participants (65.0%) were classified as exposed.

### Statistical analysis

Logistic regression (Python 3.12; statsmodels) was performed with IPF case/control status as outcome. Haplogroup counts were entered additively (0/1/2 copies); all models adjusted for age and ever-smoking (risk factors for IPF, included to improve precision of the case–control association rather than to control confounding of the germline genotype, which they cannot affect), PC1–5 (genetic ancestry), and recruitment centre (regression N=798).

**_MUC5AC_ models** (H1 as reference in joint model): (1a) *MUC5AC* H3 main effect; (1b) H3 + imputed rs35705950 T count; (1c) joint model: H2 and H3 entered simultaneously (H1 reference, the most common haplogroup at ~57%); (1d) additionally adjusted for rs35705950.

**_MUC5B_ model** (binary haplogroup; no joint model required): (2a) *MUC5B* H1 main effect; (2b) *MUC5B* H1 + rs35705950 T count. rs35705950 was also modelled alone for reference.

**Sensitivity analyses**: rs35705950-stratified analyses (non-carriers vs T-allele carriers) were performed for *MUC5AC* H2, H3 and *MUC5B* H1 (Supplementary Table S1). To confirm the H2 estimate does not depend on the choice of tag, the joint model was refitted with H2 assigned by the diplotype complement (2 − H1 − H3) and by the three panel H2 tSNPs (r²=0.862) as alternatives to the primary rs879201833 tag (Supplementary Table S3). A fifth Plender H2 tag, rs1201597847 (H2 r²=0.859), is present in the imputed data but adds nothing, being in complete LD with the panel tags. Nested models were compared by likelihood-ratio test, and rs35705950 was additionally estimated within H3-carrier strata. Because Plender et al. infer heterozygote advantage at H3 from balancing-selection signatures in 1000 Genomes Europeans, *MUC5AC* H3 was also fitted as a two-degree-of-freedom genotypic term. Exploratory haplogroup × asbestos and haplogroup × smoking interaction models were also fitted, and the same interaction tests were applied to rs35705950 (T-allele carriage) itself. As a further sensitivity analysis, *MUC5AC* VNTR repeat length was imputed into the 798 from the Bai et al. long-read reference panel using the ImputeSV server (SHAPEIT5 phasing, Minimac4), and both the 751 bp exon-31 length allele and a composite of total exon-31 repeat length were tested against IPF with and without adjustment for rs35705950 (Supplementary Table S5; `code/12_imputed_vntr_length.py`). As an agnostic test of the wider locus, every common (MAF≥0.01) well-imputed (imputation r²≥0.8) VNTR and structural-variant allele across the mucin cluster (*MUC6*–*MUC5B*; chr11:1.01–1.263 Mb; 109 variants after de-duplicating the VNTR and SV panels) was tested against IPF, marginally and conditional on rs35705950, with the cluster-wide threshold set from the effective number of independent tests (Li and Ji, 2005; Supplementary Table S6; `code/13_mucin_cluster_scan.py`). Separately, to test whether the *MUC5AC* H2 residual is distinct from the other independent signals at 11p15.5, the H2 tag was additionally conditioned on the eight conditionally-independent IPF signals reported at this locus by Kousathanas et al. (2025) (one in *MUC6*, three in *MUC2*, four in *MUC5B*; `code/14_h2_conditional_on_locus_signals.py`). Four primary haplogroup coefficients were tested (*MUC5AC* H3; H2 and H3 from the joint model; *MUC5B* H1); p-values are two-sided and uncorrected, and nominal results are described as such.

---

## Results

### Participant characteristics

Table 1 shows the participant characteristics. Cases were slightly older (median 77 vs 76 years, p=0.003) and more commonly ever-smokers (75.6% vs 68.8%, p=0.041). Occupational asbestos exposure was similar between cases and controls (67% vs 63.1%, p=0.278). The rs35705950 T allele was carried by 58.9% of cases versus 22.9% of controls (p=9.4×10⁻²⁵), consistent with the established association of this variant with IPF (Seibold et al., 2011).

### Haplogroup frequencies

In the IPFJES controls (n=401), which approximate a population reference since cases are enriched for risk haplotypes, *MUC5AC* haplogroup frequencies (H1 ~57%, H2 ~24%, H3 ~20% of haplotypes) and *MUC5B* frequencies (H2 ~85%, H1 ~15%) were consistent with European reference data (Plender et al., 2024): H1 was the most common *MUC5AC* haplogroup, H3 matched closely (~20% vs the reported European frequency of ~21%), and *MUC5B* H2 was consistent (~85% vs Plender's overall estimate of ~82%). *MUC5AC* H2 is tagged directly (rs879201833, r²=0.888); in our controls it had the intermediate frequency (~24%) expected for the middle *MUC5AC* clade, consistent with the reference data. The primary analysis relies on directly-tagged H3 (r²=1.0) and H2 (r²=0.888), with H1 as the reference baseline.

All three haplogroups were in Hardy-Weinberg equilibrium in controls (*MUC5AC* H3 p=0.65, *MUC5AC* H1 p=0.76, *MUC5B* H1 p=0.83), indicating no genotyping artefact in the tSNP proxies; the homozygote depletion Plender et al. reported in the general population was not replicated, likely reflecting limited power in the controls (n=401).

### Linkage disequilibrium between *MUC5AC* H3 and rs35705950

*MUC5AC* H3 count was correlated with rs35705950 T count (r=0.359, p=1.2×10⁻²⁵), sharing only ~13% of its variance with the risk allele. Among H3 carriers, 60.9% also carried the rs35705950 T allele. *MUC5B* H1, *MUC5AC* H1, and *MUC5AC* H2 were inversely correlated with rs35705950 T count (r=−0.100, −0.125, and −0.240), so the risk allele instead co-occurs with the *MUC5B* H2 and *MUC5AC* H3 backgrounds. Despite this correlation, H3-carrier and rs35705950-carrier status were discordant in 250 of 798 participants (31%; cases 123/397, controls 127/401). The corresponding figure for *MUC5B* H1 was 390 (49%); because *MUC5B* H1 is inversely correlated with the T allele, carrier discordance is not comparable across the two loci, and the ranking rests on |r| (0.359 vs 0.100). These LD values are compared with the 1000 Genomes European reference in Supplementary Table S7: only H3 has appreciable LD with rs35705950, and its estimate in the IPFJES controls (r²=0.09) matches the reference (0.10), so the LD is representative rather than a feature of this cohort (the full-sample estimates above are modestly larger because cases are enriched for risk haplotypes).

### Haplogroup associations with IPF

Full estimates are in Table 2. *MUC5AC* H3 was positively associated with IPF (OR 1.69, 95% CI 1.32–2.17; p=4.0×10⁻⁵), but this association was not independent of rs35705950. Adjusting for the promoter variant attenuated the H3 estimate to the null (OR 1.17, 95% CI 0.89–1.54; p=0.27), while rs35705950 itself retained a strong independent effect (OR 3.77, 95% CI 2.75–5.16; p=1.4×10⁻¹⁶).

The three-haplotype joint model (H1 reference) showed a similar pattern: H3 was enriched and H2 depleted among cases, but after additional adjustment for rs35705950 the H3 signal was no longer detectable while an inverse H2 association persisted (OR 0.66, 95% CI 0.49–0.89; p=0.006). This H2 estimate was stable across tag choices, agreeing with the diplotype complement and the three panel H2 tSNPs (Supplementary Table S3). Fitting H2 genotypically rather than additively was consistent with a simple dose-response (one copy OR 0.65, two copies OR 0.47), with no departure from the additive model (likelihood-ratio p=0.83). The H3 main-effect estimate uses the perfect-concordance H3 tag (r²=1.0) and so does not depend on the imperfectly tagged H1 baseline.

*MUC5B* H1 showed an inverse association with IPF (OR 0.72, 95% CI 0.52–0.99; p=0.041) that attenuated after rs35705950 adjustment (OR 0.81, 95% CI 0.58–1.13; p=0.22).

In rs35705950 non-carriers the H3 estimate was directionally consistent but underpowered to exclude a modest independent effect (N=468; OR 1.23, 95% CI 0.85–1.79; p=0.27; Supplementary Table S1). Fitting H3 genotypically rather than additively gave ORs of 1.95 (95% CI 1.43–2.66) for one copy and 2.03 (1.06–3.91) for two, but only 46 participants carry two copies and the genotypic model did not improve fit (likelihood-ratio p=0.13); the additive model predicts a two-copy OR of 2.86, which that interval contains. These data therefore do not distinguish an additive from a dominant model. Both terms attenuated after rs35705950 adjustment (1.26 and 1.13). Together with the non-replicated homozygote depletion (H3 Hardy-Weinberg p=0.65 in controls), these data neither support nor exclude the heterozygote-advantage model Plender et al. infer from balancing-selection signatures.

### Cluster-wide scan of imputed structural variants

The haplogroup analysis covers only *MUC5AC* and *MUC5B*, but 11p15.5 is a whole mucin cluster in which a much larger meta-analysis places independent IPF signals near *MUC2* and *MUC6* as well, not only *MUC5B* (Kousathanas et al., 2025). To ask whether any mucin structural feature beyond rs35705950 therefore carries IPF risk, we tested every common, well-imputed VNTR and structural-variant allele across the cluster (*MUC6*–*MUC5B*; 109 variants at imputation r²≥0.8), including the large *MUC2* VNTR domain that no haplogroup analysis reaches. Twenty-two reached cluster-wide significance marginally (the strongest by p-value a *MUC2* allele, OR 2.76, p=3×10⁻¹²), but every one tracked linkage disequilibrium with rs35705950 (that *MUC2* allele at r=0.53 with the T allele, exceeding even *MUC5AC* H3 at r=0.36), and none survived adjustment for it (smallest adjusted p=0.003, driven by a low-frequency variant with an unstable estimate and within chance expectation across the 54 effective independent tests; Supplementary Table S6). In this sample, then, no imputed structural variant in the cluster reached significance independently of the promoter variant, though this is a power-limited negative: the *MUC2*, *MUC6* and *MUC5B* signals noted above are ones 397 cases cannot detect (see Supplementary Table S6).

---

## Discussion

We tested whether *MUC5AC* or *MUC5B* VNTR haplogroup is associated with IPF independently of the *MUC5B* promoter variant rs35705950. To our knowledge this is the first such test in a dedicated IPF case-control study. We found that *MUC5AC* H3 is associated with IPF risk, but that this association is largely explained by its linkage disequilibrium with rs35705950, leaving no detectable residual H3 signal, as was also true of *MUC5B* H1. A residual *MUC5AC* H2 association did survive adjustment.

Adjusting rs35705950 for *MUC5AC* H3 barely changes it (OR 3.98 to 3.77), whereas adjusting H3 for rs35705950 abolishes it (OR 1.69 to 1.17). That is the expected signature of one variant carrying the association and the other tagging it. Formally, adding H3 to a model already containing rs35705950 does not improve fit (likelihood-ratio χ²₁=1.21, p=0.27), whereas adding rs35705950 to a model containing H3 improves it decisively (χ²₁=76.8, p=1.9×10⁻¹⁸); and rs35705950 retains a similar effect in both H3 strata: OR 3.29 (95% CI 2.11–5.15) among the 450 participants carrying no H3 copy, and 4.46 (2.77–7.18) among the 348 who do, with no interaction (p=0.41). *MUC5AC* H2 behaves differently, and the contrast is informative. Adjusting H2 for rs35705950 does not abolish it (OR 0.56 to 0.66); adding H2 to a model that already contains rs35705950 still improves fit (χ²₁=7.60, p=0.006), while rs35705950 itself barely changes (OR 3.77 to 3.58). Both are retained, the signature of two partially-independent haplogroup effects rather than one tagging the other.

The LD between *MUC5AC* H3 and rs35705950 (r=0.36) is the key interpretive constraint. *MUC5AC* and *MUC5B* lie within an ~100 kb mucin cluster at 11p15.5, where regional LD is well established in Europeans. Array-based tSNP proxies cannot resolve whether the H3 association reflects an independent VNTR length effect, rs35705950-mediated transcriptional dysregulation, or co-inheritance with another functional variant. The joint-model ordering (H2 depleted, H1 intermediate, H3 enriched for the rs35705950 T allele) is what would be expected if risk tracks rs35705950 frequency by haplogroup rather than VNTR length, and after adjustment the H3 signal was no longer detectable (OR 1.07, p=0.63; 1.17 when H3 is adjusted for rs35705950 without H2 in the model, as in Table 2).

H2, however, retained significance after adjustment (OR 0.66, 95% CI 0.49–0.89; p=0.006). H2 is directly tagged here (rs879201833, r²=0.888), and the estimate is stable across the alternative panel H2 tSNPs and the diplotype complement (Supplementary Table S3), so it is not an artefact of how H2 is assigned. Estimates are of similar direction and magnitude in both rs35705950 strata, with no heterogeneity between them (interaction p=0.93; Supplementary Table S1). The residual is also not accounted for by the other independent signals at 11p15.5: Kousathanas et al. resolved several conditionally-independent IPF signals across the cluster (in *MUC6*, *MUC2* and *MUC5B*) in a much larger meta-analysis (Kousathanas et al., 2025), and the H2 tag is only weakly correlated with each (|r|≤0.33). Conditioning the H2 tag on all eight alongside rs35705950 leaves the point estimate essentially unchanged (OR 0.66 to 0.69), the higher p-value reflecting the precision lost to eight added covariates in this sample (p=0.029). The estimate is thus robust to every independent signal so far identified at the locus, though its significance falls to nominal under this full conditioning.

Two cautions follow. The residual emerged from a single ~400-case cohort within an extensive set of haplogroup tests, so its magnitude may be inflated, and, by the same identifiability argument, robustness to the *known* signals does not exclude co-inheritance with a variant not yet resolved. It should be read as hypothesis-generating, requiring an individual-level conditional replication to confirm. The primary and panel-tag estimates (p=0.003–0.006) clear a Bonferroni threshold across the four primary haplogroup coefficients tested, though the diplotype complement (p=0.020) does not, and under full conditioning on the other locus signals the estimate falls to nominal (p=0.029); whether it survives correction therefore depends on the measure and the model. The stronger caution is mechanistic rather than statistical. Plender et al.'s external validation against Locityper in 2,600 genomes found 91% complete concordance, and three-quarters of the single-haplotype discordance fell on the H1-versus-H2 distinction: roughly 6% of genomes (≈3% of chromosomes) misclassified precisely along the axis this residual describes. The H2 tag cannot exclude that error, since it sits in the same linkage disequilibrium block. Residual confounding by rs35705950 is possible in principle, since it is imputed rather than directly typed; but its imputation quality at this locus is high, so any such effect is small.

If the H2 depletion is real, its shape argues against a simple repeat-count mechanism. The three clades are ordered by protein length (H1 5,654 aa, H2 5,742 aa, H3 6,325 aa), but the risk ordering is not. After adjustment for rs35705950, H3 is indistinguishable from the H1 reference (OR 1.07, p=0.63) while H2 remains distinctly protective (0.66); the intermediate-length clade, not an end of the range, is the outlier. A monotone effect of repeat number would give a gradient across the three; this is not one. That reading depends, however, on the clade least reliably tagged. Plender et al.'s validation places most tag-versus-sequence discordance on the H1-versus-H2 boundary, so even the direct H2 tag cannot exclude a systematic misassignment along that axis, and the non-monotonic shape is precisely what such an error would distort first.

Repeat number itself can also be looked up directly. Bai et al. imputed VNTR repeat counts genome-wide from long-read-derived reference panels into UK Biobank and released the association statistics; in those statistics the *MUC5AC* exon-31 repeat is not associated with the pulmonary-fibrosis phecode (p=0.65 and 0.84 at the two loci with the 24-bp *MUC5AC* repeat unit; 347 cases, N=349,660) (Bai et al., 2026; their genome-wide scan does not single out this locus). That analysis uses a heterogeneous phecode rather than adjudicated IPF and does not condition on rs35705950, so it does not answer the question posed here; but it is an independent and far larger test, and it agrees with the clade ordering above in suggesting that whatever distinguishes these haplotypes is not simply how many repeats they carry. That null should not be read as closing the question: the one study to size this VNTR directly (Guo et al., 2011, by Southern blot in cystic fibrosis) did associate a 6.4 kbp allele with severe lung disease while finding no SNP proxy that reproduced it, a reminder that a measured and an inferred length need not agree, though cystic fibrosis is not IPF.

We were able to test imputed length directly in IPFJES, rather than resting on inference. Imputing repeat number from the Bai et al. long-read panel into our data recovered the *MUC5AC* exon-31 length variants at high fidelity (a 751 bp exon-31 insertion imputed at r²=0.98). However, imputed length is equivalent to the clade tag: it correlates with H3 count at r=0.98 and carries the same weak linkage with rs35705950 (r=0.37). As a single variable it behaves identically to the H3 tag (unadjusted OR 1.67, 95% CI 1.29–2.15; adjusted for rs35705950 OR 1.13, 0.85–1.49), and neither adds to the other (likelihood-ratio p=0.76 and 0.24). This equivalence is not incidental to our sample. A single SNP tags H3 at r²=1.0 in Plender's panel (rs28542750, the concordance SNP from which we assign H3), so little haplotype information remains for imputation to add beyond that tag, and it recovers essentially the tag. Imputation can outperform single-SNP tagging where no SNP is in high LD with the target; where a near-perfect tag already exists, as at H3, the two coincide. A composite of total exon-31 repeat length gave the same result: OR 1.46 per additional kilobase unadjusted (1.12–1.90; p=0.006), abolished by rs35705950 (0.99, 0.73–1.32; p=0.92), and stable whether length was counted in base pairs or 24 bp repeat units (Supplementary Table S5). It too is null, and a higher-fidelity length cannot resolve the identifiability problem, because the imputed length carries the same LD with rs35705950 as the tags do: the barrier is co-inheritance, not crude measurement.

Several features of these haplotypes vary independently of repeat count and would be invisible to a length-based model. Plender et al. report that *MUC5AC* haplotypes differ in cysteine-rich domain copy number as well as VNTR copy number, and suggest that additional cys domains may enhance the non-covalent self-interactions that make the mucus gel impermeable. They also find that, normalised for VNTR length, the two shortest protein groups carry the highest serine and threonine density, so glycosylation density and absolute length dissociate. Motif usage varies as well, with 211 distinct 24-mers across their haplotypes. Any of these could differ between clades without tracking length, and none is measured by a tagging SNP. A distinct and arguably stronger functional candidate is the set of predicted protein-truncating VNTR alleles Plender et al. describe, which terminate the protein early to yield shorter isoforms. These recur on more than one haplogroup (most often *MUC5AC* H1) rather than defining a clade, so a haplogroup tag cannot represent them, and establishing whether they carry the risk would require genotyping and analysing the truncating alleles directly rather than the clade they sit on (Plender et al., 2024). Such an effect might also be contingent rather than constitutive, mattering only under the inflammatory load that drives *MUC5AC* expression, which the gene–environment tests reported here were far too small to address. The alternative remains that the residual reflects nothing architectural at all, but regulatory or other variation co-inherited on the H2 clade, a distinction these data cannot make.

The *MUC5B* H1 inverse association (OR 0.72) is largely explained the same way: H2, the predominant *MUC5B* haplogroup, co-occurs with rs35705950 T, making H1 appear protective by contrast, and the estimate attenuates on adjustment (0.72 to 0.81, p=0.22). As with the other haplogroups, however, this is an underpowered term rather than a demonstrated null: its estimate in rs35705950 non-carriers (0.71; Supplementary Table S1) is little changed from the marginal, so a modest residual cannot be excluded.

At *MUC5B*, as noted in the Introduction, the haplogroups index lineage rather than length, so a *MUC5B* haplogroup effect would more likely reflect co-inherited regulatory variation. The repeat-length imputation reported above accordingly has no *MUC5B* counterpart: the only *MUC5B* length variant Plender et al. resolve, the expanded central-exon haplotypes, is confined to individuals of African descent and near-absent here (imputed allele frequency <0.001), so neither a tagging SNP nor an imputed length indexes repeat number there.

Exploratory haplogroup × environment and rs35705950 × environment interactions were null on the multiplicative scale (Supplementary Table S4). On the additive scale, with each genetic factor oriented to its risk-increasing direction as additive interaction requires (Knol and VanderWeele, 2012), most contrasts were also null, but a nominal positive interaction with ever-smoking appeared for the risk-increasing background at two factors: rs35705950 (RERI +2.78, 95% CI +0.37 to +5.19) and non-carriage of the protective *MUC5AC* H2 haplogroup (+1.03, +0.10 to +1.97). With two nominal signals among many underpowered contrasts, and none on the multiplicative scale, these are hypothesis-generating at most and reported for completeness and future meta-analysis rather than as evidence of interaction.

### Strengths and limitations

Strengths include the use of validated tSNPs from a long-read sequencing study and well-characterised occupational exposure data from a study designed for gene×environment analysis. The principal limitation is that array-based tSNPs index a phylogenetic clade within an LD block: even where a haplogroup association survives adjustment, as *MUC5AC* H2 does, it identifies the clade but not which co-inherited feature drives it, and resolving that would need phased long-read sequence. Statistical power is a second limitation: where a haplogroup association attenuates, this sample cannot exclude a modest residual effect (the adjusted *MUC5AC* H3 interval, 0.89–1.54, remains compatible with a per-allele OR up to about 1.5). A further limitation applies at both loci: rs35705950 is not wholly external to the sequence defining the haplogroups. For *MUC5B*, the promoter variant lies within the gene under test and sits on the H2 background; for *MUC5AC*, the 21 kbp window Plender et al. used to build the clades spans rs35705950, which sits 18.9 kbp beyond the *MUC5AC* 3′ end. The effect is slight for *MUC5AC* (one common variant contributes little to a tree built from 27 kbp), but conditioning on rs35705950 asks the narrower question of whether clade membership matters among chromosomes matched on promoter genotype. Two further points bear on generalisability. IPFJES recruited men only, so applicability to women is unknown. And the sample was not restricted by genetic ancestry: it is predominantly of white ethnicity (96.6%, 771/798), with ancestry handled by principal-component adjustment; restricting to the 771 white participants left every estimate essentially unchanged (*MUC5AC* H3 1.11 adjusted for rs35705950, H2 0.69, *MUC5B* H1 0.71; `code/16_ancestry_and_promoter_gxe.py`), so residual population structure is unlikely to drive the results. One caveat for tagging remains: Plender's r² values are averaged across a 206-haplotype panel spanning five super-populations, so they carry no guarantee of transferring to the European majority here, and resolving *MUC5B* length variation, which is largely African-specific (above), would need trans-ancestry, isoform-level data rather than this cohort's array genotypes.

### Conclusions

*MUC5AC* H3 and *MUC5B* H1 are associated with IPF risk in directions consistent with their respective LD with rs35705950, and neither association is separable from it. A *MUC5AC* H2 haplogroup association does survive adjustment, including for every other independent signal so far identified at the locus, and is stable across the H2 tagging SNPs, though a systematic H1-versus-H2 tagging error shared across tags cannot be fully excluded. This residual is a single-cohort, nominally significant, as-yet-unreplicated finding whose magnitude may be inflated; it is hypothesis-generating and needs independent individual-level replication before any co-inherited feature of the clade could be held responsible.

---

## Tables

**Table 1. Participant characteristics**

| Characteristic | Cases (n=397) | Controls (n=401) | p-value |
|----------------|:-------------:|:----------------:|:-------:|
| Age, median (IQR), years | 77.0 (73.0–82.0) | 76.0 (70.0–81.0) | 0.003 |
| Ever-smoked, n (%) | 300 (75.6%) | 276 (68.8%) | 0.041 |
| Pack-years (smokers), median (IQR) | 20.0 (9.0–39.0) | 18.5 (6.0–35.0) | 0.069 |
| Asbestos-exposed (Peto), n (%) | 266 (67.0%) | 253 (63.1%) | 0.278 |
| **_MUC5AC_ haplogroup** | | | |
| H1 carrier, n (%) | 316 (79.6%) | 327 (81.5%) | 0.544 |
| H2 carrier, n (%) | 103 (25.9%) | 168 (41.9%) | 2.0×10⁻⁶ |
| H3 carrier, n (%) | 207 (52.1%) | 141 (35.2%) | 1.9×10⁻⁶ |
| **_MUC5B_ haplogroup** | | | |
| H1 carrier, n (%) | 81 (20.4%) | 109 (27.2%) | 0.030 |
| **rs35705950 (*MUC5B* promoter)** | | | |
| T carrier, n (%) | 234 (58.9%) | 92 (22.9%) | 9.4×10⁻²⁵ |
| TT homozygote, n (%) | 26 (6.5%) | 4 (1.0%) | 8.3×10⁻⁵ |

p-values: Mann-Whitney U (continuous); chi-squared (categorical). All statistics use the 798-participant analytical sample.

---

**Table 2. *MUC5AC* and *MUC5B* haplogroup associations with IPF**

| Model | Term | OR | 95% CI | p |
|-------|------|----|--------|:-:|
| MUC5AC H3, unadjusted | H3 count | 1.69 | 1.32–2.17 | 4.0×10⁻⁵ |
| MUC5AC H3, adjusted for rs35705950 | H3 count | 1.17 | 0.89–1.54 | 0.27 |
| | rs35705950 T count | 3.77 | 2.75–5.16 | 1.4×10⁻¹⁶ |
| MUC5AC H2, adjusted for H3 | H2 count | 0.56 | 0.42–0.75 | 8.3×10⁻⁵ |
| MUC5AC H2, adjusted for H3 + rs35705950 | H2 count | 0.66 | 0.49–0.89 | 0.006 |
| MUC5B H1, unadjusted | H1 count | 0.72 | 0.52–0.99 | 0.041 |
| MUC5B H1, adjusted for rs35705950 | H1 count | 0.81 | 0.58–1.13 | 0.22 |
| | rs35705950 T count | 3.93 | 2.91–5.31 | 6.0×10⁻¹⁹ |
| rs35705950 alone (reference) | rs35705950 T count | 3.98 | 2.95–5.37 | 2.2×10⁻¹⁹ |

All models adjusted for age, ever-smoking, PC1–5, and recruitment centre (N=798; one cases-only centre excluded). The "*MUC5AC* H2, adjusted for H3" rows are the joint model (H2 and H3 entered together, H1 reference), reporting the H2 term; in that model H3 is enriched and attenuates on rs35705950 adjustment as in the H3 main-effect rows above. rs35705950 odds ratios differ between rows because each is conditioned on the haplogroup term(s) in that model. Every haplogroup term is an allele-dosage count (0/1/2): the H2 term is the count of the direct H2 tag rs879201833, and the H3 and *MUC5B* H1 terms are the counts of their concordance SNPs.

---

## References

Abrami M, Biasin A, Tescione F, et al. Mucus structure, viscoelastic properties, and composition in chronic respiratory diseases. *Int J Mol Sci*. 2024;25(4):1933. https://doi.org/10.3390/ijms25031933

Allen RJ, et al. Genetic variants associated with susceptibility to idiopathic pulmonary fibrosis in people of European ancestry: a genome-wide association study. *Lancet Respir Med*. 2017;5(11):869–880.

Bai WY, Liu S, Duan Z, et al. Genome-wide associations of structural variants with human traits through imputation from long-read assemblies. *Nat Genet*. 2026;58(6):1258–1267. https://doi.org/10.1038/s41588-026-02612-z

Bridges JP, Vladar EK, Kurche JS, et al. Progressive lung fibrosis: reprogramming a genetically vulnerable bronchoalveolar epithelium. *J Clin Invest*. 2025;135(1):e183836. https://doi.org/10.1172/JCI183836

Chin D, Hernandez-Beeftink T, Donoghue L, et al. Genome-wide association study of idiopathic pulmonary fibrosis susceptibility using clinically-curated European-ancestry datasets. *Eur Respir J*. 2026; in press.

Dobrinskikh E, Hennessy CE, Kurche JS, et al. Epithelial endoplasmic reticulum stress enhances the risk of Muc5b-associated lung fibrosis. *Am J Respir Cell Mol Biol*. 2023;68(1):62–74. https://doi.org/10.1165/rcmb.2022-0252OC

Guo X, Pace RG, Stonebraker JR, et al. Mucin variable number tandem repeat polymorphisms and severity of cystic fibrosis lung disease: significant association with *MUC5AC*. *PLoS One*. 2011;6(10):e25452. https://doi.org/10.1371/journal.pone.0025452

Hancock LA, et al. *MUC5B* overexpression causes mucociliary dysfunction and enhances lung fibrosis in mice. *Nat Commun*. 2018;9:5363.

Hill DB, Button B, Rubinstein M, Boucher RC. Physiology and pathophysiology of human airway mucus. *Physiol Rev*. 2022;102(4):1757–1836. https://doi.org/10.1152/physrev.00004.2021

Knol MJ, VanderWeele TJ. Recommendations for presenting analyses of effect modification and interaction. *Int J Epidemiol*. 2012;41(2):514–520.

Kousathanas A, Odhams CA, Cook J, et al. Common and rare variant analyses reveal genetic factors underlying idiopathic pulmonary fibrosis and its shared aetiology with severe COVID-19. *eBioMedicine*. 2025;122:106048. https://doi.org/10.1016/j.ebiom.2025.106048

Li J, Ji L. Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix. *Heredity*. 2005;95(3):221–227. https://doi.org/10.1038/sj.hdy.6800717

Moll M, Peljto AL, Kim JS, et al. A polygenic risk score for idiopathic pulmonary fibrosis and interstitial lung abnormalities. *Am J Respir Crit Care Med*. 2023;208(7):791–801.

Oketch JW, et al. Genome-wide analysis of tandem repeat variation identifies *SLC15A4* as a susceptibility gene for idiopathic pulmonary fibrosis. *medRxiv* (preprint). 2025. https://doi.org/10.1101/2025.11.25.25340948

Okuda K, Chen G, Subramani DB, et al. Localization of secretory mucins *MUC5AC* and *MUC5B* in normal/healthy human airways. *Am J Respir Crit Care Med*. 2019;199(6):715–727. https://doi.org/10.1164/rccm.201804-0734OC

Plender EG, Prodanov T, Hsieh P, et al. Structural and genetic diversity in the secreted mucins *MUC5AC* and *MUC5B*. *Am J Hum Genet*. 2024;111:1700–1716. https://doi.org/10.1016/j.ajhg.2024.06.007

Raghu G, et al. Idiopathic pulmonary fibrosis (an update) and progressive pulmonary fibrosis in adults: an official ATS/ERS/JRS/ALAT clinical practice guideline. *Am J Respir Crit Care Med*. 2022;205(9):e18–e47.

Reynolds CJ, Sisodia R, Barber C, et al. What role for asbestos in idiopathic pulmonary fibrosis? Findings from the IPF job exposures case-control study. *Occup Environ Med*. 2023;80:97–103. https://doi.org/10.1136/oemed-2022-108404

Seibold MA, et al. A common *MUC5B* promoter polymorphism and pulmonary fibrosis. *N Engl J Med*. 2011;364(16):1503–1512.

---

**Word count (main text: Introduction through Conclusions, including Methods; excluding abstract, tables, and references):** ~4,950 words, ~1,950 over the ERJ Open Research 3,000-word limit (whitespace-split; token-based counts give ~3,300–3,400, so treat as advisory). The Discussion has been tightened during review; the remaining trim to the limit is best done in the shared Google Doc. The Discussion still carries most of the excess.

**Abstract word count:** ~272 words (structured), over the 250-word limit; trim in final edit.

*Analysis code, the exact genotype-handling commands (`code/REPRODUCE.md`), and four independent verification layers are available at https://github.com/drcjar/ipfjes-mucin-vntr. Individual-level participant data are held under IPFJES governance arrangements and are available to qualified researchers upon application to the corresponding author and completion of a data access agreement.*

---

## Supplementary methods: software and reproducibility

Analyses used PLINK v2.00a3.7 for LD-pruning, PCA, dosage export and Hardy-Weinberg testing, and PLINK v1.90b6.26 for single-variant text extractions, with Python 3.12 (numpy, scipy, pandas, statsmodels, pandas-plink) for statistical analysis. Genetic principal components were computed from LD-pruned typed-array data (`--indep-pairwise 500 50 0.2`, then `--pca`); haplogroup-tagging SNP dosages were extracted from imputed chromosome 11 (`--export A`); and tag-SNP Hardy-Weinberg equilibrium and missingness were tested in controls (`--hardy`). The exact commands, software versions, and the analysis and numerical-verification scripts are provided in the code repository (`code/REPRODUCE.md` and `code/README.md`).

---

## Supplementary Tables

**Supplementary Table S1. *MUC5AC* H3 and *MUC5B* H1 associations stratified by rs35705950 carrier status**

| Haplogroup | Stratum | N | OR | 95% CI | p |
|------------|---------|---|----|--------|---|
| MUC5AC H3 | Non-carrier (T=0) | 468 | 1.23 | 0.85–1.79 | 0.27 |
| MUC5AC H3 | Carrier (T≥1) | 320 | 1.08 | 0.70–1.67 | 0.72 |
| MUC5AC H3 | Overall | 798 | 1.69 | 1.32–2.17 | 4.0×10⁻⁵ |
| MUC5B H1 | Non-carrier (T=0) | 468 | 0.71 | 0.46–1.10 | 0.12 |
| MUC5B H1 | Carrier (T≥1) | 320 | 0.98 | 0.53–1.79 | 0.94 |
| MUC5B H1 | Overall | 798 | 0.72 | 0.52–0.99 | 0.041 |
| MUC5AC H2 | Non-carrier (T=0) | 468 | 0.65 | 0.45–0.94 | 0.021 |
| MUC5AC H2 | Carrier (T≥1) | 320 | 0.56 | 0.32–0.97 | 0.040 |
| MUC5AC H2 | Overall | 798 | 0.66 | 0.49–0.89 | 0.006 |

All models adjusted for age, ever-smoking, PC1–5, and recruitment centre. Non-carrier stratum excludes centre 18 (n=4, cases only); carrier stratum excludes centres 17 (n=2) and 21 (n=4) (cases only), to avoid complete separation. H3 OR in non-carriers (1.23) is consistent with confounding but the analysis is underpowered to exclude an independent effect of modest magnitude. *MUC5AC* H2 rows are from the joint model (H2 tag and H3 entered together, H1 reference). The "Overall" column is not uniform: for H3 and *MUC5B* H1 it is the unadjusted single-term estimate, whereas for H2 it is the rs35705950-adjusted joint estimate (0.66; Table 2), shown because H2's relevance is precisely that it persists after adjustment. The stratified estimates control for rs35705950 by design in every row.

---

**Supplementary Table S2. Haplogroup-tagging SNP panel: alleles, frequencies, and tagging quality**

\scriptsize

| tSNP | Position | Alleles | Counted | Freq | MAF | Tags | r² |
|--------------|---------|:---:|:---:|:--:|:--:|-----------|:--:|
| rs28542750 | chr11:1165501 | T/A | T | 0.803 | 0.197 | MUC5AC H3 conc. | 1.0 |
| rs769768817 | chr11:1195214 | G/C | C | 0.566 | 0.434 | MUC5AC H1 conc. | 0.899 |
| rs2075842 | chr11:1193830 | G/C | G | 0.569 | 0.431 | MUC5AC H1 | 0.923 |
| rs1132433 | chr11:1194354 | A/G | A | 0.569 | 0.431 | MUC5AC H1 | 0.923 |
| rs1132434 | chr11:1195265 | G/C | C | 0.567 | 0.433 | MUC5AC H1 | 0.923 |
| rs28652890 | chr11:1195858 | T/C | T | 0.567 | 0.433 | MUC5AC H1 | 0.923 |
| rs879136008 | chr11:1182048 | A/T | T | 0.569 | 0.431 | MUC5AC H1 | 0.923 |
| rs28519516 | chr11:1202343 | T/C | C | 0.759 | 0.241 | MUC5AC H2 | 0.862 |
| rs28558973 | chr11:1202859 | A/G | G | 0.759 | 0.241 | MUC5AC H2 | 0.862 |
| rs28368633 | chr11:1202914 | T/C | C | 0.759 | 0.241 | MUC5AC H2 | 0.862 |
| rs36154966 | chr11:1178216 | C/G | C | 0.809 | 0.191 | MUC5AC H3 | 1.0 |
| rs1004828576 | chr11:1178023 | A/G | A | 0.803 | 0.197 | MUC5AC H3 | 1.0 |
| rs940158763 | chr11:1177688 | A/T | A | 0.803 | 0.197 | MUC5AC H3 | 1.0 |
| rs36151150 | chr11:1177356 | T/C | T | 0.803 | 0.197 | MUC5AC H3 | 1.0 |
| rs36132281 | chr11:1177269 | T/C | T | 0.803 | 0.197 | MUC5AC H3 | 1.0 |
| rs35779873 | chr11:1197482 | C/T | C | 0.805 | 0.195 | MUC5AC H3 | 0.805 |
| rs35705950 | chr11:1219991 | G/T | G | 0.880 | 0.120 | *MUC5B* promoter | n/a |
| rs879201833 | chr11:1182279 | A/G | A | 0.763 | 0.237 | MUC5AC H2 direct | 0.888 |
| MUC5B discriminator | chr11:1244757 | C/G | C | 0.854 | 0.146 | MUC5B H2 | 1.0 |

\normalsize

Freq = frequency of the counted (effect) allele in the IPFJES controls (n=401), which approximate a population reference; MAF = minor-allele frequency. r² (final column) = haplogroup-tagging r² from Plender et al. Supplementary Table S4, except rs35779873, whose r²=0.805 is from their Supplementary Table S6 (it is absent from Table S4, being below their stated r²>0.85 discovery threshold, but Plender designate it a *MUC5AC* H3 tSNP in Table S6 and in their PheWAS). The *MUC5B* discriminator r² is from the Plender main text; rs35705950 is the promoter variant, not a haplogroup tag. Strand-ambiguous (palindromic) C/G or A/T tSNPs had their orientation resolved by allele frequency (Methods), unambiguous where the MAF is away from 0.5, notably the *MUC5B* C/G discriminator (MAF 0.146) and the A/T H3 concordance tag rs28542750 (MAF 0.197, matching the ~21% European H3 frequency). Orientation is in any case fixed by the reference genome, since both the TOPMed output and Plender's allele designations are reported on the GRCh38 forward strand. The H1 concordance tag (rs769768817; r²=0.899, MAF 0.43) is used only as a cross-check; in the primary analysis H1 is the reference baseline and is not entered as a dosage term. Computed by `code/06_tsnp_table.py`.

---

**Supplementary Table S3. *MUC5AC* H2 association across tagging measures**

\scriptsize

| H2 measure | Plender r² | Ctrl freq | r (compl.) | Unadjusted OR (95% CI) | +rs35705950 OR (95% CI) | p |
|--------------|:---:|:---:|:----:|----------|----------|:-:|
| rs879201833 (primary) | 0.888 | 0.237 | 0.973 | 0.56 (0.42–0.75) | 0.66 (0.49–0.89) | 0.006 |
| rs28519516 | 0.862 | 0.241 | 0.963 | 0.54 (0.40–0.72) | 0.64 (0.47–0.86) | 0.003 |
| rs28558973 | 0.862 | 0.241 | 0.963 | 0.54 (0.40–0.72) | 0.64 (0.47–0.86) | 0.003 |
| rs28368633 | 0.862 | 0.241 | 0.963 | 0.54 (0.40–0.72) | 0.64 (0.47–0.86) | 0.003 |
| Diplotype complement (2 − H1 − H3) | n/a | 0.237 | 1.000 | 0.62 (0.47–0.81) | 0.71 (0.53–0.95) | 0.020 |

\normalsize

Ctrl freq = control-allele frequency; r (compl.) = correlation with the diplotype complement. The primary H2 measure is the direct tag rs879201833 (r²=0.888). The three r²=0.862 panel tags are in complete mutual linkage disequilibrium (r=1.000) and give identical rows; the diplotype complement (2 − H1 − H3) is shown as a cross-check. All measures agree in direction and magnitude (adjusted OR 0.64–0.71), so the H2 residual is not an artefact of how H2 is assigned; the tags and the complement have comparable tagging quality (r²≈0.86–0.90), and the small differences between them should not be over-read. All models are the *MUC5AC* joint model (H2 and H3 entered together, H1 reference), adjusted for age, ever-smoking, PC1–5 and recruitment centre (N=798). Computed by `code/11_h2_direct_and_lrt.py`.

---

**Supplementary Table S4. Gene×environment interaction analyses (exploratory)**

\scriptsize

| Genetic factor | Exposure | Mult. OR (95% CI) | p | RERI (95% CI) |
|-------------|---------|-----------|:-:|--------------|
| MUC5AC H3 carrier | Asbestos (Peto) | 0.64 (0.34–1.21) | 0.17 | −0.69 (−2.04 to +0.65) |
| MUC5AC H3 carrier | Ever-smoked | 1.19 (0.60–2.34) | 0.62 | +0.71 (−0.41 to +1.84) |
| MUC5AC H2 non-carrier | Asbestos (Peto) | 0.96 (0.49–1.89) | 0.91 | +0.09 (−0.96 to +1.13) |
| MUC5AC H2 non-carrier | Ever-smoked | 1.44 (0.71–2.93) | 0.31 | +1.03 (+0.10 to +1.97) |
| MUC5B H1 non-carrier | Asbestos (Peto) | 1.45 (0.71–2.96) | 0.30 | +0.38 (−0.22 to +0.97) |
| MUC5B H1 non-carrier | Ever-smoked | 1.10 (0.49–2.51) | 0.81 | +0.32 (−0.55 to +1.20) |
| rs35705950 T carrier | Asbestos (Peto) | 0.85 (0.44–1.66) | 0.64 | −0.04 (−2.70 to +2.63) |
| rs35705950 T carrier | Ever-smoked | 1.73 (0.84–3.58) | 0.14 | +2.78 (+0.37 to +5.19) |

\normalsize

Each genetic factor and each exposure is entered as a binary indicator oriented to its risk-increasing direction, as additive-scale interaction (the relative excess risk due to interaction, RERI) requires (Knol and VanderWeele, 2012); for the protective MUC5AC H2 and MUC5B H1 haplogroups the risk-increasing indicator is non-carriage, so the RERI is not distorted by entering a protective factor against a harmful one. RERI confidence intervals are by the delta method. Models adjust for age, ever-smoking, PC1–5 and recruitment centre (N=798); the smoking models omit the duplicate ever-smoking covariate. Post-hoc power figures are not reported; the confidence intervals convey the precision. Every multiplicative-scale interaction is null. On the additive scale, the only RERIs whose confidence interval excludes zero are the two ever-smoking rows for the risk-increasing background, rs35705950 × ever-smoking (+2.78) and non-carriage of MUC5AC H2 × ever-smoking (+1.03); these are two nominal results among many underpowered contrasts and should not be over-interpreted. Computed by `code/16_ancestry_and_promoter_gxe.py`.

---

**Supplementary Table S5. Imputed *MUC5AC* exon-31 repeat length versus IPF**

\scriptsize

| Length measure | Imp. r² | r (H3) | Unadj. OR (95% CI) | +rs35705950 OR (95% CI) |
|--------------|:--:|:--:|--------------|-------------|
| 751 bp exon-31 insertion (chr11:1,188,533) | 0.98 | r=0.98 | 1.67 (1.29–2.15); p=8×10⁻⁵ | 1.13 (0.85–1.49); p=0.40 |
| MUC5AC H3 clade tag (for comparison) | n/a | n/a | 1.69 (1.32–2.17); p=4×10⁻⁵ | 1.17 (0.89–1.54); p=0.27 |
| Composite exon-31 length, per +1 kb | n/a | r=0.92 | 1.46 (1.12–1.90); p=0.006 | 0.99 (0.73–1.32); p=0.92 |
| Composite exon-31 length, per +1 repeat (24 bp) | n/a | r=0.92 | 1.009 (1.003–1.016); p=0.006 | 1.000 (0.993–1.007); p=0.92 |

\normalsize

Imp. r² = imputation r²; r (H3) = correlation with imputed H3 count. Repeat lengths were imputed into the 798 analytical participants from the Bai et al. (2026) long-read VNTR reference panel using the ImputeSV server (BCFtools QC → SHAPEIT5 v5.0.1 phasing → Minimac4 v4.1.6); the 751 bp allele is the exon-31 length variant that tracks the H3 clade in long-read assemblies. Imputation r² is the server's estimate. The composite is total inserted exon-31 sequence per participant relative to GRCh38, summed across imputed alleles (r²≥0.30); results are unchanged if built from the two common well-imputed alleles alone (r=0.96 between definitions). Adding H3 to a model containing the imputed length, or vice versa, does not improve fit (likelihood-ratio p=0.24 and 0.76): the imputed length and the clade tag are the same information. The imputed length therefore behaves as the H3 tag does, associated unadjusted and abolished by rs35705950, consistent with the paper's argument that an SNP-imputed length inherits the tags' linkage disequilibrium and cannot isolate the repeat. Models adjust for age, ever-smoking, PC1–5 and recruitment centre. Computed by `code/12_imputed_vntr_length.py`; submission and provenance in `code/REPRODUCE.md` (step 6).

---

**Supplementary Table S6. Agnostic conditional scan of imputed structural variants across the 11p15.5 mucin cluster**

| | |
|---|---|
| Variants tested (MUC6–MUC5B, chr11:1.01–1.263 Mb) | 109 common (MAF≥0.01), well-imputed (r²≥0.8) VNTR + SV alleles |
| By gene | MUC2 73, MUC5AC 15, MUC5B 8, intergenic/flanking 13, MUC6 0 (imputes poorly) |
| Effective independent tests (Li & Ji, 2005) | M_eff = 54; cluster-wide threshold 0.05/M_eff = 9.3×10⁻⁴ |
| Reaching threshold **marginally** | 22 (20 new: mostly MUC2; 2 the already-reported MUC5AC exon-31) |
| Strongest marginal | MUC2 +1479 bp @1,096,722: OR 2.76, p=3.1×10⁻¹² (r=0.53 with rs35705950 T) |
| Surviving adjustment for rs35705950 | **0**; smallest adjusted p=0.003, expected P(min<0.003) over M_eff=54 = 0.13 |

Every common, well-imputed VNTR and structural-variant allele in the mucin cluster was tested against IPF by logistic regression (adjusted for age, ever-smoking, PC1–5 and recruitment centre), marginally and with rs35705950 T-allele count added; loci present in both the VNTR and SV panels were de-duplicated and kept at their higher imputation r². The strong marginal associations, concentrated in the *MUC2* VNTR domain that no haplogroup analysis reaches, all reflect linkage disequilibrium with rs35705950 (top hits r=0.24–0.53 with the T allele) and are markedly attenuated by adjustment for it (the strongest, a *MUC2* allele, from OR 2.76 to 1.64); none survives the cluster-wide threshold, and the smallest adjusted p is within chance expectation across the correlated tests. Survival of adjustment would in any case not establish a causal repeat effect, since every variant lies in the same linkage-disequilibrium block (Discussion). Within this sample, then, no imputed mucin structural variant reached significance independently of the promoter variant. (The *MUC5AC* H2 haplogroup tags are SNPs, not structural variants, and so are not among the 109 alleles scanned here; the smallest adjusted p in this table, 0.003, is a low-frequency *MUC2* allele, whereas the residual H2 association discussed in the main text is assessed against a different, less stringent multiplicity family, the four primary haplogroup coefficients, not this 54-test cluster-wide threshold.) This is a power-limited negative, not a claim that no such variant exists: in a far larger IPF meta-analysis, GCTA-COJO resolves conditionally-independent signals annotated by nearest gene to *MUC2* (three) and *MUC6* (one), as well as several in *MUC5B* (Kousathanas et al., 2025), which 397 cases are underpowered to detect (*MUC2* and *MUC6* are canonically gastrointestinal mucins, so these labels are positional and need not imply an airway role for those proteins). The locus is therefore genuinely multi-signal; what the haplogroup analysis shows is more specifically that the *MUC5AC* and *MUC5B* *haplogroups* are not separable from rs35705950. Computed by `code/13_mucin_cluster_scan.py`; provenance in `code/REPRODUCE.md` (step 6).

---

**Supplementary Table S7. Linkage disequilibrium between the *MUC5AC*/*MUC5B* haplogroups and rs35705950**

| Haplogroup (concordance tSNP) | IPFJES controls (n=401), r (r²) | IPFJES full sample (n=798), r (r²) | 1000G EUR reference, r² |
|-------------------------------|:--:|:--:|:--:|
| MUC5AC H1 (rs769768817) | −0.11 (0.01) | −0.13 (0.02) | <0.05 |
| MUC5AC H2 (rs879201833) | −0.16 (0.03) | −0.24 (0.06) | <0.05 |
| MUC5AC H3 (rs28542750) | +0.30 (0.09) | +0.36 (0.13) | 0.10 |
| MUC5B H1 (chr11:1,244,757) | −0.08 (0.01) | −0.10 (0.01) | <0.05 |

Each haplogroup is represented by its concordance tSNP (Methods); r is the Pearson correlation of the tag allele-dosage count with the rs35705950 T-allele count, r² its square. The reference column is the r² between rs35705950 and each tSNP in the 1000 Genomes European panel (phase 3; Ensembl LD API), the panel Plender et al. used for their locus-wide LD analysis; Plender et al. do not themselves tabulate haplogroup–rs35705950 LD. Only *MUC5AC* H3 shows appreciable LD with rs35705950, and the IPFJES-controls value (r²=0.09) agrees closely with the reference (0.10), indicating the LD is representative rather than a feature of this cohort; the full-sample estimates are modestly larger because cases are enriched for risk haplotypes. H1, H2 and *MUC5B* H1 have negligible LD (r²<0.05) in both IPFJES and the reference (the Ensembl API does not return pairs below ~0.05). Computed by `code/15_ld_haplogroup_rs35705950.py`.
