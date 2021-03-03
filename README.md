# targetage-pipeline

An analysis of age-related disease targets using Open Targets Platform and Open Targets Genetics.

[Open Targets Platform](https://www.targetvalidation.org/) <sup id="a1">[1](#f1)</sup> aggregates data from a large range of sources, using "evidence" to make connections between "target" entities (a protein, protein complex, or RNA molecule) and "disease" entities (EFO terms). The evidence linking a target and disease is summarised as a target-disease association, with a score ranging from 0-1 to assist with prioritisation. Targets, diseases, and evidence are comprehensively annotated. 


[Open Targets Genetics](https://genetics.opentargets.org/) <sup id="a2">[2](#f2)</sup> aggregates human GWAS and functional genetics data in a variant-centric manner to enable large-scale exploration and prioritisation of potential causal variants and genes. It provides a disease-agnostic Variant to Gene (V2G) mapping and a disease-specific Locus to Gene (L2G) mapping for trait-associated loci, with both using a score ranging from 0-1. The L2G score  It also enables systematic comparison between studies, providing the number of shared independently-associated loci between studies, and performing pairwise colocalisation analysis. 
 

** Data
disease_list.csv - A curated list of age-related diseases, and the corresponding EFO-codes. Commented lines are ignored.
OT_platform (parquet format)
diseases - annotation information for disease entities
targets - annotation information for target entitites
evidences - all evidences used to make associations between targets and diseases
evidences/succeeded/sourceId\=ot_genetics_portal/ - evidence taken from the genetics portal (human genome-wide association data)
knownDrugs - evidence from drugs with a known mechanism of action (MOA) and indication, that links a disease to a target



OT_genetics (json format)
v2d_coloc - colocalisation analysis result for each variant+studyId (GWAS-GWAS and GWAS-xQTL)



** Analysis


<b id="f1">1</b> Ochoa et al., Open Targets Platform: supporting systematic drug–target identification and prioritisation, Nucleic Acids Research (2021) https://doi.org/10.1093/nar/gkaa1027 [↩](#a1)
<b id="f2">2</b> Ghoussaini et al., Open Targets Genetics: systematic identification of trait-associated genes using large-scale genetics and functional genomics, Nucleic Acids Research (2021) https://doi.org/10.1093/nar/gkaa840. [↩](#a2)
