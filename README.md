# targetage-pipeline

An analysis of age-related disease targets using Open Targets Platform and Open Targets Genetics.

Open Targets Platform aggregates data from a large range of sources, using "evidence" to make connections between "target" entities (a protein, protein complex, or RNA molecule) and "disease" entities (EFO terms). The evidence linking a target and disease is summarised as a target-disease association, with a score ranging from 0-1 to assist with prioritisation. Targets, diseases, and evidence are comprehensively annotated. 

Open Targets Genetics aggregates human GWAS and functional genetics data in a variant-centric manner to enable large-scale exploration and prioritisation of potential causal variants and genes. It provides a disease-agnostic Variant to Gene (V2G) mapping and a disease-specific Locus to Gene (L2G) mapping for trait-associated loci, with both using a score ranging from 0-1. The L2G score  It also enables systematic comparison between studies, providing the number of shared independently-associated loci between studies, and performing pairwise colocalisation analysis. 
 

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


Currently, the way this is done is the following: To define overlap for a given lead variant (e.g. SNP1), the LD-defined tag variants of SNP1 are cross-referenced to the tag variants of all lead variants within 5MB of SNP1. In any case where a tag variant of SNP1 is shared with another lead variant, that lead is considered part of the same signal as SNP1. (edited) 
