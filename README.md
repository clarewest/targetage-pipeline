# targetage-pipeline

An analysis of age-related disease targets using Open Targets Platform and Open Targets Genetics.

[Open Targets Platform](https://www.targetvalidation.org/) \[<span id="a1">[1](#f1)</span>\] aggregates data from a large range of sources, using "evidence" to make connections between "target" entities (a protein, protein complex, or RNA molecule) and "disease" entities (EFO terms). The evidence linking a target and disease is summarised as a target-disease association, with a score ranging from 0-1 to assist with prioritisation. Targets, diseases, and evidence are comprehensively annotated. 


[Open Targets Genetics](https://genetics.opentargets.org/) \[<span id="a2">[2](#f2)</span>\] aggregates human GWAS and functional genetics data in a variant-centric manner to enable large-scale exploration and prioritisation of potential causal variants and genes. It provides a disease-agnostic Variant to Gene (V2G) mapping and a disease-specific Locus to Gene (L2G) mapping for trait-associated loci, with both using a score ranging from 0-1. It also enables systematic comparison between studies, providing the number of shared independently-associated loci between studies, and performing pairwise colocalisation analysis. 

Open Targets Genetics [expands each lead variant](https://genetics-docs.opentargets.org/our-approach/assigning-traits-to-loci) (the variant with the most significant p-value) to include tag variants, representing a more complete set of potentially causal variants at a trait-associated locus. Where summary statistics are available, this expansion is performed using fine-mapping and credible-set analysis. Where only the lead variant is reported per locus, expansion is performed using Linkage Disequilibium (LD) analysis, with the 1000 Genomes Phase 3 (1KG) haplotype panel as a reference population. 
 

## Data

- `disease_list.csv` - A curated list of age-related diseases (ARDs), and the corresponding EFO-codes. Commented lines are ignored.

##### OT_platform (parquet format)
- `diseases` - annotation information for disease entities
- `targets` - annotation information for target entitites
- `evidences` - all evidences used to make associations between targets and diseases
- `evidences/succeeded/sourceId\=ot_genetics_portal/` - evidence taken from the genetics portal (human genome-wide association data)
knownDrugs - evidence from drugs with a known mechanism of action (MOA) and indication, that links a disease to a target



##### OT_genetics (json format)
- `v2d_coloc` - colocalisation analysis result for each variant+studyId (GWAS-GWAS and GWAS-xQTL)



## Analysis

#### Expansion to specific diseases 
First, get the disease annotations for the ARDs (`ardiseases`), including the EFO label (`diseaseName`), Therapeutic Areas, and description. 

Open Targets (optionally) propagates evidence from specific, lower-level EFO terms, to support associations for more general, higher-level EFO terms. For example, evidence to support associations for `osteoarthritis` will include evidence linked to its child term, `osteoarthritis, knee`.

Get EFO terms and names (`specificDiseaseId` and `specificDiseaseName`) for each ARD and its descendant terms. We will use this to retrieve the evidence supporting the associations, but we will also need this information later to decide whether the evidence is relevant for our purposes.

#### Get target-disease associations and associated target annotations

`overall_associations` is all target-disease associations for our ARDs, considering [all data sources](https://docs.targetvalidation.org/data-sources/data-sources) (including genetic, literature, drugs, and other types of evidence). 

`gen_associations` is all target-disease associations based on genetic evidence. 

`ard_targets` is all the target annotations for targets linked to the ARDs through genetic evidence. 

#### Get Open Target Genetics evidence

`otg_evidences` is all the evidence from Open Targets Genetics (i.e. human genome-wide association evidence)

`ard_otg_evidences` is all the Open Targets Genetics evidence for our ARDs. Each "evidence" is a genome-wide significant trait-associated locus mapped to a target (gene) with an L2G score of at least 0.05. Note that a single trait-associated locus may be mapped to more than one target, so an individual GWAS associations may be represented more than once as evidence for different targets. 

`ard_studies` is all the unique GWAS studies (`studyId`) that we are interested in for our ARDs, the trait studied (`trait_reported`), the ARD (`diseaseName`) and specific disease (`specificDiseaseName`) to which this trait was mapped, and whether the study has summary statistics available (`has_sum_stats`).

#### Establishing independent genetic signals

We are interested in the amount of genetic overlap between different ARDs. 

In order to determine how many genetic signals are shared between ARDs, we have to establish the number of independent trait-associated signals for each ARD (i.e. the overlap between GWAS studies for an ARD).

As we don't know the causal variant for each trait-associated locus, we have to use approximations to decide whether the associations from two different GWAS studies represent the same signal. To do this, we combine two approaches from Open Targets Genetics:

1. Colocalisation

Where study 1 and study 2 both have summary statistics available, they can be compared using colocalisation analysis. For each loci, this method integrates over evidence from all variants in each study to evaluate which of these four hypothesis is most likely: no association with either trait (H0), association only with trait 1 (H1), association only with trait 2 (H2), association with both traits via two independent SNPs (H3), or association with both traits through a shared causal SNP (H4). We use a cut-off of H4>0.8 to define the association for study 1 and study 2 as sharing the same causal variant. Note that this methodology is used to compare two disease traits, or to compare a disease study with a molecular trait (e.g. eQTL or pQTL). 

2. Tag variant overlap

Where only the lead variants (not full summary statistics) are available for one or both studies, colocalisation analysis cannot be performed, and tag variant overlap is used instead. A trait-associated lead variant, SNP1, is considered to be part of the same signal as another lead variant, SNP2, if they are within 5MB of each other and any of the LD-defined tag variants are shared. 



\[<b id="f1">1</b>\] Ochoa et al., Open Targets Platform: supporting systematic drug–target identification and prioritisation, Nucleic Acids Research (2021) https://doi.org/10.1093/nar/gkaa1027 [↩](#a1)

\[<b id="f2">2</b>\] Ghoussaini et al., Open Targets Genetics: systematic identification of trait-associated genes using large-scale genetics and functional genomics, Nucleic Acids Research (2021) https://doi.org/10.1093/nar/gkaa840. [↩](#a2)
