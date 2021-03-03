#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 20:07:01 2021

@author: clarewest
"""

import pyspark
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, lit, explode, concat_ws

## Spark
sc = pyspark.SparkContext()
spark = SparkSession.builder.getOrCreate()

## Data paths
datapath = "data/"
ot_platform = datapath+"OT_platform/21.02/"
ot_genetics = datapath+"data/OT_genetics/ftp.ebi.ac.uk/pub/databases/opentargets/genetics/20022712/"

## Read Open Targets data 
diseases = (spark.read.parquet(ot_platform+"ETL_parquet/diseases/", header=True)
            .withColumnRenamed("id","diseaseId")
            .withColumnRenamed("name","diseaseName")
            )
targets = (spark.read.parquet(ot_platform+"ETL_parquet/targets/")
           .withColumnRenamed("id","targetId")
           .withColumnRenamed("approvedSymbol", "targetSymbol")
           .withColumnRenamed("approvedName","targetName")
           )
evidences = spark.read.parquet(ot_platform+"ETL_parquet/evidences/succeeded")
knowndrugs = spark.read.parquet(ot_platform+"ETL_parquet/knownDrugs")
otg_evidences = spark.read.parquet(ot_platform+"ETL_parquet/evidences/succeeded/sourceId\=ot_genetics_portal/")

## Genetic data
coloc = spark.read.json(ot_genetics+"v2d_coloc/")
#variants = spark.read.json("lut/variant-index/")
studies = spark.read.json(ot_genetics+"lut/study-index/")

## Age-related diseases (ARDs)
ards = spark.read.csv(datapath+"disease_list.csv")
ards = ards.toDF(*["diseaseId", "morbidity"]).filter(~col("diseaseId").contains("#"))
ardiseases = (ards.join(diseases, "diseaseId", "left")
              .select("morbidity","diseaseId","children", "description", "diseaseName", "therapeuticAreas", "descendants")
              )

# Get EFO codes etc for descendant diseases
descendant_ardiseases = (
    ardiseases
    .select("diseaseId", "diseaseName", explode("descendants").alias("specificDiseaseId"))
    .join(diseases.withColumnRenamed("diseaseId", "specificDiseaseId").withColumnRenamed("diseaseName", "specificDiseaseName"), ["specificDiseaseId"])   
    .join(ards, "diseaseId")
    .select("morbidity", "diseaseId", "diseaseName", "specificDiseaseId", "specificDiseaseName", "therapeuticAreas", "description")
)

parent_ardiseases = (ardiseases
                     .select("morbidity",
                             "diseaseId", 
                             "diseaseName", 
                             col("diseaseId").alias("specificDiseaseId"), 
                             col("diseaseName").alias("specificDiseaseName"),
                             "therapeuticAreas",
                             "description")
                     )

all_ardiseases = parent_ardiseases.union(descendant_ardiseases)


## Association data 
overall_associations = spark.read.parquet(ot_platform+"ETL_parquet/associations/indirect/byOverall/")
associations = spark.read.parquet(ot_platform+"ETL_parquet/associations/indirect/byDatatype/")
gen_associations = associations.filter(col("datatypeId")=="genetic_association")


## Targets with genetic associations 
ard_associations = (ards.join(gen_associations,"diseaseId", "inner"))
ard_associations.groupBy("morbidity").count().show()
# are direct associations included in indirect associations?

## Get annotations for all targets implicated
ard_targets = ard_associations.join(targets, ["targetId", "targetSymbol", "targetName"])


otg_evidence_cols = ["morbidity",
                     "specificDiseaseId", "specificDiseaseName", 
                     "diseaseId", "diseaseName", 
                     "targetId", 
                     "variantId",
                     "diseaseFromSourceId",
                     "variantFunctionalConsequenceId",
                     "publicationYear",
                     "resourceScore",
                     "diseaseFromSource",
                     "oddsRatio",
                     "confidenceIntervalUpper", "confidenceIntervalLower",
                     "variantRsId",
                     "studyId",
                     "studyCases", "studySampleSize",
                     "score"
                     ]

ard_otg_evidences = (all_ardiseases.drop("therapeuticAreas", "description")
                    .join(otg_evidences
                             .withColumnRenamed("diseaseId", "specificDiseaseId")
                             .withColumnRenamed("diseaseLabel", "specificDiseaseName"), 
                             ["specificDiseaseId", "specificDiseaseName"], "inner")
                    .select(otg_evidence_cols)
                    )

ards_studies = (ard_otg_evidences
                .select("morbidity", "studyId", "diseaseName")
                .distinct()
                .join(studies.select(col("study_id").alias("studyId"), "has_sumstats","trait_reported"), "studyId")
                )
                
# variant ID = chromosome_position_reference_alternative
coloc_studies = (coloc
                 .join(ards_studies, "left_study")
                 .withColumn("variantId", concat_ws("_", col("left_chrom"), col("left_pos"), col("left_ref"), col("left_alt")))
                 .withColumn("right_variantId", concat_ws("_", col("right_chrom"), col("right_pos"), col("right_ref"), col("right_alt")))
                 )


