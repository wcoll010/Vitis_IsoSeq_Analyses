### Identify genes under differential isoform usage in the V. riparia X V. rupestris hybrids (3309, SW, 101-14)
### Nine pairwise contrasts: within genotype (Infested vs Control); between genotypes within Control; between genotypes within Infested
### Previous anaylses (Differential_isoform_usage.txt) compared infested vs control for each genotype separately 


# ============================================================
# Inputs pre-filtered: SQANTI filter_result == "Isoform"; FPM > 0.1 in >=2 of all samples; genes with >=2 retained isoforms
# DIU criteria: FDR < 0.05; |dIF| >= 0.10
# ============================================================



#!/usr/bin/env Rscript


# ============================================================
# 1. Read contrast argument
# ============================================================

args <- commandArgs(
    trailingOnly = TRUE
)

if (length(args) != 1) {

    stop(
        "Usage:\n",
        "Rscript ISA_crossgenotype_array.R <contrast_name>"
    )
}

contrast_name <- args[1]


# ============================================================
# 2. R library
# ============================================================

my_lib <- "/rhome/wcoll010/bigdata_groenlab/R_libs"

.libPaths(
    c(
        my_lib,
        .libPaths()
    )
)

library(IsoformSwitchAnalyzeR)


# ============================================================
# 3. Paths
# ============================================================

BASE <- "/rhome/wcoll010/shared/pacbio/IsoSeq_analyses_2"

ISA_BASE <- file.path(
    BASE,
    "IsoformSwitchAnalyzeR_cross_genotype"
)

INPUT_RDATA <- file.path(
    ISA_BASE,
    "Vrip_X_Vrup_ISA_common_input_objects.RData"
)

COMMON_GTF <- file.path(
    ISA_BASE,
    "Vrip_X_Vrup_common_DTU_filtered.gtf"
)

OUT_BASE <- file.path(
    ISA_BASE,
    "results"
)

dir.create(
    OUT_BASE,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 4. Define pairwise contrasts
#
# IMPORTANT:
# Use R-safe condition labels.
#
# Condition 1 = reference
# Condition 2 = comparison
#
# Positive dIF =
# increased relative isoform usage in condition 2.
# ============================================================

comparisons <- list(

    "101-14_Inf_vs_Con" = c(
        "G101_Control",
        "G101_Infested"
    ),

    "3309_Inf_vs_Con" = c(
        "G3309_Control",
        "G3309_Infested"
    ),

    "SW_Inf_vs_Con" = c(
        "SW_Control",
        "SW_Infested"
    ),

    "3309C_vs_101-14C" = c(
        "G101_Control",
        "G3309_Control"
    ),

    "SWC_vs_101-14C" = c(
        "G101_Control",
        "SW_Control"
    ),

    "3309C_vs_SWC" = c(
        "SW_Control",
        "G3309_Control"
    ),

    "3309I_vs_101-14I" = c(
        "G101_Infested",
        "G3309_Infested"
    ),

    "SWI_vs_101-14I" = c(
        "G101_Infested",
        "SW_Infested"
    ),

    "3309I_vs_SWI" = c(
        "SW_Infested",
        "G3309_Infested"
    )
)


if (!contrast_name %in% names(comparisons)) {

    stop(
        "Unknown contrast: ",
        contrast_name
    )
}


condition_1 <- comparisons[[contrast_name]][1]
condition_2 <- comparisons[[contrast_name]][2]


# ============================================================
# 5. Check required files
# ============================================================

if (!file.exists(INPUT_RDATA)) {

    stop(
        "Input RData not found:\n",
        INPUT_RDATA
    )
}


if (!file.exists(COMMON_GTF)) {

    stop(
        "Common filtered GTF not found:\n",
        COMMON_GTF
    )
}


# ============================================================
# 6. Load common globally filtered inputs
#
# Expected objects:
#
#   isoformCountMatrix
#   designMatrix
#   isoformAnnotation
# ============================================================

load(
    INPUT_RDATA
)


required_objects <- c(
    "isoformCountMatrix",
    "designMatrix",
    "isoformAnnotation"
)


missing_objects <- required_objects[
    !vapply(
        required_objects,
        exists,
        logical(1)
    )
]


if (length(missing_objects) > 0) {

    stop(
        "Missing objects in input RData:\n",
        paste(
            missing_objects,
            collapse = "\n"
        )
    )
}


# ============================================================
# 7. Recode condition labels to R-safe names
#
# Original biological labels are preserved by the contrast
# names and sample IDs. Only the model-facing condition names
# are changed.
# ============================================================

designMatrix$condition <- as.character(
    designMatrix$condition
)


condition_map <- c(

    "101-14_Control" =
        "G101_Control",

    "101-14_Infested" =
        "G101_Infested",

    "3309_Control" =
        "G3309_Control",

    "3309_Infested" =
        "G3309_Infested",

    "SW_Control" =
        "SW_Control",

    "SW_Infested" =
        "SW_Infested"
)


unknown_conditions <- setdiff(
    unique(designMatrix$condition),
    names(condition_map)
)


if (length(unknown_conditions) > 0) {

    stop(
        "Unexpected condition labels in designMatrix:\n",
        paste(
            unknown_conditions,
            collapse = "\n"
        )
    )
}


designMatrix$condition <- unname(
    condition_map[
        designMatrix$condition
    ]
)


cat("\nRecoded condition labels:\n")

print(
    table(
        designMatrix$condition
    )
)


# ============================================================
# 8. Basic input checks
# ============================================================

if (!"isoform_id" %in% colnames(isoformCountMatrix)) {

    stop(
        "isoformCountMatrix does not contain isoform_id."
    )
}


if (
    !all(
        c(
            "sampleID",
            "condition"
        ) %in%
        colnames(designMatrix)
    )
) {

    stop(
        "designMatrix must contain sampleID and condition."
    )
}


if (
    !all(
        c(
            "isoform_id",
            "gene_id"
        ) %in%
        colnames(isoformAnnotation)
    )
) {

    stop(
        "isoformAnnotation must contain isoform_id and gene_id."
    )
}


if (anyDuplicated(isoformCountMatrix$isoform_id)) {

    stop(
        "Duplicated isoform IDs in isoformCountMatrix."
    )
}


if (anyDuplicated(isoformAnnotation$isoform_id)) {

    stop(
        "Duplicated isoform IDs in isoformAnnotation."
    )
}


# ============================================================
# 9. Output directory
# ============================================================

OUT_DIR <- file.path(
    OUT_BASE,
    contrast_name
)

dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


cat("\n")
cat("====================================================\n")
cat("CROSS-GENOTYPE DIFFERENTIAL ISOFORM USAGE\n")
cat("====================================================\n")
cat("Contrast:    ", contrast_name, "\n")
cat("Condition 1: ", condition_1, "\n")
cat("Condition 2: ", condition_2, "\n")
cat("====================================================\n")


# ============================================================
# 10. Select samples for this contrast
# ============================================================

keep_samples <- designMatrix$condition %in% c(
    condition_1,
    condition_2
)


design_sub <- designMatrix[
    keep_samples,
    ,
    drop = FALSE
]


if (nrow(design_sub) == 0) {

    stop(
        "No samples found for contrast: ",
        contrast_name
    )
}


design_sub$condition <- factor(
    design_sub$condition,
    levels = c(
        condition_1,
        condition_2
    )
)


cat("\nDesign:\n")

print(
    design_sub
)


cat("\nReplicates per condition:\n")

print(
    table(
        design_sub$condition
    )
)


if (
    any(
        table(
            design_sub$condition
        ) < 2
    )
) {

    stop(
        "Each condition requires >=2 biological replicates."
    )
}


# Extra sanity check:
# exactly two modeled conditions must be present.

if (nlevels(design_sub$condition) != 2) {

    stop(
        "design_sub$condition does not contain exactly two levels."
    )
}


if (length(unique(design_sub$condition)) != 2) {

    stop(
        "Only one observed condition remains after subsetting."
    )
}


sample_names <- design_sub$sampleID


missing_samples <- setdiff(
    sample_names,
    colnames(
        isoformCountMatrix
    )
)


if (length(missing_samples) > 0) {

    stop(
        "Samples missing from count matrix:\n",
        paste(
            missing_samples,
            collapse = "\n"
        )
    )
}


# ============================================================
# 11. Subset samples only
#
# IMPORTANT:
#
# No transcript-level expression filtering is performed here.
#
# importRdata() will automatically discard transcripts that
# are not expressed in any selected sample.
# ============================================================

counts_sub <- isoformCountMatrix[
    ,
    c(
        "isoform_id",
        sample_names
    ),
    drop = FALSE
]


cat(
    "\nCommon globally filtered isoforms:",
    nrow(
        counts_sub
    ),
    "\n"
)


cat(
    "Common globally filtered genes:",
    length(
        unique(
            isoformAnnotation$gene_id
        )
    ),
    "\n"
)


# ============================================================
# 12. Save design used for this contrast
# ============================================================

write.table(
    design_sub,
    file = file.path(
        OUT_DIR,
        paste0(
            contrast_name,
            "_design.tsv"
        )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# 13. Import into IsoformSwitchAnalyzeR
#
# For this diagnostic run:
#
#   detectUnwantedEffects = FALSE
#
# We are testing whether R-safe condition names resolve the
# 3309 vs 101-14 model-matrix failure.
#
# No preFilter() is used.
# ============================================================

cat(
    "\nImporting into IsoformSwitchAnalyzeR...\n"
)


aSwitchList <- importRdata(

    isoformCountMatrix =
        counts_sub,

    designMatrix =
        design_sub,

    isoformExonAnnoation =
        COMMON_GTF,

    detectUnwantedEffects =
        FALSE,

    fixStringTieAnnotationProblem =
        FALSE,

    ignoreAfterBar =
        FALSE,

    showProgress =
        TRUE
)


# ============================================================
# 14. Save imported object
# ============================================================

saveRDS(
    aSwitchList,
    file = file.path(
        OUT_DIR,
        paste0(
            contrast_name,
            "_DIU_imported.rds"
        )
    )
)


cat(
    "\nImported isoforms:",
    nrow(
        aSwitchList$isoformFeatures
    ),
    "\n"
)


cat(
    "Imported genes:",
    length(
        unique(
            aSwitchList$isoformFeatures$gene_id
        )
    ),
    "\n"
)


# ============================================================
# 15. Run DEXSeq differential isoform usage
#
# FDR < 0.05
# |dIF| >= 0.10
#
# Keep complete tested object.
# ============================================================

cat(
    "\nRunning DEXSeq differential isoform usage...\n"
)


aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(

    switchAnalyzeRlist =
        aSwitchList,

    alpha =
        0.05,

    dIFcutoff =
        0.10,

    reduceToSwitchingGenes =
        FALSE,

    reduceFurtherToGenesWithConsequencePotential =
        FALSE,

    onlySigIsoforms =
        FALSE,

    keepIsoformInAllConditions =
        TRUE,

    showProgress =
        TRUE
)


# ============================================================
# 16. Save complete analyzed object
# ============================================================

saveRDS(
    aSwitchListAnalyzed,
    file = file.path(
        OUT_DIR,
        paste0(
            contrast_name,
            "_DIU_results_final.rds"
        )
    )
)


# ============================================================
# 17. Overall switch summary
# ============================================================

switch_summary <- extractSwitchSummary(
    aSwitchListAnalyzed
)


write.csv(
    switch_summary,
    file = file.path(
        OUT_DIR,
        paste0(
            contrast_name,
            "_DIU_summary.csv"
        )
    ),
    row.names = FALSE
)


# ============================================================
# 18. Significant switching genes
# ============================================================

switching_genes <- extractTopSwitches(

    switchAnalyzeRlist =
        aSwitchListAnalyzed,

    filterForConsequences =
        FALSE,

    extractGenes =
        TRUE,

    alpha =
        0.05,

    dIFcutoff =
        0.10,

    n =
        Inf,

    sortByQvals =
        TRUE
)


write.csv(
    switching_genes,
    file = file.path(
        OUT_DIR,
        paste0(
            contrast_name,
            "_significant_switching_genes.csv"
        )
    ),
    row.names = FALSE
)


# ============================================================
# 19. Significant switching isoforms
# ============================================================

switching_isoforms <- extractTopSwitches(

    switchAnalyzeRlist =
        aSwitchListAnalyzed,

    filterForConsequences =
        FALSE,

    extractGenes =
        FALSE,

    alpha =
        0.05,

    dIFcutoff =
        0.10,

    n =
        Inf,

    sortByQvals =
        TRUE
)


write.csv(
    switching_isoforms,
    file = file.path(
        OUT_DIR,
        paste0(
            contrast_name,
            "_significant_switching_isoforms.csv"
        )
    ),
    row.names = FALSE
)


# ============================================================
# 20. Export full isoform feature table
# ============================================================

write.csv(
    aSwitchListAnalyzed$isoformFeatures,
    file = file.path(
        OUT_DIR,
        paste0(
            contrast_name,
            "_DIU_isoform_features_full.csv"
        )
    ),
    row.names = FALSE
)


# ============================================================
# 21. Export full switch-analysis table if present
# ============================================================

if (
    "isoformSwitchAnalysis" %in%
    names(
        aSwitchListAnalyzed
    )
) {

    write.csv(
        aSwitchListAnalyzed$isoformSwitchAnalysis,
        file = file.path(
            OUT_DIR,
            paste0(
                contrast_name,
                "_DIU_isoformSwitchAnalysis_full.csv"
            )
        ),
        row.names = FALSE
    )
}


# ============================================================
# 22. Compact run summary
# ============================================================

run_summary <- data.frame(

    contrast =
        contrast_name,

    condition_1 =
        condition_1,

    condition_2 =
        condition_2,

    samples_condition_1 =
        sum(
            design_sub$condition ==
                condition_1
        ),

    samples_condition_2 =
        sum(
            design_sub$condition ==
                condition_2
        ),

    common_input_isoforms =
        nrow(
            counts_sub
        ),

    imported_isoforms =
        nrow(
            aSwitchList$isoformFeatures
        ),

    imported_genes =
        length(
            unique(
                aSwitchList$isoformFeatures$gene_id
            )
        ),

    significant_switching_genes =
        nrow(
            switching_genes
        ),

    significant_switching_isoforms =
        nrow(
            switching_isoforms
        )
)


write.csv(
    run_summary,
    file = file.path(
        OUT_DIR,
        paste0(
            contrast_name,
            "_run_summary.csv"
        )
    ),
    row.names = FALSE
)


# ============================================================
# 23. Finish
# ============================================================

cat("\n")
cat("====================================================\n")
cat("DIU COMPLETE\n")
cat("====================================================\n")

print(
    run_summary
)

cat(
    "\nResults directory:",
    OUT_DIR,
    "\n"
)

cat("====================================================\n")
