### Merge significant cross-genotype DIU results with top-50 WGCNA hub genes within the target WGCNA modules midnightblue, white, lightcyan


#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)



# ============================================================
# 1. Paths
# ============================================================

BASE <- "/rhome/wcoll010/shared/pacbio/IsoSeq_analyses_2"


HUB_FILE <- file.path(
    BASE,
    "WGCNA_Vrip_X_Vrup",
    "signed_network",
    "target_module_hubs",
    "Top_50_Hub_Genes_All_Target_Modules.csv"
)


DIU_DIR <- file.path(
    BASE,
    "IsoformSwitchAnalyzeR_cross_genotype",
    "results"
)


SQANTI_FILE <- file.path(
    BASE,
    "sqanti_results",
    "Vrip_X_Vrup_sqanti",
    "filter_results",
    "Vrip_X_Vrup_filtered_RulesFilter_result_classification.txt"
)


OUT_DIR <- file.path(
    BASE,
    "WGCNA_Vrip_X_Vrup",
    "signed_network",
    "DIU_overlap_cross_genotype",
    "top50_hubs_corrected"
)


dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 2. Contrasts
# ============================================================

contrasts <- c(
    "101-14_Inf_vs_Con",
    "3309_Inf_vs_Con",
    "SW_Inf_vs_Con",
    "3309C_vs_101-14C",
    "SWC_vs_101-14C",
    "3309C_vs_SWC",
    "3309I_vs_101-14I",
    "SWI_vs_101-14I",
    "3309I_vs_SWI"
)


contrast_labels <- c(

    "101-14_Inf_vs_Con" =
        "101-14 Infested vs Control",

    "3309_Inf_vs_Con" =
        "3309 Infested vs Control",

    "SW_Inf_vs_Con" =
        "SW Infested vs Control",

    "3309C_vs_101-14C" =
        "3309 Control vs 101-14 Control",

    "SWC_vs_101-14C" =
        "SW Control vs 101-14 Control",

    "3309C_vs_SWC" =
        "3309 Control vs SW Control",

    "3309I_vs_101-14I" =
        "3309 Infested vs 101-14 Infested",

    "SWI_vs_101-14I" =
        "SW Infested vs 101-14 Infested",

    "3309I_vs_SWI" =
        "3309 Infested vs SW Infested"
)


# ============================================================
# 3. Check required main files
# ============================================================

if (!file.exists(HUB_FILE)) {

    stop(
        "Hub file not found:\n",
        HUB_FILE
    )
}


if (!file.exists(SQANTI_FILE)) {

    stop(
        "SQANTI file not found:\n",
        SQANTI_FILE
    )
}


# ============================================================
# 4. Read Top50 hub table
# ============================================================

hubs <- read.csv(
    HUB_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


if (!"Gene_ID" %in% colnames(hubs)) {

    stop(
        "Gene_ID column missing from Top50 hub file."
    )
}


hubs$Gene_ID <- as.character(
    hubs$Gene_ID
)


cat("\n")
cat("====================================================\n")
cat("TOP50 HUB TABLE\n")
cat("====================================================\n")

cat(
    "Rows:",
    nrow(hubs),
    "\n"
)

cat(
    "Unique Gene_IDs:",
    n_distinct(hubs$Gene_ID),
    "\n"
)


if ("Module_Color" %in% colnames(hubs)) {

    cat("\nGenes per module:\n")

    print(
        table(
            hubs$Module_Color
        )
    )
}


# ------------------------------------------------------------
# Check whether a Gene_ID occurs multiple times in hub table
# ------------------------------------------------------------

duplicate_hubs <- hubs %>%
    count(
        Gene_ID,
        name = "N_Rows"
    ) %>%
    filter(
        N_Rows > 1
    )


write.csv(
    duplicate_hubs,
    file.path(
        OUT_DIR,
        "QC_Duplicate_Gene_IDs_In_Top50_File.csv"
    ),
    row.names = FALSE
)


cat(
    "\nDuplicated Top50 Gene_IDs:",
    nrow(duplicate_hubs),
    "\n"
)


# ============================================================
# 5. Optional module-membership column
# Re-computes the MM value (module membership) for each gene in the modules
# ============================================================

required_mm_cols <- c(
    "Module_Color",
    "MM.lightcyan",
    "MM.white",
    "MM.midnightblue"
)


if (
    all(
        required_mm_cols %in%
        colnames(hubs)
    )
) {

    hubs <- hubs %>%
        mutate(

            Module_Membership =
                case_when(

                    Module_Color == "lightcyan" ~
                        MM.lightcyan,

                    Module_Color == "white" ~
                        MM.white,

                    Module_Color == "midnightblue" ~
                        MM.midnightblue,

                    TRUE ~
                        NA_real_
                )
        )


    hubs <- hubs %>%
        arrange(
            Module_Color,
            desc(
                Module_Membership
            )
        ) %>%
        group_by(
            Module_Color
        ) %>%
        mutate(
            Hub_Rank_By_MM =
                row_number()
        ) %>%
        ungroup()


    cat(
        "Module_Membership and Hub_Rank_By_MM created.\n"
    )

} else {

    cat(
        "MM columns not all present; ",
        "existing hub columns will be retained unchanged.\n"
    )
}


# ============================================================
# 6. Read SQANTI classification
# ============================================================

sqanti <- read.delim(
    SQANTI_FILE,
    header = TRUE,
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)


required_sqanti <- c(
    "isoform",
    "associated_gene",
    "associated_transcript",
    "structural_category"
)


missing_sqanti <- setdiff(
    required_sqanti,
    colnames(sqanti)
)


if (length(missing_sqanti) > 0) {

    stop(
        "Missing SQANTI columns:\n",
        paste(
            missing_sqanti,
            collapse = "\n"
        )
    )
}


cat("\n")
cat("====================================================\n")
cat("SQANTI CLASSIFICATION\n")
cat("====================================================\n")

cat(
    "Rows:",
    nrow(sqanti),
    "\n"
)

cat(
    "Unique isoform IDs:",
    n_distinct(sqanti$isoform),
    "\n"
)


# ============================================================
# 7. Check exact SQANTI isoform mapping
# ============================================================

sqanti_isoform_gene_pairs <- sqanti %>%

    transmute(

        Isoform_ID =
            as.character(
                isoform
            ),

        Gene_ID =
            as.character(
                associated_gene
            )

    ) %>%

    filter(
        !is.na(Isoform_ID),
        Isoform_ID != ""
    ) %>%

    distinct()


sqanti_isoform_mapping_qc <- sqanti_isoform_gene_pairs %>%

    group_by(
        Isoform_ID
    ) %>%

    summarise(

        N_Associated_Genes =
            n_distinct(
                Gene_ID[
                    !is.na(Gene_ID) &
                    Gene_ID != ""
                ]
            ),

        Associated_Genes =
            paste(
                sort(
                    unique(
                        Gene_ID[
                            !is.na(Gene_ID) &
                            Gene_ID != ""
                        ]
                    )
                ),
                collapse = " | "
            ),

        .groups = "drop"
    )


ambiguous_isoforms <- sqanti_isoform_mapping_qc %>%
    filter(
        N_Associated_Genes > 1
    )


write.csv(
    ambiguous_isoforms,
    file.path(
        OUT_DIR,
        "QC_SQANTI_Isoforms_With_Multiple_Associated_Genes.csv"
    ),
    row.names = FALSE
)


cat(
    "Exact SQANTI isoforms mapping to >1 gene:",
    nrow(ambiguous_isoforms),
    "\n"
)


# This should be zero; if not zero, stop 

if (nrow(ambiguous_isoforms) > 0) {

    stop(
        "\nERROR:\n",
        nrow(ambiguous_isoforms),
        " exact SQANTI isoforms map to multiple associated genes.\n",
        "Inspect:\n",
        file.path(
            OUT_DIR,
            "QC_SQANTI_Isoforms_With_Multiple_Associated_Genes.csv"
        )
    )
}


# ============================================================
# 8. Build exact isoform -> SQANTI gene map
# ============================================================

sqanti_map <- sqanti %>%

    transmute(

        Isoform_ID =
            as.character(
                isoform
            ),

        Gene_ID =
            as.character(
                associated_gene
            ),

        SQANTI_Associated_Transcript =
            as.character(
                associated_transcript
            ),

        SQANTI_Structural_Category =
            as.character(
                structural_category
            ),

        SQANTI_Filter_Result =
            if (
                "filter_result" %in%
                colnames(sqanti)
            ) {
                as.character(
                    filter_result
                )
            } else {
                NA_character_
            }

    ) %>%

    distinct(
        Isoform_ID,
        .keep_all = TRUE
    )


if (anyDuplicated(sqanti_map$Isoform_ID)) {

    stop(
        "Isoform_ID is not unique in sqanti_map."
    )
}


cat(
    "Final unique SQANTI isoform mappings:",
    nrow(sqanti_map),
    "\n"
)


cat(
    "SQANTI isoforms without associated_gene:",
    sum(
        is.na(sqanti_map$Gene_ID) |
        sqanti_map$Gene_ID == ""
    ),
    "\n"
)


# ============================================================
# 9. Read significant switching genes
# ============================================================

gene_list <- list()


for (contrast in contrasts) {

    gene_file <- file.path(
        DIU_DIR,
        contrast,
        paste0(
            contrast,
            "_significant_switching_genes.csv"
        )
    )


    if (!file.exists(gene_file)) {

        stop(
            "Missing DIU gene file:\n",
            gene_file
        )
    }


    x <- read.csv(
        gene_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )


    required_gene_cols <- c(
        "gene_id",
        "gene_switch_q_value"
    )


    missing_cols <- setdiff(
        required_gene_cols,
        colnames(x)
    )


    if (length(missing_cols) > 0) {

        stop(
            "Missing columns from:\n",
            gene_file,
            "\n",
            paste(
                missing_cols,
                collapse = ", "
            )
        )
    }


    x <- x %>%

        transmute(

            Contrast =
                contrast,

            Contrast_Label =
                unname(
                    contrast_labels[
                        contrast
                    ]
                ),

            PB_Gene_ID =
                as.character(
                    gene_id
                ),

            Gene_Switch_Q_Value =
                as.numeric(
                    gene_switch_q_value
                )
        ) %>%

        distinct(
            Contrast,
            PB_Gene_ID,
            .keep_all = TRUE
        )


    gene_list[[contrast]] <- x
}


significant_genes <- bind_rows(
    gene_list
)


cat("\n")
cat("====================================================\n")
cat("SIGNIFICANT GENE CALLS\n")
cat("====================================================\n")

cat(
    "Total gene-contrast calls:",
    nrow(significant_genes),
    "\n"
)

cat(
    "Unique PB genes:",
    n_distinct(significant_genes$PB_Gene_ID),
    "\n"
)


# ============================================================
# 10. Read significant switching isoforms
# ============================================================

isoform_list <- list()


for (contrast in contrasts) {

    cat("\n")
    cat("====================================================\n")
    cat("READING:", contrast, "\n")
    cat("====================================================\n")


    isoform_file <- file.path(
        DIU_DIR,
        contrast,
        paste0(
            contrast,
            "_significant_switching_isoforms.csv"
        )
    )


    if (!file.exists(isoform_file)) {

        stop(
            "Missing switching isoform file:\n",
            isoform_file
        )
    }


    x <- read.csv(
        isoform_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )


    required_isoform_cols <- c(
        "isoform_id",
        "gene_id",
        "condition_1",
        "condition_2",
        "IF1",
        "IF2",
        "dIF",
        "isoform_switch_q_value"
    )


    missing_cols <- setdiff(
        required_isoform_cols,
        colnames(x)
    )


    if (length(missing_cols) > 0) {

        stop(
            "Missing required columns from:\n",
            isoform_file,
            "\n",
            paste(
                missing_cols,
                collapse = ", "
            )
        )
    }


    x <- x %>%

        mutate(

            Contrast =
                contrast,

            Contrast_Label =
                unname(
                    contrast_labels[
                        contrast
                    ]
                ),

            Isoform_ID =
                as.character(
                    isoform_id
                ),

            PB_Gene_ID =
                as.character(
                    gene_id
                ),

            IF1 =
                as.numeric(
                    IF1
                ),

            IF2 =
                as.numeric(
                    IF2
                ),

            dIF =
                as.numeric(
                    dIF
                ),

            Isoform_Switch_Q_Value =
                as.numeric(
                    isoform_switch_q_value
                )

        )


    cat(
        "Significant switching isoforms:",
        n_distinct(x$Isoform_ID),
        "\n"
    )


    cat(
        "PB genes represented:",
        n_distinct(x$PB_Gene_ID),
        "\n"
    )


    isoform_list[[contrast]] <- x
}


significant_isoforms <- bind_rows(
    isoform_list
)


cat("\n")
cat("====================================================\n")
cat("ALL SIGNIFICANT SWITCHING ISOFORMS\n")
cat("====================================================\n")

cat(
    "Total isoform-contrast rows:",
    nrow(significant_isoforms),
    "\n"
)

cat(
    "Unique switching isoforms:",
    n_distinct(significant_isoforms$Isoform_ID),
    "\n"
)


# ============================================================
# 11. QC: Does every significant isoform belong to a PB gene that is also in the significant switching gene table?
# ============================================================

isoform_gene_qc <- significant_isoforms %>%

    distinct(
        Contrast,
        PB_Gene_ID
    ) %>%

    anti_join(
        significant_genes %>%
            distinct(
                Contrast,
                PB_Gene_ID
            ),
        by = c(
            "Contrast",
            "PB_Gene_ID"
        )
    )


write.csv(
    isoform_gene_qc,
    file.path(
        OUT_DIR,
        "QC_Significant_Isoforms_Without_Significant_Gene_Call.csv"
    ),
    row.names = FALSE
)


cat(
    "\nSignificant isoform PB genes absent from ",
    "significant gene table:",
    nrow(isoform_gene_qc),
    "\n"
)


# ============================================================
# 12. Attach gene-level q-value to significant isoforms
# ============================================================

significant_isoforms <- significant_isoforms %>%

    left_join(

        significant_genes %>%
            dplyr::select(
                Contrast,
                PB_Gene_ID,
                Gene_Switch_Q_Value
            ),

        by = c(
            "Contrast",
            "PB_Gene_ID"
        )

    )


# ============================================================
# 13. Map each switching isoform through SQANTI
# ============================================================

significant_isoforms_mapped <- significant_isoforms %>%

    mutate(

        SQANTI_Isoform_Found =
            Isoform_ID %in%
            sqanti_map$Isoform_ID

    ) %>%

    left_join(
        sqanti_map,
        by = "Isoform_ID"
    ) %>%

    mutate(

        SQANTI_Associated_Gene_Found =
            !is.na(Gene_ID) &
            Gene_ID != ""

    )


# ============================================================
# 14. Mapping diagnostics
# ============================================================

unmapped_isoforms <- significant_isoforms_mapped %>%

    filter(
        !SQANTI_Isoform_Found
    )


missing_associated_gene <- significant_isoforms_mapped %>%

    filter(
        SQANTI_Isoform_Found,
        !SQANTI_Associated_Gene_Found
    )


write.csv(
    unmapped_isoforms,
    file.path(
        OUT_DIR,
        "QC_Significant_DIU_Isoforms_Not_Found_In_SQANTI.csv"
    ),
    row.names = FALSE
)


write.csv(
    missing_associated_gene,
    file.path(
        OUT_DIR,
        "QC_Significant_DIU_Isoforms_Missing_Associated_Gene.csv"
    ),
    row.names = FALSE
)


mapping_summary <- data.frame(

    Metric = c(
        "Significant_isoform_contrast_rows",
        "Exact_SQANTI_isoform_matches",
        "SQANTI_isoform_not_found",
        "SQANTI_match_missing_associated_gene",
        "Successfully_mapped_rows"
    ),

    N = c(

        nrow(
            significant_isoforms_mapped
        ),

        sum(
            significant_isoforms_mapped$
                SQANTI_Isoform_Found
        ),

        nrow(
            unmapped_isoforms
        ),

        nrow(
            missing_associated_gene
        ),

        sum(
            significant_isoforms_mapped$
                SQANTI_Associated_Gene_Found
        )
    )
)


write.csv(
    mapping_summary,
    file.path(
        OUT_DIR,
        "QC_SQANTI_Mapping_Summary.csv"
    ),
    row.names = FALSE
)


cat("\n")
cat("====================================================\n")
cat("SQANTI MAPPING QC\n")
cat("====================================================\n")

print(
    mapping_summary,
    row.names = FALSE
)


# ============================================================
# 15. Save all mapped switching isoforms
# ============================================================

mapped_diu_isoforms <- significant_isoforms_mapped %>%

    filter(
        SQANTI_Associated_Gene_Found
    )


write.csv(
    mapped_diu_isoforms,
    file.path(
        OUT_DIR,
        "All_DIU_Switching_Isoforms_Mapped_To_SQANTI_Genes.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 16. Build ALL mapped DIU genes by contrast
# ============================================================

all_diu_gene_contrast <- mapped_diu_isoforms %>%

    group_by(
        Gene_ID,
        Contrast,
        Contrast_Label,
        condition_1,
        condition_2
    ) %>%

    summarise(

        N_Significant_Switching_Isoforms =
            n_distinct(
                Isoform_ID
            ),

        Switching_Isoforms =
            paste(
                sort(
                    unique(
                        Isoform_ID
                    )
                ),
                collapse = " | "
            ),

        PB_Genes =
            paste(
                sort(
                    unique(
                        PB_Gene_ID
                    )
                ),
                collapse = " | "
            ),

        Min_Gene_Switch_Q =
            if (
                all(
                    is.na(
                        Gene_Switch_Q_Value
                    )
                )
            ) {

                NA_real_

            } else {

                min(
                    Gene_Switch_Q_Value,
                    na.rm = TRUE
                )

            },

        Min_Isoform_Switch_Q =
            if (
                all(
                    is.na(
                        Isoform_Switch_Q_Value
                    )
                )
            ) {

                NA_real_

            } else {

                min(
                    Isoform_Switch_Q_Value,
                    na.rm = TRUE
                )

            },

        Max_Abs_dIF =
            if (
                all(
                    is.na(
                        dIF
                    )
                )
            ) {

                NA_real_

            } else {

                max(
                    abs(dIF),
                    na.rm = TRUE
                )

            },

        .groups = "drop"
    ) %>%

    arrange(
        Gene_ID,
        Contrast
    )


write.csv(
    all_diu_gene_contrast,
    file.path(
        OUT_DIR,
        "All_DIU_Genes_By_Contrast.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 17. Collapse ALL DIU genes across contrasts
# ============================================================

all_diu_any_contrast <- all_diu_gene_contrast %>%

    group_by(
        Gene_ID
    ) %>%

    summarise(

        DIU_Any_Contrast =
            TRUE,

        N_Contrasts_DIU =
            n_distinct(
                Contrast
            ),

        Contrasts_DIU =
            paste(
                contrasts[
                    contrasts %in%
                    unique(
                        Contrast
                    )
                ],
                collapse = " | "
            ),

        Contrast_Labels_DIU =
            paste(
                unname(
                    contrast_labels[
                        contrasts[
                            contrasts %in%
                            unique(
                                Contrast
                            )
                        ]
                    ]
                ),
                collapse = " | "
            ),

        N_Total_Unique_Switching_Isoforms =
            n_distinct(
                unlist(
                    strsplit(
                        Switching_Isoforms,
                        " \\| "
                    )
                )
            ),

        Min_Gene_Switch_Q =
            if (
                all(
                    is.na(
                        Min_Gene_Switch_Q
                    )
                )
            ) {

                NA_real_

            } else {

                min(
                    Min_Gene_Switch_Q,
                    na.rm = TRUE
                )

            },

        Min_Isoform_Switch_Q =
            if (
                all(
                    is.na(
                        Min_Isoform_Switch_Q
                    )
                )
            ) {

                NA_real_

            } else {

                min(
                    Min_Isoform_Switch_Q,
                    na.rm = TRUE
                )

            },

        Max_Abs_dIF =
            if (
                all(
                    is.na(
                        Max_Abs_dIF
                    )
                )
            ) {

                NA_real_

            } else {

                max(
                    Max_Abs_dIF,
                    na.rm = TRUE
                )

            },

        .groups = "drop"
    ) %>%

    arrange(
        desc(
            N_Contrasts_DIU
        ),
        Gene_ID
    )


write.csv(
    all_diu_any_contrast,
    file.path(
        OUT_DIR,
        "All_DIU_Genes_Any_Contrast.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 18. overlay top 50 hub genes with significant switching isoforms
# ============================================================

top50_isoform_long <- mapped_diu_isoforms %>%

    inner_join(
        hubs,
        by = "Gene_ID"
    ) %>%

    arrange(
        Module_Color,
        Gene_ID,
        Contrast,
        desc(
            abs(
                dIF
            )
        )
    )


write.csv(
    top50_isoform_long,
    file.path(
        OUT_DIR,
        "Top50_Hub_DIU_Switching_Isoforms_Long.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 19. make contrast table
# ============================================================

top50_gene_contrast <- all_diu_gene_contrast %>%

    inner_join(
        hubs,
        by = "Gene_ID"
    ) %>%

    arrange(
        Module_Color,
        Gene_ID,
        Contrast
    )


write.csv(
    top50_gene_contrast,
    file.path(
        OUT_DIR,
        "Top50_Hub_Genes_DIU_By_Contrast_Long.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 20. Collapse Top50 DIU information to one row per Gene_ID
# ============================================================

top50_diu_summary <- top50_gene_contrast %>%

    group_by(
        Gene_ID
    ) %>%

    summarise(

        DIU_Any_Contrast =
            TRUE,

        N_Contrasts_DIU =
            n_distinct(
                Contrast
            ),

        Contrasts_DIU =
            paste(
                contrasts[
                    contrasts %in%
                    unique(
                        Contrast
                    )
                ],
                collapse = " | "
            ),

        Contrast_Labels_DIU =
            paste(
                unname(
                    contrast_labels[
                        contrasts[
                            contrasts %in%
                            unique(
                                Contrast
                            )
                        ]
                    ]
                ),
                collapse = " | "
            ),

        PB_Genes_DIU =
            paste(
                sort(
                    unique(
                        unlist(
                            strsplit(
                                PB_Genes,
                                " \\| "
                            )
                        )
                    )
                ),
                collapse = " | "
            ),

        N_Total_Unique_Switching_Isoforms =
            n_distinct(
                unlist(
                    strsplit(
                        Switching_Isoforms,
                        " \\| "
                    )
                )
            ),

        Min_Gene_Switch_Q =
            if (
                all(
                    is.na(
                        Min_Gene_Switch_Q
                    )
                )
            ) {

                NA_real_

            } else {

                min(
                    Min_Gene_Switch_Q,
                    na.rm = TRUE
                )

            },

        Min_Isoform_Switch_Q =
            if (
                all(
                    is.na(
                        Min_Isoform_Switch_Q
                    )
                )
            ) {

                NA_real_

            } else {

                min(
                    Min_Isoform_Switch_Q,
                    na.rm = TRUE
                )

            },

        Max_Abs_dIF =
            if (
                all(
                    is.na(
                        Max_Abs_dIF
                    )
                )
            ) {

                NA_real_

            } else {

                max(
                    Max_Abs_dIF,
                    na.rm = TRUE
                )

            },

        .groups = "drop"
    )


# ============================================================
# 21. Complete Top50 master table
# ============================================================

top50_master <- hubs %>%

    left_join(
        top50_diu_summary,
        by = "Gene_ID"
    ) %>%

    mutate(

        DIU_Any_Contrast =
            ifelse(
                is.na(
                    DIU_Any_Contrast
                ),
                FALSE,
                DIU_Any_Contrast
            ),

        N_Contrasts_DIU =
            ifelse(
                is.na(
                    N_Contrasts_DIU
                ),
                0L,
                N_Contrasts_DIU
            ),

        Contrasts_DIU =
            ifelse(
                is.na(
                    Contrasts_DIU
                ),
                "",
                Contrasts_DIU
            ),

        Contrast_Labels_DIU =
            ifelse(
                is.na(
                    Contrast_Labels_DIU
                ),
                "",
                Contrast_Labels_DIU
            ),

        PB_Genes_DIU =
            ifelse(
                is.na(
                    PB_Genes_DIU
                ),
                "",
                PB_Genes_DIU
            ),

        N_Total_Unique_Switching_Isoforms =
            ifelse(
                is.na(
                    N_Total_Unique_Switching_Isoforms
                ),
                0L,
                N_Total_Unique_Switching_Isoforms
            )

    )


write.csv(
    top50_master,
    file.path(
        OUT_DIR,
        "Top50_Hub_Genes_x_All_DIU_Contrasts_Master.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 22. Top50 genes with DIU only
# ============================================================

top50_diu_only <- top50_master %>%

    filter(
        DIU_Any_Contrast
    ) %>%

    arrange(
        Module_Color,
        desc(
            N_Contrasts_DIU
        ),
        Gene_ID
    )


write.csv(
    top50_diu_only,
    file.path(
        OUT_DIR,
        "Top50_Hub_Genes_With_DIU.csv"
    ),
    row.names = FALSE
)


# ============================================================
# 23. Top50 gene x contrast TRUE/FALSE matrix
# ============================================================

contrast_matrix <- top50_gene_contrast %>%

    distinct(
        Gene_ID,
        Contrast
    ) %>%

    mutate(
        DIU = TRUE
    ) %>%

    pivot_wider(
        names_from = Contrast,
        values_from = DIU,
        values_fill = FALSE
    )


# ------------------------------------------------------------
# Make sure all 9 contrast columns exist
# ------------------------------------------------------------

for (contrast in contrasts) {

    if (
        !contrast %in%
        colnames(contrast_matrix)
    ) {

        contrast_matrix[[contrast]] <- FALSE
    }
}


# ------------------------------------------------------------
# Put contrast columns in defined order
# ------------------------------------------------------------

contrast_matrix <- contrast_matrix %>%

    dplyr::select(
        Gene_ID,
        all_of(contrasts)
    )


# ------------------------------------------------------------
# Join back onto all 150 Top50 hub genes
# ------------------------------------------------------------

top50_contrast_matrix <- hubs %>%

    left_join(
        contrast_matrix,
        by = "Gene_ID"
    )


# ------------------------------------------------------------
# Hub genes with no DIU will have NA after left_join, change those values to FALSE.
# ------------------------------------------------------------

for (contrast in contrasts) {

    top50_contrast_matrix[[contrast]][
        is.na(
            top50_contrast_matrix[[contrast]]
        )
    ] <- FALSE
}


# ------------------------------------------------------------
# Save matrix
# ------------------------------------------------------------

write.csv(
    top50_contrast_matrix,
    file.path(
        OUT_DIR,
        "Top50_Hub_Genes_DIU_Contrast_Matrix.csv"
    ),
    row.names = FALSE
)


cat(
    "\nWritten contrast matrix:\n",
    file.path(
        OUT_DIR,
        "Top50_Hub_Genes_DIU_Contrast_Matrix.csv"
    ),
    "\n"
)


# ============================================================
# 24. Summary by WGCNA module
# ============================================================

if ("Module_Color" %in% colnames(top50_master)) {

    module_summary <- top50_master %>%

        group_by(
            Module_Color
        ) %>%

        summarise(

            Top50_Hub_Genes =
                n_distinct(
                    Gene_ID
                ),

            Hub_Genes_DIU_Any_Contrast =
                n_distinct(
                    Gene_ID[
                        DIU_Any_Contrast
                    ]
                ),

            Percent_Hub_Genes_DIU =
                100 *
                Hub_Genes_DIU_Any_Contrast /
                Top50_Hub_Genes,

            .groups = "drop"
        )


    write.csv(
        module_summary,
        file.path(
            OUT_DIR,
            "Top50_Hub_DIU_Summary_By_Module.csv"
        ),
        row.names = FALSE
    )


    cat("\n")
    cat("====================================================\n")
    cat("DIU HUB GENES BY MODULE\n")
    cat("====================================================\n")

    print(
        module_summary
    )
}


# ============================================================
# 25. Summary by module x contrast
# ============================================================

if ("Module_Color" %in% colnames(hubs)) {

    module_totals <- hubs %>%

        group_by(
            Module_Color
        ) %>%

        summarise(

            Top50_Hub_Genes =
                n_distinct(
                    Gene_ID
                ),

            .groups = "drop"
        )


    module_contrast_counts <- top50_gene_contrast %>%

        group_by(
            Module_Color,
            Contrast
        ) %>%

        summarise(

            Hub_Genes_DIU =
                n_distinct(
                    Gene_ID
                ),

            .groups = "drop"
        )


    module_contrast_summary <- tidyr::crossing(

        Module_Color =
            unique(
                hubs$Module_Color
            ),

        Contrast =
            contrasts

    ) %>%

        left_join(
            module_totals,
            by = "Module_Color"
        ) %>%

        left_join(
            module_contrast_counts,
            by = c(
                "Module_Color",
                "Contrast"
            )
        ) %>%

        mutate(

            Contrast_Label =
                unname(
                    contrast_labels[
                        Contrast
                    ]
                ),

            Hub_Genes_DIU =
                ifelse(
                    is.na(
                        Hub_Genes_DIU
                    ),
                    0L,
                    Hub_Genes_DIU
                ),

            Percent_Hub_Genes_DIU =
                100 *
                Hub_Genes_DIU /
                Top50_Hub_Genes

        )


    write.csv(
        module_contrast_summary,
        file.path(
            OUT_DIR,
            "Top50_Hub_DIU_Summary_By_Module_And_Contrast.csv"
        ),
        row.names = FALSE
    )
}


# ============================================================
# 26. Final summary
# ============================================================

cat("\n")
cat("====================================================\n")
cat("CORRECTED TOP50 HUB x DIU ANALYSIS COMPLETE\n")
cat("====================================================\n")


cat(
    "Significant isoform-contrast rows:",
    nrow(
        significant_isoforms
    ),
    "\n"
)


cat(
    "Successfully SQANTI-mapped rows:",
    nrow(
        mapped_diu_isoforms
    ),
    "\n"
)


cat(
    "Unique SQANTI-associated genes with DIU:",
    n_distinct(
        mapped_diu_isoforms$Gene_ID
    ),
    "\n"
)


cat(
    "Top50 hub genes with DIU in >=1 contrast:",
    n_distinct(
        top50_gene_contrast$Gene_ID
    ),
    "\n"
)


cat(
    "Top50 hub switching-isoform rows:",
    nrow(
        top50_isoform_long
    ),
    "\n"
)


cat("\nMain output directory:\n")
cat(
    OUT_DIR,
    "\n"
)


cat("\nMain output files:\n")

cat(
    "  All_DIU_Genes_By_Contrast.csv\n"
)

cat(
    "  All_DIU_Genes_Any_Contrast.csv\n"
)

cat(
    "  Top50_Hub_Genes_With_DIU.csv\n"
)

cat(
    "  Top50_Hub_Genes_x_All_DIU_Contrasts_Master.csv\n"
)

cat(
    "  Top50_Hub_Genes_DIU_By_Contrast_Long.csv\n"
)

cat(
    "  Top50_Hub_Genes_DIU_Contrast_Matrix.csv\n"
)

cat(
    "  Top50_Hub_DIU_Switching_Isoforms_Long.csv\n"
)

cat(
    "\nQC files begin with QC_ in the same directory.\n"
)

cat("====================================================\n")


