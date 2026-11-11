# Map significant DIE isoforms to shared V. riparia consensus genes and summarize overlap among resistant and susceptible genotypes
# What genes encode differentially expressed isoforms (infested minus control) across all resistant hosts?

#!/usr/bin/env Rscript


# ============================================================
# Resistant genotypes:
#
#   3309   -> Vrip_X_Vrup
#   SW     -> Vrip_X_Vrup
#   5BB    -> Vrip_X_Vber
#   1616   -> Vrip_X_Vace

# DIE significance was already defined in the script 'Differential_isoform_expression_PCA.R', defined as:
#
#   padj < 0.05
#   |log2FC| > 1
# ============================================================


# ============================================================
# 1. LIBRARIES
# ============================================================

suppressPackageStartupMessages({

    library(tidyverse)
})


# ============================================================
# 2. PATHS
# ============================================================

BASE <- "/rhome/wcoll010/shared/pacbio/IsoSeq_analyses_2"


DIE_DIR <- file.path(
    BASE,
    "DE_isoform_analysis"
)


CONV_DIR <- file.path(
    BASE,
    "Vrip_resistance_convergence"
)


OUT_DIR <- file.path(
    CONV_DIR,
    "DIE_gene_convergence"
)


dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


CONSENSUS_MAP_FILE <- file.path(
    CONV_DIR,
    "Vrip_gene_consensus_map.tsv"
)


# ============================================================
# 3. GENOTYPE / REFERENCE DEFINITIONS
# ============================================================

genotype_info <- tribble(

    ~Genotype, ~DIE_Name, ~Reference,      ~Phenotype,

    "3309",    "3309",    "Vrip_X_Vrup", "Resistant",
    "SW",      "sw",      "Vrip_X_Vrup", "Resistant",
    "5BB",     "5BB",     "Vrip_X_Vber", "Resistant",
    "1616",    "1616",    "Vrip_X_Vace", "Resistant",

    "101-14",  "101-14",  "Vrip_X_Vrup", "Susceptible",
    "SO4",     "SO4",     "Vrip_X_Vber", "Susceptible"
)


resistant_genotypes <- c(
    "3309",
    "SW",
    "5BB",
    "1616"
)


susceptible_genotypes <- c(
    "101-14",
    "SO4"
)


# ============================================================
# 4. READ CONSENSUS GENE MAP
# ============================================================

consensus_map <- read.delim(

    CONSENSUS_MAP_FILE,

    stringsAsFactors = FALSE,

    check.names = FALSE
)


required_map_cols <- c(
    "Reference",
    "Original_Gene_ID",
    "Consensus_Gene_ID"
)


missing_map_cols <- setdiff(
    required_map_cols,
    colnames(consensus_map)
)


if (
    length(missing_map_cols) > 0
) {

    stop(
        "Missing columns in consensus map:\n",
        paste(
            missing_map_cols,
            collapse = "\n"
        )
    )
}


gene_map <- consensus_map %>%

    select(
        Reference,
        Original_Gene_ID,
        Consensus_Gene_ID
    ) %>%

    distinct()


# ============================================================
# 5. READ SIGNIFICANT DIE ISOFORMS FOR ONE GENOTYPE
# ============================================================

read_genotype_DIE <- function(
    genotype,
    die_name,
    reference,
    phenotype
) {

    cat(
        "\n============================================================\n"
    )

    cat(
        genotype,
        "\n"
    )

    cat(
        "============================================================\n"
    )


    result_dir <- file.path(
        DIE_DIR,
        die_name,
        "results"
    )


    up_file <- file.path(
        result_dir,
        paste0(
            die_name,
            "_upregulated_isoforms.csv"
        )
    )


    down_file <- file.path(
        result_dir,
        paste0(
            die_name,
            "_downregulated_isoforms.csv"
        )
    )


    if (!file.exists(up_file)) {

        stop(
            "Missing UP file:\n",
            up_file
        )
    }


    if (!file.exists(down_file)) {

        stop(
            "Missing DOWN file:\n",
            down_file
        )
    }


    up <- read.csv(
        up_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    ) %>%

        mutate(
            DIE_Direction =
                "UP"
        )


    down <- read.csv(
        down_file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    ) %>%

        mutate(
            DIE_Direction =
                "DOWN"
        )


    die <- bind_rows(
        up,
        down
    )


    cat(
        "Significant DE isoforms: ",
        nrow(die),
        "\n",
        sep = ""
    )


    required_cols <- c(
        "isoform",
        "associated_gene",
        "log2FoldChange",
        "padj"
    )


    missing_cols <- setdiff(
        required_cols,
        colnames(die)
    )


    if (
        length(missing_cols) > 0
    ) {

        stop(
            genotype,
            ": missing columns:\n",
            paste(
                missing_cols,
                collapse = "\n"
            )
        )
    }


    local_map <- gene_map %>%

        filter(
            Reference ==
                reference
        )


    mapped <- die %>%

        rename(
            Original_Gene_ID =
                associated_gene
        ) %>%

        inner_join(
            local_map,
            by =
                "Original_Gene_ID"
        ) %>%

        mutate(

            Genotype =
                genotype,

            Reference =
                reference,

            Phenotype =
                phenotype
        )


    cat(
        "Mapped significant isoforms: ",
        nrow(mapped),
        "\n",
        sep = ""
    )


    cat(
        "Consensus genes encoding >=1 DE isoform: ",
        n_distinct(
            mapped$Consensus_Gene_ID
        ),
        "\n",
        sep = ""
    )


    mapped
}


# ============================================================
# 6. LOAD ALL SIX GENOTYPES
# ============================================================

all_DIE_isoforms <- bind_rows(

    lapply(

        seq_len(
            nrow(genotype_info)
        ),

        function(i) {

            read_genotype_DIE(

                genotype =
                    genotype_info$Genotype[i],

                die_name =
                    genotype_info$DIE_Name[i],

                reference =
                    genotype_info$Reference[i],

                phenotype =
                    genotype_info$Phenotype[i]
            )
        }
    )
)


write.csv(

    all_DIE_isoforms,

    file.path(
        OUT_DIR,
        "All_genotypes_significant_DIE_isoforms_consensus_mapped.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 7. COLLAPSE TO GENE x GENOTYPE
#
# A gene can encode:
#
#   only UP DE isoforms
#   only DOWN DE isoforms
#   both UP and DOWN DE isoforms
# ============================================================

gene_genotype <- all_DIE_isoforms %>%

    group_by(
        Consensus_Gene_ID,
        Genotype,
        Phenotype
    ) %>%

    summarise(

        N_DE_Isoforms =
            n_distinct(
                isoform
            ),

        N_UP_Isoforms =
            n_distinct(
                isoform[
                    DIE_Direction ==
                        "UP"
                ]
            ),

        N_DOWN_Isoforms =
            n_distinct(
                isoform[
                    DIE_Direction ==
                        "DOWN"
                ]
            ),

        Has_UP =
            N_UP_Isoforms > 0,

        Has_DOWN =
            N_DOWN_Isoforms > 0,

        DIE_Pattern =
            case_when(

                Has_UP &
                Has_DOWN ~
                    "UP_and_DOWN",

                Has_UP ~
                    "UP_only",

                Has_DOWN ~
                    "DOWN_only",

                TRUE ~
                    "None"
            ),

        .groups =
            "drop"
    )


write.csv(

    gene_genotype,

    file.path(
        OUT_DIR,
        "DIE_gene_by_genotype.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 8. RESISTANT-HOST GENE OVERLAP
# ============================================================

resistant_gene_matrix <- gene_genotype %>%

    filter(
        Genotype %in%
            resistant_genotypes
    ) %>%

    mutate(
        Has_DIE =
            TRUE
    ) %>%

    select(
        Consensus_Gene_ID,
        Genotype,
        Has_DIE
    ) %>%

    pivot_wider(

        names_from =
            Genotype,

        values_from =
            Has_DIE,

        values_fill =
            FALSE
    )


# Make sure all columns exist

for (
    g in resistant_genotypes
) {

    if (
        !g %in%
        colnames(
            resistant_gene_matrix
        )
    ) {

        resistant_gene_matrix[[g]] <-
            FALSE
    }
}


resistant_gene_matrix <- resistant_gene_matrix %>%

    rowwise() %>%

    mutate(

        N_Resistant_With_DIE =
            sum(
                c_across(
                    all_of(
                        resistant_genotypes
                    )
                )
            )
    ) %>%

    ungroup() %>%

    arrange(
        desc(
            N_Resistant_With_DIE
        ),
        Consensus_Gene_ID
    )


write.csv(

    resistant_gene_matrix,

    file.path(
        OUT_DIR,
        "Resistant_DIE_gene_overlap_matrix.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 9. OVERLAP SUMMARY
# ============================================================

overlap_summary <- resistant_gene_matrix %>%

    count(
        N_Resistant_With_DIE,
        name =
            "N_Genes"
    ) %>%

    arrange(
        desc(
            N_Resistant_With_DIE
        )
    )


write.csv(

    overlap_summary,

    file.path(
        OUT_DIR,
        "Resistant_DIE_gene_overlap_summary.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 10. SHARED GENE SETS
# ============================================================

shared_all4 <- resistant_gene_matrix %>%

    filter(
        N_Resistant_With_DIE == 4
    )


shared_atleast3 <- resistant_gene_matrix %>%

    filter(
        N_Resistant_With_DIE >= 3
    )


shared_atleast2 <- resistant_gene_matrix %>%

    filter(
        N_Resistant_With_DIE >= 2
    )


write.csv(

    shared_all4,

    file.path(
        OUT_DIR,
        "DIE_genes_shared_ALL4_resistant.csv"
    ),

    row.names = FALSE
)


write.csv(

    shared_atleast3,

    file.path(
        OUT_DIR,
        "DIE_genes_shared_atleast3_resistant.csv"
    ),

    row.names = FALSE
)


write.csv(

    shared_atleast2,

    file.path(
        OUT_DIR,
        "DIE_genes_shared_atleast2_resistant.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 11. DIRECTIONAL CONVERGENCE
#
# Ask whether the same consensus gene encodes an UP or DOWN
# DE isoform across multiple resistant genotypes.
#
# This still does NOT imply that the exact same PB transcript
# is shared across catalogs.
# ============================================================

resistant_direction <- gene_genotype %>%

    filter(
        Genotype %in%
            resistant_genotypes
    ) %>%

    select(
        Consensus_Gene_ID,
        Genotype,
        Has_UP,
        Has_DOWN,
        DIE_Pattern
    )


up_matrix <- resistant_direction %>%

    select(
        Consensus_Gene_ID,
        Genotype,
        Has_UP
    ) %>%

    pivot_wider(
        names_from =
            Genotype,
        values_from =
            Has_UP,
        values_fill =
            FALSE,
        names_prefix =
            "UP_"
    )


down_matrix <- resistant_direction %>%

    select(
        Consensus_Gene_ID,
        Genotype,
        Has_DOWN
    ) %>%

    pivot_wider(
        names_from =
            Genotype,
        values_from =
            Has_DOWN,
        values_fill =
            FALSE,
        names_prefix =
            "DOWN_"
    )


directional <- full_join(
    up_matrix,
    down_matrix,
    by =
        "Consensus_Gene_ID"
)


for (
    g in resistant_genotypes
) {

    up_col <- paste0(
        "UP_",
        g
    )

    down_col <- paste0(
        "DOWN_",
        g
    )


    if (
        !up_col %in%
        colnames(
            directional
        )
    ) {

        directional[[up_col]] <-
            FALSE
    }


    if (
        !down_col %in%
        colnames(
            directional
        )
    ) {

        directional[[down_col]] <-
            FALSE
    }


    directional[[up_col]][
        is.na(
            directional[[up_col]]
        )
    ] <- FALSE


    directional[[down_col]][
        is.na(
            directional[[down_col]]
        )
    ] <- FALSE
}


directional <- directional %>%

    rowwise() %>%

    mutate(

        N_Resistant_With_UP_DIE =
            sum(
                c_across(
                    all_of(
                        paste0(
                            "UP_",
                            resistant_genotypes
                        )
                    )
                )
            ),

        N_Resistant_With_DOWN_DIE =
            sum(
                c_across(
                    all_of(
                        paste0(
                            "DOWN_",
                            resistant_genotypes
                        )
                    )
                )
            )
    ) %>%

    ungroup() %>%

    arrange(
        desc(
            N_Resistant_With_UP_DIE
        ),
        desc(
            N_Resistant_With_DOWN_DIE
        )
    )


write.csv(

    directional,

    file.path(
        OUT_DIR,
        "Resistant_DIE_directional_convergence.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 12. SUSCEPTIBLE OVERLAP
#
# This lets us determine whether genes shared among resistant
# hosts are also DIE-responsive in susceptible hosts.
# ============================================================

susceptible_gene_matrix <- gene_genotype %>%

    filter(
        Genotype %in%
            susceptible_genotypes
    ) %>%

    mutate(
        Has_DIE =
            TRUE
    ) %>%

    select(
        Consensus_Gene_ID,
        Genotype,
        Has_DIE
    ) %>%

    pivot_wider(

        names_from =
            Genotype,

        values_from =
            Has_DIE,

        values_fill =
            FALSE,

        names_prefix =
            "Susceptible_"
    )


combined <- resistant_gene_matrix %>%

    left_join(
        susceptible_gene_matrix,
        by =
            "Consensus_Gene_ID"
    )


for (
    g in susceptible_genotypes
) {

    col <- paste0(
        "Susceptible_",
        g
    )


    if (
        !col %in%
        colnames(
            combined
        )
    ) {

        combined[[col]] <-
            FALSE
    }


    combined[[col]][
        is.na(
            combined[[col]]
        )
    ] <- FALSE
}


combined <- combined %>%

    rowwise() %>%

    mutate(

        N_Susceptible_With_DIE =
            sum(
                c_across(
                    all_of(
                        paste0(
                            "Susceptible_",
                            susceptible_genotypes
                        )
                    )
                )
            ),

        Resistant_Specific_All4 =
            N_Resistant_With_DIE == 4 &
            N_Susceptible_With_DIE == 0,

        Resistant_Specific_AtLeast3 =
            N_Resistant_With_DIE >= 3 &
            N_Susceptible_With_DIE == 0
    ) %>%

    ungroup()


write.csv(

    combined,

    file.path(
        OUT_DIR,
        "DIE_resistant_vs_susceptible_gene_convergence.csv"
    ),

    row.names = FALSE
)


write.csv(

    combined %>%

        filter(
            Resistant_Specific_All4
        ),

    file.path(
        OUT_DIR,
        "DIE_genes_ALL4_resistant_ZERO_susceptible.csv"
    ),

    row.names = FALSE
)


write.csv(

    combined %>%

        filter(
            Resistant_Specific_AtLeast3
        ),

    file.path(
        OUT_DIR,
        "DIE_genes_atleast3_resistant_ZERO_susceptible.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 13. SIMPLE SUMMARY FIGURE
# ============================================================

plot_data <- overlap_summary %>%

    mutate(

        Category =
            factor(

                N_Resistant_With_DIE,

                levels =
                    c(
                        4,
                        3,
                        2,
                        1
                    ),

                labels =
                    c(
                        "4 of 4",
                        "3 of 4",
                        "2 of 4",
                        "1 of 4"
                    )
            )
    )


p <- ggplot(

    plot_data,

    aes(
        x =
            Category,
        y =
            N_Genes
    )
) +

    geom_col(
        width =
            0.7
    ) +

    geom_text(

        aes(
            label =
                N_Genes
        ),

        vjust =
            -0.3
    ) +

    labs(

        title =
            paste0(
                "Shared V. riparia genes encoding ",
                "differentially expressed isoforms"
            ),

        x =
            "Number of resistant genotypes with DIE",

        y =
            "Consensus genes"
    ) +

    theme_bw(
        base_size =
            12
    ) +

    theme(
        panel.grid =
            element_blank()
    )


ggsave(

    file.path(
        OUT_DIR,
        "Resistant_DIE_gene_overlap.pdf"
    ),

    p,

    width =
        7,

    height =
        5
)


# ============================================================
# 14. PRINT SUMMARY
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "DIE GENE CONVERGENCE RESULTS\n"
)

cat(
    "============================================================\n\n"
)


print(
    overlap_summary,
    n = Inf,
    width = Inf
)


cat(
    "\nGenes shared by all 4 resistant genotypes: ",
    nrow(
        shared_all4
    ),
    "\n",
    sep = ""
)


cat(
    "Genes shared by >=3 resistant genotypes: ",
    nrow(
        shared_atleast3
    ),
    "\n",
    sep = ""
)


cat(
    "Genes shared by >=2 resistant genotypes: ",
    nrow(
        shared_atleast2
    ),
    "\n",
    sep = ""
)


cat(
    "\nAll-4 resistant / zero-susceptible genes: ",
    sum(
        combined$Resistant_Specific_All4
    ),
    "\n",
    sep = ""
)


cat(
    ">=3 resistant / zero-susceptible genes: ",
    sum(
        combined$Resistant_Specific_AtLeast3
    ),
    "\n",
    sep = ""
)


cat(
    "\nOutput:\n",
    OUT_DIR,
    "\n"
)












### Test for statistical significance of gene overlap /convergence

# Perform genotype-specific permutation test of DIE convergence with 100,000 permutations


#!/usr/bin/env Rscript


# ============================================================
# 1. LIBRARIES
# ============================================================

suppressPackageStartupMessages({

    library(tidyverse)
})


# ============================================================
# 2. PATHS
# ============================================================

BASE <- "/rhome/wcoll010/shared/pacbio/IsoSeq_analyses_2"


CONV_DIR <- file.path(
    BASE,
    "Vrip_resistance_convergence"
)


STEP10_DIR <- file.path(
    CONV_DIR,
    "DIE_gene_convergence"
)


DIE_BASE <- file.path(
    BASE,
    "DE_isoform_analysis"
)


OUT_DIR <- file.path(
    CONV_DIR,
    "DIE_gene_convergence_test"
)


dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


GENE_GENOTYPE_FILE <- file.path(
    STEP10_DIR,
    "DIE_gene_by_genotype.csv"
)


COMBINED_FILE <- file.path(
    STEP10_DIR,
    "DIE_resistant_vs_susceptible_gene_convergence.csv"
)


ISOFORM_FILE <- file.path(
    STEP10_DIR,
    "All_genotypes_significant_DIE_isoforms_consensus_mapped.csv"
)


CONSENSUS_MAP_FILE <- file.path(
    CONV_DIR,
    "Vrip_gene_consensus_map.tsv"
)


SQANTI_FILE <- file.path(
    BASE,
    "sqanti_results",
    "Vrip_X_Vrup_sqanti",
    "filter_results",
    "Vrip_X_Vrup_filtered_RulesFilter_result_classification.txt"
)


EGGNOG_FILE <- file.path(
    BASE,
    "eggnog_results",
    "Vrip_X_Vrup_eggnog.emapper.annotations"
)


# ============================================================
# 3. PARAMETERS
# ============================================================

N_PERMUTATIONS <- 100000

SEED <- 12345


resistant_genotypes <- c(
    "3309",
    "SW",
    "5BB",
    "1616"
)


susceptible_genotypes <- c(
    "101-14",
    "SO4"
)


# DIE directory names differ only for SW

die_names <- c(

    "3309" =
        "3309",

    "SW" =
        "sw",

    "5BB" =
        "5BB",

    "1616" =
        "1616",

    "101-14" =
        "101-14",

    "SO4" =
        "SO4"
)


reference_names <- c(

    "3309" =
        "Vrip_X_Vrup",

    "SW" =
        "Vrip_X_Vrup",

    "101-14" =
        "Vrip_X_Vrup",

    "5BB" =
        "Vrip_X_Vber",

    "SO4" =
        "Vrip_X_Vber",

    "1616" =
        "Vrip_X_Vace"
)


# ============================================================
# 4. CHECK CORE INPUTS
# ============================================================

required_files <- c(

    GENE_GENOTYPE_FILE,
    COMBINED_FILE,
    ISOFORM_FILE,
    CONSENSUS_MAP_FILE,
    SQANTI_FILE,
    EGGNOG_FILE
)


missing_files <- required_files[
    !file.exists(
        required_files
    )
]


if (
    length(
        missing_files
    ) > 0
) {

    stop(
        "\nMissing required file(s):\n\n",
        paste(
            missing_files,
            collapse = "\n"
        ),
        "\n"
    )
}


# ============================================================
# 5. READ STEP-10 RESULTS
# ============================================================

gene_genotype <- read.csv(

    GENE_GENOTYPE_FILE,

    stringsAsFactors = FALSE,

    check.names = FALSE
)


combined <- read.csv(

    COMBINED_FILE,

    stringsAsFactors = FALSE,

    check.names = FALSE
)


all_sig_isoforms <- read.csv(

    ISOFORM_FILE,

    stringsAsFactors = FALSE,

    check.names = FALSE
)


# ============================================================
# 6. READ CONSENSUS GENE MAP
# ============================================================

consensus_map <- read.delim(

    CONSENSUS_MAP_FILE,

    stringsAsFactors = FALSE,

    check.names = FALSE
)


gene_map <- consensus_map %>%

    dplyr::select(
        Reference,
        Original_Gene_ID,
        Consensus_Gene_ID
    ) %>%

    distinct()


# ============================================================
# 7. BUILD GENOTYPE-SPECIFIC DIE-TESTABLE GENE UNIVERSES
#
# The *_all_isoform_DE_results.csv tables contain the
# expression-filtered isoforms that entered DESeq2.
#
# Therefore these are a direct record of DIE testability.
# ============================================================

get_testable_genes <- function(
    genotype
) {

    die_name <- die_names[[
        genotype
    ]]


    reference <- reference_names[[
        genotype
    ]]


    f <- file.path(

        DIE_BASE,

        die_name,

        "results",

        paste0(
            die_name,
            "_all_isoform_DE_results.csv"
        )
    )


    if (
        !file.exists(
            f
        )
    ) {

        stop(
            "\nMissing DIE all-results file for ",
            genotype,
            ":\n",
            f,
            "\n"
        )
    }


    x <- read.csv(

        f,

        stringsAsFactors = FALSE,

        check.names = FALSE
    )


    if (
        !"associated_gene" %in%
        colnames(
            x
        )
    ) {

        stop(
            genotype,
            ": associated_gene column not found in:\n",
            f
        )
    }


    local_map <- gene_map %>%

        filter(
            Reference ==
                reference
        )


    genes <- x %>%

        filter(
            !is.na(
                associated_gene
            ),
            associated_gene != ""
        ) %>%

        dplyr::select(
            associated_gene
        ) %>%

        distinct() %>%

        inner_join(

            local_map,

            by = c(
                "associated_gene" =
                    "Original_Gene_ID"
            )
        ) %>%

        pull(
            Consensus_Gene_ID
        ) %>%

        unique()


    genes <- genes[
        !is.na(
            genes
        ) &
        genes != ""
    ]


    cat(
        genotype,
        ": ",
        length(
            genes
        ),
        " DIE-testable Vrip consensus genes\n",
        sep = ""
    )


    genes
}


cat(
    "\n============================================================\n"
)

cat(
    "BUILDING DIE-TESTABLE GENE UNIVERSES\n"
)

cat(
    "============================================================\n\n"
)


testable_sets <- lapply(

    resistant_genotypes,

    get_testable_genes
)


names(
    testable_sets
) <- resistant_genotypes


# ============================================================
# 8. OBSERVED RESISTANT DIE GENE SETS
# ============================================================

observed_sets <- list()


cat(
    "\n============================================================\n"
)

cat(
    "OBSERVED DIE GENE SETS\n"
)

cat(
    "============================================================\n\n"
)


for (
    g in resistant_genotypes
) {

    genes <- gene_genotype %>%

        filter(
            Genotype ==
                g
        ) %>%

        pull(
            Consensus_Gene_ID
        ) %>%

        unique()


    # Every observed DIE gene should have been testable.

    not_testable <- setdiff(

        genes,

        testable_sets[[
            g
        ]]
    )


    if (
        length(
            not_testable
        ) > 0
    ) {

        stop(
            "\n",
            g,
            ": ",
            length(
                not_testable
            ),
            " observed DIE genes are absent from its ",
            "testable universe.\n",
            "This indicates an ID-mapping inconsistency.\n"
        )
    }


    observed_sets[[
        g
    ]] <- genes


    cat(
        g,
        ": ",
        length(
            genes
        ),
        " observed DIE genes\n",
        sep = ""
    )
}


# ============================================================
# 9. OBSERVED OVERLAP COUNTS
# ============================================================

all_observed <- unique(
    unlist(
        observed_sets
    )
)


observed_counts <- integer(
    length(
        all_observed
    )
)


names(
    observed_counts
) <- all_observed


for (
    g in resistant_genotypes
) {

    observed_counts[
        observed_sets[[
            g
        ]]
    ] <-
        observed_counts[
            observed_sets[[
                g
            ]]
        ] +
        1L
}


OBS_ALL4 <- sum(
    observed_counts ==
        4
)


OBS_ATLEAST3 <- sum(
    observed_counts >=
        3
)


OBS_ATLEAST2 <- sum(
    observed_counts >=
        2
)


cat(
    "\nObserved genes shared by all 4: ",
    OBS_ALL4,
    "\n",
    sep = ""
)


cat(
    "Observed genes shared by >=3: ",
    OBS_ATLEAST3,
    "\n",
    sep = ""
)


cat(
    "Observed genes shared by >=2: ",
    OBS_ATLEAST2,
    "\n",
    sep = ""
)


# ============================================================
# 10. BUILD INTEGER REPRESENTATION FOR FAST PERMUTATION
# ============================================================

master_universe <- sort(

    unique(

        unlist(
            testable_sets
        )
    )
)


gene_to_index <- seq_along(
    master_universe
)


names(
    gene_to_index
) <- master_universe


testable_indices <- lapply(

    testable_sets,

    function(x) {

        unname(
            gene_to_index[
                x
            ]
        )
    }
)


observed_sizes <- sapply(

    observed_sets,

    length
)


testable_sizes <- sapply(

    testable_sets,

    length
)


qc_table <- tibble(

    Genotype =
        resistant_genotypes,

    DIE_Testable_Vrip_Genes =
        unname(
            testable_sizes[
                resistant_genotypes
            ]
        ),

    Observed_DIE_Genes =
        unname(
            observed_sizes[
                resistant_genotypes
            ]
        ),

    Percent_Testable_Genes_DIE =
        100 *
        Observed_DIE_Genes /
        DIE_Testable_Vrip_Genes
)


write.csv(

    qc_table,

    file.path(
        OUT_DIR,
        "DIE_genotype_testable_universe_summary.csv"
    ),

    row.names = FALSE
)


cat(
    "\nGenotype-specific testable universes:\n\n"
)


print(
    qc_table,
    n = Inf,
    width = Inf
)


# ============================================================
# 11. PERMUTATION
#
# Each genotype is sampled from ITS OWN testable gene universe
# while preserving its observed DIE-set size.
# ============================================================

set.seed(
    SEED
)


null_all4 <- integer(
    N_PERMUTATIONS
)


null_atleast3 <- integer(
    N_PERMUTATIONS
)


null_atleast2 <- integer(
    N_PERMUTATIONS
)


N_MASTER <- length(
    master_universe
)


cat(
    "\nRunning ",
    format(
        N_PERMUTATIONS,
        big.mark = ","
    ),
    " permutations...\n",
    sep = ""
)


for (
    i in seq_len(
        N_PERMUTATIONS
    )
) {

    counts <- integer(
        N_MASTER
    )


    for (
        g in resistant_genotypes
    ) {

        sampled <- sample(

            testable_indices[[
                g
            ]],

            size =
                observed_sizes[[
                    g
                ]],

            replace =
                FALSE
        )


        counts[
            sampled
        ] <-
            counts[
                sampled
            ] +
            1L
    }


    null_all4[i] <-
        sum(
            counts ==
                4L
        )


    null_atleast3[i] <-
        sum(
            counts >=
                3L
        )


    null_atleast2[i] <-
        sum(
            counts >=
                2L
        )
}


# ============================================================
# 12. EMPIRICAL P-VALUES
# ============================================================

empirical_p <- function(
    observed,
    null
) {

    (
        sum(
            null >=
                observed
        ) +
        1
    ) /
    (
        length(
            null
        ) +
        1
    )
}


summarize_null <- function(
    criterion,
    observed,
    null
) {

    null_mean <- mean(
        null
    )


    tibble(

        Criterion =
            criterion,

        Observed =
            observed,

        Null_Mean =
            null_mean,

        Null_Median =
            median(
                null
            ),

        Null_SD =
            sd(
                null
            ),

        Null_2.5pct =
            as.numeric(
                quantile(
                    null,
                    0.025
                )
            ),

        Null_97.5pct =
            as.numeric(
                quantile(
                    null,
                    0.975
                )
            ),

        Fold_Over_Null =
            ifelse(
                null_mean >
                    0,
                observed /
                    null_mean,
                NA_real_
            ),

        Empirical_P =
            empirical_p(
                observed,
                null
            )
    )
}


permutation_summary <- bind_rows(

    summarize_null(
        "All_4_resistant",
        OBS_ALL4,
        null_all4
    ),

    summarize_null(
        "At_least_3_of_4_resistant",
        OBS_ATLEAST3,
        null_atleast3
    ),

    summarize_null(
        "At_least_2_of_4_resistant",
        OBS_ATLEAST2,
        null_atleast2
    )
)


write.csv(

    permutation_summary,

    file.path(
        OUT_DIR,
        "DIE_gene_convergence_permutation_summary.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 13. SAVE NULL DISTRIBUTIONS
# ============================================================

null_df <- tibble(

    Permutation =
        seq_len(
            N_PERMUTATIONS
        ),

    Shared_All4 =
        null_all4,

    Shared_AtLeast3 =
        null_atleast3,

    Shared_AtLeast2 =
        null_atleast2
)


write.csv(

    null_df,

    file.path(
        OUT_DIR,
        "DIE_gene_convergence_null_distributions.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 14. READ EGGNOG
# ============================================================

read_eggnog <- function(
    file
) {

    lines <- readLines(
        file,
        warn = FALSE
    )


    header_line <- grep(
        "^#query\t",
        lines
    )


    if (
        length(
            header_line
        ) != 1
    ) {

        stop(
            "Could not identify eggNOG #query header."
        )
    }


    x <- read.delim(

        file,

        skip =
            header_line - 1,

        header =
            TRUE,

        sep =
            "\t",

        quote =
            "",

        comment.char =
            "",

        stringsAsFactors =
            FALSE,

        check.names =
            FALSE
    )


    colnames(
        x
    )[1] <-
        "isoform"


    x
}


eggnog <- read_eggnog(
    EGGNOG_FILE
)


# ============================================================
# 15. MAP EGGNOG ISOFORMS TO CONSENSUS GENES
# ============================================================

sqanti <- read.delim(

    SQANTI_FILE,

    header = TRUE,

    sep = "\t",

    quote = "",

    stringsAsFactors = FALSE,

    check.names = FALSE
)


isoform_gene <- sqanti %>%

    filter(
        filter_result ==
            "Isoform"
    ) %>%

    dplyr::select(
        isoform,
        associated_gene
    ) %>%

    distinct()


vrup_map <- gene_map %>%

    filter(
        Reference ==
            "Vrip_X_Vrup"
    ) %>%

    dplyr::select(
        Original_Gene_ID,
        Consensus_Gene_ID
    ) %>%

    distinct()


annotation_rows <- isoform_gene %>%

    inner_join(
        eggnog,
        by =
            "isoform"
    ) %>%

    inner_join(

        vrup_map,

        by = c(
            "associated_gene" =
                "Original_Gene_ID"
        )
    )


# ============================================================
# 16. COLLAPSE FUNCTIONAL ANNOTATIONS BY GENE
# ============================================================

collapse_values <- function(
    x
) {

    x <- as.character(
        x
    )


    x <- x[
        !is.na(
            x
        ) &
        x != "" &
        x != "-"
    ]


    x <- unique(
        x
    )


    if (
        length(
            x
        ) ==
        0
    ) {

        return(
            NA_character_
        )
    }


    paste(
        x,
        collapse =
            " | "
    )
}


gene_annotation <- annotation_rows %>%

    group_by(
        Consensus_Gene_ID
    ) %>%

    summarise(

        Preferred_Name =
            collapse_values(
                Preferred_name
            ),

        Description =
            collapse_values(
                Description
            ),

        GO_Terms =
            collapse_values(
                GOs
            ),

        KEGG_Pathways =
            collapse_values(
                KEGG_Pathway
            ),

        PFAMs =
            collapse_values(
                PFAMs
            ),

        .groups =
            "drop"
    )


# ============================================================
# 17. DIE DIRECTION / NUMBER OF ISOFORMS BY GENOTYPE
# ============================================================

direction_table <- gene_genotype %>%

    dplyr::select(

        Consensus_Gene_ID,
        Genotype,
        N_DE_Isoforms,
        N_UP_Isoforms,
        N_DOWN_Isoforms,
        DIE_Pattern
    ) %>%

    pivot_wider(

        names_from =
            Genotype,

        values_from =
            c(
                N_DE_Isoforms,
                N_UP_Isoforms,
                N_DOWN_Isoforms,
                DIE_Pattern
            ),

        names_sep =
            "_"
    )


# ============================================================
# 18. EXACT SIGNIFICANT ISOFORM IDs
# ============================================================

isoform_table <- all_sig_isoforms %>%

    group_by(
        Consensus_Gene_ID,
        Genotype
    ) %>%

    summarise(

        Significant_Isoforms =
            paste(

                sort(
                    unique(
                        isoform
                    )
                ),

                collapse =
                    ";"
            ),

        .groups =
            "drop"
    ) %>%

    pivot_wider(

        names_from =
            Genotype,

        values_from =
            Significant_Isoforms,

        names_prefix =
            "Significant_Isoforms_"
    )


# ============================================================
# 19. BUILD ANNOTATED CONVERGENCE TABLE
# ============================================================

convergence_table <- combined %>%

    left_join(
        direction_table,
        by =
            "Consensus_Gene_ID"
    ) %>%

    left_join(
        isoform_table,
        by =
            "Consensus_Gene_ID"
    ) %>%

    left_join(
        gene_annotation,
        by =
            "Consensus_Gene_ID"
    ) %>%

    mutate(

        Convergence_Class =
            case_when(

                N_Resistant_With_DIE ==
                    4 &
                N_Susceptible_With_DIE ==
                    0 ~
                    "4R_0S",

                N_Resistant_With_DIE ==
                    4 ~
                    "4R",

                N_Resistant_With_DIE ==
                    3 &
                N_Susceptible_With_DIE ==
                    0 ~
                    "3R_0S",

                N_Resistant_With_DIE ==
                    3 ~
                    "3R",

                N_Resistant_With_DIE ==
                    2 ~
                    "2R",

                TRUE ~
                    "1R"
            )
    ) %>%

    arrange(

        desc(
            N_Resistant_With_DIE
        ),

        N_Susceptible_With_DIE,

        Consensus_Gene_ID
    )


write.csv(

    convergence_table,

    file.path(
        OUT_DIR,
        "All_DIE_convergence_genes_annotated.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 20. WRITE KEY CONVERGENT GENE TABLES
# ============================================================

all4 <- convergence_table %>%

    filter(
        N_Resistant_With_DIE ==
            4
    )


all4_0S <- convergence_table %>%

    filter(

        N_Resistant_With_DIE ==
            4,

        N_Susceptible_With_DIE ==
            0
    )


atleast3 <- convergence_table %>%

    filter(
        N_Resistant_With_DIE >=
            3
    )


atleast3_0S <- convergence_table %>%

    filter(

        N_Resistant_With_DIE >=
            3,

        N_Susceptible_With_DIE ==
            0
    )


write.csv(

    all4,

    file.path(
        OUT_DIR,
        "Convergent_DIE_genes_ALL4_resistant_annotated.csv"
    ),

    row.names = FALSE
)


write.csv(

    all4_0S,

    file.path(
        OUT_DIR,
        "Convergent_DIE_genes_ALL4_resistant_ZERO_susceptible_annotated.csv"
    ),

    row.names = FALSE
)


write.csv(

    atleast3,

    file.path(
        OUT_DIR,
        "Convergent_DIE_genes_ATLEAST3_resistant_annotated.csv"
    ),

    row.names = FALSE
)


write.csv(

    atleast3_0S,

    file.path(
        OUT_DIR,
        "Convergent_DIE_genes_ATLEAST3_resistant_ZERO_susceptible_annotated.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 21. NULL-DISTRIBUTION PLOTS
# ============================================================

plot_null <- function(
    null,
    observed,
    title,
    filename
) {

    df <- tibble(
        Null_Overlap =
            null
    )


    p <- ggplot(

        df,

        aes(
            x =
                Null_Overlap
        )
    ) +

        geom_histogram(
            binwidth =
                1,
            boundary =
                -0.5
        ) +

        geom_vline(
            xintercept =
                observed,
            linewidth =
                1
        ) +

        labs(

            title =
                title,

            subtitle =
                paste0(
                    "Observed = ",
                    observed
                ),

            x =
                "Number of overlapping genes under null",

            y =
                "Permutations"
        ) +

        theme_bw(
            base_size =
                12
        )


    ggsave(

        file.path(
            OUT_DIR,
            filename
        ),

        p,

        width =
            7,

        height =
            5
    )
}


plot_null(

    null_all4,

    OBS_ALL4,

    "DIE convergence across all four resistant genotypes",

    "Permutation_null_ALL4_resistant.pdf"
)


plot_null(

    null_atleast3,

    OBS_ATLEAST3,

    "DIE convergence across at least three resistant genotypes",

    "Permutation_null_ATLEAST3_resistant.pdf"
)


plot_null(

    null_atleast2,

    OBS_ATLEAST2,

    "DIE convergence across at least two resistant genotypes",

    "Permutation_null_ATLEAST2_resistant.pdf"
)


# ============================================================
# 22. PRINT RESULTS
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "DIE CONVERGENCE PERMUTATION RESULTS\n"
)

cat(
    "============================================================\n\n"
)


print(
    permutation_summary,
    n = Inf,
    width = Inf
)


cat(
    "\nConvergent gene tables:\n\n"
)


cat(
    "All 4 resistant: ",
    nrow(
        all4
    ),
    "\n",
    sep = ""
)


cat(
    "All 4 resistant / zero susceptible: ",
    nrow(
        all4_0S
    ),
    "\n",
    sep = ""
)


cat(
    "At least 3 resistant: ",
    nrow(
        atleast3
    ),
    "\n",
    sep = ""
)


cat(
    "At least 3 resistant / zero susceptible: ",
    nrow(
        atleast3_0S
    ),
    "\n",
    sep = ""
)


cat(
    "\nOutput directory:\n",
    OUT_DIR,
    "\n"
)








### Assess the same question for DIU genes : do the same V. riparia genes undergo DIU across all resistant hosts? 
### Also test for overlap between shared/convergent DIU genes and shared DIE genes and perform permutation tests

#!/usr/bin/env Rscript


# ============================================================
# 1. Do the same shared V. riparia genes show significant
#    differential isoform usage (DIU) across resistant hosts?
#
# 2. Is the observed overlap greater than expected by chance,
#    after accounting for each genotype's own DIU-testable
#    gene universe and total number of significant DIU genes?
#
# 3. Which convergent DIU genes are absent from susceptible
#    hosts?
#
# 4. How does DIU convergence compare with differential
#    isoform expression (DIE) convergence?
#
#
# Resistant genotypes, V. riparia hybrids:
#
#   3309
#   SW
#   5BB
#   1616
#
# Susceptible V. riparia hybrids:
#
#   101-14
#   SO4
#
# ============================================================


suppressPackageStartupMessages({

    library(tidyverse)
})


# ============================================================
# 1. PATHS
# ============================================================

BASE <- "/rhome/wcoll010/shared/pacbio/IsoSeq_analyses_2"


DIU_BASE <- file.path(
    BASE,
    "DIU_IsoformSwitchAnalyzeR",
    "results"
)


CONV_DIR <- file.path(
    BASE,
    "Vrip_resistance_convergence"
)


DIE_DIR <- file.path(
    CONV_DIR,
    "DIE_gene_convergence"
)


DIE_TEST_DIR <- file.path(
    CONV_DIR,
    "DIE_gene_convergence_test"
)


OUT_DIR <- file.path(
    CONV_DIR,
    "DIU_gene_convergence_test"
)


dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


CONSENSUS_MAP_FILE <- file.path(
    CONV_DIR,
    "Vrip_gene_consensus_map.tsv"
)


# ============================================================
# 2. PARAMETERS
# ============================================================

N_PERMUTATIONS <- 100000L

SEED <- 12345


resistant <- c(
    "3309",
    "SW",
    "5BB",
    "1616"
)


susceptible <- c(
    "101-14",
    "SO4"
)


all_genotypes <- c(
    resistant,
    susceptible
)


# Directory names in DIU output

dir_names <- c(

    "3309" =
        "3309",

    "SW" =
        "sw",

    "5BB" =
        "5BB",

    "1616" =
        "1616",

    "101-14" =
        "101-14",

    "SO4" =
        "SO4"
)


# Reference used by each genotype

reference_names <- c(

    "3309" =
        "Vrip_X_Vrup",

    "SW" =
        "Vrip_X_Vrup",

    "101-14" =
        "Vrip_X_Vrup",

    "5BB" =
        "Vrip_X_Vber",

    "SO4" =
        "Vrip_X_Vber",

    "1616" =
        "Vrip_X_Vace"
)


# ============================================================
# 3. READ CONSENSUS MAP
# ============================================================

consensus_map <- read.delim(

    CONSENSUS_MAP_FILE,

    stringsAsFactors = FALSE,

    check.names = FALSE
)


required_map_cols <- c(
    "Reference",
    "Original_Gene_ID",
    "Consensus_Gene_ID"
)


missing_map_cols <- setdiff(
    required_map_cols,
    colnames(consensus_map)
)


if (
    length(missing_map_cols) > 0
) {

    stop(
        "Missing consensus-map columns:\n",
        paste(
            missing_map_cols,
            collapse = "\n"
        )
    )
}


gene_map <- consensus_map %>%

    dplyr::select(
        Reference,
        Original_Gene_ID,
        Consensus_Gene_ID
    ) %>%

    distinct()


# ============================================================
# 4. HELPER: PATHS FOR ONE GENOTYPE
# ============================================================

get_paths <- function(
    genotype
) {

    d <- dir_names[[
        genotype
    ]]


    genotype_dir <- file.path(
        DIU_BASE,
        d
    )


    list(

        significant_genes =
            file.path(
                genotype_dir,
                paste0(
                    d,
                    "_significant_switching_genes.csv"
                )
            ),

        significant_isoforms =
            file.path(
                genotype_dir,
                paste0(
                    d,
                    "_significant_switching_isoforms.csv"
                )
            ),

        full_features =
            file.path(
                genotype_dir,
                paste0(
                    d,
                    "_DIU_isoform_features_full.csv"
                )
            )
    )
}


# ============================================================
# 5. CHECK INPUTS
# ============================================================

for (
    g in all_genotypes
) {

    p <- get_paths(
        g
    )


    files <- unlist(
        p
    )


    missing <- files[
        !file.exists(
            files
        )
    ]


    if (
        length(missing) > 0
    ) {

        stop(
            "\nMissing DIU files for ",
            g,
            ":\n",
            paste(
                missing,
                collapse = "\n"
            ),
            "\n"
        )
    }
}


# ============================================================
# 6. MAP ONE LOCAL GENE SET TO CONSENSUS IDS
# ============================================================

map_local_genes <- function(
    genes,
    genotype
) {

    reference <- reference_names[[
        genotype
    ]]


    local_map <- gene_map %>%

        filter(
            Reference ==
                reference
        )


    mapped <- tibble(
        Original_Gene_ID =
            unique(
                genes
            )
    ) %>%

        filter(
            !is.na(
                Original_Gene_ID
            ),
            Original_Gene_ID != ""
        ) %>%

        inner_join(
            local_map,
            by =
                "Original_Gene_ID"
        ) %>%

        distinct(
            Consensus_Gene_ID
        )


    mapped$Consensus_Gene_ID
}


# ============================================================
# 7. READ SIGNIFICANT DIU GENES FOR EACH GENOTYPE
# ============================================================

sig_DIU_sets <- list()


sig_gene_rows <- list()


cat(
    "\n============================================================\n"
)

cat(
    "SIGNIFICANT DIU GENES\n"
)

cat(
    "============================================================\n\n"
)


for (
    g in all_genotypes
) {

    p <- get_paths(
        g
    )


    sig <- read.csv(

        p$significant_genes,

        stringsAsFactors = FALSE,

        check.names = FALSE
    )


    if (
        !"gene_id" %in%
        colnames(sig)
    ) {

        stop(
            g,
            ": gene_id column missing from significant gene file."
        )
    }


    n_original <- n_distinct(
        sig$gene_id
    )


    mapped_ids <- map_local_genes(
        sig$gene_id,
        g
    )


    sig_DIU_sets[[
        g
    ]] <- unique(
        mapped_ids
    )


    mapped_table <- sig %>%

        rename(
            Original_Gene_ID =
                gene_id
        ) %>%

        inner_join(

            gene_map %>%
                filter(
                    Reference ==
                        reference_names[[g]]
                ),

            by =
                "Original_Gene_ID"
        ) %>%

        mutate(
            Genotype =
                g
        )


    sig_gene_rows[[
        g
    ]] <- mapped_table


    cat(
        g,
        ": ",
        n_original,
        " significant DIU genes; ",
        length(
            sig_DIU_sets[[g]]
        ),
        " mapped Vrip consensus genes\n",
        sep = ""
    )
}


all_sig_gene_rows <- bind_rows(
    sig_gene_rows
)


write.csv(

    all_sig_gene_rows,

    file.path(
        OUT_DIR,
        "All_genotypes_significant_DIU_genes_consensus_mapped.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 8. BUILD DIU-TESTABLE GENE UNIVERSES
#
# Use DIU_isoform_features_full.csv.
#
# Every gene represented here entered the IsoformSwitchAnalyzeR
# analysis after expression filtering / multi-isoform filtering.
# ============================================================

testable_sets <- list()


cat(
    "\n============================================================\n"
)

cat(
    "DIU-TESTABLE GENE UNIVERSES\n"
)

cat(
    "============================================================\n\n"
)


for (
    g in resistant
) {

    p <- get_paths(
        g
    )


    full <- read.csv(

        p$full_features,

        stringsAsFactors = FALSE,

        check.names = FALSE
    )


    if (
        !"gene_id" %in%
        colnames(full)
    ) {

        stop(
            g,
            ": gene_id column missing from full DIU feature table."
        )
    }


    local_genes <- unique(
        full$gene_id
    )


    mapped <- map_local_genes(
        local_genes,
        g
    )


    testable_sets[[
        g
    ]] <- unique(
        mapped
    )


    cat(
        g,
        ": ",
        length(
            testable_sets[[g]]
        ),
        " DIU-testable Vrip consensus genes\n",
        sep = ""
    )
}


# ============================================================
# 9. QC: SIGNIFICANT GENES MUST BE TESTABLE
# ============================================================

for (
    g in resistant
) {

    missing_from_universe <- setdiff(

        sig_DIU_sets[[g]],

        testable_sets[[g]]
    )


    if (
        length(
            missing_from_universe
        ) > 0
    ) {

        stop(
            "\n",
            g,
            ": ",
            length(
                missing_from_universe
            ),
            " significant DIU genes are absent from the ",
            "DIU-testable universe.\n"
        )
    }
}


# ============================================================
# 10. TESTABLE-UNIVERSE SUMMARY
# ============================================================

testability_summary <- tibble(

    Genotype =
        resistant,

    DIU_Testable_Vrip_Genes =
        vapply(
            resistant,
            function(g) {
                length(
                    testable_sets[[g]]
                )
            },
            integer(1)
        ),

    Significant_DIU_Genes =
        vapply(
            resistant,
            function(g) {
                length(
                    sig_DIU_sets[[g]]
                )
            },
            integer(1)
        )
) %>%

    mutate(

        Percent_Testable_Genes_DIU =
            100 *
            Significant_DIU_Genes /
            DIU_Testable_Vrip_Genes
    )


write.csv(

    testability_summary,

    file.path(
        OUT_DIR,
        "DIU_genotype_testable_universe_summary.csv"
    ),

    row.names = FALSE
)


cat(
    "\n"
)


print(
    testability_summary,
    n = Inf,
    width = Inf
)


# ============================================================
# 11. RESISTANT DIU OVERLAP MATRIX
# ============================================================

all_resistant_DIU_genes <- sort(

    unique(
        unlist(
            sig_DIU_sets[
                resistant
            ]
        )
    )
)


resistant_matrix <- tibble(

    Consensus_Gene_ID =
        all_resistant_DIU_genes
)


for (
    g in resistant
) {

    resistant_matrix[[
        g
    ]] <-
        resistant_matrix$Consensus_Gene_ID %in%
        sig_DIU_sets[[g]]
}


resistant_matrix <- resistant_matrix %>%

    rowwise() %>%

    mutate(

        N_Resistant_With_DIU =
            sum(
                c_across(
                    all_of(
                        resistant
                    )
                )
            )
    ) %>%

    ungroup() %>%

    arrange(
        desc(
            N_Resistant_With_DIU
        ),
        Consensus_Gene_ID
    )


write.csv(

    resistant_matrix,

    file.path(
        OUT_DIR,
        "Resistant_DIU_gene_overlap_matrix.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 12. ADD SUSCEPTIBLE DIU STATUS
# ============================================================

combined <- resistant_matrix


for (
    g in susceptible
) {

    col <- paste0(
        "Susceptible_",
        g
    )


    combined[[
        col
    ]] <-
        combined$Consensus_Gene_ID %in%
        sig_DIU_sets[[g]]
}


combined <- combined %>%

    rowwise() %>%

    mutate(

        N_Susceptible_With_DIU =
            sum(
                c_across(
                    all_of(
                        paste0(
                            "Susceptible_",
                            susceptible
                        )
                    )
                )
            ),

        Resistant_Specific_All4 =
            N_Resistant_With_DIU == 4 &
            N_Susceptible_With_DIU == 0,

        Resistant_Specific_AtLeast3 =
            N_Resistant_With_DIU >= 3 &
            N_Susceptible_With_DIU == 0
    ) %>%

    ungroup()


write.csv(

    combined,

    file.path(
        OUT_DIR,
        "DIU_resistant_vs_susceptible_gene_convergence.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 13. OBSERVED OVERLAP SUMMARY
# ============================================================

overlap_summary <- combined %>%

    count(
        N_Resistant_With_DIU,
        name =
            "N_Genes"
    ) %>%

    arrange(
        desc(
            N_Resistant_With_DIU
        )
    )


write.csv(

    overlap_summary,

    file.path(
        OUT_DIR,
        "Resistant_DIU_gene_overlap_summary.csv"
    ),

    row.names = FALSE
)


OBS_ALL4 <- sum(
    combined$N_Resistant_With_DIU ==
        4
)


OBS_ATLEAST3 <- sum(
    combined$N_Resistant_With_DIU >=
        3
)


OBS_ATLEAST2 <- sum(
    combined$N_Resistant_With_DIU >=
        2
)


cat(
    "\nObserved DIU overlap:\n"
)


cat(
    "All 4 resistant: ",
    OBS_ALL4,
    "\n",
    sep = ""
)


cat(
    "At least 3 resistant: ",
    OBS_ATLEAST3,
    "\n",
    sep = ""
)


cat(
    "At least 2 resistant: ",
    OBS_ATLEAST2,
    "\n",
    sep = ""
)


# ============================================================
# 14. PREPARE PERMUTATION UNIVERSE
#
# Each genotype retains:
#
#   its own DIU-testable genes
#   its own number of significant DIU genes
# ============================================================

master_universe <- sort(

    unique(
        unlist(
            testable_sets
        )
    )
)


gene_index <- seq_along(
    master_universe
)


names(
    gene_index
) <- master_universe


testable_indices <- lapply(

    testable_sets,

    function(x) {

        unname(
            gene_index[
                x
            ]
        )
    }
)


observed_sizes <- vapply(

    resistant,

    function(g) {

        length(
            sig_DIU_sets[[g]]
        )
    },

    integer(1)
)


names(
    observed_sizes
) <- resistant


# ============================================================
# 15. PERMUTATION
# ============================================================

set.seed(
    SEED
)


N_MASTER <- length(
    master_universe
)


null_all4 <- integer(
    N_PERMUTATIONS
)


null_atleast3 <- integer(
    N_PERMUTATIONS
)


null_atleast2 <- integer(
    N_PERMUTATIONS
)


cat(
    "\nRunning ",
    format(
        N_PERMUTATIONS,
        big.mark = ","
    ),
    " DIU permutations...\n",
    sep = ""
)


for (
    i in seq_len(
        N_PERMUTATIONS
    )
) {

    counts <- integer(
        N_MASTER
    )


    for (
        g in resistant
    ) {

        sampled <- sample(

            testable_indices[[g]],

            size =
                observed_sizes[[g]],

            replace =
                FALSE
        )


        counts[
            sampled
        ] <-
            counts[
                sampled
            ] +
            1L
    }


    null_all4[i] <-
        sum(
            counts ==
                4L
        )


    null_atleast3[i] <-
        sum(
            counts >=
                3L
        )


    null_atleast2[i] <-
        sum(
            counts >=
                2L
        )
}


# ============================================================
# 16. PERMUTATION SUMMARY
# ============================================================

empirical_p <- function(
    observed,
    null
) {

    (
        sum(
            null >=
                observed
        ) +
        1
    ) /
    (
        length(
            null
        ) +
        1
    )
}


summarize_null <- function(
    criterion,
    observed,
    null
) {

    null_mean <- mean(
        null
    )


    tibble(

        Criterion =
            criterion,

        Observed =
            observed,

        Null_Mean =
            null_mean,

        Null_Median =
            median(
                null
            ),

        Null_SD =
            sd(
                null
            ),

        Null_2.5pct =
            as.numeric(
                quantile(
                    null,
                    0.025
                )
            ),

        Null_97.5pct =
            as.numeric(
                quantile(
                    null,
                    0.975
                )
            ),

        Fold_Over_Null =
            ifelse(
                null_mean > 0,
                observed /
                    null_mean,
                NA_real_
            ),

        Empirical_P =
            empirical_p(
                observed,
                null
            )
    )
}


perm_summary <- bind_rows(

    summarize_null(
        "All_4_resistant",
        OBS_ALL4,
        null_all4
    ),

    summarize_null(
        "At_least_3_of_4_resistant",
        OBS_ATLEAST3,
        null_atleast3
    ),

    summarize_null(
        "At_least_2_of_4_resistant",
        OBS_ATLEAST2,
        null_atleast2
    )
)


write.csv(

    perm_summary,

    file.path(
        OUT_DIR,
        "DIU_gene_convergence_permutation_summary.csv"
    ),

    row.names = FALSE
)


write.csv(

    tibble(

        Permutation =
            seq_len(
                N_PERMUTATIONS
            ),

        Shared_All4 =
            null_all4,

        Shared_AtLeast3 =
            null_atleast3,

        Shared_AtLeast2 =
            null_atleast2
    ),

    file.path(
        OUT_DIR,
        "DIU_gene_convergence_null_distributions.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 17. READ SIGNIFICANT DIU ISOFORM DETAILS
# ============================================================

isoform_details <- list()


for (
    g in all_genotypes
) {

    p <- get_paths(
        g
    )


    x <- read.csv(

        p$significant_isoforms,

        stringsAsFactors = FALSE,

        check.names = FALSE
    )


    required <- c(
        "isoform_id",
        "gene_id",
        "IF1",
        "IF2",
        "dIF",
        "isoform_switch_q_value"
    )


    missing <- setdiff(
        required,
        colnames(x)
    )


    if (
        length(missing) > 0
    ) {

        stop(
            g,
            ": missing DIU isoform columns:\n",
            paste(
                missing,
                collapse = "\n"
            )
        )
    }


    local_map <- gene_map %>%

        filter(
            Reference ==
                reference_names[[g]]
        )


    x <- x %>%

        rename(
            Original_Gene_ID =
                gene_id
        ) %>%

        inner_join(

            local_map,

            by =
                "Original_Gene_ID"
        ) %>%

        mutate(
            Genotype =
                g
        )


    isoform_details[[
        g
    ]] <- x
}


all_DIU_isoform_details <- bind_rows(
    isoform_details
)


write.csv(

    all_DIU_isoform_details,

    file.path(
        OUT_DIR,
        "All_genotypes_significant_DIU_isoforms_consensus_mapped.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 18. COLLAPSE DIU DETAILS TO GENE x GENOTYPE
# ============================================================

DIU_gene_detail <- all_DIU_isoform_details %>%

    group_by(
        Consensus_Gene_ID,
        Genotype
    ) %>%

    summarise(

        N_Significant_DIU_Isoforms =
            n_distinct(
                isoform_id
            ),

        Max_Abs_dIF =
            max(
                abs(
                    dIF
                ),
                na.rm = TRUE
            ),

        Positive_dIF_Isoforms =
            sum(
                dIF > 0,
                na.rm = TRUE
            ),

        Negative_dIF_Isoforms =
            sum(
                dIF < 0,
                na.rm = TRUE
            ),

        Significant_DIU_Isoforms =
            paste(
                sort(
                    unique(
                        isoform_id
                    )
                ),
                collapse = ";"
            ),

        dIF_Values =
            paste(
                signif(
                    dIF,
                    4
                ),
                collapse = ";"
            ),

        .groups =
            "drop"
    )


DIU_detail_wide <- DIU_gene_detail %>%

    pivot_wider(

        names_from =
            Genotype,

        values_from =
            c(
                N_Significant_DIU_Isoforms,
                Max_Abs_dIF,
                Positive_dIF_Isoforms,
                Negative_dIF_Isoforms,
                Significant_DIU_Isoforms,
                dIF_Values
            ),

        names_sep =
            "_"
    )


# ============================================================
# 19. WRITE CONVERGENT DIU GENE TABLES
# ============================================================

DIU_annotated <- combined %>%

    left_join(
        DIU_detail_wide,
        by =
            "Consensus_Gene_ID"
    )


# If the annotated DIE table exists, borrow its functional
# annotation columns for the same consensus genes.

ANNOTATED_DIE_FILE <- file.path(
    DIE_TEST_DIR,
    "All_DIE_convergence_genes_annotated.csv"
)


if (
    file.exists(
        ANNOTATED_DIE_FILE
    )
) {

    annotations <- read.csv(

        ANNOTATED_DIE_FILE,

        stringsAsFactors = FALSE,

        check.names = FALSE
    ) %>%

        dplyr::select(
            any_of(
                c(
                    "Consensus_Gene_ID",
                    "Preferred_Name",
                    "Description",
                    "GO_Terms",
                    "KEGG_Pathways",
                    "PFAMs"
                )
            )
        ) %>%

        distinct(
            Consensus_Gene_ID,
            .keep_all = TRUE
        )


    DIU_annotated <- DIU_annotated %>%

        left_join(
            annotations,
            by =
                "Consensus_Gene_ID"
        )
}


write.csv(

    DIU_annotated,

    file.path(
        OUT_DIR,
        "All_DIU_convergence_genes_annotated.csv"
    ),

    row.names = FALSE
)


write.csv(

    DIU_annotated %>%

        filter(
            N_Resistant_With_DIU ==
                4
        ),

    file.path(
        OUT_DIR,
        "Convergent_DIU_genes_ALL4_resistant_annotated.csv"
    ),

    row.names = FALSE
)


write.csv(

    DIU_annotated %>%

        filter(
            N_Resistant_With_DIU ==
                4,
            N_Susceptible_With_DIU ==
                0
        ),

    file.path(
        OUT_DIR,
        "Convergent_DIU_genes_ALL4_resistant_ZERO_susceptible_annotated.csv"
    ),

    row.names = FALSE
)


write.csv(

    DIU_annotated %>%

        filter(
            N_Resistant_With_DIU >=
                3
        ),

    file.path(
        OUT_DIR,
        "Convergent_DIU_genes_ATLEAST3_resistant_annotated.csv"
    ),

    row.names = FALSE
)


write.csv(

    DIU_annotated %>%

        filter(
            N_Resistant_With_DIU >=
                3,
            N_Susceptible_With_DIU ==
                0
        ),

    file.path(
        OUT_DIR,
        "Convergent_DIU_genes_ATLEAST3_resistant_ZERO_susceptible_annotated.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 20. DIE vs DIU — PER GENOTYPE
# ============================================================

DIE_GENE_FILE <- file.path(
    DIE_DIR,
    "DIE_gene_by_genotype.csv"
)


if (
    !file.exists(
        DIE_GENE_FILE
    )
) {

    stop(
        "\nMissing DIE gene-by-genotype file:\n",
        DIE_GENE_FILE
    )
}


DIE_gene_genotype <- read.csv(

    DIE_GENE_FILE,

    stringsAsFactors = FALSE,

    check.names = FALSE
)


DIE_sets <- list()


for (
    g in all_genotypes
) {

    DIE_sets[[
        g
    ]] <- DIE_gene_genotype %>%

        filter(
            Genotype ==
                g
        ) %>%

        pull(
            Consensus_Gene_ID
        ) %>%

        unique()
}


comparison_by_genotype <- list()


for (
    g in all_genotypes
) {

    union_genes <- union(
        DIE_sets[[g]],
        sig_DIU_sets[[g]]
    )


    x <- tibble(

        Genotype =
            g,

        Consensus_Gene_ID =
            union_genes,

        DIE =
            union_genes %in%
            DIE_sets[[g]],

        DIU =
            union_genes %in%
            sig_DIU_sets[[g]]
    ) %>%

        mutate(

            Category =
                case_when(

                    DIE & DIU ~
                        "DIE_and_DIU",

                    DIE ~
                        "DIE_only",

                    DIU ~
                        "DIU_only",

                    TRUE ~
                        "Neither"
                )
        )


    comparison_by_genotype[[
        g
    ]] <- x
}


DIE_DIU_genotype_table <- bind_rows(
    comparison_by_genotype
)


write.csv(

    DIE_DIU_genotype_table,

    file.path(
        OUT_DIR,
        "DIE_vs_DIU_gene_status_by_genotype.csv"
    ),

    row.names = FALSE
)


DIE_DIU_genotype_summary <- DIE_DIU_genotype_table %>%

    count(
        Genotype,
        Category,
        name =
            "N_Genes"
    )


write.csv(

    DIE_DIU_genotype_summary,

    file.path(
        OUT_DIR,
        "DIE_vs_DIU_gene_status_by_genotype_summary.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 21. DIE vs DIU — CONVERGENCE LEVEL
# ============================================================

DIE_CONVERGENCE_FILE <- file.path(
    DIE_DIR,
    "DIE_resistant_vs_susceptible_gene_convergence.csv"
)


DIE_conv <- read.csv(

    DIE_CONVERGENCE_FILE,

    stringsAsFactors = FALSE,

    check.names = FALSE
)


DIU_conv <- combined


conv_union <- union(

    DIE_conv$Consensus_Gene_ID,

    DIU_conv$Consensus_Gene_ID
)


conv_compare <- tibble(

    Consensus_Gene_ID =
        conv_union
) %>%

    left_join(

        DIE_conv %>%

            dplyr::select(

                Consensus_Gene_ID,

                N_Resistant_With_DIE,

                N_Susceptible_With_DIE
            ),

        by =
            "Consensus_Gene_ID"
    ) %>%

    left_join(

        DIU_conv %>%

            dplyr::select(

                Consensus_Gene_ID,

                N_Resistant_With_DIU,

                N_Susceptible_With_DIU
            ),

        by =
            "Consensus_Gene_ID"
    ) %>%

    mutate(

        across(
            c(
                N_Resistant_With_DIE,
                N_Susceptible_With_DIE,
                N_Resistant_With_DIU,
                N_Susceptible_With_DIU
            ),
            ~ replace_na(
                .x,
                0
            )
        ),

        DIE_All4 =
            N_Resistant_With_DIE ==
            4,

        DIU_All4 =
            N_Resistant_With_DIU ==
            4,

        DIE_AtLeast3 =
            N_Resistant_With_DIE >=
            3,

        DIU_AtLeast3 =
            N_Resistant_With_DIU >=
            3,

        DIE_All4_0S =
            DIE_All4 &
            N_Susceptible_With_DIE ==
            0,

        DIU_All4_0S =
            DIU_All4 &
            N_Susceptible_With_DIU ==
            0,

        DIE_DIU_Convergence_Class =
            case_when(

                DIE_All4 &
                DIU_All4 ~
                    "All4_DIE_and_DIU",

                DIE_All4 ~
                    "All4_DIE_only",

                DIU_All4 ~
                    "All4_DIU_only",

                DIE_AtLeast3 &
                DIU_AtLeast3 ~
                    "AtLeast3_DIE_and_DIU",

                DIE_AtLeast3 ~
                    "AtLeast3_DIE_only",

                DIU_AtLeast3 ~
                    "AtLeast3_DIU_only",

                TRUE ~
                    "Lower_overlap"
            )
    )


# Add annotation if available

if (
    exists(
        "annotations"
    )
) {

    conv_compare <- conv_compare %>%

        left_join(
            annotations,
            by =
                "Consensus_Gene_ID"
        )
}


write.csv(

    conv_compare,

    file.path(
        OUT_DIR,
        "DIE_vs_DIU_convergence_gene_table.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 22. CONVERGENCE SUMMARY
# ============================================================

convergence_summary <- tibble(

    Metric = c(

        "DIE genes shared by all 4 resistant",

        "DIU genes shared by all 4 resistant",

        "Genes shared by all 4 in BOTH DIE and DIU",

        "DIE genes shared by >=3 resistant",

        "DIU genes shared by >=3 resistant",

        "Genes shared by >=3 in BOTH DIE and DIU",

        "DIE all4 / zero susceptible",

        "DIU all4 / zero susceptible",

        "All4/zero-susceptible in BOTH DIE and DIU"
    ),

    N_Genes = c(

        sum(
            conv_compare$DIE_All4
        ),

        sum(
            conv_compare$DIU_All4
        ),

        sum(
            conv_compare$DIE_All4 &
            conv_compare$DIU_All4
        ),

        sum(
            conv_compare$DIE_AtLeast3
        ),

        sum(
            conv_compare$DIU_AtLeast3
        ),

        sum(
            conv_compare$DIE_AtLeast3 &
            conv_compare$DIU_AtLeast3
        ),

        sum(
            conv_compare$DIE_All4_0S
        ),

        sum(
            conv_compare$DIU_All4_0S
        ),

        sum(
            conv_compare$DIE_All4_0S &
            conv_compare$DIU_All4_0S
        )
    )
)


write.csv(

    convergence_summary,

    file.path(
        OUT_DIR,
        "DIE_vs_DIU_convergence_summary.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 23. FOCAL IMMUNE GENES
# ============================================================

focal_genes <- c(

    EDS1 =
        "VITVri588271_v1.0.hap1.chr17.ver1.0.g351240",

    SAG101 =
        "VITVri588271_v1.0.hap1.chr14.ver1.0.g309600",

    NIMIN2 =
        "VITVri588271_v1.0.hap1.chr01.ver1.0.g071050"
)


focal_table <- conv_compare %>%

    filter(
        Consensus_Gene_ID %in%
            focal_genes
    ) %>%

    mutate(

        Focal_Gene =
            names(
                focal_genes
            )[
                match(
                    Consensus_Gene_ID,
                    focal_genes
                )
            ]
    ) %>%

    relocate(
        Focal_Gene
    )


write.csv(

    focal_table,

    file.path(
        OUT_DIR,
        "EDS1_SAG101_NIMIN2_DIE_vs_DIU.csv"
    ),

    row.names = FALSE
)


# ============================================================
# 24. PERMUTATION PLOTS
# ============================================================

plot_null <- function(
    null,
    observed,
    title,
    filename
) {

    p <- ggplot(

        tibble(
            Null =
                null
        ),

        aes(
            x =
                Null
        )
    ) +

        geom_histogram(
            binwidth = 1,
            boundary = -0.5
        ) +

        geom_vline(
            xintercept =
                observed,
            linewidth =
                1
        ) +

        labs(

            title =
                title,

            subtitle =
                paste0(
                    "Observed = ",
                    observed
                ),

            x =
                "Overlapping genes under null",

            y =
                "Permutations"
        ) +

        theme_bw(
            base_size =
                12
        )


    ggsave(

        file.path(
            OUT_DIR,
            filename
        ),

        p,

        width = 7,
        height = 5
    )
}


plot_null(

    null_all4,

    OBS_ALL4,

    "DIU convergence across all four resistant genotypes",

    "Permutation_null_DIU_ALL4.pdf"
)


plot_null(

    null_atleast3,

    OBS_ATLEAST3,

    "DIU convergence across at least three resistant genotypes",

    "Permutation_null_DIU_ATLEAST3.pdf"
)


plot_null(

    null_atleast2,

    OBS_ATLEAST2,

    "DIU convergence across at least two resistant genotypes",

    "Permutation_null_DIU_ATLEAST2.pdf"
)


# ============================================================
# 25. PRINT FINAL RESULTS
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "DIU CONVERGENCE PERMUTATION RESULTS\n"
)

cat(
    "============================================================\n\n"
)


print(
    perm_summary,
    n = Inf,
    width = Inf
)


cat(
    "\n============================================================\n"
)

cat(
    "DIE VS DIU CONVERGENCE\n"
)

cat(
    "============================================================\n\n"
)


print(
    convergence_summary,
    n = Inf,
    width = Inf
)


cat(
    "\n============================================================\n"
)

cat(
    "FOCAL IMMUNE GENES\n"
)

cat(
    "============================================================\n\n"
)


print(
    focal_table,
    n = Inf,
    width = Inf
)


cat(
    "\nOutput directory:\n",
    OUT_DIR,
    "\n"
)




