### Create a gene map so that V. riparia genes expressed across different hybrid genotypes (sqanti3 catalogs / reference transcriptomes) are harmonized /equivilent 
### Novel genes are given arbitrary names , e.g. novelGene_1 , not mapped to a known, annotated gene
### Novel genes must be assigned new names based on chromosome,strand, and overlap (>= 80%)

#!/usr/bin/env Rscript


# ============================================================
# 1. LIBRARIES
# ============================================================

my_lib <- "/rhome/wcoll010/bigdata_groenlab/R_libs"

.libPaths(
    c(
        my_lib,
        .libPaths()
    )
)


suppressPackageStartupMessages({

    library(tidyverse)
    library(GenomicRanges)
    library(rtracklayer)
    library(igraph)
})


# ============================================================
# 2. PATHS
# ============================================================

BASE <- "/rhome/wcoll010/shared/pacbio/IsoSeq_analyses_2"


OUT_DIR <- file.path(
    BASE,
    "Vrip_resistance_convergence"
)


dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


OVERLAP_THRESHOLD <- 0.80


# ============================================================
# 3. INPUT DATASETS
# ============================================================

datasets <- list(

    Vrip_X_Vrup = list(

        classification = file.path(
            BASE,
            "sqanti_results",
            "Vrip_X_Vrup_sqanti",
            "filter_results",
            "Vrip_X_Vrup_filtered_RulesFilter_result_classification.txt"
        ),

        gtf = file.path(
            BASE,
            "sqanti_results",
            "Vrip_X_Vrup_sqanti",
            "filter_results",
            "Vrip_X_Vrup_filtered.filtered.gtf"
        )
    ),


    Vrip_X_Vber = list(

        classification = file.path(
            BASE,
            "sqanti_results",
            "Vrip_X_Vber_sqanti",
            "filter_results",
            "Vrip_X_Vber_filtered_RulesFilter_result_classification.txt"
        ),

        gtf = file.path(
            BASE,
            "sqanti_results",
            "Vrip_X_Vber_sqanti",
            "filter_results",
            "Vrip_X_Vber_filtered.filtered.gtf"
        )
    ),


    Vrip_X_Vace = list(

        classification = file.path(
            BASE,
            "sqanti_results",
            "Vrip_X_Vace_sqanti",
            "filter_results",
            "Vrip_X_Vace_filtered_RulesFilter_result_classification.txt"
        ),

        gtf = file.path(
            BASE,
            "sqanti_results",
            "Vrip_X_Vace_sqanti",
            "filter_results",
            "Vrip_X_Vace_filtered.filtered.gtf"
        )
    )
)


# ============================================================
# 4. V. RIPARIA CHROMOSOME IDENTIFIER
# Restrict to the common Vrip subgenome.
# ============================================================

VRIP_CHR_PATTERN <- "^VITVri588271_v1\\.0\\.hap1\\.chr"


# ============================================================
# 5. EXTRACT TRANSCRIPT COORDINATES FROM GTF

# Uses transcript records if available; if the GTF does not contain explicit transcript features, exon coordinates are collapsed instead.
# ============================================================

get_transcript_coords <- function(
    gtf_file
) {

    message(
        "  Reading GTF:"
    )

    message(
        "    ",
        gtf_file
    )


    gtf <- import(
        gtf_file
    )


    gtf_df <- as.data.frame(
        gtf
    )


    if (
        !"transcript_id" %in%
        colnames(
            gtf_df
        )
    ) {

        stop(
            "GTF does not contain transcript_id:\n",
            gtf_file
        )
    }


    # --------------------------------------------------------
    # Preferred:
    # explicit transcript records
    # --------------------------------------------------------

    tx <- gtf_df %>%

        filter(
            type == "transcript",
            !is.na(
                transcript_id
            )
        )


    if (
        nrow(
            tx
        ) > 0
    ) {

        coords <- tx %>%

            transmute(

                isoform =
                    as.character(
                        transcript_id
                    ),

                chrom =
                    as.character(
                        seqnames
                    ),

                strand =
                    as.character(
                        strand
                    ),

                start =
                    as.numeric(
                        start
                    ),

                end =
                    as.numeric(
                        end
                    )
            )


    } else {


        # ----------------------------------------------------
        # Fallback:
        # derive transcript span from exon coordinates
        # ----------------------------------------------------

        message(
            "  No transcript features found; ",
            "deriving transcript spans from exons."
        )


        coords <- gtf_df %>%

            filter(
                type == "exon",
                !is.na(
                    transcript_id
                )
            ) %>%

            transmute(

                isoform =
                    as.character(
                        transcript_id
                    ),

                chrom =
                    as.character(
                        seqnames
                    ),

                strand =
                    as.character(
                        strand
                    ),

                start =
                    as.numeric(
                        start
                    ),

                end =
                    as.numeric(
                        end
                    )
            ) %>%

            group_by(
                isoform,
                chrom,
                strand
            ) %>%

            summarise(

                start =
                    min(
                        start,
                        na.rm = TRUE
                    ),

                end =
                    max(
                        end,
                        na.rm = TRUE
                    ),

                .groups = "drop"
            )
    }


    coords %>%

        distinct(
            isoform,
            .keep_all = TRUE
        )
}


# ============================================================
# 6. LOAD EACH SQANTI CATALOG
# ============================================================

load_reference <- function(
    reference
) {

    info <- datasets[[
        reference
    ]]


    message(
        "\n============================================================"
    )

    message(
        "REFERENCE: ",
        reference
    )

    message(
        "============================================================"
    )


    required <- c(
        info$classification,
        info$gtf
    )


    missing <- required[
        !file.exists(
            required
        )
    ]


    if (
        length(
            missing
        ) > 0
    ) {

        stop(
            "Missing file(s):\n",
            paste(
                missing,
                collapse = "\n"
            )
        )
    }


    classification <- read_tsv(

        info$classification,

        show_col_types = FALSE,

        progress = FALSE
    )


    required_cols <- c(
        "isoform",
        "associated_gene"
    )


    missing_cols <- setdiff(
        required_cols,
        colnames(
            classification
        )
    )


    if (
        length(
            missing_cols
        ) > 0
    ) {

        stop(
            reference,
            " classification missing column(s): ",
            paste(
                missing_cols,
                collapse = ", "
            )
        )
    }


    coords <- get_transcript_coords(
        info$gtf
    )


    dat <- classification %>%

        select(
            isoform,
            associated_gene
        ) %>%

        filter(

            !is.na(
                associated_gene
            ),

            associated_gene != "",

            associated_gene !=
                "SQANTI_Filtered_Artifact"
        ) %>%

        inner_join(
            coords,
            by = "isoform"
        ) %>%

        filter(
            grepl(
                VRIP_CHR_PATTERN,
                chrom
            )
        ) %>%

        mutate(
            Reference =
                reference
        )


    message(
        "  Vrip transcripts retained: ",
        nrow(
            dat
        )
    )


    message(
        "  Vrip associated genes: ",
        n_distinct(
            dat$associated_gene
        )
    )


    dat
}


# ============================================================
# 7. COMBINE THREE REFERENCES
# ============================================================

transcript_master <- bind_rows(

    lapply(
        names(
            datasets
        ),
        load_reference
    )
)


# ============================================================
# 8. COLLAPSE PB TRANSCRIPTS TO GENE LOCI
#
# One row per gene: reference, associated_gene, chromosome, strand
# ============================================================

gene_loci <- transcript_master %>%

    group_by(

        Reference,

        Original_Gene_ID =
            associated_gene,

        Chromosome =
            chrom,

        Strand =
            strand
    ) %>%

    summarise(

        Gene_Start =
            min(
                start,
                na.rm = TRUE
            ),

        Gene_End =
            max(
                end,
                na.rm = TRUE
            ),

        N_Transcripts =
            n_distinct(
                isoform
            ),

        .groups = "drop"
    ) %>%

    mutate(

        Gene_Width =
            Gene_End -
            Gene_Start +
            1,

        Is_Novel =
            grepl(
                "^novelGene",
                Original_Gene_ID
            ),

        Is_Annotated =
            grepl(
                "^VITVri",
                Original_Gene_ID
            )
    )


# ============================================================
# 9. QC: ONE ORIGINAL GENE SHOULD NOT FORM MULTIPLE LOCI
# ============================================================

duplicate_gene_loci <- gene_loci %>%

    count(
        Reference,
        Original_Gene_ID
    ) %>%

    filter(
        n > 1
    )


write_csv(

    duplicate_gene_loci,

    file.path(
        OUT_DIR,
        "QC_genes_with_multiple_genomic_loci.csv"
    )
)


if (
    nrow(
        duplicate_gene_loci
    ) > 0
) {

    warning(
        nrow(
            duplicate_gene_loci
        ),
        " Reference × gene IDs occur at >1 genomic locus. ",
        "Inspect QC_genes_with_multiple_genomic_loci.csv"
    )
}


# ============================================================
# 10. ANNOTATED GENES: Preserve exact VITVri ID.
# ============================================================

annotated_map <- gene_loci %>%

    filter(
        Is_Annotated
    ) %>%

    mutate(

        Consensus_Gene_ID =
            Original_Gene_ID,

        Consensus_Type =
            "Annotated_VITVri"
    )


# ============================================================
# 11. NOVEL LOCI
# ============================================================

novel_loci <- gene_loci %>%

    filter(
        Is_Novel
    ) %>%

    arrange(
        Chromosome,
        Gene_Start,
        Gene_End,
        Strand,
        Reference,
        Original_Gene_ID
    ) %>%

    mutate(
        Novel_Index =
            row_number()
    )


message(
    "\nNovel gene loci to harmonize: ",
    nrow(
        novel_loci
    )
)


# ============================================================
# 12. BUILD GRANGES FOR NOVEL GENE LOCI
# ============================================================

if (
    nrow(
        novel_loci
    ) > 0
) {

    novel_gr <- GRanges(

        seqnames =
            novel_loci$Chromosome,

        ranges =
            IRanges(
                start =
                    novel_loci$Gene_Start,
                end =
                    novel_loci$Gene_End
            ),

        strand =
            novel_loci$Strand
    )


    # --------------------------------------------------------
    # Candidate overlaps already require chromosome; ignore.strand = FALSE additionally requires strand.
    # --------------------------------------------------------

    hits <- findOverlaps(

        novel_gr,

        novel_gr,

        ignore.strand = FALSE
    )


    hit_df <- as.data.frame(
        hits
    ) %>%

        filter(
            queryHits < subjectHits
        )


    if (
        nrow(
            hit_df
        ) > 0
    ) {

        overlap_width <- width(

            pintersect(

                novel_gr[
                    hit_df$queryHits
                ],

                novel_gr[
                    hit_df$subjectHits
                ]
            )
        )


        hit_df <- hit_df %>%

            mutate(

                Overlap_bp =
                    overlap_width,

                Query_Width =
                    width(
                        novel_gr[
                            queryHits
                        ]
                    ),

                Subject_Width =
                    width(
                        novel_gr[
                            subjectHits
                        ]
                    ),

                Query_Fraction =
                    Overlap_bp /
                    Query_Width,

                Subject_Fraction =
                    Overlap_bp /
                    Subject_Width,

                Query_Reference =
                    novel_loci$Reference[
                        queryHits
                    ],

                Subject_Reference =
                    novel_loci$Reference[
                        subjectHits
                    ]
            ) %>%

            # NOTE: only use matches between DIFFERENT transcript catalogs to define a cross-reference consensus
            filter(
                Query_Reference !=
                    Subject_Reference
            )


        valid_edges <- hit_df %>%

            filter(

                Query_Fraction >=
                    OVERLAP_THRESHOLD,

                Subject_Fraction >=
                    OVERLAP_THRESHOLD
            )


    } else {

        valid_edges <- tibble(

            queryHits =
                integer(),

            subjectHits =
                integer()
        )
    }


    # ========================================================
    # 13. INTEGRATE DATA: Include every novel locus as a graph vertex so singleton genes receive a component as well.
    # ========================================================

    vertices <- tibble(
        name =
            as.character(
                seq_len(
                    nrow(
                        novel_loci
                    )
                )
            )
    )


    if (
        nrow(
            valid_edges
        ) > 0
    ) {

        edges <- valid_edges %>%

            transmute(

                from =
                    as.character(
                        queryHits
                    ),

                to =
                    as.character(
                        subjectHits
                    )
            )


        graph <- graph_from_data_frame(

            d =
                edges,

            directed =
                FALSE,

            vertices =
                vertices
        )


    } else {

        graph <- make_empty_graph(

            n =
                nrow(
                    novel_loci
                ),

            directed =
                FALSE
        )


        V(
            graph
        )$name <- vertices$name
    }


    component_membership <- components(
        graph
    )$membership


    novel_loci$Novel_Cluster <- as.integer(

        component_membership[
            as.character(
                seq_len(
                    nrow(
                        novel_loci
                    )
                )
            )
        ]
    )


    # ========================================================
    # 14. ASSIGN NOVEL IDs
    # ========================================================

    novel_cluster_info <- novel_loci %>%

        group_by(
            Novel_Cluster
        ) %>%

        summarise(

            Chromosome =
                dplyr::first(
                    Chromosome
                ),

            Strand =
                dplyr::first(
                    Strand
                ),

            Consensus_Start =
                min(
                    Gene_Start
                ),

            Consensus_End =
                max(
                    Gene_End
                ),

            N_References =
                n_distinct(
                    Reference
                ),

            References =
                paste(
                    sort(
                        unique(
                            Reference
                        )
                    ),
                    collapse = ";"
                ),

            N_Local_Genes =
                n(),

            .groups = "drop"
        ) %>%

        mutate(

            Chr_Short =
                str_extract(
                    Chromosome,
                    "chr[^.]+$"
                ),

            Strand_Label =
                case_when(

                    Strand == "+" ~
                        "plus",

                    Strand == "-" ~
                        "minus",

                    TRUE ~
                        "unknown"
                ),

            Consensus_Gene_ID =
                paste0(

                    "VripNovel_",

                    Chr_Short,

                    "_",

                    Strand_Label,

                    "_",

                    Consensus_Start,

                    "_",

                    Consensus_End
                ),

            Consensus_Type =
                case_when(

                    N_References == 3 ~
                        "Novel_shared_3",

                    N_References == 2 ~
                        "Novel_shared_2",

                    TRUE ~
                        "Novel_single_reference"
                )
        )


    novel_map <- novel_loci %>%

        left_join(

            novel_cluster_info %>%

                select(
                    Novel_Cluster,
                    Consensus_Gene_ID,
                    Consensus_Type,
                    N_References,
                    References
                ),

            by =
                "Novel_Cluster"
        )


} else {


    valid_edges <-
        tibble()


    novel_cluster_info <-
        tibble()


    novel_map <-
        novel_loci
}


# ============================================================
# 15. ADD REFERENCE COVERAGE FOR ANNOTATED GENES
# ============================================================

annotated_coverage <- annotated_map %>%

    group_by(
        Consensus_Gene_ID
    ) %>%

    summarise(

        N_References =
            n_distinct(
                Reference
            ),

        References =
            paste(
                sort(
                    unique(
                        Reference
                    )
                ),
                collapse = ";"
            ),

        .groups = "drop"
    )


annotated_map <- annotated_map %>%

    left_join(
        annotated_coverage,
        by =
            "Consensus_Gene_ID"
    )


# ============================================================
# 16. FINAL GENE CONSENSUS MAP
# ============================================================

annotated_final <- annotated_map %>%

    transmute(

        Reference,

        Original_Gene_ID,

        Consensus_Gene_ID,

        Consensus_Type,

        Chromosome,

        Strand,

        Gene_Start,

        Gene_End,

        Gene_Width,

        N_Transcripts,

        N_References,

        References
    )


novel_final <- novel_map %>%

    transmute(

        Reference,

        Original_Gene_ID,

        Consensus_Gene_ID,

        Consensus_Type,

        Chromosome,

        Strand,

        Gene_Start,

        Gene_End,

        Gene_Width,

        N_Transcripts,

        N_References,

        References
    )


consensus_map <- bind_rows(

    annotated_final,

    novel_final
) %>%

    arrange(
        Chromosome,
        Gene_Start,
        Reference,
        Original_Gene_ID
    )


# ============================================================
# 17. TRANSCRIPT-LEVEL CROSS-REFERENCE; Keeps the local PB ID, but adds the common gene ID.
# ============================================================

transcript_crosswalk <- transcript_master %>%

    transmute(

        Reference,

        Isoform_ID =
            isoform,

        Original_Gene_ID =
            associated_gene,

        Chromosome =
            chrom,

        Strand =
            strand,

        Transcript_Start =
            start,

        Transcript_End =
            end
    ) %>%

    left_join(

        consensus_map %>%

            select(

                Reference,

                Original_Gene_ID,

                Consensus_Gene_ID,

                Consensus_Type,

                N_References,

                References
            ),

        by = c(
            "Reference",
            "Original_Gene_ID"
        )
    )


if (
    any(
        is.na(
            transcript_crosswalk$Consensus_Gene_ID
        )
    )
) {

    stop(
        "ERROR: some Vrip transcripts failed to receive ",
        "a Consensus_Gene_ID."
    )
}


# ============================================================
# 18. CONSENSUS-LEVEL SUMMARY
# ============================================================

consensus_summary <- consensus_map %>%

    group_by(
        Consensus_Gene_ID,
        Consensus_Type
    ) %>%

    summarise(

        N_References =
            n_distinct(
                Reference
            ),

        References =
            paste(
                sort(
                    unique(
                        Reference
                    )
                ),
                collapse = ";"
            ),

        N_Local_Gene_IDs =
            n_distinct(
                paste(
                    Reference,
                    Original_Gene_ID
                )
            ),

        .groups = "drop"
    ) %>%

    arrange(
        Consensus_Type,
        desc(
            N_References
        ),
        Consensus_Gene_ID
    )


# ============================================================
# 19. WRITE OUTPUTS
# ============================================================

write_tsv(

    consensus_map,

    file.path(
        OUT_DIR,
        "Vrip_gene_consensus_map.tsv"
    )
)


write_tsv(

    transcript_crosswalk,

    file.path(
        OUT_DIR,
        "Vrip_transcript_consensus_crosswalk.tsv"
    )
)


write_tsv(

    consensus_summary,

    file.path(
        OUT_DIR,
        "Vrip_consensus_gene_summary.tsv"
    )
)


write_tsv(

    valid_edges,

    file.path(
        OUT_DIR,
        "Vrip_novel_gene_valid_overlap_edges.tsv"
    )
)


# ============================================================
# 20. SUMMARY
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "VRIP GENE HARMONIZATION COMPLETE\n"
)

cat(
    "============================================================\n"
)


cat(
    "\nTranscript crosswalk rows: ",
    nrow(
        transcript_crosswalk
    ),
    "\n",
    sep = ""
)


cat(
    "Local gene loci: ",
    nrow(
        consensus_map
    ),
    "\n",
    sep = ""
)


cat(
    "Unique consensus genes: ",
    n_distinct(
        consensus_map$Consensus_Gene_ID
    ),
    "\n\n",
    sep = ""
)


summary_counts <- consensus_summary %>%

    count(
        Consensus_Type,
        N_References,
        name = "N_Consensus_Genes"
    ) %>%

    arrange(
        Consensus_Type,
        N_References
    )


print(
    summary_counts,
    n = Inf,
    width = Inf
)


cat(
    "\nOutput directory:\n",
    OUT_DIR,
    "\n"
)


