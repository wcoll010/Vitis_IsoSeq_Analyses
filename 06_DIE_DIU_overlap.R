################################################################################
# Summarize DIE genes, DIU genes, and within-genotype gene-level overlap.
# DIE means a gene encodes >=1 significant DE isoform based on DIE analysis in 'Differential_isoform_expression_PCA.R'
################################################################################

#!/usr/bin/env Rscript


args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args)) args[1] else
  "/rhome/wcoll010/shared/pacbio/IsoSeq_analyses_2"
DE_DIR <- file.path(BASE, "DE_isoform_analysis")
DIU_DIR <- file.path(BASE, "DIU_IsoformSwitchAnalyzeR", "results")
OUT_DIR <- file.path(BASE, "DE_DIU_summary", "DIE_DIU_overlap")
genotypes <- c("101-14", "1103", "1616", "3309", "5BB",
               "99R", "cab33", "ripa", "SO4", "sw")

paths <- data.frame(
  Genotype = genotypes,
  DE_Up = file.path(DE_DIR, genotypes, "results",
                   paste0(genotypes, "_upregulated_isoforms.csv")),
  DE_Down = file.path(DE_DIR, genotypes, "results",
                     paste0(genotypes, "_downregulated_isoforms.csv")),
  DIU = file.path(DIU_DIR, genotypes,
                 paste0(genotypes, "_significant_switching_isoforms.csv")),
  stringsAsFactors = FALSE
)
required <- unlist(paths[c("DE_Up", "DE_Down", "DIU")], use.names = FALSE)
missing <- required[!file.exists(required)]
if (length(missing)) stop(
  "Missing input files (not counted as zero):\n",
  paste(missing, collapse = "\n"), call. = FALSE
)

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | tolower(x) %in%
      c("", "na", "nan", "n/a", "none", "null", "-")] <- NA_character_
  x
}
ids <- function(x) unique(x[!is.na(x)])

read_pairs <- function(path, iso_col, gene_col) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c(iso_col, gene_col)
  if (!all(needed %in% names(x))) stop(
    "Missing required columns in ", path, ": ",
    paste(setdiff(needed, names(x)), collapse = ", "), call. = FALSE
  )
  out <- unique(data.frame(
    isoform = clean_id(x[[iso_col]]),
    gene = clean_id(x[[gene_col]]),
    stringsAsFactors = FALSE
  ))
  if (anyNA(out$isoform)) stop("Missing isoform IDs in ", path, call. = FALSE)
  mapped <- out[!is.na(out$gene), , drop = FALSE]
  if (anyDuplicated(mapped$isoform)) stop(
    "An isoform maps to multiple gene IDs in ", path,
    "; resolve the mapping before summarizing.", call. = FALSE
  )
  out
}

rows <- qc <- overlap_rows <- vector("list", length(genotypes))
for (i in seq_along(genotypes)) {
  g <- genotypes[i]
  up <- read_pairs(paths$DE_Up[i], "isoform", "associated_gene")
  down <- read_pairs(paths$DE_Down[i], "isoform", "associated_gene")
  diu <- read_pairs(paths$DIU[i], "isoform_id", "gene_id")

  if (length(intersect(up$isoform, down$isoform))) stop(
    g, ": an isoform appears in BOTH up and down lists.", call. = FALSE
  )

  # Check consistency of mappings across the three input lists.
  mapped <- unique(rbind(up, down, diu))
  mapped <- mapped[!is.na(mapped$gene), , drop = FALSE]
  if (anyDuplicated(mapped$isoform)) stop(
    g, ": conflicting gene assignments across DIE and DIU lists.",
    " Resolve mappings before calculating overlap.", call. = FALSE
  )

  up_genes <- ids(up$gene)
  down_genes <- ids(down$gene)
  die_genes <- union(up_genes, down_genes)
  diu_genes <- ids(diu$gene)
  overlap_genes <- sort(intersect(die_genes, diu_genes))

  rows[[i]] <- data.frame(
    Genotype = g,
    DE_Isoforms_Up = length(ids(up$isoform)),
    DE_Isoforms_Down = length(ids(down$isoform)),
    DE_Isoforms_Total = length(union(up$isoform, down$isoform)),
    DE_Genes_Up = length(up_genes),
    DE_Genes_Down = length(down_genes),
    DE_Genes_Total = length(die_genes),
    DIU_Isoforms_Total = length(ids(diu$isoform)),
    DIU_Genes_Total = length(diu_genes),
    DIE_DIU_Overlap_Genes = length(overlap_genes)
  )
  # rep() preserves a valid, zero-row table when there is no overlap.
  overlap_rows[[i]] <- data.frame(
    Genotype = rep(g, length(overlap_genes)),
    gene_id = overlap_genes,
    stringsAsFactors = FALSE
  )
  qc[[i]] <- data.frame(
    Genotype = g,
    DE_Genes_With_Both_Up_And_Down_Isoforms =
      length(intersect(up_genes, down_genes)),
    DE_Up_Isoforms_Without_Gene =
      length(setdiff(up$isoform, up$isoform[!is.na(up$gene)])),
    DE_Down_Isoforms_Without_Gene =
      length(setdiff(down$isoform, down$isoform[!is.na(down$gene)])),
    DIU_Isoforms_Without_Gene =
      length(setdiff(diu$isoform, diu$isoform[!is.na(diu$gene)]))
  )
}

summary_table <- do.call(rbind, rows)
overlap_table <- do.call(rbind, overlap_rows)
qc_table <- do.call(rbind, qc)
if (any(as.matrix(qc_table[3:5]) > 0)) warning(
  "Some isoforms lack gene IDs: counted as isoforms but excluded ",
  "from gene totals and overlap. See QC output."
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(summary_table,
          file.path(OUT_DIR, "DIE_DIU_overlap_summary.csv"), row.names = FALSE)
write.table(summary_table,
            file.path(OUT_DIR, "DIE_DIU_overlap_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(overlap_table,
          file.path(OUT_DIR, "DIE_DIU_overlapping_gene_ids.csv"),
          row.names = FALSE)
write.csv(qc_table, file.path(OUT_DIR, "DIE_DIU_overlap_QC.csv"),
          row.names = FALSE)
write.csv(paths, file.path(OUT_DIR, "DIE_DIU_overlap_input_paths.csv"),
          row.names = FALSE)
print(summary_table, row.names = FALSE)
cat("\nOutputs saved to:", OUT_DIR, "\n")














### Test whether overlap between DIE genes and DIU genes is statistically significant based on Fisher exact tests
### P values are adjusted across genotypes by Benjamini-Hochberg.

#   DIE: >=1 isoform with finite padj and log2FoldChange
#   DIU: >=2 isoforms in full results, and >=1 isoform with significant DIU, based on analysis in 'Differential_isoform_usage.txt'
#   Gene universe/background: intersection of these gene sets (DIE and DIU genes)



args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args)) args[1] else
  "/rhome/wcoll010/shared/pacbio/IsoSeq_analyses_2"
OUT <- file.path(BASE, "DE_DIU_summary", "DIE_DIU_Fisher")
genotypes <- c("101-14", "1103", "1616", "3309", "5BB",
               "99R", "cab33", "ripa", "SO4", "sw")

paths <- data.frame(
  Genotype = genotypes,
  DIE_Full = file.path(BASE, "DE_isoform_analysis", genotypes, "results",
                      paste0(genotypes, "_all_isoform_DE_results.csv")),
  DIE_Up = file.path(BASE, "DE_isoform_analysis", genotypes, "results",
                    paste0(genotypes, "_upregulated_isoforms.csv")),
  DIE_Down = file.path(BASE, "DE_isoform_analysis", genotypes, "results",
                      paste0(genotypes, "_downregulated_isoforms.csv")),
  DIU_Full = file.path(BASE, "DIU_IsoformSwitchAnalyzeR", "results", genotypes,
                      paste0(genotypes, "_DIU_isoform_features_full.csv")),
  DIU_Sig = file.path(BASE, "DIU_IsoformSwitchAnalyzeR", "results", genotypes,
                     paste0(genotypes, "_significant_switching_isoforms.csv")),
  stringsAsFactors = FALSE
)
required <- unlist(paths[-1], use.names = FALSE)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing files:\n", paste(missing, collapse = "\n"),
                          call. = FALSE)

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | tolower(x) %in%
      c("", "na", "nan", "n/a", "none", "null", "-")] <- NA_character_
  x
}
ids <- function(x) unique(x[!is.na(x)])
read_result <- function(path, required, iso_col, gene_col) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  absent <- setdiff(required, names(x))
  if (length(absent)) stop(
    "Missing columns in ", path, ": ", paste(absent, collapse = ", "),
    "\nAvailable columns: ", paste(names(x), collapse = ", "), call. = FALSE
  )
  x[[iso_col]] <- clean_id(x[[iso_col]])
  x[[gene_col]] <- clean_id(x[[gene_col]])
  if (anyNA(x[[iso_col]])) stop("Missing isoform ID in ", path, call. = FALSE)
  # Restrict explicitly to the comparison in the supplied analysis scripts.
  if (all(c("condition_1", "condition_2") %in% names(x))) {
    keep <- !is.na(x$condition_1) & !is.na(x$condition_2) &
      x$condition_1 == "Control" & x$condition_2 == "Infested"
    if (nrow(x) > 0 && !any(keep)) stop(
      "No Control -> Infested comparison in ", path, call. = FALSE
    )
    x <- x[keep, , drop = FALSE]
  }
  x <- unique(x)
  if (anyDuplicated(x[[iso_col]])) stop(
    "Multiple rows per isoform in ", path,
    ". Inspect duplicate rows/comparisons before testing.", call. = FALSE
  )
  x
}
as_number <- function(x, label) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  if (any(!is.na(x) & is.na(y))) stop("Non-numeric values in ", label,
                                      call. = FALSE)
  y
}
check_probability <- function(x, label) {
  if (any(!is.na(x) & (!is.finite(x) | x < 0 | x > 1)))
    stop("Invalid probability in ", label, call. = FALSE)
}
check_subset <- function(x, y, label) {
  absent <- setdiff(x, y)
  if (length(absent)) stop(
    label, ": IDs absent from eligible full results (first 10): ",
    paste(head(absent, 10), collapse = ", "), call. = FALSE
  )
}
pct <- function(n, d) if (d > 0) 100 * n / d else NA_real_

results <- membership <- exclusions <- qc <- vector("list", length(genotypes))
for (i in seq_along(genotypes)) {
  g <- genotypes[i]
  message("Processing ", g)
  de <- read_result(paths$DIE_Full[i],
                    c("isoform", "associated_gene", "padj", "log2FoldChange"),
                    "isoform", "associated_gene")
  up <- read_result(paths$DIE_Up[i], c("isoform", "associated_gene"),
                    "isoform", "associated_gene")
  down <- read_result(paths$DIE_Down[i], c("isoform", "associated_gene"),
                      "isoform", "associated_gene")
  du <- read_result(paths$DIU_Full[i],
                    c("isoform_id", "gene_id", "isoform_switch_q_value", "dIF"),
                    "isoform_id", "gene_id")
  ds <- read_result(paths$DIU_Sig[i], c("isoform_id", "gene_id"),
                    "isoform_id", "gene_id")

  # Consistent mappings are essential for a gene-set comparison.
  pair <- unique(rbind(
    data.frame(isoform = de$isoform, gene = de$associated_gene),
    data.frame(isoform = up$isoform, gene = up$associated_gene),
    data.frame(isoform = down$isoform, gene = down$associated_gene),
    data.frame(isoform = du$isoform_id, gene = du$gene_id),
    data.frame(isoform = ds$isoform_id, gene = ds$gene_id)
  ))
  pair <- pair[!is.na(pair$gene), , drop = FALSE]
  if (anyDuplicated(pair$isoform)) stop(
    g, ": conflicting isoform-to-gene mappings across files.", call. = FALSE
  )
  if (length(intersect(up$isoform, down$isoform))) stop(
    g, ": an isoform is in both UP and DOWN lists.", call. = FALSE
  )

  de$padj <- as_number(de$padj, paste(g, "DIE padj"))
  de$log2FoldChange <- as_number(de$log2FoldChange, paste(g, "DIE LFC"))
  du$isoform_switch_q_value <- as_number(
    du$isoform_switch_q_value, paste(g, "DIU q-value"))
  du$dIF <- as_number(du$dIF, paste(g, "DIU dIF"))
  check_probability(de$padj, paste(g, "DIE padj"))
  check_probability(du$isoform_switch_q_value, paste(g, "DIU q-value"))

  de_ok <- is.finite(de$padj) & is.finite(de$log2FoldChange)
  du_ok <- is.finite(du$isoform_switch_q_value) & is.finite(du$dIF)
  de_eligible <- ids(de$associated_gene[de_ok])
  n_iso <- table(du$gene_id)
  du_multi <- names(n_iso[n_iso >= 2])
  du_eligible <- intersect(ids(du$gene_id[du_ok]), du_multi)

  # Check DIE lists match the supplied pipeline thresholds exactly.
  expected_up <- de$isoform[de_ok & de$padj < 0.05 & de$log2FoldChange > 1]
  expected_down <- de$isoform[de_ok & de$padj < 0.05 & de$log2FoldChange < -1]
  if (!setequal(up$isoform, expected_up) ||
      !setequal(down$isoform, expected_down)) stop(
    g, ": significant DIE lists disagree with full results at ",
    "padj < 0.05 and abs(log2FoldChange) > 1. Check file versions.",
    call. = FALSE
  )
  # Preserve extractTopSwitches() calls rather than redefining DIU significance.
  check_subset(ds$isoform_id, du$isoform_id[du_ok],
               paste(g, "DIU significant isoforms"))
  die_sig <- union(ids(up$associated_gene), ids(down$associated_gene))
  diu_sig <- ids(ds$gene_id)
  check_subset(die_sig, de_eligible, paste(g, "DIE significant genes"))
  check_subset(diu_sig, du_eligible, paste(g, "DIU significant genes"))

  universe <- sort(intersect(de_eligible, du_eligible))
  N <- length(universe)
  die <- universe %in% die_sig
  diu <- universe %in% diu_sig
  a <- sum(die & diu)
  b <- sum(die & !diu)
  c <- sum(!die & diu)
  d <- sum(!die & !diu)
  contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                        dimnames = list(DIE = c("Yes", "No"),
                                        DIU = c("Yes", "No")))
  valid <- N > 0 && all(rowSums(contingency) > 0) &&
    all(colSums(contingency) > 0)
  test <- if (valid) fisher.test(contingency, alternative = "two.sided") else NULL
  expected <- if (N > 0) (a + b) * (a + c) / N else NA_real_
  results[[i]] <- data.frame(
    Genotype = g,
    DIE_Eligible_Genes = length(de_eligible),
    DIU_Eligible_Genes = length(du_eligible),
    Common_Background_Genes = N,
    DIE_Significant_Original = length(die_sig),
    DIU_Significant_Original = length(diu_sig),
    Overlap_Original = length(intersect(die_sig, diu_sig)),
    DIE_Significant_Excluded = length(setdiff(die_sig, universe)),
    DIU_Significant_Excluded = length(setdiff(diu_sig, universe)),
    DIE_Significant_Common = a + b,
    DIU_Significant_Common = a + c,
    Both = a, DIE_Only = b, DIU_Only = c, Neither = d,
    Expected_Overlap = expected,
    Observed_Expected_Ratio = if (is.finite(expected) && expected > 0)
      a / expected else NA_real_,
    Percent_DIE_Also_DIU = pct(a, a + b),
    Percent_DIU_Also_DIE = pct(a, a + c),
    Percent_Shared_Of_Union = pct(a, a + b + c),
    Odds_Ratio = if (valid) unname(test$estimate) else NA_real_,
    OR_CI_Lower = if (valid) test$conf.int[1] else NA_real_,
    OR_CI_Upper = if (valid) test$conf.int[2] else NA_real_,
    Fisher_P = if (valid) test$p.value else NA_real_,
    Status = if (valid) "Tested" else "Not testable: empty/constant category"
  )
  membership[[i]] <- data.frame(
    Genotype = rep(g, N), gene_id = universe,
    DIE_Significant = die, DIU_Significant = diu,
    Category = ifelse(die & diu, "Both",
                      ifelse(die, "DIE_only", ifelse(diu, "DIU_only", "Neither")))
  )
  excluded <- sort(setdiff(union(die_sig, diu_sig), universe))
  exclusions[[i]] <- data.frame(
    Genotype = rep(g, length(excluded)), gene_id = excluded,
    DIE_Significant = excluded %in% die_sig,
    DIU_Significant = excluded %in% diu_sig,
    DIE_Eligible = excluded %in% de_eligible,
    DIU_Eligible = excluded %in% du_eligible
  )
  qc[[i]] <- data.frame(
    Genotype = g,
    DIE_Full_Isoforms = nrow(de),
    DIE_Usable_Isoforms = sum(de_ok),
    DIE_Full_Isoforms_Missing_Gene = sum(is.na(de$associated_gene)),
    DIE_Sig_Isoforms_Missing_Gene =
      sum(is.na(up$associated_gene)) + sum(is.na(down$associated_gene)),
    DIU_Full_Isoforms = nrow(du),
    DIU_Usable_Isoforms = sum(du_ok),
    DIU_Full_Isoforms_Missing_Gene = sum(is.na(du$gene_id)),
    DIU_Sig_Isoforms_Missing_Gene = sum(is.na(ds$gene_id))
  )
}

summary <- do.call(rbind, results)
summary$Fisher_FDR_BH <- NA_real_
valid <- which(is.finite(summary$Fisher_P))
summary$Fisher_FDR_BH[valid] <- p.adjust(summary$Fisher_P[valid], method = "BH")
summary$Association <- "Not testable"
for (i in valid) {
  summary$Association[i] <- if (summary$Fisher_FDR_BH[i] >= 0.05) {
    "No significant association"
  } else if (summary$Odds_Ratio[i] > 1) {
    "Enriched overlap"
  } else {
    "Depleted overlap"
  }
}
qc_table <- do.call(rbind, qc)
missing_cols <- grepl("Missing_Gene$", names(qc_table))
if (any(as.matrix(qc_table[, missing_cols, drop = FALSE]) > 0))
  warning("Missing gene IDs were excluded from gene sets; inspect QC.")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write.csv(summary, file.path(OUT, "DIE_DIU_Fisher_summary.csv"), row.names = FALSE)
write.table(summary, file.path(OUT, "DIE_DIU_Fisher_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(do.call(rbind, membership),
          file.path(OUT, "DIE_DIU_common_background_gene_membership.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, exclusions),
          file.path(OUT, "DIE_DIU_significant_genes_excluded_from_background.csv"),
          row.names = FALSE)
write.csv(qc_table, file.path(OUT, "DIE_DIU_Fisher_QC.csv"), row.names = FALSE)
write.csv(paths, file.path(OUT, "DIE_DIU_Fisher_input_paths.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
print(summary[, c("Genotype", "Common_Background_Genes", "Both",
                  "Expected_Overlap", "Odds_Ratio", "Fisher_FDR_BH",
                  "Association")], row.names = FALSE)
cat("\nOutputs saved to:", OUT, "\n")






################################################################################
# ARCHIVAL STATUS: FINAL / EXACT EXTRACT
# Harmonized V. riparia consensus-gene namespace across:
#   Vrip_X_Vrup
#   Vrip_X_Vber
#   Vrip_X_Vace
#
# Annotated VITVri gene IDs are preserved.
# Novel loci are linked by chromosome + strand + >=80% reciprocal overlap.
################################################################################

#!/usr/bin/env Rscript


# ============================================================
# 01. BUILD COMMON V. RIPARIA GENE NAMESPACE
#
# References:
#   Vrip_X_Vrup
#   Vrip_X_Vber
#   Vrip_X_Vace
#
# Purpose:
#
#   Annotated VITVri genes:
#       retain exact reference gene ID
#
#   novelGene_* genes:
#       harmonize across transcript catalogs using
#       chromosome + strand + >=80% reciprocal genomic overlap
#
# IMPORTANT:
#
#   PB transcript IDs are NOT compared across catalogs.
#
#   Novel and annotated genes are NOT collapsed together.
#
#   Coordinate matching is performed at the GENE LOCUS level,
#   after first collapsing all PB transcripts assigned to the
#   same associated_gene.
# ============================================================


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
#
# Restrict to the common Vrip subgenome.
# ============================================================

VRIP_CHR_PATTERN <- "^VITVri588271_v1\\.0\\.hap1\\.chr"


# ============================================================
# 5. EXTRACT TRANSCRIPT COORDINATES FROM GTF
#
# Uses transcript records if available.
#
# If the GTF does not contain explicit transcript features,
# exon coordinates are collapsed instead.
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
# One row per:
#
#   Reference × associated_gene × chromosome × strand
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
# 9. QC:
# ONE ORIGINAL GENE SHOULD NOT FORM MULTIPLE LOCI
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
# 10. ANNOTATED GENES
#
# Preserve exact VITVri ID.
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
    # Candidate overlaps already require chromosome.
    # ignore.strand = FALSE additionally requires strand.
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

            # CRITICAL:
            # only use matches between DIFFERENT transcript
            # catalogs to define a cross-reference consensus.
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
    # 13. CONNECTED COMPONENTS
    #
    # Include every novel locus as a graph vertex so singleton
    # genes receive a component as well.
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
    # 14. DETERMINISTIC COORDINATE-BASED NOVEL IDs
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
# 17. TRANSCRIPT-LEVEL CROSSWALK
#
# Keeps the local PB ID, but adds the common gene ID.
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














