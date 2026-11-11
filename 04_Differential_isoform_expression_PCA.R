#!/usr/bin/env Rscript

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(ggplot2)


# print out all the R packages and version loaded for version control

> sessionInfo()
# R version 4.5.2 (2025-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Rocky Linux 8.10 (Green Obsidian)

# Matrix products: default
# BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.15.so;  LAPACK version 3.9.0

# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

# time zone: America/Los_Angeles
# tzcode source: system (glibc)

# attached base packages:
#  [1] stats4    stats     graphics  utils     datasets  grDevices methods  
# [8] base     

# other attached packages:
#  [1] KEGGREST_1.50.0             enrichplot_1.30.3          
# [3] GO.db_3.22.0                AnnotationDbi_1.72.0       
# [5] clusterProfiler_4.18.2      lubridate_1.9.4            
# [7] forcats_1.0.1               stringr_1.6.0              
# [9] dplyr_1.2.1                 purrr_1.2.2                
# [11] readr_2.1.6                 tidyr_1.3.2                
# [13] tibble_3.3.0                ggplot2_4.0.2              
# [15] tidyverse_2.0.0             DESeq2_1.50.2              
# [17] SummarizedExperiment_1.40.0 Biobase_2.70.0             
# [19] MatrixGenerics_1.22.0       matrixStats_1.5.0          
# [21] GenomicRanges_1.62.0        Seqinfo_1.0.0              
# [23] IRanges_2.44.0              S4Vectors_0.48.0           
 #[25] BiocGenerics_0.56.0         generics_0.1.4             
 #[27] setwidth_1.1.0              colorout_1.3-3             

# loaded via a namespace (and not attached):
#  [1] DBI_1.2.3               gson_0.1.0              rlang_1.2.0            
# [4] magrittr_2.0.4          DOSE_4.4.0              compiler_4.5.2         
# [7] RSQLite_2.4.4           systemfonts_1.3.1       png_0.1-8              
# [10] vctrs_0.7.3             reshape2_1.4.5          pkgconfig_2.0.3        
# [13] crayon_1.5.3            fastmap_1.2.0           XVector_0.50.0         
# [16] tzdb_0.5.0              bit_4.6.0               cachem_1.1.0           
# [19] aplot_0.2.9             jsonlite_2.0.0          blob_1.2.4             
# [22] DelayedArray_0.36.0     BiocParallel_1.44.0     parallel_4.5.2         
# [25] R6_2.6.1                stringi_1.8.7           RColorBrewer_1.1-3     
# [28] GOSemSim_2.36.0         Rcpp_1.1.1-1            ggtangle_0.0.8         
# [31] R.utils_2.13.0          Matrix_1.7-4            splines_4.5.2          
# [34] igraph_2.2.1            timechange_0.3.0        tidyselect_1.2.1       
# [37] qvalue_2.42.0           dichromat_2.0-0.1       abind_1.4-8            
# [40] codetools_0.2-20        lattice_0.22-7          plyr_1.8.9             
# [43] treeio_1.34.0           withr_3.0.2             S7_0.2.1               
# [46] gridGraphics_0.5-1      Biostrings_2.78.0       pillar_1.11.1          
# [49] ggtree_4.0.1            ggfun_0.2.0             hms_1.1.4              
# [52] scales_1.4.0            tidytree_0.4.6          glue_1.8.0             
# [55] gdtools_0.4.4           lazyeval_0.2.2          tools_4.5.2            
# [61] ggiraph_0.9.2           fs_1.6.6                fastmatch_1.1-6        
 #[64] cowplot_1.2.0           grid_4.5.2              ape_5.8-1              
 #[67] nlme_3.1-168            patchwork_1.3.2         cli_3.6.5              
# [70] rappdirs_0.3.3          fontBitstreamVera_0.1.1 S4Arrays_1.10.0        
 #[73] gtable_0.3.6            R.methodsS3_1.8.2       yulab.utils_0.2.1      
# [76] fontquiver_0.2.1        digest_0.6.39           SparseArray_1.10.2     
# [79] ggrepel_0.9.8           ggplotify_0.1.3         htmlwidgets_1.6.4      
# [82] farver_2.1.2            htmltools_0.5.8.1       memoise_2.0.1          
# [85] R.oo_1.27.1             lifecycle_1.0.5         httr_1.4.7             
 #[88] fontLiberation_0.1.0    bit64_4.6.0-1 







# ============================================================
# 1. Paths
# ============================================================

BASE_DIR <- "/rhome/wcoll010/shared/pacbio/IsoSeq_analyses_2/DE_isoform_analysis"

SUMMARY_FILE <- file.path(
  BASE_DIR,
  "DE_summary_10_genotypes.tsv"
)


# ============================================================
# 2. Genotypes and panel order
# ============================================================

GENOTYPES <- c(
  "101-14", "99R", "1103", "SO4", "ripa",
  "cab33", "1616", "3309", "sw", "5BB"
)

GENOTYPE_LABELS <- c(
  "101-14" = "101-14", "99R" = "99R", "1103" = "1103",
  "SO4" = "SO4", "ripa" = "V. riparia", "cab33" = "Cab33",
  "1616" = "1616", "3309" = "3309", "sw" = "SW", "5BB" = "5BB"
)
PCA_DIR <- file.path(BASE_DIR, "PCA_all_10_genotypes")
dir.create(PCA_DIR, recursive = TRUE, showWarnings = FALSE)

# For PCA plot, use 500 most variable isoforms in each genotype, after expression filtering and VST.
PCA_NTOP <- 500L

# Explicit input selection, this allows the genotypes with duplicated libraries (sequenced twice due to low output) to be merged together
COUNT_FILES <- setNames(
  file.path(BASE_DIR, GENOTYPES,
            paste0(GENOTYPES, "_sqanti_filtered_counts_DESeq2.tsv")),
  GENOTYPES
)
for (g in c("99R", "cab33", "1616")) {
  COUNT_FILES[g] <- file.path(
    BASE_DIR, g, paste0(g, "_sqanti_filtered_counts_DESeq2_merged.tsv")
  )
}
CAB33_UNMERGED <- file.path(
  BASE_DIR, "cab33", "cab33_sqanti_filtered_counts_DESeq2.tsv"
)
required_inputs <- c(
  unname(COUNT_FILES), CAB33_UNMERGED,
  file.path(BASE_DIR, GENOTYPES,
            paste0(GENOTYPES, "_sqanti_retained_annotation.tsv"))
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing input files:\n", paste(missing_inputs, collapse = "\n"))
}



# ============================================================
# 4. Function to process one genotype
# ============================================================

process_genotype <- function(genotype) {

  cat("\n")
  cat("============================================================\n")
  cat("GENOTYPE:", genotype, "\n")
  cat("============================================================\n")


  # ==========================================================
  # Paths
  # ==========================================================

  genotype_dir <- file.path(
    BASE_DIR,
    genotype
  )

  counts_file <- COUNT_FILES[[genotype]]
  cat("\nUsing count file:\n", counts_file, "\n")

  class_file <- file.path(
    genotype_dir,
    paste0(
      genotype,
      "_sqanti_retained_annotation.tsv"
    )
  )

  output_dir <- file.path(
    genotype_dir,
    "results"
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  # ==========================================================
  # Check inputs
  # ==========================================================

  if (!file.exists(counts_file)) {
    stop(
      genotype,
      ": count file not found:\n",
      counts_file
    )
  }

  if (!file.exists(class_file)) {
    stop(
      genotype,
      ": annotation file not found:\n",
      class_file
    )
  }


  # ==========================================================
  # Load counts
  # ==========================================================

  counts_data <- read.table(
    counts_file,
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
  )

  cat(
    "\nInitial count matrix:",
    nrow(counts_data),
    "isoforms x",
    ncol(counts_data),
    "sequencing libraries\n"
  )

  cat("\nInitial samples:\n")
  print(colnames(counts_data))


  # ==========================================================
  # Load SQANTI annotation
  # ==========================================================

  class_data <- read.table(
    class_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE
  ) %>%
    select(
      isoform,
      associated_gene,
      structural_category
    )


  # ==========================================================
  # Counts are already at biological-sample level; do not merge again.
  # Confirm the six names observed in the inspected input headers.
  # Column order is immaterial: treatment is assigned from each name.
  
  if (genotype %in% c("99R", "cab33", "1103", "1616")) {
    expected_samples <- paste0(
      genotype, "_", rep(c("con", "inf"), each = 3),
      "_R", rep(1:3, 2)
    )
    if (!setequal(colnames(counts_data), expected_samples)) {
      stop(genotype, ": unexpected biological-sample names in selected count file.")
    }
  }

  # ==========================================================
  # Final sample checks
  # ==========================================================

  if (anyDuplicated(colnames(counts_data))) {

    duplicated_names <- unique(
      colnames(counts_data)[
        duplicated(colnames(counts_data))
      ]
    )

    stop(
      genotype,
      ": duplicate sample names after preprocessing: ",
      paste(duplicated_names, collapse = ", ")
    )
  }


  cat("\nFinal biological samples:\n")
  print(colnames(counts_data))

  cat(
    "\nFinal matrix:",
    nrow(counts_data),
    "isoforms x",
    ncol(counts_data),
    "biological samples\n"
  )


  # ==========================================================
  # Metadata
  # ==========================================================

  sample_names <- colnames(counts_data)

  meta_data <- data.frame(
    row.names = sample_names,

    condition = case_when(
      grepl("_con_", sample_names) ~ "Control",
      grepl("_inf_", sample_names) ~ "Infested",
      TRUE ~ NA_character_
    )
  )


  if (any(is.na(meta_data$condition))) {

    bad_samples <- rownames(meta_data)[
      is.na(meta_data$condition)
    ]

    stop(
      genotype,
      ": unable to assign treatment for: ",
      paste(bad_samples, collapse = ", ")
    )
  }


  meta_data$condition <- factor(
    meta_data$condition,
    levels = c(
      "Control",
      "Infested"
    )
  )


  cat("\nMetadata:\n")
  print(meta_data)

  cat("\nBiological replicates per condition:\n")
  print(table(meta_data$condition))


  # ==========================================================
  # Check biological replication
  # ==========================================================

  n_control <- sum(
    meta_data$condition == "Control"
  )

  n_infested <- sum(
    meta_data$condition == "Infested"
  )


  if (n_control < 2 || n_infested < 2) {
    stop(
      genotype,
      ": fewer than 2 biological replicates in a treatment."
    )
  }


  # ==========================================================
  # DESeq2 object
  # ==========================================================

  dds <- DESeqDataSetFromMatrix(
    countData = counts_data,
    colData = meta_data,
    design = ~ condition
  )


  dds <- estimateSizeFactors(dds)


  cat("\nSize factors:\n")
  print(sizeFactors(dds))


  # ==========================================================
  # Expression filtering
  #
  # Retain isoforms with FPM > 0.1 in at least two
  # biological replicates of EITHER treatment.
  # ==========================================================

  fpm_mat <- fpm(dds)


  control_samples <- rownames(meta_data)[
    meta_data$condition == "Control"
  ]

  infested_samples <- rownames(meta_data)[
    meta_data$condition == "Infested"
  ]


  control_keep <- rowSums(
    fpm_mat[
      ,
      control_samples,
      drop = FALSE
    ] > 0.1
  ) >= 2


  infested_keep <- rowSums(
    fpm_mat[
      ,
      infested_samples,
      drop = FALSE
    ] > 0.1
  ) >= 2


  keep <- control_keep | infested_keep


  n_before_filter <- nrow(dds)
  n_after_filter <- sum(keep)


  cat(
    "\nIsoforms before expression filtering:",
    n_before_filter,
    "\n"
  )

  cat(
    "Isoforms retained:",
    n_after_filter,
    "\n"
  )

  cat(
    "Percent retained:",
    round(
      100 * n_after_filter / n_before_filter,
      2
    ),
    "%\n"
  )


  dds <- dds[
    keep,
  ]


  # ==========================================================
  # Differential isoform expression
  # ==========================================================

  cat("\nRunning DESeq2...\n")

  dds <- DESeq(dds)


  res <- results(
    dds,
    contrast = c(
      "condition",
      "Infested",
      "Control"
    ),
    alpha = 0.05
  )


  cat("\nDESeq2 summary:\n")
  print(summary(res))


  # ==========================================================
  # Full DE result table
  # ==========================================================

  res_df <- as.data.frame(res) %>%
    rownames_to_column(
      var = "isoform"
    ) %>%
    left_join(
      class_data,
      by = "isoform"
    )


  all_results_file <- file.path(
    output_dir,
    paste0(
      genotype,
      "_all_isoform_DE_results.csv"
    )
  )


  write.csv(
    res_df,
    all_results_file,
    row.names = FALSE
  )


  # ==========================================================
  # Significant UP isoforms
  #
  # Infested > Control
  # ==========================================================

  up_isoforms <- res_df %>%
    filter(
      !is.na(padj),
      padj < 0.05,
      log2FoldChange > 1
    ) %>%
    arrange(
      desc(log2FoldChange)
    )


  n_up <- nrow(up_isoforms)


  cat(
    "\nSignificantly UP isoforms:",
    n_up,
    "\n"
  )


  write.csv(
    up_isoforms,
    file.path(
      output_dir,
      paste0(
        genotype,
        "_upregulated_isoforms.csv"
      )
    ),
    row.names = FALSE
  )


  # ==========================================================
  # Significant DOWN isoforms
  #
  # Infested < Control
  # ==========================================================

  down_isoforms <- res_df %>%
    filter(
      !is.na(padj),
      padj < 0.05,
      log2FoldChange < -1
    ) %>%
    arrange(
      log2FoldChange
    )


  n_down <- nrow(down_isoforms)


  cat(
    "Significantly DOWN isoforms:",
    n_down,
    "\n"
  )


  write.csv(
    down_isoforms,
    file.path(
      output_dir,
      paste0(
        genotype,
        "_downregulated_isoforms.csv"
      )
    ),
    row.names = FALSE
  )


  # ==========================================================
  # PCA transformation
  # ==========================================================

  cat("\nCalculating VST for PCA...\n")

  vsd <- vst(
    dds,
    blind = FALSE
  )


  pcaData <- plotPCA(
    vsd,
    intgroup = "condition",
    ntop = PCA_NTOP,
    returnData = TRUE
  )


  percentVar <- round(
    100 * attr(
      pcaData,
      "percentVar"
    ),
    1
  )


  # ==========================================================
  # PCA 1: WITH sample labels
  # ==========================================================

  pca_labeled <- ggplot(
    pcaData,
    aes(
      x = PC1,
      y = PC2,
      color = condition,
      label = name
    )
  ) +
    geom_point(
      size = 4
    ) +
    geom_text_repel(
      size = 3,
      max.overlaps = Inf
    ) +
    xlab(
      paste0(
        "PC1: ",
        percentVar[1],
        "% variance"
      )
    ) +
    ylab(
      paste0(
        "PC2: ",
        percentVar[2],
        "% variance"
      )
    ) +
    ggtitle(
      paste0(
        genotype,
        ": Infested vs Control"
      )
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5
      )
    )


  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        genotype,
        "_PCA_with_sample_labels.pdf"
      )
    ),
    plot = pca_labeled,
    width = 7,
    height = 6
  )


  # ==========================================================
  # PCA 2: WITHOUT sample labels
  # ==========================================================

  pca_unlabeled <- ggplot(
    pcaData,
    aes(
      x = PC1,
      y = PC2,
      color = condition
    )
  ) +
    geom_point(
      size = 4
    ) +
    xlab(
      paste0(
        "PC1: ",
        percentVar[1],
        "% variance"
      )
    ) +
    ylab(
      paste0(
        "PC2: ",
        percentVar[2],
        "% variance"
      )
    ) +
    ggtitle(
      paste0(
        "PCA: ",
        genotype,
        " Infested vs Control"
      )
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5
      )
    )


  # Save directly.
  # Do NOT print() on the cluster because there is no X11
  # graphics display.

  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        genotype,
        "_PCA_no_sample_labels.pdf"
      )
    ),
    plot = pca_unlabeled,
    width = 7,
    height = 6
  )


  # ==========================================================
  # Volcano data
  # ==========================================================

  volcano_data <- res_df %>%
    mutate(

      diffexpressed = case_when(

        !is.na(padj) &
          padj < 0.05 &
          log2FoldChange > 1 ~ "UP",

        !is.na(padj) &
          padj < 0.05 &
          log2FoldChange < -1 ~ "DOWN",

        TRUE ~ "NO"
      ),

      plot_padj = ifelse(
        is.na(padj),
        NA_real_,
        pmax(
          padj,
          1e-300
        )
      )
    )


  # ==========================================================
  # Volcano plot
  # ==========================================================

  volcano_plot <- ggplot(
    volcano_data,
    aes(
      x = log2FoldChange,
      y = -log10(plot_padj),
      color = diffexpressed
    )
  ) +
    geom_point(
      alpha = 0.6,
      size = 1.5
    ) +
    geom_vline(
      xintercept = c(
        -1,
        1
      ),
      color = "black",
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      color = "black",
      linetype = "dashed"
    ) +
    scale_color_manual(
      values = c(
        "DOWN" = "blue",
        "NO" = "grey",
        "UP" = "red"
      )
    ) +
    labs(
      title = paste0(
        genotype,
        ": Infested vs Control"
      ),
      x = "Log2 fold change (Infested / Control)",
      y = "-Log10 adjusted P-value",
      color = NULL
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5
      )
    )


  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        genotype,
        "_volcano.pdf"
      )
    ),
    plot = volcano_plot,
    width = 7,
    height = 6
  )


  # ==========================================================
  # Return summary statistics for this genotype
  # ==========================================================

  result_summary <- data.frame(

    Genotype = genotype,

    Control_samples = n_control,

    Infested_samples = n_infested,

    SQANTI_retained_isoforms =
      n_before_filter,

    Expression_filtered_isoforms =
      n_after_filter,

    Percent_retained =
      round(
        100 * n_after_filter /
          n_before_filter,
        2
      ),

    Significant_UP = n_up,

    Significant_DOWN = n_down,

    Significant_total =
      n_up + n_down,

    PC1_variance =
      percentVar[1],

    PC2_variance =
      percentVar[2]
  )


  cat("\nCompleted:", genotype, "\n")

  # Keep scores and variance metadata together when leaving this function.
  pcaData$genotype <- genotype
  pcaData$PC1_variance <- percentVar[1]
  pcaData$PC2_variance <- percentVar[2]
  pcaData$isoforms_used <- min(PCA_NTOP, nrow(vsd))
  write.csv(
    pcaData,
    file.path(output_dir, paste0(genotype, "_PCA_scores.csv")),
    row.names = FALSE
  )
  return(list(summary = result_summary, pca = pcaData))
}


# ============================================================
# 5. Run all ten genotypes
# ============================================================

all_summaries <- list()
all_pca <- list()


for (genotype in GENOTYPES) {

  result <- process_genotype(
    genotype
  )

  all_summaries[[
    genotype
  ]] <- result$summary
  all_pca[[genotype]] <- result$pca
}


# ============================================================
# 6. Combined summary table
# ============================================================

summary_df <- bind_rows(
  all_summaries
)


write.table(
  summary_df,
  SUMMARY_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("ALL 10 GENOTYPES COMPLETED\n")
cat("============================================================\n")
cat("\n")

print(summary_df)

cat("\nSummary written to:\n")
cat(SUMMARY_FILE, "\n")


# ============================================================
# 7. One faceted PCA figure: two rows x five columns
# ============================================================
# These are ten independent PCAs. PC axes/loadings differ among
# genotypes, so cross-panel positions are not directly comparable.
# facet_wrap has shared axis titles; genotype-specific variance
# percentages are therefore shown in each facet heading.

pca_combined <- bind_rows(all_pca)
pca_combined$genotype <- factor(pca_combined$genotype, levels = GENOTYPES)
pca_combined$condition <- factor(
  pca_combined$condition, levels = c("Control", "Infested")
)

facet_labels <- setNames(
  vapply(GENOTYPES, function(g) {
    x <- all_pca[[g]]
    sprintf("%s\nPC1: %.1f%% | PC2: %.1f%%",
            GENOTYPE_LABELS[[g]], x$PC1_variance[1], x$PC2_variance[1])
  }, character(1)),
  GENOTYPES
)

pca_faceted <- ggplot(
  pca_combined,
  aes(x = PC1, y = PC2, color = condition, shape = condition)
) +
  geom_point(size = 3.2, alpha = 0.95) +
  facet_wrap(
    ~ genotype,
    nrow = 2, ncol = 5,
    scales = "free",
    labeller = labeller(genotype = as_labeller(facet_labels))
  ) +
  scale_color_manual(
    name = "Treatment",
    values = c("Control" = "#F8766D", "Infested" = "#00BFC4"),
    breaks = c("Control", "Infested"), drop = FALSE
  ) +
  scale_shape_manual(
    name = "Treatment",
    values = c("Control" = 16, "Infested" = 17),
    breaks = c("Control", "Infested"), drop = FALSE
  ) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(x = "PC1", y = "PC2") +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "grey65"),
    strip.text = element_text(size = 11, face = "bold", lineheight = 1.1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.spacing = grid::unit(1, "lines"),
    axis.text = element_text(size = 9, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# PDF is vector output; the PNG is useful for quick viewing or PowerPoint.
ggsave(
  file.path(PCA_DIR, "PCA_all_10_genotypes_faceted.pdf"),
  plot = pca_faceted, width = 16, height = 7.5, units = "in"
)
ggsave(
  file.path(PCA_DIR, "PCA_all_10_genotypes_faceted.png"),
  plot = pca_faceted, width = 16, height = 7.5, units = "in",
  dpi = 300, bg = "white"
)

# A second version identifies replicates with short R1/R2/R3 labels.
pca_combined$replicate_label <- sub("^.*_(R[0-9]+)$", "\\1", pca_combined$name)
pca_faceted_labeled <- pca_faceted +
  geom_text_repel(
    data = pca_combined,
    aes(label = replicate_label),
    size = 3, max.overlaps = Inf, seed = 123,
    show.legend = FALSE
  )
ggsave(
  file.path(PCA_DIR, "PCA_all_10_genotypes_faceted_with_labels.pdf"),
  plot = pca_faceted_labeled, width = 16, height = 7.5, units = "in"
)

# Scores and variance percentages allow later restyling without rerunning DE.
write.csv(pca_combined, file.path(PCA_DIR, "PCA_all_10_genotypes_scores.csv"),
          row.names = FALSE)
write.csv(
  summary_df[, c("Genotype", "PC1_variance", "PC2_variance")],
  file.path(PCA_DIR, "PCA_all_10_genotypes_variance.csv"), row.names = FALSE
)
saveRDS(pca_faceted, file.path(PCA_DIR, "PCA_all_10_genotypes_faceted.rds"))
cat("\nFaceted PCA figures and scores written to:\n", PCA_DIR, "\n")







