### Perform functional enrichment analyses 

### what KEGG and GO terms are overrepresented among genes encoding differentially expressed isoforms across geneotypes 
### gene-sets (DIE (differential isoform expression) genes or DIU (differential isoform usage) genes were tested relative to the genes in the SQANTI3-filtered reference transcriptomes)



#!/usr/bin/env Rscript

library(tidyverse)
library(clusterProfiler)
library(GO.db)
library(enrichplot)
library(ggplot2)
library(KEGGREST)


# ============================================================
# 1. Global setup
# ============================================================

base_dir <- "/rhome/wcoll010/shared/pacbio/IsoSeq_analyses_2"

de_dir <- file.path(
  base_dir,
  "DE_isoform_analysis"
)

sqanti_dir <- file.path(
  base_dir,
  "sqanti_results"
)

egg_dir <- file.path(
  base_dir,
  "eggnog_results"
)

# NEW output directory -- does not overwrite other analyses
out_root <- file.path(
  base_dir,
  "GO_KEGG_enrichment_SQANTI_background"
)

dir.create(
  out_root,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Genotype -> current reference catalog mapping
# ============================================================

genotype_reference <- c(

  "101-14" = "Vrip_X_Vrup",
  "3309"   = "Vrip_X_Vrup",
  "sw"     = "Vrip_X_Vrup",

  "1103"   = "Vber_X_Vrup",
  "99R"    = "Vber_X_Vrup",

  "5BB"    = "Vrip_X_Vber",
  "SO4"    = "Vrip_X_Vber",

  "1616"   = "Vrip_X_Vace",

  "cab33"  = "CabSauv08",

  "ripa"   = "Vripa"
)

genotypes <- names(genotype_reference)


# ============================================================
# 3. GO master map
# ============================================================

cat("Building GO term map...\n")

go_master_map <- AnnotationDbi::select(
  GO.db,
  keys = keys(GO.db),
  columns = c(
    "GOID",
    "ONTOLOGY",
    "TERM"
  ),
  keytype = "GOID"
) %>%
  dplyr::rename(
    ID = GOID,
    Description = TERM
  ) %>%
  distinct()


# ============================================================
# 4. KEGG pathway descriptions
# ============================================================

cat("Retrieving KEGG pathway descriptions...\n")

kegg_pathway_list <- KEGGREST::keggList(
  "pathway",
  "ko"
)

KEGG_master_names <- tibble(

  ID = names(
    kegg_pathway_list
  ),

  Description = unname(
    kegg_pathway_list
  )

) %>%
  mutate(

    # path:ko00010 -> ko00010
    ID = sub(
      "^path:",
      "",
      ID
    ),

    Description = sub(
      " - Reference pathway$",
      "",
      Description
    )
  )

cat(
  "Retrieved",
  nrow(KEGG_master_names),
  "KEGG pathway names.\n"
)


# ============================================================
# 5. Robust eggNOG reader
# ============================================================

read_eggnog <- function(file) {

  lines <- readLines(
    file,
    warn = FALSE
  )

  header_line <- grep(
    "^#query\t",
    lines
  )

  if (length(header_line) != 1) {
    stop(
      "Could not uniquely identify eggNOG #query header in:\n",
      file
    )
  }

  eggnog <- read.delim(
    file,
    skip = header_line - 1,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # First column is PB query ID
  colnames(eggnog)[1] <- "isoform"

  eggnog %>%
    distinct(
      isoform,
      .keep_all = TRUE
    )
}


# ============================================================
# 6. Main genotype loop
# ============================================================

for (geno in genotypes) {

  reference <- genotype_reference[[geno]]

  cat("\n")
  cat("============================================================\n")
  cat("GENOTYPE:", geno, "\n")
  cat("REFERENCE CATALOG:", reference, "\n")
  cat("BACKGROUND: FULL SQANTI-RETAINED CATALOG\n")
  cat("============================================================\n")


  input_dir <- file.path(
    de_dir,
    geno,
    "results"
  )

  geno_out <- file.path(
    out_root,
    geno
  )

  dir.create(
    geno_out,
    recursive = TRUE,
    showWarnings = FALSE
  )


  # ==========================================================
  # 7. Current SQANTI classification
  #
  # NEVER uses old_* folders.
  # ==========================================================

  sqanti_file <- file.path(
    sqanti_dir,
    paste0(
      reference,
      "_sqanti"
    ),
    "filter_results",
    paste0(
      reference,
      "_filtered_RulesFilter_result_classification.txt"
    )
  )


  # ==========================================================
  # 8. Current reference-level eggNOG annotation
  # ==========================================================

  egg_file <- file.path(
    egg_dir,
    paste0(
      reference,
      "_eggnog.emapper.annotations"
    )
  )


  if (!file.exists(sqanti_file)) {

    cat(
      "[!] Missing SQANTI file:",
      sqanti_file,
      "\n"
    )

    next
  }


  if (!file.exists(egg_file)) {

    cat(
      "[!] Missing eggNOG file:",
      egg_file,
      "\n"
    )

    next
  }


  # ==========================================================
  # 9. Build full SQANTI-retained transcript -> gene bridge
  # ==========================================================

  sqanti <- read.delim(
    sqanti_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE
  )


  bridge <- sqanti %>%

    filter(
      filter_result == "Isoform"
    ) %>%

    dplyr::select(
      isoform,
      associated_gene
    ) %>%

    filter(
      !is.na(associated_gene),
      associated_gene != ""
    ) %>%

    distinct()


  cat(
    "SQANTI-retained isoforms with associated genes:",
    nrow(bridge),
    "\n"
  )

  cat(
    "Unique SQANTI-retained genes:",
    n_distinct(bridge$associated_gene),
    "\n"
  )


  # ==========================================================
  # 10. Add eggNOG annotations to full SQANTI background
  # ==========================================================

  full_eggnog <- read_eggnog(
    egg_file
  )


  gene_annot_map <- bridge %>%

    inner_join(
      full_eggnog,
      by = "isoform"
    )


  cat(
    "SQANTI-retained isoforms with eggNOG records:",
    nrow(gene_annot_map),
    "\n"
  )


  # ==========================================================
  # 11. Full SQANTI gene-level background
  #
  # This is intentionally NOT restricted to the isoforms that
  # entered DESeq2.
  # ==========================================================

  sqanti_background_genes <- bridge %>%

    pull(
      associated_gene
    ) %>%

    unique()


  # ==========================================================
  # 12. Gene-level GO mappings
  #
  # Multiple isoforms belonging to the same gene and carrying
  # the same GO term count only once.
  # ==========================================================

  GO_links <- gene_annot_map %>%

    dplyr::select(
      associated_gene,
      GOs
    ) %>%

    filter(
      !is.na(associated_gene),
      associated_gene != "",
      !is.na(GOs),
      GOs != "",
      GOs != "-"
    ) %>%

    separate_rows(
      GOs,
      sep = ","
    ) %>%

    mutate(
      GOs = trimws(GOs)
    ) %>%

    filter(
      GOs != ""
    ) %>%

    transmute(
      ID = GOs,
      associated_gene
    ) %>%

    distinct(
      ID,
      associated_gene
    ) %>%

    left_join(
      go_master_map,
      by = "ID"
    )


  cat(
    "Gene-GO links in SQANTI background:",
    nrow(GO_links),
    "\n"
  )


  # ==========================================================
  # 13. Gene-level KEGG pathway mappings
  #
  # Keep koXXXXX versions rather than duplicate mapXXXXX forms.
  # ==========================================================

  KEGG_links <- gene_annot_map %>%

    dplyr::select(
      associated_gene,
      KEGG_Pathway
    ) %>%

    filter(
      !is.na(associated_gene),
      associated_gene != "",
      !is.na(KEGG_Pathway),
      KEGG_Pathway != "",
      KEGG_Pathway != "-"
    ) %>%

    separate_rows(
      KEGG_Pathway,
      sep = ","
    ) %>%

    mutate(
      KEGG_Pathway = trimws(
        KEGG_Pathway
      )
    ) %>%

    filter(
      grepl(
        "^ko[0-9]+$",
        KEGG_Pathway
      )
    ) %>%

    transmute(
      ID = KEGG_Pathway,
      associated_gene
    ) %>%

    distinct(
      ID,
      associated_gene
    )


  KEGG_names <- KEGG_links %>%

    distinct(
      ID
    ) %>%

    left_join(
      KEGG_master_names,
      by = "ID"
    ) %>%

    mutate(
      Description = ifelse(
        is.na(Description),
        ID,
        Description
      )
    )


  cat(
    "Gene-KEGG links in SQANTI background:",
    nrow(KEGG_links),
    "\n"
  )


  # ==========================================================
  # 14. UP / DOWN loop
  # ==========================================================

  for (direction in c(
    "upregulated",
    "downregulated"
  )) {

    cat("\n")
    cat("------------------------------------------------------------\n")
    cat("DIRECTION:", direction, "\n")
    cat("------------------------------------------------------------\n")


    de_file <- file.path(
      input_dir,
      paste0(
        geno,
        "_",
        direction,
        "_isoforms_eggnog.csv"
      )
    )


    if (!file.exists(de_file)) {

      cat(
        "[!] Missing DE file:",
        de_file,
        "\n"
      )

      next
    }


    de_data <- read.csv(
      de_file,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )


    # ========================================================
    # 15. Collapse DE isoforms to unique genes
    # ========================================================

    target_genes <- de_data %>%

      filter(
        !is.na(associated_gene),
        associated_gene != ""
      ) %>%

      pull(
        associated_gene
      ) %>%

      unique()


    # Target must exist in this SQANTI reference catalog
    target_genes <- intersect(
      target_genes,
      sqanti_background_genes
    )


    cat(
      "DE isoforms:",
      nrow(de_data),
      "\n"
    )

    cat(
      "Unique DE genes:",
      length(target_genes),
      "\n"
    )


    if (length(target_genes) < 5) {

      cat(
        "[!] Too few target genes.\n"
      )

      next
    }


    # ========================================================
    # ========================================================
    # 16. KEGG PATHWAY ENRICHMENT
    # ========================================================
    # ========================================================

    cat("\n[+] KEGG enrichment using SQANTI background...\n")


    kegg_universe <- unique(
      KEGG_links$associated_gene
    )


    kegg_targets <- intersect(
      target_genes,
      kegg_universe
    )


    cat(
      "KEGG-annotated SQANTI background genes:",
      length(kegg_universe),
      "\n"
    )

    cat(
      "KEGG-annotated DE genes:",
      length(kegg_targets),
      "\n"
    )


    if (length(kegg_targets) >= 5) {


      k_res <- enricher(

        gene = kegg_targets,

        universe = kegg_universe,

        TERM2GENE = KEGG_links %>%
          dplyr::select(
            ID,
            associated_gene
          ),

        TERM2NAME = KEGG_names,

        pvalueCutoff = 0.05,

        pAdjustMethod = "BH",

        minGSSize = 5
      )


      if (
        !is.null(k_res) &&
        nrow(as.data.frame(k_res)) > 0
      ) {


        k_res@result <- k_res@result %>%

          filter(
            Count >= 5
          )


        if (nrow(k_res@result) > 0) {


          kegg_base <- file.path(
            geno_out,
            paste0(
              geno,
              "_",
              direction,
              "_KEGG_SQANTI_background"
            )
          )


          # Full enrichment result
          write.csv(
            as.data.frame(k_res),
            paste0(
              kegg_base,
              "_results.csv"
            ),
            row.names = FALSE
          )


          # Contributing genes
          kegg_gene_list <- k_res@result %>%

            dplyr::select(
              ID,
              Description,
              GeneRatio,
              BgRatio,
              pvalue,
              p.adjust,
              Count,
              geneID
            ) %>%

            separate_rows(
              geneID,
              sep = "/"
            ) %>%

            dplyr::rename(
              associated_gene = geneID
            )


          write.csv(
            kegg_gene_list,
            paste0(
              kegg_base,
              "_gene_list.csv"
            ),
            row.names = FALSE
          )


          # Functional pathway names on y-axis
          pdf(
            paste0(
              kegg_base,
              "_dotplot.pdf"
            ),
            width = 11,
            height = ifelse(
              nrow(k_res@result) < 10,
              7,
              9
            )
          )


          print(

            dotplot(
              k_res,
              showCategory = 20,
              label_format = 55
            ) +

              ggtitle(
                paste(
                  "KEGG:",
                  geno,
                  direction,
                  "(SQANTI background)"
                )
              ) +

              theme_bw()
          )


          dev.off()


          cat(
            "[+] Saved KEGG SQANTI-background results.\n"
          )
        }

      } else {

        cat(
          "[-] No significant KEGG pathways.\n"
        )
      }

    } else {

      cat(
        "[-] Too few KEGG-annotated target genes.\n"
      )
    }


    # ========================================================
    # ========================================================
    # 17. GO ENRICHMENT
    #
    # NO redundancy reduction.
    # All significant related GO terms are retained.
    # ========================================================
    # ========================================================

    for (ont in c(
      "BP",
      "MF",
      "CC"
    )) {


      cat(
        "\n[+] GO",
        ont,
        "using SQANTI background...\n"
      )


      ont_links <- GO_links %>%

        filter(
          ONTOLOGY == ont
        ) %>%

        dplyr::select(
          ID,
          associated_gene
        ) %>%

        distinct()


      ont_names <- go_master_map %>%

        filter(
          ONTOLOGY == ont
        ) %>%

        dplyr::select(
          ID,
          Description
        ) %>%

        distinct()


      if (nrow(ont_links) == 0) {
        next
      }


      # ======================================================
      # SQANTI catalog background for this ontology
      # ======================================================

      go_universe <- unique(
        ont_links$associated_gene
      )


      go_targets <- intersect(
        target_genes,
        go_universe
      )


      cat(
        "GO",
        ont,
        "SQANTI background genes:",
        length(go_universe),
        "\n"
      )

      cat(
        "GO",
        ont,
        "DE genes:",
        length(go_targets),
        "\n"
      )


      if (length(go_targets) < 5) {
        next
      }


      # ======================================================
      # Run enrichment
      # ======================================================

      go_res <- enricher(

        gene = go_targets,

        universe = go_universe,

        TERM2GENE = ont_links,

        TERM2NAME = ont_names,

        pvalueCutoff = 0.05,

        pAdjustMethod = "BH",

        minGSSize = 5
      )


      if (
        is.null(go_res) ||
        nrow(as.data.frame(go_res)) == 0
      ) {

        cat(
          "[-] No significant GO",
          ont,
          "terms.\n"
        )

        next
      }


      go_res@result <- go_res@result %>%

        filter(
          Count >= 5
        )


      if (nrow(go_res@result) == 0) {
        next
      }


      go_base <- file.path(
        geno_out,
        paste0(
          geno,
          "_",
          direction,
          "_GO_",
          ont,
          "_SQANTI_background"
        )
      )


      # ======================================================
      # Full GO result -- NO redundancy reduction
      # ======================================================

      write.csv(
        as.data.frame(go_res),
        paste0(
          go_base,
          "_results.csv"
        ),
        row.names = FALSE
      )


      # ======================================================
      # Contributing genes for every term
      # ======================================================

      go_gene_list <- go_res@result %>%

        dplyr::select(
          ID,
          Description,
          GeneRatio,
          BgRatio,
          pvalue,
          p.adjust,
          Count,
          geneID
        ) %>%

        separate_rows(
          geneID,
          sep = "/"
        ) %>%

        dplyr::rename(
          associated_gene = geneID
        )


      write.csv(
        go_gene_list,
        paste0(
          go_base,
          "_gene_list.csv"
        ),
        row.names = FALSE
      )


      # ======================================================
      # Original-style GO plot:
      # no redundancy collapsing
      # ======================================================

      pdf(
        paste0(
          go_base,
          "_dotplot.pdf"
        ),
        width = 11,
        height = 9
      )


      print(

        dotplot(
          go_res,
          showCategory = 20,
          label_format = 50
        ) +

          ggtitle(
            paste(
              "GO",
              ont,
              ":",
              geno,
              direction,
              "(SQANTI background)"
            )
          ) +

          theme_bw()
      )


      dev.off()


      cat(
        "[+] Saved",
        nrow(go_res@result),
        "GO",
        ont,
        "terms.\n"
      )
    }
  }
}


cat("\n")
cat("============================================================\n")
cat("SQANTI-BACKGROUND GO/KEGG ANALYSIS COMPLETE\n")
cat("============================================================\n")
