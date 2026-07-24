### Performing differential isoform expression analysis 
### Question - which RNA isoforms increase or decrease in expression during phylloxera infestation

# Note: Isoform count data files were generated at the last step of Iso-seq data processing
# The *classification.txt file was generated during the SQANTI3-based transcriptome filtering step
# The scripts/workflow for generating count files and classification files is shown in the 'IsoSeq_data_processing.txt' file


# Load R libraries

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(ggrepel)
library(ggplot2)

# Set working directory
setwd("/Users/natecollison/Desktop/isoseq_analyses_2")



# ========================================================
# 1. Load Data
# ========================================================

# SO4 isoform abundance (count) data, as an example

counts_file <- "SO4_counts_for_DE.csv" 
class_file  <- "SO4_filtered_RulesFilter_result_classification.txt"

# Read the counts matrix
counts_data <- read.csv(counts_file, header=TRUE, row.names=1, check.names=FALSE)

# Check data
head(counts_data, 5)

# Read the annotation file (to link PB.IDs (generic PacBio isoform IDs) to gene names from the Vitis reference genome)
class_data <- read.table(class_file, header=TRUE, sep="\t", quote=NULL) %>%
  select(isoform, associated_gene, structural_category)

head(class_data, 5)




# ========================================================
# 2. Create Metadata (Experimental Design)
# ========================================================
# We look at the column names to define who is "Control" and "Infection"
sample_names <- colnames(counts_data)
sample_names

# Create a simple table: SampleName -> Condition
meta_data <- data.frame(
  row.names = sample_names,
  condition = ifelse(grepl("con", sample_names), "Control", "Infested")
)

# Ensure "Control" is the baseline reference level
meta_data$condition <- factor(meta_data$condition, levels = c("Control", "Infested"))

print(meta_data) # Verify it looks correct!


# ========================================================
# 3. Run DESeq2 (With Normalized Filtering)
# ========================================================
# Create the DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = meta_data,
                              design = ~ condition)

# 1. Estimate Size Factors FIRST
# (We must do this now so DESeq2 knows how big each library is)
dds <- estimateSizeFactors(dds)

# Check the library sizes (optional, just to see)
print("Size Factors (Relative Library Depths):")
print(sizeFactors(dds))

# 2. Apply the Normalized Filter (CPM / FPM)
# Logic: Keep isoform if CPM > 0.5 in at least 3 samples
# (This accounts for sequencing depth differences!)
keep <- rowSums(fpm(dds) > 0.5) >= 3



dds <- dds[keep,]
print(paste("Remaining isoforms after CPM filtering:", sum(keep)))

#[1] "Remaining isoforms after CPM filtering: 202236" , SO4

dds <- DESeq(dds)

# Get Results
res <- results(dds)
summary(res)


##

# ========================================================
# 4. Filter and Save Results
# ========================================================

# Convert results to a standard dataframe
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = "isoform")


# --------------------------------------------------------
# Filter A: The "Resistance Candidates" (UP)
# Criteria: Adjusted P < 0.05 AND Log2FoldChange > 1 (Doubled in expression)
# --------------------------------------------------------
up_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  arrange(desc(log2FoldChange)) %>%
  left_join(class_data, by="isoform") # Add gene names

print(paste("Number of High-Confidence UPREGULATED isoforms:", nrow(up_genes)))
# [1] "Number of High-Confidence UPREGULATED isoforms: 1697". SO4

write.csv(up_genes, "S04_upreg_isoforms.csv", row.names=FALSE)


# --------------------------------------------------------
# Filter B: The "Metabolic Shutdown" (DOWN)
# Criteria: Adjusted P < 0.05 AND Log2FoldChange < -1 (Halved in expression)
# --------------------------------------------------------
down_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>%
  arrange(log2FoldChange) %>% # Sort by most negative
  left_join(class_data, by="isoform")

print(paste("Number of High-Confidence DOWNREGULATED isoforms:", nrow(down_genes)))
# [1] "Number of High-Confidence DOWNREGULATED isoforms: 15628", SO4

write.csv(down_genes, "S04_downreg_isoforms.csv", row.names=FALSE)





# ========================================================
# 5. PCA Plot (Quality Control)
# ========================================================

# 1. Transform the data for plotting
# 'blind=FALSE' means "don't hide the experimental design from the transformation"
# This is standard for Quality Control.
vsd <- vst(dds, blind=FALSE)

# 2. Basic PCA Plot (Quick Check)
plotPCA(vsd, intgroup="condition")

# 3. Fancier PCA Plot (Publication Quality using ggplot2)
# This lets you customize colors and labels easier
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=4) +                          # Big dots
  # geom_text(vjust=1.5, size=3) +              # Uncomment to see sample names on plot
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  ggtitle("PCA: SO4 Infested vs Control") +
  theme_bw()                                    # Clean white background




# print with sample names 
library(ggplot2)


# 1. Extract PCA data (using the vsd_all from the previous step)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# 2. Plot with Sample Names
ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=4) +
  geom_text_repel(size=3, max.overlaps = Inf) + # This adds the labels
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("SO4 PCA with Sample Labels") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))




#######

# Make a volcano plot
library(ggplot2)

# 1. Convert results to a Data Frame for plotting
volcano_data <- as.data.frame(res)

# 2. Add a column to classify genes (Significant vs Not)
# Criteria: P-adj < 0.05 AND |Log2FC| > 1 (2-fold change)
volcano_data$diffexpressed <- "NO"
volcano_data$diffexpressed[volcano_data$log2FoldChange > 1 & volcano_data$padj < 0.05] <- "UP"
volcano_data$diffexpressed[volcano_data$log2FoldChange < -1 & volcano_data$padj < 0.05] <- "DOWN"

# 3. Create the Plot
ggplot(volcano_data, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) +
  geom_point(alpha=0.6, size=1.5) +               # Semi-transparent dots
  theme_bw() +                                    # Clean white background
  geom_vline(xintercept=c(-1, 1), col="black", linetype="dashed") + # Threshold lines
  geom_hline(yintercept=-log10(0.05), col="black", linetype="dashed") + 
  scale_color_manual(values=c("blue", "grey", "red")) + # Blue=Down, Grey=NS, Red=Up
  labs(title="Volcano Plot: SO4 Infestastion vs Control",
       x="Log2 Fold Change (Infection / Control)",
       y="-Log10 Adjusted P-value") +
  theme(plot.title = element_text(hjust = 0.5))   # Center the title





