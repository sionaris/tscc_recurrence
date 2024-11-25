# This script focusses on importing the nextflow output in RStudio
# for preprocessing and exploratory analysis

# Libraries #####
library(dplyr)
library(DESeq2)
library(tximport)
library(tximeta)
library(readr)
library(stringr)
library(edgeR)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(cowplot)
library(pheatmap)
library(rcartocolor)
library(viridisLite)
library(vsn)
library(RColorBrewer)
library(tibble)
library(ggrepel)
library(glmpca)
library(openxlsx)
library(PoiClaClu)
library(EnhancedVolcano)
library(genefilter)
library(grid)
library(openxlsx)

# Preprocessing #####

# The nextflow output is all included in a gitignored (due to size limitations)
# directory in the repository, called "NF-processed_samples". Within this directory,
# there are 80 subdirectories, each specific to a sample, named "batch_*", where * 
# ranges from 0 to 79. From each subdirectory we fetch the salmon.merged.gene_counts file

# Import the files into R
# Load pheno data
Pheno = read.xlsx("Resources_and_environments/Pheno.xlsx")

# Extract sample names in order
batch_dirs <- paste0("NF-processed_samples/batch_", seq(0, 79, 1))

# Extract sample names from COHT?-* subdirectories in order
samplenames <- c()
for (batch in batch_dirs) {
  star_salmon_dir <- file.path(batch, "star_salmon")
  
  if (dir.exists(star_salmon_dir)) {
    subdirs <- list.dirs(star_salmon_dir, recursive = FALSE, full.names = FALSE)
    sample <- grep("^COHT.*-.*$", subdirs, value = TRUE)
    
    if (length(sample) > 0) {
      samplenames <- c(samplenames, sample[1])
    } else {
      samplenames <- c(samplenames, NA)
    }
  } else {
    samplenames <- c(samplenames, NA)
  }
}

rm(sample, subdirs, star_salmon_dir, batch); gc()

# Get paths to count files
count_files <- paste0(batch_dirs, "/star_salmon/", samplenames, "/quant.sf")

# Start building pheno data
coldata = as.data.frame(list(Directory = batch_dirs, Patient.ID = samplenames, files = count_files)) %>%
  arrange(Patient.ID) %>%
  mutate(names = paste0(Patient.ID, c("_a", "_b"))) %>%
  mutate(Recurrence = c(rep("Yes", 40), rep("No", 40)))

# Get paths to tx2gene file
tx2gene <- readr::read_tsv(paste0(batch_dirs[1], "/star_salmon/tx2gene.tsv"), col_names = FALSE)[, c(1, 3)]
txi = tximport(coldata$files, type = "salmon", tx2gene = tx2gene)
colnames(txi[[1]]) = colnames(txi[[2]]) = colnames(txi[[3]]) = coldata$names

# Create DESeq object using tximport output
dds_pre = DESeqDataSetFromTximport(txi = txi,
                               colData = coldata,
                               design = ~ Recurrence)

# How many genes in this dataset are all zeros?
table(rowSums(dds_pre@assays@data@listData[["counts"]]==0)==80)

# FALSE  TRUE 
# 57925   4585

print(paste0(round(100*4585/(4585+57925), 2), "% of genes only have zero counts."))
# 7.33% of genes only have zero counts

# We remove these genes from the experiment
toremove = which(rowSums(dds_pre@assays@data@listData[["counts"]]==0)==80)
dds_pre = dds_pre[-toremove, ]

rm(toremove); gc()

# Generate GPCA and MDS plots to check how similar samples of the same patient
# are to each other

# glmpca
RNGversion("4.2.2")
set.seed(123)
gpca = glmpca(counts(dds_pre), L=2)
gpca.dat = gpca$factors
gpca.dat$Sample.ID = rownames(gpca.dat)

gpcaplot = ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Sample.ID)) +
  geom_point(size = 0.2, alpha = 0.65) +
  ggtitle("Generalized Principal Component Analysis (glmpca, counts)") +
  theme_classic() +
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 2),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank()) +
  guides(color = "none", alpha = "none")

# MDS with Poisson distances
samplePoisDistMatrix_pre = as.matrix(PoissonDistance(t(counts(dds_pre)))$dd)
dimnames(samplePoisDistMatrix_pre) = list(colnames(counts(dds_pre)), colnames(counts(dds_pre)))

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix_pre))
poisson_mdsplot = ggplot(mdsPois, aes(x = `1`, y = `2`, color = names)) +
  geom_point(size = 0.2, alpha = 0.65) +
  ggtitle("Multidimensional scaling: Poisson distance (counts)") +
  theme_classic() +
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "none",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank())+
  labs(x = "MDS1", y = "MDS2") +
  guides(alpha = "none", color = "none")

tiff("Results/Exploratory analysis/dimred_prior_to_sum.tiff", width = 2980, height = 1380, res = 700,
     compression = "lzw")
ggarrange(gpcaplot, poisson_mdsplot,
          ncol = 2, nrow = 1, common.legend = TRUE, 
          labels = c("A", "B"), legend = "right",
          font.label = list(size = 4, color = "black", face = "bold"))
dev.off()

# Samples from the same patient cluster together
# We will therefore sum these samples and reduce the sample size from 80 -> 40
# Normalization will take this into account

# Here are some resources which justify this approach:
# https://support.bioconductor.org/p/85536/
# https://www.seqanswers.com/forum/bioinformatics/bioinformatics-aa/48596-why-sum-technical-replicates-instead-of-averaging-and-rounding-for-deseq?t=60996#4
# https://www.biostars.org/p/273421/
# https://www.biostars.org/p/9542505/#9542536

colnames(counts(dds_pre))

# Extract unique sample IDs by removing the suffix "_a" and "_b"
unique_sample_ids <- unique(sub("_[ab]$", "", colnames(counts(dds_pre))))

# Initialize a matrix to store summed counts
summed_counts <- sapply(unique_sample_ids, function(sample_id) {
  # Find the columns corresponding to this sample ID
  cols_to_sum <- grep(paste0("^", sample_id, "_[ab]$"), 
                      colnames(counts(dds_pre)), value = TRUE)
  
  # Sum the columns
  rowSums(counts(dds_pre)[, cols_to_sum, drop = FALSE])
})

# Convert summed_counts to a matrix and set proper row and column names
summed_counts <- as.matrix(summed_counts)
rownames(summed_counts) <- rownames(counts(dds_pre))
colnames(summed_counts) <- unique_sample_ids

coldata <- as.data.frame(colData(dds_pre))
coldata <- coldata[match(unique_sample_ids, sub("_[ab]$", "", rownames(coldata))), , drop = FALSE]
rownames(coldata) <- unique_sample_ids

# Create a full design matrix for the filtering
design = model.matrix(~ Pheno$Sample.type + Pheno$Stage + Pheno$Recurrence + 
                        Pheno$Sample.type:Pheno$Recurrence +
                        Pheno$Stage:Pheno$Recurrence)

# We need to filter out uninformative genes
keep.exprs = names(which(filterByExpr(summed_counts, design = design) == TRUE))

# This design results in ~20,000 genes less (57925 -> 38992)
summed_counts = summed_counts[keep.exprs, ]

# Create a new DESeqDataSet object from the summed counts matrix
dds <- DESeqDataSetFromMatrix(countData = summed_counts,
                              colData = Pheno,
                              design = design)

rm(unique_sample_ids, summed_counts, dds_pre); gc()

dds = estimateSizeFactors(dds)

# Density curve function# Density curdds_preve function
create_density_curve <- function(matrix) {
  as.data.frame(matrix) %>%
    gather(.) %>%
    ggplot(aes(x = value, fill = "Density Curve")) + 
    geom_density(alpha = 0.5, position = "identity") +
    xlim(-1, NA) +
    labs(title = "Distribution of values by sample",
         x = "Value", y = "Density") +
    theme_classic()+
    theme(legend.position = "none")
}

# Per column density plot
create_density_plot_color <- function(matrix) {
  as.data.frame(matrix) %>%
    gather(., key = "sample", value = "value") %>%
    ggplot(aes(x = value, fill = sample)) + 
    geom_density(alpha = 0.5, position = "identity") +
    xlim(-1, NA) +
    labs(title = "Distribution of values by sample",
         x = "Value", y = "Density") +
    theme_classic()+
    theme(legend.position = "none")
}

# Density plots/curves for our data
create_density_curve(counts(dds)) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "OSCC_filtered_read_counts_density_curves.tiff",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# The previous plot contains count outliers (abnormally high counts) which make 
# it difficult to see the distribution of counts near zero. We replot with 
# xlim 30

create_density_curve(counts(dds)) +
  scale_x_continuous(limits = c(-1, 30), breaks = seq(0, 30, 10)) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "OSCC_filtered_read_counts_density_curves_xlim30.tiff",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# Color density plots for counts
create_density_plot_color(counts(dds)) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "OSCC_filtered_read_counts_color_densities.tiff",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# With xlim 30
create_density_plot_color(counts(dds))+
  scale_x_continuous(limits = c(-1, 30), breaks = seq(0, 30, 10)) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "OSCC_filtered_read_counts_color_densities_xlim30.tiff",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# Now with normalized counts

# Firstly, transform the normalized counts with a data-driven pseudocount
# https://www.biorxiv.org/content/10.1101/404962v1.full -> Process description

# We need a grouping for the OSCC samples. We are going to go for 
# Recurrence vs No recurrence
rec_samples = Pheno$Patient.ID[Pheno$Recurrence == "Recurrence"]
norec_samples = Pheno$Patient.ID[Pheno$Recurrence == "No recurrence"]

# Get size factors similarly to TCGA
# Mock sample conditions
sample_conditions <- data.frame(
  row.names = colnames(counts(dds)),
  condition = factor(rep("condition", ncol(counts(dds))))
)

# Create DESeqDataSet object
dds_null <- DESeqDataSetFromMatrix(countData = counts(dds), 
                              colData = sample_conditions, 
                              design = ~ 1) 

# Estimate size factors
dds_null <- estimateSizeFactors(dds_null)

# Save size factors for later if needed
size_factors = dds_null$sizeFactor
names(size_factors) = dds_null@colData@rownames

# Get the normalized counts
normalized_counts <- counts(dds_null, normalized = TRUE)

# Assign proper row and column names to the normalized counts
rownames(normalized_counts) = rownames(counts(dds))
colnames(normalized_counts) = colnames(counts(dds))

# Estimate averages
avg_rec_size_factor = mean(size_factors[rec_samples], na.rm = TRUE)
avg_norec_size_factor = mean(size_factors[norec_samples], na.rm = TRUE)
sf = c(avg_rec_size_factor, avg_norec_size_factor)

# Lun et al. suggest pseudocount = max{1, r|1/smin - 1/smax|}, where r = 1 (suggestion)
pseudocount_RNA = max(c(1, abs(1/min(sf) - 1/max(sf)))) # 1

# We proceed with the log2(norm.counts + 1) transformation
log2.norm.counts.plus1 = log2(normalized_counts + pseudocount_RNA)

rm(rec_samples, norec_samples, sample_conditions, normalized_counts, dds_null, 
   size_factors, avg_norec_size_factor, avg_rec_size_factor, sf, pseudocount_RNA); gc()

create_density_curve(log2.norm.counts.plus1) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "OSCC_filtered_log2normcounts_density_curves.tiff",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# Color density plots for counts
create_density_plot_color(log2.norm.counts.plus1) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "OSCC_filtered_log2normcounts_color_densities.tiff",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'tiff', units = "px",
       dpi = 150, compression = "lzw")
dev.off()

# Further exploratory analysis #####
# VST and rlog transformations
vsd = vst(dds, blind = FALSE)
rld = rlog(dds, blind = FALSE)

# See https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-transformations-and-visualization
# for why blind = FALSE

# Produce a scatter plot of transformed counts for the first two samples
df = bind_rows(
  as.data.frame(log2.norm.counts.plus1[, 1:2]) %>%
    mutate(transformation = "log2(norm.counts + 1)"),
  as.data.frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as.data.frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls = c("log2(norm.counts + 1)", "vst", "rlog")
df$transformation = factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

# Produce mean - sd plots 
plot1 = meanSdPlot(log2.norm.counts.plus1)$gg + rremove("ylab") + rremove("xlab") + 
  ggtitle(bquote(bold("log2(norm.counts + 1)"))) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
plot2 = meanSdPlot(assay(rld))$gg + rremove("ylab") + rremove("xlab") + 
  ggtitle(bquote(bold("rlog"))) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
plot3 = meanSdPlot(assay(vsd))$gg + rremove("ylab") + rremove("xlab") + 
  ggtitle(bquote(bold("vst"))) +
  theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
dev.off()

# Use ggarrange to align the three plots horizontally and set the x-axis and y-axis labels
fig = ggarrange(plot1, plot2, plot3, 
                ncol = 3, nrow = 1, 
                common.legend = TRUE, legend = "bottom",
                labels = c("A", "B", "C"),
                align = "hv",
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
tiff("Results/Exploratory analysis/meanSdPlot.tiff", res = 300, width = 4500, height = 1600, compression = "lzw")
fig
annotate_figure(fig, left = textGrob(bquote(bold("sd")), rot = 90, vjust = 1, gp = gpar(cex = 1)),
                bottom = textGrob(bquote(bold("rank (mean)")), gp = gpar(cex = 1)))
dev.off()

# Heatmap of the count matrix - top 50 genes
select = order(rowMeans(log2.norm.counts.plus1),
               decreasing=TRUE)[1:50]
df = as.data.frame(colData(dds)[,c("Sample.type", "Stage", "Recurrence")])

# Annotation colors
ann_colors = list(
  Recurrence = c(`No recurrence` = "dodgerblue4", Recurrence = "deeppink4"),
  Sample.type = c(Tumor = "purple4", Normal = "orange"),
  Stage = c(Early = carto_pal(n = 7, "Emrld")[4], 
                Late = carto_pal(n = 7, "Burg")[4])
)

# log2(norm.counts + 1)
tiff("Results/Exploratory analysis/top50_heatmap_log2(norm.counts+1).tiff",
     res = 200, width = 1420, 
     height = 1080, compression = "lzw")
pheatmap(log2.norm.counts.plus1[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "log2(counts + 1)",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

# rlog
tiff("Results/Exploratory analysis/top50_heatmap_rlog.tiff",
     res = 200, width = 1420, 
     height = 1080, compression = "lzw")
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "rlog",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

# vst
tiff("Results/Exploratory analysis/top50_heatmap_vst.tiff", res = 200, width = 1420, 
     height = 1080, compression = "lzw")
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "vst",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

# We will only use the vst-transformed data for the next visualisations

# Heatmap of sample distances
sampleVSTdists = dist(t(assay(vsd)))
sampleDistMatrix = as.matrix(sampleVSTdists)
#colnames(sampleDistMatrix) = NULL
colors = colorRampPalette(viridisLite::magma(10))(255)

tiff("Results/Exploratory analysis/sample_distance_heatmap_vst.tiff", 
     res = 700, width = 6760, 
     height = 4760, compression = "lzw")
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleVSTdists,
         clustering_distance_cols = sampleVSTdists,
         col = colors,
         fontsize = 8,
         main = "Sample distance heatmap",
         annotation_col = df,
         annotation_colors = ann_colors)
dev.off()

# Heatmap of sample Poisson distances (more appropriate for count data)
poisd = PoissonDistance(t(counts(dds)))
samplePoisDistMatrix = as.matrix(poisd$dd)
dimnames(samplePoisDistMatrix) = dimnames(sampleDistMatrix)
#colnames(samplePoisDistMatrix) = NULL

tiff("Results/Exploratory analysis/sample_poisson_distance_heatmap.tiff", res = 700, width = 6760, 
     height = 4760, compression = "lzw")
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         fontsize = 8,
         main = "Sample Poisson distance heatmap",
         annotation_col = df,
         annotation_colors = ann_colors)
dev.off()

# Dimensionality reduction #####
# Principal Component Analysis
pcaData = plotPCA(vsd, intgroup = c("Recurrence", "Sample.type"), 
                  returnData=TRUE) %>%
  dplyr::rename(Patient.ID = name)
stage = as.data.frame(vsd@colData@listData) %>%
  rownames_to_column %>%
  dplyr::select(Patient.ID = rowname, Stage)
pcaData = pcaData %>%
  inner_join(stage, by = "Patient.ID")
percentVar = round(100 * attr(pcaData, "percentVar"))
pcaplot = ggplot(pcaData, aes(PC1, PC2, color = Sample.type, shape = Recurrence)) +
  geom_point(size = 0.2, alpha = 0.65) +
  geom_text_repel(aes(label = Stage), segment.linetype = 1,
                  segment.color = "black", segment.size = 0.1,
                  min.segment.length = 0,
                  size = 0.5, nudge_y = 0.03, nudge_x = 0.03, show.legend = FALSE) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("Principal Component Analysis: VST") +
  theme_classic() +
  scale_color_manual(name = "Sample.type", values = c("orange", "purple4"),
                     labels = c("Normal", "Tumor"))+
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank())

# Generalized PCA (more appropriate for count data)
RNGversion("4.2.2")
set.seed(123)
gpca = glmpca(counts(dds), L=2)
gpca.dat = gpca$factors
gpca.dat$Sample.type = dds$Sample.type
gpca.dat$Stage = dds$Stage
gpca.dat$Recurrence = dds$Recurrence

gpcaplot = ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Sample.type, shape = Recurrence)) +
  geom_point(size = 0.2, alpha = 0.65) +
  geom_text_repel(aes(label = Stage), segment.linetype = 1,
                  segment.color = "black", segment.size = 0.1,
                  min.segment.length = 0,
                  size = 0.5, nudge_y = 0.03, nudge_x = 0.03, show.legend = FALSE) +
  ggtitle("Generalized Principal Component Analysis (glmpca, counts)") +
  theme_classic() +
  scale_color_manual(name = "Sample.type", values = c("orange", "purple4"),
                     labels = c("Normal", "Tumor"))+
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 2),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank())

tiff("Results/Exploratory analysis/PCA.tiff", width = 2880, height = 1380, res = 700,
     compression = "lzw")
ggarrange(pcaplot, gpcaplot, 
          ncol = 2, nrow = 1, common.legend = TRUE, 
          labels = c("A", "B"), legend = "right",
          font.label = list(size = 4, color = "black", face = "bold"))
dev.off()

# Multidimensional scaling plots for VST and Poisson distances

# VST
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
mdsplot = ggplot(mds, aes(x = `1`, y = `2`, color = Sample.type, shape = Recurrence)) +
  geom_point(size = 0.2, alpha = 0.65) +
  geom_text_repel(aes(label = Stage, segment.linetype = 1,
                      segment.color = "black", segment.size = 0.1),
                  min.segment.length = 0,
                  size = 0.5, nudge_y = 0.03, nudge_x = 0.03) +
  ggtitle("Multidimensional scaling: VST") +
  theme_classic() +
  scale_color_manual(name = "Sample.type", values = c("orange", "purple4"),
                     labels = c("Normal", "Tumor"))+
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank())+
  labs(x = "MDS1", y = "MDS2") +
  guides(alpha = "none")

# Poisson
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
poisson_mdsplot = ggplot(mdsPois, aes(x = `1`, y = `2`, color = Sample.type, shape = Recurrence)) +
  geom_point(size = 0.2, alpha = 0.65) +
  geom_text_repel(aes(label = Stage, segment.linetype = 1,
                      segment.color = "black", segment.size = 0.1),
                  min.segment.length = 0,
                  size = 0.5, nudge_y = 0.03, nudge_x = 0.03) +
  ggtitle("Multidimensional scaling: Poisson distance (counts)") +
  theme_classic() +
  scale_color_manual(name = "Sample.type", values = c("orange", "purple4"),
                     labels = c("Normal", "Tumor"))+
  theme(plot.title = element_text(size = 4, face = "bold", vjust = 0.5, hjust = 0.5),
        axis.text = element_text(size = 3, hjust = 0.5, vjust = 0.5, 
                                 color = "black"),
        axis.title = element_text(size = 4, face = "bold"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.line = element_line(linewidth = 0.3),
        legend.position = "right",
        legend.key.size = unit(2, units = "mm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(face = "bold", size = 3.5),
        legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
        legend.spacing.y = unit(0.5, units = "mm"),
        legend.spacing.x = unit(0.5, units = "mm"),
        legend.background = element_blank())+
  labs(x = "MDS1", y = "MDS2") +
  guides(alpha = "none")

fig = ggarrange(mdsplot +rremove("ylab"), poisson_mdsplot +rremove("ylab"), 
                ncol = 2, nrow = 1, 
                common.legend = TRUE, legend = "right",
                labels = c("A", "B"),
                align = "hv",
                font.label = list(size =4, color = "black", face = "bold"))
tiff("Results/Exploratory analysis/MDS.tiff", res = 700, width = 2880, height = 1080, compression = "lzw")
fig
annotate_figure(fig, left = textGrob(bquote(bold("MDS2")), rot = 90, vjust = 1, gp = gpar(cex = 0.45)))
dev.off()

# Differential Expression #####
# We are going to use a LRT as described in the DESeq2 vignette:
# https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments
# and the RNAseq workflow by Michael Love:
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments

dds = DESeq(dds, test = "LRT", reduced = ~ batch + treatment + timepoint)
resultsNames(dds) # Names of the coefficients

contrasts = list(HMBA24h_vs_Control48h = list("treatment_HMBA_vs_control",
                                              "timepoint_48h_vs_24h"),
                 HMBA48h_vs_Control48h = list(c("treatment_HMBA_vs_control", 
                                                "treatmentHMBA.timepoint48h")),
                 HMBA72h_vs_Control48h = list(c("treatment_HMBA_vs_control",
                                                "timepoint_72h_vs_24h",
                                                "treatmentHMBA.timepoint72h"),
                                              "timepoint_48h_vs_24h"),
                 HMBA48h_vs_HMBA24h = list(c("timepoint_48h_vs_24h",
                                             "treatmentHMBA.timepoint48h")),
                 HMBA72h_vs_HMBA24h = list(c("treatmentHMBA.timepoint72h",
                                             "timepoint_72h_vs_24h")),
                 HMBA72h_vs_HMBA48h = list(c("treatmentHMBA.timepoint72h",
                                             "timepoint_72h_vs_24h"), 
                                           c("timepoint_48h_vs_24h",
                                             "treatmentHMBA.timepoint48h")))

# hmba24h vs control48h
res_HMBA24h_vs_control48h <- results(dds, contrast=contrasts[["HMBA24h_vs_Control48h"]], 
                                     listValues=c(1, -1), test = "Wald")

# hmba48h vs control48h
res_HMBA48h_vs_control48h <- results(dds, 
                                     contrast=contrasts[["HMBA48h_vs_Control48h"]], 
                                     listValues=c(1, -1), test = "Wald")

# hmba72h vs control48h
res_HMBA72h_vs_control48h <- results(dds, contrast=contrasts[["HMBA72h_vs_Control48h"]],
                                     listValues=c(1, -1), test = "Wald")

# HMBA48h vs HMBA24h
res_HMBA48h_vs_HMBA24h <- results(dds, 
                                  contrast=contrasts[["HMBA48h_vs_HMBA24h"]],
                                  listValues=c(1, -1), test = "Wald")

# HMBA72h vs HMBA24h
res_HMBA72h_vs_HMBA24h <- results(dds, 
                                  contrast=contrasts[["HMBA72h_vs_HMBA24h"]], 
                                  listValues=c(1, -1), test = "Wald")

# HMBA72h vs HMBA48h
res_HMBA72h_vs_HMBA48h <- results(dds, contrast=contrasts[["HMBA72h_vs_HMBA48h"]],
                                  listValues=c(1, -1), test = "Wald")

results = list(res_HMBA24h_vs_control48h, res_HMBA48h_vs_control48h, 
               res_HMBA72h_vs_control48h, res_HMBA48h_vs_HMBA24h,
               res_HMBA72h_vs_HMBA24h, res_HMBA72h_vs_HMBA48h)
names(results) = names(contrasts)

# Shrinking logFCs
# apeglm_shrunk_results = list() # apeglm does not work with contrasts
# normal_shrunk_results = list() # normal cannot be implemented with interactions designs
ashr_shrunk_results = list()

for (i in 1:length(results)) {
  # apeglm_shrunk_results[[i]] = lfcShrink(dds, contrast = contrasts[[i]], type = "apeglm")
  # normal_shrunk_results[[i]] = lfcShrink(dds, contrast = contrasts[[i]], type = "normal")
  ashr_shrunk_results[[i]] = lfcShrink(dds, res = results[[i]], type = "ashr")
}
names(ashr_shrunk_results) = names(results)

# Create a workbook and save data there to export later
library(openxlsx)
dgea_wb = createWorkbook()
dgea_wb_shrunk = createWorkbook()

# Create header and body styles for the output Excel Sheet
header_style = createStyle(textDecoration = "bold")
body_style = createStyle(halign = "left", valign = "center")

# Choose if you want the colored gene symbols or not (color coding is time-consuming)
color_handle = FALSE # Change that to FALSE if you don't want colored symbols

if (color_handle) {
  blue_style = createStyle(fontColour = "#279AE5", textDecoration = "bold")  # Blue and bold
  red_style = createStyle(fontColour = "#B91B60", textDecoration = "bold")  # Red and bold
  grey_style = createStyle(fontColour = "#818383", textDecoration = "bold")  # Grey and bold
  
  # Define a function to determine the color style based on the conditions
  getStyle <- function(log2foldchange, padj) {
    if (is.na(log2foldchange) | is.na(padj)) {
      return(grey_style)
    } else if (log2foldchange < 0 & padj < 0.05) {
      return(blue_style)
    } else if (log2foldchange > 0 & padj < 0.05) {
      return(red_style)
    } else {
      return(grey_style)
    }
  }
}

# Create data frames that will be exported into .xlsx files
for (i in 1:length(results)) {
  # Raw LFCs
  genes = rownames(results[[i]])
  res = results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj) %>%
    dplyr::select(Gene.Symbol, everything())
  
  addWorksheet(dgea_wb, names(results)[i])
  writeData(dgea_wb, names(results)[i], res)
  
  # Apply the styles
  addStyle(dgea_wb, names(results)[i], style = header_style, rows = 1, 
           cols = 1:ncol(res), gridExpand = TRUE) # first row only
  addStyle(dgea_wb, names(results)[i], style = body_style, rows = 2:nrow(res), 
           cols = 1:ncol(res), gridExpand = TRUE)
  
  if (color_handle) {
    color_styles = mapply(getStyle, res$log2FoldChange, res$padj)
    mapply(addStyle, MoreArgs = list(wb = dgea_wb, sheet = names(results)[i], 
                                     cols = which(colnames(res) == "Gene.Symbol")), 
           style = color_styles, rows = 2:(nrow(res)+1))
  }
  
  # ASHR LFCs
  genes = rownames(ashr_shrunk_results[[i]])
  res_ashr = ashr_shrunk_results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj) %>%
    dplyr::select(Gene.Symbol, everything())
  
  addWorksheet(dgea_wb_shrunk, names(ashr_shrunk_results)[i])
  writeData(dgea_wb_shrunk, names(ashr_shrunk_results)[i], res_ashr)
  
  # Apply the styles
  addStyle(dgea_wb_shrunk, names(ashr_shrunk_results)[i], style = header_style, 
           rows = 1, cols = 1:ncol(res_ashr), gridExpand = TRUE) # first row only
  addStyle(dgea_wb_shrunk, names(ashr_shrunk_results)[i], style = body_style, 
           rows = 2:nrow(res_ashr), cols = 1:ncol(res_ashr), gridExpand = TRUE)
  
  if (color_handle) {
    color_styles = mapply(getStyle, res_ashr$log2FoldChange, res_ashr$padj)
    mapply(addStyle, MoreArgs = list(wb = dgea_wb_shrunk, 
                                     sheet = names(ashr_shrunk_results)[i], 
                                     cols = which(colnames(res_ashr) == "Gene.Symbol")), 
           style = color_styles, rows = 2:(nrow(res_ashr)+1))
  }
}

# Create DGEA folder
if (!dir.exists("DGEA")) {
  dir.create("DGEA")
}

saveWorkbook(dgea_wb, "DGEA/DGEA_raw.xlsx",
             overwrite = TRUE)
saveWorkbook(dgea_wb_shrunk, "DGEA/DGEA_ashr.xlsx",
             overwrite = TRUE)
rm(genes, res, res_ashr, i); gc()

# MA-plots
if (!dir.exists("DGEA/MA_plots")) {
  dir.create("DGEA/MA_plots")
}

for (i in 1:length(results)) {
  tiff(paste0("DGEA/MA_plots/MA_plot_",
              names(results)[i], ".tiff"), res = 700,
       width = 8000, height = 4000, compression = "lzw")
  par(mfrow=c(1,2), mar=c(4,4,2,1))
  print(DESeq2::plotMA(results[[i]], main = "MA plot: raw LFCs",
                       xlab = substitute(paste(bold('Mean of normalized counts'))),
                       ylab = substitute(paste(bold('log fold change'))),
                       alpha = 0.05))
  print(DESeq2::plotMA(ashr_shrunk_results[[i]], main = "MA plot: ashr shrunk LFCs",
                       xlab = substitute(paste(bold('Mean of normalized counts'))),
                       ylab = substitute(paste(bold('log fold change'))),
                       alpha = 0.05))
  dev.off()
}

# Volcano plots
library(EnhancedVolcano)

if (!dir.exists("DGEA/Volcano_plots")) {
  dir.create("DGEA/Volcano_plots")
}

# Volcano loop (raw LFCs)
for (i in 1:length(results)) {
  # Construct the data frame
  genes = rownames(results[[i]])
  res = results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj)
  
  # create custom key-value pairs for stat. sig genes (p.adj < 0.05) and n.s genes
  keyvals.colour <- ifelse(
    res$log2FoldChange < -1 & res$padj < 0.05, 'royalblue',
    ifelse(res$log2FoldChange > 1 & res$padj < 0.05, 'red4',
           ifelse(abs(res$log2FoldChange) < 1 & res$padj < 0.05, 'pink',
                  'grey')))
  
  # keyvals.colour[is.na(keyvals.colour)] <- 'black'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'
  names(keyvals.colour)[keyvals.colour == 'red4'] <- 'Up-regulated'
  names(keyvals.colour)[keyvals.colour == 'pink'] <- '|log2FC| < 1'
  names(keyvals.colour)[keyvals.colour == 'grey'] <- 'p.adj > 0.05'
  
  res$aes = keyvals.colour
  
  volcano = EnhancedVolcano(res,
                            lab = res[, "Gene.Symbol"],
                            caption = NULL,
                            x = 'log2FoldChange',
                            y = 'padj',
                            pCutoff = 0.05,
                            FCcutoff = 1,
                            cutoffLineType = "dashed",
                            cutoffLineWidth = 0.3,
                            cutoffLineCol = "black",
                            colCustom = keyvals.colour,
                            colAlpha = 0.7,
                            ylab = bquote(bold(-log[10]("BH adj. p-value"))),
                            xlab = bquote(bold(log[2]("Fold Change"))),
                            pointSize = 1,
                            axisLabSize = 7,
                            subtitle = NULL,
                            labSize = 2,
                            selectLab = res[1:20, "Gene.Symbol"],
                            legendLabSize = 8,
                            legendIconSize = 4,
                            title = paste0("Volcano plot: contrast = ",
                                           names(results)[i]),
                            labFace = "bold",
                            boxedLabels = TRUE,
                            drawConnectors = TRUE,
                            typeConnectors = "closed",
                            arrowheads = FALSE,
                            widthConnectors = 0.3,
                            max.overlaps = Inf) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 10, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(face = "bold", size = 8),
          axis.line = element_line(colour = "black", linewidth = 0.4),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.ticks.length = unit(1, units = "mm"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.margin = ggplot2::margin(1, 1, 1, 1, unit = "mm"),
          legend.spacing.y = unit(1, units = "mm"),
          legend.spacing.x = unit(1, units = "mm"),
          legend.background = element_blank(),
          legend.key.size = unit(0.5, units = "mm"),
          legend.text = element_text(size = 6),
          legend.box.background = element_rect(colour = "black"))
  
  volcano
  ggsave(filename = paste0("Volcano_raw_LFCs_", names(results)[i], ".tiff"),
         path = "DGEA/Volcano_plots", 
         width = 152, height = 142, device = 'tiff', units = "mm",
         dpi = 700, compression = "lzw")
}
rm(genes, i, res, keyvals.colour, volcano)

# Volcano loop (ashr-shrunk LFCs)
for (i in 1:length(ashr_shrunk_results)) {
  # Construct the data frame
  genes = rownames(ashr_shrunk_results[[i]])
  res = ashr_shrunk_results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj)
  
  # create custom key-value pairs for stat. sig genes (p.adj < 0.05) and n.s genes
  keyvals.colour <- ifelse(
    res$log2FoldChange < -1 & res$padj < 0.05, 'royalblue',
    ifelse(res$log2FoldChange > 1 & res$padj < 0.05, 'red4',
           ifelse(abs(res$log2FoldChange) < 1 & res$padj < 0.05, 'pink',
                  'grey')))
  
  # keyvals.colour[is.na(keyvals.colour)] <- 'black'
  names(keyvals.colour)[keyvals.colour == 'royalblue'] <- 'Down-regulated'
  names(keyvals.colour)[keyvals.colour == 'red4'] <- 'Up-regulated'
  names(keyvals.colour)[keyvals.colour == 'pink'] <- '|log2FC| < 1'
  names(keyvals.colour)[keyvals.colour == 'grey'] <- 'p.adj > 0.05'
  
  res$aes = keyvals.colour
  
  volcano = EnhancedVolcano(res,
                            lab = res[, "Gene.Symbol"],
                            caption = NULL,
                            x = 'log2FoldChange',
                            y = 'padj',
                            pCutoff = 0.05,
                            FCcutoff = 1,
                            cutoffLineType = "dashed",
                            cutoffLineWidth = 0.3,
                            cutoffLineCol = "black",
                            colCustom = keyvals.colour,
                            colAlpha = 0.7,
                            ylab = bquote(bold(-log[10]("BH adj. p-value"))),
                            xlab = bquote(bold(log[2]("Fold Change"))),
                            pointSize = 1,
                            axisLabSize = 7,
                            subtitle = NULL,
                            labSize = 2,
                            selectLab = res[1:20, "Gene.Symbol"],
                            legendLabSize = 8,
                            legendIconSize = 4,
                            title = paste0("Volcano plot: contrast = ",
                                           names(results)[i]),
                            labFace = "bold",
                            boxedLabels = TRUE,
                            drawConnectors = TRUE,
                            typeConnectors = "closed",
                            arrowheads = FALSE,
                            widthConnectors = 0.3,
                            max.overlaps = Inf) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.4),
          plot.title = element_text(face = "bold", size = 10, hjust = 0.5, vjust = 0.5),
          axis.title = element_text(face = "bold", size = 8),
          axis.line = element_line(colour = "black", linewidth = 0.4),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.ticks.length = unit(1, units = "mm"),
          legend.position = "right",
          legend.title = element_blank(),
          legend.margin = ggplot2::margin(1, 1, 1, 1, unit = "mm"),
          legend.spacing.y = unit(1, units = "mm"),
          legend.spacing.x = unit(1, units = "mm"),
          legend.background = element_blank(),
          legend.key.size = unit(0.5, units = "mm"),
          legend.text = element_text(size = 6),
          legend.box.background = element_rect(colour = "black"))
  
  volcano
  ggsave(filename = paste0("Volcano_ashr_shrunk_LFCs_", names(ashr_shrunk_results)[i], ".tiff"),
         path = "DGEA/Volcano_plots", 
         width = 152, height = 142, device = 'tiff', units = "mm",
         dpi = 700, compression = "lzw")
}
rm(genes, i, res, keyvals.colour, volcano)

# Gene clustering with mean-centered values from the VST data
library(genefilter)

# Create output directories
if (!dir.exists("DGEA/Post_DGEA_heatmaps")) {
  dir.create("DGEA/Post_DGEA_heatmaps")
}

# Subdirectories
if (!dir.exists("DGEA/Post_DGEA_heatmaps/raw_lfcs")) {
  dir.create("DGEA/Post_DGEA_heatmaps/raw_lfcs")
}

if (!dir.exists("DGEA/Post_DGEA_heatmaps/ashr_shrunk")) {
  dir.create("DGEA/Post_DGEA_heatmaps/ashr_shrunk")
}

# Number of features for plots
n_degs = 20

# Top n_degs variable genes fromt the VST object
topVarGenes = head(order(rowVars(assay(vsd)), decreasing = TRUE), n_degs)
matvar  = assay(vsd)[topVarGenes, ]
matvar  = matvar - rowMeans(matvar)
colnames(matvar) = colnames(assay(vsd))
anno = as.data.frame(colData(vsd)[, c("batch", "treatment", "timepoint")])

# Create a custom palette with 20*14 = 280 color breaks (one for each tile)
colors2 = colorRampPalette(c(viridisLite::mako(280)[140:280],
                             rev(viridisLite::rocket(280)[140:280])))(280)

tiff(paste0("DGEA/Post_DGEA_heatmaps/VST_top20vargenes_samples_heatmap.tiff"),
     res = 700, width = 6760, 
     height = 4760, compression = "lzw")
pheatmap(matvar, annotation_col = anno,
         col = colors2,
         fontsize = 8,
         cutree_col = 2,
         main = "Top 20 variable VST genes (centered) heatmap",
         annotation_colors = ann_colors,
         cluster_rows = FALSE)
dev.off()

# Heatmaps for raw LFCs
for (i in 1:length(results)) {
  # Construct the data frame
  genes = rownames(results[[i]])
  res = results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj)
  
  topDEGs = res[which(res$padj < 0.05),]
  
  # Gene clustering with raw counts and Poisson distance
  topDEGs = topDEGs[1:n_degs, ]
  matcount = assay(dds)[topDEGs$Gene.Symbol,]
  DEGpoisd = PoissonDistance(t(counts(dds[topDEGs$Gene.Symbol,])))
  DEGsamplePoisDistMatrix = as.matrix(DEGpoisd$dd)
  rownames(DEGsamplePoisDistMatrix) = dds$cell_line
  colnames(DEGsamplePoisDistMatrix) = NULL
  
  # Same for genes
  DEGpoisd_genes = PoissonDistance(counts(dds[topDEGs$Gene.Symbol,]))
  DEGgenesPoisDistMatrix = as.matrix(DEGpoisd_genes$dd)
  rownames(DEGgenesPoisDistMatrix) = topDEGs$Gene.Symbol
  colnames(DEGgenesPoisDistMatrix) = NULL
  
  # Log-transformed matrix
  logmat = log2(matcount + 1)
  
  # The same matrix but order samples by metastatic status and genes by logFC (raw)
  genes_order = topDEGs %>% dplyr::filter(padj < 0.05) %>%
    arrange(log2FoldChange)
  gaps_row = length(which(genes_order$log2FoldChange < 0))
  genes_order = genes_order$Gene.Symbol
  
  sample_order = anno %>% dplyr::arrange(treatment)
  gaps_col = length(which(sample_order$treatment == "control"))
  sample_order = rownames(sample_order)
  
  matcount_ordered = matcount[genes_order, sample_order]
  
  # Mean-centering
  logmat_znorm = (logmat - rowMeans(logmat))/sqrt(rowVars(logmat))
  logmat_ordered_znorm = logmat_znorm[genes_order, sample_order]
  
  # In the plot below the cluster distances are not affected by the log-transformation
  # of the matrix. It is only used to make a smoother color scale
  tiff(paste0("DGEA/Post_DGEA_heatmaps/raw_lfcs/",
              names(results)[i], "_DEG_Poisson_clustered_samples_and_genes.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
  pheatmap(logmat_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           clustering_distance_rows = DEGpoisd_genes$dd,
           clustering_distance_cols = DEGpoisd$dd,
           col = colors2,
           cutree_cols = 2,
           cutree_rows = 2,
           fontsize = 8,
           fontsize_row = 8,
           fontsize_col = 6,
           main = paste0(names(results)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  
  # Ordered genes and samples
  tiff(paste0("DGEA/Post_DGEA_heatmaps/raw_lfcs/",
              names(results)[i], "_DEG_samples_ordered.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
  pheatmap(logmat_ordered_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           col = colors2,
           gaps_row = gaps_row,
           gaps_col = gaps_col,
           cutree_cols = 2,
           fontsize = 8,
           fontsize_row = 8,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(results)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
}
rm(genes, i, res, topDEGs, matcount, matcount_ordered, DEGgenesPoisDistMatrix,
   DEGpoisd, DEGpoisd_genes, DEGsamplePoisDistMatrix, logmat, logmat_ordered_znorm,
   logmat_znorm, sample_order, genes_order, gaps_col)

# Heatmaps for ashr-shrunk LFCs
for (i in 1:length(ashr_shrunk_results)) {
  # Construct the data frame
  genes = rownames(ashr_shrunk_results[[i]])
  res = ashr_shrunk_results[[i]]@listData %>%
    as.data.frame() %>% 
    mutate(Gene.Symbol = genes) %>%
    dplyr::arrange(padj)
  
  topDEGs = res[which(res$padj < 0.05),]
  
  # Gene clustering with raw counts and Poisson distance
  topDEGs = topDEGs[1:n_degs, ]
  matcount = assay(dds)[topDEGs$Gene.Symbol,]
  DEGpoisd = PoissonDistance(t(counts(dds[topDEGs$Gene.Symbol,])))
  DEGsamplePoisDistMatrix = as.matrix(DEGpoisd$dd)
  rownames(DEGsamplePoisDistMatrix) = dds$cell_line
  colnames(DEGsamplePoisDistMatrix) = NULL
  
  # Same for genes
  DEGpoisd_genes = PoissonDistance(counts(dds[topDEGs$Gene.Symbol,]))
  DEGgenesPoisDistMatrix = as.matrix(DEGpoisd_genes$dd)
  rownames(DEGgenesPoisDistMatrix) = topDEGs$Gene.Symbol
  colnames(DEGgenesPoisDistMatrix) = NULL
  
  # Log-transformed matrix
  logmat = log2(matcount + 1)
  
  # The same matrix but order samples by metastatic status and genes by logFC (raw)
  genes_order = topDEGs %>% dplyr::filter(padj < 0.05) %>%
    arrange(log2FoldChange)
  gaps_row = length(which(genes_order$log2FoldChange < 0))
  genes_order = genes_order$Gene.Symbol
  
  sample_order = anno %>% dplyr::arrange(treatment)
  gaps_col = length(which(sample_order$treatment == "control"))
  sample_order = rownames(sample_order)
  
  matcount_ordered = matcount[genes_order, sample_order]
  
  # Mean-centering
  logmat_znorm = (logmat - rowMeans(logmat))/sqrt(rowVars(logmat))
  logmat_ordered_znorm = logmat_znorm[genes_order, sample_order]
  
  # In the plot below the cluster distances are not affected by the log-transformation
  # of the matrix. It is only used to make a smoother color scale
  tiff(paste0("DGEA/Post_DGEA_heatmaps/ashr_shrunk/",
              names(ashr_shrunk_results)[i], "_DEG_Poisson_clustered_samples_and_genes.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
  pheatmap(logmat_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           clustering_distance_rows = DEGpoisd_genes$dd,
           clustering_distance_cols = DEGpoisd$dd,
           col = colors2,
           cutree_cols = 2,
           cutree_rows = 2,
           fontsize = 8,
           fontsize_row = 8,
           fontsize_col = 6,
           main = paste0(names(ashr_shrunk_results)[i], 
                         " - DEGs vs. samples: Poisson heatmap, cluster cut = 2"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
  
  
  # Ordered genes and samples
  tiff(paste0("DGEA/Post_DGEA_heatmaps/ashr_shrunk/",
              names(ashr_shrunk_results)[i], "_DEG_samples_ordered.tiff"),
       res = 700, width = 7260, 
       height = 5760, compression = "lzw")
  pheatmap(logmat_ordered_znorm, annotation_col = anno,
           breaks = seq(-4, 4, 8/length(colors2)),
           annotation_colors = ann_colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           col = colors2,
           gaps_row = gaps_row,
           gaps_col = gaps_col,
           cutree_cols = 2,
           fontsize = 8,
           fontsize_row = 8,
           fontsize_col = 6,
           border_color = "grey60",
           main = paste0(names(ashr_shrunk_results)[i], 
                         " - DEGs vs. samples heatmap (ordered)"),
           legend_breaks = c(-4, 4),
           legend_labels = c("lower expression", "higher expression"))
  dev.off()
}

rm(genes, i, res, topDEGs, matcount, matcount_ordered, DEGgenesPoisDistMatrix,
   DEGpoisd, DEGpoisd_genes, DEGsamplePoisDistMatrix, logmat, logmat_ordered_znorm,
   logmat_znorm, sample_order, genes_order, gaps_col)

gc()

# Session info
sessionInfo()
