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
library(mclust)
library(aricode)
library(EnhancedVolcano)
library(genefilter)
library(grid)

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
# 57925   3585

print(paste0(round(100*3585/(3585+57925), 2), "% of genes only have zero counts."))
# 5.83% of genes only have zero counts

# We remove these genes from the experiment
toremove = which(rowSums(dds_pre@assays@data@listData[["counts"]]==0)==80)
dds_pre = dds_pre[-toremove, ]

rm(toremove); gc()

# Generate GPCA and MDS plots to check how similar samples of the same patient
# are to each other

# glmpca
RNGversion("4.2.2")
set.seed(123)
gpca = glmpca(counts(dds_pre), L=2, fam = "nb", minibatch = "none",
              optimizer = "fisher",
              ctl = list(maxIter = 10000, tol = 1e-6))
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

mdsPois <- as.data.frame(colData(dds_pre)) %>%
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

png("Results/Exploratory analysis/dimred_prior_to_sum.png", width = 2980, height = 1380, res = 700)
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

coldata = coldata %>% 
  dplyr::select(-Recurrence, -names) %>%
  inner_join(Pheno, by = "Patient.ID")

# Remove genes that were not mapped to gene symbols and still have ENSG# ID's
# which map to RF####, AL####.# sequences
keep.symbols = setdiff(rownames(dds_pre), 
                       rownames(dds_pre)[grep("ENSG00", rownames(dds_pre))])
dds_pre = dds_pre[keep.symbols, ]
summed_counts = summed_counts[keep.symbols, ]

# Create a full design matrix for the filtering
design = model.matrix(~ coldata$Sample.type*coldata$Stage*coldata$Recurrence)
rownames(design) = rownames(coldata)

# # We need to filter out uninformative genes
# keep.exprs = names(which(filterByExpr(summed_counts, design = design) == TRUE))
# 
# # This design results in ~4,000 genes less (38339 -> 34049)
# summed_counts = summed_counts[keep.exprs, ]

# Converting to factors
coldata$Stage = factor(coldata$Stage, levels = c("Early", "Late"),
                       labels = c("Early", "Late"))
coldata$Sample.type = factor(coldata$Sample.type, levels = c("Normal", "Tumor"),
                             labels = c("Normal", "Tumor"))
coldata$Recurrence = factor(coldata$Recurrence, levels = c("No recurrence", "Recurrence"),
                            labels = c("No recurrence", "Recurrence"))

# Create a new DESeqDataSet object from the summed counts matrix
dds <- DESeqDataSetFromMatrix(countData = summed_counts,
                              colData = coldata,
                              design = ~Sample.type*Stage*Recurrence)

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
ggsave(filename = "TSCC_filtered_read_counts_density_curves.png",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
dev.off()

# The previous plot contains count outliers (abnormally high counts) which make 
# it difficult to see the distribution of counts near zero. We replot with 
# xlim 30

create_density_curve(counts(dds)) +
  scale_x_continuous(limits = c(-1, 30), breaks = seq(0, 30, 10)) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "TSCC_filtered_read_counts_density_curves_xlim30.png",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
dev.off()

# Color density plots for counts
create_density_plot_color(counts(dds)) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "TSCC_filtered_read_counts_color_densities.png",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
dev.off()

# With xlim 30
create_density_plot_color(counts(dds))+
  scale_x_continuous(limits = c(-1, 30), breaks = seq(0, 30, 10)) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "TSCC_filtered_read_counts_color_densities_xlim30.png",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
dev.off()

# Now with normalized counts

# Firstly, transform the normalized counts with a data-driven pseudocount
# https://www.biorxiv.org/content/10.1101/404962v1.full -> Process description

# We need a grouping for the TSCC samples. We are going to go for 
# Recurrence vs No recurrence
rec_samples = Pheno$Patient.ID[Pheno$Recurrence == "Recurrence"]
norec_samples = Pheno$Patient.ID[Pheno$Recurrence == "No recurrence"]

# Get size factors
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
ggsave(filename = "TSCC_filtered_log2normcounts_density_curves.png",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
dev.off()

# Color density plots for counts
create_density_plot_color(log2.norm.counts.plus1) +
  theme(plot.title = element_text(face = "bold"))
ggsave(filename = "TSCC_filtered_log2normcounts_color_densities.png",
       path = "Results/Exploratory analysis",
       width = 1920, height = 1080, device = 'png', units = "px",
       dpi = 150)
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
png("Results/Exploratory analysis/meanSdPlot.png", res = 300, width = 4500, height = 1600)
fig
annotate_figure(fig, left = textGrob(bquote(bold("sd")), rot = 90, vjust = 1, gp = gpar(cex = 1)),
                bottom = textGrob(bquote(bold("rank (mean)")), gp = gpar(cex = 1)))
dev.off()

# Heatmap of the count matrix - top 50 genes
select1 = order(rowMeans(log2.norm.counts.plus1),
               decreasing=TRUE)[1:50]
select2 = order(rowMeans(assay(rld)),
                decreasing=TRUE)[1:50]
select3 = order(rowMeans(assay(vsd)),
                decreasing=TRUE)[1:50]
df = as.data.frame(colData(dds)[,c("Recurrence", "Sample.type", "Stage")]) %>%
  dplyr::arrange(Recurrence, Sample.type, Stage)

# Annotation colors
ann_colors = list(
  Recurrence = c(`No recurrence` = "dodgerblue4", Recurrence = "deeppink4"),
  Sample.type = c(Tumor = "purple4", Normal = "orange"),
  Stage = c(Early = carto_pal(n = 7, "Emrld")[4], 
                Late = carto_pal(n = 7, "Burg")[4])
)

# log2(norm.counts + 1)
png("Results/Exploratory analysis/top50_heatmap_log2(norm.counts+1).png",
     res = 200, width = 1420, 
     height = 1080)
pheatmap(log2.norm.counts.plus1[select1, rownames(df)], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "log2(norm.counts + 1)",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

# rlog
png("Results/Exploratory analysis/top50_heatmap_rlog.png",
     res = 200, width = 1420, 
     height = 1080)
pheatmap(assay(rld)[select2, rownames(df)], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "rlog",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

# vst
png("Results/Exploratory analysis/top50_heatmap_vst.png", res = 200, width = 1420, 
     height = 1080)
pheatmap(assay(vsd)[select3, rownames(df)], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df, main = "vst",
         annotation_colors = ann_colors, fontsize = 6)
dev.off()

# We will only use the vst-transformed data for the next visualisations

# Heatmap of sample distances
sampleVSTdists = dist(t(assay(vsd)))
sampleDistMatrix = as.matrix(sampleVSTdists)
#colnames(sampleDistMatrix) = NULL
colors = colorRampPalette(viridisLite::magma(10))(255)

png("Results/Exploratory analysis/sample_distance_heatmap_vst.png", 
     res = 700, width = 6760, 
     height = 4760)
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

png("Results/Exploratory analysis/sample_poisson_distance_heatmap.png", res = 700, width = 6760, 
     height = 4760)
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         fontsize = 8,
         main = "Sample Poisson distance heatmap",
         annotation_col = df,
         annotation_colors = ann_colors)
dev.off()

# Check the ARI of the two main clusters
vst_clust2 = c("COHTC-37", "COHTC-16", "COHTC-10", "COHTN-12",
               "COHTC-68", "COHTC-61", "COHTC-8")
vst_clust1 = setdiff(colnames(sampleDistMatrix), vst_clust2)
vst_clust = c(rep("clust1", length(vst_clust1)), rep("clust2", length(vst_clust2)))
names(vst_clust) = c(vst_clust1, vst_clust2)

pois_clust2 = c("COHTC-37", "COHTC-16", "COHTC-10", "COHTN-12",
                "COHTN-13", "COHTC-61", "COHTC-68", "COHTC-8")
pois_clust1 = setdiff(colnames(sampleDistMatrix), pois_clust2)
pois_clust = c(rep("clust1", length(pois_clust1)), rep("clust2", length(pois_clust2)))
names(pois_clust) = c(pois_clust1, pois_clust2)

pois_clust = pois_clust[names(vst_clust)]

# ARI
ari = adjustedRandIndex(as.factor(vst_clust), as.factor(pois_clust))
cat("Adjusted Rand Index (ARI):", ari, "\n")

# Dimensionality reduction #####
# Principal Component Analysis
pcaData = plotPCA(vsd, intgroup = c("Recurrence", "Sample.type"), 
                  returnData=TRUE)
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
gpca = glmpca(counts(dds), L=2, fam = "nb", minibatch = "none",
              optimizer = "fisher",
              ctl = list(maxIter = 10000, tol = 1e-6))
gpca.dat = gpca$factors
gpca.dat$Sample.type = dds$Sample.type
gpca.dat$Stage = dds$Stage
gpca.dat$Recurrence = dds$Recurrence

gpcaplot = ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Sample.type, shape = Recurrence)) +
  geom_point(size = 0.2, alpha = 0.65) +
  geom_text_repel(aes(label = Stage), segment.linetype = 1,
                  segment.color = "black", segment.size = 0.1,
                  min.segment.length = 0, max.overlaps = Inf,
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

png("Results/Exploratory analysis/PCA.png", width = 2880, height = 1380, res = 700)
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
                  size = 0.5, nudge_y = 0.03, nudge_x = 0.03, max.overlaps = 50) +
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
                  size = 0.5, nudge_y = 0.03, nudge_x = 0.03, max.overlaps = 50) +
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
png("Results/Exploratory analysis/MDS.png", res = 700, width = 2880, height = 1080)
fig
annotate_figure(fig, left = textGrob(bquote(bold("MDS2")), rot = 90, vjust = 1, gp = gpar(cex = 0.45)))
dev.off()

# Export summed counts
write.xlsx(as.data.frame(counts(dds)), "Resources_and_environments/summed_counts.xlsx",
           rowNames = TRUE)

# Export a session info text file for reproducibility
writeLines(capture.output(sessionInfo()), 
           paste0("Resources_and_environments/Preprocessing_and_EA_sessionInfo.txt"))

# Save environment for easy loading for DGEA later
save.image("Resources_and_environments/Preprocessing_and_exploratory_analysis.RData")