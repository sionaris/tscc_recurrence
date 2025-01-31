# This script uses the lists produced by DGEA for the purpose of ORA
# (Over-representation analysis) and GSEA
# (Gene Set Enrichment Analysis)
# and subsequently Active Subnetwork Analysis with pathfindR

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(openxlsx)
library(ReactomePA)
library(org.Hs.eg.db)
library(pathfindR)
library(cowplot)
library(DESeq2)
library(fastcluster)
library(parallel)
library(foreach)
library(doParallel)

# Import data (4 lists from DGEA)
ashr_dgea = list()
for (i in 1:4) {
  ashr_dgea[[i]] = openxlsx::read.xlsx("Results/DGEA/DGEA_ashr.xlsx", sheet = i)
}
rm(i); gc()
 
names(ashr_dgea) = paste0("ashr_", 
                          c("TR_vs_TNR", "LTR_vs_ETR", 
                            "ETR_vs_ETNR", "LTR_vs_LTNR"))

# Import dds object
dds = readRDS("Resources_and_environments/final_dds.rds")

# Create output folders
for (i in 1:length(ashr_dgea)) {
  if (!dir.exists(paste0("Results/Pathways/", names(ashr_dgea)[i]))) {
    dir.create(paste0("Results/Pathways/", names(ashr_dgea)[i]))
  }
}

pathway_approaches = c("ORA", "GSEA", "pathfindR")

for (i in 1:length(pathway_approaches)) {
  for (j in 1:length(ashr_dgea)) {
    if (!dir.exists(paste0("Results/Pathways/", names(ashr_dgea)[j], "/",
                           pathway_approaches[[i]]))) {
      dir.create(paste0("Results/Pathways/", names(ashr_dgea)[j], "/",
                        pathway_approaches[[i]]))
    }
  }
}

rm(i, j); gc()

# Create lists of workbooks
ora_wb_list = list()
gsea_wb_list = list()
pathfindR_wb_list = list()

for (i in 1:length(ashr_dgea)) {
  ora_wb_list[[i]] = createWorkbook()
  gsea_wb_list[[i]] = createWorkbook()
  pathfindR_wb_list[[i]] = createWorkbook()
}
names(ora_wb_list) = names(gsea_wb_list) = names(pathfindR_wb_list) =
  names(ashr_dgea)
rm(i); gc()

# Create lists of output
ora_results = list()
gsea_results = list()
pathfindR_results = list()

# Set download option to libcurl for clusterProfiler
R.utils::setOption("clusterProfiler.download.method","libcurl")

# Objects to keep after every loop
keepers = c("ashr_dgea", "gsea_wb_list", "ora_wb_list", "pathfindR_wb_list",
            "gsea_results", "ora_results", "pathfindR_results", "pathway_approaches",
            "i", "keepers")

# ORA and GSEA #####
for (i in 1:length(ashr_dgea)) {
  
  # Define the cutoff for significant genes
  sig_genes = ashr_dgea[[i]]$padj < 0.05
  
  # Extract the log2 fold changes and gene names of the significant genes
  log2fc = ashr_dgea[[i]]$log2FoldChange[sig_genes]
  genes = ashr_dgea[[i]]$Gene.Symbol[sig_genes]
  gene_list = log2fc
  names(gene_list) = genes
  rm(log2fc, genes)
  gene_list = sort(gene_list, decreasing = TRUE) # symbols
  gene_list_conv = suppressWarnings(bitr(names(gene_list), fromType = "SYMBOL",
                                         toType = "ENTREZID", OrgDb = org.Hs.eg.db))
  gene_list_entrez = gene_list
  names(gene_list_entrez) = as.character(gene_list_conv$ENTREZID)
  gene_list_entrez = sort(gene_list_entrez, decreasing = TRUE) # Entrez ids
  entrez_universe = suppressWarnings(bitr(ashr_dgea[[i]]$Gene.Symbol, fromType = "SYMBOL",
                                          toType = "ENTREZID", OrgDb = org.Hs.eg.db))$ENTREZID
  
  # Gene Set Enrichment Analysis (GSEA) #####
  
  # GSEA temporary list
  gsea_output = list()
  
  # Perform GO gene set enrichment analysis using clusterProfiler
  RNGversion("4.2.2")
  set.seed(123)
  gseaGO = suppressWarnings(gseGO(geneList = gene_list,
                                  ont = "ALL",
                                  OrgDb = "org.Hs.eg.db",
                                  keyType = "SYMBOL",
                                  pvalueCutoff = 0.1,
                                  pAdjustMethod = "BH",
                                  minGSSize = 3,
                                  maxGSSize = 800,
                                  seed = TRUE))
  
  if (exists("gseaGO")) {
    if (nrow(gseaGO) > 0) {
      # gseaGO = clusterProfiler::setReadable(gseaGO, 'org.Hs.eg.db')
      gsea_output[["Gene Ontology"]] = gseaGO@result
      
      addWorksheet(gsea_wb_list[[i]], "GSEA_GO")
      writeData(gsea_wb_list[[i]], "GSEA_GO", as.data.frame(gseaGO@result))
    }
  }
  
  # KEGG
  RNGversion("4.2.2")
  set.seed(123)
  gseaKEGG = suppressWarnings(gseKEGG(geneList = gene_list_entrez,
                                      organism = "hsa",
                                      keyType = "ncbi-geneid",
                                      pvalueCutoff = 0.1,
                                      pAdjustMethod = "BH",
                                      minGSSize = 3,
                                      maxGSSize = 800,
                                      seed = TRUE))
  
  if (exists("gseaKEGG")) {
    if (nrow(gseaKEGG) > 0) {
      gseaKEGG = clusterProfiler::setReadable(gseaKEGG, 'org.Hs.eg.db', keyType = "ENTREZID")
      gsea_output[["KEGG"]] = gseaKEGG@result
      
      addWorksheet(gsea_wb_list[[i]], "GSEA_KEGG")
      writeData(gsea_wb_list[[i]], "GSEA_KEGG", as.data.frame(gseaKEGG@result))
    }
  }
  
  # Reactome
  RNGversion("4.2.2")
  set.seed(123)
  gseaReactome = suppressWarnings(gsePathway(geneList = gene_list_entrez,
                                             organism = "human",
                                             pvalueCutoff = 0.1,
                                             pAdjustMethod = "BH",
                                             minGSSize = 3,
                                             maxGSSize = 800,
                                             seed = TRUE))
  
  if (exists("gseaReactome")) {
    if (nrow(gseaReactome) > 0) {
      gseaReactome = clusterProfiler::setReadable(gseaReactome, 'org.Hs.eg.db',
                                                  keyType = "ENTREZID")
      gsea_output[["Reactome"]] = gseaReactome@result
      
      addWorksheet(gsea_wb_list[[i]], "GSEA_Reactome")
      writeData(gsea_wb_list[[i]], "GSEA_Reactome", as.data.frame(gseaReactome@result))
    }
  }
  
  # WikiPathways
  RNGversion("4.2.2")
  set.seed(123)
  gseaWP = suppressWarnings(gseWP(geneList = gene_list_entrez,
                                  organism = "Homo sapiens",
                                  pvalueCutoff = 0.1,
                                  pAdjustMethod = "BH",
                                  minGSSize = 3,
                                  maxGSSize = 800,
                                  seed = TRUE))
  
  if (exists("gseaWP")) {
    if (nrow(gseaWP) > 0) {
      gseaWP = clusterProfiler::setReadable(gseaWP, 'org.Hs.eg.db',
                                            keyType = "ENTREZID")
      gsea_output[["WikiPathways"]] = gseaWP@result
      
      addWorksheet(gsea_wb_list[[i]], "GSEA_WP")
      writeData(gsea_wb_list[[i]], "GSEA_WP", as.data.frame(gseaWP@result))
    }
  }
  
  # DO analysis (Disease Enrichment)
  RNGversion("4.2.2")
  set.seed(123)
  gseaDO = suppressWarnings(DOSE::gseDO(geneList = gene_list_entrez,
                                        pvalueCutoff = 0.1,
                                        pAdjustMethod = "BH",
                                        minGSSize = 3,
                                        maxGSSize = 800,
                                        seed = TRUE))
  
  if (exists("gseaDO")) {
    if (nrow(gseaDO) > 0) {
      gseaDO = clusterProfiler::setReadable(gseaDO, 'org.Hs.eg.db',
                                            keyType = "ENTREZID")
      gsea_output[["Disease Ontology"]] = gseaDO@result
      
      addWorksheet(gsea_wb_list[[i]], "GSEA_DO")
      writeData(gsea_wb_list[[i]], "GSEA_DO", as.data.frame(gseaDO@result))
    }
  }
  
  # DO in the Network Cancer Gene
  # http://ncg.kcl.ac.uk/
  RNGversion("4.2.2")
  set.seed(123)
  gseaNCG = suppressWarnings(DOSE::gseNCG(geneList = gene_list_entrez,
                        pvalueCutoff = 0.1,
                        pAdjustMethod = "BH",
                        minGSSize = 3,
                        maxGSSize = 800,
                        seed = TRUE))

  if (exists("gseaNCG")) {
   if (nrow(gseaNCG) > 0) {
     gseaNCG = clusterProfiler::setReadable(gseaNCG, 'org.Hs.eg.db',
                                            keyType = "ENTREZID")
     gsea_output[["Network Cancer Gene (NCG)"]] = gseaNCG@result

     addWorksheet(gsea_wb_list[[i]], "GSEA_NCG")
     writeData(gsea_wb_list[[i]], "GSEA_NCG", as.data.frame(gseaNCG@result))
   }
  }
  
  # DO in the Disease Gene Network (DisGeNET)
  # http://disgenet.org/
  RNGversion("4.2.2")
  set.seed(123)
  gseaDGN = suppressWarnings(DOSE::gseDGN(geneList = gene_list_entrez,
                       pvalueCutoff = 0.1,
                       pAdjustMethod = "BH",
                       minGSSize = 3,
                       maxGSSize = 800,
                       seed = TRUE))
  if (exists("gseaDGN")) {
  if (nrow(gseaDGN) > 0) {
   gseaDGN = clusterProfiler::setReadable(gseaDGN, 'org.Hs.eg.db',
                                          keyType = "ENTREZID")
   gsea_output[["Disease Gene Network (DisGeNET)"]] = gseaDGN@result

    addWorksheet(gsea_wb_list[[i]], "GSEA_DGN")
    writeData(gsea_wb_list[[i]], "GSEA_DGN", as.data.frame(gseaDGN@result))
  }
  }
  
  gsea_results[[i]] = gsea_output
  
  # Over-representation analysis (ORA) #####
  
  # ORA temporary list
  ora_output = list()
  
  # GO
  oraGO = suppressWarnings(enrichGO(gene = names(gene_list),
                                    ont = "ALL",
                                    universe = ashr_dgea[[i]]$Gene.Symbol,
                                    OrgDb = "org.Hs.eg.db",
                                    keyType = "SYMBOL",
                                    pvalueCutoff = 0.1,
                                    qvalueCutoff = 0.2,
                                    pAdjustMethod = "BH",
                                    minGSSize = 3,
                                    maxGSSize = 800,
                                    pool = TRUE))
  
  if (exists("oraGO")) {
    if (nrow(oraGO) > 0) {
      oraGO_q0.2 = oraGO@result %>% dplyr::filter(qvalue <= 0.2)
      ora_output[["Gene Ontology"]] = oraGO@result
      
      addWorksheet(ora_wb_list[[i]], "ORA_GO")
      writeData(ora_wb_list[[i]], "ORA_GO", as.data.frame(oraGO_q0.2))
    }
  }
  
  # KEGG
  oraKEGG = suppressWarnings(enrichKEGG(gene = names(gene_list_entrez),
                                        universe = entrez_universe,
                                        organism = "hsa",
                                        keyType = "ncbi-geneid",
                                        pvalueCutoff = 0.1,
                                        qvalueCutoff = 0.2,
                                        pAdjustMethod = "BH",
                                        minGSSize = 3,
                                        maxGSSize = 800))
  
  if (exists("oraKEGG")) {
    if (nrow(oraKEGG) > 0) {
      oraKEGG = clusterProfiler::setReadable(oraKEGG, 'org.Hs.eg.db', keyType = "ENTREZID")
      oraKEGG_q0.2 = oraKEGG@result %>% dplyr::filter(qvalue <= 0.2)
      ora_output[["KEGG"]] = oraKEGG@result
      
      addWorksheet(ora_wb_list[[i]], "ORA_KEGG")
      writeData(ora_wb_list[[i]], "ORA_KEGG", as.data.frame(oraKEGG_q0.2))
    }
  }
  
  # WikiPathways
  oraWP = suppressWarnings(enrichWP(gene = names(gene_list_entrez),
                                    universe = entrez_universe,
                                    organism = "Homo sapiens",
                                    pvalueCutoff = 0.1,
                                    qvalueCutoff = 0.2,
                                    pAdjustMethod = "BH",
                                    minGSSize = 3,
                                    maxGSSize = 800))
  
  if (exists("oraWP")) {
    if (nrow(oraWP) > 0) {
      oraWP = clusterProfiler::setReadable(oraWP, 'org.Hs.eg.db', keyType = "ENTREZID")
      oraWP_q0.2 = oraWP@result %>% dplyr::filter(qvalue <= 0.2)
      ora_output[["WikiPathways"]] = oraWP@result
      
      addWorksheet(ora_wb_list[[i]], "ORA_WP")
      writeData(ora_wb_list[[i]], "ORA_WP", as.data.frame(oraWP_q0.2))
    }
  }
  
  # Reactome
  oraReactome = suppressWarnings(enrichPathway(gene = names(gene_list_entrez),
                                               universe = entrez_universe,
                                               organism = "human",
                                               pvalueCutoff = 0.1,
                                               qvalueCutoff = 0.2,
                                               pAdjustMethod = "BH",
                                               minGSSize = 3,
                                               maxGSSize = 800,
                                               readable = TRUE))
  
  if (exists("oraReactome")) {
    if (nrow(oraReactome) > 0) {
      oraReactome = clusterProfiler::setReadable(oraReactome, 'org.Hs.eg.db',
                                   keyType = "ENTREZID")
      ora_output[["Reactome"]] = oraReactome@result
      oraReactome_q0.2 = oraReactome@result %>% dplyr::filter(qvalue <= 0.2)
      
      addWorksheet(ora_wb_list[[i]], "ORA_Reactome")
      writeData(ora_wb_list[[i]], "ORA_Reactome", as.data.frame(oraReactome_q0.2))
    }
  }
  
  # DO
  oraDO = suppressWarnings(DOSE::enrichDO(gene = names(gene_list_entrez),
                        universe = entrez_universe,
                        pvalueCutoff = 0.1,
                        qvalueCutoff = 0.2,
                        pAdjustMethod = "BH",
                        minGSSize = 3,
                        maxGSSize = 800,
                        readable = TRUE))

  if (exists("oraDO")) {
   if (nrow(oraDO) > 0) {
     oraDO = clusterProfiler::setReadable(oraDO, 'org.Hs.eg.db',
                                          keyType = "ENTREZID")
     ora_output[["Disease Ontology"]] = oraDO@result
     oraDO_q0.2 = oraDO@result %>% dplyr::filter(qvalue <= 0.2)

     addWorksheet(ora_wb_list[[i]], "ORA_DO")
     writeData(ora_wb_list[[i]], "ORA_DO", as.data.frame(oraDO_q0.2))
   }
  }

  # DO: NCG
  oraNCG = suppressWarnings(DOSE::enrichNCG(gene = names(gene_list_entrez),
                          universe = entrez_universe,
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.2,
                          pAdjustMethod = "BH",
                          minGSSize = 3,
                          maxGSSize = 800,
                          readable = TRUE))

  if (exists("oraNCG")) {
     if (nrow(oraNCG) > 0) {
       oraNCG = clusterProfiler::setReadable(oraNCG, 'org.Hs.eg.db',
                                            keyType = "ENTREZID")
       ora_output[["Network Cancer Gene (NCG)"]] = oraNCG@result
       oraNCG_q0.2 = oraNCG@result %>% dplyr::filter(qvalue <= 0.2)

       addWorksheet(ora_wb_list[[i]], "ORA_NCG")
       writeData(ora_wb_list[[i]], "ORA_NCG", as.data.frame(oraNCG_q0.2))
     }
  }
  
  # DO: DGN
  oraDGN = suppressWarnings(DOSE::enrichDGN(gene = names(gene_list_entrez),
                          universe = entrez_universe,
                          pvalueCutoff = 0.1,
                          qvalueCutoff = 0.2,
                          pAdjustMethod = "BH",
                          minGSSize = 3,
                          maxGSSize = 800,
                          readable = TRUE))

  if (exists("oraDGN")) {
   if (nrow(oraDGN) > 0) {
     oraDGN = clusterProfiler::setReadable(oraDGN, 'org.Hs.eg.db',
                                           keyType = "ENTREZID")
     ora_output[["Disease Gene Network (DisGeNET)"]] = oraDGN@result
     oraDGN_q0.2 = oraDGN@result %>% dplyr::filter(qvalue <= 0.2)

     addWorksheet(ora_wb_list[[i]], "ORA_DGN")
     writeData(ora_wb_list[[i]], "ORA_DGN", as.data.frame(oraDGN_q0.2))
   }
  }
  
  ora_results[[i]] = ora_output
  
  # Export the workbooks
  saveWorkbook(ora_wb_list[[i]], paste0("Results/Pathways/", names(ashr_dgea)[i],
                                        "/ORA/ORA_output_", names(ashr_dgea)[i],
                                        ".xlsx"),
               overwrite = TRUE)
  saveWorkbook(gsea_wb_list[[i]], paste0("Results/Pathways/", names(ashr_dgea)[i],
                                         "/GSEA/GSEA_output_", names(ashr_dgea)[i],
                                         ".xlsx"),
               overwrite = TRUE)
  
  # Keep track of the loop's progress
  cat(paste("Done with", names(ashr_dgea)[i], "\n"))
  
  # Remove garbage
  rm(list=setdiff(ls(), keepers))
}
names(ora_results) = names(gsea_results) = names(ashr_dgea)

rm(i); gc()

# Pathway analysis with pathfindR #####
niter = 50
n_processes = 10
min_gset_size = 10
max_gset_size = 300
# score_quan_thr = 0.8
# sig_gene_thr = 0.02
enrichment_threshold = 0.05
p_val_threshold = 0.05
search_method = "GR"
pin_name_path = "Biogrid"
adj_method = 'fdr'
# combined_res = NULL

# ashr_TR_vs_TNR #####
# Preparing a pathfindR loop for enrichment analysis
subdirs = c("Reactome", "BioCarta", "GO-BP", "GO-CC", "GO-MF", "KEGG")
pathf_input_ashr_TR_vs_TNR = ashr_dgea$ashr_TR_vs_TNR %>%
  dplyr::select(Gene.Symbol, log2FoldChange, padj) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_ashr_TR_vs_TNR = paste0("Results/Pathways/ashr_TR_vs_TNR/pathfindR/", subdirs)
pathfindR_outputs_ashr_TR_vs_TNR = list()

RNGversion("4.2.2")
set.seed(123)
for (i in 1:length(dirs_ashr_TR_vs_TNR)){
  pathfindR_outputs_ashr_TR_vs_TNR[[i]] = run_pathfindR(pathf_input_ashr_TR_vs_TNR, gene_sets = subdirs[i],
                                                         p_val_threshold = p_val_threshold, convert2alias = FALSE,
                                                         output_dir = dirs_ashr_TR_vs_TNR[i], min_gset_size = min_gset_size,
                                                         max_gset_size = max_gset_size, adj_method = adj_method,
                                                         enrichment_threshold = enrichment_threshold,
                                                         pin_name_path = pin_name_path, search_method = search_method,
                                                         iterations = niter, n_processes = n_processes)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_ashr_TR_vs_TNR) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

# Use a custom fast approach
# Set up the number of cores to use
nonzero_res = which(lapply(pathfindR_outputs_ashr_TR_vs_TNR, nrow) != 0)
cl <- makeCluster(length(nonzero_res))
registerDoParallel(cl)

timestamp()
clustered_results_ashr_TR_vs_TNR <- foreach(i = nonzero_res, .packages = c("pathfindR", "fastcluster")) %dopar% {
  RNGversion("4.2.2")
  set.seed(123)
  source("R_Scripts/fast_pathfindR_hclust.R")
  
  # Perform clustering, handle errors
  result <- cluster_enriched_terms_fast(pathfindR_outputs_ashr_TR_vs_TNR[[i]], method = "hierarchical", plot_clusters_graph = FALSE,
                                        use_description = TRUE, use_active_snw_genes = FALSE)
  if (is.character(result) && result == "hclust impossible") {
    return("hclust impossible")
  } else {
    return(result)
  }
}

timestamp() # <0.5 mins
stopCluster(cl)
gc()

# Unlist
for (i in 1:length(clustered_results_ashr_TR_vs_TNR)) {
  clustered_results_ashr_TR_vs_TNR[[i]] = clustered_results_ashr_TR_vs_TNR[[i]]$clustered_df
}
names(clustered_results_ashr_TR_vs_TNR) = names(pathfindR_outputs_ashr_TR_vs_TNR)[nonzero_res]

# Are there any null sets?
which(clustered_results_ashr_TR_vs_TNR == "hclust impossible")

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_ashr_TR_vs_TNR = pathfindR_outputs_ashr_TR_vs_TNR
for (i in nonzero_res){
  wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_ashr_TR_vs_TNR[[i]]$Term_Description, 
                                                                                     width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR = clustered_results_ashr_TR_vs_TNR
for (i in 1:length(wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR)){
  wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR[[i]]$Term_Description = stringr::str_wrap(clustered_results_ashr_TR_vs_TNR[[i]]$Term_Description, 
                                                                                               width = 41)
}
rm(i)

comparisons = c("TR vs TNR", "LTR vs ETR", 
                "ETR vs ETNR", "LTR vs LTNR")
names(comparisons) = c("ashr_TR_vs_TNR", "ashr_LTR_vs_ETR",
                       "ashr_ETR_vs_ETNR", "ashr_LTR_vs_LTNR")

enrichment_dotplots_ashr_TR_vs_TNR = list()
cluster_enrichment_dotplots_ashr_TR_vs_TNR = list()

# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_TR_vs_TNR)[i]
  # unclustered results
  enrichment_dotplots_ashr_TR_vs_TNR[[names]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[names]],
                                                             top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    scale_size(range = c(0.1, 2)) +
    theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 3.5),
          axis.title.x = element_text(size = 4.5, face = "bold"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 3.5, face = "bold"),
          legend.text = element_text(size = 2.75))+
    labs(title = paste0("Top 10 ", names,
                        " terms enrichment dotplot - (", comparisons["ashr_TR_vs_TNR"], ")"))
  print(enrichment_dotplots_ashr_TR_vs_TNR[[names]])
  ggsave(filename = paste0(names, "_top10_dotplot.tiff"),
         path = paste0("Results/Pathways/ashr_TR_vs_TNR/pathfindR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  if (names %in% names(clustered_results_ashr_TR_vs_TNR)) {
    cluster_enrichment_dotplots_ashr_TR_vs_TNR[[names]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR[[names]][clustered_results_ashr_TR_vs_TNR[[names]]$Status
                                                                                                                                                                                    == "Representative", ][1:10,],
                                                                       top_terms = NULL,
                                                                       plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      scale_size(range = c(0.1, 2)) +
      theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
            axis.text.y = element_text(color = "black", size = 3.5),
            axis.text.x = element_text(color = "black", size = 3.5),
            axis.title.x = element_text(size = 4.5, face = "bold"),
            legend.key.size = unit(1.5, units = "mm"),
            legend.spacing.y = unit(0.5, units = "mm"),
            legend.spacing.x = unit(0.5, units = "mm"),
            legend.title = element_text(size = 3.5, face = "bold"),
            legend.text = element_text(size = 2.75))+
      labs(title = paste0("Top 10 clustered ", names,
                          " terms enrichment dotplot - (", comparisons["ashr_TR_vs_TNR"], ")"))
    print(cluster_enrichment_dotplots_ashr_TR_vs_TNR[[names]])
    ggsave(filename = paste0(names, "_top10_dotplot_clustered.tiff"),
           path = paste0("Results/Pathways/ashr_TR_vs_TNR/pathfindR/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'tiff', units = "px",
           dpi = 700, compression = "lzw")
    dev.off()
  }
}

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(clustered_results_ashr_TR_vs_TNR)) {
  addWorksheet(wb, names(clustered_results_ashr_TR_vs_TNR)[i])
  writeData(wb, names(clustered_results_ashr_TR_vs_TNR)[i], clustered_results_ashr_TR_vs_TNR[[i]])
}
saveWorkbook(wb, file = "Results/Pathways/ashr_TR_vs_TNR/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2_names = names(clustered_results_ashr_TR_vs_TNR)
wb2 = createWorkbook()
for(i in 1:length(clustered_results_ashr_TR_vs_TNR)) {
  addWorksheet(wb2, paste0(wb2_names[i], " - rep"))
  writeData(wb2, paste0(wb2_names[i], " - rep"), 
            clustered_results_ashr_TR_vs_TNR[[i]][clustered_results_ashr_TR_vs_TNR[[i]]$Status
                                                  == "Representative", ])
}
saveWorkbook(wb2, file = "Results/Pathways/ashr_TR_vs_TNR/pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####

# # Defining a legend alignment function 
# align_legend <- function(p, hjust = 0.5) {
#   g <- cowplot::plot_to_gtable(p)
#   grobs <- g$grobs
#   legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
#   if (length(legend_index) == 0) {
#     warning("No legend found for this plot.")
#     return(g)
#   }
#   legend <- grobs[[legend_index]]
#   guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
#   if (length(guides_index) == 0) {
#     warning("No guides layout found in legend.")
#     return(g)
#   }
#   for (gi in guides_index) {
#     guides <- legend$grobs[[gi]]
#     if (length(guides$width) < 5) {
#       warning("Insufficient guide width elements for alignment.")
#       return(g)
#     }
#     spacing <- guides$width[5]
#     guides <- gtable::gtable_add_cols(guides, hjust * spacing, 1)
#     guides$widths[6] <- (1 - hjust) * spacing
#     title_index <- guides$layout$name == "title"
#     guides$layout$l[title_index] <- 2
#     legend$grobs[[gi]] <- guides
#   }
#   g$grobs[[legend_index]] <- legend
#   g
# }

term_gene_heatmaps_ashr_TR_vs_TNR = list()
term_gene_graphs_ashr_TR_vs_TNR = list()

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_TR_vs_TNR)[i]
  # term-gene heatmaps
  term_gene_heatmaps_ashr_TR_vs_TNR[[names]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[names]],
                                                             genes_df = pathf_input_ashr_TR_vs_TNR,
                                                             num_terms = 10,
                                                             use_description = TRUE,
                                                             low = "darkblue",
                                                             high = "red",
                                                             mid = "white",
                                                             pin_name_path = 'Biogrid')+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 4),
          axis.text.x = element_text(color = "black", size = 4,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names, 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_TR_vs_TNR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  #print(ggdraw(align_legend(term_gene_heatmaps_ashr_TR_vs_TNR[[names]], hjust = 0.5)))
  print(term_gene_heatmaps_ashr_TR_vs_TNR[[names]])
  ggsave(filename = paste0(names,
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Results/Pathways/ashr_TR_vs_TNR/pathfindR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_TR_vs_TNR[[names]] = term_gene_graph(result_df = pathfindR_outputs_ashr_TR_vs_TNR[[names]],
                                                         num_terms = 5,
                                                         use_description = TRUE,
                                                         node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 2, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names, 
                        " term - gene graph (",
                        comparisons["ashr_TR_vs_TNR"], ")"))
  term_gene_graphs_ashr_TR_vs_TNR[[names]]$layers = list(term_gene_graphs_ashr_TR_vs_TNR[[names]]$layers[[1]], 
                                                     term_gene_graphs_ashr_TR_vs_TNR[[names]]$layers[[2]], 
                                                     term_gene_graphs_ashr_TR_vs_TNR[[names]]$layers[[4]])
  print(term_gene_graphs_ashr_TR_vs_TNR[[names]])
  ggsave(filename = paste0(names,
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Results/Pathways/ashr_TR_vs_TNR/pathfindR/", 
                       names, "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_TR_vs_TNR = list()
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_TR_vs_TNR)[i]
  # UpSet plot
  UpSet_plots_ashr_TR_vs_TNR[[names]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[names]],
                                               genes_df = pathf_input_ashr_TR_vs_TNR,
                                               num_terms = 5,
                                               use_description = TRUE,
                                               low = "darkgreen",
                                               high = "darkred",
                                               mid = "white")+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 3),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                              combmatrix.panel.point.size = 0.8,
                              combmatrix.panel.line.size = 0.5)+
    labs(fill = "logFC",
         title = paste0("Top 5 ", names, 
                        " term - UpSet plot (",
                        comparisons["ashr_TR_vs_TNR"], ")"))
  tiff(paste0("Results/Pathways/ashr_TR_vs_TNR/pathfindR/", names, "/",
              names, "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_ashr_TR_vs_TNR[[names]])
  dev.off()
}

# ashr_LTR_vs_ETR #####
# Preparing a pathfindR loop for enrichment analysis
pathf_input_ashr_LTR_vs_ETR = ashr_dgea$ashr_LTR_vs_ETR %>%
  dplyr::select(Gene.Symbol, log2FoldChange, padj) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_ashr_LTR_vs_ETR = paste0("Results/Pathways/ashr_LTR_vs_ETR/pathfindR/", subdirs)
pathfindR_outputs_ashr_LTR_vs_ETR = list()

RNGversion("4.2.2")
set.seed(123)
for (i in 1:length(dirs_ashr_LTR_vs_ETR)){
  pathfindR_outputs_ashr_LTR_vs_ETR[[i]] = run_pathfindR(pathf_input_ashr_LTR_vs_ETR, gene_sets = subdirs[i],
                                                         p_val_threshold = p_val_threshold, convert2alias = FALSE,
                                                         output_dir = dirs_ashr_LTR_vs_ETR[i], min_gset_size = min_gset_size,
                                                         max_gset_size = max_gset_size, adj_method = adj_method,
                                                         enrichment_threshold = enrichment_threshold,
                                                         pin_name_path = pin_name_path, search_method = search_method,
                                                         iterations = niter, n_processes = n_processes)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_ashr_LTR_vs_ETR) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

# Use a custom fast approach
# Set up the number of cores to use
nonzero_res = which(lapply(pathfindR_outputs_ashr_LTR_vs_ETR, nrow) != 0)
cl <- makeCluster(length(nonzero_res))
registerDoParallel(cl)

timestamp()
clustered_results_ashr_LTR_vs_ETR <- foreach(i = nonzero_res, .packages = c("pathfindR", "fastcluster")) %dopar% {
  RNGversion("4.2.2")
  set.seed(123)
  source("R_Scripts/fast_pathfindR_hclust.R")
  
  # Perform clustering, handle errors
  result <- cluster_enriched_terms_fast(pathfindR_outputs_ashr_LTR_vs_ETR[[i]], method = "hierarchical", plot_clusters_graph = FALSE,
                                        use_description = TRUE, use_active_snw_genes = FALSE)
  if (is.character(result) && result == "hclust impossible") {
    return("hclust impossible")
  } else {
    return(result)
  }
}

timestamp() # ~1.5 mins
stopCluster(cl)
gc()

# Unlist
for (i in 1:length(clustered_results_ashr_LTR_vs_ETR)) {
  clustered_results_ashr_LTR_vs_ETR[[i]] = clustered_results_ashr_LTR_vs_ETR[[i]]$clustered_df
}
names(clustered_results_ashr_LTR_vs_ETR) = names(pathfindR_outputs_ashr_LTR_vs_ETR)[nonzero_res]

# Are there any null sets?
which(clustered_results_ashr_LTR_vs_ETR == "hclust impossible")

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_ashr_LTR_vs_ETR = pathfindR_outputs_ashr_LTR_vs_ETR
for (i in nonzero_res){
  wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_ashr_LTR_vs_ETR[[i]]$Term_Description, 
                                                                                      width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR = clustered_results_ashr_LTR_vs_ETR
for (i in 1:length(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR)){
  wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR[[i]]$Term_Description = stringr::str_wrap(clustered_results_ashr_LTR_vs_ETR[[i]]$Term_Description, 
                                                                                                width = 41)
}
rm(i)

enrichment_dotplots_ashr_LTR_vs_ETR = list()
cluster_enrichment_dotplots_ashr_LTR_vs_ETR = list()

# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_ETR)[i]
  # unclustered results
  enrichment_dotplots_ashr_LTR_vs_ETR[[names]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[names]],
                                                                  top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    scale_size(range = c(0.1, 2)) +
    theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 3.5),
          axis.title.x = element_text(size = 4.5, face = "bold"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 3.5, face = "bold"),
          legend.text = element_text(size = 2.75))+
    labs(title = paste0("Top 10 ", names,
                        " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_ETR"], ")"))
  print(enrichment_dotplots_ashr_LTR_vs_ETR[[names]])
  ggsave(filename = paste0(names, "_top10_dotplot.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_ETR/pathfindR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  if (names %in% names(clustered_results_ashr_LTR_vs_ETR)) {
    cluster_enrichment_dotplots_ashr_LTR_vs_ETR[[names]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR[[names]][clustered_results_ashr_LTR_vs_ETR[[names]]$Status
                                                                                                                                                     == "Representative", ][1:10,],
                                                                            top_terms = NULL,
                                                                            plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      scale_size(range = c(0.1, 2)) +
      theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
            axis.text.y = element_text(color = "black", size = 3.5),
            axis.text.x = element_text(color = "black", size = 3.5),
            axis.title.x = element_text(size = 4.5, face = "bold"),
            legend.key.size = unit(1.5, units = "mm"),
            legend.spacing.y = unit(0.5, units = "mm"),
            legend.spacing.x = unit(0.5, units = "mm"),
            legend.title = element_text(size = 3.5, face = "bold"),
            legend.text = element_text(size = 2.75))+
      labs(title = paste0("Top 10 clustered ", names,
                          " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_ETR"], ")"))
    print(cluster_enrichment_dotplots_ashr_LTR_vs_ETR[[names]])
    ggsave(filename = paste0(names, "_top10_dotplot_clustered.tiff"),
           path = paste0("Results/Pathways/ashr_LTR_vs_ETR/pathfindR/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'tiff', units = "px",
           dpi = 700, compression = "lzw")
    dev.off()
  }
}

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(clustered_results_ashr_LTR_vs_ETR)) {
  addWorksheet(wb, names(clustered_results_ashr_LTR_vs_ETR)[i])
  writeData(wb, names(clustered_results_ashr_LTR_vs_ETR)[i], clustered_results_ashr_LTR_vs_ETR[[i]])
}
saveWorkbook(wb, file = "Results/Pathways/ashr_LTR_vs_ETR/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2_names = names(clustered_results_ashr_LTR_vs_ETR)
wb2 = createWorkbook()
for(i in 1:length(clustered_results_ashr_LTR_vs_ETR)) {
  addWorksheet(wb2, paste0(wb2_names[i], " - rep"))
  writeData(wb2, paste0(wb2_names[i], " - rep"), 
            clustered_results_ashr_LTR_vs_ETR[[i]][clustered_results_ashr_LTR_vs_ETR[[i]]$Status
                                                   == "Representative", ])
}
saveWorkbook(wb2, file = "Results/Pathways/ashr_LTR_vs_ETR/pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_ashr_LTR_vs_ETR = list()
term_gene_graphs_ashr_LTR_vs_ETR = list()

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_ETR)[i]
  # term-gene heatmaps
  term_gene_heatmaps_ashr_LTR_vs_ETR[[names]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[names]],
                                                                  genes_df = pathf_input_ashr_LTR_vs_ETR,
                                                                  num_terms = 10,
                                                                  use_description = TRUE,
                                                                  low = "darkblue",
                                                                  high = "red",
                                                                  mid = "white",
                                                                  pin_name_path = 'Biogrid')+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 4),
          axis.text.x = element_text(color = "black", size = 4,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names, 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_LTR_vs_ETR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  #print(ggdraw(align_legend(term_gene_heatmaps_ashr_LTR_vs_ETR[[names]], hjust = 0.5)))
  print(term_gene_heatmaps_ashr_LTR_vs_ETR[[names]])
  ggsave(filename = paste0(names,
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_ETR/pathfindR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_LTR_vs_ETR[[names]] = term_gene_graph(result_df = pathfindR_outputs_ashr_LTR_vs_ETR[[names]],
                                                              num_terms = 5,
                                                              use_description = TRUE,
                                                              node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 2, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names, 
                        " term - gene graph (",
                        comparisons["ashr_LTR_vs_ETR"], ")"))
  term_gene_graphs_ashr_LTR_vs_ETR[[names]]$layers = list(term_gene_graphs_ashr_LTR_vs_ETR[[names]]$layers[[1]], 
                                                          term_gene_graphs_ashr_LTR_vs_ETR[[names]]$layers[[2]], 
                                                          term_gene_graphs_ashr_LTR_vs_ETR[[names]]$layers[[4]])
  print(term_gene_graphs_ashr_LTR_vs_ETR[[names]])
  ggsave(filename = paste0(names,
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_ETR/pathfindR/", 
                       names, "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_LTR_vs_ETR = list()
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_ETR)[i]
  # UpSet plot
  UpSet_plots_ashr_LTR_vs_ETR[[names]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[names]],
                                                    genes_df = pathf_input_ashr_LTR_vs_ETR,
                                                    num_terms = 5,
                                                    use_description = TRUE,
                                                    low = "darkgreen",
                                                    high = "darkred",
                                                    mid = "white")+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 3),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                              combmatrix.panel.point.size = 0.8,
                              combmatrix.panel.line.size = 0.5)+
    labs(fill = "logFC",
         title = paste0("Top 5 ", names, 
                        " term - UpSet plot (",
                        comparisons["ashr_LTR_vs_ETR"], ")"))
  tiff(paste0("Results/Pathways/ashr_LTR_vs_ETR/pathfindR/", names, "/",
              names, "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_ashr_LTR_vs_ETR[[names]])
  dev.off()
}

# ashr_LTR_vs_LTNR #####
# Preparing a pathfindR loop for enrichment analysis
pathf_input_ashr_LTR_vs_LTNR = ashr_dgea$ashr_LTR_vs_LTNR %>%
  dplyr::select(Gene.Symbol, log2FoldChange, padj) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_ashr_LTR_vs_LTNR = paste0("Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/", subdirs)
pathfindR_outputs_ashr_LTR_vs_LTNR = list()

RNGversion("4.2.2")
set.seed(123)
for (i in 1:length(dirs_ashr_LTR_vs_LTNR)){
  pathfindR_outputs_ashr_LTR_vs_LTNR[[i]] = run_pathfindR(pathf_input_ashr_LTR_vs_LTNR, gene_sets = subdirs[i],
                                                          p_val_threshold = p_val_threshold, convert2alias = FALSE,
                                                          output_dir = dirs_ashr_LTR_vs_LTNR[i], min_gset_size = min_gset_size,
                                                          max_gset_size = max_gset_size, adj_method = adj_method,
                                                          enrichment_threshold = enrichment_threshold,
                                                          pin_name_path = pin_name_path, search_method = search_method,
                                                          iterations = niter, n_processes = n_processes)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_ashr_LTR_vs_LTNR) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

# Use a custom fast approach
# Set up the number of cores to use
nonzero_res = which(lapply(pathfindR_outputs_ashr_LTR_vs_LTNR, nrow) != 0)
cl <- makeCluster(length(nonzero_res))
registerDoParallel(cl)

at1 = Sys.time()
clustered_results_ashr_LTR_vs_LTNR <- foreach(i = nonzero_res, .packages = c("pathfindR", "fastcluster")) %dopar% {
  RNGversion("4.2.2")
  set.seed(123)
  source("R_Scripts/fast_pathfindR_hclust.R")
  
  # Perform clustering, handle errors
  result <- cluster_enriched_terms_fast(pathfindR_outputs_ashr_LTR_vs_LTNR[[i]], method = "hierarchical", plot_clusters_graph = FALSE,
                                        use_description = TRUE, use_active_snw_genes = FALSE)
  if (is.character(result) && result == "hclust impossible") {
    return("hclust impossible")
  } else {
    return(result)
  }
}

at2 = Sys.time() # <0.5 min
dat2 = at2-at1
stopCluster(cl)
gc()

# Unlist
for (i in 1:length(clustered_results_ashr_LTR_vs_LTNR)) {
  clustered_results_ashr_LTR_vs_LTNR[[i]] = clustered_results_ashr_LTR_vs_LTNR[[i]]$clustered_df
}
names(clustered_results_ashr_LTR_vs_LTNR) = names(pathfindR_outputs_ashr_LTR_vs_LTNR)[nonzero_res]

# Are there any null sets?
which(clustered_results_ashr_LTR_vs_LTNR == "hclust impossible")

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR = pathfindR_outputs_ashr_LTR_vs_LTNR
for (i in nonzero_res){
  wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_ashr_LTR_vs_LTNR[[i]]$Term_Description, 
                                                                                       width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR = clustered_results_ashr_LTR_vs_LTNR
for (i in 1:length(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR)){
  wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR[[i]]$Term_Description = stringr::str_wrap(clustered_results_ashr_LTR_vs_LTNR[[i]]$Term_Description, 
                                                                                                 width = 41)
}
rm(i)

enrichment_dotplots_ashr_LTR_vs_LTNR = list()
cluster_enrichment_dotplots_ashr_LTR_vs_LTNR = list()

# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i]
  # unclustered results
  enrichment_dotplots_ashr_LTR_vs_LTNR[[names]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]],
                                                                   top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    scale_size(range = c(0.1, 2)) +
    theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 3.5),
          axis.title.x = element_text(size = 4.5, face = "bold"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 3.5, face = "bold"),
          legend.text = element_text(size = 2.75))+
    labs(title = paste0("Top 10 ", names,
                        " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_LTNR"], ")"))
  print(enrichment_dotplots_ashr_LTR_vs_LTNR[[names]])
  ggsave(filename = paste0(names, "_top10_dotplot.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  if (names %in% names(clustered_results_ashr_LTR_vs_LTNR)) {
    cluster_enrichment_dotplots_ashr_LTR_vs_LTNR[[names]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]][clustered_results_ashr_LTR_vs_LTNR[[names]]$Status
                                                                                                                                                       == "Representative", ][1:10,],
                                                                             top_terms = NULL,
                                                                             plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      scale_size(range = c(0.1, 2)) +
      theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
            axis.text.y = element_text(color = "black", size = 3.5),
            axis.text.x = element_text(color = "black", size = 3.5),
            axis.title.x = element_text(size = 4.5, face = "bold"),
            legend.key.size = unit(1.5, units = "mm"),
            legend.spacing.y = unit(0.5, units = "mm"),
            legend.spacing.x = unit(0.5, units = "mm"),
            legend.title = element_text(size = 3.5, face = "bold"),
            legend.text = element_text(size = 2.75))+
      labs(title = paste0("Top 10 clustered ", names,
                          " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_LTNR"], ")"))
    print(cluster_enrichment_dotplots_ashr_LTR_vs_LTNR[[names]])
    ggsave(filename = paste0(names, "_top10_dotplot_clustered.tiff"),
           path = paste0("Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'tiff', units = "px",
           dpi = 700, compression = "lzw")
    dev.off()
  }
}

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(clustered_results_ashr_LTR_vs_LTNR)) {
  addWorksheet(wb, names(clustered_results_ashr_LTR_vs_LTNR)[i])
  writeData(wb, names(clustered_results_ashr_LTR_vs_LTNR)[i], clustered_results_ashr_LTR_vs_LTNR[[i]])
}
saveWorkbook(wb, file = "Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2_names = names(clustered_results_ashr_LTR_vs_LTNR)
wb2 = createWorkbook()
for(i in 1:length(clustered_results_ashr_LTR_vs_LTNR)) {
  addWorksheet(wb2, paste0(wb2_names[i], " - rep"))
  writeData(wb2, paste0(wb2_names[i], " - rep"), 
            clustered_results_ashr_LTR_vs_LTNR[[i]][clustered_results_ashr_LTR_vs_LTNR[[i]]$Status
                                                    == "Representative", ])
}
saveWorkbook(wb2, file = "Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_ashr_LTR_vs_LTNR = list()
term_gene_graphs_ashr_LTR_vs_LTNR = list()

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i]
  # term-gene heatmaps
  term_gene_heatmaps_ashr_LTR_vs_LTNR[[names]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]],
                                                                   genes_df = pathf_input_ashr_LTR_vs_LTNR,
                                                                   num_terms = 10,
                                                                   use_description = TRUE,
                                                                   low = "darkblue",
                                                                   high = "red",
                                                                   mid = "white",
                                                                   pin_name_path = 'Biogrid')+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 4),
          axis.text.x = element_text(color = "black", size = 4,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names, 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_LTR_vs_LTNR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  #print(ggdraw(align_legend(term_gene_heatmaps_ashr_LTR_vs_LTNR[[names]], hjust = 0.5)))
  print(term_gene_heatmaps_ashr_LTR_vs_LTNR[[names]])
  ggsave(filename = paste0(names,
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_LTR_vs_LTNR[[names]] = term_gene_graph(result_df = pathfindR_outputs_ashr_LTR_vs_LTNR[[names]],
                                                               num_terms = 5,
                                                               use_description = TRUE,
                                                               node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 2, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names, 
                        " term - gene graph (",
                        comparisons["ashr_LTR_vs_LTNR"], ")"))
  term_gene_graphs_ashr_LTR_vs_LTNR[[names]]$layers = list(term_gene_graphs_ashr_LTR_vs_LTNR[[names]]$layers[[1]], 
                                                           term_gene_graphs_ashr_LTR_vs_LTNR[[names]]$layers[[2]], 
                                                           term_gene_graphs_ashr_LTR_vs_LTNR[[names]]$layers[[4]])
  print(term_gene_graphs_ashr_LTR_vs_LTNR[[names]])
  ggsave(filename = paste0(names,
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/", 
                       names, "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_LTR_vs_LTNR = list()
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i]
  # UpSet plot
  UpSet_plots_ashr_LTR_vs_LTNR[[names]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]],
                                                     genes_df = pathf_input_ashr_LTR_vs_LTNR,
                                                     num_terms = 5,
                                                     use_description = TRUE,
                                                     low = "darkgreen",
                                                     high = "darkred",
                                                     mid = "white")+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 3),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                              combmatrix.panel.point.size = 0.8,
                              combmatrix.panel.line.size = 0.5)+
    labs(fill = "logFC",
         title = paste0("Top 5 ", names, 
                        " term - UpSet plot (",
                        comparisons["ashr_LTR_vs_LTNR"], ")"))
  tiff(paste0("Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/", names, "/",
              names, "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_ashr_LTR_vs_LTNR[[names]])
  dev.off()
}

# ashr_ETR_vs_ETNR #####
# Preparing a pathfindR loop for enrichment analysis
pathf_input_ashr_ETR_vs_ETNR = ashr_dgea$ashr_ETR_vs_ETNR %>%
  dplyr::select(Gene.Symbol, log2FoldChange, padj) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_ashr_ETR_vs_ETNR = paste0("Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/", subdirs)
pathfindR_outputs_ashr_ETR_vs_ETNR = list()

RNGversion("4.2.2")
set.seed(123)
for (i in 1:length(dirs_ashr_ETR_vs_ETNR)){
  pathfindR_outputs_ashr_ETR_vs_ETNR[[i]] = run_pathfindR(pathf_input_ashr_ETR_vs_ETNR, gene_sets = subdirs[i],
                                                          p_val_threshold = p_val_threshold, convert2alias = FALSE,
                                                          output_dir = dirs_ashr_ETR_vs_ETNR[i], min_gset_size = min_gset_size,
                                                          max_gset_size = max_gset_size, adj_method = adj_method,
                                                          enrichment_threshold = enrichment_threshold,
                                                          pin_name_path = pin_name_path, search_method = search_method,
                                                          iterations = niter, n_processes = n_processes)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_ashr_ETR_vs_ETNR) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

# Use a custom fast approach
# Set up the number of cores to use
nonzero_res = which(lapply(pathfindR_outputs_ashr_ETR_vs_ETNR, nrow) != 0)
cl <- makeCluster(length(nonzero_res))
registerDoParallel(cl)

bt1 = Sys.time()
clustered_results_ashr_ETR_vs_ETNR <- foreach(i = nonzero_res, .packages = c("pathfindR", "fastcluster")) %dopar% {
  RNGversion("4.2.2")
  set.seed(123)
  source("R_Scripts/fast_pathfindR_hclust.R")
  
  # Perform clustering, handle errors
  result <- cluster_enriched_terms_fast(pathfindR_outputs_ashr_ETR_vs_ETNR[[i]], method = "hierarchical", plot_clusters_graph = FALSE,
                                        use_description = TRUE, use_active_snw_genes = FALSE)
  if (is.character(result) && result == "hclust impossible") {
    return("hclust impossible")
  } else {
    return(result)
  }
}

bt2 = Sys.time()
dbt2 = bt2-bt1 # <0.5 mins
stopCluster(cl)
gc()

# Unlist
for (i in 1:length(clustered_results_ashr_ETR_vs_ETNR)) {
  clustered_results_ashr_ETR_vs_ETNR[[i]] = clustered_results_ashr_ETR_vs_ETNR[[i]]$clustered_df
}
names(clustered_results_ashr_ETR_vs_ETNR) = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[nonzero_res]

# Are there any null sets?
which(clustered_results_ashr_ETR_vs_ETNR == "hclust impossible")

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR = pathfindR_outputs_ashr_ETR_vs_ETNR
for (i in nonzero_res){
  wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_ashr_ETR_vs_ETNR[[i]]$Term_Description, 
                                                                                       width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR = clustered_results_ashr_ETR_vs_ETNR
for (i in 1:length(wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR)){
  wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[i]]$Term_Description = stringr::str_wrap(clustered_results_ashr_ETR_vs_ETNR[[i]]$Term_Description, 
                                                                                                 width = 41)
}
rm(i)

enrichment_dotplots_ashr_ETR_vs_ETNR = list()
cluster_enrichment_dotplots_ashr_ETR_vs_ETNR = list()

# Producing dotplots with the results
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i]
  # unclustered results
  enrichment_dotplots_ashr_ETR_vs_ETNR[[names]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]],
                                                                   top_terms = 10)+
    scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
    scale_size(range = c(0.1, 2)) +
    theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 3.5),
          axis.title.x = element_text(size = 4.5, face = "bold"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 3.5, face = "bold"),
          legend.text = element_text(size = 2.75))+
    labs(title = paste0("Top 10 ", names,
                        " terms enrichment dotplot - (", comparisons["ashr_ETR_vs_ETNR"], ")"))
  print(enrichment_dotplots_ashr_ETR_vs_ETNR[[names]])
  ggsave(filename = paste0(names, "_top10_dotplot.tiff"),
         path = paste0("Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  if (names %in% names(clustered_results_ashr_ETR_vs_ETNR)) {
    cluster_enrichment_dotplots_ashr_ETR_vs_ETNR[[names]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]][clustered_results_ashr_ETR_vs_ETNR[[names]]$Status
                                                                                                                                                       == "Representative", ][1:10,],
                                                                             top_terms = NULL,
                                                                             plot_by_cluster = TRUE)+
      scale_color_gradient(low = "#fca4a4", high = "#fc0303")+
      scale_size(range = c(0.1, 2)) +
      theme(plot.title = element_text(size = 5, face = "bold", vjust = 2, hjust = 0.5),
            axis.text.y = element_text(color = "black", size = 3.5),
            axis.text.x = element_text(color = "black", size = 3.5),
            axis.title.x = element_text(size = 4.5, face = "bold"),
            legend.key.size = unit(1.5, units = "mm"),
            legend.spacing.y = unit(0.5, units = "mm"),
            legend.spacing.x = unit(0.5, units = "mm"),
            legend.title = element_text(size = 3.5, face = "bold"),
            legend.text = element_text(size = 2.75))+
      labs(title = paste0("Top 10 clustered ", names,
                          " terms enrichment dotplot - (", comparisons["ashr_ETR_vs_ETNR"], ")"))
    print(cluster_enrichment_dotplots_ashr_ETR_vs_ETNR[[names]])
    ggsave(filename = paste0(names, "_top10_dotplot_clustered.tiff"),
           path = paste0("Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/", 
                         names, "/"),
           width = 2880, height = 1620, device = 'tiff', units = "px",
           dpi = 700, compression = "lzw")
    dev.off()
  }
}

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(clustered_results_ashr_ETR_vs_ETNR)) {
  addWorksheet(wb, names(clustered_results_ashr_ETR_vs_ETNR)[i])
  writeData(wb, names(clustered_results_ashr_ETR_vs_ETNR)[i], clustered_results_ashr_ETR_vs_ETNR[[i]])
}
saveWorkbook(wb, file = "Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2_names = names(clustered_results_ashr_ETR_vs_ETNR)
wb2 = createWorkbook()
for(i in 1:length(clustered_results_ashr_ETR_vs_ETNR)) {
  addWorksheet(wb2, paste0(wb2_names[i], " - rep"))
  writeData(wb2, paste0(wb2_names[i], " - rep"), 
            clustered_results_ashr_ETR_vs_ETNR[[i]][clustered_results_ashr_ETR_vs_ETNR[[i]]$Status
                                                    == "Representative", ])
}
saveWorkbook(wb2, file = "Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_ashr_ETR_vs_ETNR = list()
term_gene_graphs_ashr_ETR_vs_ETNR = list()

for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i]
  # term-gene heatmaps
  term_gene_heatmaps_ashr_ETR_vs_ETNR[[names]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]],
                                                                   genes_df = pathf_input_ashr_ETR_vs_ETNR,
                                                                   num_terms = 10,
                                                                   use_description = TRUE,
                                                                   low = "darkblue",
                                                                   high = "red",
                                                                   mid = "white",
                                                                   pin_name_path = 'Biogrid')+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 4),
          axis.text.x = element_text(color = "black", size = 4,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names, 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_ETR_vs_ETNR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  #print(ggdraw(align_legend(term_gene_heatmaps_ashr_ETR_vs_ETNR[[names]], hjust = 0.5)))
  print(term_gene_heatmaps_ashr_ETR_vs_ETNR[[names]])
  ggsave(filename = paste0(names,
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/", 
                       names, "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_ETR_vs_ETNR[[names]] = term_gene_graph(result_df = pathfindR_outputs_ashr_ETR_vs_ETNR[[names]],
                                                               num_terms = 5,
                                                               use_description = TRUE,
                                                               node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 2, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 5, face = "bold", hjust = 0.5),
          legend.text = element_text(size = 4),
          #legend.title = element_text(hjust = 0.5),
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names, 
                        " term - gene graph (",
                        comparisons["ashr_ETR_vs_ETNR"], ")"))
  term_gene_graphs_ashr_ETR_vs_ETNR[[names]]$layers = list(term_gene_graphs_ashr_ETR_vs_ETNR[[names]]$layers[[1]], 
                                                           term_gene_graphs_ashr_ETR_vs_ETNR[[names]]$layers[[2]], 
                                                           term_gene_graphs_ashr_ETR_vs_ETNR[[names]]$layers[[4]])
  print(term_gene_graphs_ashr_ETR_vs_ETNR[[names]])
  ggsave(filename = paste0(names,
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/", 
                       names, "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_ETR_vs_ETNR = list()
for (i in nonzero_res){
  names = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i]
  # UpSet plot
  UpSet_plots_ashr_ETR_vs_ETNR[[names]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]],
                                                     genes_df = pathf_input_ashr_ETR_vs_ETNR,
                                                     num_terms = 5,
                                                     use_description = TRUE,
                                                     low = "darkgreen",
                                                     high = "darkred",
                                                     mid = "white")+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 3),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    ggupset::theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                              combmatrix.panel.point.size = 0.8,
                              combmatrix.panel.line.size = 0.5)+
    labs(fill = "logFC",
         title = paste0("Top 5 ", names, 
                        " term - UpSet plot (",
                        comparisons["ashr_ETR_vs_ETNR"], ")"))
  tiff(paste0("Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/", names, "/",
              names, "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_ashr_ETR_vs_ETNR[[names]])
  dev.off()
}

# Per sample pathway scores #####
dds = readRDS("Resources_and_environments/final_dds.rds")
coldata = colData(dds)
recurrent = coldata$Patient.ID[coldata$Recurrence == "Recurrence"]
non_recurrent = coldata$Patient.ID[coldata$Recurrence == "No recurrence"]
z_exprs_gs = readRDS("Resources_and_environments/log2.norm.counts.plus1.rds")

# ashr_TR_vs_TNR #####
nonzero_res = which(lapply(pathfindR_outputs_ashr_TR_vs_TNR, nrow) != 0)
# Raw output
raw_score_matrices_ashr_TR_vs_TNR = list()
for (p in nonzero_res){
  names = names(pathfindR_outputs_ashr_TR_vs_TNR)[p]
  if (nrow(wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[p]]) == 0) {
    next
  } else if (nrow(wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[p]]) < 10) {
    rows = 1:nrow(wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[p]])
  } else {
    rows = c(1:10)
  }
  tiff(paste0(dirs_ashr_TR_vs_TNR[p], "/",
              names, "_per_sample_pathway_scores.tiff"), 
       width = 5760*1.5, height = 2160*1.5, res = 700, compression = "lzw")
  raw_score_matrices_ashr_TR_vs_TNR[[names]] = score_terms(enrichment_table = wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[names]][rows,],
                                        exp_mat = z_exprs_gs,
                                        cases = recurrent,
                                        use_description = TRUE, # default FALSE
                                        label_samples = TRUE, # default = TRUE
                                        case_title = "Recurrent",  # default = "Case"
                                        control_title = "Non-recurrent", # default = "Control"
                                        low = "darkblue", # default = "green"
                                        mid = "white", # default = "black"
                                        high = "darkred") # default = "red"
  dev.off()
}

# Clustered terms
clustered_score_matrices_ashr_TR_vs_TNR = list()
for (p in 1:length(wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR)){
  if (nrow(wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR[[p]]) == 0){
    next
  } else if (nrow(wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR[[p]][wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR[[p]]$Status == "Representative",]) < 10) {
    rows = 1:nrow(wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR[[p]][wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR[[p]]$Status == "Representative",])
  } else {
    rows = c(1:10)
  }
  names = names(wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR)[p]
  tiff(paste0(dirs_ashr_TR_vs_TNR[p], "/",
              names, "_per_sample_pathway_scores_clustered.tiff"), 
       width = 5760*1.5, height = 2160*1.5, res = 700, compression = "lzw")
  clustered_score_matrices_ashr_TR_vs_TNR[[names]] = score_terms(enrichment_table = wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR[[names]][wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR[[names]]$Status == "Representative",][rows,], 
                                              exp_mat = z_exprs_gs,
                                              cases = recurrent,
                                              use_description = TRUE, # default FALSE
                                              label_samples = TRUE, # default = TRUE
                                              case_title = "Recurrent",  # default = "Case"
                                              control_title = "Non-recurrent", # default = "Control"
                                              low = "darkblue", # default = "green"
                                              mid = "white", # default = "black"
                                              high = "darkred") # default = "red"
  dev.off()
}

# ashr_LTR_vs_ETR #####
nonzero_res = which(lapply(pathfindR_outputs_ashr_LTR_vs_ETR, nrow) != 0)
# Raw output
raw_score_matrices_ashr_LTR_vs_ETR = list()
for (p in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_ETR)[p]
  if (nrow(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[p]]) == 0) {
    next
  } else if (nrow(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[p]]) < 10) {
    rows = 1:nrow(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[p]])
  } else {
    rows = c(1:10)
  }
  tiff(paste0(dirs_ashr_LTR_vs_ETR[p], "/",
              names, "_per_sample_pathway_scores.tiff"), 
       width = 5760*1.5, height = 2160*1.5, res = 700, compression = "lzw")
  raw_score_matrices_ashr_LTR_vs_ETR[[names]] = score_terms(enrichment_table = wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[names]][rows,],
                                                            exp_mat = z_exprs_gs,
                                                            cases = recurrent,
                                                            use_description = TRUE, # default FALSE
                                                            label_samples = TRUE, # default = TRUE
                                                            case_title = "Recurrent",  # default = "Case"
                                                            control_title = "Non-recurrent", # default = "Control"
                                                            low = "darkblue", # default = "green"
                                                            mid = "white", # default = "black"
                                                            high = "darkred") # default = "red"
  dev.off()
}

# Clustered terms
clustered_score_matrices_ashr_LTR_vs_ETR = list()
for (p in 1:length(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR)){
  if (nrow(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR[[p]]) == 0){
    next
  } else if (nrow(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR[[p]][wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR[[p]]$Status == "Representative",]) < 10) {
    rows = 1:nrow(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR[[p]][wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR[[p]]$Status == "Representative",])
  } else {
    rows = c(1:10)
  }
  names = names(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR)[p]
  tiff(paste0(dirs_ashr_LTR_vs_ETR[p], "/",
              names, "_per_sample_pathway_scores_clustered.tiff"), 
       width = 5760*1.5, height = 2160*1.5, res = 700, compression = "lzw")
  clustered_score_matrices_ashr_LTR_vs_ETR[[names]] = score_terms(enrichment_table = wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR[[names]][wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR[[names]]$Status == "Representative",][rows,], 
                                                                  exp_mat = z_exprs_gs,
                                                                  cases = recurrent,
                                                                  use_description = TRUE, # default FALSE
                                                                  label_samples = TRUE, # default = TRUE
                                                                  case_title = "Recurrent",  # default = "Case"
                                                                  control_title = "Non-recurrent", # default = "Control"
                                                                  low = "darkblue", # default = "green"
                                                                  mid = "white", # default = "black"
                                                                  high = "darkred") # default = "red"
  dev.off()
}

# ashr_ETR_vs_ETNR #####
nonzero_res = which(lapply(pathfindR_outputs_ashr_ETR_vs_ETNR, nrow) != 0)
# Raw output
raw_score_matrices_ashr_ETR_vs_ETNR = list()
for (p in nonzero_res){
  names = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[p]
  if (nrow(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[p]]) == 0) {
    next
  } else if (nrow(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[p]]) < 10) {
    rows = 1:nrow(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[p]])
  } else {
    rows = c(1:10)
  }
  tiff(paste0(dirs_ashr_ETR_vs_ETNR[p], "/",
              names, "_per_sample_pathway_scores.tiff"), 
       width = 5760*1.5, height = 2160*1.5, res = 700, compression = "lzw")
  raw_score_matrices_ashr_ETR_vs_ETNR[[names]] = score_terms(enrichment_table = wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]][rows,],
                                                             exp_mat = z_exprs_gs,
                                                             cases = recurrent,
                                                             use_description = TRUE, # default FALSE
                                                             label_samples = TRUE, # default = TRUE
                                                             case_title = "Recurrent",  # default = "Case"
                                                             control_title = "Non-recurrent", # default = "Control"
                                                             low = "darkblue", # default = "green"
                                                             mid = "white", # default = "black"
                                                             high = "darkred") # default = "red"
  dev.off()
}

# Clustered terms
clustered_score_matrices_ashr_ETR_vs_ETNR = list()
for (p in 1:length(wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR)){
  if (nrow(wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[p]]) == 0){
    next
  } else if (nrow(wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[p]][wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[p]]$Status == "Representative",]) < 10) {
    rows = 1:nrow(wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[p]][wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[p]]$Status == "Representative",])
  } else {
    rows = c(1:10)
  }
  names = names(wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR)[p]
  tiff(paste0(dirs_ashr_ETR_vs_ETNR[p], "/",
              names, "_per_sample_pathway_scores_clustered.tiff"), 
       width = 5760*1.5, height = 2160*1.5, res = 700, compression = "lzw")
  clustered_score_matrices_ashr_ETR_vs_ETNR[[names]] = score_terms(enrichment_table = wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]][wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[names]]$Status == "Representative",][rows,], 
                                                                   exp_mat = z_exprs_gs,
                                                                   cases = recurrent,
                                                                   use_description = TRUE, # default FALSE
                                                                   label_samples = TRUE, # default = TRUE
                                                                   case_title = "Recurrent",  # default = "Case"
                                                                   control_title = "Non-recurrent", # default = "Control"
                                                                   low = "darkblue", # default = "green"
                                                                   mid = "white", # default = "black"
                                                                   high = "darkred") # default = "red"
  dev.off()
}

# ashr_LTR_vs_LTNR #####
nonzero_res = which(lapply(pathfindR_outputs_ashr_LTR_vs_LTNR, nrow) != 0)
# Raw output
raw_score_matrices_ashr_LTR_vs_LTNR = list()
for (p in nonzero_res){
  names = names(pathfindR_outputs_ashr_LTR_vs_LTNR)[p]
  if (nrow(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[p]]) == 0) {
    next
  } else if (nrow(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[p]]) < 10) {
    rows = 1:nrow(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[p]])
  } else {
    rows = c(1:10)
  }
  tiff(paste0(dirs_ashr_LTR_vs_LTNR[p], "/",
              names, "_per_sample_pathway_scores.tiff"), 
       width = 5760*1.5, height = 2160*1.5, res = 700, compression = "lzw")
  raw_score_matrices_ashr_LTR_vs_LTNR[[names]] = score_terms(enrichment_table = wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]][rows,],
                                                             exp_mat = z_exprs_gs,
                                                             cases = recurrent,
                                                             use_description = TRUE, # default FALSE
                                                             label_samples = TRUE, # default = TRUE
                                                             case_title = "Recurrent",  # default = "Case"
                                                             control_title = "Non-recurrent", # default = "Control"
                                                             low = "darkblue", # default = "green"
                                                             mid = "white", # default = "black"
                                                             high = "darkred") # default = "red"
  dev.off()
}

# Clustered terms
clustered_score_matrices_ashr_LTR_vs_LTNR = list()
for (p in 1:length(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR)){
  if (nrow(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR[[p]]) == 0){
    next
  } else if (nrow(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR[[p]][wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR[[p]]$Status == "Representative",]) < 10) {
    rows = 1:nrow(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR[[p]][wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR[[p]]$Status == "Representative",])
  } else {
    rows = c(1:10)
  }
  names = names(wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR)[p]
  tiff(paste0(dirs_ashr_LTR_vs_LTNR[p], "/",
              names, "_per_sample_pathway_scores_clustered.tiff"), 
       width = 5760*1.5, height = 2160*1.5, res = 700, compression = "lzw")
  clustered_score_matrices_ashr_LTR_vs_LTNR[[names]] = score_terms(enrichment_table = wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]][wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR[[names]]$Status == "Representative",][rows,], 
                                                                   exp_mat = z_exprs_gs,
                                                                   cases = recurrent,
                                                                   use_description = TRUE, # default FALSE
                                                                   label_samples = TRUE, # default = TRUE
                                                                   case_title = "Recurrent",  # default = "Case"
                                                                   control_title = "Non-recurrent", # default = "Control"
                                                                   low = "darkblue", # default = "green"
                                                                   mid = "white", # default = "black"
                                                                   high = "darkred") # default = "red"
  dev.off()
}

# Pathway results processing for tabular presentation #####
# ORA
ora_sheets_names = c("ORA_GO", "ORA_KEGG", "ORA_WP", "ORA_Reactome",
                     "ORA_DO", "ORA_NCG", "ORA_DGN")
ora_names = c("Gene Ontology", "KEGG", "WikiPathways", "Reactome",
              "Disease Ontology", "Network Cancer Gene (NCG)", "Disease Gene Network (DisGeNET)")
for (j in 1:length(ora_results)) {
  for (i in 1:length(ora_names)) {
    if (ora_names[i] %in% names(ora_results[[j]])) {
      raw = read.xlsx(paste0("Results/Pathways/", names(ora_results)[j],
                             "/ORA/ORA_output_", names(ora_results)[j], ".xlsx"), 
                      sheet = ora_sheets_names[i])
      raw = tidyr::separate(raw, col = GeneRatio, sep = "/", into = c("num", "den"))
      raw$GeneRatio = as.numeric(raw$num)/as.numeric(raw$den)
      top10max = raw %>% dplyr::filter(p.adjust < 0.05) %>% 
        dplyr::arrange(p.adjust, desc(GeneRatio)) %>% 
        slice_head(n = 10)
      colnames(top10max)[1] = "Category"
      top10max = top10max %>% dplyr::select(Category, Description, GeneRatio, p.adjust, geneID)
      top10max$geneID = stringr::str_wrap(gsub("/", ", ", top10max$geneID), 
                                          width = 40)
      write.xlsx(top10max, paste0("Results/Pathways/", names(ora_results)[j],
                                  "/ORA/top_ORA_output_", names(ora_results)[j], "_",
                                  ora_names[i], ".xlsx"),
                 overwrite =TRUE)
    } else {
      next
    }
  }
}

# GSEA
gsea_sheets_names = c("GSEA_GO", "GSEA_KEGG", "GSEA_WP", "GSEA_Reactome",
                      "GSEA_DO", "GSEA_NCG", "GSEA_DGN")
gsea_names = c("Gene Ontology", "KEGG", "WikiPathways", "Reactome",
               "Disease Ontology", "Network Cancer Gene (NCG)", "Disease Gene Network (DisGeNET)")
for (j in 1:length(gsea_results)) {
  for (i in 1:length(gsea_names)) {
    if (gsea_names[i] %in% names(gsea_results[[j]])) {
      raw = read.xlsx(paste0("Results/Pathways/", names(gsea_results)[j],
                             "/GSEA/GSEA_output_", names(gsea_results)[j], ".xlsx"), 
                      sheet = gsea_sheets_names[i])
      top10max = raw %>% dplyr::filter(p.adjust < 0.05) %>% 
        dplyr::arrange(p.adjust, desc(abs(NES))) %>% 
        slice_head(n = 10)
      colnames(top10max)[1] = "Category"
      top10max = top10max %>% dplyr::select(Category, Description, 
                                  `Normalized Enrichment Score` = NES, p.adjust, 
                                  geneID = core_enrichment)
      top10max$geneID = stringr::str_wrap(gsub("/", ", ", top10max$geneID), 
                                          width = 40)
      write.xlsx(top10max, paste0("Results/Pathways/", names(gsea_results)[j],
                                  "/GSEA/top_GSEA_output_", names(gsea_results)[j], "_",
                                  gsea_names[i], ".xlsx"),
                 overwrite =TRUE)
    } else {
      next
    }
  }
}

# pathfindR
pathfindR_sheets_names = c("Reactome", "BioCarta", "GO-BP", 
                           "GO-CC", "GO-MF", "KEGG")
pathfindR_results = list(pathfindR_outputs_ashr_TR_vs_TNR,
                         pathfindR_outputs_ashr_LTR_vs_ETR,
                         pathfindR_outputs_ashr_ETR_vs_ETNR,
                         pathfindR_outputs_ashr_LTR_vs_LTNR)
names(pathfindR_results) = names(ashr_dgea)

for (j in 1:length(pathfindR_results)) {
  for (i in 1:length(pathfindR_sheets_names)) {
    if (nrow(pathfindR_results[[j]][[pathfindR_sheets_names[[i]]]]) > 0) {
      raw = read.xlsx(paste0("Results/Pathways/", names(pathfindR_results)[j],
                             "/pathfindR/Comprehensive_pathfindR_output.xlsx"), 
                      sheet = pathfindR_sheets_names[i])
      top10max_comprehensive = raw %>% dplyr::filter(lowest_p < 0.05) %>% 
        # dplyr::arrange(lowest_p, desc(abs(Fold_Enrichment))) %>% 
        slice_head(n = 100) %>% 
        dplyr::select(Description = Term_Description, `Fold enrichment` = Fold_Enrichment,
                      p.adjust = lowest_p, `Up-regulated` = Up_regulated, 
                      `Down-regulated` = Down_regulated, Status)
      top10max_representative = raw %>% dplyr::filter(lowest_p < 0.05 & Status == "Representative") %>% 
        dplyr::arrange(lowest_p, desc(abs(Fold_Enrichment))) %>% 
        slice_head(n = 10) %>% 
        dplyr::select(Description = Term_Description, `Fold enrichment` = Fold_Enrichment,
                      p.adjust = lowest_p, `Up-regulated` = Up_regulated, 
                      `Down-regulated` = Down_regulated, Status)
      write.xlsx(top10max_comprehensive, paste0("Results/Pathways/", names(pathfindR_results)[j],
                                                "/pathfindR/top_Comprehensive_pathfindR_output_",
                                                pathfindR_sheets_names[i], ".xlsx"),
                 overwrite =TRUE)
      write.xlsx(top10max_comprehensive, paste0("Results/Pathways/", names(pathfindR_results)[j],
                                                "/pathfindR/top_Representative_pathfindR_output_",
                                                pathfindR_sheets_names[i], ".xlsx"),
                 overwrite =TRUE)
    } else {
      next
    }
  }
}

rm(raw, top10max, top10max_comprehensive, top10max_representative, i, j); gc()

# Export a session info text file for reproducibility
writeLines(capture.output(sessionInfo()), 
           paste0("Resources_and_environments/Pathways_sessionInfo.txt"))

# Save environment for easy loading for DGEA later
save.image("Resources_and_environments/Pathways.RData")
