# This script uses the lists produced by DGEA for the purpose of ORA
# (Over-representation analysis) and GSEA
# (Gene Set Enrichment Analysis)

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(openxlsx)
library(ReactomePA)
library(org.Hs.eg.db)
library(pathfindR)
library(cowplot)

# Import data (6 lists from DGEA)
ashr_dgea = list()
for (i in 1:4) {
  ashr_dgea[[i]] = openxlsx::read.xlsx("Results/DGEA/DGEA_ashr.xlsx", sheet = i)
}
rm(i); gc()
 
names(ashr_dgea) = paste0("ashr_", 
                          c("TR_vs_TNR", "LTR_vs_ETR", 
                            "ETR_vs_ETNR", "LTR_vs_LTNR"))

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

rm(i); gc()

# Pathway analysis with pathfindR #####

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
                                                        p_val_threshold = 0.05, convert2alias = FALSE,
                                                        output_dir = dirs_ashr_TR_vs_TNR[i], min_gset_size = 10,
                                                        max_gset_size = 300, adj_method = 'fdr',
                                                        enrichment_threshold = 0.05,
                                                        pin_name_path = 'Biogrid', search_method = 'GR',
                                                        iterations = 10, n_processes = 10)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_ashr_TR_vs_TNR) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

# Clustering only for Reactome and GO-BP

cluster_names = subdirs
RNGversion("4.2.2")
set.seed(123)
clustered_results_ashr_TR_vs_TNR = list()
for (i in c(1, 3)){
  clustered_results_ashr_TR_vs_TNR[[i]] = cluster_enriched_terms(pathfindR_outputs_ashr_TR_vs_TNR[[cluster_names[[i]]]],
                                                                 method = "hierarchical",
                                                                 use_description = TRUE,
                                                                 plot_clusters_graph = FALSE)
}

# Add the original info to the non-clustered slots
clustered_results_ashr_TR_vs_TNR[[2]] = pathfindR_outputs_ashr_TR_vs_TNR[["BioCarta"]]
clustered_results_ashr_TR_vs_TNR[[4]] = pathfindR_outputs_ashr_TR_vs_TNR[["GO-CC"]]
clustered_results_ashr_TR_vs_TNR[[4]]$Cluster = 1
clustered_results_ashr_TR_vs_TNR[[4]]$Status = "Member"
clustered_results_ashr_TR_vs_TNR[[5]] = pathfindR_outputs_ashr_TR_vs_TNR[["GO-MF"]]
clustered_results_ashr_TR_vs_TNR[[5]]$Cluster = 1
clustered_results_ashr_TR_vs_TNR[[5]]$Status = "Member"
clustered_results_ashr_TR_vs_TNR[[6]] = pathfindR_outputs_ashr_TR_vs_TNR[["KEGG"]]
clustered_results_ashr_TR_vs_TNR[[6]]$Cluster = 1
clustered_results_ashr_TR_vs_TNR[[6]]$Status = "Member"
names(clustered_results_ashr_TR_vs_TNR) = cluster_names

# Reactome     : The maximum average silhouette width was 0.96 for k = 5
# BioCarta     : The maximum average silhouette width was NA
# GO-BP        : The maximum average silhouette width was 0.8 for k = 5   
# GO-CC        : The maximum average silhouette width was NA
# GO-MF        : The maximum average silhouette width was NA
# KEGG         : The maximum average silhouette width was NA

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_ashr_TR_vs_TNR = pathfindR_outputs_ashr_TR_vs_TNR
for (i in c(1, 3:6)){
  wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_ashr_TR_vs_TNR[[i]]$Term_Description, 
                                                                                     width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR = clustered_results_ashr_TR_vs_TNR
for (i in c(1, 3:6)){
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
for (i in c(1, 3:6)){
  # unclustered results
  enrichment_dotplots_ashr_TR_vs_TNR[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[i]],
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
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i],
                        " terms enrichment dotplot - (", comparisons["ashr_TR_vs_TNR"], ")"))
  print(enrichment_dotplots_ashr_TR_vs_TNR[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_ashr_TR_vs_TNR)[i], "_top10_dotplot.tiff"),
         path = paste0("Results/Pathways/ashr_TR_vs_TNR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  cluster_enrichment_dotplots_ashr_TR_vs_TNR[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_TR_vs_TNR[[names(pathfindR_outputs_ashr_TR_vs_TNR)[i]]][clustered_results_ashr_TR_vs_TNR[[names(pathfindR_outputs_ashr_TR_vs_TNR)[i]]]$Status
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
    labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i],
                        " terms enrichment dotplot - (", comparisons["ashr_TR_vs_TNR"], ")"))
  print(cluster_enrichment_dotplots_ashr_TR_vs_TNR[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_ashr_TR_vs_TNR)[i], "_top10_dotplot_clustered.tiff"),
         path = paste0("Results/Pathways/ashr_TR_vs_TNR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

names(enrichment_dotplots_ashr_TR_vs_TNR) = names(pathfindR_outputs_ashr_TR_vs_TNR)
names(cluster_enrichment_dotplots_ashr_TR_vs_TNR) = names(pathfindR_outputs_ashr_TR_vs_TNR)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb, subdirs[i])
  writeData(wb, subdirs[i], clustered_results_ashr_TR_vs_TNR[[i]])
}
saveWorkbook(wb, file = "Results/Pathways/ashr_TR_vs_TNR/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2_names = subdirs
wb2 = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb2, paste0(wb2_names[i], " - rep"))
  writeData(wb2, paste0(wb2_names[i], " - rep"), 
            clustered_results_ashr_TR_vs_TNR[[i]][clustered_results_ashr_TR_vs_TNR[[i]]$Status
                                                  == "Representative", ])
}
saveWorkbook(wb2, file = "Results/Pathways/ashr_TR_vs_TNR/pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####

# Defining a legend alignment function 
align_legend <- function(p, hjust = 0.5) {
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  if (length(legend_index) == 0) {
    warning("No legend found for this plot.")
    return(g)
  }
  legend <- grobs[[legend_index]]
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  if (length(guides_index) == 0) {
    warning("No guides layout found in legend.")
    return(g)
  }
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    if (length(guides$width) < 5) {
      warning("Insufficient guide width elements for alignment.")
      return(g)
    }
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust * spacing, 1)
    guides$widths[6] <- (1 - hjust) * spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    legend$grobs[[gi]] <- guides
  }
  g$grobs[[legend_index]] <- legend
  g
}

term_gene_heatmaps_ashr_TR_vs_TNR = list()
term_gene_graphs_ashr_TR_vs_TNR = list()

for (i in c(1, 3:6)){
  # term-gene heatmaps
  term_gene_heatmaps_ashr_TR_vs_TNR[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[i]],
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
          legend.title = element_text(size = 5, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i], 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_TR_vs_TNR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  print(ggdraw(align_legend(term_gene_heatmaps_ashr_TR_vs_TNR[[i]], hjust = 0.5)))
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i],
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Results/Pathways/ashr_TR_vs_TNR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_TR_vs_TNR[[i]] = term_gene_graph(result_df = pathfindR_outputs_ashr_TR_vs_TNR[[i]],
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
          legend.title = element_text(size = 5, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i], 
                        " term - gene graph (",
                        comparisons["ashr_TR_vs_TNR"], ")"))
  term_gene_graphs_ashr_TR_vs_TNR[[i]]$layers = list(term_gene_graphs_ashr_TR_vs_TNR[[i]]$layers[[1]], 
                                                     term_gene_graphs_ashr_TR_vs_TNR[[i]]$layers[[2]], 
                                                     term_gene_graphs_ashr_TR_vs_TNR[[i]]$layers[[4]])
  print(term_gene_graphs_ashr_TR_vs_TNR[[i]])
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i],
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Results/Pathways/ashr_TR_vs_TNR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i], "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_TR_vs_TNR = list()
for (i in c(1, 3:6)){
  # UpSet plot
  UpSet_plots_ashr_TR_vs_TNR[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_TR_vs_TNR[[i]],
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
         title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i], 
                        " term - UpSet plot (",
                        comparisons["ashr_TR_vs_TNR"], ")"))
  tiff(paste0("Results/Pathways/ashr_TR_vs_TNR/pathfindR/", names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i], "/",
              names(wrapped_pathfindR_outputs_ashr_TR_vs_TNR)[i], "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_ashr_TR_vs_TNR[[i]])
  dev.off()
}

names(term_gene_graphs_ashr_TR_vs_TNR) = names(pathfindR_outputs_ashr_TR_vs_TNR)
names(term_gene_heatmaps_ashr_TR_vs_TNR) = names(pathfindR_outputs_ashr_TR_vs_TNR)
names(UpSet_plots_ashr_TR_vs_TNR) = names(pathfindR_outputs_ashr_TR_vs_TNR)

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
                                                         p_val_threshold = 0.05, convert2alias = FALSE,
                                                         output_dir = dirs_ashr_LTR_vs_ETR[i], min_gset_size = 10,
                                                         max_gset_size = 300, adj_method = 'fdr',
                                                         enrichment_threshold = 0.05,
                                                         pin_name_path = 'Biogrid', search_method = 'GR',
                                                         iterations = 10, n_processes = 10)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_ashr_LTR_vs_ETR) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

cluster_names = subdirs
RNGversion("4.2.2")
set.seed(123)
clustered_results_ashr_LTR_vs_ETR = list()
for (i in 1:length(subdirs)){
  clustered_results_ashr_LTR_vs_ETR[[i]] = cluster_enriched_terms(pathfindR_outputs_ashr_LTR_vs_ETR[[cluster_names[[i]]]],
                                                                  method = "hierarchical",
                                                                  use_description = TRUE,
                                                                  plot_clusters_graph = FALSE)
}
names(clustered_results_ashr_LTR_vs_ETR) = cluster_names

# Reactome     : The maximum average silhouette width was 0.59 for k = 200
# BioCarta     : The maximum average silhouette width was 0.46 for k = 40
# GO-BP        : The maximum average silhouette width was 0.16 for k = 18  
# GO-CC        : The maximum average silhouette width was 0.22 for k = 20
# GO-MF        : The maximum average silhouette width was 0.19 for k = 30
# KEGG         : The maximum average silhouette width was 0.20 for k = 30

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_ashr_LTR_vs_ETR = pathfindR_outputs_ashr_LTR_vs_ETR
for (i in 1:length(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)){
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
for (i in 1:length(pathfindR_outputs_ashr_LTR_vs_ETR)){
  # unclustered results
  enrichment_dotplots_ashr_LTR_vs_ETR[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[i]],
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
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i],
                        " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_ETR"], ")"))
  print(enrichment_dotplots_ashr_LTR_vs_ETR[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_ashr_LTR_vs_ETR)[i], "_top10_dotplot.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_ETR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  cluster_enrichment_dotplots_ashr_LTR_vs_ETR[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_ETR[[names(pathfindR_outputs_ashr_LTR_vs_ETR)[i]]][clustered_results_ashr_LTR_vs_ETR[[names(pathfindR_outputs_ashr_LTR_vs_ETR)[i]]]$Status
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
    labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i],
                        " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_ETR"], ")"))
  print(cluster_enrichment_dotplots_ashr_LTR_vs_ETR[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_ashr_LTR_vs_ETR)[i], "_top10_dotplot_clustered.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_ETR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

names(enrichment_dotplots_ashr_LTR_vs_ETR) = names(pathfindR_outputs_ashr_LTR_vs_ETR)
names(cluster_enrichment_dotplots_ashr_LTR_vs_ETR) = names(pathfindR_outputs_ashr_LTR_vs_ETR)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb, subdirs[i])
  writeData(wb, subdirs[i], clustered_results_ashr_LTR_vs_ETR[[i]])
}
saveWorkbook(wb, file = "Results/Pathways/ashr_LTR_vs_ETR/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2_names = subdirs
wb2 = createWorkbook()
for(i in 1:length(subdirs)) {
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

for (i in 1:length(pathfindR_outputs_ashr_LTR_vs_ETR)){
  # term-gene heatmaps
  term_gene_heatmaps_ashr_LTR_vs_ETR[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[i]],
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
          legend.title = element_text(size = 5, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i], 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_LTR_vs_ETR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  print(ggdraw(align_legend(term_gene_heatmaps_ashr_LTR_vs_ETR[[i]], hjust = 0.5)))
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i],
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_ETR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_LTR_vs_ETR[[i]] = term_gene_graph(result_df = pathfindR_outputs_ashr_LTR_vs_ETR[[i]],
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
          legend.title = element_text(size = 5, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i], 
                        " term - gene graph (",
                        comparisons["ashr_LTR_vs_ETR"], ")"))
  term_gene_graphs_ashr_LTR_vs_ETR[[i]]$layers = list(term_gene_graphs_ashr_LTR_vs_ETR[[i]]$layers[[1]], 
                                                      term_gene_graphs_ashr_LTR_vs_ETR[[i]]$layers[[2]], 
                                                      term_gene_graphs_ashr_LTR_vs_ETR[[i]]$layers[[4]])
  print(term_gene_graphs_ashr_LTR_vs_ETR[[i]])
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i],
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_ETR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i], "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_LTR_vs_ETR = list()
for (i in 1:length(pathfindR_outputs_ashr_LTR_vs_ETR)){
  # UpSet plot
  UpSet_plots_ashr_LTR_vs_ETR[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_ETR[[i]],
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
         title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i], 
                        " term - UpSet plot (",
                        comparisons["ashr_LTR_vs_ETR"], ")"))
  tiff(paste0("Results/Pathways/ashr_LTR_vs_ETR/pathfindR/", names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i], "/",
              names(wrapped_pathfindR_outputs_ashr_LTR_vs_ETR)[i], "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_ashr_LTR_vs_ETR[[i]])
  dev.off()
}

names(term_gene_graphs_ashr_LTR_vs_ETR) = names(pathfindR_outputs_ashr_LTR_vs_ETR)
names(term_gene_heatmaps_ashr_LTR_vs_ETR) = names(pathfindR_outputs_ashr_LTR_vs_ETR)
names(UpSet_plots_ashr_LTR_vs_ETR) = names(pathfindR_outputs_ashr_LTR_vs_ETR)

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
                                                          p_val_threshold = 0.05, convert2alias = FALSE,
                                                          output_dir = dirs_ashr_LTR_vs_LTNR[i], min_gset_size = 10,
                                                          max_gset_size = 300, adj_method = 'fdr',
                                                          enrichment_threshold = 0.05,
                                                          pin_name_path = 'Biogrid', search_method = 'GR',
                                                          iterations = 10, n_processes = 10)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_ashr_LTR_vs_LTNR) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

cluster_names = subdirs
RNGversion("4.2.2")
set.seed(123)
clustered_results_ashr_LTR_vs_LTNR = list()
for (i in 1:length(subdirs)){
  clustered_results_ashr_LTR_vs_LTNR[[i]] = cluster_enriched_terms(pathfindR_outputs_ashr_LTR_vs_LTNR[[cluster_names[[i]]]],
                                                                   method = "hierarchical",
                                                                   use_description = TRUE,
                                                                   plot_clusters_graph = FALSE)
}
names(clustered_results_ashr_LTR_vs_LTNR) = cluster_names

# Reactome     : The maximum average silhouette width was 0.66 for k = 200
# BioCarta     : The maximum average silhouette width was 0.61 for k = 14
# GO-BP        : The maximum average silhouette width was 0.30 for k = 30  
# GO-CC        : The maximum average silhouette width was 0.24 for k = 16
# GO-MF        : The maximum average silhouette width was 0.51 for k = 2
# KEGG         : The maximum average silhouette width was 0.28 for k = 16

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR = pathfindR_outputs_ashr_LTR_vs_LTNR
for (i in 1:length(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)){
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
for (i in 1:length(pathfindR_outputs_ashr_LTR_vs_LTNR)){
  # unclustered results
  enrichment_dotplots_ashr_LTR_vs_LTNR[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[i]],
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
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i],
                        " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_LTNR"], ")"))
  print(enrichment_dotplots_ashr_LTR_vs_LTNR[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i], "_top10_dotplot.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  cluster_enrichment_dotplots_ashr_LTR_vs_LTNR[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_LTR_vs_LTNR[[names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i]]][clustered_results_ashr_LTR_vs_LTNR[[names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i]]]$Status
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
    labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i],
                        " terms enrichment dotplot - (", comparisons["ashr_LTR_vs_LTNR"], ")"))
  print(cluster_enrichment_dotplots_ashr_LTR_vs_LTNR[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_ashr_LTR_vs_LTNR)[i], "_top10_dotplot_clustered.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

names(enrichment_dotplots_ashr_LTR_vs_LTNR) = names(pathfindR_outputs_ashr_LTR_vs_LTNR)
names(cluster_enrichment_dotplots_ashr_LTR_vs_LTNR) = names(pathfindR_outputs_ashr_LTR_vs_LTNR)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb, subdirs[i])
  writeData(wb, subdirs[i], clustered_results_ashr_LTR_vs_LTNR[[i]])
}
saveWorkbook(wb, file = "Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2_names = subdirs
wb2 = createWorkbook()
for(i in 1:length(subdirs)) {
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

for (i in 1:length(pathfindR_outputs_ashr_LTR_vs_LTNR)){
  # term-gene heatmaps
  term_gene_heatmaps_ashr_LTR_vs_LTNR[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[i]],
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
          legend.title = element_text(size = 5, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i], 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_LTR_vs_LTNR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  print(ggdraw(align_legend(term_gene_heatmaps_ashr_LTR_vs_LTNR[[i]], hjust = 0.5)))
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i],
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_LTR_vs_LTNR[[i]] = term_gene_graph(result_df = pathfindR_outputs_ashr_LTR_vs_LTNR[[i]],
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
          legend.title = element_text(size = 5, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i], 
                        " term - gene graph (",
                        comparisons["ashr_LTR_vs_LTNR"], ")"))
  term_gene_graphs_ashr_LTR_vs_LTNR[[i]]$layers = list(term_gene_graphs_ashr_LTR_vs_LTNR[[i]]$layers[[1]], 
                                                       term_gene_graphs_ashr_LTR_vs_LTNR[[i]]$layers[[2]], 
                                                       term_gene_graphs_ashr_LTR_vs_LTNR[[i]]$layers[[4]])
  print(term_gene_graphs_ashr_LTR_vs_LTNR[[i]])
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i],
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i], "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_LTR_vs_LTNR = list()
for (i in 1:length(pathfindR_outputs_ashr_LTR_vs_LTNR)){
  # UpSet plot
  UpSet_plots_ashr_LTR_vs_LTNR[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR[[i]],
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
         title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i], 
                        " term - UpSet plot (",
                        comparisons["ashr_LTR_vs_LTNR"], ")"))
  tiff(paste0("Results/Pathways/ashr_LTR_vs_LTNR/pathfindR/", names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i], "/",
              names(wrapped_pathfindR_outputs_ashr_LTR_vs_LTNR)[i], "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_ashr_LTR_vs_LTNR[[i]])
  dev.off()
}

names(term_gene_graphs_ashr_LTR_vs_LTNR) = names(pathfindR_outputs_ashr_LTR_vs_LTNR)
names(term_gene_heatmaps_ashr_LTR_vs_LTNR) = names(pathfindR_outputs_ashr_LTR_vs_LTNR)
names(UpSet_plots_ashr_LTR_vs_LTNR) = names(pathfindR_outputs_ashr_LTR_vs_LTNR)

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
                                                          p_val_threshold = 0.05, convert2alias = FALSE,
                                                          output_dir = dirs_ashr_ETR_vs_ETNR[i], min_gset_size = 10,
                                                          max_gset_size = 300, adj_method = 'fdr',
                                                          enrichment_threshold = 0.05,
                                                          pin_name_path = 'Biogrid', search_method = 'GR',
                                                          iterations = 10, n_processes = 10)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_ashr_ETR_vs_ETNR) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

cluster_names = subdirs
RNGversion("4.2.2")
set.seed(123)
clustered_results_ashr_ETR_vs_ETNR = list()
for (i in 1:3){
  clustered_results_ashr_ETR_vs_ETNR[[i]] = cluster_enriched_terms(pathfindR_outputs_ashr_ETR_vs_ETNR[[cluster_names[[i]]]],
                                                                   method = "hierarchical",
                                                                   use_description = TRUE,
                                                                   plot_clusters_graph = FALSE)
}
clustered_results_ashr_ETR_vs_ETNR[[4]] = pathfindR_outputs_ashr_ETR_vs_ETNR[["GO-CC"]]
clustered_results_ashr_ETR_vs_ETNR[[4]]$Cluster = 1
clustered_results_ashr_ETR_vs_ETNR[[4]]$Status = "Member"
clustered_results_ashr_ETR_vs_ETNR[[5]] = pathfindR_outputs_ashr_ETR_vs_ETNR[["GO-MF"]]
clustered_results_ashr_ETR_vs_ETNR[[5]]$Cluster = 1
clustered_results_ashr_ETR_vs_ETNR[[5]]$Status = "Member"
clustered_results_ashr_ETR_vs_ETNR[[6]] = pathfindR_outputs_ashr_ETR_vs_ETNR[["KEGG"]]
clustered_results_ashr_ETR_vs_ETNR[[6]]$Cluster = 1
clustered_results_ashr_ETR_vs_ETNR[[6]]$Status = "Member"
names(clustered_results_ashr_ETR_vs_ETNR) = cluster_names

# Reactome     : The maximum average silhouette width was 0.89 for k = 40
# BioCarta     : The maximum average silhouette width was 0.82 for k = 14
# GO-BP        : The maximum average silhouette width was 0.29 for k = 4  
# GO-CC        : The maximum average silhouette width was NA
# GO-MF        : The maximum average silhouette width was NA
# KEGG         : The maximum average silhouette width was NA

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR = pathfindR_outputs_ashr_ETR_vs_ETNR
for (i in 1:length(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)){
  wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_ashr_ETR_vs_ETNR[[i]]$Term_Description, 
                                                                                       width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR = clustered_results_ashr_ETR_vs_ETNR
for (i in 1:3){
  wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[i]]$Term_Description = stringr::str_wrap(clustered_results_ashr_ETR_vs_ETNR[[i]]$Term_Description, 
                                                                                                 width = 41)
}
rm(i)

enrichment_dotplots_ashr_ETR_vs_ETNR = list()
cluster_enrichment_dotplots_ashr_ETR_vs_ETNR = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_ashr_ETR_vs_ETNR)){
  # unclustered results
  enrichment_dotplots_ashr_ETR_vs_ETNR[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[i]],
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
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i],
                        " terms enrichment dotplot - (", comparisons["ashr_ETR_vs_ETNR"], ")"))
  print(enrichment_dotplots_ashr_ETR_vs_ETNR[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i], "_top10_dotplot.tiff"),
         path = paste0("Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  if (i %in% 1:3) {
    cluster_enrichment_dotplots_ashr_ETR_vs_ETNR[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_ashr_ETR_vs_ETNR[[names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i]]][clustered_results_ashr_ETR_vs_ETNR[[names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i]]]$Status
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
      labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i],
                          " terms enrichment dotplot - (", comparisons["ashr_ETR_vs_ETNR"], ")"))
    print(cluster_enrichment_dotplots_ashr_ETR_vs_ETNR[[i]])
    ggsave(filename = paste0(names(pathfindR_outputs_ashr_ETR_vs_ETNR)[i], "_top10_dotplot_clustered.tiff"),
           path = paste0("Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/", 
                         names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i], "/"),
           width = 2880, height = 1620, device = 'tiff', units = "px",
           dpi = 700, compression = "lzw")
    dev.off()
  }
}

names(enrichment_dotplots_ashr_ETR_vs_ETNR) = names(pathfindR_outputs_ashr_ETR_vs_ETNR)
names(cluster_enrichment_dotplots_ashr_ETR_vs_ETNR) = names(pathfindR_outputs_ashr_ETR_vs_ETNR)[1:3]

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb, subdirs[i])
  writeData(wb, subdirs[i], clustered_results_ashr_ETR_vs_ETNR[[i]])
}
saveWorkbook(wb, file = "Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2_names = subdirs
wb2 = createWorkbook()
for(i in 1:length(subdirs)) {
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

for (i in 1:length(pathfindR_outputs_ashr_ETR_vs_ETNR)){
  # term-gene heatmaps
  term_gene_heatmaps_ashr_ETR_vs_ETNR[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[i]],
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
          legend.title = element_text(size = 5, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i], 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["ashr_ETR_vs_ETNR"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  print(ggdraw(align_legend(term_gene_heatmaps_ashr_ETR_vs_ETNR[[i]], hjust = 0.5)))
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i],
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_ashr_ETR_vs_ETNR[[i]] = term_gene_graph(result_df = pathfindR_outputs_ashr_ETR_vs_ETNR[[i]],
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
          legend.title = element_text(size = 5, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i], 
                        " term - gene graph (",
                        comparisons["ashr_ETR_vs_ETNR"], ")"))
  term_gene_graphs_ashr_ETR_vs_ETNR[[i]]$layers = list(term_gene_graphs_ashr_ETR_vs_ETNR[[i]]$layers[[1]], 
                                                       term_gene_graphs_ashr_ETR_vs_ETNR[[i]]$layers[[2]], 
                                                       term_gene_graphs_ashr_ETR_vs_ETNR[[i]]$layers[[4]])
  print(term_gene_graphs_ashr_ETR_vs_ETNR[[i]])
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i],
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/", 
                       names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i], "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_ashr_ETR_vs_ETNR = list()
for (i in 1:length(pathfindR_outputs_ashr_ETR_vs_ETNR)){
  # UpSet plot
  UpSet_plots_ashr_ETR_vs_ETNR[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR[[i]],
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
         title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i], 
                        " term - UpSet plot (",
                        comparisons["ashr_ETR_vs_ETNR"], ")"))
  tiff(paste0("Results/Pathways/ashr_ETR_vs_ETNR/pathfindR/", names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i], "/",
              names(wrapped_pathfindR_outputs_ashr_ETR_vs_ETNR)[i], "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_ashr_ETR_vs_ETNR[[i]])
  dev.off()
}

names(term_gene_graphs_ashr_ETR_vs_ETNR) = names(pathfindR_outputs_ashr_ETR_vs_ETNR)
names(term_gene_heatmaps_ashr_ETR_vs_ETNR) = names(pathfindR_outputs_ashr_ETR_vs_ETNR)
names(UpSet_plots_ashr_ETR_vs_ETNR) = names(pathfindR_outputs_ashr_ETR_vs_ETNR)

# Export a session info text file for reproducibility
writeLines(capture.output(sessionInfo()), 
           paste0("Resources_and_environments/Pathways_sessionInfo.txt"))

# Save environment for easy loading for DGEA later
save.image("Resources_and_environments/Pathways.RData")