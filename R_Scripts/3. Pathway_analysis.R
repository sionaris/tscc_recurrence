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
    if (!dir.exists(paste0("Pathways/", names(ashr_dgea)[j], "/",
                           pathway_approaches[[i]]))) {
      dir.create(paste0("Pathways/", names(ashr_dgea)[j], "/",
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

# HMBA24h_vs_Control48h #####
# Loading the input to pathfindR (the stage 1 vs normal topTable output):
pathf_input_HMBA24h_vs_Control48h = ashr_dgea$HMBA24h_vs_Control48h %>%
  dplyr::select(Gene.Symbol, log2FoldChange, padj) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_HMBA24h_vs_Control48h = paste0("Pathways/HMBA24h_vs_Control48h/pathfindR/", subdirs)
pathfindR_outputs_HMBA24h_vs_Control48h = list()

RNGversion("4.2.2")
set.seed(123)
for (i in 1:length(dirs_HMBA24h_vs_Control48h)){
  pathfindR_outputs_HMBA24h_vs_Control48h[[i]] = run_pathfindR(pathf_input_HMBA24h_vs_Control48h, gene_sets = "Custom",
                                                               p_val_threshold = 0.05, convert2alias = FALSE,
                                                               custom_genes = genes_of_gene_sets[[i]],
                                                               custom_descriptions = descriptions_of_gene_sets[[i]],
                                                               output_dir = dirs_HMBA24h_vs_Control48h[i], min_gset_size = 10,
                                                               max_gset_size = 300, adj_method = 'fdr',
                                                               enrichment_threshold = 0.05,
                                                               pin_name_path = path2SIF, search_method = 'GR',
                                                               grMaxDepth = 1, grSearchDepth = 1,
                                                               iterations = 10, n_processes = 10)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_HMBA24h_vs_Control48h) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

cluster_names = subdirs
RNGversion("4.2.2")
set.seed(123)
clustered_results_HMBA24h_vs_Control48h = list()
for (i in 1:length(cluster_names)){
  clustered_results_HMBA24h_vs_Control48h[[i]] = cluster_enriched_terms(pathfindR_outputs_HMBA24h_vs_Control48h[[cluster_names[[i]]]],
                                                                        method = "hierarchical")
}
names(clustered_results_HMBA24h_vs_Control48h) = cluster_names

# MSigDB Reactome     : The maximum average silhouette width was 0.35 for k = 100
# MSigDB WikiPathways : The maximum average silhouette width was 0.12 for k = 90
# MSigDB GO-BP        : The maximum average silhouette width was 0.4 for k = 400  
# MSigDB GO-CC        : The maximum average silhouette width was 0.29 for k = 20
# MSigDB GO-MF        : The maximum average silhouette width was 0.25 for k = 30
# MSigDB BioCarta     : The maximum average silhouette width was 0.17 for k = 14 
# KEGG                : The maximum average silhouette width was 0.15 for k = 2

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_HMBA24h_vs_Control48h = pathfindR_outputs_HMBA24h_vs_Control48h
for (i in 1:length(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)){
  wrapped_pathfindR_outputs_HMBA24h_vs_Control48h[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_HMBA24h_vs_Control48h[[i]]$Term_Description, 
                                                                                            width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_HMBA24h_vs_Control48h = clustered_results_HMBA24h_vs_Control48h
for (i in 1:length(wrapped_clustered_pathfindR_outputs_HMBA24h_vs_Control48h)){
  wrapped_clustered_pathfindR_outputs_HMBA24h_vs_Control48h[[i]]$Term_Description = stringr::str_wrap(clustered_results_HMBA24h_vs_Control48h[[i]]$Term_Description, 
                                                                                                      width = 41)
}
rm(i)

comparisons = c("HMBA24h vs. Control48h", "HMBA48h vs. Control48h", 
                "HMBA48h vs. HMBA24h", "HMBA72h vs. Control48h", 
                "HMBA72h vs. HMBA24h", "HMBA72h vs. HMBA48h")
names(comparisons) = c("HMBA24h_vs_Control48h", "HMBA48h_vs_Control48h", 
                       "HMBA48h_vs_HMBA24h", "HMBA72h_vs_Control48h", 
                       "HMBA72h_vs_HMBA24h", "HMBA72h_vs_HMBA48h")

enrichment_dotplots_HMBA24h_vs_Control48h = list()
cluster_enrichment_dotplots_HMBA24h_vs_Control48h = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_HMBA24h_vs_Control48h)){
  # unclustered results
  enrichment_dotplots_HMBA24h_vs_Control48h[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_HMBA24h_vs_Control48h[[i]],
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
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA24h_vs_Control48h"], ")"))
  print(enrichment_dotplots_HMBA24h_vs_Control48h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA24h_vs_Control48h)[i], "_top10_dotplot.tiff"),
         path = paste0("Pathways/HMBA24h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  cluster_enrichment_dotplots_HMBA24h_vs_Control48h[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_HMBA24h_vs_Control48h[[names(pathfindR_outputs_HMBA24h_vs_Control48h)[i]]][clustered_results_HMBA24h_vs_Control48h[[names(pathfindR_outputs_HMBA24h_vs_Control48h)[i]]]$Status
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
    labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA24h_vs_Control48h"], ")"))
  print(cluster_enrichment_dotplots_HMBA24h_vs_Control48h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA24h_vs_Control48h)[i], "_top10_dotplot_clustered.tiff"),
         path = paste0("Pathways/HMBA24h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

names(enrichment_dotplots_HMBA24h_vs_Control48h) = names(pathfindR_outputs_HMBA24h_vs_Control48h)
names(cluster_enrichment_dotplots_HMBA24h_vs_Control48h) = names(pathfindR_outputs_HMBA24h_vs_Control48h)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb, subdirs[i])
  writeData(wb, subdirs[i], clustered_results_HMBA24h_vs_Control48h[[i]])
}
saveWorkbook(wb, file = "Pathways/HMBA24h_vs_Control48h/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2_names = c("Reactome", "WikiPathways", "GO-Biological Processes", 
              "GO-Cellular Components", "GO-Molecular Functions",
              "BioCarta", "KEGG")
wb2 = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb2, paste0(wb2_names[i], " - rep"))
  writeData(wb2, paste0(wb2_names[i], " - rep"), 
            clustered_results_HMBA24h_vs_Control48h[[i]][clustered_results_HMBA24h_vs_Control48h[[i]]$Status
                                                         == "Representative", ])
}
saveWorkbook(wb2, file = "Pathways/HMBA24h_vs_Control48h/pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####

# Defining a legend alignment function 
align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
}

term_gene_heatmaps_HMBA24h_vs_Control48h = list()
term_gene_graphs_HMBA24h_vs_Control48h = list()

for (i in 1:length(pathfindR_outputs_HMBA24h_vs_Control48h)){
  # term-gene heatmaps
  term_gene_heatmaps_HMBA24h_vs_Control48h[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_HMBA24h_vs_Control48h[[i]],
                                                                    genes_df = pathf_input_HMBA24h_vs_Control48h,
                                                                    num_terms = 10,
                                                                    use_description = TRUE,
                                                                    low = "darkblue",
                                                                    high = "red",
                                                                    mid = "white",
                                                                    pin_name_path = path2SIF)+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 2.5,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i], 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["HMBA24h_vs_Control48h"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  print(ggdraw(align_legend(term_gene_heatmaps_HMBA24h_vs_Control48h[[i]], hjust = 0.5)))
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i],
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Pathways/HMBA24h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_HMBA24h_vs_Control48h[[i]] = term_gene_graph(result_df = pathfindR_outputs_HMBA24h_vs_Control48h[[i]],
                                                                num_terms = 5,
                                                                use_description = TRUE,
                                                                node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 1, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i], 
                        " term - gene graph (",
                        comparisons["HMBA24h_vs_Control48h"], ")"))
  term_gene_graphs_HMBA24h_vs_Control48h[[i]]$layers = list(term_gene_graphs_HMBA24h_vs_Control48h[[i]]$layers[[1]], 
                                                            term_gene_graphs_HMBA24h_vs_Control48h[[i]]$layers[[2]], 
                                                            term_gene_graphs_HMBA24h_vs_Control48h[[i]]$layers[[4]])
  print(term_gene_graphs_HMBA24h_vs_Control48h[[i]])
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i],
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Pathways/HMBA24h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i], "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_HMBA24h_vs_Control48h = list()
for (i in 1:length(pathfindR_outputs_HMBA24h_vs_Control48h)){
  # UpSet plot
  UpSet_plots_HMBA24h_vs_Control48h[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_HMBA24h_vs_Control48h[[i]],
                                                      genes_df = pathf_input_HMBA24h_vs_Control48h,
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
         title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i], 
                        " term - UpSet plot (",
                        comparisons["HMBA24h_vs_Control48h"], ")"))
  tiff(paste0("Pathways/HMBA24h_vs_Control48h/pathfindR/", names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i], "/",
              names(wrapped_pathfindR_outputs_HMBA24h_vs_Control48h)[i], "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_HMBA24h_vs_Control48h[[i]])
  dev.off()
}

names(term_gene_graphs_HMBA24h_vs_Control48h) = names(pathfindR_outputs_HMBA24h_vs_Control48h)
names(term_gene_heatmaps_HMBA24h_vs_Control48h) = names(pathfindR_outputs_HMBA24h_vs_Control48h)
names(UpSet_plots_HMBA24h_vs_Control48h) = names(pathfindR_outputs_HMBA24h_vs_Control48h)

# HMBA48h_vs_Control48h #####
# Loading the input to pathfindR (the stage 1 vs normal topTable output):
pathf_input_HMBA48h_vs_Control48h = ashr_dgea$HMBA48h_vs_Control48h %>%
  dplyr::select(Gene.Symbol, log2FoldChange, padj) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_HMBA48h_vs_Control48h = paste0("Pathways/HMBA48h_vs_Control48h/pathfindR/", subdirs)
pathfindR_outputs_HMBA48h_vs_Control48h = list()

RNGversion("4.2.2")
set.seed(123)
for (i in 1:length(dirs_HMBA48h_vs_Control48h)){
  pathfindR_outputs_HMBA48h_vs_Control48h[[i]] = run_pathfindR(pathf_input_HMBA48h_vs_Control48h, gene_sets = "Custom",
                                                               p_val_threshold = 0.05, convert2alias = FALSE,
                                                               custom_genes = genes_of_gene_sets[[i]],
                                                               custom_descriptions = descriptions_of_gene_sets[[i]],
                                                               output_dir = dirs_HMBA48h_vs_Control48h[i], min_gset_size = 10,
                                                               max_gset_size = 300, adj_method = 'fdr',
                                                               enrichment_threshold = 0.05,
                                                               pin_name_path = path2SIF, search_method = 'GR',
                                                               grMaxDepth = 1, grSearchDepth = 1,
                                                               iterations = 10, n_processes = 10)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_HMBA48h_vs_Control48h) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

cluster_names = subdirs
RNGversion("4.2.2")
set.seed(123)
clustered_results_HMBA48h_vs_Control48h = list()
for (i in 1:length(cluster_names)){
  clustered_results_HMBA48h_vs_Control48h[[i]] = cluster_enriched_terms(pathfindR_outputs_HMBA48h_vs_Control48h[[cluster_names[[i]]]],
                                                                        method = "hierarchical")
}
names(clustered_results_HMBA48h_vs_Control48h) = cluster_names

# MSigDB Reactome     : The maximum average silhouette width was 0.35 for k = 100
# MSigDB WikiPathways : The maximum average silhouette width was 0.11 for k = 100
# MSigDB GO-BP        : The maximum average silhouette width was 0.39 for k = 400 
# MSigDB GO-CC        : The maximum average silhouette width was 0.34 for k = 30
# MSigDB GO-MF        : The maximum average silhouette width was 0.34 for k = 50
# MSigDB BioCarta     : The maximum average silhouette width was 0.2 for k = 30
# KEGG                : The maximum average silhouette width was 0.14 for k = 2

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_HMBA48h_vs_Control48h = pathfindR_outputs_HMBA48h_vs_Control48h
for (i in 1:length(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)){
  wrapped_pathfindR_outputs_HMBA48h_vs_Control48h[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_HMBA48h_vs_Control48h[[i]]$Term_Description, 
                                                                                            width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_HMBA48h_vs_Control48h = clustered_results_HMBA48h_vs_Control48h
for (i in 1:length(wrapped_clustered_pathfindR_outputs_HMBA48h_vs_Control48h)){
  wrapped_clustered_pathfindR_outputs_HMBA48h_vs_Control48h[[i]]$Term_Description = stringr::str_wrap(clustered_results_HMBA48h_vs_Control48h[[i]]$Term_Description, 
                                                                                                      width = 41)
}
rm(i)

enrichment_dotplots_HMBA48h_vs_Control48h = list()
cluster_enrichment_dotplots_HMBA48h_vs_Control48h = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_HMBA48h_vs_Control48h)){
  # unclustered results
  enrichment_dotplots_HMBA48h_vs_Control48h[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_HMBA48h_vs_Control48h[[i]],
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
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA48h_vs_Control48h"], ")"))
  print(enrichment_dotplots_HMBA48h_vs_Control48h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA48h_vs_Control48h)[i], "_top10_dotplot.tiff"),
         path = paste0("Pathways/HMBA48h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  cluster_enrichment_dotplots_HMBA48h_vs_Control48h[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_HMBA48h_vs_Control48h[[names(pathfindR_outputs_HMBA48h_vs_Control48h)[i]]][clustered_results_HMBA48h_vs_Control48h[[names(pathfindR_outputs_HMBA48h_vs_Control48h)[i]]]$Status
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
    labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA48h_vs_Control48h"], ")"))
  print(cluster_enrichment_dotplots_HMBA48h_vs_Control48h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA48h_vs_Control48h)[i], "_top10_dotplot_clustered.tiff"),
         path = paste0("Pathways/HMBA48h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

names(enrichment_dotplots_HMBA48h_vs_Control48h) = names(pathfindR_outputs_HMBA48h_vs_Control48h)
names(cluster_enrichment_dotplots_HMBA48h_vs_Control48h) = names(pathfindR_outputs_HMBA48h_vs_Control48h)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb, subdirs[i])
  writeData(wb, subdirs[i], clustered_results_HMBA48h_vs_Control48h[[i]])
}
saveWorkbook(wb, file = "Pathways/HMBA48h_vs_Control48h/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2 = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb2, paste0(wb2_names[i], " - rep"))
  writeData(wb2, paste0(wb2_names[i], " - rep"), 
            clustered_results_HMBA48h_vs_Control48h[[i]][clustered_results_HMBA48h_vs_Control48h[[i]]$Status
                                                         == "Representative", ])
}
saveWorkbook(wb2, file = "Pathways/HMBA48h_vs_Control48h/pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_HMBA48h_vs_Control48h = list()
term_gene_graphs_HMBA48h_vs_Control48h = list()

for (i in 1:length(pathfindR_outputs_HMBA48h_vs_Control48h)){
  # term-gene heatmaps
  term_gene_heatmaps_HMBA48h_vs_Control48h[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_HMBA48h_vs_Control48h[[i]],
                                                                    genes_df = pathf_input_HMBA48h_vs_Control48h,
                                                                    num_terms = 10,
                                                                    use_description = TRUE,
                                                                    low = "darkblue",
                                                                    high = "red",
                                                                    mid = "white",
                                                                    pin_name_path = path2SIF)+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 2.5,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i], 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["HMBA48h_vs_Control48h"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  print(ggdraw(align_legend(term_gene_heatmaps_HMBA48h_vs_Control48h[[i]], hjust = 0.5)))
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i],
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Pathways/HMBA48h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_HMBA48h_vs_Control48h[[i]] = term_gene_graph(result_df = pathfindR_outputs_HMBA48h_vs_Control48h[[i]],
                                                                num_terms = 5,
                                                                use_description = TRUE,
                                                                node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 1, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i], 
                        " term - gene graph (",
                        comparisons["HMBA48h_vs_Control48h"], ")"))
  term_gene_graphs_HMBA48h_vs_Control48h[[i]]$layers = list(term_gene_graphs_HMBA48h_vs_Control48h[[i]]$layers[[1]], 
                                                            term_gene_graphs_HMBA48h_vs_Control48h[[i]]$layers[[2]], 
                                                            term_gene_graphs_HMBA48h_vs_Control48h[[i]]$layers[[4]])
  print(term_gene_graphs_HMBA48h_vs_Control48h[[i]])
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i],
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Pathways/HMBA48h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i], "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_HMBA48h_vs_Control48h = list()
for (i in 1:length(pathfindR_outputs_HMBA48h_vs_Control48h)){
  # UpSet plot
  UpSet_plots_HMBA48h_vs_Control48h[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_HMBA48h_vs_Control48h[[i]],
                                                      genes_df = pathf_input_HMBA48h_vs_Control48h,
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
         title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i], 
                        " term - UpSet plot (",
                        comparisons["HMBA48h_vs_Control48h"], ")"))
  tiff(paste0("Pathways/HMBA48h_vs_Control48h/pathfindR/", names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i], "/",
              names(wrapped_pathfindR_outputs_HMBA48h_vs_Control48h)[i], "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_HMBA48h_vs_Control48h[[i]])
  dev.off()
}

names(term_gene_graphs_HMBA48h_vs_Control48h) = names(pathfindR_outputs_HMBA48h_vs_Control48h)
names(term_gene_heatmaps_HMBA48h_vs_Control48h) = names(pathfindR_outputs_HMBA48h_vs_Control48h)
names(UpSet_plots_HMBA48h_vs_Control48h) = names(pathfindR_outputs_HMBA48h_vs_Control48h)

# HMBA48h_vs_HMBA24h #####
# Loading the input to pathfindR (the stage 1 vs normal topTable output):
pathf_input_HMBA48h_vs_HMBA24h = ashr_dgea$HMBA48h_vs_HMBA24h %>%
  dplyr::select(Gene.Symbol, log2FoldChange, padj) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_HMBA48h_vs_HMBA24h = paste0("Pathways/HMBA48h_vs_HMBA24h/pathfindR/", subdirs)
pathfindR_outputs_HMBA48h_vs_HMBA24h = list()

RNGversion("4.2.2")
set.seed(123)
for (i in 1:length(dirs_HMBA48h_vs_HMBA24h)){
  pathfindR_outputs_HMBA48h_vs_HMBA24h[[i]] = run_pathfindR(pathf_input_HMBA48h_vs_HMBA24h, gene_sets = "Custom",
                                                            p_val_threshold = 0.05, convert2alias = FALSE,
                                                            custom_genes = genes_of_gene_sets[[i]],
                                                            custom_descriptions = descriptions_of_gene_sets[[i]],
                                                            output_dir = dirs_HMBA48h_vs_HMBA24h[i], min_gset_size = 10,
                                                            max_gset_size = 300, adj_method = 'fdr',
                                                            enrichment_threshold = 0.05,
                                                            pin_name_path = path2SIF, search_method = 'GR',
                                                            grMaxDepth = 1, grSearchDepth = 1,
                                                            iterations = 10, n_processes = 10)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_HMBA48h_vs_HMBA24h) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

cluster_names = subdirs
RNGversion("4.2.2")
set.seed(123)
clustered_results_HMBA48h_vs_HMBA24h = list()
for (i in 1:length(cluster_names)){
  clustered_results_HMBA48h_vs_HMBA24h[[i]] = cluster_enriched_terms(pathfindR_outputs_HMBA48h_vs_HMBA24h[[cluster_names[[i]]]],
                                                                     method = "hierarchical")
}
names(clustered_results_HMBA48h_vs_HMBA24h) = cluster_names

# MSigDB Reactome     : The maximum average silhouette width was 0.38 for k = 150
# MSigDB WikiPathways : The maximum average silhouette width was 0.17 for k = 100
# MSigDB GO-BP        : The maximum average silhouette width was 0.41 for k = 450 
# MSigDB GO-CC        : The maximum average silhouette width was 0.42 for k = 40
# MSigDB GO-MF        : The maximum average silhouette width was 0.32 for k = 50
# MSigDB BioCarta     : The maximum average silhouette width was 0.32 for k = 30
# KEGG                : The maximum average silhouette width was 0.13 for k = 90

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h = pathfindR_outputs_HMBA48h_vs_HMBA24h
for (i in 1:length(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)){
  wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_HMBA48h_vs_HMBA24h[[i]]$Term_Description, 
                                                                                         width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_HMBA48h_vs_HMBA24h = clustered_results_HMBA48h_vs_HMBA24h
for (i in 1:length(wrapped_clustered_pathfindR_outputs_HMBA48h_vs_HMBA24h)){
  wrapped_clustered_pathfindR_outputs_HMBA48h_vs_HMBA24h[[i]]$Term_Description = stringr::str_wrap(clustered_results_HMBA48h_vs_HMBA24h[[i]]$Term_Description, 
                                                                                                   width = 41)
}
rm(i)

enrichment_dotplots_HMBA48h_vs_HMBA24h = list()
cluster_enrichment_dotplots_HMBA48h_vs_HMBA24h = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_HMBA48h_vs_HMBA24h)){
  # unclustered results
  enrichment_dotplots_HMBA48h_vs_HMBA24h[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h[[i]],
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
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA48h_vs_HMBA24h"], ")"))
  print(enrichment_dotplots_HMBA48h_vs_HMBA24h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA48h_vs_HMBA24h)[i], "_top10_dotplot.tiff"),
         path = paste0("Pathways/HMBA48h_vs_HMBA24h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  cluster_enrichment_dotplots_HMBA48h_vs_HMBA24h[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_HMBA48h_vs_HMBA24h[[names(pathfindR_outputs_HMBA48h_vs_HMBA24h)[i]]][clustered_results_HMBA48h_vs_HMBA24h[[names(pathfindR_outputs_HMBA48h_vs_HMBA24h)[i]]]$Status
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
    labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA48h_vs_HMBA24h"], ")"))
  print(cluster_enrichment_dotplots_HMBA48h_vs_HMBA24h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA48h_vs_HMBA24h)[i], "_top10_dotplot_clustered.tiff"),
         path = paste0("Pathways/HMBA48h_vs_HMBA24h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

names(enrichment_dotplots_HMBA48h_vs_HMBA24h) = names(pathfindR_outputs_HMBA48h_vs_HMBA24h)
names(cluster_enrichment_dotplots_HMBA48h_vs_HMBA24h) = names(pathfindR_outputs_HMBA48h_vs_HMBA24h)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb, subdirs[i])
  writeData(wb, subdirs[i], clustered_results_HMBA48h_vs_HMBA24h[[i]])
}
saveWorkbook(wb, file = "Pathways/HMBA48h_vs_HMBA24h/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2 = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb2, paste0(wb2_names[i], " - rep"))
  writeData(wb2, paste0(wb2_names[i], " - rep"), 
            clustered_results_HMBA48h_vs_HMBA24h[[i]][clustered_results_HMBA48h_vs_HMBA24h[[i]]$Status
                                                      == "Representative", ])
}
saveWorkbook(wb2, file = "Pathways/HMBA48h_vs_HMBA24h/pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_HMBA48h_vs_HMBA24h = list()
term_gene_graphs_HMBA48h_vs_HMBA24h = list()

for (i in 1:length(pathfindR_outputs_HMBA48h_vs_HMBA24h)){
  # term-gene heatmaps
  term_gene_heatmaps_HMBA48h_vs_HMBA24h[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h[[i]],
                                                                 genes_df = pathf_input_HMBA48h_vs_HMBA24h,
                                                                 num_terms = 10,
                                                                 use_description = TRUE,
                                                                 low = "darkblue",
                                                                 high = "red",
                                                                 mid = "white",
                                                                 pin_name_path = path2SIF)+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 2.5,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i], 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["HMBA48h_vs_HMBA24h"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  print(ggdraw(align_legend(term_gene_heatmaps_HMBA48h_vs_HMBA24h[[i]], hjust = 0.5)))
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i],
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Pathways/HMBA48h_vs_HMBA24h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_HMBA48h_vs_HMBA24h[[i]] = term_gene_graph(result_df = pathfindR_outputs_HMBA48h_vs_HMBA24h[[i]],
                                                             num_terms = 5,
                                                             use_description = TRUE,
                                                             node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 1, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i], 
                        " term - gene graph (",
                        comparisons["HMBA48h_vs_HMBA24h"], ")"))
  term_gene_graphs_HMBA48h_vs_HMBA24h[[i]]$layers = list(term_gene_graphs_HMBA48h_vs_HMBA24h[[i]]$layers[[1]], 
                                                         term_gene_graphs_HMBA48h_vs_HMBA24h[[i]]$layers[[2]], 
                                                         term_gene_graphs_HMBA48h_vs_HMBA24h[[i]]$layers[[4]])
  print(term_gene_graphs_HMBA48h_vs_HMBA24h[[i]])
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i],
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Pathways/HMBA48h_vs_HMBA24h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i], "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_HMBA48h_vs_HMBA24h = list()
for (i in 1:length(pathfindR_outputs_HMBA48h_vs_HMBA24h)){
  # UpSet plot
  UpSet_plots_HMBA48h_vs_HMBA24h[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h[[i]],
                                                   genes_df = pathf_input_HMBA48h_vs_HMBA24h,
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
         title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i], 
                        " term - UpSet plot (",
                        comparisons["HMBA48h_vs_HMBA24h"], ")"))
  tiff(paste0("Pathways/HMBA48h_vs_HMBA24h/pathfindR/", names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i], "/",
              names(wrapped_pathfindR_outputs_HMBA48h_vs_HMBA24h)[i], "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_HMBA48h_vs_HMBA24h[[i]])
  dev.off()
}

names(term_gene_graphs_HMBA48h_vs_HMBA24h) = names(pathfindR_outputs_HMBA48h_vs_HMBA24h)
names(term_gene_heatmaps_HMBA48h_vs_HMBA24h) = names(pathfindR_outputs_HMBA48h_vs_HMBA24h)
names(UpSet_plots_HMBA48h_vs_HMBA24h) = names(pathfindR_outputs_HMBA48h_vs_HMBA24h)

# HMBA72h_vs_Control48h #####
# Loading the input to pathfindR (the stage 1 vs normal topTable output):
pathf_input_HMBA72h_vs_Control48h = ashr_dgea$HMBA72h_vs_Control48h %>%
  dplyr::select(Gene.Symbol, log2FoldChange, padj) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_HMBA72h_vs_Control48h = paste0("Pathways/HMBA72h_vs_Control48h/pathfindR/", subdirs)
pathfindR_outputs_HMBA72h_vs_Control48h = list()

RNGversion("4.2.2")
set.seed(123)
for (i in 1:length(dirs_HMBA72h_vs_Control48h)){
  pathfindR_outputs_HMBA72h_vs_Control48h[[i]] = run_pathfindR(pathf_input_HMBA72h_vs_Control48h, gene_sets = "Custom",
                                                               p_val_threshold = 0.05, convert2alias = FALSE,
                                                               custom_genes = genes_of_gene_sets[[i]],
                                                               custom_descriptions = descriptions_of_gene_sets[[i]],
                                                               output_dir = dirs_HMBA72h_vs_Control48h[i], min_gset_size = 10,
                                                               max_gset_size = 300, adj_method = 'fdr',
                                                               enrichment_threshold = 0.05,
                                                               pin_name_path = path2SIF, search_method = 'GR',
                                                               grMaxDepth = 1, grSearchDepth = 1,
                                                               iterations = 10, n_processes = 10)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_HMBA72h_vs_Control48h) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

cluster_names = subdirs
RNGversion("4.2.2")
set.seed(123)
clustered_results_HMBA72h_vs_Control48h = list()
for (i in 1:length(cluster_names)){
  clustered_results_HMBA72h_vs_Control48h[[i]] = cluster_enriched_terms(pathfindR_outputs_HMBA72h_vs_Control48h[[cluster_names[[i]]]],
                                                                        method = "hierarchical")
}
names(clustered_results_HMBA72h_vs_Control48h) = cluster_names

# MSigDB Reactome     : The maximum average silhouette width was 0.35 for k = 100
# MSigDB WikiPathways : The maximum average silhouette width was 0.12 for k = 90
# MSigDB GO-BP        : The maximum average silhouette width was 0.39 for k = 450 
# MSigDB GO-CC        : The maximum average silhouette width was 0.32 for k = 40
# MSigDB GO-MF        : The maximum average silhouette width was 0.29 for k = 50
# MSigDB BioCarta     : The maximum average silhouette width was 0.17 for k = 50
# KEGG                : The maximum average silhouette width was 0.14 for k = 2

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_HMBA72h_vs_Control48h = pathfindR_outputs_HMBA72h_vs_Control48h
for (i in 1:length(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)){
  wrapped_pathfindR_outputs_HMBA72h_vs_Control48h[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_HMBA72h_vs_Control48h[[i]]$Term_Description, 
                                                                                            width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_HMBA72h_vs_Control48h = clustered_results_HMBA72h_vs_Control48h
for (i in 1:length(wrapped_clustered_pathfindR_outputs_HMBA72h_vs_Control48h)){
  wrapped_clustered_pathfindR_outputs_HMBA72h_vs_Control48h[[i]]$Term_Description = stringr::str_wrap(clustered_results_HMBA72h_vs_Control48h[[i]]$Term_Description, 
                                                                                                      width = 41)
}
rm(i)

enrichment_dotplots_HMBA72h_vs_Control48h = list()
cluster_enrichment_dotplots_HMBA72h_vs_Control48h = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_HMBA72h_vs_Control48h)){
  # unclustered results
  enrichment_dotplots_HMBA72h_vs_Control48h[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_HMBA72h_vs_Control48h[[i]],
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
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA72h_vs_Control48h"], ")"))
  print(enrichment_dotplots_HMBA72h_vs_Control48h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA72h_vs_Control48h)[i], "_top10_dotplot.tiff"),
         path = paste0("Pathways/HMBA72h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  cluster_enrichment_dotplots_HMBA72h_vs_Control48h[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_HMBA72h_vs_Control48h[[names(pathfindR_outputs_HMBA72h_vs_Control48h)[i]]][clustered_results_HMBA72h_vs_Control48h[[names(pathfindR_outputs_HMBA72h_vs_Control48h)[i]]]$Status
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
    labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA72h_vs_Control48h"], ")"))
  print(cluster_enrichment_dotplots_HMBA72h_vs_Control48h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA72h_vs_Control48h)[i], "_top10_dotplot_clustered.tiff"),
         path = paste0("Pathways/HMBA72h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

names(enrichment_dotplots_HMBA72h_vs_Control48h) = names(pathfindR_outputs_HMBA72h_vs_Control48h)
names(cluster_enrichment_dotplots_HMBA72h_vs_Control48h) = names(pathfindR_outputs_HMBA72h_vs_Control48h)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb, subdirs[i])
  writeData(wb, subdirs[i], clustered_results_HMBA72h_vs_Control48h[[i]])
}
saveWorkbook(wb, file = "Pathways/HMBA72h_vs_Control48h/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2 = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb2, paste0(wb2_names[i], " - rep"))
  writeData(wb2, paste0(wb2_names[i], " - rep"), 
            clustered_results_HMBA72h_vs_Control48h[[i]][clustered_results_HMBA72h_vs_Control48h[[i]]$Status
                                                         == "Representative", ])
}
saveWorkbook(wb2, file = "Pathways/HMBA72h_vs_Control48h/pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_HMBA72h_vs_Control48h = list()
term_gene_graphs_HMBA72h_vs_Control48h = list()

for (i in 1:length(pathfindR_outputs_HMBA72h_vs_Control48h)){
  # term-gene heatmaps
  term_gene_heatmaps_HMBA72h_vs_Control48h[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_HMBA72h_vs_Control48h[[i]],
                                                                    genes_df = pathf_input_HMBA72h_vs_Control48h,
                                                                    num_terms = 10,
                                                                    use_description = TRUE,
                                                                    low = "darkblue",
                                                                    high = "red",
                                                                    mid = "white",
                                                                    pin_name_path = path2SIF)+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 2.5,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i], 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["HMBA72h_vs_Control48h"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  print(ggdraw(align_legend(term_gene_heatmaps_HMBA72h_vs_Control48h[[i]], hjust = 0.5)))
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i],
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Pathways/HMBA72h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_HMBA72h_vs_Control48h[[i]] = term_gene_graph(result_df = pathfindR_outputs_HMBA72h_vs_Control48h[[i]],
                                                                num_terms = 5,
                                                                use_description = TRUE,
                                                                node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 1, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i], 
                        " term - gene graph (",
                        comparisons["HMBA72h_vs_Control48h"], ")"))
  term_gene_graphs_HMBA72h_vs_Control48h[[i]]$layers = list(term_gene_graphs_HMBA72h_vs_Control48h[[i]]$layers[[1]], 
                                                            term_gene_graphs_HMBA72h_vs_Control48h[[i]]$layers[[2]], 
                                                            term_gene_graphs_HMBA72h_vs_Control48h[[i]]$layers[[4]])
  print(term_gene_graphs_HMBA72h_vs_Control48h[[i]])
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i],
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Pathways/HMBA72h_vs_Control48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i], "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_HMBA72h_vs_Control48h = list()
for (i in 1:length(pathfindR_outputs_HMBA72h_vs_Control48h)){
  # UpSet plot
  UpSet_plots_HMBA72h_vs_Control48h[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_HMBA72h_vs_Control48h[[i]],
                                                      genes_df = pathf_input_HMBA72h_vs_Control48h,
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
         title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i], 
                        " term - UpSet plot (",
                        comparisons["HMBA72h_vs_Control48h"], ")"))
  tiff(paste0("Pathways/HMBA72h_vs_Control48h/pathfindR/", names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i], "/",
              names(wrapped_pathfindR_outputs_HMBA72h_vs_Control48h)[i], "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_HMBA72h_vs_Control48h[[i]])
  dev.off()
}

names(term_gene_graphs_HMBA72h_vs_Control48h) = names(pathfindR_outputs_HMBA72h_vs_Control48h)
names(term_gene_heatmaps_HMBA72h_vs_Control48h) = names(pathfindR_outputs_HMBA72h_vs_Control48h)
names(UpSet_plots_HMBA72h_vs_Control48h) = names(pathfindR_outputs_HMBA72h_vs_Control48h)

# HMBA72h_vs_HMBA24h #####
# Loading the input to pathfindR (the stage 1 vs normal topTable output):
pathf_input_HMBA72h_vs_HMBA24h = ashr_dgea$HMBA72h_vs_HMBA24h %>%
  dplyr::select(Gene.Symbol, log2FoldChange, padj) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_HMBA72h_vs_HMBA24h = paste0("Pathways/HMBA72h_vs_HMBA24h/pathfindR/", subdirs)
pathfindR_outputs_HMBA72h_vs_HMBA24h = list()

RNGversion("4.2.2")
set.seed(123)
for (i in 1:length(dirs_HMBA72h_vs_HMBA24h)){
  pathfindR_outputs_HMBA72h_vs_HMBA24h[[i]] = run_pathfindR(pathf_input_HMBA72h_vs_HMBA24h, gene_sets = "Custom",
                                                            p_val_threshold = 0.05, convert2alias = FALSE,
                                                            custom_genes = genes_of_gene_sets[[i]],
                                                            custom_descriptions = descriptions_of_gene_sets[[i]],
                                                            output_dir = dirs_HMBA72h_vs_HMBA24h[i], min_gset_size = 10,
                                                            max_gset_size = 300, adj_method = 'fdr',
                                                            enrichment_threshold = 0.05,
                                                            pin_name_path = path2SIF, search_method = 'GR',
                                                            grMaxDepth = 1, grSearchDepth = 1,
                                                            iterations = 10, n_processes = 10)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_HMBA72h_vs_HMBA24h) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

cluster_names = subdirs
RNGversion("4.2.2")
set.seed(123)
clustered_results_HMBA72h_vs_HMBA24h = list()
for (i in 1:length(cluster_names)){
  clustered_results_HMBA72h_vs_HMBA24h[[i]] = cluster_enriched_terms(pathfindR_outputs_HMBA72h_vs_HMBA24h[[cluster_names[[i]]]],
                                                                     method = "hierarchical")
}
names(clustered_results_HMBA72h_vs_HMBA24h) = cluster_names

# MSigDB Reactome     : The maximum average silhouette width was 0.33 for k = 150
# MSigDB WikiPathways : The maximum average silhouette width was 0.11 for k = 80
# MSigDB GO-BP        : The maximum average silhouette width was 0.39 for k = 350 
# MSigDB GO-CC        : The maximum average silhouette width was 0.28 for k = 40
# MSigDB GO-MF        : The maximum average silhouette width was 0.33 for k = 50
# MSigDB BioCarta     : The maximum average silhouette width was 0.20 for k = 50
# KEGG                : The maximum average silhouette width was 0.15 for k = 2

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h = pathfindR_outputs_HMBA72h_vs_HMBA24h
for (i in 1:length(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)){
  wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_HMBA72h_vs_HMBA24h[[i]]$Term_Description, 
                                                                                         width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_HMBA72h_vs_HMBA24h = clustered_results_HMBA72h_vs_HMBA24h
for (i in 1:length(wrapped_clustered_pathfindR_outputs_HMBA72h_vs_HMBA24h)){
  wrapped_clustered_pathfindR_outputs_HMBA72h_vs_HMBA24h[[i]]$Term_Description = stringr::str_wrap(clustered_results_HMBA72h_vs_HMBA24h[[i]]$Term_Description, 
                                                                                                   width = 41)
}
rm(i)

enrichment_dotplots_HMBA72h_vs_HMBA24h = list()
cluster_enrichment_dotplots_HMBA72h_vs_HMBA24h = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_HMBA72h_vs_HMBA24h)){
  # unclustered results
  enrichment_dotplots_HMBA72h_vs_HMBA24h[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h[[i]],
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
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA72h_vs_HMBA24h"], ")"))
  print(enrichment_dotplots_HMBA72h_vs_HMBA24h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA72h_vs_HMBA24h)[i], "_top10_dotplot.tiff"),
         path = paste0("Pathways/HMBA72h_vs_HMBA24h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  cluster_enrichment_dotplots_HMBA72h_vs_HMBA24h[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_HMBA72h_vs_HMBA24h[[names(pathfindR_outputs_HMBA72h_vs_HMBA24h)[i]]][clustered_results_HMBA72h_vs_HMBA24h[[names(pathfindR_outputs_HMBA72h_vs_HMBA24h)[i]]]$Status
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
    labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA72h_vs_HMBA24h"], ")"))
  print(cluster_enrichment_dotplots_HMBA72h_vs_HMBA24h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA72h_vs_HMBA24h)[i], "_top10_dotplot_clustered.tiff"),
         path = paste0("Pathways/HMBA72h_vs_HMBA24h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

names(enrichment_dotplots_HMBA72h_vs_HMBA24h) = names(pathfindR_outputs_HMBA72h_vs_HMBA24h)
names(cluster_enrichment_dotplots_HMBA72h_vs_HMBA24h) = names(pathfindR_outputs_HMBA72h_vs_HMBA24h)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb, subdirs[i])
  writeData(wb, subdirs[i], clustered_results_HMBA72h_vs_HMBA24h[[i]])
}
saveWorkbook(wb, file = "Pathways/HMBA72h_vs_HMBA24h/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2 = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb2, paste0(wb2_names[i], " - rep"))
  writeData(wb2, paste0(wb2_names[i], " - rep"), 
            clustered_results_HMBA72h_vs_HMBA24h[[i]][clustered_results_HMBA72h_vs_HMBA24h[[i]]$Status
                                                      == "Representative", ])
}
saveWorkbook(wb2, file = "Pathways/HMBA72h_vs_HMBA24h/pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_HMBA72h_vs_HMBA24h = list()
term_gene_graphs_HMBA72h_vs_HMBA24h = list()

for (i in 1:length(pathfindR_outputs_HMBA72h_vs_HMBA24h)){
  # term-gene heatmaps
  term_gene_heatmaps_HMBA72h_vs_HMBA24h[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h[[i]],
                                                                 genes_df = pathf_input_HMBA72h_vs_HMBA24h,
                                                                 num_terms = 10,
                                                                 use_description = TRUE,
                                                                 low = "darkblue",
                                                                 high = "red",
                                                                 mid = "white",
                                                                 pin_name_path = path2SIF)+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 2.5,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i], 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["HMBA72h_vs_HMBA24h"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  print(ggdraw(align_legend(term_gene_heatmaps_HMBA72h_vs_HMBA24h[[i]], hjust = 0.5)))
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i],
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Pathways/HMBA72h_vs_HMBA24h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_HMBA72h_vs_HMBA24h[[i]] = term_gene_graph(result_df = pathfindR_outputs_HMBA72h_vs_HMBA24h[[i]],
                                                             num_terms = 5,
                                                             use_description = TRUE,
                                                             node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 1, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i], 
                        " term - gene graph (",
                        comparisons["HMBA72h_vs_HMBA24h"], ")"))
  term_gene_graphs_HMBA72h_vs_HMBA24h[[i]]$layers = list(term_gene_graphs_HMBA72h_vs_HMBA24h[[i]]$layers[[1]], 
                                                         term_gene_graphs_HMBA72h_vs_HMBA24h[[i]]$layers[[2]], 
                                                         term_gene_graphs_HMBA72h_vs_HMBA24h[[i]]$layers[[4]])
  print(term_gene_graphs_HMBA72h_vs_HMBA24h[[i]])
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i],
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Pathways/HMBA72h_vs_HMBA24h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i], "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_HMBA72h_vs_HMBA24h = list()
for (i in 1:length(pathfindR_outputs_HMBA72h_vs_HMBA24h)){
  # UpSet plot
  UpSet_plots_HMBA72h_vs_HMBA24h[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h[[i]],
                                                   genes_df = pathf_input_HMBA72h_vs_HMBA24h,
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
         title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i], 
                        " term - UpSet plot (",
                        comparisons["HMBA72h_vs_HMBA24h"], ")"))
  tiff(paste0("Pathways/HMBA72h_vs_HMBA24h/pathfindR/", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i], "/",
              names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA24h)[i], "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_HMBA72h_vs_HMBA24h[[i]])
  dev.off()
}

names(term_gene_graphs_HMBA72h_vs_HMBA24h) = names(pathfindR_outputs_HMBA72h_vs_HMBA24h)
names(term_gene_heatmaps_HMBA72h_vs_HMBA24h) = names(pathfindR_outputs_HMBA72h_vs_HMBA24h)
names(UpSet_plots_HMBA72h_vs_HMBA24h) = names(pathfindR_outputs_HMBA72h_vs_HMBA24h)

# HMBA72h_vs_HMBA48h #####
# Loading the input to pathfindR (the stage 1 vs normal topTable output):
pathf_input_HMBA72h_vs_HMBA48h = ashr_dgea$HMBA72h_vs_HMBA48h %>%
  dplyr::select(Gene.Symbol, log2FoldChange, padj) %>%
  na.omit()

# Preparing a pathfindR loop for enrichment analysis
dirs_HMBA72h_vs_HMBA48h = paste0("Pathways/HMBA72h_vs_HMBA48h/pathfindR/", subdirs)
pathfindR_outputs_HMBA72h_vs_HMBA48h = list()

RNGversion("4.2.2")
set.seed(123)
for (i in 1:length(dirs_HMBA72h_vs_HMBA48h)){
  pathfindR_outputs_HMBA72h_vs_HMBA48h[[i]] = run_pathfindR(pathf_input_HMBA72h_vs_HMBA48h, gene_sets = "Custom",
                                                            p_val_threshold = 0.05, convert2alias = FALSE,
                                                            custom_genes = genes_of_gene_sets[[i]],
                                                            custom_descriptions = descriptions_of_gene_sets[[i]],
                                                            output_dir = dirs_HMBA72h_vs_HMBA48h[i], min_gset_size = 10,
                                                            max_gset_size = 300, adj_method = 'fdr',
                                                            enrichment_threshold = 0.05,
                                                            pin_name_path = path2SIF, search_method = 'GR',
                                                            grMaxDepth = 1, grSearchDepth = 1,
                                                            iterations = 10, n_processes = 10)
  cat(paste0("Done with ", subdirs[i], "\n"))
}
names(pathfindR_outputs_HMBA72h_vs_HMBA48h) = subdirs

# Perform hierarchical clustering on the results (average distance metric)
# Bear in mind the algorithm time complexity is O(n^3)

cluster_names = subdirs
RNGversion("4.2.2")
set.seed(123)
clustered_results_HMBA72h_vs_HMBA48h = list()
for (i in 1:length(cluster_names)){
  clustered_results_HMBA72h_vs_HMBA48h[[i]] = cluster_enriched_terms(pathfindR_outputs_HMBA72h_vs_HMBA48h[[cluster_names[[i]]]],
                                                                     method = "hierarchical")
}
names(clustered_results_HMBA72h_vs_HMBA48h) = cluster_names

# MSigDB Reactome     : The maximum average silhouette width was 0.35 for k = 250
# MSigDB WikiPathways : The maximum average silhouette width was 0.14 for k = 150
# MSigDB GO-BP        : The maximum average silhouette width was 0.41 for k = 550 
# MSigDB GO-CC        : The maximum average silhouette width was 0.3 for k = 80
# MSigDB GO-MF        : The maximum average silhouette width was 0.35 for k = 80
# MSigDB BioCarta     : The maximum average silhouette width was 0.23 for k = 30
# KEGG                : The maximum average silhouette width was 0.13 for k = 5

# Wrapping the text of terms with too many characters in their description
wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h = pathfindR_outputs_HMBA72h_vs_HMBA48h
for (i in 1:length(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)){
  wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h[[i]]$Term_Description = stringr::str_wrap(pathfindR_outputs_HMBA72h_vs_HMBA48h[[i]]$Term_Description, 
                                                                                         width = 41)
}
rm(i)

wrapped_clustered_pathfindR_outputs_HMBA72h_vs_HMBA48h = clustered_results_HMBA72h_vs_HMBA48h
for (i in 1:length(wrapped_clustered_pathfindR_outputs_HMBA72h_vs_HMBA48h)){
  wrapped_clustered_pathfindR_outputs_HMBA72h_vs_HMBA48h[[i]]$Term_Description = stringr::str_wrap(clustered_results_HMBA72h_vs_HMBA48h[[i]]$Term_Description, 
                                                                                                   width = 41)
}
rm(i)

enrichment_dotplots_HMBA72h_vs_HMBA48h = list()
cluster_enrichment_dotplots_HMBA72h_vs_HMBA48h = list()

# Producing dotplots with the results
for (i in 1:length(pathfindR_outputs_HMBA72h_vs_HMBA48h)){
  # unclustered results
  enrichment_dotplots_HMBA72h_vs_HMBA48h[[i]] = enrichment_chart(result_df = wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h[[i]],
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
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA72h_vs_HMBA48h"], ")"))
  print(enrichment_dotplots_HMBA72h_vs_HMBA48h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA72h_vs_HMBA48h)[i], "_top10_dotplot.tiff"),
         path = paste0("Pathways/HMBA72h_vs_HMBA48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # clustered results
  cluster_enrichment_dotplots_HMBA72h_vs_HMBA48h[[i]] = enrichment_chart(result_df = wrapped_clustered_pathfindR_outputs_HMBA72h_vs_HMBA48h[[names(pathfindR_outputs_HMBA72h_vs_HMBA48h)[i]]][clustered_results_HMBA72h_vs_HMBA48h[[names(pathfindR_outputs_HMBA72h_vs_HMBA48h)[i]]]$Status
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
    labs(title = paste0("Top 10 clustered ", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i],
                        " terms enrichment dotplot - (", comparisons["HMBA72h_vs_HMBA48h"], ")"))
  print(cluster_enrichment_dotplots_HMBA72h_vs_HMBA48h[[i]])
  ggsave(filename = paste0(names(pathfindR_outputs_HMBA72h_vs_HMBA48h)[i], "_top10_dotplot_clustered.tiff"),
         path = paste0("Pathways/HMBA72h_vs_HMBA48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

names(enrichment_dotplots_HMBA72h_vs_HMBA48h) = names(pathfindR_outputs_HMBA72h_vs_HMBA48h)
names(cluster_enrichment_dotplots_HMBA72h_vs_HMBA48h) = names(pathfindR_outputs_HMBA72h_vs_HMBA48h)

# Write out results in a comprehensive .xlsx file
wb = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb, subdirs[i])
  writeData(wb, subdirs[i], clustered_results_HMBA72h_vs_HMBA48h[[i]])
}
saveWorkbook(wb, file = "Pathways/HMBA72h_vs_HMBA48h/pathfindR/Comprehensive_pathfindR_output.xlsx",
             overwrite = TRUE); rm(wb)

# Representative terms
wb2 = createWorkbook()
for(i in 1:length(subdirs)) {
  addWorksheet(wb2, paste0(wb2_names[i], " - rep"))
  writeData(wb2, paste0(wb2_names[i], " - rep"), 
            clustered_results_HMBA72h_vs_HMBA48h[[i]][clustered_results_HMBA72h_vs_HMBA48h[[i]]$Status
                                                      == "Representative", ])
}
saveWorkbook(wb2, file = "Pathways/HMBA72h_vs_HMBA48h/pathfindR/Representative_terms.xlsx",
             overwrite = TRUE); rm(wb2)

# Term-gene heatmaps and term-gene graphs #####
term_gene_heatmaps_HMBA72h_vs_HMBA48h = list()
term_gene_graphs_HMBA72h_vs_HMBA48h = list()

for (i in 1:length(pathfindR_outputs_HMBA72h_vs_HMBA48h)){
  # term-gene heatmaps
  term_gene_heatmaps_HMBA72h_vs_HMBA48h[[i]] = term_gene_heatmap(result_df = wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h[[i]],
                                                                 genes_df = pathf_input_HMBA72h_vs_HMBA48h,
                                                                 num_terms = 10,
                                                                 use_description = TRUE,
                                                                 low = "darkblue",
                                                                 high = "red",
                                                                 mid = "white",
                                                                 pin_name_path = path2SIF)+
    theme(plot.title = element_text(size = 4.5, face = "bold", vjust = 2, hjust = 0.5),
          axis.text.y = element_text(color = "black", size = 3.5),
          axis.text.x = element_text(color = "black", size = 2.5,  vjust = 0.5),
          axis.title.x = element_text(size = 5, face = "bold"),
          legend.key.size = unit(3, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical") +
    scale_fill_gradient2(name = "logFC",
                         low = "darkblue", mid = "white", na.value = "white",
                         high = "red") +
    labs(title = paste0("Top 10 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i], 
                        " terms - differentially expressed genes heatmap (",
                        comparisons["HMBA72h_vs_HMBA48h"], ")"),
         fill = expression(log[2] ~ "FoldChange"))
  print(ggdraw(align_legend(term_gene_heatmaps_HMBA72h_vs_HMBA48h[[i]], hjust = 0.5)))
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i],
                           "_top10_term_gene_heatmap.tiff"),
         path = paste0("Pathways/HMBA72h_vs_HMBA48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i], "/"),
         width = 2880, height = 1620, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
  
  # term-gene graphs
  term_gene_graphs_HMBA72h_vs_HMBA48h[[i]] = term_gene_graph(result_df = pathfindR_outputs_HMBA72h_vs_HMBA48h[[i]],
                                                             num_terms = 5,
                                                             use_description = TRUE,
                                                             node_size = "p_val")+
    aes(max.overlaps = 5)+
    scale_size(range = c(1, 3)) +
    suppressWarnings(ggraph::geom_node_text(ggplot2::aes_(label = ~name),  nudge_y = .1,
                                            repel = TRUE, size = 1, max.overlaps = 10, check_overlap = T))+
    theme(plot.title = element_text(size = 6, face = "bold", hjust = 0.5, vjust = 0.5),
          plot.background = element_rect(fill = "white"),
          plot.subtitle = element_text(size = 4.5, face = "italic", hjust = 0.5, vjust = 1.5),
          legend.key.size = unit(2, units = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.title = element_text(size = 4, face = "bold"),
          legend.text = element_text(size = 4),
          legend.title.align = 0.5,
          legend.direction = "vertical")+
    labs(title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i], 
                        " term - gene graph (",
                        comparisons["HMBA72h_vs_HMBA48h"], ")"))
  term_gene_graphs_HMBA72h_vs_HMBA48h[[i]]$layers = list(term_gene_graphs_HMBA72h_vs_HMBA48h[[i]]$layers[[1]], 
                                                         term_gene_graphs_HMBA72h_vs_HMBA48h[[i]]$layers[[2]], 
                                                         term_gene_graphs_HMBA72h_vs_HMBA48h[[i]]$layers[[4]])
  print(term_gene_graphs_HMBA72h_vs_HMBA48h[[i]])
  ggsave(filename = paste0(names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i],
                           "_top5_term_gene_graph.tiff"),
         path = paste0("Pathways/HMBA72h_vs_HMBA48h/pathfindR/", 
                       names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i], "/"),
         width = 1920*2.5, height = 1080*2.5, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# UpSet plots #####
UpSet_plots_HMBA72h_vs_HMBA48h = list()
for (i in 1:length(pathfindR_outputs_HMBA72h_vs_HMBA48h)){
  # UpSet plot
  UpSet_plots_HMBA72h_vs_HMBA48h[[i]] = UpSet_plot(result_df = wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h[[i]],
                                                   genes_df = pathf_input_HMBA72h_vs_HMBA48h,
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
         title = paste0("Top 5 ", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i], 
                        " term - UpSet plot (",
                        comparisons["HMBA72h_vs_HMBA48h"], ")"))
  tiff(paste0("Pathways/HMBA72h_vs_HMBA48h/pathfindR/", names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i], "/",
              names(wrapped_pathfindR_outputs_HMBA72h_vs_HMBA48h)[i], "_top5_UpSet_plot.tiff"), 
       width = 2880, height = 1620*4, res = 700, compression = "lzw")
  print(UpSet_plots_HMBA72h_vs_HMBA48h[[i]])
  dev.off()
}

names(term_gene_graphs_HMBA72h_vs_HMBA48h) = names(pathfindR_outputs_HMBA72h_vs_HMBA48h)
names(term_gene_heatmaps_HMBA72h_vs_HMBA48h) = names(pathfindR_outputs_HMBA72h_vs_HMBA48h)
names(UpSet_plots_HMBA72h_vs_HMBA48h) = names(pathfindR_outputs_HMBA72h_vs_HMBA48h)