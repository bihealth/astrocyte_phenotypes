library(tidyr)
library(dplyr)
library(DESeq2)
library(factoextra)
library(readxl)
library(tmod)
library(tidyr)
library(tibble)
library(pca3d)


sample_annot<-read.table("../additional_data/sample_parameter.csv", sep=",", header=T)
load("../additional_data/feature_annot.rdata")
load("../additional_data/deseq2_results.rda")

 synapse_pr <- readxl::read_excel("../additional_data/GO_term_summary_20200131_074208.xlsx") %>% tbl_df()
 gap_junction <- readxl::read_excel("../additional_data/GO_term_summary_20200130_073756.xlsx") %>% tbl_df()
 glial_cell_migration <- readxl::read_excel("../additional_data/GO_term_summary_20200130_064709.xlsx") %>% tbl_df()
 lactate <- readxl::read_excel("../additional_data/GO_term_summary_20200123_114803.xlsx") %>% tbl_df()
 calcium <- readxl::read_excel("../additional_data/GO_term_summary_20200123_115022.xlsx") %>% tbl_df()
 glial_cell_prolif <- readxl::read_excel("../additional_data/GO_term_summary_20200211_074931.xlsx") %>% tbl_df()
 glucose<- readxl::read_excel("../additional_data/GO_term_summary_20200123_114921.xlsx") %>% tbl_df()
 astrocyte_markers <- readxl::read_xls("../additional_data/Astrocyte_markers.xls", col_names = "gene_name") %>% tbl_df()




get_reg_log_counts <- function(dds, blind = T, ...){

  rld <- DESeq2::rlog(dds, blind = blind, ...)
  sample_annot <- get_sample_table(dds)
  df_reg_log_count <-
  SummarizedExperiment::assay(rld) %>%
    as.data.frame() %>%
    tibble::rownames_to_column('feature_id') %>%
    tibble::as_data_frame() %>%
    tidyr::gather(sample_id, reg_log_count, -feature_id)

  if (sample_annot %has_name% "sample_name") {
    df_reg_log_count <-
    dplyr::left_join(df_reg_log_count,
                     sample_annot[, c("sample_id", "sample_name")],
                     by = "sample_id") %>%
    dplyr::select(feature_id, sample_id, sample_name, reg_log_count)
  }

  return(df_reg_log_count)
}

get_sample_table <- function(dds){
  sample_table <-
    SummarizedExperiment::colData(dds) %>%
    SummarizedExperiment::as.data.frame() %>%
    tibble::as_data_frame()
  return(sample_table)
}

rld <- DESeq2::rlog(dds, blind = TRUE)

rld_df <- get_reg_log_counts(dds, blind = T)

rld_df1 <- rld_df %>% dplyr::select(-sample_name) %>% spread(sample_id, reg_log_count)

glucose_rld_pca <- rld_df1  %>% left_join(feature_annot, by=c("feature_id"="gene_id")) %>% filter(gene_name %in% glucose$Symbol)  %>% dplyr::select(-gene_name) %>% tibble::column_to_rownames("feature_id") %>% as.data.frame() %>% as.matrix() %>% t() %>% prcomp()
fviz_pca_ind(glucose_rld_pca, axes=c(1,2), repel=T, palette = "Dark2", legend="none")
  

lactate_rld_pca <- rld_df1  %>% left_join(feature_annot, by=c("feature_id"="gene_id")) %>% filter(gene_name %in% lactate$Symbol)  %>% dplyr::select(-gene_name) %>% tibble::column_to_rownames("feature_id") %>% as.data.frame() %>% as.matrix() %>% t() %>% prcomp()
fviz_pca_ind(lactate_rld_pca, axes=c(1,2), repel=T, palette = "Dark2", legend="none")


calcium_rld_pca <- rld_df1  %>% left_join(feature_annot, by=c("feature_id"="gene_id")) %>% filter(gene_name %in% calcium$Symbol)  %>% dplyr::select(-gene_name) %>% tibble::column_to_rownames("feature_id") %>% as.data.frame() %>% as.matrix() %>% t() %>% prcomp()
fviz_pca_ind(calcium_rld_pca, axes=c(1,2), repel=T, palette = "Dark2", legend="none")

synapse_pr_rld_pca <- rld_df1  %>% left_join(feature_annot, by=c("feature_id"="gene_id")) %>% filter(gene_name %in% synapse_pr$Symbol)  %>% dplyr::select(-gene_name) %>% tibble::column_to_rownames("feature_id") %>% as.data.frame() %>% as.matrix() %>% t() %>% prcomp()
fviz_pca_ind(synapse_pr_rld_pca, axes=c(1,2), repel=T, palette = "Dark2", legend="none")


gap_junction_rld_pca <- rld_df1  %>% left_join(feature_annot, by=c("feature_id"="gene_id")) %>% filter(gene_name %in% gap_junction$Symbol)  %>% dplyr::select(-gene_name) %>% tibble::column_to_rownames("feature_id") %>% as.data.frame() %>% as.matrix() %>% t() %>% prcomp()
fviz_pca_ind(gap_junction_rld_pca, axes=c(1,2), repel=T, palette = "Dark2", legend="none")

glial_cell_migration_rld_pca <- rld_df1  %>% left_join(feature_annot, by=c("feature_id"="gene_id")) %>% filter(gene_name %in% glial_cell_migration$Symbol)  %>% dplyr::select(-gene_name) %>% tibble::column_to_rownames("feature_id") %>% as.data.frame() %>% as.matrix() %>% t() %>% prcomp()

fviz_pca_ind(glial_cell_migration_rld_pca, axes=c(1,2), repel=T, palette = "Dark2", legend="none")  

glial_cell_prolif_rld_pca <- rld_df1  %>% left_join(feature_annot, by=c("feature_id"="gene_id")) %>% filter(gene_name %in% glial_cell_prolif$Symbol)  %>% dplyr::select(-gene_name) %>% tibble::column_to_rownames("feature_id") %>% as.data.frame() %>% as.matrix() %>% t() %>% prcomp()

fviz_pca_ind(glial_cell_prolif_rld_pca, axes=c(1,2), repel=T, palette = "Dark2", legend="none")

astro_m_rld_pca <- rld_df1  %>% left_join(feature_annot, by=c("feature_id"="gene_id")) %>% filter(gene_name %in% astrocyte_markers$gene_name)  %>% dplyr::select(-gene_name, -description, -gene_type, -feature_length) %>% tibble::column_to_rownames("feature_id") %>% as.data.frame() %>% as.matrix() %>% t() %>% prcomp()
fviz_pca_ind(astro_m_rld_pca, axes=c(1,2), repel=T, palette = "Dark2", legend="none")




rld_forTMOD<-assay(rld) %>% as.data.frame() %>% rownames_to_column("gene_id") %>% left_join(feature_annot) %>% dplyr::select(-gene_id) %>% dplyr::select(gene_name, everything())
E <- as.matrix(rld_forTMOD[,-c(1)])
pca <- prcomp(t(E), scale.=TRUE)
group <- rep(c("AdAC_CellC", "AdAC_direct", "AGES", "neoAC_CellC", "neoAC_direct", "NSC"), each=4)



pdf("pca_var2.pdf")
par(mfrow=c(2,2))
f1<-fviz_eig(pca)
f2<-fviz_pca_ind(pca, axes=c(1,2), repel=T,habillage=group,addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2", legend="none")
f3<-fviz_pca_ind(pca, axes=c(3,4), repel=T,habillage=group,addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2", legend="none")
f4<-fviz_pca_ind(pca, axes=c(5,6), repel=T,habillage=group,addEllipses=TRUE, ellipse.level=0.95, palette = "Dark2", legend="none")
cowplot::plot_grid(f1, f2, f3, f4, labels = c('A', 'B', 'C', 'D'), label_size = 12)
dev.off()

data(tmod)
source("tmod_functions.R")


tmod_db_path <- "../additional_data/tmod_db"
tmod_sort_order <- "pvalue"
tmod_qval_max <- .01
tmod_auc_min  <- .75
tmod_n_max    <- 25
tmod_func <- tmodCERNOtest
tmod_params <- list()
tmod_dbs <- list(MSigDB=readRDS(file.path(tmod_db_path, "msigdb_HS_v6.2.rds")),BTM=tmod)


tmod_dbs_desc <- list(
  MSigDB="MSigDB is the database of gene sets from the Broad Institute,
          used by default by the GSEA program. The database includes a
          manually curated set of 50 'hallmark' gene sets, the Gene
          Ontology (GO) gene sets and many others. For a full list of
          gene sets included here, see the [MSigDB website](http://software.broadinstitute.org/gsea/msigdb/index.jsp).
          This version contains only gene sets not contained in other
          databases (such as GO, REACTOME, Hallmark etc.)",
  BTM="BTM are blood transcriptional modules derived from human blood gene
       expression by means of co-expression analysis. That is, genes
       present in a gene set show a similar transcriptional behavior under
       a variety of conditions and diseases in human peripherial blood
       cells. The BTMs used here have been collected from two sources; see
       tmod user manual for details.")

tmod_mapping <- list(
  default=feature_annot %>% dplyr::select(from_gene_id=gene_id, to_gene_id=gene_name)
  )
tmod_mapping$default$to_gene_id = toupper(tmod_mapping$default$to_gene_id)

## Additional entries may be created on the go.
msigdb <- tmod_dbs$MSigDB
tmod_dbs$Hallmark <- msigdb[ msigdb$MODULES$Category == "H" ]
tmod_dbs_desc$Hallmark <- "The Hallmark gene sets from MSigDB"

tmod_dbs$GO <- msigdb[ msigdb$MODULES$Category == "C5" ]
tmod_dbs_desc$GO <- "Gene Ontology data sets from MSigDB"

for(db in c("CP:REACTOME", "CP:KEGG", "CP:BIOCARTA", "CP")) {
  dbname <- gsub("CP:", "", db)
  tmod_dbs[[dbname]] <- msigdb[ msigdb$MODULES$Subcategory == db ]
  tmod_dbs_desc[[dbname]] <- paste0(dbname, " canonical pathway gene sets from MSigDB")
}

## remove these data sets from msigdb
tmod_dbs$MSigDB <- msigdb[ ! (msigdb$MODULES$Category %in% c("H", "C5", "C2")) ]


## Load the tmod databases
.tmod_load_db <- function(x, y) {
  if(class(x) == "tmod") return(x)
  if(is.character(x))    return(readRDS(x))
  stop(sprintf("tmod: Cannot recognize database %s", y))
}

tmod_dbs <- imap(tmod_dbs, .tmod_load_db)

go <- msigdb[ msigdb$MODULES$Category == "C5" ]



tmod_all_results_12<-list()
tmod_all_results_34<-list()


for (i in names(tmod_dbs)){
  print (i)
  tmod_all_results_12[[i]]<-tmodPCA(pca, genes=toupper(rld_forTMOD$gene_name), components=1:2, plot.params=list(group=group, legend="topright"), mode="leftbottom", plotfunc = pca2d, mset=tmod_dbs[[i]], plot=F, qval=0.05)
  tmod_all_results_34[[i]]<-tmodPCA(pca, genes=toupper(rld_forTMOD$gene_name), components=3:4, plot.params=list(group=group, legend="topright"), mode="leftbottom", plotfunc = pca2d, mset=tmod_dbs[[i]], plot=F, qval=0.05)
}

# only GO

tmodPCA_results <- tmodPCA(pca, genes=toupper(rld_forTMOD$gene_name), components=1:2, plot.params=list(group=group, legend="topright"), mode="leftbottom", plotfunc = pca2d, mset=go, plot=F, qval=0.05)

pca_weights <- pca$rotation %>% tbl_df() %>% mutate(gene_name=rld_forTMOD$gene_name)

evidencePlot(m='M13456', l=toupper((pca_weights %>% arrange(-PC1))$gene_name) , mset=go, gene.labels=T, main="DNA replication initiation\n Dim1:right", cex.main=0.7, cex.lab=0.7, gl.cex=0.7)
evidencePlot(m='M13478', l=toupper((pca_weights %>% arrange(PC1))$gene_name), mset=go, gene.labels=T, main="Cilium or flagellum dependent cell motility\n Dim1:left", cex.main=0.7, cex.lab=0.7, gl.cex=0.7)
evidencePlot(m='M10697', l=toupper((pca_weights %>% arrange(-PC2))$gene_name), mset=go, gene.labels=T, main="Cell cycle dna replication\n Dim2:top", cex.main=0.8, cex.lab=0.7, gl.cex=0.7)
evidencePlot(m='M13707', l=toupper((pca_weights %>% arrange(PC2))$gene_name), mset=go, gene.labels=T, main="Rregulation of synaptic vesicle exocytosis\n Dim2:bottom", cex.main=0.8, cex.lab=0.7, gl.cex=0.7)
evidencePlot(m='M12468', l=toupper((pca_weights %>% arrange(-PC3))$gene_name), mset=go, gene.labels=T, main="Regulation of macrophage chemotaxis\n Dim3:right", cex.main=0.8, cex.lab=0.7, gl.cex=0.7)
evidencePlot(m='M10697', l=toupper((pca_weights %>% arrange(PC3))$gene_name), mset=go, gene.labels=T, main="Cell cycle dna replication\n Dim3:left", cex.main=0.8, cex.lab=0.7, gl.cex=0.7)
evidencePlot(m='M12248', l=toupper((pca_weights %>% arrange(-PC4))$gene_name), mset=go, gene.labels=T, main="Retrograde protein transport er to cytosol\n Dim4:top", cex.main=0.7, cex.lab=0.7, gl.cex=0.7)
evidencePlot(m='M10159', l=toupper((pca_weights %>% arrange(PC4))$gene_name), mset=go, gene.labels=T, main="Neuronal stem cell division\nDim4:bottom", cex.main=0.8, cex.lab=0.7, gl.cex=0.7)



tm_go12 <-list()


tm_go12$Component1.left <- tmod_all_results_12$GO$enrichments$Component1.left %>% filter(AUC>0.75)%>%  head(10)
tm_go12$Component1.left$Title <- sub("Go ","", tm_go12$Component1.left$Title )

tm_go12$Component1.right <- tmod_all_results_12$GO$enrichments$Component1.right %>% filter(AUC>0.75)%>% head(10)
tm_go12$Component1.right$Title <- sub("Go ","", tm_go12$Component1.right$Title )


tm_go12$Component2.top <- tmod_all_results_12$GO$enrichments$Component2.top %>% filter(AUC>0.75)%>% head(10)
tm_go12$Component2.top$Title <- sub("Go ","", tm_go12$Component2.top$Title )

tm_go12$Component2.bottom <- tmod_all_results_12$GO$enrichments$Component2.bottom %>% filter(AUC>0.75)%>% head(10)
tm_go12$Component2.bottom$Title <- sub("Go ","", tm_go12$Component2.bottom$Title )



tm_go34 <-list()


tm_go34$Component3.left <- tmod_all_results_34$GO$enrichments$Component3.left %>% filter(AUC>0.75)%>%  head(10) 
tm_go34$Component3.left$Title <- sub("Go ","", tm_go34$Component3.left$Title )
tm_go34$Component3.right <- tmod_all_results_34$GO$enrichments$Component3.right %>% filter(AUC>0.75)%>% head(10)
tm_go34$Component3.right$Title <- sub("Go ","", tm_go34$Component3.right$Title )
tm_go34$Component4.top <- tmod_all_results_34$GO$enrichments$Component4.top %>% filter(AUC>0.75)%>% head(10)
tm_go34$Component4.top$Title <- sub("Go ","", tm_go34$Component4.top$Title )
tm_go34$Component4.bottom <- tmod_all_results_34$GO$enrichments$Component4.bottom %>% filter(AUC>0.75)%>% head(10)
tm_go34$Component4.bottom$Title <- sub("Go ","", tm_go34$Component4.bottom$Title )

tm_go34$Component4.bottom$Title[5]="reg. of AMPA select. glutamate receptor activity"






tmodPanelPlot(tm_go12)
tmodPanelPlot(tm_go34)

contrasts_g<-c("AdAC_CellC_vs_AdAC_direct", "AGES_vs_AdAC_CellC", "AGES_vs_neoAC_CellC", "AGES_vs_NSC", "neoAC_CellC_vs_AdAC", "neoAC_CellC_vs_AdAC", "neoAC_CellC_vs_neoAC_direct", "neoAC_direct_vs_AdAC_direct")


tmod_ab_res_GO <- list()
tmod_ab_pie_GO <- list()

for (i in contrasts_g){
    res <- readRDS(paste("../additional_data/",i,"_results.rds", sep=""))
    res <- res %>% as.data.frame() %>%  tibble::rownames_to_column("gene_id") %>% tbl_df() %>% arrange(padj) %>% left_join(feature_annot)
    res_genes <- toupper(res$gene_name)
    res_lfc <- res$log2FoldChange
    res_pval <- res$padj
    tmod_ab_pie_GO[[i]] <- tmodDecideTests(res_genes, res_lfc, res_pval, mset=go) 
    tmod_ab_res_GO[[i]] <- tmodCERNOtest(res_genes, mset=go)
}




# or read the objects from supplementary data

tmod_ab_res_GO <- readRDS("../additional_data/tmod_results_GO.rds")
tmod_ab_pie_GO <- readRDS("../additional_data/tmod_pie_GO.rds")


go_modules<- unique(c("M13456","M13478","M10902", "M11222", "M16412", "M13478", "M13478", "M16412", "M14460","M12441", "M12183", "M13469", "M13169","M13810", "M14007", "M11199", "M10697", "M13169","M12526","M15091"))



tmodPanelPlot(tmod_ab_res_GO, pie = tmod_ab_pie_GO, select = go_modules)


