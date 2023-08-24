################################################################################
# Name: plotHeatmap                                                            #
#                                                                              #
# Author: jean resende (jean.s.s.resende@gmail.com)                            #
# Date: 2023/08/24                                                             #
# Project: lncRNA_TCGA                                                         #
################################################################################

# -- Funcoes -- ################################################################
download.TPM.TCGA <- function(project){
  require(TCGAbiolinks)
  query <- GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "STAR - Counts")
  
  GDCdownload(query)
  
  tcgaProject <- GDCprepare(query, summarizedExperiment = TRUE) # 60660 genes
  save(tcgaProject, file = paste(project,".RData",sep = ""))
  
  rm(list = ls())
}

# -- baixar os dados -- ########################################################

# -- download dos dados
download.TPM.TCGA(project = "TCGA-PRAD")

# -- automatizando para os canceres do TCGA
project_tcga <- c("TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
                  "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
                  "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LGG",
                  "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV",
                  "TCGA-PAAD", "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC",
                  "TCGA-SKCM", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM",
                  "TCGA-UCEC", "TCGA-UCS", "TCGA-UVM")

for (i in project_tcga) {
  download.TPM.TCGA(project = i)
}

# -- carregando objetos
load("TCGA-PRAD.RData")

lnc_bruno <- c("ENSG00000224592.5","ENSG00000286552.1","ENSG00000227290.2",
               "ENSG00000237505.8","ENSG00000232184.1","ENSG00000270087.7",
               "ENSG00000230967.2","ENSG00000254789.4","ENSG00000254946.2",
               "ENSG00000228061.7","ENSG00000257379.1","ENSG00000282977.1",
               "ENSG00000273046.1","ENSG00000257545.5","ENSG00000234377.8",
               "ENSG00000287635.1","ENSG00000224243.2","ENSG00000257986.3",
               "ENSG00000257126.6","ENSG00000257523.3","ENSG00000258038.6",
               "ENSG00000258107.3","ENSG00000248975.2","ENSG00000287153.1",
               "ENSG00000257585.1","ENSG00000248079.4","ENSG00000259673.6",
               "ENSG00000259768.8","ENSG00000233101.10","ENSG00000260372.7",
               "ENSG00000228835.1","ENSG00000285940.4","ENSG00000286716.1",
               "ENSG00000238133.6","ENSG00000224577.4","ENSG00000233723.11",
               "ENSG00000233891.8","ENSG00000271955.2","ENSG00000285611.2",
               "ENSG00000231609.7","ENSG00000204929.13","ENSG00000228655.7",
               "ENSG00000257277.2","ENSG00000232377.1","ENSG00000226383.7",
               "ENSG00000283436.1","ENSG00000286889.1","ENSG00000287048.1",
               "ENSG00000286530.1","ENSG00000235724.10","ENSG00000237844.2",
               "ENSG00000234350.7","ENSG00000289413.1","ENSG00000228956.9",
               "ENSG00000240405.7","ENSG00000243620.2","ENSG00000242536.2",
               "ENSG00000241479.1","ENSG00000145075.13","ENSG00000285336.1",
               "ENSG00000250467.1","ENSG00000288692.1","ENSG00000248118.2",
               "ENSG00000249483.1","ENSG00000287862.1","ENSG00000247828.10",
               "ENSG00000245526.13","ENSG00000248309.9","ENSG00000237187.9",
               "ENSG00000271860.9","ENSG00000232790.3","ENSG00000285592.1",
               "ENSG00000254369.6","ENSG00000273433.1","ENSG00000224223.1",
               "ENSG00000231764.11","ENSG00000224595.1","ENSG00000240973.1",
               "ENSG00000230316.8","ENSG00000254102.2","ENSG00000205636.3",
               "ENSG00000235298.1","ENSG00000234323.10","ENSG00000233033.1",
               "ENSG00000288597.1")

idx <- match(gsub("\\..*","",lnc_bruno),
             gsub("\\..*","",tcgaProject@rowRanges$gene_id))
idx <- idx[!is.na(idx)]

barc_tpnt <- tcgaProject$barcode[tcgaProject$shortLetterCode == "TP" |
                                   tcgaProject$shortLetterCode == "NT"]

tpm <- tcgaProject@assays@data$tpm_unstrand[idx, barc_tpnt]

rownames(tpm) <- tcgaProject@rowRanges$gene_id[idx]
colnames(tpm) <- barc_tpnt

coldata <- tcgaProject@colData[tcgaProject$barcode %in% barc_tpnt,
                               c("barcode","shortLetterCode")]

coldata$shortLetterCode <- as.factor(coldata$shortLetterCode)

rm(tcgaProject,barc_tpnt,idx,idxPac,lnc_bruno)
################################################################################








