################################################################################
# Name: downloadTCGA-BRCA                                                      #
#                                                                              #
# Author: Jean Resende (jean.s.s.resende@gmail.com)                            #
# Creation date: 2023/07/06                                                    #
# Last update date: 2023/07/06                                                 #
#                                                                              #
# Descrption: Download of TCGA-BRCA cohort counts                              #
################################################################################

# -- download of samples -- ####################################################
library(TCGAbiolinks)

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")

GDCdownload(query)
tcgaBRCA <- GDCprepare(query, summarizedExperiment = TRUE) # 60660 genes
save(tcgaBRCA, file = "tcgaBRCA_2023.RData")

rm(list = ls())

# -- Filter for lncRNA -- ######################################################

load("tcgaBRCA_2023.RData")

table(tcgaBRCA@rowRanges$gene_type)

lncRNA <- tcgaBRCA[tcgaBRCA@rowRanges$gene_type=="lncRNA",]
save(lncRNA, file = "lncRNA_tcgaBRCA_2023.RData")
################################################################################
