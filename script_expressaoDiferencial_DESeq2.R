library(DESeq2)

load("lncRNA_tcgaBRCA_2023.RData")

nrow(lncRNA)
ncol(lncRNA)

# -- preparando o objeto
## -- permanecer somente coom TP e NT
table(lncRNA$shortLetterCode)
lncRNA_TPNT <- lncRNA[,lncRNA$shortLetterCode != "TM"]
table(lncRNA_TPNT$shortLetterCode)

## -- coldata
class(lncRNA_TPNT@colData)
coldata <- as.data.frame(lncRNA_TPNT@colData)
coldata$shortLetterCode <- factor(coldata$shortLetterCode, levels = c("NT","TP"))

## -- matriz
mat_TPNT <- lncRNA_TPNT@assays@data$unstranded
rownames(mat_TPNT) <- rownames(lncRNA_TPNT)
colnames(mat_TPNT) <- lncRNA_TPNT$barcode

## -- deseq
dds_TPNT <- DESeqDataSetFromMatrix(countData = mat_TPNT,
                                   colData = coldata,
                                   design = ~ shortLetterCode)

# -- pre-filtragem
#keep <- rowSums(counts(dds_TPNT)) >= 10
#dds <- dds_TPNT[keep,]

dds <- dds_TPNT

# -- setando a refencia
dds$shortLetterCode <- relevel(dds$shortLetterCode, ref = "NT")

# -- expressao diferencial
dds <- DESeq(dds)
res_padj05 <- results(dds, contrast = c("shortLetterCode","TP","NT"), alpha = 0.05)
head(res_padj05)
write.csv(res_padj05, file = "differentialExpression_TPNT_res05.csv")

lnc_bruno <- c("ENSG00000224592.5","ENSG00000286552.1","ENSG00000227290.2","ENSG00000237505.8","ENSG00000232184.1","ENSG00000270087.7","ENSG00000230967.2","ENSG00000254789.4","ENSG00000254946.2","ENSG00000228061.7","ENSG00000257379.1","ENSG00000282977.1","ENSG00000273046.1","ENSG00000257545.5","ENSG00000234377.8","ENSG00000287635.1","ENSG00000224243.2","ENSG00000257986.3","ENSG00000257126.6","ENSG00000257523.3","ENSG00000258038.6","ENSG00000258107.3","ENSG00000248975.2","ENSG00000287153.1","ENSG00000257585.1","ENSG00000248079.4","ENSG00000259673.6","ENSG00000259768.8","ENSG00000233101.10","ENSG00000260372.7","ENSG00000228835.1","ENSG00000285940.4","ENSG00000286716.1","ENSG00000238133.6","ENSG00000224577.4","ENSG00000233723.11","ENSG00000233891.8","ENSG00000271955.2","ENSG00000285611.2","ENSG00000231609.7","ENSG00000204929.13","ENSG00000228655.7","ENSG00000257277.2","ENSG00000232377.1","ENSG00000226383.7","ENSG00000283436.1","ENSG00000286889.1","ENSG00000287048.1","ENSG00000286530.1","ENSG00000235724.10","ENSG00000237844.2","ENSG00000234350.7","ENSG00000289413.1","ENSG00000228956.9","ENSG00000240405.7","ENSG00000243620.2","ENSG00000242536.2","ENSG00000241479.1","ENSG00000145075.13","ENSG00000285336.1","ENSG00000250467.1","ENSG00000288692.1","ENSG00000248118.2","ENSG00000249483.1","ENSG00000287862.1","ENSG00000247828.10","ENSG00000245526.13","ENSG00000248309.9","ENSG00000237187.9","ENSG00000271860.9","ENSG00000232790.3","ENSG00000285592.1","ENSG00000254369.6","ENSG00000273433.1","ENSG00000224223.1","ENSG00000231764.11","ENSG00000224595.1","ENSG00000240973.1","ENSG00000230316.8","ENSG00000254102.2","ENSG00000205636.3","ENSG00000235298.1","ENSG00000234323.10","ENSG00000233033.1","ENSG00000288597.1")

table(lnc_bruno %in% rownames(res_padj05))
table(nchar(lnc_bruno))
table(gsub("\\..*","", lnc_bruno) %in% gsub("\\..*","", rownames(res_padj05)))
table(substr(lnc_bruno,1,15) %in% substr(rownames(res_padj05),1,15))

lnc_bruno[gsub("\\..*","", lnc_bruno) %in% gsub("\\..*","", rownames(res_padj05)) == FALSE]
lnc_bruno[substr(lnc_bruno,1,15) %in% substr(rownames(res_padj05),1,15) == FALSE]

res_padj05_lncRNASelect <- res_padj05[gsub("\\..*","", rownames(res_padj05)) %in% 
                                        gsub("\\..*","", lnc_bruno),]
write.csv(res_padj05_lncRNASelect, file = "differentialExpression_lncRNASelect_TPNT_res05.csv")

################################################################################
library(ggplot2)
library(ggrepel)

data_vulc <- as.data.frame(res_padj05) 

data_vulc <- data_vulc[!is.na(data_vulc$padj),]
data_vulc <- data_vulc[data_vulc$padj < 0.05, ]

data_vulc$diffexpressed <- "NO"
data_vulc$diffexpressed[data_vulc$log2FoldChange < (-2)] <- "DOWN"
data_vulc$diffexpressed[data_vulc$log2FoldChange > (2)] <- "UP"

data_vulc$delabel <- NA
#idx <- match(gsub("\\..*","", lnc_bruno), gsub("\\..*","", rownames(data_vulc)))
#idx <- idx[!is.na(idx)]
#data_vulc$delabel[idx] <- rownames(data_vulc[idx,])

mycolors <- c("#00AFBB","#bb0c00","grey")
names(mycolors) <- c("DOWN","UP","NO")

graph <- ggplot(data_vulc, aes(x=log2FoldChange, y= -log10(padj),
                               col=diffexpressed, label=delabel))+
  geom_point(size=0.5)+
  labs(color = 'DEGs', x=expression("log"[2]*"FC"), y=expression("-log"[10]*"padj"))+
  coord_cartesian(ylim = c(0,100), xlim = c(-5, 5)) +
  #theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values = mycolors,
                     labels=c(paste("DOWN\n (",
                                    length(data_vulc$diffexpressed[data_vulc$diffexpressed=="DOWN"]),
                                    ")", sep = ""),
                              paste("NO\n (", length(data_vulc$diffexpressed[data_vulc$diffexpressed=="NO"]),
                                    ")", sep = ""),
                              paste("UP\n (",
                                    length(data_vulc$diffexpressed[data_vulc$diffexpressed=="UP"]),
                                    ")", sep = "")))+
  ggtitle("lncRNA All")+
  guides(color = guide_legend(label.position = "bottom",
                              title.position="top",
                              title.hjust=0.5,
                              override.aes=list(size=3,
                                                shape=21,
                                                fill=c("#00AFBB","grey","#bb0c00"))))+
  theme(
    plot.title=element_text(hjust=0.5),
    legend.position="bottom",
    legend.margin = margin(-5, -5, 5, -5),
    legend.spacing.y = unit(0, 'cm'), 
    title = element_text(size=10), 
    axis.title.y = element_text(size=10), 
    axis.title.x = element_text(size=10), 
    axis.text.y=element_text(size=8), 
    axis.text.x=element_text(size=8),
    axis.ticks = element_line(linetype=1, color="grey"),
    legend.text=element_text(size=10),
    legend.title=element_text(size=10 ),
    legend.key.size = unit(0,"mm"),
    legend.background = element_rect(fill=NULL, color=NULL),
    axis.line = element_blank(),
    panel.grid.major.x = element_line(linetype=111, color="grey80", size = 0.4),
    panel.grid.major = element_line(linetype=3, color="grey", size = 0.2),
    panel.background = element_rect(fill = "grey98", colour = "grey50"),
    panel.border = element_rect(colour = "grey", fill=NA, size=1))

ggsave("vulcanoplot.pdf",graph, width = 7, height = 5)
