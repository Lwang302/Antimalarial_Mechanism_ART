library(ggplot2)
library(readxl)
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(DMwR2)
library(patchwork)
library(impute)
library(DMwR2)


data = read_excel("inputfile",sheet = 1)
data<-data.frame(data)
head(data)

index<-duplicated(data$Accession)
print(table(index))
data<-data[!index,]
rownames(data)<-data$Accession
data<-data[,-1]
head(data)

meta<-data.frame(type=colnames(data),group=c(rep("case",3),rep("control",3)))
print(meta)

tmp_data_log<-log2(data+1)
samps<-factor(meta$group)
print(samps)
design <- model.matrix(~0+factor(samps))
colnames(design) <- levels(factor(samps))
rownames(design) <- meta$type
print(design)

contrast.matrix <- makeContrasts(contrasts="case-control",levels = design)

fit1 <- lmFit(tmp_data_log,design)                 #拟合模型
fit2 <- contrasts.fit(fit1, contrast.matrix) #统计检验
efit <- eBayes(fit2)                         #修正
tempOutput <- topTable(efit, coef="case-control", n=Inf)
tempOutput<-tempOutput[order(-tempOutput$logFC),]

head(tempOutput)
data["PF3D7_0526600",]

tempOutput$type = ifelse(tempOutput$adj.P.Val < 0.05 & abs(tempOutput$logFC) >= 1.5,
                             ifelse(tempOutput$logFC > 1.5 ,'Up','Down'),'Stable')
table(tempOutput$type)
head(tempOutput)
plot.volc <- ggplot(tempOutput) +
  geom_point(aes(x = logFC , y = -log10(adj.P.Val), color = type),alpha=0.6, size=0.5) +
  ggtitle(paste0("AHA ","vs control diff genes")) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value")+
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    #panel.grid.major = element_blank(), # get rid of major grid
    #panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    #legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ) +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  scale_color_manual(values=c("#546de5","#DDDDDD","#DDDDDD"))+
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.8,alpha=0.8)+
  geom_vline(xintercept = c(-1.5,1.5),lty=4,lwd=0.8,alpha=0.8)+
  scale_x_continuous(limits = c(-18,18))


print(plot.volc)


data_dn<-tempOutput[tempOutput$type=="Down",]
dim(data_dn)
data_dn_heatmap<-data[rownames(data_dn),]

head(data_dn_heatmap)
colnames(data_dn_heatmap)
data_dn_heatmap<-data_dn_heatmap[,c("DMSO..1","DMSO..2","DMSO..3","ATS_AHA..1","ATS_AHA..2","ATS_AHA..3")]

pheatmap(data_dn_heatmap,scale = "row",show_rownames = F,color=viridis(20),cluster_cols = F)


library(org.Pf.plasmo.db)
library(clusterProfiler)

dn_go <- enrichGO(gene = rownames(data_dn_heatmap),
                  keyType = "SYMBOL",
                  OrgDb = org.Pf.plasmo.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.5,
                  qvalueCutoff  = 0.5)


dotplot(dn_go,showCategory = 10)




