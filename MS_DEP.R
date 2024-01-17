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

data_list<-readRDS(input)


diff_protein<-function(sample_names){
  tmp_data<-data_list[[sample_names]]
  print(colnames(tmp_data))

  meta<-data.frame(type=colnames(tmp_data),group=c(rep("case",3),rep("control",3)))
  #print(meta)
  tmp_data_log<-log2(tmp_data+1)
  samps<-factor(meta$group)
  #print(samps)
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
  return(tempOutput)
}

names(data_list)
diff_pro_list<-list()
for (i in names(data_list)) {
  print(i)
  diff_pro_list[[i]]<-diff_protein(sample_names=i)
}


plot.volc<-function(data_in,title){
  #data statistics
  print(title)
  #diff_statistics<-subset(data_in,P.Value<=0.05)
  #diff_statistics<-subset(diff_statistics,logFC >= 0.263)
  #print(dim(diff_statistics))

  diff_tmp<-subset(data_in,logFC >= -0.263)
  diff_tmp$type = ifelse(diff_tmp$P.Value < 0.05 & abs(diff_tmp$logFC) >= 0.263,
                         ifelse(diff_tmp$logFC >= 0.263 ,'Up','Down'),'Stable')
  diff_tmp<-na.omit(diff_tmp)
  stat<-table(diff_tmp$type)
  print(stat)
  a<-stat[["Stable"]]
  b<-stat[["Up"]]
  plot <- ggplot(diff_tmp) +
    geom_point(aes(x = logFC, y = -log10(P.Value), color = type),alpha=0.6, size=1.5) +
    ggtitle(paste0(title,"vs Dmso ","stable=",a,",Up=",b)) +
    xlab("log2 fold change") +
    ylab("-log10 p-value") +
    scale_y_continuous(limits = c(0,9)) +scale_x_continuous(limits = c(-2,8))+
    theme(
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      #panel.grid.major = element_blank(), # get rid of major grid
      #panel.grid.minor = element_blank(), # get rid of minor grid
      legend.background = element_rect(fill = "transparent"), # get rid of legend bg
      #legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    ) +
    theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
    scale_color_manual(values=c("#DDDDDD","#ff4757"))+
    geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.8,alpha=0.8)+
    geom_vline(xintercept = c(0.263),lty=4,lwd=0.8,alpha=0.8)

  return(plot)
}


###############
str(diff_pro_list)
data_in<-diff_pro_list[["AP2_R_UV.1"]]

pro_pick<-function(data_in){
  diff_statistics<-subset(data_in,P.Value<=0.05)
  diff_statistics<-subset(diff_statistics,logFC >= 0.263)
  pro_name<-rownames(diff_statistics)
  print(length(pro_name))
  return(pro_name)
}

pro_list<-list()

for (i in names(diff_pro_list)) {
  print(i)
  pro_list[[i]]<-pro_pick(data_in=diff_pro_list[[i]])
}

str(pro_list)

pro_list_a<-unique(unlist(pro_list))
length(pro_list_a)

sub_data<-data[rownames(data) %in% pro_list_a,]


pheatmap(sub_data,scale = "row",show_rownames = F,cluster_cols = F,cluster_rows = F,color=plasma(20))


