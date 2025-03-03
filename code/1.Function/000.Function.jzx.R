#### 2024/9/14
ggDBC_GSVA <- function(group=df_sample$group,group1='group1',group2='group2',gsva.res=gsva.res,output){
  ## limma差异通路分析
  
  #BiocManager::install('limma')
  library(limma)
  library(stringr)
  library(ggplot2)
  library(ggthemes)
  library(ggprism)
  library(dplyr)
  # 设置或导入分组
  design <- model.matrix(~0+group)
  colnames(design) = levels(factor(group))
  rownames(design) = colnames(gsva.res)
  design
  # Tunor VS Normal
  compareName <- paste0(group1, "-", group2)
  compare <- makeContrasts(contrasts = compareName, levels=design)
  fit <- lmFit(gsva.res, design)
  fit2 <- contrasts.fit(fit, compare)
  fit3 <- eBayes(fit2)
  Diff <- topTable(fit3, coef=1, number=200)
  head(Diff)
  
  ## 发散条形图绘制
  ## barplot
  dat_plot <- data.frame(id = row.names(Diff),
                         t = Diff$t)
  # 去掉"HALLMARK_"

  dat_plot$id <- str_replace(dat_plot$id , "HALLMARK_","")
  # 新增一列 根据t阈值分类
  dat_plot$threshold = factor(ifelse(dat_plot$t  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
  # 排序
  dat_plot <- dat_plot %>% arrange(t)
  # 变成因子类型
  dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
  # 绘制
  ##install.packages("ggthemes")
  #install.packages("ggprism")
  p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
    geom_col()+
    coord_flip() +
    scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
    geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
    xlab('') + 
    ylab(paste("t value of GSVA score", group1, "VS", group2)) + #注意坐标轴旋转了
    guides(fill=F)+ # 不显示图例
    theme_prism(border = T) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  p
  # 添加标签
  
  # 小于-1的数量
  low1 <- dat_plot %>% filter(t < -1) %>% nrow()
  # 小于0总数量
  low0 <- dat_plot %>% filter( t < 0) %>% nrow()
  # 小于1总数量
  high0 <- dat_plot %>% filter(t < 1) %>% nrow()
  # 总的柱子数量
  high1 <- nrow(dat_plot)
  
  # 依次从下到上添加标签
  p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                     hjust = 0,color = 'black',size=0.75) + # 小于-1的为黑色标签
    geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
              hjust = 0,color = 'grey',size=0.75) + # 灰色标签
    geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
              hjust = 1,color = 'grey',size=0.75) + # 灰色标签
    geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
              hjust = 1,color = 'black',size=0.75) # 大于1的为黑色标签
  
  p
  
  ggsave(filename = paste0(output,'.pdf'),p,width = 10,height  = 10)
  ggsave(filename = paste0(output,'.tiff'), plot = p, width = 8, height = 6, dpi = 600) 
}

gsva_yingbai <- function(expr=df_used,phenotype=df_sample,group.by='group',plotid='11',species = "human",group=df_sample$group,group1='Lac1',group2='Lac2',output){
  library(GSVA)
  library(msigdbr)
  #BiocManager::install('msigdbr')
  library(pheatmap)
  
  if(species=='human')
  {
    DBkeyset='org.Hs.eg.db'
    org='Homo sapiens'
    category='H'
    #library(org.Hs.eg.db)
    
  }else if(species=='mouse')
  {
    category='H'
    #library(org.Mm.eg.db)
    DBkeyset='org.Mm.eg.db'
    org='Mus musculus'
  }else if(species=='rat')
  {
    #library(org.Rn.eg.db)
    #BiocManager::install('org.Rn.eg.db')
    DBkeyset='org.Rn.eg.db'
    org='Mus musculus'
    category='H'
  }else{
    DBkeyset='org.Hs.eg.db'
    #library(org.Hs.eg.db)
    org='Homo sapiens'
    category='H'
  }
  
  expr<-as.matrix(expr)
  genesets <- msigdbr(species = org, category = category) 
  genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
  genesets <- split(genesets$gene_symbol, genesets$gs_name)
  
  gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
  gsva.df <- data.frame(gsva.res, check.names = F)
  gsva.df[1:3,1:3]
  write.table(gsva.df, paste0("fig",plotid,"a.",group.by,".gsva.hallmark.tsv"),sep='\t',row.names = T,quote = FALSE)
  
  
  ggheatmap_yingbio_2groups_for_gsva(gsva.df,phenotype,geneid='genesets',pheight=6,
                                     output=paste0('fig',plotid,'a.gsva.hallmark'),save.data=T)
  
  
  genesets <- msigdbr(species = org, category = "C2") 
  genesets <- subset(genesets, gs_subcat=="CP:KEGG", select = c("gs_name", "gene_symbol")) %>% as.data.frame()
  genesets <- split(genesets$gene_symbol, genesets$gs_name)
  
  gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
  gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
  write.table(gsva.df, paste0("fig",plotid,"b.",group.by,".gsva.kegg.tsv"),sep='\t',row.names = F,quote = FALSE)
  
  if(0)
  {
    gsva.df<-read.delim(paste0("fig",plotid,"b.",group.by,".gsva.kegg.tsv"),row.names = 1)
    gsva.df[1:3,1:3]
    
  }
  ggheatmap_yingbio_2groups_for_gsva(gsva.df,phenotype,geneid='genesets',pheight=16,
                                     output=paste0('fig',plotid,'b.gsva.kegg'),save.data=T)
  
  ggDBC_GSVA(group=phenotype$group,group1='Lac1',group2='Lac2',gsva.res=gsva.res,output)
}

out_tsv<-function(mx=mx,my_dir=my_dir)
{
  write.table(mx,file = paste0('./',my_dir,'/fig.result.txt'),sep = '\t',row.names = T,quote = F)
}

ggplot2_boxplot_Pvalue_stat <- function(df_plot, group = '', score = '', comparisons = my_comparisons, method = "wilcox.test", pd = FALSE, pw = 4, ph = 6) 
{
  library(ggplot2)
  library(ggpubr)
  
  df_plot$id <- factor(rep(1:(nrow(df_plot)/length(unique(df_plot[,group]))), length(unique(df_plot[,group]))))
  
  plota <- ggplot(df_plot, aes(x = .data[[group]], y = .data[[score]], color = .data[[group]])) + 
    stat_boxplot(geom = "errorbar", width = 0.35, lwd = 0.9) +
    geom_boxplot(width = 0.35, outlier.shape = NA, lwd = 0.75)
  
  if(method == "anova" || method == "kruskal.test"){
    plota <- plota + stat_compare_means(method = method, paired = pd, hjust = 0.5, vjust = 0.5, show.legend = T,label.y = max(df_plot[,score]) * 1.1)
  }else if(method == "friedman.test"){
    plota <- plota + stat_friedman_test(aes(wid = id),label.y = max(df_plot[,score]) * 1.1)
  }else{
    plota <- plota + stat_compare_means(comparisons = comparisons, method = method, paired = pd, label = 'p.format', hide.ns = FALSE, show.legend = FALSE)
  }
  
  
  plota <- plota + 
    theme(panel.background = element_blank(),  
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),  
          axis.line = element_line(colour = "black"),  
          axis.text.x = element_text(angle = 90, hjust = 0.6, vjust = 0.6, size = 10),  
          axis.title.x.bottom = element_blank(),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 20, hjust = 0.5, vjust = 0.5),
          legend.position = "none",
          panel.border = element_blank())  
  
  print('good good!!!')
  return(plota)
}

stat_jzx<-function(mx ,pd = F ,my_dir = my_dir)
{
  library(reshape2)
  library(car)
  library(nortest)
  library(PMCMRplus)
  library(dplyr)
  library(stats)
  
  #数据整理
  target<-'value'
  group<-'variable'
  mx_raw<-mx
  
  mean_list<-list()
  for (i in colnames(mx_raw)) {
    mean_list[[i]]<-mean(mx_raw[,i])
  }
  mean_df<-data.frame(mean_list)
  
  mx<-reshape2::melt(mx)
  mx <- mx %>%
    dplyr::select(-variable, variable)
  
  #定义参数
  nfactor<-length(colnames(mx))
  factor<-colnames(mx)
  
  plota <- '666'
  
  #正态分布检验
  st_tf  <- T
  st_result = list()
  for (i in colnames(mx_raw)) {
    st<-shapiro.test(as.matrix(mx_raw[[i]]))
    st_result[[i]] <- st$p.value
  }
  
  ano = ''
  for (i in 1:length(st_result)) {
    if (st_result[[i]] > 0.05) {
      status <- TRUE
    } else {
      status <- FALSE
      st_tf  <- F
    }
    ano <- paste0(ano,paste0(names(st_result)[i],':',status),' , ')
  }
  
  #是否为配对样本
  if(pd){
    annotation=paste0('正态检验:',ano)
    
    if(length(unique(mx[,2])) == 2){
      if(st_tf){
        #两配对样本t测验
        tt <- t.test(mx_raw[,1], mx_raw[,2], paired = TRUE)
        result <- data.frame(
          mean_df,
          method = 'Paired t.test',
          Pvalue = tt$p.value,
          FC     = (mean(mx_raw[,2])/mean(mx_raw[,1])),
          '正态检验'=ano
        )
        
      }else{
        #两配对样本秩和检验
        wt<-wilcox.test(mx_raw[,1], mx_raw[,2], paired = TRUE)
        
        result <- data.frame(
          mean_df,
          method = 'Paired wilcox.test',
          Pvalue = wt$p.value,
          FC     = (mean(mx_raw[,2])/mean(mx_raw[,1])),
          '正态检验'=ano
        )
        
      }
    }else{
      #多配对样本弗里德曼检验
      ft <- friedman.test(as.matrix(mx_raw))
      
      result <- data.frame(
        mean_df,
        method = 'friedman.test',
        Pvalue = ft$p.value,
        '正态检验'=ano
      )
    }
    
    for (i in 1:length(colnames(mx_raw))) {
      colnames(result)[i] <- paste0('mean_',colnames(mx_raw)[i])
    }
    out_tsv(result ,my_dir = my_dir)
    return(result)
  }
  
  
  if(nfactor == 2){
    formula <- as.formula(paste(target, "~", group))
    mx[,2]<-as.factor(mx[,2])
  }else{
    factor_list<-list()
    for (i in 2:length(factor)) {
      if(i != length(factor)){
        factor_temp<-paste0(factor[i],' * ')
        factor_list<-append(factor_list,factor_temp)
        mx[,i]<-as.factor(mx[,i])
      }else{
        factor_list<-append(factor_list,factor[i])
        mx[,i]<-as.factor(mx[,i])
      }
    }
    factor_used<-paste(unlist(factor_list), collapse = "")
    formula <- as.formula(paste(target, "~", factor_used))
  }
  
  #方差齐性检验
  lt<-leveneTest(formula, data = mx)
  
  if(st_tf && lt$`Pr(>F)`[1]>0.05){
    if(nfactor > 2 || length(unique(mx[,2])) > 2){
      #多组比较
      av<-aov(formula,data = mx)
      
      av_sr<-summary(av)
      av_sr<-data.frame(av_sr[[1]])
      
      result <- data.frame(
        mean_df,
        method = 'anova',
        Pvalue = av_sr[1,5],
        '正态检验'=ano,
        '方差齐性检验'=lt$`Pr(>F)`[1]>0.05
      )
      
    }else{
      #独立样本t测验
      tt <- t.test(formula, mx)
      
      result <- data.frame(
        mean_df,
        method = 't.test',
        Pvalue = tt$p.value,
        FC     = (mean(mx_raw[,2])/mean(mx_raw[,1])),
        '正态检验'=ano,
        '方差齐性检验'=lt$`Pr(>F)`[1]>0.05
      )
      
    }
  }else{
    if(nfactor > 2 || length(levels(mx[,2])) > 2){
      #多组比较
      kt<-kruskal.test(formula, data = mx)
      
      result <- data.frame(
        mean_df,
        method = 'kruskal.test',
        Pvalue = kt$p.value,
        '正态检验'=ano,
        '方差齐性检验'=lt$`Pr(>F)`[1]>0.05
      )
      
    }else{
      #两组比较
      if((nrow(mx_raw) * ncol(mx_raw)) > 6){
        wt<-wilcox.test(formula, data = mx)
        
        result <- data.frame(
          mean_df,
          method = 'wilcox.test',
          Pvalue = wt$p.value,
          FC     = (mean(mx_raw[,2])/mean(mx_raw[,1])),
          '正态检验'=ano,
          '方差齐性检验'=lt$`Pr(>F)`[1]>0.05
        )
        
      }else{
        #独立样本t测验
        tt <- t.test(formula, mx)
        
        result <- data.frame(
          mean_df,
          method = 't.test',
          Pvalue = tt$p.value,
          FC     = (mean(mx_raw[,2])/mean(mx_raw[,1])),
          '正态检验'=ano,
          '方差齐性检验'=lt$`Pr(>F)`[1]>0.05,
          Warning= 'The sample size is too small, and the P-value effect is low!!!'
        )
      }
    }
  }
  
  for (i in 1:length(colnames(mx_raw))) {
    colnames(result)[i] <- paste0('mean_',colnames(mx_raw)[i])
  }
  
  if(my_dir == ''){
    return(result)
  }else{
    out_tsv(result ,my_dir = my_dir)
    return(result)
  }
} 

box_plot_stat<-function(mx ,result = result ,my_dir = my_dir ,plot_only = F)
{
  library(dplyr)
  library(reshape2)
  
  #数据整理
  target<-'value'
  group<-'variable'
  
  mx<-reshape2::melt(mx)
  mx <- mx %>%
    dplyr::select(-variable, variable)
  
  my_comparisons <- list(levels(mx[,group]))
  pw<-(length(levels(mx[,group])))*2
  ph<-6
  
  method <- result$method
  my_method <- gsub('Paired ','',method)
  pd = F
  
  if(method == 'Paired t.test' || method == 'Paired wilcox.test' || method == 'friedman.test'){
    pd        = T
  }
  
  plota<-ggplot2_boxplot_Pvalue_stat(mx,group=group,score=target,comparisons=my_comparisons 
                                     ,method = my_method ,pd = pd ,pw = pw ,ph = ph)
  
  if(plot_only){
    return(plota)
  }else{
    ggsave(paste0("fig.box_plot.tiff"), plot = plota, path = paste0('./',my_dir), width = pw, height = ph,compression='lzw')
    ggsave(paste0("fig.box_plot.pdf"), plot = plota, path = paste0('./',my_dir), width = pw, height = ph)
    return(plota)
  }
}

#从fa文件中提取想要的基因生成新fa
fasta_subset<-function(file,target,out='output')
{
  library(seqinr)
  library(dplyr)
  
  fa<-read.fasta(file = file,forceDNAtolower = F)
  genes<-read.delim(target, header = F)
  genes<-genes[,1]
  
  # 处理列表元素名称
  names(fa) <- sub("\\.\\d+$", "", names(fa))  # 去掉.及其后的数字
  
  # 去重并保留后面的名称
  fa <- as.list(
    fa[!duplicated(names(fa), fromLast = TRUE)]
  )
  
  fa_out<-list()
  for (i in genes) {
    fa_out[i]<-fa[i]
  }
  write.fasta(sequences = fa_out,
              names=names(fa_out),
              file.out = paste0(out,".fa"))
}

#读取geo芯片数据
geo_parser <- function(GSE) {
  series_matrix_file = paste0(GSE,'_series_matrix.txt.gz')
  lines <- readLines(series_matrix_file)
  pd <- stringr::str_remove_all(lines, "^!|\"") %>% 
    str_split("\t", simplify = TRUE)
  k <- pd[,1] %>% str_starts("Sample_")
  pd <- pd[k, ]
  rownames(pd) <- pd[,1] %>% 
    str_remove("Sample_")
  pd <- t(pd[,-1])
  exp <- read.delim(series_matrix_file, 
                    comment.char = "!", 
                    row.names = 1) %>% 
    as.matrix()
  return(list(pd = pd, exp= exp))
}

seurat_cell_annotation_Neuron_PMID38981007<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
  #参考来源：PMID: 38981007
  Idents(pbmc) <- group.by
  nc<-length(levels(pbmc))
  ph<-max(nc/2,6)
  th=theme(axis.text.x = element_text(angle = 45,
                                      vjust = 1, hjust=1))
  markers_first = list(
    Radial_glial=c("PAX6", "SOX2",'MKI67','PCNA'),
    Olig_Precursor=c("CSPG4", "CNP", "SOX10"),
    Pre_Ast=c("APOE", "GJA1", "PHGDH",'THBS1','TNC'),
    Ast=c("GFAP", "CLU", "SOX9",'GLUD1'),
    MG=c("CD68", "AIF1","C1QA",'CST3','CX3CR1','P2RY12'),
    Neurons=c("NEUROD6", "DCX", "MAP2",'TUBB3','TBR1','SOX4')
  )
  if(species =='mouse')
  {
    library(stringr)
    for(i in 1:length(markers_first))
    {
      #i=1
      markers_first[[i]]<-str_to_title(markers_first[[i]])
    }
  }
  p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
  ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}

seurat_cell_annotation_Neuron_ExIn<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
  #参考来源：PMID: 32989152 ,PMID: 30385464,PMID：32101709
  Idents(pbmc) <- group.by
  nc<-length(levels(pbmc))
  ph<-max(nc/2,6)
  th=theme(axis.text.x = element_text(angle = 45,
                                      vjust = 1, hjust=1))
  markers_first = list(
     astrocytes = c("AQP4", "ADGRV1", "GPC5", "RYR3"),
     microglia = c("C3", "LRMDA", "DOCK8"),
     excitatory = c("CAMK2A", "CBLN2", "LDB2","SLC17A6","TAC1","CALB2"),
     inhibitory = c("GAD1","GAD2","LHFPL3", "PCDH15","SLC32A1","NOS1","GAL","VIP","NPY"),
     dopaminergic = c("DDC", "SLC6A3", "SLC18A2"),
     oligodendrocytes = c("MBP", "PLP1", "ST18")
   )
  if(species =='mouse')
  {
    library(stringr)
    for(i in 1:length(markers_first))
    {
      #i=1
      markers_first[[i]]<-str_to_title(markers_first[[i]])
    }
  }
  p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
  ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}

seurat_cell_annotation_Neuron_first<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
  #参考来源：https://zhuanlan.zhihu.com/p/685904134
  #PMID: 38245794
  Idents(pbmc) <- group.by
  nc<-length(levels(pbmc))
  ph<-max(nc/2,6)
  th=theme(axis.text.x = element_text(angle = 45,
                                      vjust = 1, hjust=1))
  markers_first = list(
    Fib=c('LUM'),
    Vascular=c('FLT1','CLDN5'),
    Ependymal=c('FOXJ1'),
    Choroid_plexus=c('TTR'),
    Neurons=c("DCX", "MAP2",'SNAP25'),
    OPCs=c('PDGFRA','VCAN','SOX10'),
    Olig =  c('PLP1','SLC44A1','MOG'),
    Ast =  c('GFAP','AQP4','SLC1A2'),
    Immune=c('PTPRC'),
    Microglia =  c("C1QC",'BHLHE41','CD68'),
    Tcells=c('CD44','THEMIS')
  )
  if(species =='mouse')
  {
    library(stringr)
    for(i in 1:length(markers_first))
    {
      #i=1
      markers_first[[i]]<-str_to_title(markers_first[[i]])
    }
  }
  p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
  ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}

seurat_cell_annotation_Neuron_sub<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human'){
  #参考来源：PMID：30385464 https://zhuanlan.zhihu.com/p/46201084 PMID：38981007
  #https://www.cellsignal.cn/pathways/neuronal-and-glial-cell-markers
  #http://xteam.xbio.top/CellMarker/
  Idents(pbmc) <- group.by
  nc<-length(levels(pbmc))
  ph<-max(nc/2,6)
  th=theme(axis.text.x = element_text(angle = 45,
                                      vjust = 1, hjust=1))
  markers_first = list(
     Neu_Stem=c('SOX2','MSI1','ID1','ID3','PROM1'),
     Neuroepithelial=c('HES1','NES'),
     Radial_glial=c('PAX6','HOPX'),
     OPCs=c('PDGFRA','VCAN'),
     Olig = c("OLIG1", "OLIG2",'CNP'),
     Astrocytes = c("GFAP", "ALDH1L1", "SLC1A2", "SLC1A3",'S100B'),
     Microglia = c("AIF1", "CD40LG", "CD68","SLC2A5"),
     Neurons=c('DCX'),
     Glu_Neurons = c('VGLUT1','VGLUT2','NEUROD2','NEUROD6','TBR1',"SLC17A6","SLC17A7","SLC17A8",
                     "GLS","GOT1","GRIN1", "GRIN2B"),
     GABA_Neurons = c("GAD1", "GAD2", "DLX1", "DLX5","ABAT", "GABBR1", "GABBR2"),
     DA_Neurons = c("DBH","FOXA2","LMX1B", "NR4A2","SLC6A2", "SLC6A3", "TH"),
     Serotonergic_Neurons = c("FEV","SLC6A4","TPH2"),
     Cholinergic_Neurons = c("ACHE", "CHAT", "SLC5A7", "SLC18A3")
   )
  if(species =='mouse')
  {
    library(stringr)
    for(i in 1:length(markers_first))
    {
      #i=1
      markers_first[[i]]<-str_to_title(markers_first[[i]])
    }
  }
  p <- DotPlot(pbmc, features = markers_first,assay='RNA' ,group.by = group.by )  +th
  ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
}

seurat_cell_annotation_MCs<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F,assay='RNA'){
  ####  https://pmc.ncbi.nlm.nih.gov/articles/PMC11647153/
  Idents(pbmc) <- group.by
  nc<-length(levels(pbmc))
  ph<-max(nc/2,6)
  th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
  markers_first = list(
    # 定义细胞类型及其对应的 marker
    MCs = c("KRT8", "TM4F1", "KRT19", "KRT18", "ITLN1"),
    Fibroblast = c("COL1A1", "COL3A1", "COL1A2",'SFRP4','IGFBP6'),
    EMCs = c("WFDC2", "CRCT1", "KRT13",'LCN2','SPRR3')
  )
  if(species =='mouse')
  {
    library(stringr)
    for(i in 1:length(markers_first))
    {
      #i=1
      markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
    }
    #markers_first[[21]][1]<-'Fcgr4'
    #markers_first[[21]][1]<-'Siglecf'
  }
  p <- DotPlot(pbmc, features = markers_first,assay=assay ,group.by = group.by )  +th
  ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.first.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
  if(plot.full)
  {
    markers<-as.character(unlist(markers_first))
    seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,ngene_cutoff=20000,dotplot.only=F,plot.heatmap=T)
  }
  
}

seurat_cell_annotation_eye<- function(pbmc,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F,assay='RNA'){
  ####  PMID：21091424  https://www.nature.com/articles/s41467-019-12780-8/figures/1
  Idents(pbmc) <- group.by
  nc<-length(levels(pbmc))
  ph<-max(nc/2,6)
  th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
  markers_first = list(
    # 定义细胞类型及其对应的 marker
    RPE = c('RPE65','RLBP1','RDH5','STRA6'),             #retinal pigmented epithelial cell
    vascular = c("ADAMTS9", "CD34", "CDH5",'RGS5'),      #vascular cells
    Rods = c("PDE6A", "PPEF2", "NR2E3"),
    Cones = c("GNAT2", "OPN1MW", "OPN1LW",'OPN1SW'),
    RGCs = c("NEFM", "SNCG", "SLC17A6"),                 #retinal ganglion cells
    Bipolar = c("CAMK2B", "GRM6", "TMEM215", "TRPM1"),   #bipolar cells
    Amacrine = c("GAD1", "C1QL2"),                       #amacrine cells
    Horizontal = c("LHX1", "ONECUT1", "ONECUT2"),        #horizontal cells
    Macroglia = c("CLU", "GLUL", "APOE"),
    microglia = c("C1QA", "TMEM119", "AIF1",'CD163')
  )
  if(species =='mouse')
  {
    library(stringr)
    for(i in 1:length(markers_first))
    {
      #i=1
      markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
    }
    #markers_first[[21]][1]<-'Fcgr4'
    #markers_first[[21]][1]<-'Siglecf'
  }
  p <- DotPlot(pbmc, features = markers_first,assay=assay ,group.by = group.by )  +th
  ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.first.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
  if(plot.full)
  {
    markers<-as.character(unlist(markers_first))
    seurat_markers_plot(pbmc,markers,group.by=group.by,plotid=plotid,ngene_cutoff=20000,dotplot.only=F,plot.heatmap=T)
  }
  
}


# BlockCounts Counts
#1           1   1275
#2           2   1398
ggplot2_count_barplot<-function(value,sample='BlockCounts',count='Counts',color='brown1',y_range=1500,ouput=paste0('fig.',plotid,'.barPlot'))
{
  library(ggplot2)
  col=sample
  plota <- ggplot(data = value, aes(fill=.data[[sample]],x = factor(.data[[sample]]), y = .data[[count]])) +
    geom_col(width = 0.5) +  # 使用 geom_col() 绘制柱状图，fill 设置柱子的填充颜色
    #geom_hline(yintercept = median_value, linetype = "dashed", color = "red") +  # 绘制虚线
    #geom_text(aes(label = round(count, 2)),size = 1.75, vjust = -0.5, color = "black") +  # 在每个柱子上显示值
    scale_y_continuous(expand = c(0, 0))+  # 让柱子紧贴 x 轴
    labs(title = sample)+
    xlab('')+
    ylab(count) +
    theme_minimal() +  # 使用简洁的主题，去掉背景网格线
    theme(
      axis.line = element_line(color = "black", size = 0.5),  # 设置 x 和 y 轴为实线
      plot.background = element_rect(fill = "white", color = "white"),  # 设置整个图的背景为白色
      panel.grid = element_blank(), # 去掉背景网格线
      axis.title.x = element_text(vjust = 2),  # 设置 x 轴标签移到顶部
      axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴标签，并使标签对齐
      plot.title = element_text(size = 10, hjust = 0.5)  # 调整标题的字体大小和水平居中
    )+
    coord_cartesian(ylim = c(0, y_range))  # 设置 y 轴的显示范围
  #coord_cartesian(ylim = c(0, max(value[,2]) + y_range))  # 设置 y 轴的显示范围
  ggsave(paste0(ouput,".tiff"), plot = plota, width = 8, height = 5,compression='lzw')
  ggsave(paste0(ouput,".pdf"), plot = plota, width = 8, height = 5)
}

#####  数据格式为ggplot2的标准格式，数据在一列，group在一列。
# my_comparisons <- list( c("Normal", "LGSC"), c("Normal", "HGSOC"), c("LGSC", "HGSOC"))
# ggplot2_boxplot_Pvalue(df_plot=df_box,group='group',score='Inflammatory_score',output='fig30.Inflammatory_score',comparisons=my_comparisons,pw=4,ph=6)
ggplot2_boxplot_Pvalue<-function(df_plot,group='group',score='score',output='fig8a.ALYREF_boxplot',
                                 comparisons=my_comparisons,method="wilcox.test",pd=F,pw=4,ph=6,P=T)
{
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  if(P)
  {label='p.format'}
  else
  {
   label='p.signif' 
      }
  plota<-NULL
  plota<-ggplot(df_plot, aes(x=.data[[group]], y=.data[[score]],color=.data[[group]])) + 
  stat_boxplot(geom = "errorbar", width=0.35,lwd = 0.9) +
  geom_boxplot(width=0.35,outlier.shape = NA,lwd = 0.75) +
  stat_compare_means(comparisons = comparisons,method = method, paired = pd,label = label,hide.ns = F,show.legend = FALSE)+
  #labs(title = colnames(df_plot[group]))+
  theme(panel.background = element_blank(),  # 去除背景
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank(),  # 去除次要网格线
        axis.line = element_line(colour = "black"),  # 显示xy轴
        axis.text.x = element_text(angle = 90, hjust = 0.6, vjust = 0.6,size = 10),  # 旋转x轴标签
        axis.title.x.bottom = element_blank(),
        axis.title.y = element_text(size = 20),
        plot.title=element_text(size=20,hjust=0.5,vjust=0.5),
        legend.position = "none",
        panel.border = element_blank())  # 去除边框
 print('good good!!!')
 ggsave(paste0(output,".tiff"), plot = plota, width = pw, height = ph,compression='lzw')
 ggsave(paste0(output,".pdf"), plot = plota, width = pw, height = ph)
}

#使用单细胞metadata画带标准误的柱状图
#x轴是celltype，柱子是grop时tf = T，反之为F
#ggplot2_barplot_paired_with_errorbar_jzx(metadata,x='celltype_sub',group='group',tf = T,sample='orig.ident',plotid=plotid,pw=8,ph=6)
ggplot2_barplot_paired_with_errorbar_jzx<-function(metadata,x='celltype_manual',group='group',sample='orig.ident',tf=F,plotid=plotid,pw=8,ph=6)
{
  if(tf)
  {
  df_data<-as.data.frame(table(metadata[,c(sample,x)]))
  df_group<-unique(metadata[,c(sample,group)])
  data<-merge(df_data,df_group,by=sample)
  }else{
  df_data<-as.data.frame(table(metadata[,c(sample,group)]))
  df_group<-unique(metadata[,c(sample,x)])
  data<-merge(df_data,df_group,by=sample)
  }
  ggplot2_barplot_paired_with_errorbar(data,value='Freq',x=group,y=x,y_title='Cell count',output=paste0("fig",plotid,".b.",x,'-',group,"_with_errorbar"),pw=8,ph=6)
  
  # 生成数据摘要
  data_summary <- summarySE(data, measurevar='Freq', groupvars=c(x,group))
  
  write.table(data_summary,file=paste0("fig",plotid,".a.",x,'_',group,"_with_errorbar.tsv"),quote=F,sep="\t", row.names=F)
  
  plota<-ggplot(data_summary, aes(x = .data[[x]], y = Freq, fill = .data[[group]],color = .data[[group]]), color = .data[[group]]) +
    geom_errorbar(aes(ymin = Freq - se, ymax = Freq + se), 
                  position = position_dodge(width = 0.83), 
                  width = 0.4,
                  show.legend = FALSE) + 
    geom_bar(stat = "identity", 
             position = position_dodge(width = 0.83), 
             width = 0.8) + 
    labs(y = NULL) +            
    theme_bw(base_size = 10) + 
    theme(
      axis.title = element_blank(),       # 隐藏坐标轴标题
      panel.border = element_blank(),     # 去掉包围线
      axis.line = element_line(),         # 显示 x 和 y 轴线
      panel.grid = element_blank(),       # 去掉网格线
      axis.text.x = element_text(colour = 'black',angle = 45,hjust = 1),  # 刻度字体颜色为黑色
      legend.title=element_blank()
    ) +
    scale_y_continuous(expand = c(0, 0),limits = c(0,max((data_summary$Freq) + max(data_summary$sd))*1.1))   # 柱子底部贴紧 x 轴
  ggsave(paste0("fig",plotid,".a.",x,'-',group,"_with_errorbar.tiff"), plot = plota, width = pw, height = ph,compression='lzw')
  ggsave(paste0("fig",plotid,".a.",x,'-',group,"_with_errorbar.pdf"), plot = plota, width = pw, height = ph)
}

cibersort_yingbio<-function(normalized_count,df_sample,sig_matrix='',plotid=6)
{
  library(CIBERSORT)
  source('/data/panzhong/genome_sequence/000cibersort_function.R')
  getwd()
  
  if(1)  ###使用cibersort函数，计算LM22中的22种细胞的肿瘤浸润得分
  {
    # 分别定义signature矩阵LM22和我的数据（演示）矩阵mixed_expr
    # results <- cibersort(sig_matrix = ref_custom, mixture_file = df_data,QN = F,perm=100)
    #normalized_count[is.na(normalized_count)] <- '0'
    if(sig_matrix=='')
    {
      sig_matrix=LM22
    }
    results <- cibersort(sig_matrix = sig_matrix, mixture_file = normalized_count,QN = F,perm=10)
    df_cibersort<-results[,1:(ncol(results)-3)]
    gc()
    
    class(results)
    dim(results)
    df_output<-as.data.frame(results)
    df_output<-cbind(rownames(df_output),df_output)
    colnames(df_output)[1] <-'sampleID'
    write.table(df_output,file=paste0('fig.',plotid,'.',projectid,'.cibersort.results.tsv'),row.names = F,sep='\t',quote = FALSE)
    
    # 理解一下results的结果
    # 你可以理解为返回了一个列名为细胞类型、行名为样本名的细胞浸润程度（占比）的矩阵
    # 此外result中还会多出三列：
    # P-value: 用来展示去卷积的结果在所有细胞类群中是否具有差异
    # Correlation:参考矩阵与输入矩阵的特征基因相关性
    # RMSE: Root mean squared error，参考矩阵与输入矩阵的特征基因标准差
  }
  #devtools::install_github("JanCoUnchained/ggunchained")
  
  if(1)   #### cibersort绘图
  {
    phenotype<-df_sample[,'group',drop=F]
    cibersort_plot(df_cibersort,phenotype,groupby='group',plotid=plotid)
  }
}

#快速ID转换
ID_convert<-function(df,fromType = "ENSEMBL",toType="SYMBOL",species='human')
{
  library(clusterProfiler)
  if(species=='human')
  {
    #BiocManager::install('org.Hs.eg.db')
    OrgDb='org.Hs.eg.db'
    library(org.Hs.eg.db)
  }else if(species=='mouse')
  {
    #BiocManager::install('org.Mm.eg.db')
    library(org.Mm.eg.db)
    OrgDb='org.Mm.eg.db'
  }else if(species=='rat')
  {
    library(org.Rn.eg.db)
    #BiocManager::install('org.Rn.eg.db')
    OrgDb='org.Rn.eg.db'
  }else{
    OrgDb='org.Hs.eg.db'
    library(org.Hs.eg.db)
  }
  df_used<-df
  df_used$gene<-rownames(df_used)
  genes_ensembl_id<-df_used$gene
  geneid_to_symbol<- bitr(genes_ensembl_id, fromType = fromType, toType=toType,OrgDb = OrgDb) #ENSEMBL  ENTREZID
  colnames(geneid_to_symbol)<-c('gene','gene_name')
  head(geneid_to_symbol)
  dim(geneid_to_symbol)
  
  geneid_to_symbol <- distinct(geneid_to_symbol)
  
  df_used <- merge(df_used,geneid_to_symbol,by='gene')
  
  duplicated_row_names <- duplicated(df_used$gene_name)
  df_used <- df_used[!duplicated_row_names, ]
  rownames(df_used) <- df_used$gene_name
  
  df_used$gene <- NULL
  df_used$gene_name <- NULL
  
  return(df_used)
}

#快速行名取交集两个数据框
df_common_merge<-function(df_1,df_2,plotid='1',projectid,rm_low=F,exp_tsv=F)
{
  if(rm_low)
  {
  df_1 <- df_1[rowMeans(df_1)>1,]
  df_2 <- df_2[rowMeans(df_2)>1,]
  }
  
  df_merge<-merge(df_1,df_2,by=0)
  
  if(exp_tsv)
  {
  df_row<-data.frame(df_merge[,1],row.names = df_merge[,1])
  df_1_m<-merge(df_row,df_1,by=0)
  rownames(df_1_m)<-df_1_m[,1]
  df_1_m<-df_1_m[,-1:-2]
  write.table(df_1_m, file= paste0("fig",plotid,"A.train.",projectid,".used.tsv"),sep="\t",row.names=T)
  
  df_2_m<-merge(df_row,df_2,by=0)
  rownames(df_2_m)<-df_2_m[,1]
  df_2_m<-df_2_m[,-1:-2]
  write.table(df_2_m, file= paste0("fig",plotid,"B.test.",projectid,".used.tsv"),sep="\t",row.names=T)
  }
  
  rownames(df_merge)<-df_merge[,1]
  df_merge[,1]<-NULL
  return(df_merge)
}

#快速去批次
combat_2df<-function(combi_data,test_phenotype,ref_phenotype,test_project,ref_project,plotid=10)
{
  test_phenotype$project=test_project
  ref_phenotype$project=ref_project
  combi_phenotype<-rbind(test_phenotype,ref_phenotype)
  output_table_with_rowname(combi_phenotype,name='sample',output=paste0('fig',plotid,'.',ref_project,'-',test_project,'-combat.phenotype'))
  
  if(1)  #########  绘制PCA图
  {
    sds<-apply(combi_data,1,sd)
    df_data_pca<-combi_data[sds>0,]
    ggplot2_pca_2factor(df_data_pca,combi_phenotype,output=paste0('fig',plotid,'a.before',".PCA-2factor"))
  }
  
  library(sva)
  model <- model.matrix(~group, data = combi_phenotype)
  df_expr_combat <- ComBat(dat = combi_data, batch = combi_phenotype$project, mod = model,ref.batch=ref_project)
  
  if(1)  #########  绘制PCA图
  {
    sds<-apply(df_expr_combat,1,sd)
    df_data_pca<-df_expr_combat[sds>0,]
    ggplot2_pca_2factor(df_data_pca,combi_phenotype,output=paste0('fig',plotid,'b.after',".PCA-2factor"))
  }
  
  df_expr_combat<-as.data.frame(df_expr_combat)
  write.table(df_expr_combat[,rownames(test_phenotype)], file= paste0("fig",plotid,".",test_project,".used.tsv"),sep="\t",row.names=T)
  output_table_with_rowname(df_expr_combat,name='gene',output=paste0('fig',plotid,'.',ref_project,'-',test_project,'-combat.data'))
  return(df_expr_combat)
}

#快速排序
rownames_match<-function(df1,df2,column='')
{
  if(column!='')
  {
    df2_order <- order(df2[column])
    df1_sorted <- df1[match(df2_order, df2[column]), ]
    return(df1_sorted)
  }else{
    df1_sorted <- df1[match(rownames(df2), rownames(df1)), ]
    return(df1_sorted)
  }
}

#快速韦恩图与热图
venn_and_2heat<-function(DEG_sig,gene_col='gene',target='',genes_use='',plotid=7,gs=rds_gene_set,
                         df_count,df_sample,comparison='',species='human',ph=2,pw=4)
{
  if(length(genes_use)==1)  ######  读入关注基因
  {
    gene_set <- readRDS(gs)
    names(gene_set)
    genes<-gene_set[[target]]
    
    if(species=='mouse')
    {    
     genes<-gene_human2mouse_vector(genes)
    }
  }else{
    genes=genes_use
    if(species=='mouse')
    {    
      genes<-gene_human2mouse_vector(genes)
    }
  }
  
  if(1)   ####绘制韦恩图
  {
    DEG_genes<-DEG_sig[[gene_col]]
    venn_genes<-ggplot2_venn2(DEG_genes,genes,name1='DEG',name2=target,output=paste0('fig',plotid,'.venn.DEGs_venn_',target))
    print(venn_genes)
  }
  
  if(1)   ######  绘制差异表达基因的热图（不区分PCG和LNC）
  {
    phenotype<-df_sample[,'group',drop=F]
    DEG_sig_filtered <- DEG_sig[venn_genes, ]
    df_count_heatmap1<-df_count[rownames(df_count) %in% DEG_genes,]
    df_count_heatmap2<-df_count[rownames(df_count) %in% venn_genes,]
    
    ggheatmap_yingbio_2groups(df_count_heatmap1,phenotype,geneid=gene_col,show_rowname = F,pheight=10,output=paste0('fig',plotid,".a.Heatmap.",comparison),save.data=T)
    ggheatmap_yingbio_2groups(df_count_heatmap2,phenotype,geneid=gene_col,show_rowname = T,pheight=ph,pw=pw,output=paste0('fig',plotid,".b.Heatmap.",comparison),save.data=T)
  }
  return(venn_genes)
}

#data是一个表达矩阵，横坐标基因名称纵坐标样本名，group是一个二分类向量，必须与data表达矩阵纵坐标一一对应
#参考：https://blog.csdn.net/weixin_55842556/article/details/128828895?spm=1001.2014.3001.8078#comments_35678317
SVM_RFE_2factor<-function(data=data,group,projectid,plotid='6b',nfold=10,nf=20)
{
  library(future)
  library(e1071)
  library(caret)
  plan("multicore", workers = 10)
  options(future.globals.maxSize = 8000 * 1024^2)
  
  svmRFE.wrap <- function(test.fold, X, ...) {
    # Wrapper to run svmRFE function while omitting a given test fold
    train.data = X[-test.fold, ]
    test.data  = X[test.fold, ]
    
    # Rank the features
    features.ranked = svmRFE(train.data, ...)
    
    return(list(feature.ids=features.ranked, train.data.ids=row.names(train.data), test.data.ids=row.names(test.data)))
  }
  
  svmRFE <- function(X, k=1, halve.above=5000) {
    # Feature selection with Multiple SVM Recursive Feature Elimination (RFE) algorithm
    n = ncol(X) - 1
    
    # Scale data up front so it doesn't have to be redone each pass
    cat('Scaling data...')
    X[, -1] = scale(X[, -1])
    cat('Done!\n')
    flush.console()
    
    pb = txtProgressBar(1, n, 1, style=3)
    
    i.surviving = 1:n
    i.ranked    = n
    ranked.list = vector(length=n)
    
    # Recurse through all the features
    while(length(i.surviving) > 0) {
      if(k > 1) {
        # Subsample to obtain multiple weights vectors (i.e. mSVM-RFE)            
        folds = rep(1:k, len=nrow(X))[sample(nrow(X))]
        folds = lapply(1:k, function(x) which(folds == x))
        
        # Obtain weights for each training set
        w = lapply(folds, getWeights, X[, c(1, 1+i.surviving)])
        w = do.call(rbind, w)
        
        # Normalize each weights vector
        w = t(apply(w, 1, function(x) x / sqrt(sum(x^2))))
        
        # Compute ranking criteria
        v    = w * w
        vbar = apply(v, 2, mean)
        vsd  = apply(v, 2, sd)
        c    = vbar / vsd
      } else {
        # Only do 1 pass (i.e. regular SVM-RFE)
        w = getWeights(NULL, X[, c(1, 1+i.surviving)])
        c = w * w
      }
      
      # Rank the features
      ranking = sort(c, index.return=T)$ix
      if(length(i.surviving) == 1) {
        ranking = 1
      }
      
      if(length(i.surviving) > halve.above) {
        # Cut features in half until less than halve.above
        nfeat = length(i.surviving)
        ncut  = round(nfeat / 2)
        n     = nfeat - ncut
        
        cat('Features halved from', nfeat, 'to', n, '\n')
        flush.console()
        
        pb = txtProgressBar(1, n, 1, style=3)
        
      } else ncut = 1
      
      # Update feature list
      ranked.list[i.ranked:(i.ranked-ncut+1)] = i.surviving[ranking[1:ncut]]
      i.ranked    = i.ranked - ncut
      i.surviving = i.surviving[-ranking[1:ncut]]
      
      setTxtProgressBar(pb, n-length(i.surviving))
      flush.console()
    }
    
    close(pb)
    
    return (ranked.list)
  }
  
  getWeights <- function(test.fold, X) {
    # Fit a linear SVM model and obtain feature weights
    train.data = X
    if(!is.null(test.fold)) train.data = X[-test.fold, ]
    
    svmModel = svm(train.data[, -1], train.data[, 1], cost=10, cachesize=500,
                   scale=F, type="C-classification", kernel="linear")
    
    t(svmModel$coefs) %*% svmModel$SV
  }
  
  WriteFeatures <- function(results, input, save=T, file='features_ranked.txt') {
    # Compile feature rankings across multiple folds
    featureID = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$ix
    avg.rank  = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$x
    feature.name = colnames(input[, -1])[featureID]
    features.ranked = data.frame(FeatureName=feature.name, FeatureID=featureID, AvgRank=avg.rank)
    if(save==T) {
      write.table(features.ranked, file=file, quote=F, row.names=F)
    } else {
      features.ranked
    }
  }
  
  FeatSweep.wrap <- function(i, results, input) {
    # Wrapper to estimate generalization error across all hold-out folds, for a given number of top features
    svm.list = lapply(results, function(x) tune(svm,
                                                train.x      = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                                                train.y      = input[x$train.data.ids, 1],
                                                validation.x = input[x$test.data.ids, 1+x$feature.ids[1:i]],
                                                validation.y = input[x$test.data.ids, 1],
                                                # Optimize SVM hyperparamters
                                                ranges       = tune(svm,
                                                                    train.x = input[x$train.data.ids, 1+x$feature.ids[1:i]],
                                                                    train.y = input[x$train.data.ids, 1],
                                                                    ranges  = list(gamma=2^(-12:0), cost=2^(-6:6)))$best.par,
                                                tunecontrol  = tune.control(sampling='fix'))$perf)
    
    error = mean(sapply(svm.list, function(x) x$error))
    return(list(svm.list=svm.list, error=error))
  }
  
  PlotErrors <- function(errors, errors2=NULL, no.info=0.5, ylim=range(c(errors, errors2), na.rm=T), xlab='Number of Features',  ylab='10x CV Error') {
    # Makes a plot of average generalization error vs. number of top features
    AddLine <- function(x, col='black') {
      lines(which(!is.na(errors)), na.omit(x), col=col)
      points(which.min(x), min(x, na.rm=T), col='red')
      text(which.min(x), min(x, na.rm=T), paste(which.min(x), '-', format(min(x, na.rm=T), dig=3)), pos=4, col='red', cex=0.75)
    }
    
    plot(errors, type='n', ylim=ylim, xlab=xlab, ylab=ylab)
    AddLine(errors)
    if(!is.null(errors2)) AddLine(errors2, 'gray30')
    abline(h=no.info, lty=3)
  }
  
  # 确保数据集是数据框格式，并且目标变量是因子类型
  data<-data.frame(data)
  data$group<-group
  data <- data %>% dplyr::select(group, everything())
  
  data <- as.data.frame((data))
  data$group <- factor(group)  # 确保 group 是因子类型
  
  nfold = nfold #10倍交叉验证
  nrows = nrow(data)
  folds = rep(1:nfold, len=nrows)[sample(nrows)]
  folds = lapply(1:nfold, function(x) which(folds == x))
  
  results = lapply(folds, svmRFE.wrap, data, k=10, halve.above=100)
  top.features = WriteFeatures(results, data, save=F)

  featsweep = lapply(1:nf, FeatSweep.wrap, results, data)
  save(featsweep,file = paste0('fig',plotid,'.',projectid,".svm_featsweep.Rdata"))
  
  no.info = min(prop.table(table(data[,1])))
  errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
  
  pdf(file=paste0('fig',plotid,'.',projectid,'.svm_rfe.pdf'), height = 8, width = 10)
  PlotErrors(errors, no.info=no.info)
  dev.off()
  tiff(file=paste0('fig',plotid,'.',projectid,'.svm_rfe.tiff'),width=5*300,height=5*300,res=300)
  PlotErrors(errors, no.info=no.info)
  dev.off()
  
  write.table(top.features,file=paste0('fig',plotid,'.svm_features_ALL.tsv'),sep = '\t',row.names = F)
}
  
#df_lasso是一个表达矩阵，横坐标基因名称纵坐标样本名，group是一个二分类向量，必须与df_lasso表达矩阵纵坐标一一对应
lasso_2factor<-function(df_lasso,group,plotid=6,projectid,nlambda=150,ROC_only=F)
{
  if(!ROC_only)
  {
    library(glmnet)
    #BiocManager::install('glmnet')
    
    fit <- glmnet(df_lasso, group, alpha = 1,family = "binomial", nlambda = nlambda) #glmnet是构建模型的
    
    pdf(paste0("fig",plotid,".",projectid,".fit.pdf"),height=5,width=5)
    plot(fit) 
    dev.off()
    tiff(paste0("fig",plotid,".",projectid,".fit.tiff"),height=5*300,width=5*300,res=300)
    plot(fit) 
    dev.off()
    
    cvfit=cv.glmnet(df_lasso, group, family = "binomial", type.measure = "deviance")  #交叉验证
    
    #绘图展示
    
    pdf(paste0("fig",plotid,".",projectid,".cvfit.pdf"),height=5,width=5)
    plot(cvfit)  
    dev.off()
    tiff(paste0("fig",plotid,".",projectid,".cvfit.tiff"),height=5*300,width=5*300,res=300)
    plot(cvfit)  
    dev.off()
    
    save(cvfit,file = paste0("fig",plotid,".",projectid,".cvfit.Rdata"))
    
    feature_all <- as.data.frame(as.matrix(coef(cvfit, s = cvfit$lambda.min)))
    colnames(feature_all) <- "coff"
    library(dplyr)
    feature_opt <-  feature_all %>% filter(abs(coff) > 0)
    feature_opt <- cbind(row.names(feature_opt), feature_opt)
    feature_opt<-feature_opt[-1,]
    colnames(feature_opt)[1]<-'gene'
    write.table(feature_opt, file= paste0("fig",plotid,".",projectid,".feature_opt.tsv"),sep="\t",row.names=F)
  }
  
  if(ROC_only)
  {
    library(dplyr)
    load(paste0("fig",plotid,".",projectid,".cvfit.Rdata"))
    feature_all <- as.data.frame(as.matrix(coef(cvfit, s = cvfit$lambda.min)))
    colnames(feature_all) <- "coff"
    feature_opt <-  feature_all %>% filter(abs(coff) > 0)
  }
  
  lassogenes<-feature_opt$gene
  
  if(0)  ## 构建列线图
  {
    library(rms)
    df_lasso_lxt <- as.data.frame(df_lasso)
    # 只保留feature_opt中的基因
    df_lasso_lxt <- df_lasso_lxt[, lassogenes, drop = FALSE]
    df_lasso_lxt$group <- as.factor(group)
    
    # 创建一个数据描述符（Data Descriptor）
    dd <- datadist(df_lasso_lxt)
    options(datadist = dd)
    
    # 拟合逻辑回归模型
    fit_lrm <- lrm(group ~ ., data = df_lasso_lxt)
    summary(fit_lrm) 
    
    # 绘制列线图
    nom <- nomogram(fit_lrm, funlabel = "Probability")
    pdf(paste0("fig",plotid,"a.",projectid,".nomogram.pdf"), height = 6, width = 8)
    plot(nom)
    dev.off()
    
    tiff(paste0("fig",plotid,"a.",projectid,".nomogram.tiff"), height = 6*300, width = 8*300, res = 300)
    plot(nom)
    dev.off()
  }
  
  if(1) #绘制ROC曲线
  {
    ## 计算riskscore
    library(dplyr)
    risk_score_table <- df_lasso[,lassogenes]
    for(each_sig_gene in colnames(risk_score_table)){
      risk_score_table[,each_sig_gene] <- risk_score_table[,each_sig_gene]*(feature_opt[each_sig_gene,2])
    }
    head(risk_score_table)
    head(rowSums(risk_score_table))
    head(exp(rowSums(risk_score_table)))
    risk_score_table <- cbind(risk_score_table, 'risk_score'=rowSums(risk_score_table))
    
    
    library(pROC)
    for(i in 1:length(lassogenes))
    {
      rocobj1=roc(group,risk_score_table[,lassogenes[i]],smooth=F)
      pdf(file=paste0("fig",plotid+1,".",projectid,".",lassogenes[i],".ROC.pdf"),width=5,height=5)
      plot(rocobj1,legacy.axes = TRUE,print.auc=TRUE, col="red", main = paste0(lassogenes[i]," ROC Curve"))
      dev.off()
      tiff(file=paste0("fig",plotid+1,".",projectid,".",lassogenes[i],".ROC.tiff"),width=5*300,height=5*300,res=300)
      plot(rocobj1,legacy.axes = TRUE,print.auc=TRUE, col="red", main = paste0(lassogenes[i]," ROC Curve"))
      dev.off()
    }
  }
}

#绘制ROC曲线
ROC_genes<-function(df_lasso,feature_all,lassogenes='',group,plotid,projectid)
{
  rownames(feature_all)<-feature_all$gene
  feature_opt <-  feature_all %>% filter(abs(coff) > 0)
  
  if(length(lassogenes)==1){
    lassogenes<-feature_opt$gene
  }

  ## 计算riskscore
  library(dplyr)
  risk_score_table <- df_lasso[,lassogenes]
  for(each_sig_gene in colnames(risk_score_table)){
    risk_score_table[,each_sig_gene] <- risk_score_table[,each_sig_gene]*(feature_opt[each_sig_gene,2])
  }
  print(head(risk_score_table))

  risk_score_table <- cbind(risk_score_table, 'risk_score'=rowSums(risk_score_table))
  
  if(1)  ##总体 ROC 验证曲线绘制
  {
    library(pROC)
    rocobj1=roc(group,risk_score_table[,"risk_score"])
    pdf(file=paste0("fig",plotid,"A.",projectid,".ALL.ROC.pdf"),width=5,height=5)
    plot(rocobj1,legacy.axes = TRUE,print.auc=TRUE, col="red", main = paste0("ROC Curve"))
    dev.off()
    tiff(file=paste0("fig",plotid,"A.",projectid,".ALL.ROC.tiff"),width=5*300,height=5*300,res=300)
    plot(rocobj1,legacy.axes = TRUE,print.auc=TRUE, col="red", main = paste0("ROC Curve"))
    dev.off()
  }
  
  for(i in 1:length(lassogenes))##基因 ROC 验证曲线绘制
  {
    rocobj1=roc(group,risk_score_table[,lassogenes[i]],smooth=F)
    pdf(file=paste0("fig",plotid,"B.",projectid,".",lassogenes[i],".ROC.pdf"),width=5,height=5)
    plot(rocobj1,legacy.axes = TRUE,print.auc=TRUE, col="red", main = paste0(lassogenes[i]," ROC Curve"))
    dev.off()
    tiff(file=paste0("fig",plotid,"B.",projectid,".",lassogenes[i],".ROC.tiff"),width=5*300,height=5*300,res=300)
    plot(rocobj1,legacy.axes = TRUE,print.auc=TRUE, col="red", main = paste0(lassogenes[i]," ROC Curve"))
    dev.off()
  }
  
}

#绘制列线图
nom_2factor<-function(df_lasso,lassogenes,group,projectid,plotid)
{
  library(rms)
  df_lasso_lxt <- as.data.frame(df_lasso)
  # 只保留feature_opt中的基因
  df_lasso_lxt <- df_lasso_lxt[, lassogenes, drop = FALSE]
  df_lasso_lxt$group <- as.factor(group)
  
  # 创建一个数据描述符（Data Descriptor）
  dd <- datadist(df_lasso_lxt)
  options(datadist = dd)
  
  # 拟合逻辑回归模型
  fit_lrm <- lrm(group ~ ., data = df_lasso_lxt)
  #summary(fit_lrm) 
  
  # 绘制列线图
  nom <- rms::nomogram(fit_lrm, funlabel = "Probability")
  pdf(paste0("fig",plotid,"a.",projectid,".nomogram.pdf"), height = 6, width = 8)
  plot(nom)
  dev.off()
  
  tiff(paste0("fig",plotid,"a.",projectid,".nomogram.tiff"), height = 6*300, width = 8*300, res = 300)
  plot(nom)
  dev.off()
}

#绘制DCA曲线
fit_2factor_DCA<-function(fit=fit,Train_set=Train_set,Type=Train_data$Type,projectid='',plotid='15')
{
  # 计算预测概率 pred_prob
  pred_prob <- predict(fit, newx = as.matrix(Train_set), type = "response")
  #pred_prob <- predict(fit, newx = as.matrix(Train_set), type = "response")
  
  dcadata<-data.frame(Type,pred_prob[,1])
  colnames(dcadata)<-c('Type','pred_prob')
  
  # 执行决策曲线分析
  dca_result <- decision_curve(Type ~ pred_prob, 
                               data = dcadata,
                               study.design = c("case-control"))
  
  # 绘制决策曲线
  pdf(paste0("fig",plotid,"c.",projectid,".DCA.pdf"),height=5,width=6)
  plot_decision_curve(dca_result, curve.names = projectid,
                      cost.benefit.axis = T, # 是否需要损失：获益比 轴
                      confidence.intervals = "none" ,# 不画可信区间
                      legend.position = c("bottomleft")
  )
  dev.off()
  
  tiff(paste0("fig",plotid,"c.",projectid,".DCA.tiff"),height=5*300,width=6*300,compression = 'lzw',res=300)
  plot_decision_curve(dca_result, curve.names = projectid,
                      cost.benefit.axis = T, # 是否需要损失：获益比 轴
                      confidence.intervals = "none" ,# 不画可信区间
                      legend.position = c("bottomleft")
  )
  dev.off()
}

#绘制校准曲线
fit_2factor_jiaozhun<-function(fit,data=Train_set,Type,projectid,plotid)
{
  data<-as.data.frame(data)
  pred_prob <- predict(fit, newdata = data, type = "response")
  #pred_prob <- predict(fit, newx = as.matrix(Train_set), type = "response")
  pred_prob<-pred_prob[,1]
  #建立模型公式
  form.bestglm<-as.formula(Type~pred_prob)
  #form.all<-as.formula(Type~.)
  #打包
  library(rms)
  dd=datadist(data)
  options(datadist=dd)
  #建立logistic模型
  fit.glm<- lrm(formula=form.bestglm,data=data,x=TRUE,y=TRUE) 
  cal.glm<-calibrate(fit.glm,method = "boot",B=10000)
  
  #par(mar=c(6,6,6,6))
  pdf(paste0("fig", plotid, "b.", projectid, ".cal.pdf"), height = 8, width = 8)
  plot(cal.glm,
       xlim = c(0,1),
       ylim = c(0,1),
       xlab = "Predicted Probability",
       ylab="Observed  Probability",
       xaxs = "i",
       yaxs = "i",
       legend =FALSE,
       subtitles = FALSE,#不显示副标题
       cex=1.5,
       cex.axis=1.5,
       cex.lab=1.5)
  abline(0,1,col="blue",lty=2,lwd=2) #绘制参考线
  #调用cal.glm中的“calibrated.orig",即未校正的曲线，设置为实线，红色
  lines(cal.glm[,c("predy","calibrated.orig")],type="l",lwd=2,col="red")
  #调用cal.glm中的"calibrated.corrected",即校准的曲线，实线绿色
  lines(cal.glm[,c("predy","calibrated.corrected")],type="l",lwd=2,col="green")
  legend(x=0.55,y=0.35,
         legend=c("Ideal","Apparent","Bias-corrected"),
         lty = c(2,1,1),
         lwd = c(2,1,1),
         col = c("blue","red","green"),
         bty="n",
         cex=1.5)
  dev.off()
  
  #par(mar=c(6,6,6,6))
  tiff(paste0("fig", plotid, "b.", projectid, ".cal.tiff"), height = 8 * 300, width = 8 * 300, compression = 'lzw', res = 300)
  plot(cal.glm,
       xlim = c(0,1),
       ylim = c(0,1),
       xlab = "Predicted Probability",
       ylab="Observed  Probability",
       xaxs = "i",
       yaxs = "i",
       legend =FALSE,
       subtitles = FALSE,#不显示副标题
       cex=1.5,
       cex.axis=1.5,
       cex.lab=1.5)
  abline(0,1,col="blue",lty=2,lwd=2) #绘制参考线
  #调用cal.glm中的“calibrated.orig",即未校正的曲线，设置为实线，红色
  lines(cal.glm[,c("predy","calibrated.orig")],type="l",lwd=2,col="red")
  #调用cal.glm中的"calibrated.corrected",即校准的曲线，实线绿色
  lines(cal.glm[,c("predy","calibrated.corrected")],type="l",lwd=2,col="green")
  legend(x=0.55,y=0.35,
         legend=c("Ideal","Apparent","Bias-corrected"),
         lty = c(2,1,1),
         lwd = c(2,1,1),
         col = c("blue","red","green"),
         bty="n",
         cex=1.5)
  dev.off()
}

# 创建一个函数来绘制临床影响曲线
fit_2factor_CIC <- function(fit, data, Type, projectid = 'GSE160306', plotid = '15', population_prevalence = 0.4) 
{
  # 导入 rmda 包
  library(rmda)
  
  # 计算预测概率 pred_prob
  pred_prob <- predict(fit, newdata = as.data.frame(data), type = "response")
  
  # 将预测概率和实际类别合并成数据框
  dcadata <- data.frame(Type, pred_prob)
  colnames(dcadata) <- c('Type', 'pred_prob')
  
  # 执行决策曲线分析 (DCA) 
  # 使用 logistic 回归的模型和不同的阈值来计算决策曲线
  simple <- decision_curve(
    Type ~ pred_prob,                         # Type 是目标变量，pred_prob 是预测概率
    data = dcadata,
    family = binomial(link = 'logit'),        # 使用逻辑回归
    thresholds = seq(0, 1, by = 0.01),       # 设置阈值范围
    confidence.intervals = 0.95,              # 设置置信区间
    study.design = 'case-control',           # 研究设计类型：case-control 或 cohort
    population.prevalence = population_prevalence  # 设置患病率（这个值可以根据具体的研究领域调整）
  )
  
  # 绘制临床影响曲线
  pdf(paste0("fig", plotid, "d.", projectid, ".CIC.pdf"), height = 5, width = 6)
  plot_clinical_impact(
    simple,
    population.size = 1000,                # 样本量
    cost.benefit.axis = TRUE,               # 显示损失/收益坐标轴
    n.cost.benefits = 8,                    # 设置损失/收益的数量
    col = c('red', 'blue'),                 # 设置曲线颜色
    confidence.intervals = F             # 不显示置信区间
  )
  dev.off()
  
  # 保存为 tiff 格式
  tiff(paste0("fig", plotid, "d.", projectid, ".CIC.tiff"), height = 5 * 300, width = 6 * 300, compression = 'lzw', res = 300)
  plot_clinical_impact(
    simple,
    population.size = 1000,
    cost.benefit.axis = TRUE,
    n.cost.benefits = 8,
    col = c('red', 'blue'),
    confidence.intervals = F
  )
  dev.off()
}

seurat_integrate_null_sketch <- function(pbmc,ncells=50000,group.by='orig.ident',species = "human",dim.use=50,nfeatures=2000,res=0.5,plotid='02',
                                         mt.pattern="^MT-",mt.list=NULL,mt.cutoff=5,nf.low=500,nf.high=6000) 
{
  library(BPCells)
  Idents(pbmc) <- group.by
  options(future.globals.maxSize = 1e9*80)

  pbmc <- NormalizeData(pbmc)
  pbmc <- FindVariableFeatures(pbmc)
  pbmc <- SketchData(
    object = pbmc,
    ncells = 50000,
    method = "LeverageScore",
    sketched.assay = "sketch"
  )
  pbmc
  
  DefaultAssay(pbmc) <- "sketch"
  pbmc <- FindVariableFeatures(pbmc)
  pbmc <- ScaleData(pbmc, verbose = FALSE)
  pbmc <- RunPCA(pbmc,npcs = 50,pc.genes=VariableFeatures(object = pbmc))
  pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:dim.use)
  pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:dim.use)
  pbmc <- FindClusters(pbmc, resolution = res)
  
  pbmc <- ProjectData(
    object = pbmc,
    assay = "RNA",
    full.reduction = "pca.full",
    sketched.assay = "sketch",
    sketched.reduction = "pca",
    umap.model = "umap",
    dims = 1:dim.use,
    refdata = list(cluster_full = "seurat_clusters")
  )
  # now that we have projected the full dataset, switch back to analyzing all cells
  DefaultAssay(pbmc) <- "RNA"
  
  #输出特征方差图
  top10 <- head(x = VariableFeatures(object = pbmc), 10)
  plot1 <- VariableFeaturePlot(object = pbmc)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plota<-CombinePlots(plots = list(plot1, plot2))
  ggsave(paste0("fig",plotid,"a.VariableFeatures.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"a.VariableFeatures.pdf"), plot = plota, width = 12, height = 6)
  
  plota<-NULL
  plota<-VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
  ggsave(paste0("fig",plotid,"b.pcaGene.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
  ggsave(paste0("fig",plotid,"b.pcaGene.pdf"), plot = plota, width = 10, height = 10)
  plota<-NULL
  plota<-ElbowPlot(pbmc, ndims = 50)
  ggsave(paste0("fig",plotid,"c.pca50.tiff"), plot = plota, width = 10, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"c.pca50.pdf"), plot = plota, width = 10, height = 6)
  
  #绘制主成分分析图形
  plota<-NULL
  plota<-DimPlot(object = pbmc, reduction = "pca")
  ggsave(paste0("fig",plotid,"d.PCA.tiff"), plot = plota, width = 10, height = 10,compression='lzw')
  ggsave(paste0("fig",plotid,"d.PCA.pdf"), plot = plota, width = 10, height = 10)
  
  plota<-NULL
  plota<-DimPlot(pbmc, reduction = "umap", group.by = 'seurat_clusters',pt.size = 0.5,label=T)
  ggsave(paste0("fig",plotid,"e0.umap.cluster.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
  ggsave(paste0("fig",plotid,"e0.umap.cluster.pdf"), plot = plota, width = 8, height = 8)
  
  meta<-pbmc@meta.data
  if('group' %in% colnames(meta))
  {
    plota<-NULL
    plota<-DimPlot(pbmc, reduction = "umap", group.by = 'group',pt.size = 0.5,label=T)
    ggsave(paste0("fig",plotid,"e1.umap.group.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
    ggsave(paste0("fig",plotid,"e1.umap.group.pdf"), plot = plota, width = 8, height = 8)
    
    ngroup<-length(unique(meta$group))
    
    plotc<-NULL
    plotc <- DimPlot(pbmc, reduction = "umap", split.by = "group")
    ggsave(paste0("fig",plotid,"e2.umap.group.split.tiff"), plot = plotc, width = ngroup*6, height = 6,compression='lzw')
    ggsave(paste0("fig",plotid,"e2.umap.group.split.pdf"), plot = plotc, width = ngroup*6, height = 6)
  }
  plota<-NULL
  plota<-DimPlot(pbmc, reduction = "umap", group.by = 'orig.ident',pt.size = 0.5,label=T)
  ggsave(paste0("fig",plotid,"e3.umap.sample.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
  ggsave(paste0("fig",plotid,"e3.umap.sample.pdf"), plot = plota, width = 8, height = 8)
  
  nsample<-length(unique(meta$orig.ident))
  nh<-ceiling(nsample/4)
  
  plotc<-NULL
  plotc <- DimPlot(pbmc, reduction = "umap", split.by = "orig.ident",ncol=4)
  ggsave(paste0("fig",plotid,"e4.umap.sample.split.tiff"), plot = plotc, width = 12, height = nh*3,compression='lzw')
  ggsave(paste0("fig",plotid,"e4.umap.sample.split.pdf"), plot = plotc, width = 12, height = nh*3)
  plota<-NULL
  plota<-FeaturePlot(object = pbmc, features = c('percent.mt','percent.HB','percent.rb'), cols = c("green", "red"),ncol=3)
  ggsave(paste0("fig",plotid,"e5.umap.percentmt.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"e5.umap.percentmt.pdf"), plot = plota, width = 18, height = 6)
  
  if(0)
  {
    #细胞周期回归：上一步找到的高变基因，常常会包含一些细胞周期相关基因。
    #它们会导致细胞聚类发生一定的偏移，即相同类型的细胞在聚类时会因为细胞周期的不同而分开。
    ?CaseMatch
    cc.genes
    CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(scRNA1))
    #细胞周期评分
    g2m_genes = cc.genes$g2m.genes
    g2m_genes = CaseMatch(search = g2m_genes, match = rownames(scRNA1))
    s_genes = cc.genes$s.genes
    s_genes = CaseMatch(search = s_genes, match = rownames(scRNA1))
    scRNA1 <- CellCycleScoring(object=scRNA1,  g2m.features=g2m_genes,  s.features=s_genes)
    scRNAb <- ScaleData(scRNA1, vars.to.regress = c("S.Score", "G2M.Score"))
    #     然后重新PCA、UMAP等计算。
    if(species=='human')
    {
      pbmc <- CellCycleScoring(
        object = pbmc,
        g2m.features = cc.genes$g2m.genes,
        s.features = cc.genes$s.genes
      )
    }else{
      library(stringr)
      pbmc <- CellCycleScoring(
        object = pbmc,
        g2m.features = str_to_title(cc.genes$g2m.genes),
        s.features = str_to_title(cc.genes$s.genes)
      )
    }
    
    plota<-NULL
    plota<-VlnPlot(pbmc, features = c("S.Score","G2M.Score"),group.by="orig.ident")
    ggsave(paste0("fig",plotid,"f1.umap.cellCycle.tiff"), plot = plota, width = 8, height = 8,compression='lzw')
    ggsave(paste0("fig",plotid,"f1.umap.cellCycle.pdf"), plot = plota, width = 8, height = 8)
    #head(pbmc@meta.data)
    #table(pbmc$Phase)
    plota<-NULL
    plota<-DimPlot(pbmc,group.by='Phase',split.by = "Phase")
    ggsave(paste0("fig",plotid,"f2.umap.Phase.tiff"), plot = plota, width = 18, height = 6,compression='lzw')
    ggsave(paste0("fig",plotid,"f2.umap.Phase.pdf"), plot = plota, width = 18, height = 6)
    DefaultAssay(pbmc) <- "RNA"
  }
  
  DefaultAssay(pbmc) <- "RNA"
  all.genes<-rownames(pbmc)
  pbmc<-ScaleData(pbmc,features=all.genes)
  
  return(pbmc)
}




