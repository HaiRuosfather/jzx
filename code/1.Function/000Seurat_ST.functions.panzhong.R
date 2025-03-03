###  pbmc<-seurat_vertion4to5(pbmc)
seurat_vertion4to5<-function(pbmc)
  {
    
    pbmc[["RNA5"]] <- as(object = pbmc[["RNA"]], Class = "Assay5")
    DefaultAssay(pbmc)<-'RNA5'
    pbmc[['RNA']]=NULL  
    pbmc[["RNA"]] <- as(object = pbmc[["RNA5"]], Class = "Assay5")
    DefaultAssay(pbmc)<-'RNA'
    pbmc[['RNA5']]=NULL
    return(pbmc)
}

seurat_ST_add_celltype<-function(obj_st) 
{
  metadata<-obj_st[["predictions"]]@data
  class(metadata)
  metadata<-as.data.frame(t(metadata))
  dim(metadata)
  head(metadata)
  metadata$celltype<-'unknown'
  
  for(i in 1:nrow(metadata))
  {
    for(j in 1:(ncol(metadata)-2))
    {
      if(metadata[i,j]==metadata[i,'max'])
      {
        metadata[i,'celltype']<-colnames(metadata)[j]
      }
    }
  }
  obj_st@meta.data$celltype<-metadata$celltype
  return(obj_st)
}
seurat_ST_add_celltypes<-function(pbmc_st) 
{
  decon_mtrx = t(pbmc_st@assays$predictions@data)
  head(decon_mtrx)
  cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "max")]
  
  decon_df <- decon_mtrx %>%
    data.frame(check.names = F) %>%
    tibble::rownames_to_column("barcodes")
  pbmc_st@meta.data <- pbmc_st@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  return(pbmc_st)
}
  
  
########################################## 空间转录组相关函数
Load10X_Spatial_with_h5<-function(sample){
  obj_st = Load10X_Spatial(data.dir = sample,
                           filename = "filtered_feature_bc_matrix.h5",
                           assay = "Spatial", 
                           slice = "Slice",project=sample)
  head(obj_st@meta.data)
  obj_st@project.name <- sample
  #Idents(obj_st) <- sample
  obj_st$orig.ident <- sample
  return(obj_st)
}

Load10X_Spatial_without_h5<-function(sample){
  dir_mtx<-paste0(sample,'/filtered_feature_bc_matrix')
  dir_image<-paste0(sample,'/spatial')
  
  pbmc <- CreateSeuratObject(
    counts = Read10X( data.dir = dir_mtx ,  gene.column = 2),assay = 'Spatial',project=sample)
  my_image <- Read10X_Image( image.dir = dir_image,filter.matrix = F)
  image <- my_image[Cells(x = pbmc)]
  DefaultAssay(object = image) <- 'Spatial'
  pbmc[['Slice']] <- image
  pbmc@meta.data$orig.ident<-sample
  pbmc@meta.data$tissue<-pbmc@images$Slice@coordinates$tissue
  pbmc <- subset(pbmc,subset = tissue > 0)
  return(pbmc)
}

seurat_st_qc<-function(pbmc_st,plotid=0,species='human'){
  
  ##UMI统计画图
  plot1 <- VlnPlot(pbmc_st, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(pbmc_st, features = "nCount_Spatial") + theme(legend.position = "right")
  plota<-wrap_plots(plot1, plot2)
  ggsave(paste0("fig",plotid,"1.nCount_Spatial.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"1.nCount_Spatial.pdf"), plot = plota, width = 12, height = 6)
  
  plot1 <- VlnPlot(pbmc_st, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(pbmc_st, features = "nFeature_Spatial") + theme(legend.position = "right")
  plota<-wrap_plots(plot1, plot2)
  ggsave(paste0("fig",plotid,"2.nFeature_Spatial.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"2.nFeature_Spatial.pdf"), plot = plota, width = 12, height = 6)  
  
  mt.genes <- grep(pattern = "^MT-", x = rownames(pbmc_st), value = TRUE)
  string_pattern='^MT-'
  if(species=='mouse')
  {
    mt.genes <- grep(pattern = "^mt-", x = rownames(pbmc_st), value = TRUE)
    string_pattern='^mt-'
  }else if(species=='rat')
  {
    mt.genes <- grep(pattern = "^mt-", x = rownames(pbmc_st), value = TRUE)
    string_pattern='^mt-'
  }
  #pbmc[["percent.mt"]] <- (Matrix::colSums(pbmc[["Spatial"]]$counts[mt.genes, ])/Matrix::colSums(pbmc[["Spatial"]]$counts))*100
  pbmc_st[["percent.mt"]] <- PercentageFeatureSet(object = pbmc_st, pattern = string_pattern)
  
  plot1 <- VlnPlot(pbmc_st, features = "percent.mt", pt.size = 0.1) + NoLegend()
  dim(pbmc_st)
  plot2 <- SpatialFeaturePlot(pbmc_st, features = "percent.mt") + theme(legend.position = "right")
  plota<-wrap_plots(plot1, plot2)
  ggsave(paste0("fig",plotid,"3.percent.mt.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"3.percent.mt.pdf"), plot = plota, width = 12, height = 6) 
  
  genes_to_keep <- setdiff(names(which(Matrix::rowSums(pbmc_st[["Spatial"]]$counts)>5)),mt.genes)
  pbmc_st <- subset(pbmc_st,
                 features =genes_to_keep,  #针对基因
                 subset = nFeature_Spatial > 300 & percent.mt < 30 #针对spots
  )
  
  plot1 <- VlnPlot(pbmc_st, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(pbmc_st, features = "nCount_Spatial") + theme(legend.position = "right")
  plota<-wrap_plots(plot1, plot2)
  ggsave(paste0("fig",plotid,"4.nCount_Spatial.after_filter.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"4.nCount_Spatial.after_filter.pdf"), plot = plota, width = 12, height = 6)
  
  plot1 <- VlnPlot(pbmc_st, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(pbmc_st, features = "nFeature_Spatial") + theme(legend.position = "right")
  plota<-wrap_plots(plot1, plot2)
  ggsave(paste0("fig",plotid,"5.nFeature_Spatial.after_filter.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"5.nFeature_Spatial.after_filter.pdf"), plot = plota, width = 12, height = 6) 
  
  plot1 <- VlnPlot(pbmc_st, features = "percent.mt", pt.size = 0.1) + NoLegend()
  dim(pbmc_st)
  plot2 <- SpatialFeaturePlot(pbmc_st, features = "percent.mt") + theme(legend.position = "right")
  plota<-wrap_plots(plot1, plot2)
  ggsave(paste0("fig",plotid,"6.percent.mt.after_filter.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"6.percent.mt.after_filter.pdf"), plot = plota, width = 12, height = 6) 
  
  return(pbmc_st)
}

movefiles<-function(pattern='fig',res_dir)
{
  files<-list.files()
  for(fig in files)
  {
    if(grepl(pattern,fig))
    {
      file.rename(fig,paste0(res_dir,'/',fig))
    }
  }
}

seurat_st_normalization<-function(obj_st,plotid=2){
  obj_st <- SCTransform(obj_st, assay = "Spatial", verbose = FALSE)
  obj_st <- RunPCA(obj_st, assay = "SCT", verbose = FALSE)
  obj_st <- FindNeighbors(obj_st, reduction = "pca", dims = 1:30)
  obj_st <- FindClusters(obj_st, verbose = FALSE)
  obj_st <- RunUMAP(obj_st, reduction = "pca", dims = 1:30)
  
  p1 <- DimPlot(obj_st, reduction = "umap", label = TRUE,group.by='seurat_clusters')
  p2 <- SpatialDimPlot(obj_st, label = TRUE, label.size = 3,group.by='seurat_clusters')
  plota<-NULL
  plota<-wrap_plots(p1, p2)
  ggsave(paste0("fig",plotid,"a.umap_Spatial.clusters.tiff"), plot = plota, width = 12, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,"a.umap_Spatial.clusters.pdf"), plot = plota, width = 12, height = 6) 
  
  clusters<-unique(obj_st@meta.data$seurat_clusters)
  for(i in 1:length(clusters))
  {
    #i=1
    cluster<-clusters[i]
    plota<-NULL
    plota<-SpatialDimPlot(obj_st, cells.highlight = CellsByIdentities(object = obj_st, idents = cluster), 
                          facet.highlight = TRUE, pt.size.factor = 1.5 )
    ggsave(paste0("fig",plotid,"b.Spatial.cluster",cluster,".tiff"), plot = plota, width = 6, height = 6,compression='lzw')
    ggsave(paste0("fig",plotid,"b.Spatial.cluster",cluster,".pdf"), plot = plota, width = 6, height = 6) 
  }
  
  return(obj_st)
}

seurat_st_markers<-function(obj_st,markers,plotid=0){
  #markers<-m6A
  markers<-unique(markers)
  markers_used<-markers[markers %in% rownames(obj_st)]
  ngene<-length(markers_used)
  nsample<-length(pbmc_st@images)
  
  if(ngene==1)
  {
    plota<-NULL
    plota<-SpatialFeaturePlot(obj_st, features = markers,pt.size.factor = 1.5,ncol=nsample)
    ggsave(paste0("fig",plotid,".Spatial.markers.tiff"), plot = plota, width = 6*nsample, height = 6,compression='lzw',limitsize=F)
    ggsave(paste0("fig",plotid,".Spatial.markers.pdf"), plot = plota, width = 6*nsample, height = 6,limitsize=F)
    
  }else if(ngene<=6)
  {
    plota<-NULL
    plota<-SpatialFeaturePlot(obj_st, features = markers,pt.size.factor = 1.5,ncol=nsample)
    ggsave(paste0("fig",plotid,".Spatial.markers.tiff"), plot = plota, width = 6*nsample, height = ngene*6,compression='lzw',limitsize=F)
    ggsave(paste0("fig",plotid,".Spatial.markers.pdf"), plot = plota, width = 6*nsample, height = ngene*6,limitsize=F)
  }else
  {
    ncut<-ceiling(ngene/6)
    i<-0
    while(i<ncut)
    {
      markersi<-markers_used[(i*6+1):(i*6+6)]
      markersi<-na.omit(markersi)
      ngenei<-length(markersi)
      seurat_st_markers(obj_st,markersi,plotid=paste0(plotid,i))
      i=i+1
    }
  }
}

seurat_st_deg_cluster<-function(pbmc_st,plot.only=F,min.pct=0.25,logFCfilter=1,adjPvalFilter=0.05,plotid=0)
{
		library(ggplot2)
		library(presto)
		if(!plot.only)
		{
		DefaultAssay(pbmc_st) <- "SCT"
		Idents(pbmc_st) <- "seurat_clusters"
		scttemp <- SCTResults(pbmc_st, slot="umi.assay")
		for(i in 1:length(scttemp)) {
		  # i=1
		  slot(object = pbmc_st@assays$SCT@SCTModel.list[[i]], name="umi.assay")<-"Spatial"
		}
		pbmc_st<-PrepSCTFindMarkers(pbmc_st,assay = "SCT", verbose = TRUE)
		deg.markers <- FindAllMarkers(object = pbmc_st,only.pos = FALSE,min.pct = min.pct,logfc.threshold = 0,test.use = "wilcox",return.thresh=1)
		write.table(deg.markers,file=paste0("fig",plotid,".deg_cluster.all.tsv"),sep="\t",row.names=F,quote=F)
		sig.markers=deg.markers[(abs(as.numeric(as.vector(deg.markers$avg_log2FC)))>=logFCfilter & as.numeric(as.vector(deg.markers$p_val_adj))<=adjPvalFilter),]
		write.table(sig.markers,file=paste0("fig",plotid,".deg_cluster.sig.tsv"),sep="\t",row.names=F,quote=F)
		}else
		{
		sig.markers<-read.table(paste0("fig",plotid,".deg_cluster.sig.tsv"),header=T,sep='\t')
		}
		ncluster<-length(levels(pbmc_st))
		top5 <- sig.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
		top5.head<-top5$gene
		top5.head<-unique(top5.head)
		seurat_markers_plot(pbmc_st,top5.head,assay='SCT',group.by="seurat_clusters",
						  plotid=plotid,dotplot.only=F,plot.heatmap=T)
		seurat_st_markers(pbmc_st,top5.head,plotid=plotid)
}
    
seurat_st_deg_group<-function(obj_st,group.by='group',group_test='PE',group_control='control',plot.only=F,logFCfilter=1,adjPvalFilter=0.05,plotid=0)
{
  library(ggplot2)
  library(presto)
  if(!plot.only)
  {
    deg_group <- FindMarkers(obj_st,ident.1 = group_test, 
                             ident.2 = group_control,group_by=group.by, 
                             verbose = FALSE,only.pos = FALSE,min.pct = 0.1,logfc.threshold = 0)
    write.table(deg_group,file=paste0("fig",plotid,".deg_group.all.tsv"),sep="\t",row.names=F,quote=F)
    deg_group_sig=deg_group[(abs(as.numeric(as.vector(deg_group$avg_log2FC)))>logFCfilter & as.numeric(as.vector(deg_group$p_val_adj))<adjPvalFilter),]
    write.table(deg_group_sig,file=paste0("fig",plotid,".deg_group.sig.tsv"),sep="\t",row.names=F,quote=F)
  }else
  {
    deg_group_sig<-read.table(paste0("fig",plotid,".deg_group.sig.tsv"),header=T,sep='\t')
  }
  comparison=paste0(group_test,'_vs_',group_control)
  ggplot2_volcano(deg_group[,c('p_val_adj','avg_log2FC')],comparison=comparison,output=paste0('fig',plotid,'.volcano.',comparison),fc_threshold=2^logFCfilter,ylab='p_val_adj')
  
  ngene<-nrow(deg_group_sig)
  if(ngene>0)
  {
    top5_deg_group_sig <- deg_group_sig %>% top_n(n = 5, wt = avg_log2FC)
    top5<-top5_deg_group_sig$gene
    top5_deg_group_sig <- deg_group_sig %>% top_n(n = -5, wt = avg_log2FC)
    tail5<-top5_deg_group_sig$gene
    genes_used<-unique(append(top5,tail5))
    seurat_st_markers(obj_st,genes_used,plotid=plotid)
  }
}

seurat_st_cell_annotation_first<- function(pbmc_st,plotid=20,group.by = 'seurat_clusters',species='human',plot.full=F){
  ####  此函数中的marker的来自潘重2023年6月之前的单细胞数据分析的总结,后续再不停更新
  Idents(pbmc_st) <- group.by
  DefaultAssay(pbmc_st) <- "SCT"
  nc<-length(levels(pbmc_st))
  ph<-max(nc/2,6)
  th=theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
  markers_first = list(
    Epithelial=c('WFDC2','PAX8','EPCAM','GPC3','MMP7','KRT19','PROM1','ALDH1A1','CD24','KRT14','KRT5','JUP'),
    Fibroblasts=c('PDPN','FN1','VIM','DCN','GSN','COL1A1','COL1A2','COL4A1','COL4A2','PDGFRA'),
    SMC =  c('DES','ACTA2','TAGLN','MYH11','MYLK','ACTG2'),
    Pericyte=c('BCAM','PDGFRB','KCNJ8','RGS5','ABCC9'),
    Endothelial =  c('MME','PECAM1','THBD','CDH5','ENG','VWF','CCL21'),
    #mesothelial =  c('LY6G','MSIN'),
    #mesenchymal =  c('PDGFRA','TCF21','CDH11'),                   ###  32
    Neuronal=c('NRXN1'),
    Immune = c("PTPRC"),
    Lymphoid = c("CD69",'MKI67'),
    Tcells = c('CD3D','CD3E','CD3G','CD247','CD4','CD8A','CD8B'),
    NK=c('NCAM1','CD19','GNLY','KLRD1','NKG7','FCGR3A','FGFBP2','CX3CR1','KLRF1'),           ####20
    ILCs=c('AREG','TNFRSF18','TNFRSF25'),
    Bcells=c( 'CD79A','MS4A1'),
    Plasma=c( 'MZB1'),
    myeloid =  c("CD14",'LYZ','CD300E'),  
    Macrophages=c('CD68','CD86','CD163','FCER1G','MARCO'),
    Monocytes=c('FCN1','APOBEC3A','THBS1'),
    DC = c('HLA-DQB2','HLA-DPB1','BIRC3'),                                                 
    Mast =  c('MS4A2','TPSB2',"TPSAB1",'CST3','KIT'),                                       ##24
    #Neutrophils =  c("FCGR3B","CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2",'CSF3R'),           ###15
    #####     CD11b-->ITGAM   CD66b-->CEACAM8
    Neutrophils =  c("ITGAM","CEACAM8",'LY6G','Retnlg'),           ###15
    #https://mp.weixin.qq.com/s/lUZeQEPUSaoDumMnpspMQg
    #https://pubmed.ncbi.nlm.nih.gov/32719519/
    Basophils =  c('CD34','CD200R3','FCER1A','CD9',"CPA3",'ITGB2','IL3RA'),
    Eosinophils =  c('SIGLEC5', 'IL5RA', 'CCR3', 'EPX')
  )
  if(species =='mouse')
  {
    library(stringr)
    for(i in 1:length(markers_first))
    {
      #i=1
      markers_first[[i]]<-gene_human2mouse_vector(markers_first[[i]])
    }
  }
  p <- DotPlot(pbmc_st, features = markers_first,group.by = group.by )  +th
  ggsave(plot=p, filename=paste0("fig",plotid,".dotplot.markers.first.tiff"),compression='lzw',limitsize=F,width = 25, height =ph)
  if(plot.full)
  {
    plotid_new<-paste0(plotid,'a')
    markers<-as.character(unlist(markers_first))
    seurat_markers_plot(pbmc_st,markers,group.by=group.by,assay='SCT',plotid=plotid_new,dotplot.only=F,plot.heatmap=T)
    for(i in 1:length(markers_first))
    {
			markers<-markers_first[[i]]
			if(species =='mouse')
			  {
				library(stringr)
				  markers<-gene_human2mouse_vector(markers)
			  }  
			cellname<-names(markers_first)[i]
			  plotid_new<-paste0(plotid,'b.',i,'.',cellname)
			seurat_st_markers(pbmc_st,markers,plotid=plotid_new)
		}
  }
}


## object<-pbmc_st
#### top.features <-  rownames(pbmc_st[["predictions"]]@meta.features)[which(pbmc_st[["predictions"]]@meta.features$moransi.spatially.variable.rank < 6)]
#### top.clusters <- head(SpatiallyVariableFeatures_workaround(pbmc_st, assay="predictions",selection.method = "moransi"), 4)
SpatiallyVariableFeatures_workaround <- function(object, assay="SCT", selection.method = "moransi") {
  #' This is work around function to replace SeuratObject::SpatiallyVariableFeatures function.
  #' return ranked list of Spatially Variable Features
  
  # Check if object is a Seurat object
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  # Check if assay is a valid assay
  if (!assay %in% names(object@assays)) {
    stop("assay must be a valid assay")
  }
  
  # Extract meta.features from the specified object and assay
  data <- object@assays[[assay]]@meta.features
  
  # Select columns starting with the provided col_prefix
  moransi_cols <- grep(paste0("^", selection.method), colnames(data), value = TRUE)
  
  # Filter rows where "moransi.spatially.variable" is TRUE
  ## filtered_data <- data[data[[paste0(selection.method, ".spatially.variable")]], moransi_cols]
  filtered_data <- data[
    data[[paste0(selection.method, ".spatially.variable")]] & (!is.na(data[[paste0(selection.method, ".spatially.variable")]])), 
    moransi_cols
  ]
  filtered_data<-filtered_data[rownames(filtered_data)!='max',]
  
  # Sort filtered data by "moransi.spatially.variable.rank" column in ascending order
  sorted_data <- filtered_data[order(filtered_data[[paste0(selection.method, ".spatially.variable.rank")]]), ]
  
  # Return row names of the sorted data frame
  rownames(sorted_data)
}



## seurat_gsva_custom_violin(pbmc,genesets=gs,group.by='celltype_manual',projectid=projectid,plotid=plotid,species = "human",title="m6A score")
seurat_st_gsva_custom_violin <- function(pbmc_st,assay="SCT",genesets=gmt,group.by='celltype_manual',projectid='xxxxxx',plotid='14',species = "human",title="Ferroptosis score",plot.only=F){
  library(GSVA)
  library(msigdbr)
  #BiocManager::install('msigdbr')
  
  if(species=='human')
  {
    DBkeyset='org.Hs.eg.db'
    org='Homo sapiens'
    #library(org.Hs.eg.db)
    
  }else if(species=='mouse')
  {
    #library(org.Mm.eg.db)
    DBkeyset='org.Mm.eg.db'
    org='Mus musculus'
  }else if(species=='rat')
  {
    #library(org.Rn.eg.db)
    #BiocManager::install('org.Rn.eg.db')
    DBkeyset='org.Rn.eg.db'
    org='Mus musculus'
  }else{
    DBkeyset='org.Hs.eg.db'
    #library(org.Hs.eg.db)
    org='Homo sapiens'
  }
  if(plot.only)
  {
    df_plot<-read.delim(paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.tsv"),header = T)
  }else{
    DefaultAssay(pbmc) <- "RNA"
    Idents(pbmc_st) <- group.by
    levels(pbmc_st)
    expr <- GetAssayData(object = pbmc_st,assay="SCT",layer = "data")
    expr <- expr[rowSums(expr)>0,]  #选取非零基因
    expr <- as.matrix(expr)		
    gsvaPar <- ssgseaParam(expr, genesets,minSize=2)
    gsva.res <- gsva(gsvaPar, verbose=FALSE)
    #gsva.res <- gsva(expr, genesets, method="ssgsea",parallel.sz=10,min.sz=1) 
    #gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
    metadata<-pbmc_st@meta.data
    df_plot<-as.data.frame(t(gsva.res))
    dim(df_plot) 
    dim(metadata) 
    df_plot[[group.by]]<-metadata[[group.by]]
    head(df_plot)
    write.table(df_plot, paste0("fig",plotid,".",projectid,'.',group.by,".cell.gsva.tsv"),sep='\t',row.names = T,quote = FALSE)
  }
  colnames(df_plot)[1]<-'score'
  library(ggplot2)  
  ggplot2_boxplot_gsva(df_plot,groupby=group.by,score='score',y_title=title,pw=7,ph=4,output=paste0('fig',plotid,'.',projectid,".ssgsea")) 
  
  score_name=gsub(' ','_',title)
  pbmc_st@meta.data[[score_name]]<-df_plot$score
  
  plota<-FeaturePlot(pbmc_st, features =score_name,reduction = "umap", label = TRUE)
  ggsave(paste0("fig",plotid,".umap.",score_name,".tiff"), plot = plota, width = 8, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,".umap.",score_name,".pdf"), plot = plota, width = 8, height = 6)
  
  plota<-SpatialFeaturePlot(pbmc_st, features =score_name, pt.size.factor = 1.6, crop = TRUE)
  ggsave(paste0("fig",plotid,".SpatialFeaturePlot.",score_name,".tiff"), plot = plota, width = 6, height = 6,compression='lzw')
  ggsave(paste0("fig",plotid,".SpatialFeaturePlot.",score_name,".pdf"), plot = plota, width = 6, height = 6)
}



### pair=c("Epi-SRGN","Epi-LEPR")
### seurat_ST_plot_colocation(spatial_coord,pair,plotid='31',pw = 12, ph = 6) 
seurat_ST_plot_colocation<-function(spatial_coord,pair,plotid='31',pw = 8, ph = 6) 
{
  knn = 6
  pt.size=2
  alpha.min=0.1
  max.cut=0.95
  ####选择两种细胞类型
  LRpair = pair
  
  location = spatial_coord[,c('x','y')]
  topn=floor(0.2*dim(location)[1])
  expr = spatial_coord[,LRpair]
  ncell<-dim(expr)[1]
  nnmatrix<-RANN::nn2(location,k=knn)$nn.idx
  countsum<-Matrix::colSums(expr)
  ####normalize
  expr<-Matrix::t(log(Matrix::t(expr)/countsum*median(countsum)+1))
  
  ligand<-expr[,LRpair[1]]
  receptor<-expr[,LRpair[2]]
  
  LRexp<-rbind(ligand,receptor)
  neighexp<-apply(nnmatrix,1,function(x){apply(LRexp[,x[2:knn]],1,max)})
  
  LRadd<-pmax(LRexp[1,]*neighexp[2,],LRexp[2,]*neighexp[1,])
  LRadd_max<-quantile(LRadd,probs=max.cut)
  LRadd[LRadd>LRadd_max]<-LRadd_max
  if(sum(ligand>0)>topn){n1<-order(ligand,sample(ncell,ncell),decreasing=T)[1:topn]}else{n1<-which(ligand>0)}
  if(sum(receptor>0)>topn){n2<-order(receptor,sample(ncell,ncell),decreasing=T)[1:topn]}else{n2<-which(receptor>0)}
  expcol<-rep(0,ncell)
  expcol[n1]<-1
  expcol[n2]<-2
  expcol[intersect(n1,n2)]<-3
  tmp<-data.frame(x=location[,1],y=location[,2],Exp=as.factor(expcol))
  tmpLRadd<-data.frame(x=location[,1],y=location[,2],LR=LRadd)
  alpha=(LRadd-min(LRadd))/(max(LRadd)-min(LRadd))*(1-alpha.min)+alpha.min
  head(tmp)
  min(tmp$x)
  max(tmp$x)
  
  min(tmp$y)
  max(tmp$y)
  
  
  min(tmpLRadd$x)
  max(tmpLRadd$x)
  
  min(tmpLRadd$y)
  max(tmpLRadd$y)
  
  head(tmpLRadd)
  
  p1<-ggplot(tmp,aes(x=y,y=x,col=Exp))+geom_point(size=pt.size)+
    scale_color_manual(values=c("gray","blue","green","red"),labels=c("Bothlow",paste0(LRpair[1],"_high"),paste0(LRpair[2],"_High"),"BothHigh"))+
    ggtitle(paste0(LRpair,collapse="_n_"))+xlab("")+ylab("")+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),axis.ticks.y=element_blank()) + theme_minimal() + 
    theme(axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank())
  
  p2<-ggplot(tmpLRadd,aes(x=y,y=x,col=LR))+geom_point(size=pt.size,alpha=alpha)+
    scale_color_gradient2(midpoint=quantile(LRadd,probs=0.5),low="gray",high="red",mid="gray")+xlab("")+ylab("")+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),
          axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+labs(color="colocalization") + 
    theme_minimal() + theme(axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank())
  plota<-p1+p2&scale_y_reverse()
  ggsave(paste0("fig",plotid,".colocation.tiff"), plot = plota, width = pw, height = ph,compression='lzw')
  ggsave(paste0("fig",plotid,".colocation.pdf"), plot = plota, width = pw, height = ph)
  
}
#





###  pbmc@meta.data$celltype_manual是factor
seurat_ST_RCTD_calculation<-function(pbmc_st,pbmc,celltype='celltype_manual',mode='full',rds_output='output.RCTD.rds')
{
	if(1)  #####  RCTD准备reference
	{
		library(spacexr)	
		sc_counts <- GetAssayData(object = pbmc,assay="RNA",layer = "counts")	
		metadata <- pbmc@meta.data
		head(metadata)
		celltypes <- metadata[[celltype]]
		head(celltypes)
		names(celltypes) <- rownames(metadata)
		# celltypes <- factor(celltypes,levels=level) # convert to factor data type
		table(celltypes)	
		nUMI <- metadata$nCount_RNA
		names(nUMI) <- rownames(metadata)	
		reference <- Reference(sc_counts, celltypes, nUMI)
	}
	st_counts <- GetAssayData(object = pbmc_st,assay="SCT",layer = "counts")	
	coords<-seurat_ST_fetch_coords(pbmc_st)
	## rownames(coords) <- gsub("^.+\\.", "", rownames(coords))
	head(coords)
	dim(coords)
	query <- SpatialRNA(coords, st_counts, colSums(st_counts))
	RCTD <- create.RCTD(query, reference, max_cores = 32)
	RCTD_res <- run.RCTD(RCTD, doublet_mode = mode)   #####  full; doublet; multi	
	saveRDS(RCTD_res,rds_output)
	return(RCTD_res)
}

seurat_ST_fetch_coords<-function(pbmc_st)
{
	samples<-names(pbmc_st@images)
	#counts <- GetAssayData(object = pbmc_st,assay="Spatial",layer = "counts")        
	coordinates_list <- lapply(samples, function(image_name) {
	pos <- GetTissueCoordinates(pbmc_st, image = image_name)
	colnames(pos) <- c('x','y')
	return(pos)
	})
	coords <- do.call(rbind, coordinates_list)
	return(coords)
}

seurat_ST_add_celltypes_rctd<-function(pbmc_st,RCTD_res,celltype='celltype_manual',level=c('xxx','yy')) 
{
		df_norm_weights<-seurat_ST_fetch_rctd_celltypes(RCTD_res,celltype=celltype,level=level)
		celltypes <- unique(df_norm_weights[[celltype]])
		level_used<-level[level %in% celltypes]
		df_norm_weights[[celltype]]<-factor(df_norm_weights[[celltype]],levels=level_used)		
		pbmc_st <- AddMetaData(pbmc_st, metadata = df_norm_weights)
		return(pbmc_st)
}

seurat_ST_fetch_rctd_celltypes<-function(RCTD_res,celltype='celltype_manual',level=c('xxx','yy'))
{
	
	library(MatrixGenerics)
		results <- RCTD_res@results
		barcodes <- colnames(RCTD_res@spatialRNA@counts)
		weights <- RCTD_res@results$weights
		norm_weights <- normalize_weights(weights)
		
		head(norm_weights)
		class(norm_weights)
		df_norm_weights<-as.data.frame(norm_weights)
		df_norm_weights$max<-rowMaxs(as.matrix(df_norm_weights))
		head(df_norm_weights)
	
		df_norm_weights[[celltype]]<-'unknown'
		for(i in 1:nrow(df_norm_weights))
		{
		for(j in 1:(ncol(df_norm_weights)-2))
		{
			if(df_norm_weights[i,j]==df_norm_weights[i,'max'])
			{
			df_norm_weights[i,celltype]<-colnames(df_norm_weights)[j]
			}
		}
		}
		df_norm_weights$max<-NULL
		return(df_norm_weights)
}

seurat_ST_plot_celltypes<-function(pbmc_st,df_norm_weights,celltype='celltype_manual',plotid='04',pie_scale = 0.4) 
{
		samples<-names(pbmc_st@images)
		nsample<-length(samples)
		
		plota<-DimPlot(pbmc_st, reduction = "umap", label = TRUE,group.by = celltype)
		ggsave(paste0("fig",plotid,"c1.umap.",'celltype.nominal',".tiff"), plot = plota, width = 8, height = 6,compression='lzw',limitsize=F)
		ggsave(paste0("fig",plotid,"c1.umap.",'celltype.nominal',".pdf"), plot = plota, width = 8, height = 6,limitsize=F)
		
		plota<-NULL
		plota<-SpatialDimPlot(pbmc_st,group.by = celltype, label = TRUE, label.size = 3,ncol=nsample)
		ggsave(paste0("fig",plotid,"c2.SpatialPlot.celltype.nominal",".tiff"), plot = plota, width = 6*nsample, height = 4,compression='lzw',limitsize=F)
		ggsave(paste0("fig",plotid,"c2.SpatialPlot.celltype.nominal",".pdf"), plot = plota, width = 6*nsample, height = 4,limitsize=F)

		celltypes<-colnames(df_norm_weights)[1:(length(df_norm_weights)-1)]
		ncelltype<-length(celltypes)
		
		plota<-SpatialFeaturePlot(pbmc_st, features =celltypes, pt.size.factor = 1.6, ncol = nsample, crop = TRUE)
		ggsave(paste0("fig",plotid,"c3.SpatialFeaturePlot.",'celltypes',".tiff"), plot = plota, width = nsample*4, height = ncelltype*4,compression='lzw',limitsize=F)
		ggsave(paste0("fig",plotid,"c3.SpatialFeaturePlot.",'celltypes',".pdf"), plot = plota, width = nsample*4, height = ncelltype*4,limitsize=F)
		
		library(STdeconvolve)
		library(ggplot2)
		library(ggsci)
		packageVersion("STdeconvolve")
		for(i in 1:length(samples))
		{
		  # i =1
		  sample<-samples[i]
		pos <- GetTissueCoordinates(pbmc_st, image = sample)
		colnames(pos) <- c('x','y')
		head(pos)

		norm_weights<-df_norm_weights[rownames(pos),]
		norm_weights<-norm_weights[,-ncol(norm_weights)]
		norm_weights <- as.matrix(norm_weights)

		plt <- vizAllTopics(theta = norm_weights,
							pos = pos,
							topicOrder=seq(ncol(norm_weights)),
							topicCols=rainbow(ncol(norm_weights)),
							groups = NA,
							group_cols = NA,
							r = 3, # size of scatterpies; adjust depending on the coordinates of the pixels
							lwd = 0.3,
							showLegend = TRUE,
							plotTitle = "scatterpies")

		## function returns a `ggplot2` object, so other aesthetics can be added on:
		plt <- plt + ggplot2::guides(fill=ggplot2::guide_legend(ncol=2));plt
		ggsave(paste0("fig",plotid,"c4.Spaital_scatterpies.",sample,".pdf"), width=9, height=6, plot=plt, bg="white",limitsize=F)
		ggsave(paste0("fig",plotid,"c4.Spaital_scatterpies.",sample,".tiff"), width=9, height=6, plot=plt, bg="white",compression='lzw',limitsize=F)
		}
		
		for(i in 1:length(samples))
		{
		  # i =1
		  sample<-samples[i]
		  pos <- GetTissueCoordinates(pbmc_st, image = sample)
		  colnames(pos) <- c('x','y')
		  head(pos)
		  
		  norm_weights<-df_norm_weights[rownames(pos),]
		  norm_weights<-norm_weights[,-ncol(norm_weights)]
		  norm_weights <- as.matrix(norm_weights)
		  
		  library(SPOTlight)
		  class(norm_weights)
		  cell_types_all <- colnames(norm_weights)
		  ncelltype<-dim(norm_weights)[2]
		  
		  paletteMartin <- c(
			"#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
			"#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
			"#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
		  
		  pal <- colorRampPalette(paletteMartin)(ncelltype)
		  names(pal) <- cell_types_all
		  
		  mycolors<-get_colors(ncelltype)
		  names(mycolors) <- cell_types_all
		  

		  pos_new<-pos
		  pos_new$x<-max(pos$x)-pos$x
		  head(pos_new)
		  
		  p3 <- SPOTlight::plotSpatialScatterpie(
			x = pos_new,
			y = norm_weights,
			cell_types = colnames(norm_weights),
			img = F, #以tiff为背景
			scatterpie_alpha = 1,
			pie_scale = pie_scale)+scale_fill_manual(values = mycolors,breaks = names(mycolors))
		  ggsave(paste0("fig",plotid,"c5.Spaital_scatterpies.SPOTlight.",sample,".pdf"), width=9, height=6, plot=p3, bg="white",limitsize=F)
		  ggsave(paste0("fig",plotid,"c5.Spaital_scatterpies.SPOTlight.",sample,".tiff"), width=9, height=6, plot=p3, bg="white",compression='lzw',limitsize=F)
		  
			graph_heatmap <- SPOTlight::plotInteractions(norm_weights,"heatmap")
			ggsave(paste0("fig",plotid,"c6.plotInteractions.SPOTlight.",sample,".heatmap.tiff"), plot = graph_heatmap, width = 6, height = 6,compression='lzw',limitsize=F)
			ggsave(paste0("fig",plotid,"c6.plotInteractions.SPOTlight.",sample,".heatmap.pdf"), plot = graph_heatmap, width = 6, height = 6,limitsize=F) 
			
			pdf(file=paste0("fig",plotid,"c7.plotInteractions.SPOTlight.",sample,".network.pdf"),width=8,height=8)
			print(SPOTlight::plotInteractions(norm_weights,"network"))
			dev.off()
			tiff(file=paste0("fig",plotid,"c7.plotInteractions.SPOTlight.",sample,".network.tiff"),width=8*300,height=8*300,compression='lzw',res=300)
			print(SPOTlight::plotInteractions(norm_weights,"network"))
			dev.off()
			
			plota<-plotCorrelationMatrix(norm_weights)
			ggsave(paste0("fig",plotid,"c8.plotCorrelationMatrix.SPOTlight.",sample,".tiff"), plot = plota, width = 6, height = 6,compression='lzw',limitsize=F)
			ggsave(paste0("fig",plotid,"c8.plotCorrelationMatrix.SPOTlight.",sample,".pdf"), plot = plota, width = 6, height = 6,limitsize=F) 
		}
}


seurat_ST_plot_celltype_count<-function(pbmc_st,df_norm_weights,celltype='celltype_manual',plotid='04',level=c('xxx','yy')) 
{
		samples<-names(pbmc_st@images)
		nsample<-length(samples)
		
		celltypes<-colnames(df_norm_weights)[1:(length(df_norm_weights)-1)]
		ncelltype<-length(celltypes)
		
		for(i in 1:length(samples))
		{
			# i =1
			sample<-samples[i]
			pos <- GetTissueCoordinates(pbmc_st, image = sample)
			colnames(pos) <- c('x','y')
			head(pos)


			norm_weights<-df_norm_weights[rownames(pos),]
			norm_weights<-norm_weights[,-ncol(norm_weights)]
			norm_weights <- as.matrix(norm_weights)

			head(norm_weights)
			rowSums(norm_weights)

			df_ncell<-colSums(norm_weights)
			df_ncell<-as.data.frame(df_ncell)
			colnames(df_ncell)<-sample
			df_ncell_percentage<-df_ncell/sum(df_ncell[,1])
			ggplot2_barplot(df_ncell_percentage,group=NULL,value=sample,title='celltype',y_title='Percentage',output=paste0("fig",plotid,".",sample,".b1.barplot.celltype.count"))
			df_ncell_percentage_used<-df_ncell_percentage[df_ncell_percentage[[sample]]>0,,drop=F]
			ggplot2_pie_with_legend(df_ncell_percentage_used,group=NULL,count=sample,title='celltype',output=paste0("fig",plotid,".",sample,".b2.pie.celltype.percentage"))

			df_count<-table(df_norm_weights[[celltype]])
			df_count<-as.data.frame(df_count)
			colnames(df_count)<-c(celltype,'count')
			rownames(df_count)<-df_count[,1]
			df_count$percentage<-df_count[,2]/sum(df_count[,2])
			df_percentage_nominal<-df_count[,'percentage',drop=F]
			ggplot2_pie_with_legend(df_percentage_nominal,group=NULL,count='percentage',title='celltype',output=paste0("fig",plotid,".",sample,".b3.pie.celltype.percentage.nominal"))

		}

		if(1)  #####  总体细胞比例作图
		{
			plotid_new<-paste0(plotid,".total")
			df_nspots<-table(pbmc_st@meta.data$orig.ident)
			df_nspots<-as.data.frame(df_nspots)
			colnames(df_nspots)<-c('sample','count')
			rownames(df_nspots)<-df_nspots[,1]
			ggplot2_barplot(df_nspots,group='sample',value='count',title='Numer of Spots',y_title='Number',output=paste0("fig",plotid_new,".a.barplot.Nspots"))

			norm_weights<-df_norm_weights[,-ncol(df_norm_weights)]
			norm_weights <- as.matrix(norm_weights)
			head(norm_weights)

			df_ncell<-colSums(norm_weights)
			df_ncell<-as.data.frame(df_ncell)
			colnames(df_ncell)<-'count'
			df_ncell_percentage<-df_ncell/sum(df_ncell[,1])
			ggplot2_barplot(df_ncell_percentage,group=NULL,value='count',title='celltype',y_title='Percentage',output=paste0("fig",plotid_new,".b1.barplot.celltype.count"))
			df_ncell_percentage_used<-df_ncell_percentage[df_ncell_percentage[['count']]>0,,drop=F]
			ggplot2_pie_with_legend(df_ncell_percentage_used,group=NULL,count='count',title='celltype',output=paste0("fig",plotid_new,".b2.pie.celltype.percentage"))
			
			
			library(reshape2)
			df_ncells<-melt(as.matrix(norm_weights))
			head(df_ncells)
			colnames(df_ncells)<-c('cellid','celltype','count')
			
			df_plot<-merge(df_ncells,pbmc_st@meta.data[,c('orig.ident','group')],by.x='cellid',by.y=0)
			head(df_plot)
			dim(df_ncells)
			dim(df_plot)
			
			df_plot2<-aggregate(df_plot$count, by = list(df_plot$orig.ident,df_plot$celltype), FUN = sum)
			colnames(df_plot2)<-c('sample','celltype','count')
			head(df_plot2)
			
			df_percentage <- df_plot2 %>% group_by(sample) %>% dplyr::mutate(pct=prop.table(count))
			head(df_percentage)
			df_percentage<-as.data.frame(df_percentage)
			
			if(is.null(level))
			{
			  print(levels(df_percentage$celltype))
			}else
			{
			  df_percentage$celltype<-factor(df_percentage$celltype,levels =rev(level) )
			}
			
			df_plot3<-df_percentage[df_percentage$count>0,]
			ggplot2_barplot_strack(df_plot3,x='sample',y='pct',type='celltype',pw=8,ph=6,
								   output=paste0("fig",plotid_new,".c1.",'celltype_percentage.by_sample'))
			ggplot2_barplot_strack_with_label(df_plot3,x='sample',y='pct',type='celltype',pw=8,ph=6,
											  output=paste0("fig",plotid_new,".c1.",'celltype_percentage.by_sample.labled'))
											  
											  
			df_plot2<-aggregate(df_plot$count, by = list(df_plot$group,df_plot$celltype), FUN = sum)
			colnames(df_plot2)<-c('group','celltype','count')
			head(df_plot2)
			
			df_percentage <- df_plot2 %>% group_by(group) %>% dplyr::mutate(pct=prop.table(count))
			head(df_percentage)
			df_percentage<-as.data.frame(df_percentage)
			
			if(is.null(level))
			{
			  print(levels(df_percentage$celltype))
			}else
			{
			  df_percentage$celltype<-factor(df_percentage$celltype,levels =rev(level) )
			}
			
			df_plot3<-df_percentage[df_percentage$count>0,]
			ggplot2_barplot_strack(df_plot3,x='group',y='pct',type='celltype',pw=8,ph=6,
								   output=paste0("fig",plotid_new,".c2.",'celltype_percentage.by_group'))
			ggplot2_barplot_strack_with_label(df_plot3,x='group',y='pct',type='celltype',pw=8,ph=6,
											  output=paste0("fig",plotid_new,".c2.",'celltype_percentage.by_group.labled'))
		}

}

seurat_ST_colocation_pcc<-function(pbmc_st,df_norm_weights,plotid='04',pcc_cutoff=0.6,pvalue_cutoff=0.05) 
{
		samples<-names(pbmc_st@images)
		nsample<-length(samples)
		
		celltypes<-colnames(df_norm_weights)[1:(length(df_norm_weights)-1)]
		ncelltype<-length(celltypes)
		
		for(i in 1:length(samples))
		{
		  # i=1
		  
		  sample<-samples[i]

		  plotid_new1<-paste0(plotid,".",sample)
		  plotid_new2<-paste0(plotid+i,".",sample)
		  
		  metadata<- data.frame(pbmc_st@meta.data)
		  colnames(metadata) <- colnames(pbmc_st@meta.data)
		  head(metadata)
		  metadata_ds<-metadata[metadata$orig.ident==sample,]
		  
		  pos <- GetTissueCoordinates(pbmc_st, image = sample)
		  colnames(pos) <- c('x','y')
		  head(pos)
		  
		  spatial_coord<-merge(metadata_ds[,celltypes],pos,by=0)
		  head(spatial_coord)
		  rownames(spatial_coord)<-spatial_coord[,1]
		  spatial_coord<-spatial_coord[,-1]
		  

		  df_expr = spatial_coord[,celltypes,drop=F]
		  head(df_expr)
		  sds<-apply(df_expr,2,sd)
		  df_expr_cor<-df_expr[,sds>0]
		  ggcorrplot_1matrix(df_expr_cor,cor_method='pearson',output=paste0("fig",plotid_new1,".a.ggcorrplot"),hc.order = T,pw=6,ph=6)
		  
		  matrix_cor<-read.delim(paste0("fig",plotid_new1,".a.ggcorrplot.correlation.r.tsv"),row.names = 1,check.names=F)
		  head(matrix_cor)
		  dim(matrix_cor)
		  values<-as.numeric(unlist(matrix_cor))
		  values<-values[values!=1]
		  max(values)
		  min(values)
		  
		  matrix_cor_p<-read.delim(paste0("fig",plotid_new1,".a.ggcorrplot.correlation.pvalue.tsv"),row.names = 1,check.names=F)
		  head(matrix_cor_p)
		  dim(matrix_cor_p)
		  
		  cell_types_used<-colnames(df_expr_cor)
		  for(m in 1:(length(cell_types_used)-1))
		  {
			for(n in (m+1):length(cell_types_used))
			{
			  ### m=1
			  ### n=2
			  c1<-cell_types_used[m]
			  c2<-cell_types_used[n]
			  pair=c(c1,c2)
			  plotid_used=paste0(plotid_new2,'.b.',m,'_',n)
			  pcc<-matrix_cor[c1,c2]
			  pvalue<-matrix_cor_p[c1,c2]
			  if(abs(pcc)>0.6 & pvalue<0.05)
			  {
				if(1)
				{
				  seurat_ST_plot_colocation(spatial_coord,pair,plotid=plotid_used,pw = 12, ph = 6) 
				}
				if(1)
				{
				  plota<-SpatialPlot(object = pbmc_st,images=sample, features = pair, ncol = 2)
				  ggsave(paste0("fig",plotid_used,".SpatialPlot.",".tiff"), plot = plota, width = 8, height = 4,compression='lzw')
				  ggsave(paste0("fig",plotid_used,".SpatialPlot.",".pdf"), plot = plota, width = 8, height = 4)
				}
				if(1)
				{
				  expr = spatial_coord[,pair]
				  head(expr)
				  ggplot2_scatter_with_density(expr,output=paste0("fig",plotid_used,".scatterplot"),method="pearson")
				}
			  }
			}
		  }
		  
		}
		
		
	
}
