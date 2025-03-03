#####  行是样本，列是细胞类型
cibersort_plot<-function(df_infiltration,phenotype,groupby='group',plotid=10)
{
	library(ggplot2)
library(pheatmap)
    head(df_infiltration)
    if(1)#堆积柱状图
    {
      data=t(df_infiltration)
      col=rainbow(nrow(data),s=1,v=1)
      #绘制柱状图
      pdf(paste0('fig',plotid,"a.cibersort",".barplot2.pdf"),height=10,width=22)
      par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
      a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8,axisnames=F)
      a2=axis(2,tick=F,labels=F)
      axis(2,a2,paste0(a2*100,"%"))
      #axis(1,a1,labels=F)
      par(srt=60,xpd=T);
      #text(a1,-0.02,colnames(data),adj=1,cex=0.6);
      par(srt=0)
      ytick2 = cumsum(data[,ncol(data)])
      ytick1 = c(0,ytick2[-length(ytick2)])
      legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
      dev.off()
      
      tiff(paste0('fig',plotid,"a.cibersort",".barplot2.tiff"),height=10*300,width=22*300,units="px",res=300,compression='lzw')
      par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
      a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8,axisnames=F)
      a2=axis(2,tick=F,labels=F)
      axis(2,a2,paste0(a2*100,"%"))
      #axis(1,a1,labels=F)
      par(srt=60,xpd=T);
      #text(a1,-0.02,colnames(data),adj=1,cex=0.6);
      par(srt=0)
      ytick2 = cumsum(data[,ncol(data)])
      ytick1 = c(0,ytick2[-length(ytick2)])
      legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
      dev.off()
    }
    if(1)   ##### 使用生信自学网的脚本绘制相关性图形
    {
      library(corrplot)  
      #绘制相关性图形
      pdf(paste0('fig',plotid,"a.cibersort",".corrplot.pdf"),height=13,width=13)              #保存图片的文件名称
      par(oma=c(0.5,1,1,1.2))
      data=df_infiltration
      head(data)
      dim(data)
      data=data[,colMeans(data)>0]
      M=cor(data)
      corrplot(M,
               order="hclust",
               method = "color",
               addCoef.col = "black",
               diag = TRUE,
               tl.col="black",
               #col=colorRampPalette(c("blue", "white", "red"))(50)
      )
      dev.off()
      
      tiff(paste0('fig',plotid,"a.cibersort",".corrplot.tiff"),height=13*300,width=13*300,units="px",res=300,compression='lzw')              #保存图片的文件名称
      par(oma=c(0.5,1,1,1.2))
      data=df_infiltration
      head(data)
      dim(data)
      data=data[,colMeans(data)>0]
      M=cor(data)
      corrplot(M,
               order="hclust",
               method = "color",
               addCoef.col = "black",
               diag = TRUE,
               tl.col="black",
               #col=colorRampPalette(c("blue", "white", "red"))(50)
      )
      dev.off()
    }
    if(1)  ##### 使用生信自学网的脚本绘制热图
    {
      library(stringr)
      data=t(df_infiltration)
      phenotype_used<-phenotype[,groupby,drop=F]
      groups<-levels(phenotype_used[,groupby])
      
      annotation_col=phenotype_used
      head(annotation_col)
      
      ann_colors = list(groupby= c('a' ="#ce6700",'b'="#1b8c19") )
      class(ann_colors$groupby)
      names(ann_colors$groupby)<-groups
      
      pdf(paste0('fig',plotid,"c.cibersort",".heatmap.pdf"),height=6,width=12)
      pheatmap::pheatmap(data, 
                         annotation=annotation_col,
                         annotation_colors = ann_colors,
                         color = colorRampPalette(c('black',"red", "yellow"))(50),
                         show_colnames=F,
                         cluster_cols =F,
                         fontsize = 8,
                         fontsize_row=8,
                         fontsize_col=5)
      dev.off()
      tiff(paste0('fig',plotid,"c.cibersort",".heatmap.tiff"),height=6*300,width=12*300,units="px",res=300,compression='lzw')
      pheatmap::pheatmap(data, 
                         annotation=annotation_col,
                         annotation_colors = ann_colors,
                         color = colorRampPalette(c('black',"red", "yellow"))(50),
                         show_colnames=F,
                         cluster_cols =F,
                         fontsize = 8,
                         fontsize_row=8,
                         fontsize_col=5)
      dev.off()
    }
    if(1)   #### boxplot 像tcga一样带tumor和normal两种配对的样本。
    {
      data<-t(df_infiltration)
      phenotype_used<-phenotype[,groupby,drop=F]
      #####   data的行为celltype
      ggplot2_boxplot_paired_cibersort(data, phenotype_used,groupby=groupby,plotid=paste0('fig',plotid,'e'),ymax=0.9,label_y=0.75)
    }
}


#### 每种细胞类型分开作boxpot，并且比较两组间的差异。   
cibersort_celltype_boxplot<-function(df_infiltration,phenotype,groupby='group',plotid=10)
{
	library(ggplot2)
	library(pheatmap)
 
  data=df_infiltration
  head(data)
  dim(data)
  class(data)
  head(phenotype)
  dim(phenotype)
  for(i in 1:ncol(data))
  {
     #i<-6
     mydata<-cbind(data[,i,drop=F],phenotype[,groupby,drop=F])
     head(mydata)
     class(mydata)
     dim(mydata)
     
     groups<-levels(phenotype[,groupby])
     
     mydata_low<-mydata[mydata[,groupby]==groups[1],]
     dim(mydata_low)
     mydata_high<-mydata[mydata[,groupby]==groups[1],]
     meana<-mean(mydata_low[,1])
     meanb<-mean(mydata_high[,1])
     celltype<-colnames(data)[i]
     #celltype<-gsub('[.]',' ',celltype)
     #celltype<-gsub('CD8  ','CD8+ ',celltype)
     #celltype<-gsub('CD4  ','CD4+ ',celltype)
     ymax<-max(mydata[,1])
     ggboxplot_2groups_comparison(mydata,group=groupby,value=colnames(data)[i],title=celltype,output=paste0('fig',plotid,'.boxplot.',i),label_y=ymax+0.02,pw=4)
  }

}

