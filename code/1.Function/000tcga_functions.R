
project_all<-c('ACC','BLCA','BRCA','CESC','CHOL',
            'COAD','COADREAD','DLBC','ESCA','FPPP',
            'GBM','GBMLGG','HNSC','KICH','KIRC',
            'KIRP','LAML','LGG','LIHC','LUAD','LUNG',
            'LUSC','MESO','OV','PAAD','PCPG',
            'PRAD','READ','SARC','SKCM','STAD','TGCT',
            'THCA','THYM','UCEC','UCS','UVM',
            'PANCAN')
tcga_rna_seq_dir<-'/data/panzhong/genome_sequence/TCGA_expression_matrix'
tcga_rna_seq_workdir<-'/data/panzhong/000tcga'


## mycolors<-get_colors(2)
get_colors<-function(n)
{
  if(n<=9)
  {
    library(ggsci)
    library("scales")
    library(RColorBrewer)
    
    #mycolors= pal_lancet('lanonc')(9)
    mycolors<-brewer.pal(9,"Set1")
    mycolors<-mycolors[1:n]
    #show_col(mycolors)
  }else if(n<=20)
  {
    library(ggsci)
    library("scales")
    mycolors= pal_d3('category20')(20)
    # show_col(mycolors)
    mycolors<-mycolors[1:n]
  }else if(n<=51)
  {
    library(ggsci)
    library("scales")
    mycolors= pal_igv()(51)
    #show_col(mycolors)
    mycolors<-mycolors[1:n]
  }else
  {
    library(viridis)
    mycolors<-viridis(n, option = "D")
    #show_col(mycolors)
  }
  return(mycolors)
}


#tcga_rna_seq_dir<-'G:/rawdata/TCGA_expression_matrix'
#' get_tcga_projects_data(projects=c('COAD'),datatype='tpm',tcga_rna_seq_dir=tcga_rna_seq_dir)
get_tcga_projects_data<-function(projects=c('COAD'),datatype='tpm',tcga_rna_seq_dir='G:/rawdata/TCGA_expression_matrix')
{
  nproject<-length(projects)
  df_count<-data.frame()
  datadir<-paste0(tcga_rna_seq_dir,'/count')
  filetype='.unstranded.csv'
  if(datatype=='count')
  {
    datadir<-paste0(tcga_rna_seq_dir,'/count')
    filetype='.unstranded.csv'
  }else if(datatype=='tpm')
  {
    datadir<-paste0(tcga_rna_seq_dir,'/tpm')
    filetype='.tpm_unstranded.csv'
  }

  if(nproject>1)
  {
    project<-projects[1]
    df_count_file=paste0(datadir,'/TCGA-',project,filetype)
    df_count <- read.csv(file = df_count_file,header=T,check.names=F)

    for(i in 2:nproject)
    {
      project<-projects[i]
      df_count_file=paste0(datadir,'/TCGA-',project,filetype)
      df_count_tmp <- read.csv(file = df_count_file,header=T,check.names=F)
      df_count<-merge(df_count,df_count_tmp,by='gene_id')
      rm(df_count_tmp)
      gc()
    }
  }else{
    project<-projects[1]
    df_count_file=paste0(datadir,'/TCGA-',project,filetype)
    df_count <- read.csv(file = df_count_file,header=T,check.names=F)
  }
  return(df_count)
}

#' get_tcga_projects_data_symbol(projects=c('COAD'),datatype='tpm')
get_tcga_projects_data_symbol<-function(projects=c('COAD'),datatype='tpm',tcga_rna_seq_dir='G:/rawdata/TCGA_expression_matrix')
{
  nproject<-length(projects)
  df_count<-data.frame()
  datadir<-paste0(tcga_rna_seq_dir,'/count')
  filetype='.unstranded.csv'
  if(datatype=='count')
  {
    datadir<-paste0(tcga_rna_seq_dir,'/count')
    filetype='.unstranded.csv'
  }else if(datatype=='tpm')
  {
    datadir<-paste0(tcga_rna_seq_dir,'/tpm')
    filetype='.tpm_unstranded.csv'
  }

  if(nproject>1)
  {
    project<-projects[1]
    df_count_file=paste0(datadir,'/TCGA-',project,filetype)
    df_count <- read.csv(file = df_count_file,header=T,check.names=F)

    for(i in 2:nproject)
    {
      project<-projects[i]
      df_count_file=paste0(datadir,'/TCGA-',project,filetype)
      df_count_tmp <- read.csv(file = df_count_file,header=T,check.names=F)
      df_count<-merge(df_count,df_count_tmp,by='gene_id')
      rm(df_count_tmp)
      gc()
    }
  }else{
    project<-projects[1]
    df_count_file=paste0(datadir,'/TCGA-',project,filetype)
    df_count <- read.csv(file = df_count_file,header=T,check.names=F)
  }
  id_map_gene_file=paste0(tcga_rna_seq_dir,"/tcga.rna_seq.id.map.csv")
  id_map_gene=read.csv(id_map_gene_file,header=T)
  #head(counts[,1:3])
  if(1)     ####  将count矩阵的ID由EnsemblID转换成Symbol
  {
    df_count<-subset(df_count,grepl('^ENS',df_count$gene_id))
    #删除Y染色体上的冗余基因
    df_count<-subset(df_count,!grepl('_PAR_',df_count$gene_id))
    #添加gene_name信息，并筛选只保留编码基因
    #expr<-merge(counts,subset(id_map_gene,id_map_gene$gene_type=='protein_coding'),by.x = 'gene_id', by.y = 'gene_id')
    df_count<-merge(df_count,id_map_gene,by.x = 'gene_id', by.y = 'gene_id')
    #expr[1:3,]
    #expr[1:3,1:3]
    #dim(expr)
    #dim(expr)
    #删除指定列
    df_count$gene_id<-NULL
    df_count$gene_type<-NULL
    #df_data1[1:3,1:3]
    #dim(df_data1)
    #调整列顺序
    #cols <- colnames(expr)
    #new_cols <- c(cols[length(cols)], cols[1:(length(cols)-1)])
    #expr <- expr[, new_cols]
    #dim(df_data1)
    #删除基因名重复的行
    df_count<-df_count[!duplicated(df_count$gene_name),]
    #df_count[1:3,1:3]
    rownames(df_count)<-df_count$gene_name
    df_count$gene_name<-NULL
    #df_counts=df_counts[complete.cases(df_counts), ]
    dim(df_count)
  }
  return(df_count)
}


#### mydata的数据格式是列明：OS_STATUS,OS_MONTHS,hubgenes
univar_cox<-function(mydata,projectid='TCGA',plotid='08',sur='OS_STATUS',sur_time='OS_MONTHS')
{
  library(survival)
  library(survminer)
  #mydata<-train_data
  ngene<-ncol(mydata)
  outTab=data.frame()
  for(j in 3:ngene){
    #j=3
    gene_name=colnames(mydata)[j]
    mydata_used=mydata[,c(gene_name,sur,sur_time)]
    mydata_used[[sur]]=as.numeric(unlist(mydata_used[[sur]]))
    class(mydata_used[[sur]])
    
    cox=coxph(Surv(mydata_used[[sur_time]],mydata_used[[sur]]) ~ mydata_used[,1], data = mydata_used)
    coxSummary = summary(cox)
    outTab=rbind(outTab,
                 cbind(cancer=projectid,
                       gene_name=colnames(mydata_used)[1],
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                       sample_num=dim(mydata_used)[1]
                 ) )
  }
  outfile=paste0('fig',plotid,'.univariate_cox.result.tsv')
  write.table(outTab,file=outfile,sep="\t",row.names=F,quote=F)
  outTab$pvalue<-as.numeric(outTab$pvalue)
  cox_genes<-outTab[outTab$pvalue<0.05,'gene_name']
  cox_genes<-as.vector(na.omit(cox_genes))
  return(cox_genes)
}

lasso_cox_filter<-function(mydata,projectid='TCGA-LIHC',plotid='21',sur='OS_STATUS',sur_time='OS_MONTHS',pw=6)
{
  #mydata<-df_clinical_used
  mydata1<-mydata[,3:(ncol(mydata))]
  mydata2<-mydata[,c(sur_time,sur)]
  colnames(mydata2)<-c('time','status')
  head(mydata2)
  y <- data.matrix(Surv(mydata2$time,mydata2$status))
  head(y)
  mydata1<-as.matrix(mydata1)
  if(1)
  {
    fit <- glmnet(mydata1, y, family = "cox", alpha=1)
    pdf(paste0("fig",plotid,".lasso.lambda.pdf"), height=pw, width=pw)
    plot(fit, xvar = "lambda", label = TRUE)
    dev.off()
    tiff(paste0("fig",plotid,".lasso.lambda.tiff"), height=pw*300, width=pw*300,units="px",res=300,compression='lzw')
    plot(fit, xvar = "lambda", label = TRUE)
    dev.off()

    cvfit <- cv.glmnet(mydata1, y, family="cox", alpha=1,nfolds=10)
    pdf(paste0("fig",plotid,".lasso.cvfit.pdf"), height=pw, width=pw)
    plot(cvfit)
    abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
    dev.off()

    tiff(paste0("fig",plotid,".lasso.cvfit.tiff"), height=pw*300, width=pw*300,units="px",res=300,compression='lzw')
    plot(cvfit)
    abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
    dev.off()
	lambda_best<-cvfit$lambda.min
	print(lambda_best)
    coef <- coef(cvfit, s = cvfit$lambda.min)
    index <- which(coef != 0)
    actCoef <- coef[index]
    lassoGene=row.names(coef)[index]
    df_coef<-cbind(lassoGene,actCoef)
    print(lassoGene)
    write.table(df_coef,file=paste0("fig",plotid,"a.lasso.coef.txt"),sep="\t",row.names=F,quote=F)
    lassoGene<-lassoGene[lassoGene!='(Intercept)']
    lassoSigExp=mydata[,lassoGene]
    lassoSigExp=cbind(id=row.names(lassoSigExp),time=mydata2$time,status=mydata2$status,lassoSigExp)
    write.table(lassoSigExp,file=paste0("fig",plotid,"a.lasso.SigExp.txt"),sep="\t",row.names=F,quote=F)
    
    mydata3<-mydata1[,lassoGene]
    fit_final <- cv.glmnet(mydata3, y, family="cox", alpha=1)
    print(fit_final$lambda.min)    
    coef <- coef(fit_final, s = fit_final$lambda.min)
    print(coef)
    coef<-as.data.frame(as.matrix( coef))
    colnames(coef)<-'coefficients'
    write.table(coef,file=paste0("fig",plotid,"a.lasso.coef.final.txt"),sep="\t",row.names=T,quote=F)
    
    saveRDS(fit_final,paste0("fig",plotid,"a.lasso_model_final.rds"))
  }
  return(lassoGene)
}

#calculate the risk score of each sample.
riskscore <- function(survival_cancer_df, candidate_genes_for_cox, cox_report) {
  #survival_cancer_df<-lassoSigExp
  #candidate_genes_for_cox<-lassoGene_multi
  #cox_report<-multi_variate_cox_2
  library('dplyr')
  head(survival_cancer_df)
  risk_score_table <- survival_cancer_df[,candidate_genes_for_cox]
  head(risk_score_table)
  head(survival_cancer_df)
  coxSummary<-summary(cox_report)$coefficients
  riskscore_value<-crossprod(t(risk_score_table),coxSummary$actCoef)
  head(risk_score_table)
  risk_score_table <- cbind(risk_score_table, 'risk_score'=riskscore_value) %>%
    cbind(survival_cancer_df[,c('time','status')])
  risk_score_table <- risk_score_table[,c('time','status', candidate_genes_for_cox, 'risk_score')]
  return(risk_score_table)
}


riskscore_lasso <- function(survival_cancer_df, candidate_genes_for_cox, lasso_fit) {
  #survival_cancer_df<-lassoSigExp
  #candidate_genes_for_cox<-lassoGene_multi
  #cox_report<-multi_variate_cox_2
  library('dplyr')
  head(survival_cancer_df)
  risk_score_table <- survival_cancer_df[,candidate_genes_for_cox]
  head(risk_score_table)
  head(survival_cancer_df)
  coef <- coef(lasso_fit, s = lasso_fit$lambda.min)
  
  riskscore_value<-crossprod(t(risk_score_table),coef[,1])
  head(risk_score_table)
  risk_score_table <- cbind(risk_score_table, 'risk_score'=riskscore_value) %>%
    cbind(survival_cancer_df[,c('time','status')])
  risk_score_table <- risk_score_table[,c('time','status', candidate_genes_for_cox, 'risk_score')]
  return(risk_score_table)
}


####  data.frame中包含3列：time,status,risk_score
multi_ROC <- function(time_vector, risk_score_table){
        #install.packages('survivalROC')
        #install.packages("plotROC")
        library('survivalROC')
        library(plotROC)
        nobs <- nrow(risk_score_table)
        single_ROC <- function(single_time){
          for_ROC <- survivalROC(Stime = risk_score_table$time,
                                 status = risk_score_table$status,
                                 marker = risk_score_table$risk_score,
                                 predict.time = single_time,span = 0.25*nobs^(-0.20))
          data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
                     'Cut_values'=for_ROC$cut.values, 'Time_point'=rep(single_time, length(for_ROC$TP)),
                     'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
        }
        multi_ROC_list <- lapply(time_vector, single_ROC)
        do.call(rbind, multi_ROC_list)
}



sur_ROC_genes <- function(risk_score_table,pw=6,ph=6,output='output'){
          library(survivalROC)
          library(plotROC)
          library(ggplot2)
          nobs <- nrow(risk_score_table)
          nc<-ncol(risk_score_table)
          genes<-colnames(risk_score_table[,3:nc])
          single_ROC <- function(gene){
            for_ROC <- survivalROC(Stime = risk_score_table$time,
                                   status = risk_score_table$status,
                                   marker = risk_score_table[[gene]],
                                   predict.time = max(risk_score_table$time),span = 0.25*nobs^(-0.20))
            data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
                       'Cut_values'=for_ROC$cut.values, 'Gene'=rep(gene, length(for_ROC$TP)),
                       'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
          }
          multi_ROC_list <- lapply(genes, single_ROC)
          # do.call(rbind, multi_ROC_list)
          for_multi_ROC<-do.call(rbind, multi_ROC_list)
          AUC_max <- max(for_multi_ROC$AUC)
          AUC_max_time <- for_multi_ROC$Time_point[which(for_multi_ROC$AUC == AUC_max)]
          AUC_max_time <- AUC_max_time[!duplicated(AUC_max_time)]
          AUC_max_time <- AUC_max_time[length(AUC_max_time)]
          gene_auc<-unique(for_multi_ROC[,c('Gene','AUC')])
          head(for_multi_ROC)
          ncolor<-length(unique(for_multi_ROC$Gene))
          mycolors<-get_colors(ncolor)
          pROC<-ggplot(for_multi_ROC, aes(x = False_positive, y = True_positive, label = Cut_values, color = Gene)) + 
            geom_roc(labels = F, stat = 'identity', n.cuts = 0) + scale_color_manual(values=mycolors)+
            geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
            theme_bw()+
            theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
                  panel.grid = element_blank(),legend.position='none')
          for(i in 1:ncolor)
          {
            # i=4
            j=ncolor-i
            gene=gene_auc[i,1]
            auc=gene_auc[i,2]
            texts<-paste0("AUC of ", gene, ':')
            nc<-nchar(texts)
            nb<- 20-nc
            bks<-paste0(rep(' ',nb),sep='',collapse='')
            pROC=pROC+annotate("segment",x = 0.5,xend = 0.55, y = j*0.04+0.05, yend = j*0.04+0.05,colour = mycolors[i])
            pROC=pROC+annotate("text",x = 0.6, y = j*0.04+0.05,size=3,label = paste0(texts,bks, round(auc,2),sep=''),hjust=0)
          }
          ggsave(paste0(output,".tiff"),plot=pROC, width = pw, height = ph,compression='lzw')
          ggsave(paste0(output,".pdf"),plot=pROC, width = pw, height = ph)
          
          return(for_multi_ROC)
}
        

### input是data.frame，表头为project,gene_name,HR,HR.95L,HR.95H,pvalue,sample_num
### ggplot2_cox_forestplot(df_data,x_start=0,x_end=3.5,step=0.6,pw=9,ph=6,output='fig11.cox_forestplot')  
ggplot2_cox_forestplot<- function(df_data,geneid='miRNA',x_start=0,x_end=3.5,step=0.6,pw=9,ph=6,output='cox_forestplot')  ######  绘制hub基因的配对boxplot
{
  rownames(df_data)<-NULL
  df_data$id<-as.numeric(rownames(df_data))

  df_data[,"HR"]<-as.numeric(df_data[,"HR"])
  df_data[,"HR.95L"]<-as.numeric(df_data[,"HR.95L"])
  df_data[,"HR.95H"]<-as.numeric(df_data[,"HR.95H"])
  df_data[,"pvalue"]<-as.numeric(df_data[,"pvalue"])

  hr=sprintf("%.3f",df_data[,"HR"])
  hrLow=sprintf("%.3f",df_data[,"HR.95L"])
  hrHigh=sprintf("%.3f",df_data[,"HR.95H"])
  pVal=df_data$pvalue
  pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))

  # df_data$HR.95H[df_data$HR.95H>5] <- 5.1
  # df_data$HR.95L[df_data$HR.95L>5] <- 5.0
  # df_data$HR[df_data$HR>5] <- 5.0

  #df_data$HR.95H <- log10(df_data$HR.95H)
  #df_data$HR.95L <- log10(df_data$HR.95L)
  #df_data$HR <- log10(df_data$HR)

  len<-dim(df_data)[1]
  x1<-rep(x_start+0.3,len)
  x2<-rep(x_end-2*step,len)
  x3<-rep(x_end-1*step,len)
  x4<-rep(x_end,len)
  x<-c(x_start+0.3,1,x_end-2*step,x_end-1*step,x_end,x1,x2,x3,x4)
  y<-c(len+1,len+1,len+1,len+1,len+1,len:1,len:1,len:1,len:1)
  texts<-c(geneid,'Hazard Ratio','HR','95%CI','P Value',rev(df_data$gene_name),rev(hr),rev(paste0(hrLow,'-',hrHigh)),rev(pVal))
  annotation <- data.frame(x=x,y=y,label=texts)

  annotation1 <- data.frame(x=x[1:(5+len)],y=y[1:(5+len)],label=texts[1:(5+len)])
  annotation2 <- data.frame(x=x[(6+len):length(x)],y=y[(6+len):length(x)],label=texts[(6+len):length(x)])

  plota <- ggplot(df_data, aes(HR, id))
  plota <- plota+ geom_point(size=1.2, aes(col=gene_name)) +
    geom_segment(aes(x = 1, xend = 1, y = 0, yend = len+0.5),lwd = 0.8,linetype=6)+
    geom_segment(aes(x = x_start, xend = x_end+0.3, y = 0, yend = 0),lwd = 1.1)+
    geom_segment(aes(x = x_start, xend = x_end+0.3, y = len+0.5, yend = len+0.5),lwd = 1.1)+scale_y_discrete(expand = expansion(mult = c(0, 0.1)))+
    #geom_vline(aes(xintercept=1,lwd=0.3))+geom_hline(aes(yintercept=len+0.5,lwd=0.3))+
    geom_errorbarh(aes(xmax =HR.95H, xmin = HR.95L,col=gene_name),size=1,height=0.4)+             #height = 0.4
    scale_x_continuous(limits=c(x_start, x_end+0.3), breaks=seq(-3, 4, 1)) +
    #,labels=c('1e-03','1e-02','1e-01','1e+00','1e+01','1e+02','1e+03','1e+04')
    labs(y = "",x='') + theme_bw() + theme(legend.position ="none",panel.border=element_blank()) +
    theme(plot.margin=unit(c(1,1,1,1), 'lines'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text=element_text(size=10, face = "bold"),legend.text=element_text(size=11),
          axis.text.y=element_blank(),axis.ticks.y=element_blank())
  plota <- plota+geom_text(annotation1,mapping=aes(x=x,y=y,label=label),fontface="bold")
  plota <- plota+geom_text(annotation2,mapping=aes(x=x,y=y,label=label),size=5)
  ggsave(file=paste0(output,".tiff"),plot=plota,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".pdf"),plot=plota,width = pw, height = ph)
}


ggplot2_cox_forestplot_redblue<- function(df_data,geneid='miRNA',x_start=0,x_end=3.5,step=0.6,xstep=2,pw=9,ph=6,output='cox_forestplot',hjust=0.5)  ######  绘制hub基因的配对boxplot
{
  ngene<-nrow(df_data)
  rownames(df_data)<-NULL
  df_data$id<-ngene-as.numeric(rownames(df_data))+1
  df_data<-df_data[order(df_data$id),]
  
  df_data[,"HR"]<-as.numeric(df_data[,"HR"])
  df_data[,"HR.95L"]<-as.numeric(df_data[,"HR.95L"])
  df_data[,"HR.95H"]<-as.numeric(df_data[,"HR.95H"])
  df_data[,"pvalue"]<-as.numeric(df_data[,"pvalue"])
  
  hr=sprintf("%.3f",df_data[,"HR"])
  hrLow=sprintf("%.3f",df_data[,"HR.95L"])
  hrHigh=sprintf("%.3f",df_data[,"HR.95H"])
  pVal=df_data$pvalue
  pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
  
  # df_data$HR.95H[df_data$HR.95H>5] <- 5.1
  # df_data$HR.95L[df_data$HR.95L>5] <- 5.0
  # df_data$HR[df_data$HR>5] <- 5.0
  
  #df_data$HR.95H <- log10(df_data$HR.95H)
  #df_data$HR.95L <- log10(df_data$HR.95L)
  #df_data$HR <- log10(df_data$HR)
  
  len<-dim(df_data)[1]
  x1<-rep(x_start+0.3,len)
  x2<-rep(x_end-2*step,len)
  x3<-rep(x_end-1*step,len)
  x4<-rep(x_end,len)
  x<-c(x_start+0.3,1,x_end-2*step,x_end-1*step,x_end,x1,x2,x3,x4)
  y<-c(len+1,len+1,len+1,len+1,len+1,len:1,len:1,len:1,len:1)
  texts<-c(geneid,'Hazard Ratio','HR','95%CI','P Value',rev(df_data$gene_name),rev(hr),rev(paste0(hrLow,'-',hrHigh)),rev(pVal))
  annotation <- data.frame(x=x,y=y,label=texts)
  
  annotation1 <- data.frame(x=x[1:5],y=y[1:5],label=texts[1:5])
  annotation1b <- data.frame(x=x[6:(5+len)],y=y[6:(5+len)],label=texts[6:(5+len)])
  annotation2 <- data.frame(x=x[(6+len):length(x)],y=y[(6+len):length(x)],label=texts[(6+len):length(x)])
  
  mycolors<-get_colors(ngene)
  
  plota <- ggplot(df_data, aes(HR, id))
  plota <- plota+ geom_point(size=2, col='red') +
    geom_segment(aes(x = 1, xend = 1, y = 0, yend = len+0.5),lwd = 0.4,linetype=6)+
    geom_segment(aes(x = x_start, xend = x_end+0.3, y = 0, yend = 0),lwd = 1.1)+
    geom_segment(aes(x = x_start, xend = x_end+0.3, y = len+0.5, yend = len+0.5),lwd = 1.1)+scale_y_discrete(expand = expansion(mult = c(0, 0.1)))+
    #geom_vline(aes(xintercept=1,lwd=0.3))+geom_hline(aes(yintercept=len+0.5,lwd=0.3))+
    geom_errorbarh(aes(xmax =HR.95H, xmin = HR.95L),col='blue',size=0.5,height=0.4)+             #height = 0.4
    scale_x_continuous(limits=c(x_start, x_end+0.3), breaks=seq(1, floor(x_end), xstep)) +
    #scale_color_manual(values=mycolors)
    #,labels=c('1e-03','1e-02','1e-01','1e+00','1e+01','1e+02','1e+03','1e+04')
    labs(y = "",x='') + theme_bw() + theme(legend.position ="none",panel.border=element_blank()) +
    theme(plot.margin=unit(c(1,1,1,1), 'lines'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text=element_text(size=10, face = "bold"),legend.text=element_text(size=11),
          axis.text.y=element_blank(),axis.ticks.y=element_blank())
  plota <- plota+geom_text(annotation1,mapping=aes(x=x,y=y,label=label),fontface="bold")
  plota <- plota+geom_text(annotation1b,mapping=aes(x=x,y=y,label=label),fontface="bold",hjust=hjust)
  plota <- plota+geom_text(annotation2,mapping=aes(x=x,y=y,label=label),size=5)
  print(plota)
  ggsave(file=paste0(output,".tiff"),plot=plota,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".pdf"),plot=plota,width = pw, height = ph)
}


ggplot2_cox_forestplot_large<- function(df_data,geneid='miRNA',x_start=0,x_end=3.5,step=0.6,pw=9,ph=6,output='cox_forestplot')  ######  绘制hub基因的配对boxplot
{
  # df_data<-df_cox
  df_data$id<-as.numeric(rownames(df_data))
  
  df_data[,"HR"]<-as.numeric(df_data[,"HR"])
  df_data[,"HR.95L"]<-as.numeric(df_data[,"HR.95L"])
  df_data[,"HR.95H"]<-as.numeric(df_data[,"HR.95H"])
  df_data[,"pvalue"]<-as.numeric(df_data[,"pvalue"])
  
  hr= as.character(sprintf("%.3f",df_data[,"HR"]))
  hrLow= as.character(sprintf("%.3f",df_data[,"HR.95L"]))
  hrHigh=as.character(sprintf("%.3f",df_data[,"HR.95H"]))
  pVal=df_data$pvalue
  pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
  
  for(i in 1:nrow(df_data))
  {
    if(df_data[i,"HR.95L"]>1000)
    {
      hrLow[i]<- as.character(format(df_data[i,"HR.95L"],scientific=TRUE,digit=3))
    }
    if(df_data[i,"HR.95H"]>1000)
    {
      hrHigh[i]<-as.character(format(df_data[i,"HR.95H"],scientific=TRUE,digit=3))
    }
    if(df_data[i,"HR"]>1000)
    {
      hr[i]<-as.character(format(df_data[i,"HR"],scientific=TRUE,digit=3))
    }
  }
  
  
  df_data$HR.95H[df_data$HR.95H>5] <- 5.1
  df_data$HR.95L[df_data$HR.95L>5] <- 5.0
  df_data$HR[df_data$HR>5] <- 5.0
  
  
  
  
  
  #df_data$HR.95H <- log10(df_data$HR.95H)
  #df_data$HR.95L <- log10(df_data$HR.95L)
  #df_data$HR <- log10(df_data$HR)
  
  len<-dim(df_data)[1]
  x1<-rep(x_start+0.3,len)
  x2<-rep(x_end-2*step,len)
  x3<-rep(x_end-1*step,len)
  x4<-rep(x_end,len)
  x<-c(x_start+0.3,1,x_end-2*step,x_end-1*step,x_end,x1,x2,x3,x4)
  y<-c(len+1,len+1,len+1,len+1,len+1,len:1,len:1,len:1,len:1)
  texts<-c(geneid,'OS Hazard Ratio','HR','95%CI','P Value',rev(df_data$gene_name),rev(hr),rev(paste0(hrLow,'-',hrHigh)),rev(pVal))
  annotation <- data.frame(x=x,y=y,label=texts)
  
  annotation1 <- head(annotation,5)
  annotation3 <- annotation[6:(5+len),]
  annotation2 <- tail(annotation,length(x)-5-len)
  
  plota <- ggplot(df_data, aes(HR, id))
  plota <- plota+ geom_point(size=1.2, aes(col=gene_name)) +
    geom_segment(aes(x = 1, xend = 1, y = 0, yend = len+0.5),lwd = 0.8,linetype=6)+
    geom_segment(aes(x = x_start, xend = x_end+0.3, y = 0, yend = 0),lwd = 1.1)+
    geom_segment(aes(x = x_start, xend = x_end+0.3, y = len+0.5, yend = len+0.5),lwd = 1.1)+scale_y_discrete(expand = expansion(mult = c(0, 0.1)))+
    #geom_vline(aes(xintercept=1,lwd=0.3))+geom_hline(aes(yintercept=len+0.5,lwd=0.3))+
    geom_errorbarh(aes(xmax =HR.95H, xmin = HR.95L,col=gene_name),size=1,height=0.4)+             #height = 0.4
    scale_x_continuous(limits=c(x_start, x_end+0.3), breaks=seq(-3, x_end, 1)) +
    #,labels=c('1e-03','1e-02','1e-01','1e+00','1e+01','1e+02','1e+03','1e+04')
    labs(y = "",x='') + theme_bw() + theme(legend.position ="none",panel.border=element_blank()) +
    theme(plot.margin=unit(c(1,1,1,1), 'lines'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text=element_text(size=10, face = "bold"),legend.text=element_text(size=11),
          axis.text.y=element_blank(),axis.ticks.y=element_blank())
  plota <- plota+geom_text(annotation1,mapping=aes(x=x,y=y,label=label),fontface="bold",hjust=0.5)
  plota <- plota+geom_text(annotation3,mapping=aes(x=x,y=y,label=label),fontface="bold",hjust=0)
  plota <- plota+geom_text(annotation2,mapping=aes(x=x,y=y,label=label),size=4)
  ggsave(file=paste0(output,".tiff"),plot=plota,width = pw, height = ph,compression='lzw')
  ggsave(file=paste0(output,".pdf"),plot=plota,width = pw, height = ph)
}

##### ggplot2_forestplot与ggplot2_forestplot_log2的所有设置都可以一样，
##### zero应该都设置为1
##### x轴会自动选择xlim的范围
##### 当有较大的HR值时，选择ggplot2_forestplot_log2函数，此函数也只是将x轴按log2进行缩放，x_start始终要大于0
ggplot2_forestplot<- function(df_data,geneid='',zero=1,x_start=1,x_end=3.5,pw=9,ph=6,output='cox_forestplot')  ######  绘制hub基因的配对boxplot
{
	library(grid)
	library(forestplot)
	library(tidyverse)
  ###  df_data<-df_cox
  rownames(df_data)<-NULL
  df_data$id<-as.numeric(rownames(df_data))
  df_data[,"HR"]<-as.numeric(df_data[,"HR"])
  df_data[,"HR.95L"]<-as.numeric(df_data[,"HR.95L"])
  df_data[,"HR.95H"]<-as.numeric(df_data[,"HR.95H"])
  df_data[,"pvalue"]<-as.numeric(df_data[,"pvalue"])
  hr=sprintf("%.3f",df_data[,"HR"])
  hrLow=sprintf("%.3f",df_data[,"HR.95L"])
  hrHigh=sprintf("%.3f",df_data[,"HR.95H"])
  pVal=df_data$pvalue
  pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
  range<-paste0(hrLow,'-',hrHigh)

  df_text<-data.frame(gene_name=df_data[['gene_name']],
    hr=hr,range=range,pvalue=pVal)
  
  plota<-df_data |>
    forestplot(df_text, 
               mean = HR,
               lower = HR.95L,
               upper = HR.95H,
               zero=zero,
               graph.pos = 4,
               hrzl_lines = gpar(col = "#444444"),
               new_page = TRUE,
               clip = c(x_start,x_end), 
               xlog = F, 
               boxsize = 0.1,
               col = fpColors(box = "red",line = "blue"),
               vertices = TRUE,
               title = NA) |> 
    fp_add_lines(h_3 = gpar(lty = 2) ,
                 h_11 = gpar(lwd = 1, columns = 1:3, col = "#000044")) |>
    fp_add_header(gene_name = geneid |>
                    fp_align_center(),
                  hr = "Hazard Ratio"|>
                    fp_align_center(),
                  range = "95%CI"|>
                    fp_align_center(),
                  pvalue = "P Value" |>
                    fp_align_center())
  
  
  tiff(paste0(output,".tiff"),width=pw*300,height=ph*300,units="px",res=300,compression='lzw')
  plot(plota)
  dev.off()
  pdf(paste0(output,".pdf"),width=pw,height=ph)
  plot(plota)
  dev.off()

}


##### ggplot2_forestplot与ggplot2_forestplot_log2的所有设置都可以一样，
##### zero应该都设置为1
##### x轴会自动选择xlim的范围
##### 当有较大的HR值时，选择ggplot2_forestplot_log2函数，此函数也只是将x轴按log2进行缩放，x_start始终要大于0
ggplot2_forestplot_log2<- function(df_data,geneid='',zero=1,x_start=1,x_end=3,pw=9,ph=6,output='cox_forestplot')  ######  绘制hub基因的配对boxplot
{
	library(grid)
	library(forestplot)
	library(tidyverse)
  ###  df_data<-df_cox
  rownames(df_data)<-NULL
  df_data$id<-as.numeric(rownames(df_data))
  df_data[,"HR"]<-as.numeric(df_data[,"HR"])
  df_data[,"HR.95L"]<-as.numeric(df_data[,"HR.95L"])
  df_data[,"HR.95H"]<-as.numeric(df_data[,"HR.95H"])
  df_data[,"pvalue"]<-as.numeric(df_data[,"pvalue"])
  hr=sprintf("%.3f",df_data[,"HR"])
  hrLow=sprintf("%.3f",df_data[,"HR.95L"])
  hrHigh=sprintf("%.3f",df_data[,"HR.95H"])
  pVal=df_data$pvalue
  pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
  range<-paste0(hrLow,'-',hrHigh)
  
  
  
  df_text<-data.frame(gene_name=df_data[['gene_name']],
                      hr=hr,range=range,pvalue=pVal)
  
  df_data_used<-df_data[,c('HR','HR.95L','HR.95H')]
  df_data_used<- log10(df_data_used)
  plota<-df_data |>
    forestplot(df_text, 
               mean = HR,
               lower = HR.95L,
               upper = HR.95H,
               zero=zero,
               graph.pos = 4,
               hrzl_lines = gpar(col = "#444444"),
               new_page = TRUE,
               clip = c(x_start,x_end), 
               xlog = T, 
               boxsize = 0.1,
               col = fpColors(box = "red",line = "blue"),
               vertices = TRUE,
               title = NA) |> 
    fp_add_lines(h_3 = gpar(lty = 2) ,
                 h_11 = gpar(lwd = 1, columns = 1:3, col = "#000044")) |>
    fp_add_header(gene_name = geneid |>
                    fp_align_center(),
                  hr = "Hazard Ratio"|>
                    fp_align_center(),
                  range = "95%CI"|>
                    fp_align_center(),
                  pvalue = "P Value" |>
                    fp_align_center())
  
  
  tiff(paste0(output,".tiff"),width=pw*300,height=ph*300,units="px",res=300,compression='lzw')
  plot(plota)
  dev.off()
  pdf(paste0(output,".pdf"),width=pw,height=ph)
  plot(plota)
  dev.off()
  
}




#pv_train<-ggplot2_KM_plot(lassoSigExp_risk_score,plotid=plotid,projectid=projectid,sur='status',sur_time='time',group='group',title='TCGA-LIHC',legene_title='Risk_score',legene_labels=c('high','low'))

ggplot2_KM_plot<-function(mydata,projectid='TCGA',sur='OS_STATUS',sur_time='OS_MONTHS',group='group',title='TCGA',y_title='OS',
legene_title=NA,legene_labels=NULL,output='KM')
{
      library(survival)
      library(survminer)
      #mydata<-df_cox_used
      formula <- as.formula(paste0('Surv(',sur_time,', ',sur,')~', group, sep = '', collapse = ' + '))
      diff=survdiff(formula=formula,data = mydata)
      diff$call$formula <- formula
      pValue=1-pchisq(diff$chisq,df=1)
      print(pValue)
      pvalue_string<-pValue
      if(pValue<0.001){
        pvalue_string="Pvalue<0.001"
      }else{
        pvalue_string=paste0("Pvalue=",sprintf("%.03f",pValue))
      }
      cols<-c("blue","red",'green')
      ngroup<-length(unique(mydata[,group]))
      fit <- survfit(formula=formula, data = mydata)
      fit$call$formula <- formula
      surPlot=ggsurvplot(fit,
                         data=mydata,
                         title=title,
                         pval=pvalue_string,
                         pval.size=6,
                         legend.labs=legene_labels,
                         legend.title=NULL,
                         font.legend=12,
                         xlab="Time(Months)",
                         ylab=y_title,
                         break.time.by = 12,
                         palette=cols[1:ngroup],
                         conf.int=F,
                         fontsize=4,
                         risk.table=F,
                         risk.table.title="",
                         risk.table.height=.25)
      #tiff(file=paste0(sur,".",i,".tiff"),width = 6, height =5)
      tiff(paste0(output,".a.tiff"),width=2400,height=2000,units="px",res=300,compression='lzw')
      print(surPlot)
      dev.off()
      pdf(file=paste0(output,".a.pdf"),width=8,height =6.6)
      print(surPlot)
      dev.off()
      
      surPlot=ggsurvplot(fit,
                         data=mydata,
                         title=title,
                         pval=pvalue_string,
                         pval.size=6,
                         legend.labs=legene_labels,
                         legend.title=NULL,
                         font.legend=12,
                         xlab="Time(Months)",
                         ylab=y_title,
                         break.time.by = 12,
                         palette=cols[1:ngroup],
                         conf.int=F,
                         fontsize=4,
                         risk.table=T,
                         risk.table.title="",
                         risk.table.height=.25)
      #tiff(file=paste0(sur,".",i,".tiff"),width = 6, height =5)
      tiff(paste0(output,".b.tiff"),width=2400,height=2400,units="px",res=300,compression='lzw')
      print(surPlot)
      dev.off()
      pdf(file=paste0(output,".b.pdf"),width=8,height =8)
      print(surPlot)
      dev.off()
      return(pValue)
    }


tcga_clinical_format <-function(df_clinical)
{
      df_clinical[which(df_clinical$OS_STATUS=='0:LIVING'),'OS_STATUS']<-'0'
      df_clinical[which(df_clinical$OS_STATUS=='1:DECEASED'),'OS_STATUS']<-'1'
      df_clinical$OS_STATUS<-as.numeric(df_clinical$OS_STATUS)
      df_clinical$OS_MONTHS<-as.numeric(df_clinical$OS_MONTHS)
      
      df_clinical[which(df_clinical$DFS_STATUS=='0:DiseaseFree'),'DFS_STATUS']<-'0'
      df_clinical[which(df_clinical$DFS_STATUS=='1:Recurred/Progressed'),'DFS_STATUS']<-'1'
      df_clinical$DFS_STATUS<-as.numeric(df_clinical$DFS_STATUS)
      df_clinical$DFS_MONTHS<-as.numeric(df_clinical$DFS_MONTHS)
      
      
      df_clinical<-df_clinical[!duplicated(df_clinical$patientId),]
      dim(df_clinical)
      rownames(df_clinical)<-df_clinical$patientId
      return(df_clinical)
}

tcga_tumor_sample<-function(df_logtpm,df_sample,group='Tumor')
{
  tumor_samples<-row.names(df_sample[df_sample$group==group,,drop=F])
  df_tumor_logtpm<-df_logtpm[,tumor_samples]
  return(df_tumor_logtpm)
}
      
      
lasso_cox_validation<-function(df_cox_data,lassoGene,lasso_fit,cut_off,plotid=26,sur_time='OS_MONTHS',sur='OS_STATUS')
{
  
  df_cox_data[1:3,c(1:3,ncol(df_cox_data))]
  df_cox_data_used<-df_cox_data[,c(sur_time,sur,lassoGene)]
  colnames(df_cox_data_used)[1:2]<-c('time','status')
  total_risk_score <- riskscore_lasso(df_cox_data_used, lassoGene, lasso_fit)
  head(total_risk_score)
  
  for_multi_ROC <- multi_ROC(time_vector = c(12*seq(1,10,2)), risk_score_table = total_risk_score)
  AUC_max <- max(for_multi_ROC$AUC)
  #maybe AUCs are identical in different time points. So select the last time point indicating longer survival.
  AUC_max_time <- for_multi_ROC$Time_point[which(for_multi_ROC$AUC == AUC_max)]
  AUC_max_time <- AUC_max_time[!duplicated(AUC_max_time)]
  AUC_max_time <- AUC_max_time[length(AUC_max_time)]
  for_multi_ROC$Time_point_factor <- as.factor(for_multi_ROC$Time_point/12)
  
  head(for_multi_ROC)
  unique(for_multi_ROC$Time_point_factor)
  ncolor<-length(unique(for_multi_ROC$Time_point_factor))
  mycolors<-get_colors(ncolor)
  #visualization of the ROC curves of multiple time points.
  pROC<-ggplot(for_multi_ROC, aes(x = False_positive, y = True_positive, label = Cut_values, color = Time_point_factor)) + 
	geom_roc(labels = F, stat = 'identity', n.cuts = 0) + scale_color_manual(values=mycolors)+
	geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
	theme_bw()+
	theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
		  panel.grid = element_blank())+
	annotate("text",x = 0.75, y = 0.15,
			 label = paste("AUC max = ", round(AUC_max, 2), '\n', 'AUC max time = ', AUC_max_time/12, ' years', sep = ''))
  #pROC
  ggsave(paste0("fig",plotid,".",projectid,".pROC.tiff"),plot=pROC, width = 8, height = 8,compression='lzw')
  ggsave(paste0("fig",plotid,".",projectid,".pROC.pdf"),plot=pROC, width = 8, height = 8)
  
  total_risk_score$group<-'high'
  total_risk_score[total_risk_score$risk_score <= cut_off,'group']<-'low'
  write.table(total_risk_score,file=paste0('fig',plotid,'.',projectid,'.risk_score_group.tsv'),row.names = T,sep='\t',quote = FALSE)
  pv_test<-ggplot2_KM_plot(total_risk_score,projectid=projectid,
						   sur='status',sur_time='time',group='group',title=NULL,legene_title='Risk_score',legene_labels=c('high','low'),output=paste0('fig',plotid,'.',projectid,'.test.KM'))
  
  
  
  if(1)   ######  绘制差异表达基因的热图（PCG和LNC分开绘制）
  {
	phenotype<-data.frame(group=total_risk_score$group)
	head(phenotype)
	row.names(phenotype)<-row.names(total_risk_score)
	
	#lassoSigExp_heatmap$id<-NULL
	head(df_cox_data)
	dim(df_cox_data)
	data_heatmap<-t(df_cox_data[,lassoGene])
	ggheatmap_yingbio_2groups(data_heatmap,phenotype,show_rowname=T,pheight=2,
							  geneid='gene',output=paste0('fig',plotid,'.',projectid,'.test.hubgene.heatmap'))
  }
  return(pv_test)
}

### sur_nomogram_regplot(cox_fit,months=c(36,60),output=paste0('fig',plotid,'.',projectid,'.train.nomogram'))
sur_nomogram_regplot<-function(cox_fit,months=c(36,60),output='output')
{
library(regplot)
regplot(cox_fit,plots=c("density","boxes"),observation=FALSE,
		failtime = months,#预测3年和5年的死亡风险，此处单位是mouth
		prfail = TRUE, #cox回归中需要TRUE
		points=TRUE,showP = T,dencol="yellow",boxcol="yellow") #是否展示统计学差异
dev.copy2pdf(file=paste0(output,".pdf"), width=8, height=6, out.type="pdf")
dev.off()
}
      
### sur_nomogram_rms(cox_fit,months=c(36,60),output=paste0('fig',plotid,'.',projectid,'.train.nomogram.rms'))
sur_nomogram_rms<-function(cox_fit,months=c(36,60),output='output',pw=8,ph=6)
{
library(rms)
surv1 <- function(x) surv(36,x) # 3年OS
surv2 <- function(x) surv(60,x) # 5年OS
## nom <- nomogram(cox_fit, fun=function(x) med(lp=x),funlabel="Median Survival Time")
nom <- nomogram(cox_fit,fun = list(surv1,surv2),lp = T,
				funlabel = c('3-year OS','5-year OS'),
				maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
pdf(file=paste0(output,".pdf"),width=pw,height=ph)
plot(nom, 
	 lplabel="Linear Predictor",
	 xfrac = 0.2, # 左侧标签距离坐标轴的距离
	 #varname.label = TRUE, 
	 tcl = -0.2, # 刻度长短和方向 
	 lmgp = 0.1, # 坐标轴标签距离坐标轴远近
	 points.label ='Points', 
	 total.points.label = 'Total Points',
	 cap.labels = FALSE,
	 cex.var = 1, # 左侧标签字体大小
	 cex.axis = 1, # 坐标轴字体大小
	 col.grid = gray(c(0.8, 0.95))) # 竖线颜色
dev.off()
tiff(file=paste0(output,".tiff"),width=pw*300,height=ph*300,compression='lzw',res=300)
plot(nom, 
	 lplabel="Linear Predictor",
	 xfrac = 0.2, # 左侧标签距离坐标轴的距离
	 #varname.label = TRUE, 
	 tcl = -0.2, # 刻度长短和方向 
	 lmgp = 0.1, # 坐标轴标签距离坐标轴远近
	 points.label ='Points', 
	 total.points.label = 'Total Points',
	 cap.labels = FALSE,
	 cex.var = 1, # 左侧标签字体大小
	 cex.axis = 1, # 坐标轴字体大小
	 col.grid = gray(c(0.8, 0.95))) # 竖线颜色
dev.off()
}


wilcox.test_deg<-function(df_count,df_normalized,phenotype,group.by='group',groups=c('URG_cluster1','URG_cluster2'),logFC_cutoff=1,padj_cutoff=0.05,plotid=20)
  {
    comparison=paste0(groups[2],'_vs_',groups[1])
    group_list<-phenotype[,group.by]
    group_list<-factor(group_list,levels=groups)
    
    if(1)   ######  输出差异表达基因表格
    {
      group_mean=apply(df_count,1,function(x) aggregate(x,by=list(group_list),mean)$x)
      group_mean<-t(group_mean)
      colnames(group_mean)<-groups
      head(group_mean)

      # Run the Wilcoxon rank-sum test for each gene
      pvalues <- sapply(1:nrow(df_normalized),function(i){
        data<-cbind.data.frame(gene=as.numeric(t(df_normalized[i,])),group_list)
        p=wilcox.test(gene~group_list, data)$p.value
        return(p)
      })
      fdr=p.adjust(pvalues,method = "fdr")

      # Calculate fold-change for each gene
      dataCon1=df_normalized[,c(which(group_list==groups[1]))]
      dataCon2=df_normalized[,c(which(group_list==groups[2]))]
      foldChanges=rowMeans(dataCon2)-rowMeans(dataCon1)
      # Output results based on FDR threshold
      res<-data.frame(log2foldChange=foldChanges, pValues=pvalues, FDR=fdr)
      
      
      rownames(res)=rownames(df_count)
      head(res)
      res$change = as.factor(
        ifelse(res$FDR < padj_cutoff & abs(res$log2foldChange) > logFC_cutoff,
               ifelse(res$log2foldChange > logFC_cutoff ,'UP','DOWN'),'NOT'))
      DEG_full<-merge(res,group_mean,by='row.names')
      colnames(DEG_full)[1]<-'gene'
      write.table(DEG_full,file=paste0('fig',plotid,'a.tpm.wilcox.',comparison,'.deg.all.tsv'),row.names = F,sep='\t',quote = FALSE)
      DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]
      write.table(DEG_sig,file=paste0('fig',plotid,'a.tpm.wilcox.',comparison,'.deg.sig.tsv'),row.names = F,sep='\t',quote = FALSE)
      dim(DEG_sig)

    }
    if(1) ########## 火山图 使用自编写的函数绘制
    {
      head(DEG_full)
      ggplot2_volcano(DEG_full[,c('FDR','log2foldChange')],comparison=comparison,
                      output=paste0('fig',plotid,'b.tpm.wilcox.volcano.',comparison),fc_threshold=2^logFC_cutoff,ylab='FDR',ymax=26)
    }
    if(1)   ######  绘制差异表达基因的热图
    {
      df_expr_deg<-df_normalized[rownames(df_normalized) %in% DEG_sig$gene,]
      ggheatmap_yingbio_2groups(df_expr_deg,phenotype,geneid='gene',pheight=8,
                                output=paste0('fig',plotid,'c.tpm.wilcox.heatmap.deg.',comparison),save.data=T)
    }

    if(0)  ######  进行GO和KEGG pathway分析
    {
      library(stats)
      #DEG_full <- read.table(file = 'TCGA-KIRC.count.DESeq2.deg.all.tsv',header=T,sep='\t',check.names=FALSE)
      #DEG_sig <- read.table(file = 'TCGA-KIRC.count.DESeq2.deg.sig.tsv',header=T,sep='\t',check.names=FALSE)
      #comparison='Tumor_vs_Normal'
      go_analysis=T
      head(DEG_sig)
      #plotid<-'17'
      if(go_analysis)
      {
        #volcano_ggplot2(df_deg_all[,c('P.Value','logFC')],comparison=paste0(comparison),plotid=plotid,threshold=2,ylab='P.value',ymax=11)
        source('g:/100sop/sop_perl/001Seurat.functions.panzhong.R')
        aaa<-DEG_sig$gene[DEG_sig$log2foldChange>0]
        aaa
        gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UP'),paste0('kegg_',comparison,'_UP')),species = 'human')

        aaa<-DEG_sig$gene[DEG_sig$log2foldChange<0]
        gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_DOWN'),paste0('kegg_',comparison,'_DOWN')),species = 'human')

        aaa<-DEG_sig$gene
        length(aaa)
        #gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UPDOWN'),paste0('kegg_',comparison,'_UPDOWN')),species = 'human')

        #gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UPDOWN_nofilter'),paste0('kegg_',comparison,'_UPDOWN_nofilter')),species = 'human',go_pvalueCutoff=1,kegg_pvalueCutoff=1)


        #aaa<-DEG_full[,c('gene','avg_log2FC')]
        #gsea_yingbai(aaa,dir_output=paste0('gsea_',comparison),species = "human")
      }
    }
  }

wilcox.test_deg_reFilter<-function(DEG_full,df_count,df_expr,phenotype,group.by='group',groups=c('URG_cluster1','URG_cluster2'),logFC_cutoff=1,padj_cutoff=0.05,plotid=20)
{
  comparison=paste0(groups[2],'_vs_',groups[1])
  #DEG_full<-read.delim(paste0('fig1a.tpm.wilcox.',comparison,'.deg.all.tsv'))
  head(DEG_full)
  DEG_full$change = as.factor(
	ifelse(DEG_full$FDR < padj_cutoff & abs(DEG_full$log2foldChange) > logFC_cutoff,
		   ifelse(DEG_full$log2foldChange > logFC_cutoff ,'UP','DOWN'),'NOT'))
	DEG_sig<-DEG_full[DEG_full$change %in% c('UP','DOWN'),]
      write.table(DEG_sig,file=paste0('fig',plotid,'a.tpm.wilcox.',comparison,'.deg.sig.tsv'),row.names = F,sep='\t',quote = FALSE)

  if(1) ########## 火山图 使用自编写的函数绘制
  {
	head(DEG_full)
	ggplot2_volcano(DEG_full[,c('FDR','log2foldChange')],comparison=comparison,
					output=paste0('fig',plotid,'b.tpm.wilcox.volcano.',comparison),fc_threshold=2^logFC_cutoff,ylab='FDR',ymax=26)
  }
  if(1)   ######  绘制差异表达基因的热图
  {
	df_expr_deg<-df_expr[rownames(df_expr) %in% DEG_sig$gene,]
	ggheatmap_yingbio_2groups(df_expr_deg,phenotype,geneid='gene',pheight=8,
							  output=paste0('fig',plotid,'c.tpm.wilcox.heatmap.deg.',comparison),save.data=T)
  }

  if(0)  ######  进行GO和KEGG pathway分析
  {
	library(stats)
	#DEG_full <- read.table(file = 'TCGA-KIRC.count.DESeq2.deg.all.tsv',header=T,sep='\t',check.names=FALSE)
	#DEG_sig <- read.table(file = 'TCGA-KIRC.count.DESeq2.deg.sig.tsv',header=T,sep='\t',check.names=FALSE)
	#comparison='Tumor_vs_Normal'
	go_analysis=T
	head(DEG_sig)
	plotid<-'17'
	if(go_analysis)
	{
	  #volcano_ggplot2(df_deg_all[,c('P.Value','logFC')],comparison=paste0(comparison),plotid=plotid,threshold=2,ylab='P.value',ymax=11)
	  source('g:/100sop/sop_perl/001Seurat.functions.panzhong.R')
	  aaa<-DEG_sig$gene[DEG_sig$log2foldChange>0]
	  aaa
	  gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UP'),paste0('kegg_',comparison,'_UP')),species = 'human')

	  aaa<-DEG_sig$gene[DEG_sig$log2foldChange<0]
	  gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_DOWN'),paste0('kegg_',comparison,'_DOWN')),species = 'human')

	  aaa<-DEG_sig$gene
	  length(aaa)
	  #gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UPDOWN'),paste0('kegg_',comparison,'_UPDOWN')),species = 'human')

	  #gokegg_yingbai(aaa,dir_output=c(paste0('go_',comparison,'_UPDOWN_nofilter'),paste0('kegg_',comparison,'_UPDOWN_nofilter')),species = 'human',go_pvalueCutoff=1,kegg_pvalueCutoff=1)


	  #aaa<-DEG_full[,c('gene','avg_log2FC')]
	  #gsea_yingbai(aaa,dir_output=paste0('gsea_',comparison),species = "human")
	}
  }
}


matrix_probeid2symbol<-function(mt_data,ids){
  colnames(ids)<-c('probe_id','symbol')
  ids<-ids[,c('probe_id','symbol')]
  head(ids)
  ids=ids[ids$symbol != '',]
  ids=ids[ids$probe_id %in%  rownames(mt_data),]
  class(ids$probe_id)
  ids$probe_id<-as.character(ids$probe_id)
  
  mt_data=mt_data[ids$probe_id,] 
  ids$median=apply(mt_data,1,median)
  ids=ids[order(ids$symbol,ids$median,decreasing = T),]
  ids=ids[!duplicated(ids$symbol),]
  mt_data=mt_data[ids$probe_id,]
  rownames(mt_data)=ids$symbol
  return(mt_data)
}


###  df_sur_data<-tcga_matrix_sampleid_format(df_sur_data,length=12)
#### 行名是tcga样本名，列是基因或者别的数据
tcga_matrix_sampleid_format<-function(df_data,length=15)
{
  df_data<-as.data.frame(df_data)
  rownames(df_data)<-gsub('[.]','-',rownames(df_data))
  df_data$patientId<-substring(rownames(df_data),1,length)
  df_data<-df_data[!duplicated(df_data$patientId),]
  dim(df_data)
  rownames(df_data)<-df_data$patientId
  df_data$patientId<-NULL
  return(df_data)
}

### 数据包含3列：status、time、riskscore
### risk_score_table<-df_sur_data_os[,c('OS_STATUS','OS_MONTHS',genei)]
### ggplot2_riskscore_status(risk_score_table,output=paste0('fig',plotid,'.riskscore_status.',genei))
ggplot2_riskscore_status<-function(df_risk_score,y_title='Risk score',cutoff=0.5,labels=c('low','high'),output='output')
{
  library(survival)
  risk_score_data<-df_risk_score
  colnames(risk_score_data)<- c('OS_STATUS','OS_MONTHS','risk_score','group')
  risk_score_data <- risk_score_data[order(risk_score_data[,'risk_score']),]
  risk_score_data$patientid <- 1:length(risk_score_data$group)
  risk_score_data$OS_STATUS <- ifelse(risk_score_data$OS_STATUS == 1,"Dead","Alive")
  risk_score_data$group<-factor(risk_score_data$group,levels=labels)
  
  mycolors<-c('blue','red')
  ###第一个图
  library(ggplot2)
  p1=ggplot(risk_score_data,aes(x=patientid,y=risk_score))+geom_point(aes(color=group))+
    theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                     panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
                     plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                     legend.position=c(0.1,0.9),legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                     legend.key.height=unit(0.2,"line"),legend.background=element_blank(),axis.title.x=element_blank(),
                     axis.text.x=element_text(size=10,color='black',angle=45,hjust=1))+
    labs(y=y_title)+scale_color_manual(values=mycolors)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(risk_score_data$group==labels[1]),colour="black", linetype="dotted",size=0.8)
  
  #第二个图
  p2=ggplot(risk_score_data,aes(x=patientid,y=OS_MONTHS))+geom_point(aes(col=OS_STATUS))+
    theme_bw()+theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
                     panel.background = element_blank(),panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
                     plot.margin=unit(c(1,1,0.5,1.5), 'lines'),
                     legend.position='bottom',legend.title=element_blank(),legend.spacing=unit(0.01, 'line'),
                     legend.key.height=unit(0.2,"line"),legend.background=element_blank(),
                     axis.text.x=element_text(size=10,color='black',angle=45,hjust=1))+
    labs(x="Patient ID(increasing risk score)",y="Survival time(months)")+scale_color_manual(values=mycolors)+
    geom_vline(xintercept=sum(risk_score_data$group==labels[1]),colour="black", linetype="dotted",size=0.8)
  
  #拼图实现三图联动
  library(ggplotify)
  plots = list(p1,p2)
  library(gridExtra)
  lay1 = rbind(c(rep(1,6)),c(rep(2,6))) #布局矩阵
  pdf(file=paste0(output,".pdf"), height=5, width=5)
  print(grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,2),weights=c(10,10)))
  dev.off()
  tiff(file=paste0(output,".tiff"), height=5*300, width=5*300,compression='lzw',res=300)
  print(grid.arrange(grobs = plots, layout_matrix = lay1, heigths = c(2,2),weights=c(10,10)))
  dev.off()
}

####  df_risk_score的列必须为3列：status、time、risk_score
ggplot2_survROC<-function(df_plot_roc,time_vector=c(12,36,60),output='output')
{
  ### df_plot_roc<-df_survROC
  head(df_plot_roc)
  library(ggsci)
  library("scales")
  mycolors= pal_lancet()(10)  
  show_col(mycolors)
  mycolors<-mycolors[1:length(time_vector)]
  for_multi_ROC <- multi_ROC(time_vector = time_vector, risk_score_table = df_plot_roc)
  AUC_max <- max(for_multi_ROC$AUC)
  #maybe AUCs are identical in different time points. So select the last time point indicating longer survival.
  AUC_max_time <- for_multi_ROC$Time_point[which(for_multi_ROC$AUC == AUC_max)]
  AUC_max_time <- AUC_max_time[!duplicated(AUC_max_time)]
  AUC_max_time <- AUC_max_time[length(AUC_max_time)]
  for_multi_ROC$Time_point_factor <- as.factor(for_multi_ROC$Time_point/12)
  head(for_multi_ROC)
  unique(for_multi_ROC$Time_point_factor)
  #visualization of the ROC curves of multiple time points.
  pROC<-ggplot(for_multi_ROC, aes(x = False_positive, y = True_positive, label = Cut_values, color = Time_point_factor)) + 
    geom_roc(labels = F, stat = 'identity', n.cuts = 0) + 
    geom_abline(slope = 1, intercept = 0, color = 'grey', linetype = 2)+
    theme_bw()+scale_color_manual(values=mycolors,name='Years')+
    theme(text=element_text(size=20),
          panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
          panel.grid = element_blank(),legend.position=c(0.9,0.2),
          legend.key = element_blank(),legend.background=element_blank(),
          legend.text = element_text(size = 20))+
    annotate("text",x = 0.25, y = 0.95,
             label = paste("AUC max = ", round(AUC_max, 2), '\n', 'AUC max time = ', AUC_max_time/12, ' years', sep = ''),size=8)
  ggsave(paste0(output,".tiff"),plot=pROC, width = 8, height = 8,compression='lzw')
  ggsave(paste0(output,".pdf"),plot=pROC, width = 8, height = 8)
}
####  df_risk_score的列必须为3列：status、time、risk_score
ggplot2_KM_ROC<-function(df_risk_score,y_title='OS',output='output')
{
  
  cutoff<-median(df_risk_score[[genei]])
  df_risk_score$group<-'high'
  df_risk_score[df_risk_score[[genei]] <= cutoff,'group']<-'low'
  head(df_risk_score)
  
  library(survival)
  df_km<-df_risk_score
  colnames(df_km)<- c('status','time','risk_score','group')
  pv_train<-ggplot2_KM_plot(df_km,output=paste0(output,'.KM'),sur='status',sur_time='time',group='group',
                            title='TCGA-LIHC',y_title=y_title,legene_title='Risk_score',legene_labels=c('high','low'))
  df_plot_roc<-df_risk_score
  head(df_plot_roc)
  colnames(df_plot_roc)<- c('status','time','risk_score','group')
  ggplot2_survROC(df_plot_roc,time_vector=c(12,36,60),output=paste0(output,'.ROC'))
}


ggplot2_roc_specificity_plot<-function(obj_roc,title='Train',output='output.roc'){
library(pROC)
library(ggplot2)
auc_train <- round(obj_roc$auc[1],4)
print(auc_train)
text_add<-paste0('AUC = ', auc_train)
plota<-NULL
plota<-ggroc(obj_roc, colour = 'red', size = 1) +
  ggtitle(title) +
  #theme_bw()+
  theme(text=element_text(size=20),panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
		panel.border=element_rect(linetype='solid',fill=NA,colour ="black"),
		panel.background = element_blank(),plot.title = element_text(hjust = 0.5,size=20))
plota<-plota+annotate('text',x=0.5,y=0.9, label=text_add,size = 5)
ggsave(file=paste0(output,".tiff"),plot=plota,width = 6, height = 6,compression='lzw')
ggsave(file=paste0(output,".pdf"),plot=plota,width = 6, height = 6)
}
      
      
