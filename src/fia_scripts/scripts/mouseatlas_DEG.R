setwd("/Volumes/RBMB/MouseAtlas/")
#.libPaths("/Volumes/RBMB/library/")
library(limma)
library(nlme)
library(gplots)
library(dplyr)
library(ggrepel)
library(GSVA)
library(pheatmap)
library(ggplot2)
library(data.table)
library(ggpubr)
library(edgeR)
library(DESeq2)
library(glue)
library(ggfortify)
library(Seurat)

exp=read.delim("data/deg/newFia/pseudobulk_percelltypes_130125_kylemethod_smoke_combined.tsv", row.names=1)
#colnames(exp)=gsub("^(X)","",colnames(exp)) #remove string that's only at the beginning of the character
colnames(exp)=gsub("[.]"," ",colnames(exp))
colnames(exp)=gsub("| |","",colnames(exp))
colnames(exp)=gsub('2 W','2W',colnames(exp))
colnames(exp)=gsub('4 W','4W',colnames(exp))
colnames(exp)=gsub('6 W','6W',colnames(exp))
colnames(exp)=gsub('8 W','8W',colnames(exp))
colnames(exp)=gsub('k C','kC',colnames(exp))
colnames(exp)=gsub('k A','kA',colnames(exp))

celltypes=strsplit(colnames(exp),"_");celltypes=sapply(celltypes,`[`,1)
celltypes=unique(celltypes)
meta=read.csv("data/deg/newFia/meta_pseudobulk_percelltypes_130125_smokecombined (1).csv")
#a=meta$sex;a[grep("Male",a)]="Male"
#meta$sex=a
meta$combined_category=gsub("[.]|/|\\^| |[+]|[-]","",meta$combined_category)
#meta=meta[meta$study!="Mouse_Post_Sendai",]
#meta$age=gsub("18-20","19",meta$age);meta$age=gsub("16-20","18",meta$age);meta$age=gsub("12-13","13",meta$age);meta$age=gsub("10-12","11",meta$age)
#meta$age=gsub("16-36","26",meta$age);meta$age=gsub("Postnatal Day 7","1",meta$age);meta$age=gsub("6-20","13",meta$age);meta$age=gsub("6-8","7",meta$age);
#meta$age=gsub("5-6","6",meta$age);meta$age=gsub("\\b10\\b","14",meta$age) #use \\b to limit to only that string
#meta$age=as.numeric(meta$age)
#meta$sex=gsub('0','Female',meta$sex);meta$sex=gsub('1','Male',meta$sex)
#meta$strain[meta$study == 'Mouse_Covid'] <- 'K18-ACE2(2Prlmn/J)'
#meta$strain[meta$study == 'Mouse_Fibro_Age'] <- 'C57BL6/J-rj'
#meta$platform[meta$study == 'Mouse_Covid'] <- 'Rhapsody'
#meta$platform[meta$study == 'Mouse_Fibro_Age'] <- 'Drop-seq'
meta=meta[meta$combined_category!="Removed",]
#write.csv(meta,'data/deg/newFia/meta_pseudobulk_percelltypes_0611.csv')

exp=exp[,colnames(exp) %in% meta$combined_category]
exp=as.matrix(exp)
class(exp)<-'integer'


#colnames(meta)[1]="age"
#meta[,'sample']=c(rep("COPD",3),rep("Control",3),"COPD")


#exp=exp[,meta$Group]
#exp[, 1:7] <- sapply(exp[, 1:7], as.integer)

#meta$sample=gsub(" ","_",meta$sample)


##Code for plotting#####
plotting<-function(x,y,label="comparisons",ct="celltype"){
  temp_plot= ggplot(x, aes(x = log2FoldChange, y = -log10(pvalue))) + # This line creates the axes plus background of the graph
    geom_point(aes(color = Legend)) + # This line plus previous coding line creates the volcano plot on to the graph with values coloured in three groups
    geom_hline(yintercept =-log10(max(y$pvalue)),colour="#990000", linetype="dashed")+
    geom_vline(xintercept =-log2(2),colour="#990000", linetype="dashed")+
    geom_vline(xintercept =log2(2),colour="#990000", linetype="dashed")+
    # Creates a line at the intercept which shows significance
    theme_bw(base_size = 12) + theme(legend.position = "bottom",plot.title = element_text(hjust=0.5,face="bold")) + # Changes the legend to below the graph
    labs(title=paste0(ct," ",label))+
    geom_text_repel(data = top20,
                    aes(label =hgnc_symbol),size = 3, box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"),
                    max.overlaps = 30)     
  
  if(length(unique(x$Legend))>2){
    p=temp_plot+scale_color_manual(values = c("blue", "red","grey"))
  }else if(length(unique(x$Legend)) == 2 &&
           all(c("Decreased", "Not Sig") %in% unique(x$Legend))){
    p=temp_plot+scale_color_manual(values = c("blue", "grey"))
  }else if(length(unique(x$Legend)) == 2 &&
           all(c("Increased", "Not Sig") %in% unique(x$Legend))){
    p=temp_plot+scale_color_manual(values = c("red", "grey"))
  }else if(length(unique(x$Legend))==1){
    p=temp_plot+scale_colour_manual(values=c("grey"))
  }
  
  return(p)
}
#sex=as.data.frame(cbind(t(exp[c('Xist','Eif2s3y'),]),meta$sex))
#meta$sex=ifelse(sex$Xist<sex$Eif2s3y,"Male","Female")
#meta$sex=as.factor(meta$sex)

#Differential exp for each cell type #####
###sex#####
deg=list()
for(c in 1:length(celltypes)){
  
  skip_to_next<-FALSE
  
  diffexp<-function(x,y){
    
    exp2=x[,startsWith(colnames(x),celltypes[c])]
    meta2=y[startsWith(y$combined_category,celltypes[c]),]
    
    dds <- DESeqDataSetFromMatrix(countData = exp2,
                                  colData = meta2,
                                  design= ~sex)
    
    keep<-matrix(ncol=1,nrow=nrow(exp2));rownames(keep)=rownames(exp2)
    for (i in 1:nrow(exp2)){
      table=exp2[i,]
      if(table[order(table)][round((ncol(exp2)/4)*3)]>=1){
        keep[i,]=TRUE
      }else{keep[i,]=FALSE}
      
    }
    #row.names(keep)=row.names(exp2)
    keep2=as.logical(keep)
    
    #keep <- rowMedians(counts(dds)) >= 10
    dds <- dds[keep2,]
    dds <- DESeq(dds)
    resultsNames(dds)
    
    res <- results(dds, contrast=c("sex","Male","Female"))
    #res$log2FoldChange=res$log2FoldChange*-1
    res=as.data.frame(res)
    tT=res
    tT$padj=p.adjust(tT$pvalue,method="BH")
    tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(1.5), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(1.5), "Increased", "Not Sig"))
    tT=na.omit(tT)
    #write.csv(tT,paste0("deg/sex_",celltypes[c],".csv"))
    deg[[c]]<<-tT
    names(deg)[c]<<-celltypes[c]
    
    
    tT2=tT[which((tT$log2FoldChange>log2(1.5)|tT$log2FoldChange<(-log2(1.5)))&tT$padj<0.05),]
    tT2=tT2[order(tT2$padj),]
    tT2$hgnc_symbol=row.names(tT2)
    
    top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:20,]
    top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:20,]
    top20<<-rbind(top20up,top20down)
    
    #p=plotting(tT,tT2, label="Male v Female",ct=celltypes[c])
    
    #ggsave(paste0("plots/sex_",celltypes[c],".png"),plot=p,width=2000,height = 1500,units="px",dpi=300);print(paste0(celltypes[c]," saved"))
  }
  
  tryCatch( diffexp(exp,meta), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
}

openxlsx::write.xlsx(deg,"deg/new2024/covariates/sex_DEGs_eachcelltypes.xlsx",rowNames=T)


for(c in celltypes){
  print(paste0(c," number of samples = ",ncol(exp[,startsWith(colnames(exp),c)])))
}


###age#####
deg=list()
for(c in 1:length(celltypes)){
  
  skip_to_next<-FALSE
  
  diffexp<-function(x,y){
    exp2=x[,startsWith(colnames(x),celltypes[c])]
    meta2=y[startsWith(y$combined_category,celltypes[c]),]
    
    dds <- DESeqDataSetFromMatrix(countData = exp2,
                                  colData = meta2,
                                  design= ~age)
    
    keep<-matrix(ncol=1,nrow=nrow(exp2));rownames(keep)=rownames(exp2)
    for (i in 1:nrow(exp2)){
      table=exp2[i,]
      if(table[order(table)][(ncol(exp2)/4)*3]>=1){
        keep[i,]=TRUE
      }else{keep[i,]=FALSE}
      
    }
    #row.names(keep)=row.names(exp2)
    keep2=as.logical(keep)
    
    #keep <- rowMedians(counts(dds)) >= 10
    dds <- dds[keep2,]
    dds <- DESeq(dds)
    resultsNames(dds)
    
    res<-results(dds,name="age")
    #res <- results(dds, contrast=c("sex","Male","Female"))
    #res$log2FoldChange=res$log2FoldChange*-1
    res=as.data.frame(res)
    tT=res
    tT$padj=p.adjust(tT$pvalue,method="BH")
    tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(1), "Decreased",
                        ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(1), "Increased", "Not Sig")) #for continuous variable, logfc needs to be 0
    tT=na.omit(tT)
    #write.csv(tT,paste0("deg/age_",celltypes[c],".csv"))
    deg[[c]]<<-tT
    names(deg)[c]<<-celltypes[c]
    
    
    tT2=tT[which((tT$log2FoldChange>log2(1)|tT$log2FoldChange<(-log2(1)))&tT$padj<0.05),]
    tT2=tT2[order(tT2$padj),]
    tT2$hgnc_symbol=row.names(tT2)
    
    top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:20,]
    top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:20,]
    top20<<-rbind(top20up,top20down)
    
    plotting_age<-function(x,y,label="comparisons",ct="celltype"){
      temp_plot= ggplot(x, aes(x = log2FoldChange, y = -log10(pvalue))) + # This line creates the axes plus background of the graph
        geom_point(aes(color = Legend)) + # This line plus previous coding line creates the volcano plot on to the graph with values coloured in three groups
        geom_hline(yintercept =-log10(max(y$pvalue)),colour="#990000", linetype="dashed")+
        geom_vline(xintercept =-log2(1),colour="#990000", linetype="dashed")+
        geom_vline(xintercept =log2(1),colour="#990000", linetype="dashed")+
        # Creates a line at the intercept which shows significance
        theme_bw(base_size = 12) + theme(legend.position = "bottom",plot.title = element_text(hjust=0.5,face="bold")) + # Changes the legend to below the graph
        labs(title=paste0(ct," ",label))+
        geom_text_repel(data = top20,
                        aes(label =hgnc_symbol),size = 3, box.padding = unit(0.35, "lines"),
                        point.padding = unit(0.3, "lines"),
                        max.overlaps = 30)     
      
      if(length(unique(x$Legend))>2){
        p=temp_plot+scale_color_manual(values = c("blue", "red","grey"))
      }else if(length(unique(x$Legend))==2&
               c(unique(x$Legend=="Decreased")|
                 c(unique(x$Legend=="Not Sig")))){
        p=temp_plot+scale_color_manual(values = c("blue", "grey"))
      }else if(length(unique(x$Legend))==2&
               c(unique(x$Legend=="Increased")|
                 c(unique(x$Legend=="Not Sig")))){
        p=temp_plot+scale_color_manual(values = c("red", "grey"))
      }else if(length(unique(x$Legend))==1){
        p=temp_plot+scale_colour_manual(values=c("grey"))
      }
      
      return(p)
    }
    #p=plotting_age(tT,tT2,label="age",ct=celltypes[c])
    
    
    #ggsave(paste0("plots/age_",celltypes[c],".png"),plot=p,width=2000,height = 1500,units="px",dpi=300)
  }
  
  tryCatch( diffexp(exp,meta), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
}

openxlsx::write.xlsx(deg,"deg/new2024/covariates/age_DEGs_eachcelltypes.xlsx",rowNames=T)



###Strains#####
deg=list()
for(c in 1:(length(celltypes)-1)){ #viral monocyte were only from one strains - messed up DE analyses
  
  skip_to_next<-FALSE
  
  diffexp<-function(x,y){
    y=y[c(y$strain=="C57BL/6J"|
            y$strain=="C57BL/6NCrl or C57BL/6N"),] #exclude ACE2 strain mice - strain only bred in one lab
    
    
    x=x[,y$combined_category]
    exp2=x[,startsWith(colnames(x),celltypes[c])]
    meta2=y[startsWith(y$combined_category,celltypes[c]),]
    
    meta2$strain=gsub("C57BL/6NCrl or C57BL/6N","C57BL/6N",meta2$strain)
    
    dds <- DESeqDataSetFromMatrix(countData = exp2,
                                  colData = meta2,
                                  design= ~strain)
    print(celltypes[c])
    keep<-matrix(ncol=1,nrow=nrow(exp2));rownames(keep)=rownames(exp2)
    for (i in 1:nrow(exp2)){
      table=exp2[i,]
      if(table[order(table)][(ncol(exp2)/4)*3]>=1){
        keep[i,]=TRUE
      }else{keep[i,]=FALSE}
      
    }
    #row.names(keep)=row.names(exp2)
    keep2=as.logical(keep)
    
    #keep <- rowMedians(counts(dds)) >= 10
    dds <- dds[keep2,]
    dds <- DESeq(dds)
    resultsNames(dds)
    
    res <- results(dds, contrast=c("strain","C57BL/6J","C57BL/6N"))
    tT=as.data.frame(res)
    tT$padj=p.adjust(tT$pvalue,method="BH")
    tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(2), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(2), "Increased", "Not Sig"))
    #write.csv(tT,paste0("deg/strain_6Jv6N",celltypes[c],".csv"))
    deg[[c]]<<-tT
    names(deg)[c]<<-celltypes[c]
    
    
    tT2=tT[which((tT$log2FoldChange>log2(2)|tT$log2FoldChange<(-log2(2)))&tT$padj<0.05),]
    tT2=tT2[order(tT2$padj),]
    tT2$hgnc_symbol=row.names(tT2)
    
    top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:20,]
    top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:20,]
    top20<<-rbind(top20up,top20down)
    
    #p=plotting(tT,tT2,label='6Jv6N',ct=celltypes[c])
    
    
    #ggsave(paste0("plots/strain_",celltypes[c],"_",names(res)[r],".png"),plot=p,width=2000,height = 1500,units="px",dpi=300);print(paste0(c,"_",names(res)[r]," saved"))
    
  }
  
  tryCatch( diffexp(exp,meta), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
  #openxlsx::write.xlsx(deg,paste0("deg/new2024/strain_",names(res)[r],".xlsx"),rowNames=T)
  
  
}

openxlsx::write.xlsx(deg,paste0("deg/new2024/covariates/strain_6Nvs6J_DEGs_eachcelltypes.xlsx"),rowNames=T)



###platforms#####
#DEG was done manually for each comparison. 
deg=list()

for(c in 1:length(celltypes)){
  
  skip_to_next<-FALSE
  
  diffexp<-function(x,y){
    y=y[c(y$platform=="Rhapsody"| y$platform=="Drop-seq"),]
    y$platform=factor(y$platform,levels=c("Rhapsody","Drop-seq"))
    
    exp2=x[,startsWith(colnames(x),celltypes[c])]
    meta2=y[startsWith(y$combined_category,celltypes[c]),]
    
    unique_values <- table(meta2$platform)
    factors_to_remove <- names(unique_values[unique_values == 1])
    # Remove factors with only one unique value from 'platform' column
    if (length(factors_to_remove) > 0) {
      for (factor in factors_to_remove) {
        meta2 <- subset(meta2, !(platform == factor))
      } 
      print("Factors with only one unique value removed.")
    } else { print("No factors with only one unique value found.")}   
    
    exp2=exp2[,meta2$combined_category]
    
    dds <- DESeqDataSetFromMatrix(countData = exp2,
                                  colData = meta2,
                                  design= ~platform)
    
    keep<-matrix(ncol=1,nrow=nrow(exp2));rownames(keep)=rownames(exp2)
    for (i in 1:nrow(exp2)){
      table=exp2[i,]
      if(table[order(table)][(ncol(exp2)/4)*3]>=1){
        keep[i,]=TRUE
      }else{keep[i,]=FALSE}
      
    }
    #row.names(keep)=row.names(exp2)
    keep2=as.logical(keep)
    
    #keep <- rowMedians(counts(dds)) >= 10
    dds <- dds[keep2,]
    dds <- DESeq(dds)
    print(resultsNames(dds))
    
    res<-results(dds)
    
    # res1 <- results(dds, contrast=c("platform","10X Genomics","Rhapsody"))
    # res2 <- results(dds, contrast=c("platform","10X Genomics","Drop-seq"))
    # res3 <- results(dds, contrast=c("platform","Rhapsody","Drop-seq"))
    # 
    # res<-list(res1,res2,res3);names(res)=c("10XvRhapsody","10xvDropseq","RhapsodyvDropseq")
    # 
    # tTdf=rownames(as.data.frame(res[[1]]))
    
    res=as.data.frame(res)
    tT=res
    tT$padj=p.adjust(tT$pvalue,method="BH")
    tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(1.5), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(1.5), "Increased", "Not Sig"))
    tT=na.omit(tT)
    #write.csv(tT,paste0("deg/sex_",celltypes[c],".csv"))
    deg[[c]]<<-tT
    names(deg)[c]<<-celltypes[c]
    
    
    tT2=tT[which((tT$log2FoldChange>log2(1.5)|tT$log2FoldChange<(-log2(1.5)))&tT$padj<0.05),]
    tT2=tT2[order(tT2$padj),]
    tT2$hgnc_symbol=row.names(tT2)
    
    top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:20,]
    top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:20,]
    top20<<-rbind(top20up,top20down)
    
    #p=plotting(tT,tT2,label=names(platforms)[r],ct=celltypes[c])
    
    
    #ggsave(paste0("plots/platform_",celltypes[c],"_",names(platforms)[r],".png"),plot=p,width=2000,height = 1500,units="px",dpi=300);print(paste0(celltypes[c],"_",names(res)[r]," saved"))
    
  }
  tryCatch( diffexp(exp,meta), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  # for(p in 1:3){
  #   platforms[[p]]<-deg[paste0(celltypes[c],"_",names(res[p]))]
  #   openxlsx::write.xlsx(platforms[[p]],paste0("deg/new2024/platform_",names(res[p]),"_eachcelltype.xlsx"),rowNames=T)
  #   
  #}
}

# cells=strsplit(names(deg),"_");cells=sapply(cells,`[`,1)
# names(deg)=cells
# deg=deg[-19]

openxlsx::write.xlsx(deg,paste0("deg/new2024/covariates/Supplementary Tables/platform_RhapsodyvsDropseq_DEGs_eachcelltype.xlsx"),rowNames=T)

###Smoking Status########
meta2=meta[which(meta$study=="Mouse_Copd"|
                   meta$study=="Mouse_Copd_Covid"|
                   meta$study=="Mouse_COPD_SHAM_SARSCov2"),]
exp2=exp[,colnames(exp) %in% meta2$combined_category];dim(exp2)
meta2=meta2[meta2$combined_category%in%colnames(exp2),];dim(meta2)

smokecelltypes <- celltypes[grepl("^Smoke", celltypes)]
smokecelltypes=gsub("| |","",smokecelltypes)

#exp2=exp2+1

deg=list()
for(c in 1:length(smokecelltypes)){
  print(smokecelltypes[c])
  
  skip_to_next<-FALSE
  
  diffexp<-function(x,y){
    
    exp3=x[,startsWith(colnames(x),smokecelltypes[c])]
    meta3=y[startsWith(y$combined_category,smokecelltypes[c]),]
    #meta3$smoking_status=gsub("Non-Smoking","Control",meta3$smoking_status)
    #meta3$smoking_status[meta3$Disease == "Copd_Covid"] <- "Smoke"
    #exp3 <- exp3+1
    exp3=as.matrix(exp3)
    class(exp3)<-'integer'
    dds <- DESeqDataSetFromMatrix(countData = exp3,
                                  colData = meta3,
                                  design= ~0+smoking_status)
    design_matrix <- model.matrix(design(dds), colData(dds))
    
    keep<-matrix(ncol=1,nrow=nrow(exp3));rownames(keep)=rownames(exp3)
    for (i in 1:nrow(exp3)){
      table=exp3[i,]
      if(table[order(table)][round((ncol(exp3)/4)*3)]>=1){
        keep[i,]=TRUE
      }else{keep[i,]=FALSE}
      
    }
    #row.names(keep)=row.names(exp2)
    keep2=as.logical(keep)
    
    #keep <- rowMedians(counts(dds)) >= 10
    dds <- dds[keep2,]
    dds <- DESeq(dds)
    print("deg analysis results out")
    print(resultsNames(dds))
    
    res <- results(dds, contrast=c("smoking_status","Smoke","Non-Smoke"))
    #res$log2FoldChange=res$log2FoldChange*-1
    res=as.data.frame(res)
    tT=res
    tT$padj=p.adjust(tT$pvalue,method="BH")
    tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(1.5), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(1.5), "Increased", "Not Sig"))
    tT=na.omit(tT)
    write.csv(tT,paste0("data/deg/newFia/Smoke_",smokecelltypes[c],"_130125.csv"))
    deg[[c]]<<-tT
    names(deg)[c]<<-smokecelltypes[c]
    
    
    tT2=tT[which((tT$log2FoldChange>log2(1.5)|tT$log2FoldChange<(-log2(1.5)))&tT$padj<0.05),]
    tT2=tT2[order(tT2$padj),]
    tT2$hgnc_symbol=row.names(tT2)
    
    top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:20,]
    top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:20,]
    top20<<-rbind(top20up,top20down)
    
    p=plotting(tT,tT2, label="Smoking v Control",ct=smokecelltypes[c])
    
    ggsave(paste0("plots/newFia/smoking_",smokecelltypes[c],"_130125.png"),plot=p,width=2000,height = 1500,units="px",dpi=300);print(paste0(smokecelltypes[c]," saved"))
  }
  
  tryCatch( diffexp(exp2,meta2), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
}


openxlsx::write.xlsx(deg,"data/deg/newFia/smoking_DEGs_allsmokecelltypes_130125.xlsx",rowNames=T)

####Validating Smoking DEGs##########

deg=read.csv("deg/new2024/covariates/Supplementary Tables/S9_smoking_DEGs_allsmokingcelltypes.csv",row.names=1)
deg=deg[[1]]
deg=deg[which(deg$Legend=="Increased"|deg$Legend=="Decreased"),]
deg=deg[order(deg$padj),]

up=deg[deg$Legend=="Increased",]
down=deg[deg$Legend=="Decreased",]

#deg_eachct=as.data.frame(readxl::read_xlsx("deg/new2024/covariates/Supplementary Tables/S10_smoking_DEGs_eachsmokingcelltypes.xlsx"))
#rownames(deg_eachct)=deg_eachct[,1];deg_eachct=deg_eachct[,-1]

counts=read.csv("data/deg/newFia/GSE199853_genecounts.csv",row.names=1)
s=read.delim("data/deg/newFia/GSE199853_series_matrix_main.txt")

colnames(counts)=s$X.Sample_title

library(org.Mm.eg.db)
# Remove version numbers from ENSEMBL IDs
rownames(counts) <- gsub("\\..*", "", rownames(counts))

# Map ENSEMBL IDs to SYMBOLs
mart <- mapIds(
  org.Mm.eg.db,
  keys = rownames(counts),
  column = 'SYMBOL',
  keytype = 'ENSEMBL'
)


genenames=select(
  org.Mm.eg.db,
  keys = rownames(counts),
  column = c('SYMBOL', 'ENTREZID', 'ENSEMBL'),
  keytype = 'ENSEMBL')

# Remove duplicates in genenames based on SYMBOL, keeping the first occurrence
genenames_unique <- genenames[!duplicated(genenames$SYMBOL), ]
genenames_unique<-na.omit(genenames_unique)
# Filter counts to include only rows with ENSEMBL IDs present in genenames$ENSEMBL
counts_unique <- counts[rownames(counts) %in% genenames_unique$ENSEMBL, ]

# Match row names of counts_unique to ENSEMBL in genenames_unique
matching_indices <- match(rownames(counts_unique), genenames_unique$ENSEMBL)
genenames_unique

# Replace row names with corresponding SYMBOLs from genenames_unique
rownames(counts_unique) <- genenames_unique$SYMBOL[matching_indices]

# View the updated counts matrix
counts_unique





up_int=intersect(rownames(counts_unique),rownames(up))

down_int=intersect(rownames(counts_unique),rownames(down))


#create data frame containing all the necessary information 
gsva_res=gsva(as.matrix(counts_unique),list(up_int,down_int),mx.diff=T,verbose=F,parallel.sz=1)
rownames(gsva_res)=c("Upregulated Signature", "Downregulated Signature")


df=data.frame(groups=colnames(counts),
              t(gsva_res))
df$groups=ifelse(grepl("DEX",df$groups),"Remove",df$groups) 
df$groups=ifelse(grepl("HDM",df$groups),"Remove",df$groups) 


df=df[df$groups!='Remove',]


df$groups <- ifelse(grepl("CS", df$groups) & grepl("11", df$groups)&grepl('Parenchyma',df$groups), 
                    "Parenchyma 11 weeks Smoke", 
                    df$groups) #replace all components of smoke column that contains S with Smoke, otherwise Control

df$groups <- ifelse(grepl("CS", df$groups) & grepl("11", df$groups)&grepl('Airway',df$groups), 
                    "Airway 11 weeks Smoke", 
                    df$groups)

df$groups <- ifelse(grepl("PBS", df$groups) & grepl("11", df$groups)& grepl('Airway',df$groups), 
                    "Airway 11 weeks Air", 
                    df$groups) 

df$groups <- ifelse(grepl("PBS", df$groups) & grepl("3", df$groups)& grepl('Airway',df$groups), 
                    "Airway 3 weeks Air", 
                    df$groups) 

df$groups <- ifelse(grepl("PBS", df$groups) & grepl("3", df$groups)& grepl('Parenchyma',df$groups), 
                    "Parenchyma 3 weeks Air", 
                    df$groups) 

df$groups <- ifelse(grepl("PBS", df$groups) & grepl("11", df$groups)& grepl('Parenchyma',df$groups), 
                    "Parenchyma 11 weeks Air", 
                    df$groups) 
df$region=df$groups
df <- df %>%
  separate(groups, into = c("region", "groups"), sep = " ", extra = "merge", fill = "right")

df$treatment =df$groups
df <- df %>%
  separate(groups, into = c("groups", "treatment"), sep = "weeks ", extra = "merge", fill = "right")

df$groups=paste0(df$groups,"weeks")
colnames(df)[colnames(df) == "groups"] <- "time"

full_df=df

df=full_df[full_df$region=='Parenchyma'&full_df$time=='11 weeks',]

for(d in colnames(df)[(ncol(df)-2):ncol(df)]){
  
  p=ggplot(df, aes(x=treatment, y=as.numeric(df[,d]) , fill=treatment))+
    geom_boxplot(alpha=0.8, outlier.shape = NA)+
    geom_point(aes(y=as.numeric(df[,d]),fill=treatment),  size=3, shape=21
               ,position=position_jitterdodge()) +
    theme_classic()+
    labs(y = paste0("GSVA Enrichment Score"),title=d)+
    theme(axis.text.x = element_text(color = "black", size = 16, face = "bold"
                                     , hjust = 1, vjust = 1, angle = 45),
          axis.text.y = element_text(color = "black", size = 12, angle = 0, 
                                     hjust = 0.5, vjust = 0.5, face = "bold"),  
          axis.title.x = element_blank(),
          axis.title.y = element_text(color = "black", size = 14, angle = 90, 
                                      hjust = .5, vjust = 2, face = "bold"),
          legend.position = "none",
          plot.title = element_text(face="bold",hjust = 0.5,size=16),
          text=element_text(face = "bold"))+ 
    scale_fill_manual(values=c("#00008B", "#FF9999", "#FF6666", "#FF3333", "#CC0000"))+
    geom_signif(
      comparisons = list(c("Air", "Smoke")),
      test="t.test",fontface="bold",textsize=3.75,
      map_signif_level = function(p) sprintf("%.2g", p),
      step_increase = 0.08) 
  
  print(p)
  
  ggsave(paste0("plots/newFia/publiclyavailable_smoke_",d,"_130125.png"),width=3.5,height=5,units="in")
}

###Viral#######
meta2=meta[which(meta$study=="Mouse_Covid"|
                   meta$study=="Mouse_Copd_Covid"|
                   meta$study=="Mouse_COPD_SHAM_SARSCov2"|
                   meta$study=="Mouse_Herpesvirus"|
                   meta$study=="Mouse_influenza"),]
exp2=exp[,colnames(exp) %in% meta2$combined_category]

deg=list()
for(c in 1:length(celltypes)){
  
  skip_to_next<-FALSE
  
  diffexp<-function(x,y){
    
    exp3=x[,startsWith(colnames(x),celltypes[c])]
    meta3=y[startsWith(y$combined_category,celltypes[c]),]
    
    dds <- DESeqDataSetFromMatrix(countData = exp3,
                                  colData = meta3,
                                  design= ~Viral+age+sex)
    design_matrix <- model.matrix(design(dds), colData(dds))
    
    keep<-matrix(ncol=1,nrow=nrow(exp3));rownames(keep)=rownames(exp3)
    for (i in 1:nrow(exp3)){
      table=exp3[i,]
      if(table[order(table)][round((ncol(exp3)/4)*3)]>=1){
        keep[i,]=TRUE
      }else{keep[i,]=FALSE}
      
    }
    #row.names(keep)=row.names(exp2)
    keep2=as.logical(keep)
    
    #keep <- rowMedians(counts(dds)) >= 10
    dds <- dds[keep2,]
    dds <- DESeq(dds)
    print("deg analysis results out")
    print(resultsNames(dds))
    
    res <- results(dds, contrast=c("Viral","Viral Group","Non-viral"))
    #res$log2FoldChange=res$log2FoldChange*-1
    res=as.data.frame(res)
    tT=res
    tT$padj=p.adjust(tT$pvalue,method="BH")
    tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(2), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(2), "Increased", "Not Sig"))
    tT=na.omit(tT)
    #write.csv(tT,paste0("deg/sex_",celltypes[c],".csv"))
    deg[[c]]<<-tT
    names(deg)[c]<<-celltypes[c]
    
    
    tT2=tT[which((tT$log2FoldChange>log2(2)|tT$log2FoldChange<(-log2(2)))&tT$padj<0.05),]
    tT2=tT2[order(tT2$padj),]
    tT2$hgnc_symbol=row.names(tT2)
    
    top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:20,]
    top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:20,]
    top20<<-rbind(top20up,top20down)
    
    p=plotting(tT,tT2, label="Viral v Control",ct=celltypes[c])
    
    ggsave(paste0("plots/viral_",celltypes[c],"_0207.png"),plot=p,width=2000,height = 1500,units="px",dpi=300);print(paste0(celltypes[c]," saved"))
  }
  
  tryCatch( diffexp(exp2,meta2), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
}

#Disease across all cell types########
exp=read.delim("deg/pseudobulk_percelltypes.tsv", row.names=1)
colnames(exp)=gsub("^(X)","",colnames(exp)) #remove string that's only at the beginning of the character
colnames(exp)=gsub("[.]","",colnames(exp))
celltypes=strsplit(colnames(exp),"_");celltypes=sapply(celltypes,`[`,1)
celltypes=unique(celltypes)
meta=read.csv("deg/meta_pseudobulk_percelltypes.csv")
#a=meta$sex;a[grep("Male",a)]="Male"
#meta$sex=a
meta$combined_category=gsub("[.]|/|\\^| |[+]","",meta$combined_category)
#meta=meta[meta$study!="Mouse_Post_Sendai",]
meta$age=gsub("18-20","19",meta$age);meta$age=gsub("16-20","18",meta$age);meta$age=gsub("12-13","13",meta$age);meta$age=gsub("10-12","11",meta$age)
meta$age=gsub("16-36","26",meta$age);meta$age=gsub("Postnatal Day 7","1",meta$age);meta$age=gsub("6-20","13",meta$age);meta$age=gsub("6-8","7",meta$age);
meta$age=gsub("5-6","6",meta$age);meta$age=gsub("\\b10\\b","14",meta$age) #use \\b to limit to only that string
meta$age=as.numeric(meta$age)

#exp=exp[,meta$combined_category]
exp=as.matrix(exp)
class(exp)<-'integer'

columns=c("Age","Asthma","Asthma.w.Steroids","Bleoold","Bleoyoung","Cancer","CancerTumor",
          "Cigarette.Smoke","Copd","CopdCovid","Covid","Fibrosis","Herpes", "Hyperoxia",
          "Hypoxia","Influenza","Mir155KO","Oldmice","Postsendai","Radiation",
          "Severe.Asthma","Severe.Asthma.w.Steroids","TB")
deg=list()
for(d in 1:length(columns)){
  skip_to_next<-FALSE
  
  diffexp<-function(x,y){
    for(c in 1:length(celltypes)){
      
      exp2=x[,startsWith(colnames(x),celltypes[c])]
      meta2=y[startsWith(y$combined_category,celltypes[c]),]
      
      meta2=meta2[meta2[,columns[d]]!="",]
      meta2$DE_group=meta2[,columns[d]];print(table(meta2$DE_group))
      exp2=exp2[,meta2$combined_category]
      
      dds <- DESeqDataSetFromMatrix(countData = exp2,
                                    colData = meta2,
                                    design= ~DE_group)
      
      keep<-matrix(ncol=1,nrow=nrow(exp2));rownames(keep)=rownames(exp2)
      for (i in 1:nrow(exp2)){
        table=exp2[i,]
        if(table[order(table)][(ncol(exp2)/4)*3]>=1){
          keep[i,]=TRUE
        }else{keep[i,]=FALSE}
        
      }
      #row.names(keep)=row.names(exp2)
      keep2=as.logical(keep)
      
      #keep <- rowMedians(counts(dds)) >= 10
      dds <- dds[keep2,]
      dds <- DESeq(dds)
      resultsNames(dds)
      
      res <- results(dds, contrast=c("DE_group","Disease","Control"))
      res=as.data.frame(res)
      tT=res
      tT$padj=p.adjust(tT$pvalue,method="BH")
      tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(2), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(2), "Increased", "Not Sig"))
      tT=na.omit(tT)
      tT=tT[order(tT$Legend),]
      write.csv(tT,paste0("deg/",columns[d],"_",celltypes[c],".csv"))
      deg[[c]]<<-tT
      names(deg)[c]<<-celltypes[c]
      
      
      tT2=tT[which((tT$log2FoldChange>log2(2)|tT$log2FoldChange<(-log2(2)))&tT$padj<0.05),]
      tT2=tT2[order(tT2$padj),]
      tT2$hgnc_symbol=row.names(tT2)
      
      top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:20,]
      top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:20,]
      top20<<-rbind(top20up,top20down)
      
      p=plotting(tT,tT2,label=columns[d],ct=celltypes[c])
      
      ggsave(paste0("plots/",columns[d],"_",celltypes[c],".png"),plot=p,width=2000,height = 1500,units="px",dpi=300);print(paste0(celltypes[c],columns[d]," saved"))
      openxlsx::write.xlsx(deg,paste0("deg/",columns[d],"_eachcelltypes.xlsx"),rowNames=T)
      
    }
  }
  tryCatch( diffexp(exp,meta), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
}





##Heatmap####
#average the females and males in each cell type -> heatmap
exp2=exp+1
gset=vst(as.matrix(exp2))

males=gset[,meta[meta$sex=="Male","combined_category"]]
females=gset[,meta[meta$sex=="Female","combined_category"]]

m2=data.frame(gene=rownames(gset))
f2=data.frame(gene=rownames(gset))
df=data.frame(gene=rownames(gset)) #number of celltypesx2 for male and female
for(c in 1:length(celltypes)){
  m=as.data.frame(rowMeans(males[,startsWith(colnames(males),celltypes[c])]));colnames(m)=paste0(celltypes[c],"_Male")
  f=as.data.frame(rowMeans(females[,startsWith(colnames(females),celltypes[c])]));colnames(f)=paste0(celltypes[c],"_Female")
  m2=cbind(m2,m)
  f2=cbind(f2,f)
}
df=cbind(m2[,-1],f2[,-1])

pheatmap(df[c("Xist","Ddx3y","Eif2s3y"),],cluster_cols=F,cluster_rows = F,scale = "row")


#Check individual genes####
exp2=exp+1
gset=vst(as.matrix(exp2))

#AM= "Rabgap1l", "Sc5d"
#AT2= "Ighm","Fhad1"
#B-cells = "Zfp383","Nudt6", "Gm5617","Rbbp9","Il12a","3010003L21Rik","Ddx25","Gon7","Prr7",
#"2810001G20Rik","Tmem108","Spns2", "Ccdc88a","Pon3","Hepacam2","Sdccag8","Itgb1","Kif1b","Zfp955b",      
#"Hmgn3","Tent5c","Kcnmb4","Ptms","Nacc2","Atp1a2","Racgap1","Igha","Ptpn13"
#Capillary EC= "Serpina3m","Eif2s3x","Serpina3n","Cplane1","Sdhaf1","Syne1","Bnip3","Clec4g",
#"Igha","Tmcc1","Fam220a","Gon4l","Cemip2","Ube2v1"
#
#Granulocyte

d1=data.frame(gene=gset["Sfta3-ps",],
              groups=meta2$sex)
p<- ggplot(d1, aes(x=groups, y=gene))+
  geom_boxplot(alpha=0.8, outlier.shape = NA,aes(fill=groups))+
  geom_point(aes(y=gene,fill=groups),  size=1, shape=21
             ,position=position_jitterdodge())+theme_classic()




#DEG across all cell types#####
exp_persample=as.matrix(read.delim("data/deg/newFia/pseudobulk_permice_2910_HLCA_method_raw.tsv",row.names=1))
colnames(exp_persample)=gsub("^(X)","",colnames(exp_persample))
colnames(exp_persample)=gsub("[.]","",colnames(exp_persample))
class(exp_persample)<-'integer'

meta_persample=read.csv("data/deg/newFia/Meta_Pseudo_0611.csv")
meta_persample$Pseudo=gsub("[.]|/|\\^| |[+]","",meta_persample$Pseudo)
meta_persample=meta_persample[meta_persample$study!="Mouse_Post_Sendai",]
length(intersect(colnames(exp_persample),meta_persample$Pseudo)) #expect 184 - all samples are the same
meta_persample$age=gsub("18-20","19",meta_persample$age);meta_persample$age=gsub("16-20","18",meta_persample$age);meta_persample$age=gsub("12-13","13",meta_persample$age);meta_persample$age=gsub("10-12","11",meta_persample$age)
meta_persample$age=gsub("16-36","26",meta_persample$age);meta_persample$age=gsub("Postnatal Day 7","1",meta_persample$age);meta_persample$age=gsub("6-20","13",meta_persample$age);meta_persample$age=gsub("6-8","7",meta_persample$age);
meta_persample$age=gsub("Unknon","6",meta_persample$age);meta_persample$age=gsub("\\b10\\b","14",meta_persample$age) #use \\b to limit to only that string
meta_persample$age=as.numeric(meta_persample$age)
meta_persample$strain[meta_persample$study == 'Mouse_Covid'] <- 'K18-ACE2(2Prlmn/J)'
meta_persample$strain[meta_persample$study == 'Mouse_Fibro_Age'] <- 'C57BL6/J-rj'
meta_persample$platform[meta_persample$study == 'Mouse_Covid'] <- 'Rhapsody'
meta_persample$platform[meta_persample$study == 'Mouse_Fibro_Age'] <- 'Drop-seq'
controls=meta_persample[meta_persample$Disease=='Control',]

#meta2=meta_persample

#exclude ACE2 strain mice - strain only bred in one lab
#meta2$strain=gsub("C57BL/6NCrl or C57BL/6N","C57BL/6N",meta2$strain)
#meta2$strain=gsub("C57BL/6NTac","C57BL/6N",meta2$strain)

#exp2=exp_persample

#exp2=exp2[,meta2$Pseudo]
#exp2=as.matrix(exp2)


###strains#####
meta2=controls
meta2=meta2[c(meta2$strain!="C57BL/6J-Notch3-def"&
                meta2$strain!="C57BL/6J-LAIR1KO"),] 
exp2=exp_persample[,meta2$Pseudo]

#gene="Atg2b"
#d1=data.frame(expr=exp2[gene,],
#               group=meta2$strain)
#ggplot(d1,aes(y=expr,x=group))+geom_boxplot()+geom_point(aes(y=expr,fill=group))

dds <- DESeqDataSetFromMatrix(countData = exp2,
                              colData = meta2,
                              design= ~age+sex+strain)

keep<-matrix(ncol=1,nrow=nrow(exp2));rownames(keep)=rownames(exp2)
for (i in 1:nrow(exp2)){
  table=exp2[i,]
  if(table[order(table)][(ncol(exp2)/4)*3]>=10){
    keep[i,]=TRUE
  }else{keep[i,]=FALSE}
  
}
#row.names(keep)=row.names(exp2)
keep2=as.logical(keep)

#keep <- rowMedians(counts(dds)) >= 10
dds <- dds[keep2,]
dds <- DESeq(dds)
resultsNames(dds)

res1<- results(dds, contrast=c("strain","C57BL/6J","C57BL/6NCrl or C57BL/6N"))
res2<- results(dds, contrast=c("strain","C57BL6/J-rj","C57BL/6J"))
res3<- results(dds, contrast=c("strain","K18-ACE2(2Prlmn/J)","C57BL/6J"))
res4<- results(dds, contrast=c("strain","C57BL/6NTac","C57BL/6J"))
res5<- results(dds, contrast=c("strain","K18-ACE2(2Prlmn/J)","C57BL/6NTac"))
res6<- results(dds, contrast=c("strain","C57BL/6NTac","C57BL/6NCrl or C57BL/6N"))
res7<- results(dds, contrast=c("strain","C57BL/6NTac","C57BL6/J-rj"))
res8<- results(dds, contrast=c("strain","K18-ACE2(2Prlmn/J)","C57BL/6NCrl or C57BL/6N"))
res9<- results(dds, contrast=c("strain","K18-ACE2(2Prlmn/J)","C57BL6/J-rj"))


res<-list(res1,res2,res3,res4,res5,res6,res7,res8,res9);names(res)=c("6J vs 6N","6J-rj vs 6J","ACE2KO vs 6J","6NTac vs 6J","ACE2KO vs 6NTac","6NTac vs 6N",
                                                                     "6NTac vs 6J-rj","ACE2KO vs 6N", "ACE2KO vs 6J-rj")

for (r in 1:length(res)){
  #res$log2FoldChange=res$log2FoldChange*-1
  tT=as.data.frame(res[[r]])
  tT=na.omit(tT)
  tT$padj=p.adjust(tT$pvalue,method="BH")
  tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(2), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(2), "Increased", "Not Sig"))
  write.csv(tT,paste0("data/deg/newFia/strain_allcells_",names(res)[r],"_0912.csv"))
  
  tT2=tT[which((tT$log2FoldChange>log2(2)|tT$log2FoldChange<(-log2(2)))&tT$padj<0.05),]
  tT2=tT2[order(tT2$padj),]
  tT2$hgnc_symbol=row.names(tT2)
  
  top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:40,]
  top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:40,]
  top20<<-rbind(top20up,top20down)
  
  p=plotting(tT,tT2,label=names(res)[r],ct="All Cells")
  
  
  ggsave(paste0("plots/newFia/strain_allcells_",names(res)[r],"_0912.png"),plot=p,width=2000,height = 1500,units="px",dpi=300);print(paste0("allcells_",names(res)[r]," saved"))
}


tT=as.data.frame(res)
tT$padj=p.adjust(tT$pvalue,method="BH")
tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(2), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(2), "Increased", "Not Sig"))
write.csv(tT,paste0("data/deg/newFia/strain_allcells_6Jv6N_2102.csv"))

tT2=tT[which((tT$log2FoldChange>log2(2)|tT$log2FoldChange<(-log2(2)))&tT$padj<0.05),]
tT2=tT2[order(tT2$padj),]
tT2$hgnc_symbol=row.names(tT2)

top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:40,]
top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:40,]
top20<<-rbind(top20up,top20down)

p=plotting(tT,tT2,label="6Jv6N",ct="All Cells")


ggsave(paste0("plots/newFia/strain_allcells_6Jv6N_1911.png"),plot=p,width=2000,height = 1500,units="px",dpi=300)



###age####
meta2=controls
exp2=exp_persample
exp2=exp2[,meta2$Pseudo]
dds <- DESeqDataSetFromMatrix(countData = exp2,
                              colData = meta2,
                              design= ~age)

keep<-matrix(ncol=1,nrow=nrow(exp2));rownames(keep)=rownames(exp2)
for (i in 1:nrow(exp2)){
  table=exp2[i,]
  if(table[order(table)][(ncol(exp2)/4)*3]>=10){
    keep[i,]=TRUE
  }else{keep[i,]=FALSE}
  
}
#row.names(keep)=row.names(exp2)
keep2=as.logical(keep)

#keep <- rowMedians(counts(dds)) >= 10
dds <- dds[keep2,]
dds <- DESeq(dds)
resultsNames(dds)

res<-results(dds,name="age")
#res <- results(dds, contrast=c("sex","Male","Female"))
#res$log2FoldChange=res$log2FoldChange*-1
res=as.data.frame(res)
tT=res
tT$padj=p.adjust(tT$pvalue,method="BH")
tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(1), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(1), "Increased", "Not Sig"))
tT=na.omit(tT)
write.csv(tT,paste0("data/deg/newFia/age_covar_allcells_1911.csv"))

tT2=tT[which((tT$log2FoldChange>log2(1)|tT$log2FoldChange<(-log2(1)))&tT$padj<0.05),]
tT2=tT2[order(tT2$padj),]
tT2$hgnc_symbol=row.names(tT2)

top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:40,]
top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:40,]
top20<<-rbind(top20up,top20down)

plotting_age<-function(x,y,label="comparisons",ct="celltype"){
  temp_plot= ggplot(x, aes(x = log2FoldChange, y = -log10(pvalue))) + # This line creates the axes plus background of the graph
    geom_point(aes(color = Legend)) + # This line plus previous coding line creates the volcano plot on to the graph with values coloured in three groups
    geom_hline(yintercept =-log10(max(y$pvalue)),colour="#990000", linetype="dashed")+
    geom_vline(xintercept =-log2(1),colour="#990000", linetype="dashed")+
    geom_vline(xintercept =log2(1),colour="#990000", linetype="dashed")+
    # Creates a line at the intercept which shows significance
    theme_bw(base_size = 12) + theme(legend.position = "bottom",plot.title = element_text(hjust=0.5,face="bold")) + # Changes the legend to below the graph
    labs(title=paste0(ct," ",label))+
    geom_text_repel(data = top20,
                    aes(label =hgnc_symbol),size = 3, box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"),
                    max.overlaps = 50)     
  
  if(length(unique(x$Legend))>2){
    p=temp_plot+scale_color_manual(values = c("blue", "red","grey"))
  }else if(length(unique(x$Legend))==2&
           c(unique(x$Legend=="Decreased")|
             c(unique(x$Legend=="Not Sig")))){
    p=temp_plot+scale_color_manual(values = c("blue", "grey"))
  }else if(length(unique(x$Legend))==2&
           c(unique(x$Legend=="Increased")|
             c(unique(x$Legend=="Not Sig")))){
    p=temp_plot+scale_color_manual(values = c("red", "grey"))
  }else if(length(unique(x$Legend))==1){
    p=temp_plot+scale_colour_manual(values=c("grey"))
  }
  
  return(p)
}
p=plotting_age(tT,tT2, label="Age",ct="All Cells")


ggsave(paste0("plots/newFia/age_covar_allcells_1911.png"),plot=p,width=2000,height = 1500,units="px",dpi=300)


###sex#####
meta2=controls
exp2=exp_persample
exp2=exp2[,meta2$Pseudo]
dds <- DESeqDataSetFromMatrix(countData = exp2,
                              colData = meta2,
                              design= ~sex)

keep<-matrix(ncol=1,nrow=nrow(exp2));rownames(keep)=rownames(exp2)
for (i in 1:nrow(exp2)){
  table=exp2[i,]
  if(table[order(table)][round((ncol(exp2)/4)*3)]>=10){
    keep[i,]=TRUE
  }else{keep[i,]=FALSE}
  
}
#row.names(keep)=row.names(exp2)
keep2=as.logical(keep)

#keep <- rowMedians(counts(dds)) >= 10
dds <- dds[keep2,]
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast=c("sex","Male","Female"))
#res$log2FoldChange=res$log2FoldChange*-1
res=as.data.frame(res)
tT=res
tT$padj=p.adjust(tT$pvalue,method="BH")
tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(2), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(2), "Increased", "Not Sig"))
tT=na.omit(tT)
write.csv(tT,paste0("data/deg/newFia/sex_allcells_1911.csv"))

tT2=tT[which((tT$log2FoldChange>log2(2)|tT$log2FoldChange<(-log2(2)))&tT$padj<0.05),]
tT2=tT2[order(tT2$padj),]
tT2$hgnc_symbol=row.names(tT2)

top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:20,]
top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:20,]
top20<<-rbind(top20up,top20down)

p=plotting(tT,tT2, label="Male v Female",ct="All Cells")

ggsave(paste0("plots/newFia/sex_allcells_1911.png"),plot=p,width=2000,height = 1500,units="px",dpi=300)


###platforms#####
meta2=controls
meta2=meta2[c(meta2$platform=="10X Genomics"|
                meta2$platform=="Drop-seq"|
                meta2$platform=="Rhapsody"),] #subset only to platforms that were used by multiple studies
exp2=exp_persample
exp2=exp2[,meta2$Pseudo]

dds <- DESeqDataSetFromMatrix(countData = exp2,
                              colData = meta2,
                              design= ~platform)

keep<-matrix(ncol=1,nrow=nrow(exp2));rownames(keep)=rownames(exp2)
for (i in 1:nrow(exp2)){
  table=exp2[i,]
  if(table[order(table)][(ncol(exp2)/4)*3]>=10){
    keep[i,]=TRUE
  }else{keep[i,]=FALSE}
  
}
#row.names(keep)=row.names(exp2)
keep2=as.logical(keep)

#keep <- rowMedians(counts(dds)) >= 10
dds <- dds[keep2,]
dds <- DESeq(dds)
resultsNames(dds)

res1 <- results(dds, contrast=c("platform","10X Genomics","Rhapsody"))
res2 <- results(dds, contrast=c("platform","10X Genomics","Drop-seq"))
res3 <- results(dds, contrast=c("platform","Rhapsody","Drop-seq"))

res<-list(res1,res2,res3);names(res)=c("10XvRhapsody","10xvDropseq","RhapsodyvDropseq")

for (r in 1:length(res)){
  #res$log2FoldChange=res$log2FoldChange*-1
  tT=as.data.frame(res[[r]])
  tT$padj=p.adjust(tT$pvalue,method="BH")
  tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(2), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(2), "Increased", "Not Sig"))
  write.csv(tT,paste0("data/deg/newFia/platform_allcells_",names(res)[r],"_1911.csv"))
  
  tT2=tT[which((tT$log2FoldChange>log2(2)|tT$log2FoldChange<(-log2(2)))&tT$padj<0.05),]
  tT2=tT2[order(tT2$padj),]
  tT2$hgnc_symbol=row.names(tT2)
  
  top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:50,]
  top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:20,]
  top20<<-rbind(top20up,top20down)
  
  p=plotting(tT,tT2,label=names(res)[r],ct="All Cells")
  
  
  ggsave(paste0("plots/newFia/platform_allcells_",names(res)[r],"_1911.png"),plot=p,width=2000,height = 1500,units="px",dpi=300);print(paste0("allcells_",names(res)[r]," saved"))
}


###Disease#####
colnames(meta_persample)[colnames(meta_persample) == "Severe.Asthma"] <- "Severe.Asthma.w.Chlamydia"
colnames(meta_persample)[colnames(meta_persample) == "Asthma"] <- "Severe.Asthma"
colnames(meta_persample)[colnames(meta_persample) == "Severe.Asthma.w.Steroids"] <- "Severe.Asthma.w.Chlamydia.w.Steroids"
colnames(meta_persample)[colnames(meta_persample) == "Asthma.w.Steroids"] <- "Severe.Asthma.w.Steroids"

columns=c("Age","Severe.Asthma","Severe.Asthma.w.Steroids","Bleoold","Bleoyoung","Cancer","CancerTumor",
          "Chlamydia","Cigarette.Smoke","Copd","CopdCovid","Covid","Fibrosis","Herpes", "Hyperoxia",
          "Hypoxia","Influenza","Oldmice","Radiation",
          "Severe.Asthma.w.Chlamydia","Severe.Asthma.w.Chlamydia.w.Steroids","TB")
#meta2=meta_persample
#exp2=exp_persample
#exp2=exp2[,meta2$Pseudo]

deg=list()
for(d in 1:length(columns)){
  skip_to_next<-FALSE
  
  diffexp<-function(x,y){
    
    meta2=y[y[,columns[d]]!="",]
    meta2$DE_group=meta2[,columns[d]]
    exp2=x
    exp2=exp2[,meta2$Pseudo]
    
    ###rename to Disease and Control ONLY!!!!##################
    meta2$DE_group=ifelse(startsWith(meta2$DE_group,"Control_"), "Control", meta2$DE_group)
    meta2$DE_group=ifelse(meta2$DE_group != "Control", "Disease", meta2$DE_group)
    
    
    dds <- DESeqDataSetFromMatrix(countData = exp2,
                                  colData = meta2,
                                  design= ~DE_group)
    
    keep<-matrix(ncol=1,nrow=nrow(exp2));rownames(keep)=rownames(exp2)
    for (i in 1:nrow(exp2)){
      table=exp2[i,]
      if(table[order(table)][(ncol(exp2)/4)*3]>=10){
        keep[i,]=TRUE
      }else{keep[i,]=FALSE}
      
    }
    #row.names(keep)=row.names(exp2)
    keep2=as.logical(keep)
    
    #keep <- rowMedians(counts(dds)) >= 10
    dds <- dds[keep2,]
    dds <- DESeq(dds)
    resultsNames(dds)
    
    res <- results(dds, contrast=c("DE_group","Disease","Control"))
    res=as.data.frame(res)
    tT=res
    tT$padj=p.adjust(tT$pvalue,method="BH")
    tT$Legend <- ifelse(tT$padj < 0.05&tT$log2FoldChange < -log2(2), "Decreased",ifelse(tT$padj < 0.05&tT$log2FoldChange > log2(2), "Increased", "Not Sig"))
    tT=na.omit(tT)
    write.csv(tT,paste0("deg/new2024/",columns[d],"_allcells_2102.csv"))
    deg[[d]]<<-tT
    names(deg)[d]<<-columns[d]
    print(table(tT$Legend));print(columns[d])
    
    
    
    tT2=tT[which((tT$log2FoldChange>log2(2)|tT$log2FoldChange<(-log2(2)))&tT$padj<0.05),]
    tT2=tT2[order(tT2$padj),]
    tT2$hgnc_symbol=row.names(tT2)
    
    top20up=tT2[tT2$Legend=="Increased",];top20up=top20up[1:40,]
    top20down=tT2[tT2$Legend=="Decreased",];top20down=top20down[1:40,]
    top20<<-rbind(top20up,top20down)
    
    p=plotting(tT,tT2,label=columns[d],ct="All Cells")
    
    ggsave(paste0("plots/new2024/",columns[d],"_allcells_2102.png"),plot=p,width=2000,height = 1500,units="px",dpi=300);print(paste0("allcells_",columns[d]," saved"))
  }
  tryCatch( diffexp(exp_persample,meta_persample), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
}

names(deg)[21]="Severe.Asthma.Chlam.w.Steroids"

openxlsx::write.xlsx(deg,paste0("deg/diseases_DEGs_2102.xlsx"),rowNames=T)


#celltype composition####
ff<-list.files(path="celltypeproportion", pattern="*celltypeproportion.csv", full.names=TRUE)
each<-sapply(ff,fread,simplify=FALSE,USE.NAMES = TRUE)
filenames=names(each)
filenames=strsplit(filenames,"/")
filenames=sapply(filenames,`[`,2)
filenames=strsplit(filenames,"_")
filenames <- lapply(filenames, function(x) head(x, -1))
names(each)=filenames
names(each) <- gsub("^c\\(|\\)$|\"|\\s", "", names(each))
names(each) <- gsub(",", "_", names(each))


for (e in 1:length(each)){
  ct=each[[e]]
  ct=ct[,-1]
  ct=ct[ct$Pseudo!="Removed",]
  ct$meta=ct$Pseudo
  ct$proportion=ct$Pseudo
  #remove rows that are 0 for all groups
  rows_to_remove <- numeric(0)
  for (i in 1:nrow(ct)){
    if(sum(ct[i,3:(ncol(ct)-3)])==0){
      rows_to_remove <- c(rows_to_remove, i)
    }
  }
  ct=as.data.frame(ct[-rows_to_remove,])
  
  for(n in 1:nrow(ct)){
    renaming=as.character(colnames(ct[n,3:(ncol(ct)-3)]>0)[apply(ct[n,3:(ncol(ct)-3)]>0, 2, function(x) any(x == TRUE))])
    prop=ct[n,renaming]
    ct$meta[n]=renaming
    ct$proportion[n]=prop;ct$proportion=as.numeric(ct$proportion)
    ct$proportion[n]=ct$proportion[n]/ct$`0`[n]
    
  }
  ct$average=ct$proportion
  
  df=data.frame()
  for(m in 1:length(unique(ct$meta))){
    for (c in 1:length(unique(ct$Diff_Exp))){
      average <- mean(ct[c(ct$meta==unique(ct$meta)[m] &ct$Diff_Exp==unique(ct$Diff_Exp)[c]),'proportion'])
      ct[which(ct$meta==unique(ct$meta)[m]&ct$Diff_Exp==unique(ct$Diff_Exp)[c]),'average']=average
      test=ct[ct$meta==unique(ct$meta)[m]&ct$Diff_Exp==unique(ct$Diff_Exp)[c],]
      test2=test[1,]
      
      df=rbind(df,test2)
    }
  }
  df=na.omit(df)
  cols=c('#ff9966', '#33ffcc', '#3366ff', '#8B4513', '#cc33ff',
         '#ff66cc', '#ccff66', '#6600ff', '#ffcc00', '#00ffcc',
         '#dd5544', '#77dd66', '#4477dd', '#ffbb55', '#55bbff',
         '#ff5599', '#99ff66', '#9955ff', '#ccaa66', '#8B0000',"red")
  
  p=ggplot(df, aes(x=meta, y=average, fill=Diff_Exp))+ 
    geom_col( width=0.6,position=position_stack(vjust=1))+
    theme(axis.title.x =element_blank(),
          axis.text.x=element_text(angle=45,hjust=1,size=14),
          axis.line = element_line(colour = 'black'))+
    labs(y="Average Cell Type Proportion")+
    scale_fill_manual(values=cols)+WhiteBackground()
  
  #ggsave(paste0("plots/",names(each)[e],"_celltype.png"),plot = p,width=2500,height=2000,units="px",dpi=300)
  write.csv(ct,paste0("celltypeproportion/",names(each)[e],"_celltypeproportion2.csv"))
}

#ct$meta=factor(ct$meta,levels=c("Young","Adult","Old"))
#write.csv(ct,'celltypeproportion/age_cont_celltypeproportion2.csv')
# 
# 
# dff2 = ct %>% arrange(desc(Diff_Exp)) %>% group_by(meta) %>% mutate(newy=cumsum(proportion))
# dff2$average=dff2$proportion
# df=data.frame()
# for(m in 1:length(unique(ct$meta))){
#   for (c in 1:length(unique(ct$Diff_Exp))){
#     test=dff2[dff2$meta==unique(ct$meta)[m]&dff2$Diff_Exp==unique(ct$Diff_Exp)[c],]
#     maxval=mean(test$proportion)
#     test2=cbind(test[1,],maxval);colnames(test2)[ncol(test2)]="average"
#     test3=as.data.frame(test2[,c("meta","Diff_Exp","newy","average")])
#     print(test3)
#     df=rbind(df,test3)
#   }
# }
# 
# ggplot(df, aes(x=meta, y=max, fill=Diff_Exp))+ 
#   geom_col( width=0.6,position=position_stack(vjust=1))



##Calculate standard error for error bars
# Calculate standard deviation
# ct$se=ct$proportion
# for(m in unique(ct$meta)){
#   for (c in unique(ct$Diff_Exp)){
#       standard_deviation <- sd(ct[c(ct$meta==m &ct$Diff_Exp==c),'proportion'])
#       sample_size <- length(ct[c(ct$meta==m&ct$Diff_Exp==c),"proportion"])
#       standard_error <- standard_deviation / sqrt(sample_size)
#       ct[which(ct$meta==m&ct$Diff_Exp==c),'se']=standard_error
#   }
# }
# # # Calculate the number of observations
# sample_size <- length(ct[c(ct$meta=='Adult'&ct$Diff_Exp=='Alveolar macrophages'),"proportion"])
# # Calculate standard error
# standard_error <- standard_deviation / sqrt(sample_size)


