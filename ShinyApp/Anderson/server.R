library(shiny)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(Seurat)
library(DT)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(DESeq2)
library(qs)

#navigate to the directory that contains the datasets:
setwd(paste0(dirname(getwd()),"/dataexplorer"))
retinalmgtimecoursedata<-read.table("vetterdata/mgtimecourse.txt",sep="\t",stringsAsFactors = FALSE)
#bax<-readRDS(file = "vetterdata/bax + plx seurat.rds")
bax<-qread(file = "vetterdata/bax + plx seurat.qs")
itgax<-read_tsv("Ghena/Itgax edgeR results.txt")
p4<-readRDS(file = "Anderson/axl mertk ko P4 deseq.rds")
#p7<-readRDS(file = "Anderson/axl mertk dko P7 deseq.rds")
#res_p7<-results(p7)
#counts_p7<-as.data.frame(assay(p7))
#meta_p7<-as.data.frame(colData(p7))
res_p7<-read_tsv("Anderson/axl mertk dko P7 results.txt")
counts_p7<-read.table("Anderson/axl mertk dko P7 counts.txt")
meta_p7<-read.table("Anderson/axl mertk dko P7 meta data.txt")

pal<-c("#ede9cf","#dbe7c0","#c1e3b8","#a3dcb5","#81d3b7","#60c8bc","#48bac2","#45a9c5","#5895c5","#707ebf","#8464b2","#91489d","#942a81")

shinyServer(function(input, output) {
  
  #print dataset selection at the top of the page ####
  output$dataheader<-renderUI({
    if (input$dataset==1)
    {HTML("retinal microglia timecourse")}
    else if (input$dataset==2){HTML("retinal microglia single cell")}
    else if (input$dataset==3){HTML("Cd11c lineage bulk seq")}
    else if (input$dataset==4){HTML("Axl Mertk KO bulk seq")}
  })
  
  #onegene tab ###
  #onegene tab legends ####
  
  output$legend1<-renderUI({
    if (input$dataset==1)
    {HTML("FPKM values \u00b1 sem. Significance is padj from DESeq2. *<=0.05, **<=0.01, ***<=0.001. 
            Results from Anderson et al., Cell Reports, 2019. 
          Brain microglia Fastq files from Matcovich Natan et al., Science, 2016")}
    else if (input$dataset==2){HTML("Single cell RNA-seq results from P7 sorted microglia. In the feature plot, darker dots indicate cells with higher relative expression of the feature/gene.")}
    else if (input$dataset==3){HTML("Bulk sequencing results from P5 CD11c-Cre-gfp;Rosa-TDTomato animals sorted for active Cd11c expression (GFP+TdTom+) or Cd11c lineage (GFP-TdTom+).
                                    Very few microglia are TdTom-. 
                                    Cells were combined into a single sample for each condition and differential expression analysis is done with EdgeR.")}
    else if (input$dataset==4){HTML("For each of the triplicate samples at P4 CD11b+ CD45+ Ly6c- cells from 4 retinas (two animals) were pooled. 
                                    P7 samples are CD11b+ CD45+ cells from two retinas (one animal). Published in Anderson et al. eLife, 2022. GSE192601. 
                                    Significance is padj from DESeq2. *<=0.05, **<=0.01, ***<=0.001.")}
  })
  
  #plot function for retinal MG timecourse ####
  bargraph<-function(dataset="retmg",conditions=3)({
    
    #return the appropriate asterisk for signifigance
    star<-function(table=mytable,i=1){
      if(!table[i,4]<=0.05){return("NS")}
      else if(table[i,4]<=0.001){return("***")}
      else if(table[i,4]<=0.01){return("**")}
      else{return("*")}
    }
    
    if(dataset %in% c("retmg","wholeret","brnmg")){
      dat<-retinalmgtimecoursedata[input$goi,]
      pal<-c("#81d3b7","#5895c5","#942A81")
    }

    #generated table of values to graph
    if(dataset=="retmg"){
      mytable<-data.frame(cond=c("e16.5","P7","P60"),
                          FPKM=c(dat[[35]],dat[[36]],dat[[37]]),
                          sem=c(dat[[38]],dat[[39]],dat[[40]]),
                          sig=c(dat[[21]],dat[[23]],dat[[22]]))
      title<-"retinal microglia"}
    else if (dataset=="wholeret"){
      mytable<-data.frame(cond=c("e16.5","P7","P60"),
                          FPKM=c(dat[[41]],dat[[30]],dat[[31]]),
                          sem=c(dat[[42]],0,0),
                          sig=c(1,1,1))
      title<-"whole retina"}
    else if (dataset=="brnmg"){
      mytable<-data.frame(cond=c("e16.5","P6","P60"),
                          FPKM=c(dat[[43]],dat[[44]],dat[[45]]),
                          sem=c(dat[[46]],dat[[47]],dat[[48]]),
                          sig=c(dat[[24]],dat[[25]],dat[[26]]))
      title<-"brain microglia"}
    
    mytable$cond<-factor(mytable$cond,levels=unique(mytable$cond))
    
    ggplot(mytable, aes(x=cond, y=FPKM)) + 
      geom_bar(stat="identity",color="black",fill=pal,width=0.6) +
      geom_errorbar(aes(ymin=FPKM-sem, ymax=FPKM+sem), width=.2)+
      ggtitle(paste0(input$goi,", ",title))+
      theme_bw()+
      theme(plot.title = element_text(size = 20),
            axis.text.x = element_text(size=16),
            axis.text.y = element_text(size=12),
            aspect.ratio=0.65,
            legend.title=element_text("centered"),
            axis.title.x=element_blank(),
            axis.title.y=element_text(size=16))+
      #add 'brackets' for significance and make space for them
      {if(mytable[3,4]<=0.05&(mytable[1,4]<=0.05|mytable[2,4]<=0.05))scale_y_continuous(limits=c(0,max(rowSums(mytable[,2:3]))*1.2))}+
      {if(mytable[3,4]<=0.05&mytable[1,4]>0.05&mytable[2,4]>0.05)scale_y_continuous(limits=c(0,max(rowSums(mytable[,2:3]))*1.1))}+
      {if(mytable[3,4]>0.05&(mytable[1,4]<=0.05|mytable[2,4]<=0.05))scale_y_continuous(limits=c(0,max(rowSums(mytable[,2:3]))*1.1))}+
      {if(mytable[1,4]<=0.05)geom_signif(stat="identity",data=data.frame(x=1, xend=1.95,y=max(rowSums(mytable[1:2,2:3]))*1.05,annotation=star()),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}+
      {if(mytable[2,4]<=0.05)geom_signif(stat="identity",data=data.frame(x=2.05, xend=3,y=max(rowSums(mytable[2:3,2:3]))*1.05,annotation=star(i=2)),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}+ 
      {if(mytable[3,4]<=0.05&mytable[1,4]>0.05&mytable[2,4]>0.05)geom_signif(stat="identity",data=data.frame(x=1, xend=3,y=max(rowSums(mytable[,2:3]))*1.05,annotation=star(i=3)),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}+
      {if(mytable[3,4]<=0.05&(mytable[1,4]<=0.05|mytable[2,4]<=0.05))geom_signif(stat="identity",data=data.frame(x=1, xend=3,y=max(rowSums(mytable[,2:3]))*1.15,annotation=star(i=3)),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}
  })
  
  star<-function(q=mytable[1,4]){
    if(!q<=0.05){return("NS")}
    else if(q<=0.001){return("***")}
    else if(q<=0.01){return("**")}
    else{return("*")}
  }
  
  #bar function specifically for the itgax edgeR results. Two conditions with only one sample each.
  bar_2<-function(dds=itgax,counts_type="counts", #specific edge R dataset
                     title="in Cd11c groups",order=c("lineage","active"),pal=c("#cc2a49","#f99e4c"))({ 
                       #lineage=red,active=yellow or green
                       dat<-filter(dds,gene == input$goi)
                       
                       as_tibble(dat[,3:5]) %>% 
                         pivot_longer(-gene,names_to = "trt") %>%
                         mutate(trt = factor(trt,levels=order)) %>%
                         ggplot(aes(x=trt, y=value)) + 
                         geom_bar(aes(fill=trt),stat="identity",position=position_dodge(),color="black",width=0.6) +
                         ggtitle(paste0(input$goi,", ",title))+
                         scale_fill_manual(values=pal)+
                         theme_bw()+
                         labs(y = counts_type)+
                         theme(axis.text.x = element_text(size=16),
                               plot.title = element_text(size = 20),
                               axis.text.y = element_text(size=12),
                               aspect.ratio=0.65,
                               legend.position = ('none'),
                               axis.title.x=element_blank(),
                               axis.title.y=element_text(size=16))+
                         #add 'brackets' for significance and make space for them
                         {if(dat[1,1]<=0.05)scale_y_continuous(limits=c(0,max(dat[1,4:5])*1.1))}+
                         {if(dat[1,1]<=0.05)geom_signif(stat="identity",data=data.frame(x=1, xend=2,y=max(dat[1,4:5])*1.05,
                                         annotation=star(dat[1,1])),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}
                     })
  
  volcano_gene <- function(dds=itgax){
    ggplot(dds,aes(-logFC,-log(PValue)))+
      geom_point(color="grey")+
      geom_point(data=filter(dds,gene==input$goi),color="#7a9200")+
      geom_text_repel(data=filter(dds,gene==input$goi),aes(label=gene))+
      theme_bw()+
      xlab("logFC Cd11c active/lineage")
      
  }
  
  bar_star<-function(dds=p4, #DESeq data object with conditions as trt
                     title="P4",order=c("WT","Mertk_KO","Axl_KO"),colors=pal) ({
                       
                       t <- length(order) #number of different conditions
                       
                       counts<-t(as.data.frame(assay(dds))[input$goi,]) #pulled counts or transformed data for goi
                       dat <- cbind(colData(dds),counts) #added metadata information
                       sig<-as.data.frame(results(dds,contrast = c("trt",order[c(1,2)]))) %>%
                         rownames_to_column("gene") %>%
                         select(gene,padj)
                      if(t==3) {sig <- sig %>% mutate(two = results(dds,contrast = c("trt",order[c(3,2)]))$padj) %>%
                         mutate(three = results(dds,contrast = c("trt",order[c(3,1)]))$padj)}
                      sig <- filter(sig,gene==input$goi)
                       
                       p<-as_tibble(dat) %>% 
                         mutate(trt = factor(trt,levels=order)) %>%
                         ggplot(aes(x=trt, y=.data[[input$goi]])) + 
                         geom_bar(aes(fill=id),stat="identity",position=position_dodge(),color="black",width=0.6) +
                         ggtitle(paste0(input$goi,", ",title))+
                         scale_fill_manual(values=c(colors))+
                         theme_bw()+
                         labs(y = "counts")+
                         theme(axis.text.x = element_text(size=16),
                               plot.title = element_text(size = 20),
                               axis.text.y = element_text(size=12),
                               aspect.ratio=0.65,
                               legend.position = ('none'),
                               axis.title.x=element_blank(),
                               axis.title.y=element_text(size=16))
                         #add 'brackets' for significance and make space for them
                       if(t == 2){p+
                           {if(sig[1,2]<=0.05)scale_y_continuous(limits=c(0,max(dat[[input$goi]])*1.1))}+
                           {if(sig[1,2]<=0.05)geom_signif(stat="identity",data=data.frame(x=1, xend=1.95,y=max(dat[dat$trt%in%order[c(1,2)],input$goi])*1.05,
                                                                                          annotation=star(sig[1,2])),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}
                       }   
                       if(t == 3){p+
                         {if(sig[1,4]<=0.05&(sig[1,2]<=0.05|sig[1,3]<=0.05))scale_y_continuous(limits=c(0,max(dat[[input$goi]])*1.2))}+
                         {if(sig[1,4]<=0.05&sig[1,2]>0.05&sig[1,3]>0.05)scale_y_continuous(limits=c(0,max(dat[[input$goi]])*1.1))}+
                         {if(sig[1,4]>0.05&(sig[1,2]<=0.05|sig[1,3]<=0.05))scale_y_continuous(limits=c(0,max(dat[[input$goi]])*1.1))}+
                         
                         {if(sig[1,2]<=0.05)geom_signif(stat="identity",data=data.frame(x=1, xend=1.95,y=max(dat[dat$trt%in%order[c(1,2)],input$goi])*1.05,
                                                                                        annotation=star(sig[1,2])),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}+
                         {if(sig[1,3]<=0.05)geom_signif(stat="identity",data=data.frame(x=2.05, xend=3,y=max(dat[dat$trt%in%order[c(2,3)],input$goi])*1.05,
                                                                                        annotation=star(sig[1,3])),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}+ 
                         {if(sig[1,4]<=0.05&sig[1,2]>0.05&sig[1,3]>0.05)geom_signif(stat="identity",data=data.frame(x=1, xend=3,y=max(dat[,input$goi])*1.05,
                                                                                                                    annotation=star(sig[1,4])),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}+
                         {if(sig[1,4]<=0.05&(sig[1,2]<=0.05|sig[1,3]<=0.05))geom_signif(stat="identity",data=data.frame(x=1, xend=3,y=max(dat[[input$goi]])*1.15,
                                                                                                                        annotation=star(sig[1,4])),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}
                       }
                      })
  
  #Idk why on earth, but it takes forever to run results(p7) so I just do it up top and use that df in this specific function
  bar_star_p7<-function(res=res_p7,cts=counts_p7,meta=meta_p7, #DESeq data object with conditions as trt
                     title="P7",order=c("WT","Mertk_Axl_dKO"),colors=pal) ({
                       
                       t <- length(order) #number of different conditions
                       
                       counts<-t(cts[input$goi,]) #pulled counts or transformed data for goi
                       dat <- cbind(meta,counts) #added metadata information
                       sig<-res %>%
                         select(gene,padj) %>%
                         filter(gene==input$goi)
                       
                       as_tibble(dat) %>% 
                         mutate(trt = factor(trt,levels=order)) %>%
                         ggplot(aes(x=trt, y=.data[[input$goi]])) + 
                         geom_bar(aes(fill=id),stat="identity",position=position_dodge(),color="black",width=0.6) +
                         ggtitle(paste0(input$goi,", ",title))+
                         scale_fill_manual(values=c(colors))+
                         theme_bw()+
                         labs(y = "counts")+
                         theme(axis.text.x = element_text(size=16),
                               plot.title = element_text(size = 20),
                               axis.text.y = element_text(size=12),
                               aspect.ratio=0.65,
                               legend.position = ('none'),
                               axis.title.x=element_blank(),
                               axis.title.y=element_text(size=16))+
                       #add 'brackets' for significance and make space for them
                           {if(sig[1,2]<=0.05)scale_y_continuous(limits=c(0,max(dat[[input$goi]])*1.1))}+
                           {if(sig[1,2]<=0.05)geom_signif(stat="identity",data=data.frame(x=1, xend=1.95,y=max(dat[dat$trt%in%order[c(1,2)],input$goi])*1.05,
                                                                                          annotation=star(sig[1,2])),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}
                     })
  
  volcano_deseq <- function(dds=p4,con=c("Axl_KO","WT"),xlim=c(-1.8,1.8),ylim=c(0,60),col="#7ad151"){ #should work with any deseq object
    x<-as.data.frame(results(dds,contrast = c("trt",con))) %>%
      rownames_to_column("gene")
    ggplot(x,aes(log2FoldChange,-log(padj)))+
      geom_point(color="grey")+
      geom_point(data=filter(x,gene==input$goi),color=col)+
      geom_text_repel(data=filter(x,gene==input$goi),aes(label=gene))+
      theme_bw()+
      xlab(paste0("logFC ",con[1],"/",con[2])) +
    {if (!identical(xlim,c(NA,NA))) scale_x_continuous(limits=xlim, expand = c(0, 0))} +
    {if (!identical(ylim,c(NA,NA)))  scale_y_continuous(limits=ylim, expand = c(0, 0))}
    
  }
  volcano_p7 <- function(x=res_p7,con=c("logFC Mertk Axl dKO/WT"),xlim=c(NA,NA),ylim=c(NA,NA),col="#fbe829"){
    ggplot(x,aes(-log2FoldChange,-log(padj)))+
      geom_point(color="grey")+
      geom_point(data=filter(x,gene==input$goi),color=col)+
      geom_text_repel(data=filter(x,gene==input$goi),aes(label=gene))+
      theme_bw()+
      xlab(con) +
      {if (!identical(xlim,c(NA,NA))) scale_x_continuous(limits=xlim, expand = c(0, 0))} +
      {if (!identical(ylim,c(NA,NA)))  scale_y_continuous(limits=ylim, expand = c(0, 0))}
    
  }

  #generate boxes for gene tab depending on selected dataset ####
  output$geneboxes<-renderUI({
    if(input$dataset==1 & input$goi%in%rownames(retinalmgtimecoursedata)){
      output$plot2<-renderPlot({bargraph()})
      output$plot1<-renderPlot({bargraph("wholeret")})
      output$plot3<-renderPlot({bargraph("brnmg")})
      fluidRow(
        box(plotOutput("plot1"),width=4),
        box(plotOutput("plot2"),width=4),
        box(plotOutput("plot3"),width=4)
      )}
    else if(input$dataset==3 & input$goi%in%itgax$gene){
      output$itgax_bar<-renderPlot({bar_2()})
      output$itgax_volcano<-renderPlot({volcano_gene()})
      fluidRow(
        box(plotOutput("itgax_bar"),width=4),
        box(plotOutput("itgax_volcano"),width=6)
      )}
    else if (input$dataset==2 & (input$goi%in%(c(rownames(bax@assays$RNA@data),names(bax@meta.data))))){
      output$orig_ident<-renderPlot({DimPlot(bax, group.by = "orig.ident")})
      output$clust<-renderPlot({DimPlot(bax, label = TRUE) + NoLegend()})
      output$feature<-renderPlot({FeaturePlot(object = bax,input$goi, cols=c(pal),pt.size=1.5) + NoLegend()+theme(axis.text.x = element_text(angle = 45, hjust=1))})
      output$violin<-renderPlot({VlnPlot(object = bax,input$goi) + NoLegend()})
      output$clust_wt<-renderPlot({DimPlot(bax, label = TRUE,cells= grep("wt", Cells(bax), value = T)) + NoLegend()})
      output$clust_ko<-renderPlot({DimPlot(bax, label = TRUE,cells= grep("ko", Cells(bax), value = T)) + NoLegend()})
      output$clust_veh<-renderPlot({DimPlot(bax, label = TRUE,cells= grep("veh", Cells(bax), value = T)) + NoLegend()})
      output$clust_plx<-renderPlot({DimPlot(bax, label = TRUE,cells= grep("plx", Cells(bax), value = T)) + NoLegend()})
      output$feat_wt<-renderPlot({FeaturePlot(object = bax,input$goi, cols=c(pal),pt.size=1.5,cells= grep("wt", Cells(bax))) + NoLegend()})
      output$feat_ko<-renderPlot({FeaturePlot(object = bax,input$goi, cols=c(pal),pt.size=1.5,cells= grep("ko", Cells(bax))) + NoLegend()})
      output$feat_veh<-renderPlot({FeaturePlot(object = bax,input$goi, cols=c(pal),pt.size=1.5,cells= grep("veh", Cells(bax))) + NoLegend()})
      output$feat_plx<-renderPlot({FeaturePlot(object = bax,input$goi, cols=c(pal),pt.size=1.5,cells= grep("plx", Cells(bax))) + NoLegend()})
      fluidRow(
        box(plotOutput("orig_ident"),width=3),
        box(plotOutput("clust"),width=3),
        box(plotOutput("feature"),width=3),
        box(plotOutput("violin"),width=3),
        box("bax ko", plotOutput("feat_ko"),width=3),
        box("bax wt", plotOutput("feat_wt"),width=3),
        box("vehicle", plotOutput("feat_veh"),width=3),
        box("plx", plotOutput("feat_plx"),width=3),
        box("bax ko", plotOutput("clust_ko"),width=3),
        box("bax wt", plotOutput("clust_wt"),width=3),
        box("vehicle", plotOutput("clust_veh"),width=3),
        box("plx", plotOutput("clust_plx"),width=3)
      )}
    else if(input$dataset==4 & input$goi%in%c(rownames(p4),res_p7$gene)){
      output$p4_bar<-renderPlot({bar_star(colors=c("#22a884","#660066",rep("#22a884",2),rep("#7ad151",3),rep("#660066",2)))})
      output$p7_bar<-renderPlot({bar_star_p7(colors=c(rep("#660066",3),rep("#fbe829",3)))})
      output$dko_v<-renderPlot({volcano_p7()})
      output$axl_v1<-renderPlot({volcano_deseq(xlim=c(NA,NA),ylim=c(NA,NA))})
      output$axl_v2<-renderPlot({volcano_deseq()})
      output$mertk_v1<-renderPlot({volcano_deseq(con=c("Mertk_KO","WT"),xlim=c(NA,NA),ylim=c(NA,NA),col="#22a884")})
      output$mertk_v2<-renderPlot({volcano_deseq(con=c("Mertk_KO","WT"),col="#22a884")})
      fluidRow(
        box(plotOutput("p4_bar"),width=4),
        box(plotOutput("p7_bar"),width=4),
        box(plotOutput("dko_v"),width=4),
        box(plotOutput("mertk_v1"),width=3),
        box(plotOutput("mertk_v2"),width=3),
        box(plotOutput("axl_v1"),width=3),
        box(plotOutput("axl_v2"),width=3)
      )}
    })
  })

