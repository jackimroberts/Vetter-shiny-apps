library(shiny)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(qs)
source("../plots.R")

ug_time<-qread("../data/mgtimecourse.qs")
p7_sc<-qread("../data/bax+plx_tiny.qs")
#itgax<-read_tsv("Ghena/Itgax edgeR results.txt")
tam_p4<-qread("../data/axl mertk ko P4 deseq.qs")
tam_p7<-qread("../data/axl mertk dko P7 deseq.qs")

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
    else if (input$dataset==2){HTML("Single cell RNA-seq results from P7 sorted microglia. In the feature plot, darker dots indicate cells with higher relative expression of the feature/gene. <br> 
                                    This is a showing a VERY small dataset with ~1/5 of total cells and the more highly expressed genes. <br>
                                    Published in Anderson et al. eLife, 2022. GSE192601.")}
    else if (input$dataset==3){HTML("Bulk sequencing results from P5 CD11c-Cre-gfp;Rosa-TDTomato animals sorted for active Cd11c expression (GFP+TdTom+) or Cd11c lineage (GFP-TdTom+).
                                    Very few microglia are TdTom-. 
                                    Cells were combined into a single sample for each condition and differential expression analysis is done with EdgeR.")}
    else if (input$dataset==4){HTML("For each of the triplicate samples at P4 CD11b+ CD45+ Ly6c- cells from 4 retinas (two animals) were pooled. 
                                    P7 samples are CD11b+ CD45+ cells from two retinas (one animal). Significance is padj from DESeq2. *<=0.05, **<=0.01, ***<=0.001.")}
  })
  
  #generate boxes for gene tab depending on selected dataset ####
  output$geneboxes<-renderUI({
    if(input$dataset==1 & (input$goi %in% ug_time$gene.name)){
      output$plot2<-renderPlot({barplot_replicates(
        mutate(ug_time[,c(27:31,1)],name=gene.name),
        "whole retina",
        conditions=c("e16","P7","P60"),
        pallet=c("#81d3b7","#942A81","#5895c5"),gene=input$goi)})
      output$plot1<-renderPlot({barplot_conditions(
        mutate(ug_time[,c(21,23,22,35:40,1)],name=gene.name),
        "retinal microglia",
        conditions=c("e16","P7","P60"),
        pallet=c("#81d3b7","#942A81","#5895c5"),gene=input$goi)})
      output$rm_rep<-renderPlot({barplot_replicates(
        mutate(ug_time[,c(4:12,1)],name=gene.name),
        "retinal microglia",
        conditions=c("e16","P7","P60"),
        pallet=c("#81d3b7","#942A81","#5895c5"),gene=input$goi)})
      output$plot3<-renderPlot({barplot_conditions(
        mutate(ug_time[,c(24:26,43:48,1)],name=gene.name),
        "brain microgila",
        conditions=c("e16","P6","P60"),
        pallet=c("#81d3b7","#5895c5","#942A81"),gene=input$goi)})
      output$bm_rep<-renderPlot({barplot_replicates(
        mutate(ug_time[,c(14:19,1)],name=gene.name),
        "brain microgila",
        conditions=c("e16","P6","P60"),
        pallet=c("#81d3b7","#5895c5","#942A81"),gene=input$goi)})
      fluidRow(
        box(plotOutput("plot3"),width=4),
        box(plotOutput("plot1"),width=4),
        box(plotOutput("plot2"),width=4),
        box(plotOutput("bm_rep"),width=4),
        box(plotOutput("rm_rep"),width=4)
      )}
    else if (input$dataset==2 & (input$goi%in%(c(rownames(p7_sc@assays$SCT@data),names(p7_sc@meta.data))))){
      output$orig_ident<-renderPlot({DimPlot(p7_sc, group.by = "orig.ident")})
      output$clust<-renderPlot({DimPlot(p7_sc, label = TRUE) + NoLegend()})
      output$feature<-renderPlot({FeaturePlot(object = p7_sc,input$goi,pt.size=1.5) + 
          NoLegend()+theme(axis.text.x = element_text(angle = 45, hjust=1))+scale_color_gradientn(colors=pal)})
      output$violin<-renderPlot({VlnPlot(object = p7_sc,input$goi) + NoLegend()})
      output$clust_wt<-renderPlot({DimPlot(p7_sc, label = TRUE,cells= grep("wt", Cells(p7_sc), value = T)) + NoLegend()})
      output$clust_ko<-renderPlot({DimPlot(p7_sc, label = TRUE,cells= grep("ko", Cells(p7_sc), value = T)) + NoLegend()})
      output$clust_veh<-renderPlot({DimPlot(p7_sc, label = TRUE,cells= grep("veh", Cells(p7_sc), value = T)) + NoLegend()})
      output$clust_plx<-renderPlot({DimPlot(p7_sc, label = TRUE,cells= grep("plx", Cells(p7_sc), value = T)) + NoLegend()})
      output$feat_wt<-renderPlot({FeaturePlot(object = p7_sc,input$goi, pt.size=1.5,cells= grep("wt", Cells(p7_sc))) + NoLegend()+ scale_color_gradientn(colors=pal)})
      output$feat_ko<-renderPlot({FeaturePlot(object = p7_sc,input$goi, pt.size=1.5,cells= grep("ko", Cells(p7_sc))) + NoLegend()+ scale_color_gradientn(colors=pal)})
      output$feat_veh<-renderPlot({FeaturePlot(object = p7_sc,input$goi, pt.size=1.5,cells= grep("veh", Cells(p7_sc))) + NoLegend()+ scale_color_gradientn(colors=pal)})
      output$feat_plx<-renderPlot({FeaturePlot(object = p7_sc,input$goi, pt.size=1.5,cells= grep("plx", Cells(p7_sc))) + NoLegend()+ scale_color_gradientn(colors=pal)})
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
    else if(input$dataset==4 & input$goi%in%c(tam_p4$name,tam_p7$name)){
      output$p4_bar<-renderPlot({barplot_replicates(tam_p4[,1:10],"P4",
                                                    c("WT","Mertk KO","Axl KO"),
                                                    c("#7ad151","#22a884","#660066"),gene=input$goi)})
      output$p4_cond<-renderPlot({barplot_conditions(tam_p4[,c(20,18,1,10:16)],"P4",c("wt","mko","ako"),c("#7ad151","#22a884","#660066"),
                                                     c("WT","Mertk KO","Axl KO"),gene=input$goi)})
      output$p7_bar<-renderPlot({barplot_replicates(tam_p7[,1:7],"P7",
                                                    c("WT","Mertk_Axl_dKO"),
                                                    c("#fbe829","#660066"),gene=input$goi)})
      output$p7_cond<-renderPlot({barplot_conditions(tam_p7[,c(9:13,7)],"P7",c("wt","dko"),c("#fbe829","#660066"),c("WT","Mertk_Axl_dKO"),gene=input$goi)})
      output$dko_v<-renderPlot({volcano_gene(tam_p7,x_label="logFC Mertk Axl dKO/WT",col="#fbe829",p="padj",fc="l2fc",gene=input$goi)})
      output$axl_v1<-renderPlot({volcano_gene(tam_p4,x_label="logFC Axl KO/WT",col="#7ad151",p="q_ako_wt",fc="l2fc_ako_wt",gene=input$goi)+
          geom_hline(yintercept = 60,linetype="dashed",color="grey50")+
          geom_vline(xintercept = c(-1.8,1.8),linetype="dashed",color="grey50")})
      output$axl_v2<-renderPlot({volcano_gene(tam_p4,x_label="logFC Axl KO/WT",col="#7ad151",p="q_ako_wt",fc="l2fc_ako_wt",xlim=c(-1.8,1.8),ylim=c(0,60),gene=input$goi)})
      output$mertk_v1<-renderPlot({volcano_gene(tam_p4,x_label="logFC Mertk KO/WT",col="#22a884",p="q_mko_wt",fc="l2fc_mko_wt",gene=input$goi)+
          geom_hline(yintercept = 60,linetype="dashed",color="grey50")+
          geom_vline(xintercept = c(-1.8,1.8),linetype="dashed",color="grey50")})
      output$mertk_v2<-renderPlot({volcano_gene(tam_p4,x_label="logFC Mertk KO/WT",col="#22a884",p="q_mko_wt",fc="l2fc_mko_wt",xlim=c(-1.8,1.8),ylim=c(0,60),gene=input$goi)})
      fluidRow(
        box(plotOutput("p4_bar"),width=3),
        box(plotOutput("p4_cond"),width=3),
        box(plotOutput("p7_bar"),width=3),
        box(plotOutput("p7_cond"),width=3),
        box(plotOutput("mertk_v1"),width=3),
        box(plotOutput("axl_v1"),width=3),
        box(plotOutput("dko_v"),width=4),
        box(plotOutput("mertk_v2"),width=3),
        box(plotOutput("axl_v2"),width=3)
      )}
    })
  })

