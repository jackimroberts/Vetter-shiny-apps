library(shiny)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(qs)
source("../plots.R")

#navigate to the directory that contains the datasets:
ja2<-qread("../data/jarid2_P0_n2.qs")
ja3<-qread("../data/jarid2_P0_n3.qs")
ezh<-qread("../data/ezh2_e16.qs")
je16<-qread("../data/jarid2_e16.qs")
scj2<-qread("../data/scJarid2mini.qs")
ctx<-qread("../data/jarid2 edgeR results.qs")

j2pal<-brewer.pal(11,"Spectral")
scj2pal <- c("seagreen","yellow2","darkorange","darkorchid4","chartreuse3","midnightblue","maroon3","grey90","grey20")

shinyServer(function(input, output) {
  
  #print dataset selection at the top of the page ####
  output$dataheader<-renderUI({
    if (input$dataset==4){HTML("PRC2 conditional knock-out in developing retina")}
    else if (input$dataset==5){HTML("Jarid2 e18.5 retina, single cell")}
  })
  
  #tab legends ####
  
  output$legend1<-renderUI({
    if (input$dataset==4){HTML("Bulk RNA-seq of whole retina. Presented data are FPKM values \u00b1 sem. Significance is padj from DESeq2. *<=0.05, **<=0.01, ***<=0.001. 
    Ezh2 results from Zhang et al., Dev. Bio. 2015 and Zhang et al., Cell Reports 2023. <br><br> Jarid2 cKO cortex data is bulk RNA-seq from single samples 
                               of Jarid2 fl/fl;Emx1 cre/+ (cKO) and Jarid2 +/+;Emx1 cre/+ (control) at E15.5,
                               S1 cortical region. <br> E15.5 would be considered a late cortical progenitor stage. 
                               Differential expression analysis done with EdgeR. Unpublished.")}
    else if (input$dataset==5){HTML("Jarid2 conditional knockout. Single cell RNA-seq results from e18 whole retina. <br>
                                    These plots only show SCT data (not RNA) and represent only half of the cells to save file space. <br>
                                    Full data in Zhang et al., Cell Reports 2023.")}
  })

  #generate boxes for gene tab depending on selected dataset ####
  output$geneboxes<-renderUI({
     if (input$dataset==4){
      output$plot6<-renderPlot({barplot_conditions(je16[,c(7:13)],"e16.5, n=3",
                                          label=c("control","Jarid2 cKO"),
                                          pallet=c(j2pal[10],j2pal[2]))})
      output$plot1<-renderPlot({bar_2cond(ja2[,c(5:11)],
                                          title="P0, n=2",
                                          label=c("control","Jarid2 cKO"),
                                          pallet=c(j2pal[11],j2pal[3]))})
      output$plot2<-renderPlot({bar_2cond(je16[,c(7:13)],"P0, n=3",
                                          label=c("control","Jarid2 cKO"),
                                          pallet=c(j2pal[11],j2pal[3]))})
      output$plot3<-renderPlot({bar_2cond(ezh[,c(4,2,13:16)],"e16.5, n=4",
                                          label=c("control","Ezh2 cKO"),
                                          pallet=c(j2pal[9],j2pal[1]))})
      output$plot7<-renderPlot({barplot_replicates(je16[,c(1:6,13)],"e16.5",
                                        label=c("control","Jarid2 cKO"),
                                        pallet=c(j2pal[10],j2pal[2]))})
      output$plot4<-renderPlot({barplot_replicates(ja3[,c(1:6,13)],"P0",
                                        label=c("control","Jarid2 cKO"),
                                        pallet=c(j2pal[11],j2pal[3]))})
      output$plot5<-renderPlot({barplot_replicates(ezh[,c(2,5:12)],"e16.5",
                                        label=c("control","Ezh2 cKO"),
                                        pallet=c(j2pal[9],j2pal[1]))})
      output$plot8<-renderPlot({bar_2cond(ctx[,c(-1)]%>%mutate(cko_sem=0,con_sem=0,con_mean=con,cko_mean=cko,name=gene),
                                          "e15.5, cortex",
                                          label=c("control","Jarid2 cKO"),
                                          pallet=c("#cc2a49","#f99e4c"))}) #cortex data
      output$plot9<-renderPlot({volcano_gene()}) #cortex data
      fluidRow(
        {if (input$goi %in% ezh$name)box(plotOutput("plot3"),width=3)},#ezh2
        {if (input$goi %in% je16$name)box(plotOutput("plot6"),width=3)},#j2 e16        
        {if (input$goi %in% ja2$name)box(plotOutput("plot1"),width=3)},#j2 p0
        {if (input$goi %in% ja3$name)box(plotOutput("plot2"),width=3)},
        {if (input$goi %in% ezh$name)box(plotOutput("plot5"),width=4)},
        {if (input$goi %in% je16$name)box(plotOutput("plot7"),width=4)},
        {if (input$goi %in% ja3$name)box(plotOutput("plot4"),width=4)},
        {if (input$goi %in% ctx$gene)box(plotOutput("plot8"),width=3)},#j2 ctx
        {if (input$goi %in% ctx$gene)box(plotOutput("plot9"),width=3)}
        
      )
    }
    else if (input$dataset==5 & (input$goi%in%(c(rownames(scj2@assays$SCT@data),names(scj2@meta.data))))){
      output$clust<-renderPlot({DimPlot(scj2, label = TRUE,cols=scj2pal) + NoLegend()})
      output$con<-renderPlot({FeaturePlot(object = scj2,input$goi, pt.size=1.5,cells= grep("con", Cells(scj2))) + NoLegend()})
      output$cko<-renderPlot({FeaturePlot(object = scj2,input$goi, pt.size=1.5,cells= grep("cko", Cells(scj2))) + NoLegend()})
      fluidRow(
        {box("all",plotOutput("clust"),width=4)},
        {box("Jarid2 control", plotOutput("con"),width=4)},
        {box("Jarid2 cko", plotOutput("cko"),width=4)}
      )
    }
  })
})
