library(shiny)
library(ggplot2)
library(ggsignif)
library(reshape2)
library(Seurat)
library(DT)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)

#navigate to the directory that contains the datasets:
setwd(paste0(dirname(getwd()),"/dataexplorer"))
ja2<-readRDS(file="vetterdata/jarid2_P0_n2.RDS") %>% ungroup()
ja3<-readRDS(file="vetterdata/jarid2_P0_n3.RDS") %>% ungroup()
ezh<-readRDS(file="vetterdata/ezh2_e16.RDS") %>% ungroup()
je16<-readRDS(file="vetterdata/jarid2_e16.RDS") %>% ungroup()
scj2<-readRDS(file="vetterdata/scJarid2 v2.rds")
clark<-readRDS(file="otherdata/clark_minimal_data.rds")
ctx<-read_tsv("Roberts/jarid2 edgeR results.txt")

pal<-c("#ede9cf","#dbe7c0","#c1e3b8","#a3dcb5","#81d3b7","#60c8bc","#48bac2","#45a9c5","#5895c5","#707ebf","#8464b2","#91489d","#942a81")
j2pal<-brewer.pal(11,"Spectral")
scj2pal <- c("seagreen","yellow2","darkorange","darkorchid4","chartreuse3","midnightblue","maroon3","grey90","grey20")
clarkpal<-c("darkorange","red","chartreuse3","seagreen","skyblue","black","yellow2","maroon3","navy","darkorchid4")

shinyServer(function(input, output) {
  
  #print dataset selection at the top of the page ####
  output$dataheader<-renderUI({
    if (input$dataset==4){HTML("PRC2 conditional knock-out in developing retina")}
    else if (input$dataset==5){HTML("jarid2 e18 retina, single cell")}
  })
  
  #tab legends ####
  
  output$legend1<-renderUI({
    if (input$dataset==4){HTML("FPKM values \u00b1 sem. Significance is padj from DESeq2. *<=0.05, **<=0.01, ***<=0.001. 
    Ezh2 results from Zhang et al., Dev. Bio., 2015. <br> <br> Jarid2 cKO cortex data is bulk RNA-seq from single samples 
                               of Jarid2 fl/fl;Emx1 cre/+ (cKO) and Jarid2 +/+;Emx1 cre/+ (control) at E15.5, targeting the
                               S1 cortical region. E15.5 would be considered a late cortical progenitor stage. 
                               Differential expression analysis done with EdgeR.")}
    else if (input$dataset==5){HTML("Jarid2 conditional knockout. Single cell RNA-seq results from e18 whole retina.")}
  })
  
  #plot function three replicates ####
  bar_rep<-function(table="scrry2",title="sequencing round 2",replicates=3)({
    
    if(table=="scrry2"){dat<-scrry[scrry$name==input$goi,]
    mytable<-data.frame("type"=c("YFP","sCrry"),"one"=c(dat[[4]],dat[[7]]),"two"=c(dat[[5]],dat[[8]]),"three"=c(dat[[6]],dat[[9]]))
    pal<-c("#007284","#00546d")}
    
    else if(table=="scrry3"){dat<-scrry3 %>% filter(gene_name==input$goi)
    mytable<-data.frame("type"=c("naive","YFP","sCrry"),"one"=c(dat[[9]],dat[[6]],dat[[3]]),"two"=c(dat[[10]],dat[[7]],dat[[4]]),"three"=c(dat[[11]],dat[[8]],dat[[5]]))
    pal<-c("#f6ffaa","#007284","#00546d")}
    
    else if (table=="ja3"){dat<-filter(ja3,name==input$goi)
    mytable<-data.frame("type"=c("control","Jarid2 cKO"),"one"=c(dat[[1]],dat[[4]]),"two"=c(dat[[2]],dat[[5]]),"three"=c(dat[[3]],dat[[6]]))
    pal<-c(j2pal[11],j2pal[3])}
    
    else if (table=="ezh"){dat<-filter(ezh,name==input$goi)
    mytable<-data.frame("type"=c("control","Ezh2 cKO"),"one"=c(dat[[9]],dat[[5]]),"two"=c(dat[[10]],dat[[6]]),"three"=c(dat[[11]],dat[[7]]),"four"=c(dat[[12]],dat[[8]]))
    pal<-c(j2pal[9],j2pal[1])}

    else if (table=="je16"){dat<-filter(je16,name==input$goi)
    mytable<-data.frame("type"=c("control","Jarid2 cKO"),"one"=c(dat[[1]],dat[[3]]),"two"=c(dat[[2]],dat[[4]]),"three"=c(dat[[6]],dat[[5]]))
    pal<-c(j2pal[10],j2pal[2])}
    
    dat<-melt(mytable)
    dat$a<-letters[seq( from = 1, to = nrow(dat) )]
    dat$type<-factor(dat$type,levels=unique(dat$type))
    
    ggplot(dat, aes(x=type, y=value,fill=a)) + 
      geom_bar(stat="identity",position=position_dodge(),color="black",width=0.6) +
      ggtitle(paste0(input$goi,", ",title))+
      scale_fill_manual(values=rep(pal,replicates))+
      theme_bw()+
      labs(y = "FPKM")+
      theme(axis.text.x = element_text(size=16),
            plot.title = element_text(size = 20),
            axis.text.y = element_text(size=12),
            aspect.ratio=0.65,
            legend.position = ('none'),
            axis.title.x=element_blank(),
            axis.title.y=element_text(size=16))
  })

  bar_2cond<-function(file="ja2")({
    
    #return the appropriate asterisk for signifigance
    star<-function(table=mytable,i=1,j=4){
      if(!table[i,j]<=0.05){return("NS")}
      else if(table[i,j]<=0.001){return("***")}
      else if(table[i,j]<=0.01){return("**")}
      else{return("*")}
    }
    
    if (file=="ja2"){
      dat<-filter(ja2,name==input$goi)
      mytable<-data.frame(cond=c("control","Jarid2 cKO"),
                        FPKM=c(dat[[7]],dat[[9]]),
                        sem=c(dat[[8]],dat[[10]]),
                        sig=c(dat[[5]],0))
      pal<-c(j2pal[11],j2pal[3])
      tit<-", P0, n=2"
    }
    else if (file=="ja3"){
      dat<-filter(ja3,name==input$goi)
      mytable<-data.frame(cond=c("control","Jarid2 cKO"),
                          FPKM=c(dat[[9]],dat[[11]]),
                          sem=c(dat[[10]],dat[[12]]),
                          sig=c(dat[[7]],0))
      pal<-c(j2pal[11],j2pal[3])
      tit<-", P0, n=3"
    }
    else if (file=="ezh"){
      dat<-filter(ezh,name==input$goi)
      mytable<-data.frame(cond=c("control","Ezh2 cKO"),
                          FPKM=c(dat[[13]],dat[[15]]),
                          sem=c(dat[[14]],dat[[16]]),
                          sig=c(dat[[4]],0))
      pal<-c(j2pal[9],j2pal[1])
      tit<-", e16.5, n=4"
    }
    else if (file=="je16"){
      dat<-filter(je16,name==input$goi)
      mytable<-data.frame(cond=c("control","Jarid2 cKO"),
                          FPKM=c(dat[[9]],dat[[11]]),
                          sem=c(dat[[10]],dat[[12]]),
                          sig=c(dat[[7]],0))
      pal<-c(j2pal[10],j2pal[2])
      tit<-", e16.5, n=3"
    }

    
    mytable$cond<-factor(mytable$cond,levels=unique(mytable$cond))
    
      ggplot(mytable, aes(x=cond, y=FPKM)) + 
      geom_bar(stat="identity",color="black",fill=pal,width=0.6) +
      geom_errorbar(aes(ymin=FPKM-sem, ymax=FPKM+sem), width=.2)+
      ggtitle(paste0(input$goi,tit))+
      theme_bw()+
      theme(plot.title = element_text(size = 20),
            axis.text.x = element_text(size=16),
            axis.text.y = element_text(size=12),
            aspect.ratio=0.65,
            legend.title=element_text("centered"),
            axis.title.x=element_blank(),
            axis.title.y=element_text(size=16))+
      {if(mytable[1,4]<=0.05)scale_y_continuous(limits=c(0,max(rowSums(mytable[,2:3]))*1.1))}+
      {if(mytable[1,4]<=0.05)geom_signif(stat="identity",data=data.frame(x=1, xend=2,y=max(rowSums(mytable[,2:3]))*1.05,annotation=star()),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}
  })

  star<-function(q=mytable[1,4]){
    if(!q<=0.05){return("NS")}
    else if(q<=0.001){return("***")}
    else if(q<=0.01){return("**")}
    else{return("*")}
  }  
  #bar function specifically for the edgeR results. Two conditions with only one sample each.
  bar_2<-function(dds=ctx,counts_type="counts", #specific edge R dataset
                  title="in cortex",order=c("con","cko"),pal=c("#cc2a49","#f99e4c"))({ 
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
  
  volcano_gene <- function(dds=ctx){
    ggplot(dds,aes(-logFC,-log(PValue)))+
      geom_point(color="grey")+
      geom_point(data=filter(dds,gene==input$goi),color="#7a9200")+
      geom_text_repel(data=filter(dds,gene==input$goi),aes(label=gene))+
      theme_bw()+
      xlab("logFC cortex Jarid2 control/cKO")
  }
  
  #generate boxes for gene tab depending on selected dataset ####
  output$geneboxes<-renderUI({
     if (input$dataset==4){
      output$plot6<-renderPlot({bar_2cond("je16")})
      output$plot1<-renderPlot({bar_2cond()})
      output$plot2<-renderPlot({bar_2cond("ja3")})
      output$plot3<-renderPlot({bar_2cond("ezh")})
      output$plot7<-renderPlot({bar_rep("je16",title="e16.5")})
      output$plot4<-renderPlot({bar_rep("ja3",title="P0")})
      output$plot5<-renderPlot({bar_rep("ezh",title="e16.5",replicates=4)})
      output$plot8<-renderPlot({bar_2()}) #cortex data
      output$plot9<-renderPlot({volcano_gene()}) #cortex data
      fluidRow(
        {if (input$goi %in% ezh$name)box(plotOutput("plot3"),width=3)},#ezh2
        {if (input$goi %in% je16$name)box(plotOutput("plot6"),width=3)},#j2 e16        
        {if (input$goi %in% ja2$name)box(plotOutput("plot1"),width=3)},#j2 p0
        {if (input$goi %in% ctx$gene)box(plotOutput("plot8"),width=3)},#j2 ctx
        {if (input$goi %in% ezh$name)box(plotOutput("plot5"),width=3)},
        {if (input$goi %in% je16$name)box(plotOutput("plot7"),width=3)},
        {if (input$goi %in% ja3$name)box(plotOutput("plot4"),width=3)},
        {if (input$goi %in% ctx$gene)box(plotOutput("plot9"),width=3)},
        {if (input$goi %in% ja3$name)box(plotOutput("plot2"),width=4)}
      )
    }
    else if (input$dataset==5 & (input$goi%in%(c(rownames(scj2@assays$RNA@data),names(scj2@meta.data))))){
      output$clust<-renderPlot({DimPlot(scj2, label = TRUE,cols=scj2pal) + NoLegend()})
      output$con<-renderPlot({FeaturePlot(object = scj2,input$goi, cols=c(pal),pt.size=1.5,cells= grep("con", Cells(scj2))) + NoLegend()})
      output$cko<-renderPlot({FeaturePlot(object = scj2,input$goi, cols=c(pal),pt.size=1.5,cells= grep("cko", Cells(scj2))) + NoLegend()})
      fluidRow(
        {box("all",plotOutput("clust"),width=4)},
        {box("Jarid2 control", plotOutput("con"),width=4)},
        {box("Jarid2 cko", plotOutput("cko"),width=4)}
      )
    }
  })
})
