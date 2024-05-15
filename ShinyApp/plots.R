library(ggplot2)
library(ggsignif)
library(ggrepel)
library(tidyverse)
library(RColorBrewer)

#plot n replicates, length(conditions) ####
#assumes gene "name" column
#and ony ave columns labeled with 'conditions'
#no mean deviation or significance
barplot_replicates<-function(table=ja3[,c(1:6,13)],
                  title="P0",
                  conditions=c("con","cko"),
                  pallet=brewer.pal(11,"Spectral"),
                  label=conditions,
                  gene=input$goi)({
                    
                    #get rid of extra columns
                    table<-table[,c(grep(paste(conditions,collapse="|"),names(table),value = T),"name")]
                    #manipulate table to contain only the gene of interest
                    to_plot<- table %>% 
                      filter(name==gene) %>%
                      pivot_longer(!name, names_to = "cond", values_to = "value") %>%
                      mutate(lab=label[1],col=pallet[1])
                    #adjust layers and colors for multiple conditions
                    for(i in 2:length(conditions)){
                      idx <- grepl(conditions[i],to_plot$cond)
                      to_plot<- to_plot %>%
                        mutate( lab=replace(lab,idx,label[i]),
                                col=replace(col,idx,pallet[i]))
                    }
        
                    to_plot <- to_plot %>% mutate(lab=factor(lab,levels=label)) %>% arrange(lab)
                  #plot
                  ggplot(to_plot, aes(x=lab, y=value,fill=cond)) + 
                    geom_bar(stat="identity",position=position_dodge(),color="black",width=0.6) +
                    ggtitle(paste0(gene,", ",title))+
                    scale_fill_manual(values=to_plot$col)+
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

#This plot will expect the first columns to contain significance
#for multiple comparisons the significance it will be in the order below
#example: c(a,b,c) a-b, b-c, a-c
#currently it can only support 3 conditions
#I think that more conditions would require a vector of heights that would adjust as the graph is built...
#The other columns should have mean sem and "name" for the gene
#I could probably do the calculations here, but that's not how I'd designed the tables
barplot_conditions<-function(table=ja2[,c(5,7:11)],
                    title="P0, n=2",
                    conditions=c("con","cko"),
                    pallet=brewer.pal(11,"Spectral"),
                    label=conditions,
                    gene=input$goi)({
  
  #just building an empty table to plot                    
  to_plot<-tibble(cond=conditions,
                  lab=label,
                  pal=pallet[1:length(conditions)],
                  ave=0,
                  sem=0)
  #find gene of interest row
  goi_idx<-which(table$name==gene)
  #fill in empty table
  for(i in conditions){
    #find the name of the appropriate mean column
    a<-grep("mean|ave",grep(i,names(table),value=T),value=T) 
    s<-grep("sem|sd",grep(i,names(table),value=T),value=T)
    to_plot <- to_plot %>%
      #just in case there's more than one number only the first is going to be used ([1])
      mutate(ave=replace(ave,cond==i,as.numeric(table[goi_idx,a])[1]),
             sem=replace(sem,cond==i,as.numeric(table[goi_idx,s])[1]))
  }
  to_plot<-to_plot %>% mutate(lab=factor(lab,levels=label))
  
  #pull out the pvalues. Could be any number of samples
  sig<-as.numeric(table[goi_idx,1:sum(seq(1:length(conditions))-1)])
  
  #base height to adjust for asterisks
  h=max(rowSums(to_plot[,c("ave","sem")]))

  p<-to_plot %>% 
    ggplot(aes(x=lab, y=ave))+
    geom_bar(aes(fill=cond),stat="identity",color="black",width=0.6) +
    scale_fill_manual(values=arrange(to_plot,lab)$pal)+
    geom_errorbar(aes(ymin=ave-sem, ymax=ave+sem), width=.2)+
    ggtitle(paste0(gene,", ",title))+
    theme_bw()+
    labs(y = "FPKM")+
    theme(plot.title = element_text(size = 20),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=12),
          aspect.ratio=0.65,
          legend.position = ('none'),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=16))
  #add first layer of asterisks and adjust h
  for(i in 1:(length(conditions)-1)){
    new_h=0
    if(sig[i]<=0.05){
      new_h=max(rowSums(to_plot[c(i,i+1),c("ave","sem")]))
      p<- p+geom_signif(stat="identity",data=data.frame(x=i+.02, xend=i+.98,y=new_h*1.05,annotation=star(sig[i])),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))
    }
    h<-max(h,new_h*1.1)
  }
  if (length(conditions)>2){
    if (sig[3]<0.05){
      p<- p+geom_signif(stat="identity",data=data.frame(x=1.02, xend=2.98,y=h*1.05,annotation=star(sig[3])),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))
      h=h*1.1
    }
  }
    p+scale_y_continuous(limits=c(0,h))
})

#returns asterisk for significance
star<-function(q=sig[1,4]){
  if(!q<=0.05){return("NS")}
  else if(q<=0.001){return("***")}
  else if(q<=0.01){return("**")}
  else{return("*")}
}

#Volcano plot of fold change and p-value data. Can change the column names as necessary
volcano_gene <- function(table=ctx,x_label="logFC cortex Jarid2 control/cKO",
                         col="#7a9200",p="PValue",fc="logFC",xlim=c(NA,NA),ylim=c(NA,NA),
                         gene=input$goi){
  ggplot(table,aes(-.data[[fc]],-log(.data[[p]])))+
    geom_point(color="grey")+
    geom_point(data=filter(table,name==gene),color=col,size=2)+
    geom_text_repel(data=filter(table,name==gene),aes(label=name))+
    theme_bw()+
    xlab(x_label)+
    scale_x_continuous(limits=xlim, expand = c(0, 0)) +
    scale_y_continuous(limits=ylim, expand = c(0, 0))
}
