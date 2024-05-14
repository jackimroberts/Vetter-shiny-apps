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
                  label=c("control","cKO"),
                  pallet=brewer.pal(11,"Spectral"),
                  gene=input$goi)({
                    
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
                  
                  #plot
                  ggplot(to_plot, aes(x=lab, y=value,fill=cond)) + 
                    geom_bar(stat="identity",position=position_dodge(),color="black",width=0.6) +
                    ggtitle(paste0(gene,", ",title))+
                    scale_fill_manual(values=arrange(to_plot,lab)$col)+
                    theme_bw()+
                    labs(y = "ave")+
                    theme(axis.text.x = element_text(size=16),
                          plot.title = element_text(size = 20),
                          axis.text.y = element_text(size=12),
                          aspect.ratio=0.65,
                          legend.position = ('none'),
                          axis.title.x=element_blank(),
                          axis.title.y=element_text(size=16))
})

#This plot will expect the first columns to contain significance
#for multiple comparisons it will be in the order of the provided conditions
#example: c(a,b,c) a-b, a-c, b-c
#The other colums should have mean sem and "name" for the gene
#I could probably do the calculations here, but that's not how I'd designed the tables
barplot_conditions<-function(table=ja2[,c(5,7:11)],
                    title="P0, n=2",
                    conditions=c("con","cko"),
                    label=c("control","cKO"),
                    pallet=brewer.pal(11,"Spectral"),
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
    a<-grep("mean|ave",grep(i,names(file),value=T),value=T) 
    s<-grep("sem|sd",grep(i,names(file),value=T),value=T)
    to_plot <- to_plot %>%
      mutate(ave=replace(ave,cond==i,as.numeric(table[goi_idx,a])),
             sem=replace(sem,cond==i,as.numeric(table[goi_idx,s])))
  }
  #pull out the pvalues. Could be any number of samples
  sig<-as.numeric(table[goi_idx,1:sum(seq(1:length(conditions))-1)])
  
  #base height to adjust for asterisks
  h=max(rowSums(to_plot[,c("ave","sem")]))
  
  to_plot %>% 
    mutate(lab=factor(lab,levels=label)) %>%
    ggplot(aes(x=lab, y=ave))+
    geom_bar(aes(fill=cond),stat="identity",color="black",width=0.6) +
    scale_fill_manual(values=arrange(to_plot,lab)$pal)+
    geom_errorbar(aes(ymin=ave-sem, ymax=ave+sem), width=.2)+
    ggtitle(paste0(gene,", ",title))+
    theme_bw()+
    theme(plot.title = element_text(size = 20),
          axis.text.x = element_text(size=16),
          axis.text.y = element_text(size=12),
          aspect.ratio=0.65,
          legend.position = ('none'),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=16))+
    {if(sig[1]<=0.05)scale_y_continuous(limits=c(0,h*1.1))}+
    {if(sig[1]<=0.05)geom_signif(stat="identity",data=data.frame(x=1, xend=2,y=h*1.05,annotation=star(sig[1])),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))}
})

#returns asterisk for significance
star<-function(q=mytable[1,4]){
  if(!q<=0.05){return("NS")}
  else if(q<=0.001){return("***")}
  else if(q<=0.01){return("**")}
  else{return("*")}
}

#Volcano plot of fold change and p-value data. Can change the column names as necessary
volcano_gene <- function(table=ctx,x_label="logFC cortex Jarid2 control/cKO",col="#7a9200",p="PValue",fc="logFC",gene_name=input$goi){
  ggplot(table,aes(-.data[[fc]],-log(.data[[p]])))+
    geom_point(color="grey")+
    geom_point(data=filter(table,gene==gene_name),color=col)+
    geom_text_repel(data=filter(table,gene==gene_name),aes(label=gene))+
    theme_bw()+
    xlab(x_label)
}
