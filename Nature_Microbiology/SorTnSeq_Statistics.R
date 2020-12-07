library(edgeR)
library(tidyverse)
library(readxl)
library(writexl)

# Set working directory
setwd("~/SorTnSeq_analysis/")

# Variables to set
  read.cutoff.depleted<-2 # Only consider features that have > x insertions in all samples in the DEPLETED library.

# Read in the data file (table of counts)
  data<-read_excel("output/SorTnSeq_table_unique_insertions.xlsx", sheet = "Sheet1")
  
# Fix column names issue - sort later
  colnames(data)<-gsub("-","_r",colnames(data))
      
# Filter
  data.filtered<-data%>%filter(S_D_r1>=read.cutoff.depleted & S_D_r2>=read.cutoff.depleted & S_D_r3>=read.cutoff.depleted)
  annotations<-data.filtered%>%select(locus_tag)
  reads<-data.filtered%>%select(-locus_tag)

## Here, reads = unique insertions

# List of samples and groups/conditions - these must match the column order exactly.
    samples<-factor(c("Depleted","Depleted","Depleted","High","High","High","Input","Input","Input","Low","Low","Low"))

  # Reorder the factors
    samples<-relevel(samples,ref="Depleted")

# Initiate the edgeR object
    dataset<-DGEList(counts=reads,group=samples,genes=annotations)
  # Check
      dataset

# Look at similarity between samples using multidimensional scaling plots 
  # Standard MDS
    plotMDS.DGEList(dataset,col=rep(1:4,each=3),labels=samples,top=100)
          
  # Based on the biological coefficient of variation (top 100 genes).
    plotMDS.DGEList(dataset, method="bcv",top=100,asp=1,pch=19,col=(rep(1:4,each=3)),cex=2.5,
                xlim=c(-0.75,1),
                xlab="Biological coefficient of variation (distance 1)",
                ylab="Biological coefficient of variation (distance 2)")

## Use the Classic analysis approach

  # Estimate the dispersion - common is required first, then tagwise estimates the dispersion individually for each feature  
    dataset.classic<-estimateCommonDisp(dataset) 
    dataset.classic<-estimateTagwiseDisp(dataset.classic)  
    plotBCV(dataset.classic)
    summary(dataset.classic)
    
## Compare samples to depleted libraries
    results.exact.low<-exactTest(dataset.classic,pair=c("Depleted","Low"))
    results.exact.high<-exactTest(dataset.classic,pair=c("Depleted","High"))
    results.exact.input<-exactTest(dataset.classic,pair=c("Depleted","Input"))
    
  # View the top hits  
    topTags(results.exact.low)  # By default the FDR is the BH qvalue
    topTags(results.exact.high)
    topTags(results.exact.input)
  
  # Plots    
    plotMD(results.exact.low)
    plotMD(results.exact.high)
    plotMD(results.exact.input)

  # Add the data to the annotations table
    results.classic.low<-annotations%>%bind_cols(results.exact.low$table,q.value=p.adjust(results.exact.low$table$PValue,"BH"))
    results.classic.high<-annotations%>%bind_cols(results.exact.high$table,q.value=p.adjust(results.exact.high$table$PValue,"BH"))
    results.classic.input<-annotations%>%bind_cols(results.exact.input$table,q.value=p.adjust(results.exact.input$table$PValue,"BH"))
    
  # Write these data out to excel
    sheet.list<-list(results.classic.low,results.classic.high,results.classic.input)
    names(sheet.list)<-c("Low","High","Input")
    write_xlsx(sheet.list,"output/SorTnSeq_statistics_exact_depleted_unique_insertions.xlsx")
    
# Volcano plots with Padj < 0.05 and log2 fold change = 1
    volcano_plot<-function(to.plot,sample.name){
      plot(to.plot$logFC,-log(to.plot$q.value,base=10),xlim=range(c(-6,6)),xlab=paste("Fold-change (Log2),",sample.name),ylab="p.adjusted (-Log10)",cex=.5,pch=20)
      abline(h=-log10(0.05),col="red")
      abline(v=-1,col="red")
      abline(v=1,col="red")
    }

volcano_plot(results.classic.low,"Low") 
volcano_plot(results.classic.high,"High")    
volcano_plot(results.classic.input,"Input")  

# ggplot for volcano plots # Padj < 0.1 and log FC > 0
  volcano_plot_gg<-function(to.plot,plot.name,point.colour){
    print(ggplot(to.plot,aes(logFC,-log10(q.value))) +
      geom_point(aes(col=q.value<0.1 & logFC >0),show.legend=F,size=0.5) +
      scale_color_manual(values=c("grey 35",point.colour)) +
      geom_hline(yintercept=1,size=0.5) +
      geom_vline(xintercept=0,size=0.5) +
      xlim(-6,4) +
      ylim(-1,35)+
      theme_light()+
      xlab("Fold change from depleted (Log2)") +
      ylab("P. adjusted value (-Log10)") +
      theme(text=element_text(size=6)))
    
  # Output the plot  
    ggsave(paste0("output/volcano_",plot.name,".pdf"),width=4,height=4,units="cm")  
    
  }

# Plot text/sizes are scaled for the pdf output  
  volcano_plot_gg(results.classic.high,"High","red")
  volcano_plot_gg(results.classic.low,"Low","dodgerblue3")
  volcano_plot_gg(results.classic.input,"Input","grey 70")

# Combine with Master features Table

  # Read in the LacA feature file....
    features.in<-read_excel("Input/Serratia_master_features_table.xlsx", sheet = "Serratia39006_LacA")
    features.in[is.na(features.in)]<-""
    features.in%>%group_by(type)%>%summarise(counts=n())
  
  # Select columns that are most useful
    features<-features.in%>%select(locus_tag,type,start,end,strand,gene,product)
  
  # Rename the columns for each sample treatment and add to an output table
    classic.output<-features%>%full_join(results.classic.low%>%rename(logFC_Low=logFC,PValue_Low=PValue,padj_Low=q.value)%>%select(logCPM,everything()),by="locus_tag")
    classic.output<-classic.output%>%full_join(results.classic.high%>%rename(logFC_High=logFC,PValue_High=PValue,padj_High=q.value)%>%select(-logCPM),by="locus_tag")
    classic.output<-classic.output%>%full_join(results.classic.input%>%rename(logFC_Input=logFC,PValue_Input=PValue,padj_Input=q.value)%>%select(-logCPM),by="locus_tag")
    classic.output[is.na(classic.output)]<-""

# Format number columnns so they open correctly in Excel  
    classic.output[8:17]<-sapply(classic.output[8:17],as.numeric)
    
# Output the table  
    write_xlsx(classic.output,"output/SorTnSeq_statistics_exact_depleted_unique_insertions_features.xlsx")
      

  
