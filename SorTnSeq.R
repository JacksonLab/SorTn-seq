library(readxl)
library(writexl)
library(tidyverse)
library(iNEXT)

# Set working directory
  setwd("~/SorTnSeq_analysis/")
  
# Make an output directory and another for plots
  dir.create("output",warnings(F))
  dir.create("plots",warnings(F))

# Variables to set
  three.prime.ignore<-0.1 # Percentage of the 3' end of features to exclude for read counting
  normalisation.scale<-1000000 # Crude reads per million scale

# Automated file loading - all files must be in /plot_files
  file.list<-list.files(path="plot_files",pattern="*.insert_site_plot.gz")  
  file.names<-gsub(".gz","",file.list)
  
# Quick looping for all files - fix later....

for (file.to.analyse in 1:length(file.list)){
  print(paste("Processing sample",file.to.analyse,file.names[file.to.analyse]))
    
# Read the specified file
  data<-read_delim(gzfile(paste0("plot_files/",file.list[file.to.analyse])),delim=" ",col_names=c("plus","minus"))

# Make a new column with the genomic coordinates
  data<-data%>%mutate(position=row_number())

# Combine the strand column data into a single column called "strand" and call the values "reads"
  data.long<-data%>%gather(key ="strand",value="reads","plus","minus")
  
# Rename the reads column based on the input file number
  sample.name<-paste0("reads_S",file.to.analyse)
  data.long<-data.long%>%rename(!!sample.name:=reads) 
  
# Assign every insertion locus (strand inclusive) a unique id
  data.long<-data.long%>%mutate(id=paste0("ins.nt",position,"_",strand))
  
# Only keep the id and reads columns
  data.long<-data.long%>%select(id,!!as.name(sample.name))

# Combine the data from all files into one table by insertion site id
  if(file.to.analyse==1){data.all<-data.long} else {data.all<-data.all%>%full_join(data.long,by="id")}  
    
} # End of the file reading loop
  
# Add nt and strand columns
  data.all<-data.all%>%mutate(nt=gsub("ins.nt","",id))%>%
                        separate(nt,c("nt","strand"),sep="_")%>%
                        select(id,nt,strand,starts_with("reads"))

# Gather the data and then remove empty rows.
  data.clean<-data.all%>%gather(key="sample",value="reads",starts_with("reads_"))%>%filter(reads!=0)

# Also make a human-friendly version...
  final<-data.clean%>%spread("sample","reads",fill=0)

# Write out the tables
  write.csv(data.clean,"output/all_data_clean_reads_table.csv",row.names=F)
  write.csv(final,"output/all_data_wide_reads_table.csv",row.names=F)
   
# Also normalise the read depth per sample - not ideal, but useful for some cases
   data.norm<-data.clean%>%group_by(sample)%>%mutate(norm=reads/sum(reads)*normalisation.scale)%>%ungroup()%>%select(-reads)
   data.norm<-data.norm%>%spread("sample","norm",fill=0)
   colnames(data.norm)<-gsub("reads","norm.reads",colnames(data.norm))
   write.csv(data.norm,"output/all_data_wide_norm_read_table.csv",row.names=F)

#### Now we can run some analyses ####

# First, a rarefaction plot using iNEXT .. very slow...- include prediction of greater read depth (optional).
  plot.data<-final%>%select(starts_with("reads_S")) # Filter here for specific samples
  plot.data<-as.data.frame(plot.data[,1:ncol(plot.data)])
  plot.object<-iNEXT(plot.data,q=0,datatype="abundance",size=NULL,endpoint=NULL,knots=50,se=T,conf=0.95,nboot=1) 
 print(ggiNEXT(plot.object,type=1,se=T,facet.var="none",color.var="site",grey=F) +
    theme_bw(base_size = 10) + xlab("Reads") + ylab("Unique insertions") + ggtitle(file.names[file.to.analyse]))
 ggsave(paste("output/all_data_rarefaction.png",sep=""),dpi=600)

# clear some memory...
 
 rm(data,data.all,data.clean,data.long,final)
 gc()

###### Now we can merge the insertion site data with summary file of all Serratia features..... ######
  
# Read in the data table - start here to save time/resources
  data<-read.csv("output/all_data_wide_reads_table.csv",as.is=T)

# Change the sample names - fix this later...
  colnames(data)<-c("id","nt","strand","S_L-1","S_I-1","S_I-2","S_I-3","S_H-1","S_D-1","S_L-2","S_H-2","S_D-2","S_L-3","S_H-3","S_D-3")

# Read in the LacA feature file....
 features.in<-read_excel("input/serratia_master_features_table.xlsx",sheet="Serratia39006_LacA")
 features.in[is.na(features.in)]<-""
 features.in<-features.in%>%select(-ID) # Remove to avoid clash later
 features.in%>%group_by(type)%>%summarise(counts=n())
  
# Filter the features to only those that we want to analyse
  features<-features.in#%>%filter(type=="gene"|type=="CRISPR")#%>%mutate(locus_tag=gsub(".*;locus_tag=","",attributes))%>%mutate(locus_tag=gsub(";.*","",locus_tag))
 
# Optional (variable set above) shorten the genes by x% to exclude insertions in the 3' ends.... 
  features<-features%>%mutate(end=ifelse(type=="gene" & strand=="+",floor(end-((end-start)*three.prime.ignore)),end))  
  features<-features%>%mutate(start=ifelse(type=="gene" & strand=="-",floor(start+((end-start)*three.prime.ignore)),start))  
  
# Stretch the annotations out along their entire nt positions and number accordingly
  features<-features%>%mutate(length=(end-start+1),feature.length=length)
  expanded<-features%>%uncount(length)
  expanded<-expanded%>%group_by(locus_tag)%>%mutate(nt=start+row_number()-1)
 
# Remove the strand column - this information is in the insertion site id
  sample<-data%>%select(-strand)
  
# Join the insertion site data to the features
  sor.tn.seq<-expanded%>%left_join(sample,by="nt")
  
# Remove the empty rows (no insertion for a particular nt position)
  sor.tn.seq<-sor.tn.seq%>%filter(!is.na(id))

# Gather the data for easier grouping
  sor.tn.seq<-sor.tn.seq%>%gather(key="sample",value="reads",starts_with("S_"))  
  
# Remove entires with zero reads
  sor.tn.seq<-sor.tn.seq%>%filter(reads!=0)
    
# Summarise the insertion information for each locus tag feature - customise this to change the information you get...
  # Note that reads for overlapping features are counted for each feature
  statistics<-sor.tn.seq%>%group_by(sample,locus_tag)%>%summarise(unique.insertions=n(),sum.reads=sum(reads),feature.length=max(feature.length))

# You can also calculate additional statisitics - e.g. the insertion index here...    
  statistics<-statistics%>%mutate(ins.index=unique.insertions/feature.length,nt.per.unq.ins=round(feature.length/unique.insertions,digits=2),reads.per.nt=round(sum.reads/feature.length,digits=5))
    
## Print out seperate tables for each sample

  # First make a list of the samples
    sample.list<-statistics%>%select(sample)%>%unique()
  
  # Loop over all samples
    for(sample.number in 1:nrow(sample.list)){
    
    # Filter for each sample and output the table
      sample<-statistics%>%filter(sample==sample.list[sample.number,])
  
    # Join back to the master features table - you can filter as required later..  
      sample.out<-features%>%left_join(sample)
    
    # Fill in the blank columns with zero
      sample.out[is.na(sample.out)]<-0
  
    # Remove some columns for clarity
      sample.out<-sample.out%>%select(-seqid,-source,-phase,-score,-sample,-length)
    
    # Make a list of all the data frames to become sheets in Excel
      if(sample.number==1){sheet.list<-list(sample.out)} else {sheet.list[[sample.number]]<-sample.out}
    
    } # End of the sample loop  
    
  # Update the list names - these will become the sheet names in Excel
    names(sheet.list)<-gsub("S_","",sample.list$sample)

  # Write the table to a csv file
    write_xlsx(sheet.list,"output/SorTnSeq_all_features_mapped.xlsx")

# Make a table of read counts with all samples - this can be used for EdgeR etc...
    table.reads<-statistics%>%select(locus_tag,sample,sum.reads)
    table.reads<-table.reads%>%spread(key=sample,value=sum.reads,fill=0)
   
  # Write the table to Excel  
    write_xlsx(table.reads,"output/SorTnSeq_table_reads.xlsx")  
    
# Make a table of unique insertion counts with all samples - this can be used for EdgeR etc...
    table.unq.ins<-statistics%>%select(locus_tag,sample,unique.insertions)
    table.unq.ins<-table.unq.ins%>%spread(key=sample,value=unique.insertions,fill=0)

# Write the table to Excel  
    write_xlsx(table.unq.ins,"output/SorTnSeq_table_unique_insertions.xlsx")