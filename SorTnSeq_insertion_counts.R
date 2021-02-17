library(tidyverse)
library(readxl)
library(writexl)

# variables
  genome.prefix<-"GCF_002847015.1_ASM284701v1" # Please update for your target genome
  trim.3.prime<-0.1 # Proportion of the 3' end of features (non-intergenic) to exclude for read counting
  trim.5.prime<-0.1 # Proportion of the 5' end of features (non-intergenic) to exclude for read counting

# read in the sample data
  sample.metadata<-read_xlsx("sample_metadata.xlsx")%>%
                            mutate(sample.name=paste(sample.type,replicate,sep="_"))

for (sample in 1:nrow(sample.metadata)){
  print(paste0("Processing sample #",sample,": ",sample.metadata[sample,1]))
    
# Read the specified file
  data<-read_delim(Sys.glob(file.path(paste0("bam/",sample.metadata$plot.file.prefix[sample],"*",".mapped.bam.bed"))),
                    delim="\t", 
                    col_names=c("seqid",
                                "start",
                                "end",
                                "read.id",
                                "map.quality",
                                "strand"),
                    col_types=cols(seqid = col_character(),
                                   start = col_double(),
                                   end = col_double(),
                                   read.id = col_character(),
                                   map.quality = col_double(),
                                   strand = col_character())) 
  
  insertion.sites<-data%>%mutate(nt=ifelse(strand=="+",start+1,end))%>%
                          group_by(seqid,nt,strand)%>%
                          summarise(reads=n())%>%
                          mutate(strand=ifelse(strand=="+","plus","minus"),
                                 sample=sample.metadata$sample.name[sample])%>%
                          mutate(tn.id=paste0("ins-",seqid,"-nt",sprintf("%08d",nt),"_",strand))

  # combine the data from all files into one table
  if(sample==1){data.all<-insertion.sites} else {data.all<-data.all%>%bind_rows(insertion.sites)}  
    
} 
  
# format
  data.wide<-data.all%>%spread("sample","reads",fill=0)%>%
                        arrange(seqid,nt)
  
# merge the insertion site data with genome features
  
# read in the feature file
 features<-read_excel(paste0(genome.prefix,"_features_sortnseq.xlsx"))

# optional (variable set above) shorten the features (excluding intergenic) to exclude insertions in the 3' and 5' ends 
  features<-features%>%mutate(length=end-start+1)%>%
                      mutate(end=  case_when(type!="intergenic" & strand=="+" ~ floor(end-length*trim.3.prime),
                                             type!="intergenic" & strand=="-" ~ floor(end-length*trim.5.prime),
                                             T ~ end),
                             start=case_when(type!="intergenic" & strand=="-" ~ floor(start+length*trim.3.prime),
                                             type!="intergenic" & strand=="+" ~ floor(start+length*trim.5.prime),   
                                             T ~ start))%>%
                      # update length
                      mutate(length=end-start+1)
  
# stretch the annotations out along their entire nt positions and number accordingly
  expanded<-features%>%mutate(feature.length=length)%>%
                        uncount(length)%>%
                        group_by(seqid,id)%>%mutate(nt=start+row_number()-1)
  
# remove the strand column - clashes with the feature strands and is included in the insertion site id (tn.id)
  sample<-data.wide%>%select(-strand)
  
# join the insertion site data to the features and format
  sor.tn.seq<-expanded%>%left_join(sample,by=c("seqid","nt"))%>%
                filter(!is.na(tn.id))%>%
                gather(key="sample",value="reads",any_of(sample.metadata$sample.name))%>%
                filter(reads!=0)%>%
                select(seqid,tn.id,sample,reads,nt,everything())%>%
                # reads for overlapping features need to be split
                group_by(seqid,sample,tn.id)%>%
                mutate(overlapping.features=n())%>%
                mutate(reads.adj=reads/overlapping.features,
                       insertions.adj=1/overlapping.features)
  
# summarise the insertion information for each locus tag feature
  insertion.table<-sor.tn.seq%>%group_by(sample,seqid,id,feature.length)%>%
                                summarise(unique.insertions=sum(insertions.adj),
                                          reads=sum(reads.adj))
  
# calculate the insertion index (adjust the feature length by the excluded 3 and 5 prime ends)
  insertion.table<-insertion.table%>%mutate(ins.index=unique.insertions/(feature.length-feature.length*(trim.3.prime+trim.5.prime)))%>%
                                  select(-feature.length)

# table of all expected features
  all.features<-features%>%select(seqid,id)
  
  # make count tables for EdgeR analysis - user specified metric (e.g. unique insertions, reads, insertion index)
  format_count_table<-function(metric){
    insertion.table%>%select(seqid,id,sample,!!metric)%>%
                    spread(key=sample,value=!!metric,fill=0)%>%
                    full_join(all.features)%>%
                    replace(is.na(.),0)%>%
                    arrange(seqid,id)}
  
  format_count_table("unique.insertions")%>%write_xlsx("SorTnSeq_table_unique_insertions.xlsx") 
  format_count_table("ins.index")%>%write_xlsx("SorTnSeq_table_insertion_index.xlsx")   
  format_count_table("reads")%>%write_xlsx("SorTnSeq_table_reads.xlsx")   
  
# feature information for each sample 
  
  for(sample.number in 1:nrow(sample.metadata)){
    
    # filter for each sample and output the table
      sample<-insertion.table%>%ungroup()%>%
                                filter(sample==sample.metadata$sample.name[sample.number])%>%
                                select(-sample)
  
    # join to the feature information 
      sample.out<-features%>%replace(is.na(.),"")%>%
                             left_join(sample)%>%
                             replace(is.na(.),0)%>%
                             select(seqid,id,unique.insertions,ins.index,reads,type,start,end,strand,everything())%>%
                             arrange(seqid,start)
    
    # list data frames to become sheets in Excel
      if(sample.number==1){sheet.list<-list(sample.out)} else {sheet.list[[sample.number]]<-sample.out}
    }
    
  # update sheet names and write Excel file
    names(sheet.list)<-sample.metadata$sample.name
    write_xlsx(sheet.list,"SorTnSeq_all_features_by_sample.xlsx")    
    
