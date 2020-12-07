library(tidyverse)
library(readxl)
library(writexl)

generate_feature_table<-function(file.prefix){

# file.prefix<-"GCF_002847015.1_ASM284701v1"  
gff.refseq<-read_delim(paste0(file.prefix,"_genomic.gff"),delim="\t",comment="#",
                       col_names=c("seqid",
                                   "source",
                                   "type",
                                   "start",
                                   "end",
                                   "score",
                                   "strand",
                                   "phase",
                                   "attributes"),
                       col_types=cols(seqid = col_character(),
                                      source = col_character(),
                                      type = col_character(),
                                      start = col_double(),
                                      end = col_double(),
                                      score = col_character(),
                                      strand = col_character(),
                                      phase = col_character(),
                                      attributes = col_character()))

features.in<-gff.refseq%>%
                  # remove contig/chromosome rows
                  filter(type!="region")%>%
                  # convert attributes to columns
                  separate_rows(attributes,sep=";")%>%
                  separate(attributes,c("attribute","value"),sep="=")%>%
                  mutate(attribute=tolower(attribute))%>%
                  spread(attribute,value,fill="")

  # deal with any CRISPR arrays or repeats

  if ("rpt_family" %in% colnames(features.in)) {
      features.in<-features.in%>%
                  # merge CRISPR columns...
                  unite("crispr",rpt_family,rpt_type,rpt_unit_range,rpt_unit_seq,sep="_")%>%
                  mutate(crispr=sub("___","",crispr),
                         partial=sub("true","partial",partial),
                         pseudo=sub("true","pseudo",pseudo))%>%
                  unite("other",crispr,regulatory_class,partial,pseudo,exception,sep="")%>%
                  # assign the other column to empty locus tags
                  mutate(locus_tag=ifelse(locus_tag=="",other,locus_tag))
  
  } else {features.in<-features.in%>%mutate(crispr="",other="")}

# split the parent and children features, then remove empty columns
  features.parent<-features.in%>%filter(parent=="")%>%
                                  select_if(~!(all(.=="")))
  
  features.children<-features.in%>%filter(parent!="")%>%
                                    select_if(~!(all(.=="")))

  # find sub-children (ignore)
  #features.siblings<-features.children%>%filter(parent %in% features.children$id)%>%
  #                                        select(-protein_id,-name,-other)
  
  features.children<-features.children%>%filter(!parent %in% features.children$id)
  
  # tidy children data  
  features.children<-features.children%>%mutate(id.child=id,
                                 id=parent)%>%
                            rename(type.child=type,
                                   source.child=source,
                                   dbxref.child=dbxref,
                                   gbkey.child=gbkey,
                                   inference.child=inference,
                                   note.child=note)%>%
                            select(-score,-strand)%>%
                            # frameshifted genes have two entries... need just one
                            group_by(id)%>%
                            mutate(start=min(start),
                                   end=max(end),
                                   phase=0)%>%
                            distinct()%>%
                            ungroup()%>%
                            # select parent, product and columns with .child
                            select(id,product,ends_with(".child"))
  
  # merge children/parent data
  features.detailed<-left_join(features.parent,features.children)%>%
                                select(seqid,source,type,id,start,end,strand,locus_tag,product,everything())

  # update the type field
  features.detailed<-features.detailed%>%mutate(type=case_when(grepl("CRISPR",locus_tag) ~ "CRISPR",
                                                                                   T ~ type ))%>%
                                        filter(type!="sequence_feature")
    
  # find all the intergenic regions (any nt position with no features present)
  
    # first, we need to know how long each contig (seqid) is
    contig.lengths<-gff.refseq%>%filter(type=="region")%>%
                                  group_by(seqid)%>%
                                  summarise(length=end-start+1)

    # make a table with a row for each nt
    genome.nt<-contig.lengths%>%group_by(seqid)%>%uncount(length)%>%
                                          mutate(nt=row_number())
    
  # annotate intergenic regions
    
      # convert features to long format then keep a list of nt positions with features
      feature.coverage<-features.detailed%>%
                    select(seqid,id,start,end,strand)%>%
                    group_by(id)%>%
                    mutate(length=end-start+1)%>%
                    uncount(length)%>%
                    mutate(nt=start+row_number()-1)%>%
                    ungroup()%>%
                    select(seqid,nt)%>%
                    distinct()
  
      # match features against genome
      intergenic<-genome.nt%>%anti_join(feature.coverage,by=c("seqid","nt"))%>%
                              group_by(seqid)%>%
                              # number each intergenic region
                              mutate(intergenic.region=cumsum(c(1,abs(diff(nt))>1)))%>%
                              group_by(seqid,intergenic.region)%>%
                              mutate(start=min(nt),
                                     end=max(nt))%>%
                              select(-nt)%>%
                              distinct()%>%
                              #mutate(id=paste0("intergenic_",sprintf("%04d",intergenic.region)),
                              mutate(id=paste0("intergenic_",seqid,"-nt",sprintf("%08d",start)),            
                                      type="intergenic")%>%
                              ungroup()%>%
                              select(-intergenic.region)%>%
                              mutate(source="SorTnSeq",
                                     score=".",
                                     strand="+",
                                     phase=".",
                                     name=id)

    # merge intergenic and features, then format
    features.all<-bind_rows(features.detailed,intergenic)%>%
                        select(seqid,source,type,start,end,score,strand,phase,id,name,product,everything())%>%
                        arrange(seqid,start)
    
  # write a GFF file to check feature assignments
   features.gff<-features.all%>%gather(key="attribute",value="value",9:ncol(features.all))%>%
                      filter(!is.na(value))%>%
                      filter(value!="")%>%
                      unite("attribute",attribute,value,sep="=")%>%
                      group_by(seqid,source,type,start,end,score,strand,phase)%>%
                      mutate(attributes=paste(attribute,collapse=";"))%>%
                      ungroup()%>%
                      select(-attribute)%>%
                      distinct()%>%
                      mutate(phase=replace_na(phase,0))
    
  features.gff%>%write_delim(paste0(file.prefix,"_features_sortnseq.gff"),delim="\t",col_names=F)  

  # write the SorTnSeq format feature table
    
  features.sortnseq<-features.all%>%select(-any_of(c("phase",
                                        "score",
                                        "partial",
                                        "pseudo",
                                        "end_range",
                                        "start_range",
                                        "inference",
                                        "regulatory_class",
                                        "gene_biotype",
                                        "id.child")))
  
  features.sortnseq[is.na(features.sortnseq)]<-""

  features.sortnseq%>%write_xlsx(paste0(file.prefix,"_features_sortnseq.xlsx"))
}

generate_feature_table("GCF_002847015.1_ASM284701v1") # serratia




