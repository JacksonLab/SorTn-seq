library(tidyverse)
library(readxl)
library(writexl)
library(edgeR)
library(scales)
library(ggrepel)
library(ggiraph)

# variables
  genome.prefix<-"GCF_002847015.1_ASM284701v1"
  read.cutoff.depleted<-2 # only consider features that have > x insertions in ALL DEPLETED samples - could be an issue if the depleted pool isn't sampled deep
  bcv.features<-100
  reference.sample<-"depleted"
  # thresholds for plotting
  threshold.fc<-0.5
  threshold.p.adj<-0.05

# read in the counts table (unique insertions)
  data<-read_excel("SorTnSeq_table_unique_insertions.xlsx")
  
# require a features to have a minimum count in ALL depleted samples (see above)
  data.filtered<-data%>%filter_at(vars(starts_with("depleted")),all_vars(. >= read.cutoff.depleted))
  features.pass<-data.filtered%>%select(id)
  
  # check the excluded data
  data.excluded<-data%>%filter_at(vars(starts_with("depleted")),any_vars(. < read.cutoff.depleted))
    
  count.table<-data.filtered%>%select(-seqid,-id)

# sample names are derived from column names  
  samples<-factor(gsub("_.*","",colnames(count.table)))

# initiate an edgeR object
    dataset<-DGEList(counts=count.table,
                     group=samples,
                     genes=features.pass)
  # check
    dataset$samples

# check similarity between samples using the biological coefficient of variation (top n features is set by bcv.features)
    mds<-plotMDS.DGEList(dataset, method="bcv",top=bcv.features,asp=1,pch=19,col=as.numeric(samples),cex=2.5,
                #xlim=c(-0.75,1),
                xlab="Biological coefficient of variation (distance 1)",
                ylab="Biological coefficient of variation (distance 2)")

    # formatted plot
    mds.plot<-full_join(enframe(mds$x)%>%rename(x=value),
                        enframe(mds$y)%>%rename(y=value))%>%
                  separate(name,c("sample.type","replicate"))
    
    ggplot(mds.plot) + geom_point(aes(x,y,color=sample.type)) +
                scale_color_manual(values=c("grey 35","red","grey 70","dodgerblue3")) +
                scale_x_continuous(expand=c(0,0),limits=c(floor(min(mds.plot$x)*6)/6,ceiling(max(mds.plot$x)*6)/6)) +
                scale_y_continuous(expand=c(0,0),limits=c(floor(min(mds.plot$y)*6)/6,ceiling(max(mds.plot$y)*6)/6)) +
                theme_light() +
                geom_text_repel(aes(x,y,label=replicate)) +
                xlab("Biological coefficient of variation (distance 1)") +
                ylab("Biological coefficient of variation (distance 2)") +
        theme(text=element_text(size=8))
    
    ggsave(paste0("SorTnSeq_bcv_plot.pdf"),width=12,height=8,units="cm")  

# estimate the dispersion
    dataset.classic<-estimateCommonDisp(dataset) 
    dataset.classic<-estimateTagwiseDisp(dataset.classic)  
    plotBCV(dataset.classic)
    summary(dataset.classic)

# run the analyses    
  all.results<-tibble()
  for (sample.type in setdiff(samples,reference.sample)) {

    # compare each sample type against the specified reference type
      results<-exactTest(dataset.classic,pair=c(reference.sample,sample.type))
    # general information
      print(topTags(results))
      plotMD(results)

    # add data to annotations, calculate p.adj, update column names, and assign to object (dynamic)      
      tmp<-features.pass%>%bind_cols(results$table)%>%
                         mutate(!!paste0(sample.type,"_p.adj") :=p.adjust(PValue,"BH"))%>%
                         rename(!!paste0(sample.type,"_p.value") :=PValue,
                                !!paste0(sample.type,"_logFC") :=logFC)
  
      # combine all results into one table  
      if (nrow(all.results)==0){all.results<-tmp} else {all.results<-all.results%>%left_join(tmp)}
  }

  all.results<-all.results%>%select(id,logCPM,everything())

  # long results version for plotting
  all.results.long<-all.results%>%gather("key","value",-id,-logCPM)%>%
                    separate(key,c("sample","parameter"),sep="_")%>%
                    spread(parameter,value)

  # add feature information
    features<-read_excel(paste0(genome.prefix,"_features_sortnseq.xlsx"))
  
  all.results.out<-full_join(all.results,features)%>%
                          select(seqid,type,start,end,product,id,contains("high_"),contains("low_"),everything())%>%
                          arrange(seqid,start)
  
  all.results.out%>%write_xlsx(paste0("SorTnSeq_results_",reference.sample,"_unique_insertions.xlsx"))   
  
# volcano plots # Padj < 0.05 and log FC > 0.5

  volcano_plot_gg<-function(plot.sample,point.colour){
  
    print(ggplot(all.results.long%>%filter(sample==plot.sample),aes(logFC,-log10(p.adj))) +
      geom_point(aes(col=p.adj<threshold.p.adj & logFC >threshold.fc),show.legend=F,size=1.5) +
      scale_color_manual(values=c("grey 35",point.colour)) +
      geom_hline(yintercept=-log10(threshold.p.adj),size=0.5) +
      geom_vline(xintercept=threshold.fc,size=0.5) +
      # use the same scale limits for all plots
      scale_x_continuous(expand=c(0,0),limits=c(floor(min(all.results.long$logFC)),ceiling(max(all.results.long$logFC))+0.5)) +
      scale_y_continuous(expand=c(0,0),limits=c(0,ceiling(-log10(min(all.results.long$p.adj))+0.5))) +  
      theme_light()+
      xlab(paste0("Fold change from ",reference.sample," (Log2)")) +
      ylab("P. adjusted value (-Log10)") +
      theme(text=element_text(size=8)))
    
    ggsave(paste0("volcano_",reference.sample,"-",plot.sample,".pdf"),width=8,height=8,units="cm")  
    
  }  
  
volcano_plot_gg("high","red")    
volcano_plot_gg("low","dodgerblue3")    
volcano_plot_gg("input","grey 70")    
  
# plot the positive fold-change for low and high as a SorTnSeq enrichment plot (interactive)

# scale the point size based on the insertion index
  insertion.index.mean<-read_excel("SorTnSeq_table_insertion_index.xlsx")%>%
                      gather("sample","insertion.index",-seqid,-id)%>%
                      separate(sample,c("sample","replicate"),sep="_")%>%
                      group_by(seqid,id,sample)%>%
                      summarise(ins.per.kb.enriched=mean(insertion.index)*1000)
                      
  
  enrichment<-all.results.long%>%filter(sample %in% c("high","low") &
                                          logFC >= 0)%>%
    mutate(enrichment=ifelse(sample=="high",logFC,logFC),
           color=ifelse(sample=="high","high","low"))%>%
    mutate(color=ifelse(p.adj<threshold.p.adj & logFC >threshold.fc, color, "depleted"))
  
  to.plot<-enrichment%>%left_join(insertion.index.mean)%>%left_join(features)%>%
    mutate(display.name=ifelse(is.na(name),id,paste0(name,": ",product)))
  
  plot<-ggplot(to.plot,aes(x=enrichment,y=-log10(p.adj),color=color,size=ins.per.kb.enriched)) +
    geom_point_interactive(aes(tooltip=display.name)) +
    scale_color_manual(values=c("grey 35","red","dodgerblue3"), labels=c("not significant", "high expression", "low expression")) +
    geom_hline(yintercept=-log10(threshold.p.adj),size=0.5) +
    geom_vline(xintercept=0.5,size=0.5) +
    theme_light()+
    xlab(paste0("SorTnSeq enrichment relative to the ",reference.sample," sample (Log2)")) +
    ylab("P. adjusted value (-Log10)") +
    theme(text=element_text(size=8)) +
    guides(color=guide_legend("enriched features"))+
    facet_wrap("sample")
    

print(plot)
ggsave(paste0("SorTnSeq_enrichment_",reference.sample,".pdf"),width=12,height=8,units="cm")  

girafe(code={print(plot)})

  