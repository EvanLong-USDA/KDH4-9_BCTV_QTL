setwd("C:/Users/evan.long/USDA/REE-ARS-PWA-2054-SBGradedPool - General/BCTV_GPSeq_9-27/Analysis_result/Summary/Summary/4_GradedPool-Seq")
library(data.table)
library(dplyr)
library(ggplot2)


windowsizes <- c(400,800)
Sig_Thresholds <- c(8,12,16)
slidingType <- c("non-overlap")

SigTable <- cbind.data.frame(WindowSize=NA,SigValue=NA,Slide=NA,max_snp_ratio=NA,Mean_snp_ratio=NA,NumCandidateGenes=NA)

for(wind in windowsizes){
  for(sig in Sig_Thresholds){
    for(slide in slidingType){
    print(paste(wind,sig,slide))
     step <- ifelse(slide=="overlap",wind/2,0)
     Genes <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/Gene_in_Candidate_region.txt"))
     ncand <- nrow(Genes)
     Snpratio <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/SNP.ratio.txt"))
     mean_snp_ratio <- mean(Snpratio$Ratio)
     max_snp_ratio <- max(Snpratio$Ratio)
     SigTable <- rbind.data.frame(SigTable,c(wind,sig,slide,max_snp_ratio,mean_snp_ratio,ncand))
          }
  }
}

SigTable <- na.omit(SigTable)
SigTable$Mean_snp_ratio <- as.numeric(SigTable$Mean_snp_ratio)
SigTable$max_snp_ratio <- as.numeric(SigTable$max_snp_ratio)
SigTable$NumCandidateGenes <- as.numeric(SigTable$NumCandidateGenes)

library(RColorBrewer)

ggplot(SigTable, aes(x=max_snp_ratio, y=NumCandidateGenes)) +
  geom_point(aes(col=WindowSize, size=SigValue, shape=Slide)) +
  scale_color_brewer(palette="Set1") +
  labs(title="Relationship Between SNP Ratio and Candidate Genes",
       x="Max SNP Ratio",
       y="Number of Candidate Genes") +
  theme_minimal(base_size = 15) +
  theme(panel.grid.major = element_line(color = "grey80"), 
        panel.grid.minor = element_blank()) +
  guides(col=guide_legend(title="Window Size"),
         size=guide_legend(title="Significance Value"),
         shape=guide_legend(title="Slide Type"))


####Look at the Significant windows and see how many how many candidate genes are also differentially expressed

DE_Genes_table <-  fread("DEG/R_PVSS_P_Transcript_differential_expression.csv")
DE_Genes_table <- DE_Genes_table %>% filter(significant=="yes")
DE_Genes <- DE_Genes_table$gene_name

DE_SigTable <- cbind.data.frame(WindowSize=NA,SigValue=NA,Slide=NA,max_snp_ratio=NA,Mean_snp_ratio=NA,NumCandidateGenes=NA,CandidateGenes_DE=NA)

for(wind in windowsizes){
  for(sig in Sig_Thresholds){
    for(slide in slidingType){
      # print(paste(wind,sig,slide))
      step <- ifelse(slide=="overlap",wind/2,0)
      Genes <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/Gene_in_Candidate_region.txt"))
      Genes$`#GeneID` <- gsub("gene-","",Genes$`#GeneID`)
      ncand <- nrow(Genes)
      Snpratio <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/SNP.ratio.txt"))
      mean_snp_ratio <- mean(Snpratio$Ratio)
      max_snp_ratio <- max(Snpratio$Ratio)
      NumsigDE_in_Candidates <- sum(as.character(DE_Genes) %in% as.character(unlist(Genes[,1])))
      sigDE_in_Candidates <- DE_Genes[which(as.character(DE_Genes) %in% as.character(unlist(Genes[,1])))]
      if(length(sigDE_in_Candidates)>=1){
      DE_SigTable <- rbind.data.frame(DE_SigTable,cbind.data.frame(WindowSize=wind,SigValue=sig,Slide=slide,max_snp_ratio=max_snp_ratio,Mean_snp_ratio=mean_snp_ratio,NumCandidateGenes=ncand,CandidateGenes_DE=sigDE_in_Candidates))
      }
    }
  }
}




DE_SigTable <- na.omit(DE_SigTable)
DE_SigTable$Mean_snp_ratio <- as.numeric(DE_SigTable$Mean_snp_ratio)
DE_SigTable$max_snp_ratio <- as.numeric(DE_SigTable$max_snp_ratio)
DE_SigTable$NumCandidateGenes <- as.numeric(DE_SigTable$NumCandidateGenes)
# DE_SigTable$NumCandidateGenes_DE <- as.numeric(DE_SigTable$NumCandidateGenes_DE)

colnames(DE_Genes_table)[11] <- "CandidateGenes_DE"
DE_Genes_table <- DE_Genes_table[,c(11,12,13)]
DE_SigTable <- inner_join(DE_SigTable,DE_Genes_table)

# write.table(DE_SigTable,"DEG_in_GPwindows.txt",quote=F,row.names=F,col.names=T,sep="\t")

#######
####### More stringent method of grabbing significant windows (Only grab windows that fall in top 3-5 snp ratio)
#######

get_genes_in_windows <- function(genes,snpratio_wind){
  filtergenes <- data.frame(matrix(nrow = 0,ncol=ncol(genes)))
  colnames(filtergenes) <- colnames(genes)
  if(nrow(genes) >=1 & nrow(snpratio_wind)>=1){
    for(x in 1:nrow(genes)){
      for(y in 1:nrow(snpratio_wind)){
        if(genes[x,2]==snpratio_wind[y,1] && genes[x,3] > snpratio_wind[y,4] && genes[x,4] < snpratio_wind[y,5]){filtergenes <- rbind.data.frame(filtergenes,genes[x,])}
      }
    }
  }
  return(filtergenes)
}


DEG_comparisons <- c("H1VSM1","H1VSM1VSM2VSS","M1VSM2","M2VSS","H1VSS","R_PVSS_P")
DE_Method <- "H1VSM1"
for(DE_Method in DEG_comparisons){

          DE_Genes_table <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv"))
          DE_Genes_table <- DE_Genes_table %>% filter(significant=="yes")
          DE_Genes <- DE_Genes_table$gene_name
          
          windowsizes <- c(400,800)
          Sig_Thresholds <- c(8,12,16)
          slidingType <- c("non-overlap")
                           
          DE_SigTable <- cbind.data.frame(WindowSize=NA,SigValue=NA,Slide=NA,max_snp_ratio=NA,Mean_snp_ratio=NA,NumCandidateGenes=NA,CandidateGenes_DE=NA)
          for(wind in windowsizes){
              for(sig in Sig_Thresholds){
                  for(slide in slidingType){
                    # print(paste(wind,sig,slide))
                    step <- ifelse(slide=="overlap",wind/2,0)
                    Genes <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/Gene_in_Candidate_region.txt"))
                    Snpratio <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/SNP.ratio.txt"))
                    mean_snp_ratio <- mean(Snpratio$Ratio)
                    max_snp_ratio <- max(Snpratio$Ratio)
                    #Get genes only in top 3-5 windows
                    Snpratio <- Snpratio %>% dplyr::arrange(Ratio) %>% filter(Ratio>0.02) %>% dplyr::slice_tail(n=5)
                    Candidate_Genes <-  get_genes_in_windows(Genes,Snpratio)
                    ncand <- nrow(Candidate_Genes)
                    NumsigDE_in_Candidates <- sum(as.character(DE_Genes) %in% as.character(unlist(Candidate_Genes[,1])))
                    sigDE_in_Candidates <- DE_Genes[which(as.character(DE_Genes) %in% as.character(unlist(Candidate_Genes[,1])))]
                      if(length(sigDE_in_Candidates)>=1){
                        DE_SigTable <- rbind.data.frame(DE_SigTable,cbind.data.frame(WindowSize=wind,SigValue=sig,Slide=slide,max_snp_ratio=max_snp_ratio,Mean_snp_ratio=mean_snp_ratio,NumCandidateGenes=ncand,CandidateGenes_DE=sigDE_in_Candidates))
                      }
                  }
              }
          }
          
          
          
          
          DE_SigTable <- na.omit(DE_SigTable)
          DE_SigTable$Mean_snp_ratio <- as.numeric(DE_SigTable$Mean_snp_ratio)
          DE_SigTable$max_snp_ratio <- as.numeric(DE_SigTable$max_snp_ratio)
          DE_SigTable$NumCandidateGenes <- as.numeric(DE_SigTable$NumCandidateGenes)
          # DE_SigTable$NumCandidateGenes_DE <- as.numeric(DE_SigTable$NumCandidateGenes_DE)
          
          colnames(DE_Genes_table)[10] <- "CandidateGenes_DE"
          DE_Genes_table <- DE_Genes_table[,c(10,12,13)]
          DE_SigTable <- inner_join(DE_SigTable,DE_Genes_table)
          
          write.table(DE_SigTable,paste0(DE_Method,"_DEG_in_Top5GPwindows.txt"),quote=F,row.names=F,col.names=T,sep="\t")
          DE_SigTable <- DE_SigTable %>% select(c(CandidateGenes_DE,GO,KEGG)) %>% distinct()
      
          write.table(DE_SigTable,paste0(DE_Method,"_DEG_in_Top5GPwindows_UniqueGenes.txt"),quote=F,row.names=F,col.names=T,sep="\t")

}


####
#### Let's make a heatmap of how similar each DEG comparison is
####

DEG_comparisons <- c("H1VSM1","H1VSM1VSM2VSS","M1VSM2","M2VSS","H1VSS","R_PVSS_P")

DEG_matrix <- data.frame(matrix(ncol=1,nrow=0))
colnames(DEG_matrix) <- c("Gene")
for(DE_Method in DEG_comparisons){
  DE_Genes_table <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv"))
  if(nrow(DEG_matrix)==0){   
    DE_Genes_table <- DE_Genes_table %>% select(gene_id,significant) 
    colnames(DE_Genes_table) <- c("Gene",DE_Method)
    DEG_matrix<- DE_Genes_table
  }else{
    DE_Genes_table <-  DE_Genes_table %>% select(gene_id,significant) 
    colnames(DE_Genes_table) <- c("Gene",DE_Method)
    DEG_matrix <- inner_join(DEG_matrix,DE_Genes_table) %>% distinct()
  }


}


DEG_matrix_num <- DEG_matrix

DEG_matrix_num[,2] <- ifelse(DEG_matrix_num[,2]=="yes",1,0)
DEG_matrix_num[,3] <- ifelse(DEG_matrix_num[,3]=="yes",1,0)
DEG_matrix_num[,4] <- ifelse(DEG_matrix_num[,4]=="yes",1,0)
DEG_matrix_num[,5] <- ifelse(DEG_matrix_num[,5]=="yes",1,0)
DEG_matrix_num[,6] <- ifelse(DEG_matrix_num[,6]=="yes",1,0)
DEG_matrix_num[,7] <- ifelse(DEG_matrix_num[,7]=="yes",1,0)

#rownames(DEG_matrix_num) <- DEG_matrix_num[,1]
DEG_matrix_num <- DEG_matrix_num[,2:7]


library(pheatmap)
# Create the heatmap
DEG_matrix_cor <- cor(DEG_matrix_num)
heatmap(DEG_matrix_cor,Colv=NA,col=heat.colors(256))


# plot_pheatmap <- pheatmap(DEG_matrix_cor,
#          main = "Heatmap",
#          clustering_distance_rows = "correlation",
#          clustering_distance_cols = "correlation",
#          color= colorRampPalette(c("lightblue","white","red"))(50),
#          cluster_rows = TRUE,
#          cluster_cols = TRUE,
#          display_numbers = TRUE,
#          fontsize = 22)
# print(plot_pheatmap)


Difference_Bars <- t(colSums(DEG_matrix_num))
Difference_Bars <- t(Difference_Bars)
colnames(Difference_Bars) <- "N_DEG"
Difference_Bars <- cbind.data.frame(Comparison=rownames(Difference_Bars),Difference_Bars)
Difference_Bars$Comparison <- as.factor(Difference_Bars$Comparison)
ggplot(Difference_Bars,aes(x=Comparison,y=N_DEG))+
         geom_col(aes(fill=Comparison))+
  theme_minimal()

Most_common_DEG <- rowSums(DEG_matrix_num[,1:5])
Most_common_DEG <- cbind.data.frame(Gene=DEG_matrix[,1],DEG_Frequency=Most_common_DEG)
ggplot(Most_common_DEG,aes(x=DEG_Frequency))+
  geom_histogram(aes(),binwidth = 1,fill="lightblue")

sum(Most_common_DEG$DEG_Frequency==5)

Always_sig_Genes <- Most_common_DEG$Gene[Most_common_DEG$DEG_Frequency==5]

#####
##### Now look only at ~200 genes DE in all comparisons (Except Rmax vs R)
#####

get_genes_in_windows <- function(genes,snpratio_wind){
  filtergenes <- data.frame(matrix(nrow = 0,ncol=ncol(genes)))
  colnames(filtergenes) <- colnames(genes)
  if(nrow(genes) >=1 & nrow(snpratio_wind)>=1){
    for(x in 1:nrow(genes)){
      for(y in 1:nrow(snpratio_wind)){
        if(genes[x,2]==snpratio_wind[y,1] && genes[x,3] > snpratio_wind[y,4] && genes[x,4] < snpratio_wind[y,5]){filtergenes <- rbind.data.frame(filtergenes,genes[x,])}
      }
    }
  }
  return(filtergenes)
}


windowsizes <- c(400,800)
Sig_Thresholds <- c(8,12,16)
slidingType <- c("non-overlap")
DE_SigTable <- cbind.data.frame(WindowSize=NA,SigValue=NA,Slide=NA,max_snp_ratio=NA,Mean_snp_ratio=NA,NumCandidateGenes=NA,CandidateGenes_DE=NA)
  for(wind in windowsizes){
    for(sig in Sig_Thresholds){
      for(slide in slidingType){
        # print(paste(wind,sig,slide))
        step <- ifelse(slide=="overlap",wind/2,0)
        Genes <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/Gene_in_Candidate_region.txt"))
        Genes$`#GeneID` <- gsub("gene-","",Genes$`#GeneID`)
        Snpratio <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/SNP.ratio.txt"))
        mean_snp_ratio <- mean(Snpratio$Ratio)
        max_snp_ratio <- max(Snpratio$Ratio)
        #Get genes only in top 3-5 windows
        Snpratio <- Snpratio %>% dplyr::arrange(Ratio) %>% filter(Ratio>0.02) %>% dplyr::slice_tail(n=15)
        Candidate_Genes <-  get_genes_in_windows(Genes,Snpratio)
        ncand <- nrow(Candidate_Genes)
        NumsigDE_in_Candidates <- sum(as.character(DE_Genes) %in% as.character(unlist(Candidate_Genes[,1])))
        sigDE_in_Candidates <- DE_Genes[which(as.character(DE_Genes) %in% as.character(unlist(Candidate_Genes[,1])))]
        if(length(sigDE_in_Candidates)>=1){
          DE_SigTable <- rbind.data.frame(DE_SigTable,cbind.data.frame(WindowSize=wind,SigValue=sig,Slide=slide,max_snp_ratio=max_snp_ratio,Mean_snp_ratio=mean_snp_ratio,NumCandidateGenes=ncand,CandidateGenes_DE=sigDE_in_Candidates))
        }
      }
    }
  }
  
  DE_SigTable <- na.omit(DE_SigTable)
  DE_SigTable$Mean_snp_ratio <- as.numeric(DE_SigTable$Mean_snp_ratio)
  DE_SigTable$max_snp_ratio <- as.numeric(DE_SigTable$max_snp_ratio)
  DE_SigTable$NumCandidateGenes <- as.numeric(DE_SigTable$NumCandidateGenes)
  # DE_SigTable$NumCandidateGenes_DE <- as.numeric(DE_SigTable$NumCandidateGenes_DE)
  
  DE_Genes_table <-  fread("DEG/H1VSS_Transcript_differential_expression.csv")
  colnames(DE_Genes_table)[11] <- "CandidateGenes_DE"
  DE_Genes_table <- DE_Genes_table[,c(11,12,13)]
  DE_SigTable <- inner_join(DE_SigTable,DE_Genes_table)
  
  #write.table(DE_SigTable,"AllcomparisonDEG_in_Top15GPwindows.txt",quote=F,row.names=F,col.names=T,sep="\t")
  DE_SigTable <- DE_SigTable %>% select(c(CandidateGenes_DE,GO,KEGG)) %>% distinct()
  
  #write.table(DE_SigTable,"AllcomparisonDEG_in_Top15GPwindows_UniqueGenes.txt",quote=F,row.names=F,col.names=T,sep="\t")
  



##########
########## Next iteration. Incoportating Eggnog Annotations, including snpratio from which the Gene was found, and making a column for each DE comparison
##########
  
  get_genes_in_windows <- function(genes,snpratio_wind){
    filtergenes <- data.frame(matrix(nrow = 0,ncol=ncol(genes)+1))
    colnames(filtergenes) <- c(colnames(genes),"snpRatio")
    if(nrow(genes) >=1 & nrow(snpratio_wind)>=1){
      for(x in 1:nrow(genes)){
        for(y in 1:nrow(snpratio_wind)){
          if(genes[x,2]==snpratio_wind[y,1] && genes[x,3] > snpratio_wind[y,4] && genes[x,4] < snpratio_wind[y,5]){filtergenes <- rbind.data.frame(filtergenes,cbind.data.frame(genes[x,],snpRatio=snpratio_wind[y,8]))}
        }
      }
    }
    return(filtergenes)
  }
  
  
DE_SigTable <- cbind.data.frame(WindowSize=NA,SigValue=NA,Slide=NA,CandidateGenes_DE=NA,snpRatio=NA)
DEG_comparisons <- c("H1VSM1","H1VSM1VSM2VSS","M1VSM2","M2VSS","H1VSS","R_PVSS_P")
windowsizes <- c(400,800)
Sig_Thresholds <- c(8,12,16)
slidingType <- c("non-overlap")


for(DE_Method in DEG_comparisons){
    
  DE_Genes_table <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv"))
  DE_Genes_table <- DE_Genes_table %>% filter(significant=="yes")
  DE_Genes <- DE_Genes_table$gene_name
  for(wind in windowsizes){
    for(sig in Sig_Thresholds){
      for(slide in slidingType){
        print(paste(wind,sig,slide))
        step <- ifelse(slide=="overlap",wind/2,0)
        Genes <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/Gene_in_Candidate_region.txt"))
        Genes$`#GeneID` <- gsub("gene-","",Genes$`#GeneID`)
        Snpratio <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/SNP.ratio.txt"))
        #Get genes only in top 3-5 windows
        Snpratio <- Snpratio %>% dplyr::arrange(Ratio)  %>% dplyr::slice_tail(n=20)
        Candidate_Genes <-  get_genes_in_windows(Genes,Snpratio)
        sigDE_in_Candidates <- DE_Genes[which(as.character(DE_Genes) %in% as.character(unlist(Candidate_Genes[,1])))]
        
        if(length(sigDE_in_Candidates)>=1){
          sigDE_in_Candidates <- cbind.data.frame(CandidateGenes_DE=sigDE_in_Candidates,Sig="Significant")
          colnames(sigDE_in_Candidates)[2] <- DE_Method
        if(nrow(DE_SigTable)==0){
          DE_SigTable <- rbind.data.frame(DE_SigTable,cbind.data.frame(WindowSize=wind,SigValue=sig,Slide=slide,CandidateGenes_DE=as.character(unlist(Candidate_Genes[,1])), snpRatio=Candidate_Genes[,6]))
          DE_SigTable <- left_join(DE_SigTable,sigDE_in_Candidates)
        }
          else{
            DE_temptable <- cbind.data.frame(WindowSize=wind,SigValue=sig,Slide=slide,CandidateGenes_DE=as.character(unlist(Candidate_Genes[,1])), snpRatio=Candidate_Genes[,6])
            DE_temptable <- left_join(DE_temptable,sigDE_in_Candidates)
            DE_SigTable <- full_join(DE_SigTable,DE_temptable)
          }
        }
      }
    }
  }
  
}  

DE_SigTable <- DE_SigTable[rowSums(is.na(DE_SigTable))<=6,]
DE_SigTable[is.na(DE_SigTable)] <- 0
DE_SigTable[DE_SigTable=="Significant"] <- 1
DE_SigTable <- DE_SigTable %>% group_by(WindowSize,SigValue,Slide,CandidateGenes_DE,snpRatio) %>% summarise(H1VSM1=max(as.numeric(H1VSM1)),
                                                                                                            H1VSM1VSM2VSS =max(as.numeric(H1VSM1VSM2VSS )),
                                                                                                            R_PVSS_P=max(as.numeric(R_PVSS_P)),
                                                                                                            M1VSM2 =max(as.numeric(M1VSM2 )),
                                                                                                            M2VSS  =max(as.numeric(M2VSS  )),
                                                                                                            H1VSS=max(as.numeric(H1VSS)))



###
DE_file <- DE_SigTable
EL101_to_EL102 <- fread("EL10.1_vs_EL10.2_filter.blast")
EL101_to_EL102 <- EL101_to_EL102[,1:2]
colnames(EL101_to_EL102) <- c("EL10.1","EL10.2")
EL101_to_EL102$EL10.2 <-  unlist(lapply(EL101_to_EL102$EL10.2, function(x) paste(strsplit(x, "\\.")[[1]][1:2], collapse = ".")))
EL101_to_EL102 <- distinct(EL101_to_EL102)

Gene_ID <- fread("GeneID.txt")
colnames(Gene_ID) <- c("gene_id","EL10.1")
#Gene_ID$EL10.1 <- gsub("M","P",Gene_ID$EL10.1)
EL101_to_EL102 <- inner_join(EL101_to_EL102,Gene_ID)
DE_file$gene_id <- as.integer(gsub("LOC","",DE_file$CandidateGenes_DE))
DE_file_up <- inner_join(DE_file,EL101_to_EL102)
DE_file_up <- DE_file_up %>% select(-EL10.1) %>% distinct()

###


DE_SigTable <- cbind.data.frame(DE_file_up,Num_Sig_COmparisons=rowSums(DE_file_up[,6:11],na.rm = T))
colnames(DE_Genes_table)[11] <- "CandidateGenes_DE"
colnames(DE_Genes_table)
DE_Genes_table <- DE_Genes_table[,c(11,12,13,16)]

DE_SigTable <- inner_join(DE_SigTable,DE_Genes_table)
###Need to uplift version or get annotations for EL10.1
Eggnog_Annotations <- fread("MM_vbbc8k1v.emapper.annotations.csv")
Eggnog_Annotations <- Eggnog_Annotations[,c(1,8,9,21)]
colnames(Eggnog_Annotations) <- c("EL10.2","Eggnog_Description","Eggnog_name","PFAM_Domains")
Eggnog_Annotations$EL10.2 <- gsub(".EL10.2_2","",Eggnog_Annotations$EL10.2)
DE_SigTable <- inner_join(DE_SigTable,Eggnog_Annotations) %>% distinct()
#write.table(DE_SigTable,"DEG_in_Top20GPwindows_Eggnog.txt",quote=F,row.names=F,col.names=T,sep="\t")

DE_SigTable_Unique <- DE_SigTable %>% select(-c(WindowSize,SigValue,Slide,snpRatio)) %>% distinct()


DEG_comparisons <- c("H1VSM1","M1VSM2","M2VSS","H1VSS","R_PVSS_P")
for(DE_Method in DEG_comparisons){
  DE_Genes_table <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv"))
  DE_Genes_table <- DE_Genes_table %>% select(gene_name,`log2(fc)`)
  DE_Genes_table <- DE_Genes_table %>% group_by(gene_name) %>% summarise(enrichment=mean(`log2(fc)`,na.rm = T))
  colnames(DE_Genes_table) <- c("CandidateGenes_DE",paste0(DE_Method,"_","log2"))
  DE_SigTable_Unique <- left_join(DE_SigTable_Unique,DE_Genes_table)
}
DE_SigTable_Unique <- DE_SigTable_Unique %>% distinct()

ggplot(DE_SigTable_Unique,aes(x=H1VSS_log2,y=R_PVSS_P_log2))+
  geom_point()

colnames(DE_SigTable_Unique)
cor(DE_SigTable_Unique[,17:21])

plot_pheatmap <- pheatmap(cor(DE_SigTable_Unique[,17:21]),
                          main = "Heatmap",
                          clustering_distance_rows = "correlation",
                          clustering_distance_cols = "correlation",
                          color= colorRampPalette(c("lightblue","white","red"))(50),
                          cluster_rows = TRUE,
                          cluster_cols = TRUE,
                          display_numbers = TRUE,
                          fontsize = 22)
print(plot_pheatmap)



#write.table(DE_SigTable_Unique,"DEG_in_Top20GPwindows_Eggnog_UniqueGenes.txt",quote=F,row.names=F,col.names=T,sep="\t")
  
  
  
  
##############
############## Make Manhattan plots
##############

library(qqman)
GPSEQ_Paval <- fread('SNP.pvalue.txt')
colnames(GPSEQ_Paval)[1] <- "CHROM"
colnames(GPSEQ_Paval)[6] <- "p"
colnames(GPSEQ_Paval)
GPSEQ_Paval$p <- exp(-GPSEQ_Paval$p)

levels(as.factor(GPSEQ_Paval$CHROM))
GPSEQ_Paval$CHROM <- as.numeric(as.factor(GPSEQ_Paval$CHROM ))

GPSEQ_Paval <- na.omit(GPSEQ_Paval)
colors <- c("orange","blue")
GPSEQ_Paval <- GPSEQ_Paval %>%
  group_by(CHROM) %>%
  mutate(entry_number = row_number())
GPSEQ_Paval <- GPSEQ_Paval %>% filter(CHROM<10)

op <- par(mar = c(5, 6, 4, 2) + 0.1,  # bottom, left, top, right
          yaxs = "i")                 # use exact axis range (no extra 4% padding)

# Compute a padded ylim to avoid clipping tall peaks (optional but helpful)
ylim_max <- max(-log10(GPSEQ_Paval$p), na.rm = TRUE)
ylim_pad <- c(0, ylim_max * 1.05)

# Draw the Manhattan plot
manplot <- manhattan(
  GPSEQ_Paval,
  chr = "CHROM",
  bp  = "Pos",
  p   = "p",
  snp = "ID",
  col = colors,
  main = "Manhattan Plot of GPSEQ P-Values",
  cex.axis = 2,
  cex.lab  = 2,
  suggestiveline = FALSE,
  genomewideline = FALSE,
  ylim = ylim_pad
)
manplot
# Restore original par settings
par(op)

manplot
bonferroni_log10 <- -log10(0.05/nrow(GPSEQ_Paval))
abline(h = bonferroni_log10, col = "red", lty = 2) # Red dashed line for significance threshold
# qq(GPSEQ_Paval$p)

vcf <- fread("C:/Users/evan.long/USDA/REE-ARS-PWA-2054-SBGradedPool - General/BCTV_GPSeq_9-27/Analysis_result/Summary/Summary/3_Variant_Calling/1_SNP/vcf/SNP.vcf")
  
Sig_Sites <- GPSEQ_Paval %>% filter(p<10e-10)
colnames(vcf)[1] <- "CHROM"
vcf_clean <- vcf %>% filter(grepl("NC_0792",CHROM))
vcf_clean$CHROM <- as.numeric(as.factor(vcf_clean$CHROM))
colnames(vcf_clean)
Sig_Sites <- Sig_Sites %>% dplyr::select(CHROM,Pos,Ref,Alt,p)
colnames(Sig_Sites) <- c("CHROM","POS","REF","ALT","p")
Sig_Sites <- inner_join(Sig_Sites,vcf_clean)

colnames(Sig_Sites)[11:17] <- c("R(P)","S(P)","H1","H2","M1","M2","S")
Sig_Sites_chr3 <- Sig_Sites %>% filter(CHROM==3)

unique(vcf_clean$CHROM) %in% unique(GPSEQ_Paval$CHROM)
vcf_clean <- vcf_clean %>% dplyr::select(-ID)
colnames(GPSEQ_Paval)[2] <- "POS"
vcf_clean_p <- inner_join(GPSEQ_Paval,vcf_clean)


write.csv(Sig_Sites,"SignificantSites.csv")
# 
# 
# Get_AD <- function(x, entry){
#   x <- as.character(x)
#   AD <- strsplit(strsplit(x,":")[[1]][2],",")[[1]][[entry]]
#   return(as.numeric(AD))
# }
# 
# 
# Get_AD(as.character(test[1,14]),2)
# apply(test[,14:ncol(test)],2,function(y) sapply(y,FUN = function(z) Get_AD(z,2)))
# 
# vcf_clean_REF <- vcf_clean_p
# vcf_clean_ALT <- vcf_clean_p
# 
# vcf_clean_REF[,14:ncol(vcf_clean_REF)] <- apply(vcf_clean_REF[,14:ncol(vcf_clean_REF)],2,function(y) sapply(y,FUN = function(z) Get_AD(z,1)))
# vcf_clean_ALT[,14:ncol(vcf_clean_ALT)] <- apply(vcf_clean_ALT[,14:ncol(vcf_clean_ALT)],2,function(y) sapply(y,FUN = function(z) Get_AD(z,2)))
# 
# library(vcfR)
# vcf <- read.vcfR("C:/Users/evan.long/USDA/REE-ARS-PWA-2054-SBGradedPool - General/BCTV_GPSeq_9-27/Analysis_result/Summary/Summary/3_Variant_Calling/1_SNP/vcf/SNP.vcf")
# 
# poly <- is.polymorphic(vcf, na.omit = TRUE)  # set na.omit=TRUE to ignore missing genotypes
# 
# # Subset to keep only polymorphic rows
# vcf_poly <- vcf
# vcf_poly@fix <- vcf@fix[poly, , drop = FALSE]
# vcf_poly@gt  <- vcf@gt [poly, , drop = FALSE]
# 
# 
# 
# # Extract AD counts: assumes GT:AD format and samples named to match severity_levels
# ad <- extract.gt(vcf_poly, element = "AD", as.numeric = FALSE)
# severity_levels <- c("R(P)","S(P)","H1","H2","M1","M2","S")
# colnames(ad) <- severity_levels
# # GPSEQ_Paval$CHROM <- as.numeric(as.factor(GPSEQ_Paval$CHROM ))
# 
# # Parse "ref,alt" AD strings per sample into numeric counts
# parse_ad_ref <- function(x) {
#   if (is.na(x)) return(c(NA_integer_, NA_integer_))
#   parts <- strsplit(x, ",")[[1]]
#   if (length(parts) < 2) return(c(NA_integer_, NA_integer_))
#   return(as.integer(parts[1]))
# }
# parse_ad_alt <- function(x) {
#   if (is.na(x)) return(c(NA_integer_, NA_integer_))
#   parts <- strsplit(x, ",")[[1]]
#   if (length(parts) < 2) return(c(NA_integer_, NA_integer_))
#   return(as.integer(parts[2]))
# }
# 
# 
# ad_parsed_ref <- apply(ad, c(1, 2), parse_ad_ref)
# ad_parsed_alt <- apply(ad, c(1, 2), parse_ad_alt)
# 
# 
# 



SNPRatio <- read.table("sig8/sig8.win400kb.step0kb/SNP.ratio.txt",header=T)
colnames(SNPRatio)[1] <- "CHROM"
SNPRatio <- SNPRatio %>% filter(CHROM!="NC_002511.2", CHROM!="NC_059012.1")
levels(as.factor(SNPRatio$CHROM ))
head(SNPRatio)
SNPRatio$CHROM <- as.numeric(as.factor(SNPRatio$CHROM ))
SNPRatio <- na.omit(SNPRatio)
colors <- c("orange","blue")
SNPRatio <- SNPRatio %>%
  group_by(CHROM) %>%
  mutate(entry_number = row_number())
SNPRatio <- SNPRatio %>% filter(CHROM<10)



# Define custom colors for alternating chromosomes
chrom_colors <- c("orange", "blue") # Define your own colors as needed

# Convert Chromosome to factor and assign custom colors
SNPRatio$CHROM <- factor(SNPRatio$CHROM, levels = unique(SNPRatio$CHROM))
chrom_colors <- rep(chrom_colors, length.out = length(levels(SNPRatio$CHROM )))
names(chrom_colors) <- levels(SNPRatio$CHROM )

SNPRatio$Relative_Pos <- seq_len(nrow(SNPRatio))
chr_starts <- SNPRatio %>%
  group_by(CHROM) %>%
  summarize(Pos = min(Relative_Pos)) %>%
  ungroup()

new_breaks <- (chr_starts$Pos[-length(chr_starts$Pos)] + chr_starts$Pos[-1]) / 2
new_breaks <- c(new_breaks,(max(SNPRatio$Relative_Pos)+chr_starts$Pos[9])/2)
p1 <- ggplot(SNPRatio, aes(x = Relative_Pos, y = Ratio, col = CHROM)) +
  geom_point(size = 1) +
  scale_color_manual(values = chrom_colors, guide = "none") + 
  ylim(0,0.5) +
  ylab("SNP Ratio") +
  xlab("Chromsome") +
  theme_minimal(base_size = 24) + 
  scale_x_continuous(breaks = new_breaks, labels = chr_starts$CHROM) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = chr_starts$Pos, linetype = "solid", color = "grey", linewidth=0.25)+
  geom_hline(yintercept = 0.15, linetype = "dashed")

p1








###
### Dig into Raj Gene list
###
Raj_Genes <- cbind.data.frame(gene_id=c("Bevul.4G018700.EL10.2_2","Bevul.4G018600.EL10.2_2","Bevul.4G161100.EL10.2_2","Bevul.5G243900.EL10.2_2",
                                        "Bevul.9G034600.EL10.2_2","Bevul.3G055200.EL10.2_2","Bevul.1G087400.EL10.2_2","Bevul.6G032500.EL10.2_2",
                                        "Bevul.2G249600.EL10.2_2","Bevul.8G009000.EL10.2_2"))
DE_Genes_table <-  fread(paste0("C:/Users/evan.long/USDA/REE-ARS-PWA-2054-SBGradedPool - General/Rhizomania_RNA_DE/RVSR_maxVSS_Transcript_differential_expression.csv"))

Raj_Genes <- inner_join(Raj_Genes,DE_Genes_table)

Raj_Genes <- cbind.data.frame(gene_id=Raj_Genes$gene_id,
                              R.FPKM=rowMeans(Raj_Genes[,18:21]), 
                              R_max.FPKM=rowMeans(Raj_Genes[,22:25]),
                              S.FPKM=rowMeans(Raj_Genes[,26:29]),
                              significant=Raj_Genes$significant)



windowsizes <- c(200,400,800)
Sig_Thresholds <- c(12,16,18)
slidingType <- c("non-overlap")
  for(wind in windowsizes){
    for(sig in Sig_Thresholds){
      for(slide in slidingType){
        print(paste(wind,sig,slide))
        step <- ifelse(slide=="overlap",wind/2,0)
        Genes <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/Gene_in_Candidate_region.txt"))
        
        
        Raj_Genes <- Raj_Genes %>% mutate(tempname=gene_id %in% as.character(unlist(Genes[,1])))
        colnames(Raj_Genes)[ncol(Raj_Genes)] <- paste0(wind,"sig",sig,slide)
          }
        }
      }
 
Raj_Genes %>% mutate(tempname=gene_id %in% as.character(unlist(DE_Genes_table$gene_id)))
write.table(Raj_Genes,"SB_SNPGenes_wth_PooledResults.txt",col.names=T,row.names = F,quote=F,sep="\t")







####DEG Correlation
DEG_comparisons <- c("H1VSM1","H1VSM1VSM2VSS","M1VSM2","M2VSS","H1VSS","R_PVSS_P")
windowsizes <- c(400,800)
Sig_Thresholds <- c(8,12,16)
slidingType <- c("non-overlap")


  DE_Method <- "H1VSM1VSM2VSS"
  DE_Genes_table <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv"))
  DEG_Cortable <- DE_Genes_table %>% group_by(t_name) %>% summarise(H=(FPKM.HI_R1+FPKM.HI_R2+FPKM.HI_R3+FPKM.HI_R4)/4,M1=(FPKM.MI_R1+FPKM.MI_R2+FPKM.MI_R3+FPKM.MI_R4)/4,M2=(FPKM.M2_R1+FPKM.M2_R2+FPKM.M2_R3+FPKM.M2_R4)/4,S=(FPKM.S_R1+FPKM.S_R2+FPKM.S_R3+FPKM.S_R4)/4)
  
  DE_Method <- "R_PVSS_P"
  DE_Genes_table <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv"))
  Parents_DEG <- DE_Genes_table %>% group_by(t_name) %>% summarise(R_P=(FPKM.R_P_R1+FPKM.R_P_R2+FPKM.R_P_R3+FPKM.R_P_R4)/4,S_P=(FPKM.S_P_R1+FPKM.S_P_R2+FPKM.S_P_R3+FPKM.S_P_R4)/4)
  DEG_Cortable <- inner_join(DEG_Cortable,Parents_DEG)

  library(pheatmap)
  # Create the heatmap

  DEG_matrix_cor <- cor(DEG_Cortable[,2:7])
  heatmap(DEG_matrix_cor,Colv=NA,col=heat.colors(256))
  
  
  plot_pheatmap <- pheatmap(DEG_matrix_cor,
                            main = "Heatmap",
                            clustering_distance_rows = "correlation",
                            clustering_distance_cols = "correlation",
                            color= colorRampPalette(c("lightblue","white","red"))(50),
                            cluster_rows = TRUE,
                            cluster_cols = TRUE,
                            display_numbers = TRUE,
                            fontsize = 22)
  print(plot_pheatmap) 
  