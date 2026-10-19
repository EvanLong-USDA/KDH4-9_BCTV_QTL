setwd("C:/Users/evan.long/USDA/REE-ARS-PWA-2054-SBGradedPool - General/BCTV_GPSeq_9-27/Analysis_result/Summary/Summary/4_GradedPool-Seq")
library(data.table)
library(dplyr)
library(ggplot2)


####DEG Correlation
# DEG_comparisons <- c("H1VSM1","H1VSM1VSM2VSS","M1VSM2","M2VSS","H1VSS","R_PVSS_P")
# windowsizes <- c(400,800)
# Sig_Thresholds <- c(8,12,16)
# slidingType <- c("non-overlap")


DE_Method <- "H1VSM1VSM2VSS"
DE_Genes_table <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv"))
DEG_Cortable <- DE_Genes_table %>% group_by(t_name) %>% summarise(H=(FPKM.HI_R1+FPKM.HI_R2+FPKM.HI_R3+FPKM.HI_R4)/4,M1=(FPKM.MI_R1+FPKM.MI_R2+FPKM.MI_R3+FPKM.MI_R4)/4,M2=(FPKM.M2_R1+FPKM.M2_R2+FPKM.M2_R3+FPKM.M2_R4)/4,S=(FPKM.S_R1+FPKM.S_R2+FPKM.S_R3+FPKM.S_R4)/4)

DE_Method <- "R_PVSS_P"
DE_Genes_table <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv"))
Parents_DEG <- DE_Genes_table %>% group_by(t_name) %>% summarise(R_P=(FPKM.R_P_R1+FPKM.R_P_R2+FPKM.R_P_R3+FPKM.R_P_R4)/4,S_P=(FPKM.S_P_R1+FPKM.S_P_R2+FPKM.S_P_R3+FPKM.S_P_R4)/4)
DEG_Cortable <- inner_join(DEG_Cortable,Parents_DEG)

library(pheatmap)
library(grid)
# Create the heatmap

write.table(DEG_Cortable,"BCTV_mean_expression.txt",quote=F,row.names=F,col.names=T,sep="\t")

DEG_matrix_cor <- cor(DEG_Cortable[,2:7])
heatmap(DEG_matrix_cor,Colv=NA,col=heat.colors(256))



plot_pheatmap <- pheatmap(DEG_matrix_cor,
                          clustering_distance_rows = "correlation",
                          clustering_distance_cols = "correlation",
                          color= colorRampPalette(c("white","lightblue","red"))(50),
                          cluster_rows = F,
                          cluster_cols = F,
                          display_numbers = TRUE,
                          fontsize = 38,legend_breaks = c(0.88,0.94,1))

grid.newpage()  # Start a new page for plotting
pushViewport(viewport(width = unit(0.9, "npc"), height = unit(0.9, "npc")))  # Full size viewportprint(heatmap_plot)  # Print your heatmap with adjusted margins
print(plot_pheatmap) 
popViewport()  # Return to previous viewport




DE_Method <- "H1VSS"
DE_Genes_table_HS <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv")) %>% dplyr::select(t_name,gene_name,FPKM.HI_R1,FPKM.HI_R2,FPKM.HI_R3,FPKM.HI_R4,FPKM.S_R1,FPKM.S_R2,FPKM.S_R3,FPKM.S_R4,significant,regulation)
colnames(DE_Genes_table_HS)[11] <- "signficant_HS"
colnames(DE_Genes_table_HS)[12] <- "regulation_HS"
DE_Method <- "R_PVSS_P"
DE_Genes_table_P <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv"))  %>% dplyr::select(t_name,gene_name,FPKM.R_P_R1,FPKM.R_P_R2,FPKM.R_P_R3,FPKM.R_P_R4,FPKM.S_P_R1,FPKM.S_P_R2,FPKM.S_P_R3,FPKM.S_P_R4,FPKM.S_P_R5, significant,regulation)
colnames(DE_Genes_table_P)[12] <- "signficant_p"
colnames(DE_Genes_table_P)[13] <- "regulation_p"

DEG_HSP <- inner_join(DE_Genes_table_HS,DE_Genes_table_P) %>% filter(signficant_HS=="yes",signficant_p=="yes",regulation_p==regulation_HS)
###

sig <- 8
wind <- 400
step <- 0
Genes <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/Gene_in_Candidate_region.txt"))
Genes[,1] <- unlist(lapply(Genes[,1],function(x) gsub(pattern = "gene-",replacement = "",x)))
Genes <- Genes[,1] 
colnames(Genes) <- "gene_name"
DEG_HSP_GPseq <- inner_join(DEG_HSP,Genes)

DEG_HSP_GPseq_heat <- DEG_HSP_GPseq %>% dplyr::select(FPKM.R_P_R1,FPKM.R_P_R2,FPKM.R_P_R3,FPKM.R_P_R4,FPKM.HI_R1,FPKM.HI_R2,FPKM.HI_R3,FPKM.HI_R4,FPKM.S_P_R1,FPKM.S_P_R2,FPKM.S_P_R3,FPKM.S_P_R4,FPKM.S_P_R5,FPKM.S_R1,FPKM.S_R2,FPKM.S_R3,FPKM.S_R4)

DEG_HSP_GPseq_heat[,c(1:4,9:13)] <- as.data.frame(t(apply(DEG_HSP_GPseq_heat[,c(1:4,9:13)], 1, scale)))
DEG_HSP_GPseq_heat[,c(5:8,14:17)] <- as.data.frame(t(apply(DEG_HSP_GPseq_heat[,c(5:8,14:17)], 1, scale)))

#scaled_DEG_HSP_GPseq_heat <- t(apply(DEG_HSP_GPseq_heat, 1, scale))
scaled_DEG_HSP_GPseq_heat <- DEG_HSP_GPseq_heat
#breaks <- seq(min(DEG_HSP_GPseq_heat), max(DEG_HSP_GPseq_heat), length.out = 51)
#hist(scaled_DEG_HSP_GPseq_heat)

colnames(scaled_DEG_HSP_GPseq_heat) <- gsub("FPKM.","",colnames(DEG_HSP_GPseq_heat))
rownames(scaled_DEG_HSP_GPseq_heat) <- DEG_HSP_GPseq$t_name

scaled_DEG_HSP_GPseq_heat <- as.data.frame(scaled_DEG_HSP_GPseq_heat)
scaled_DEG_HSP_GPseq_heat <- scaled_DEG_HSP_GPseq_heat %>% mutate(Gene=DEG_HSP_GPseq$t_name,rank=(S_R1+S_R2+S_R3+S_R4+S_P_R1+S_P_R2+S_P_R3+S_P_R4+S_P_R5) - (R_P_R1+R_P_R2+R_P_R3+R_P_R4+HI_R1+HI_R2+HI_R3+HI_R4)) %>% arrange(rank) %>% dplyr::select(-rank)
scaled_DEG_HSP_GPseq_heat <- scaled_DEG_HSP_GPseq_heat %>% rowwise %>% mutate(rank=min(S_R1,S_R2,S_R3,S_R4,S_P_R1,S_P_R2,S_P_R3,S_P_R4,S_P_R5)) %>% arrange(rank) %>% dplyr::select(-rank)
scaled_DEG_HSP_GPseq_heat$Gene
colnames(scaled_DEG_HSP_GPseq_heat) <- gsub("FPKM.","",colnames(scaled_DEG_HSP_GPseq_heat))
scaled_DEG_HSP_GPseq_heat <- as.data.frame(scaled_DEG_HSP_GPseq_heat)
rownames(scaled_DEG_HSP_GPseq_heat) <- scaled_DEG_HSP_GPseq_heat$Gene
scaled_DEG_HSP_GPseq_heat <- scaled_DEG_HSP_GPseq_heat[,1:17]
#hist(as.numeric(unlist(scaled_DEG_HSP_GPseq_heat)))#hist(as.numeric(unmax()list(scaled_DEG_HSP_GPseq_heat)))
breaks <- c(min(scaled_DEG_HSP_GPseq_heat),seq(-1,0,length.out = 25),
            seq(0.075,1,length.out = 23),2,3)
plot_pheatmap <- pheatmap(scaled_DEG_HSP_GPseq_heat,
                          main = "Heatmap",
                          clustering_distance_rows = "correlation",
                          clustering_distance_cols = "correlation",
                          color= colorRampPalette(c("blue","white","red"))(50),
                          breaks=breaks,
                          cluster_rows = F,
                          cluster_cols = F,
                          display_numbers = F,
                          fontsize = 22)
print(plot_pheatmap) 



EL101_to_EL102 <- fread("EL10.1_vs_EL10.2_filter.blast")
EL101_to_EL102 <- EL101_to_EL102[,1:2]
colnames(EL101_to_EL102) <- c("EL10.1","EL10.2")
EL101_to_EL102$EL10.2 <-  unlist(lapply(EL101_to_EL102$EL10.2, function(x) paste(strsplit(x, "\\.")[[1]][1:2], collapse = ".")))
EL101_to_EL102 <- distinct(EL101_to_EL102)

Gene_ID <- fread("GeneID.txt")
colnames(Gene_ID) <- c("gene_id","EL10.1")
#Gene_ID$EL10.1 <- gsub("M","P",Gene_ID$EL10.1)
EL101_to_EL102 <- inner_join(EL101_to_EL102,Gene_ID)

Expr_ID <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv")) %>% dplyr::select(gene_id,t_name)

EL101_to_EL102 <- inner_join(EL101_to_EL102,Expr_ID)
EL101_to_EL102





EL10.2_names <- EL101_to_EL102$EL10.2[ match(rownames(scaled_DEG_HSP_GPseq_heat),EL101_to_EL102$t_name)]
#######
Expr_ID
Gene_ID
Gene_names <- inner_join(Expr_ID,Gene_ID)
# Gene_names[Gene_names$t_name=="XM_048640834.1",]
# Gene_names[Gene_names$t_name=="XM_048640634.1",]
# Gene_names[Gene_names$t_name=="XM_048640752.1",]

####### Write table of 21 genes annotations
EL10.2_names[is.na(EL10.2_names)] <- c("Bevul.3G023000","OMIT","Bevul.3G012500")
EL10.2_names_df <- cbind.data.frame(rownames(scaled_DEG_HSP_GPseq_heat),EL10.2_names)
colnames(EL10.2_names_df) <- c("t_name","EL10.2")

DEG_HSP_GPseq$t_name
DEG_HSP_GPseq_EL10.2 <- left_join(DEG_HSP_GPseq,EL10.2_names_df) %>% distinct()
DEG_HSP_GPseq_EL10.2
DEG_HSP_GPseq_EL10.2 <- DEG_HSP_GPseq_EL10.2[DEG_HSP_GPseq_EL10.2$EL10.2!="OMIT"]
DEG_HSP_GPseq_EL10.2 <- DEG_HSP_GPseq_EL10.2 %>% dplyr::select(-c(signficant_HS,signficant_p,regulation_HS,regulation_p))

Expr_ID <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv")) %>% dplyr::select(gene_id,t_name,GO,KEGG,Description)
DEG_HSP_GPseq_EL10.2 <- left_join(DEG_HSP_GPseq_EL10.2,Expr_ID)

Eggnog_Annotations <- fread("MM_vbbc8k1v.emapper.annotations.csv")
Eggnog_Annotations <- Eggnog_Annotations[,c(1,8,9,21)]
colnames(Eggnog_Annotations) <- c("EL10.2","Eggnog_Description","Eggnog_name","PFAM_Domains")
Eggnog_Annotations$EL10.2 <- gsub(".EL10.2_2","",Eggnog_Annotations$EL10.2)
DEG_HSP_GPseq_EL10.2 <- left_join(DEG_HSP_GPseq_EL10.2,Eggnog_Annotations)


write.csv(DEG_HSP_GPseq_EL10.2,"21Gene_Annotations.csv",row.names = F)
###

Genes <- fread(paste0("sig",sig,"/sig",sig,".win",wind,"kb.step",step,"kb/Gene_in_Candidate_region.txt"))
Genes[,1] <- unlist(lapply(Genes[,1],function(x) gsub(pattern = "gene-",replacement = "",x)))
Genes <- Genes[,1] 
colnames(Genes) <- "gene_name"
Expr_ID <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv")) %>% dplyr::select(t_name,gene_name,GO,KEGG,Description)
GPseq_Genes <- left_join(Genes,Expr_ID)
GPseq_Genes <- inner_join(GPseq_Genes,EL101_to_EL102)
GPseq_Genes <- inner_join(GPseq_Genes,Eggnog_Annotations)
GPseq_Genes <- GPseq_Genes %>% dplyr::select(-t_name,-gene_name,-EL10.1) %>% distinct()

write.csv(GPseq_Genes,"GPseq_Genes_Annotations.csv",row.names = F)

EL101_to_EL102
DE_Genes_table_HS <- inner_join(EL101_to_EL102,DE_Genes_table_HS) %>% dplyr::select(-EL10.1,-t_name,-gene_id,-gene_name) %>% distinct()
DE_Genes_table_P <- inner_join(EL101_to_EL102,DE_Genes_table_P) %>% dplyr::select(-EL10.1,-t_name,-gene_id,-gene_name) %>% distinct()


write.csv(DE_Genes_table_HS,"DE_Genes_table_HS_EL10.2.csv")
write.csv(DE_Genes_table_P,"DE_Genes_table_P_EL10.2.csv")










EL10.2_names_df

rownames(scaled_DEG_HSP_GPseq_heat)


scaled_DEG_HSP_GPseq_heat <- as.data.frame(scaled_DEG_HSP_GPseq_heat)
rownames(scaled_DEG_HSP_GPseq_heat)[!is.na(EL10.2_names)] <- EL10.2_names[!is.na(EL10.2_names)]
scaled_DEG_HSP_GPseq_heat <- scaled_DEG_HSP_GPseq_heat[rownames(scaled_DEG_HSP_GPseq_heat)!="OMIT",]

nrow(scaled_DEG_HSP_GPseq_heat)
plot_pheatmap <- pheatmap(scaled_DEG_HSP_GPseq_heat,
                          clustering_distance_rows = "correlation",
                          clustering_distance_cols = "correlation",
                          color= colorRampPalette(c("blue","white","red"))(50),
                          breaks=breaks,
                          cluster_rows = F,
                          cluster_cols = F,
                          display_numbers = F,
                          fontsize = 22)
print(plot_pheatmap) 

library(grid)
heatmap_plot <- pheatmap(scaled_DEG_HSP_GPseq_heat,
                         clustering_distance_rows = "correlation",
                         clustering_distance_cols = "correlation",
                         color= colorRampPalette(c("blue","white","red"))(50),
                         breaks=breaks,
                         cluster_rows = F,
                         cluster_cols = F,
                         display_numbers = F,
                         fontsize = 22,
                         silent = TRUE)  # Prevent immediate plotting
heatmap_plot
grid.newpage()  # Start a new page for plotting
pushViewport(viewport(width = unit(0.9, "npc"), height = unit(0.9, "npc")))  # Full size viewportprint(heatmap_plot)  # Print your heatmap with adjusted margins
print(heatmap_plot)  # Print your heatmap with adjusted layout
popViewport()  # Return to previous viewport

scaled_DEG_HSP_GPseq_heat


SigGenesUniqKDH <- fread("KDH49_UniqueGenes_QTL_quant.txt")
SigGenesUniqKDH$Gene <- gsub("BetaVulgaris_Chromosome","BvKDH49.",SigGenesUniqKDH$Gene)
SigGenesUniqKDH$Gene <- gsub("_","G",SigGenesUniqKDH$Gene)

SigGenesUniqKDH <- SigGenesUniqKDH %>% select(Sample,Gene,TPM) %>% distinct()
#Spread Sample by TPM
SigGenesUniqKDH <- SigGenesUniqKDH %>% group_by(Gene) %>% pivot_wider(names_from = Sample, values_from = TPM, values_fill = 0)
colnames(scaled_DEG_HSP_GPseq_heat)
SigGenesUniqKDH <- SigGenesUniqKDH %>% select(Gene,colnames(scaled_DEG_HSP_GPseq_heat))
SigGenesUniqKDH <- SigGenesUniqKDH %>% filter(Gene %in% c("BvKDH49.3G000200.1","BvKDH49.3G000201.1"))
SigGenesUniqKDH$Gene <- gsub("\\.1","",SigGenesUniqKDH$Gene)

SigGenesUniqKDH <- as.data.frame(SigGenesUniqKDH)
rownames(SigGenesUniqKDH) <- SigGenesUniqKDH$Gene
SigGenesUniqKDH <- SigGenesUniqKDH[,-1]
# SigGenesUniqKDH <- apply(SigGenesUniqKDH,1,scale)
# SigGenesUniqKDH <- t(SigGenesUniqKDH)

SigGenesUniqKDH[,c(1:4,9:13)] <- as.data.frame(t(apply(SigGenesUniqKDH[,c(1:4,9:13)], 1, scale)))
SigGenesUniqKDH[,c(5:8,14:17)] <- as.data.frame(t(apply(SigGenesUniqKDH[,c(5:8,14:17)], 1, scale)))



colnames(scaled_DEG_HSP_GPseq_heat) 
colnames(SigGenesUniqKDH) <- colnames(scaled_DEG_HSP_GPseq_heat)

scaled_DEG_HSP_GPseq_heat_wthKDH <- rbind.data.frame(SigGenesUniqKDH,scaled_DEG_HSP_GPseq_heat)

colnames(scaled_DEG_HSP_GPseq_heat_wthKDH) <- gsub(colnames(scaled_DEG_HSP_GPseq_heat_wthKDH) ,pattern = "R_P_",replacement = "R(P)_")
colnames(scaled_DEG_HSP_GPseq_heat_wthKDH) <- gsub(colnames(scaled_DEG_HSP_GPseq_heat_wthKDH) ,pattern = "S_P_",replacement = "S(P)_")
colnames(scaled_DEG_HSP_GPseq_heat_wthKDH) <- gsub(colnames(scaled_DEG_HSP_GPseq_heat_wthKDH) ,pattern = "HI",replacement = "H")


heatmap_plot <- pheatmap(scaled_DEG_HSP_GPseq_heat_wthKDH,
                         clustering_distance_rows = "correlation",
                         clustering_distance_cols = "correlation",
                         color= colorRampPalette(c("blue","white","red"))(50),
                         breaks=breaks,
                         cluster_rows = F,
                         cluster_cols = F,
                         display_numbers = F,
                         fontsize = 22,
                         silent = TRUE)  # Prevent immediate plotting
heatmap_plot
grid.newpage()  # Start a new page for plotting
pushViewport(viewport(width = unit(0.9, "npc"), height = unit(0.9, "npc")))  # Full size viewportprint(heatmap_plot)  # Print your heatmap with adjusted margins
print(heatmap_plot)  # Print your heatmap with adjusted layout
popViewport()  # Return to previous viewport
















#######

Sig_Genes <- rownames(scaled_DEG_HSP_GPseq_heat)


KEGG_Annotation <- read.csv("C:/Users/evan.long/USDA/REE-ARS-PWA-2054-SBGradedPool - General/BCTV_GPSeq_9-27/Analysis_result_new/Summary/summary/6_differential_expression/R_PVSS_P/2_KEGG_Enrichment/R_PVSS_P_KEGG_enrichment_gene.csv") %>% dplyr::select(pathway_name,Gene_ID,Gene_name)
Eggnog_Annotations <- fread("MM_vbbc8k1v.emapper.annotations.csv")
Eggnog_Annotations <- Eggnog_Annotations[,c(1,8,9,21)]
colnames(Eggnog_Annotations) <- c("EL10.2","Eggnog_Description","Eggnog_name","PFAM_Domains")
Eggnog_Annotations$EL10.2 <- gsub(".EL10.2_2","",Eggnog_Annotations$EL10.2)

missingGenes <- Gene_names %>% filter(t_name %in% c("XM_048640834.1","XM_048640752.1"))
missingGenes <- cbind.data.frame(missingGenes,EL10.2=c("Bevul.3G023000","Bevul.3G012500"))
EL101_to_EL102 <- EL101_to_EL102 %>% dplyr::select(gene_id,t_name,EL10.1,EL10.2)
EL101_to_EL102 <- rbind.data.frame(EL101_to_EL102,missingGenes)
colnames(EL101_to_EL102)[1] <- "Gene_ID"
Sig_Genes_Annotation <- full_join(EL101_to_EL102,Eggnog_Annotations) %>% filter(EL10.2 %in% Sig_Genes) %>% dplyr::select(-t_name,-Gene_ID,-EL10.1) %>% distinct()
Sig_Genes_Annotation

######
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
windowsizes <- c(400)
Sig_Thresholds <- c(8)
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
DE_file$Gene_ID <- as.integer(gsub("LOC","",DE_file$CandidateGenes_DE))
DE_file_up <- inner_join(DE_file,EL101_to_EL102)
DE_file_up <- DE_file_up %>% dplyr::select(-EL10.1,-t_name) %>% distinct()
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

DE_SigTable_Unique <- DE_SigTable %>% dplyr::select(-c(WindowSize,SigValue,Slide,snpRatio)) %>% distinct()


DEG_comparisons <- c("H1VSM1","M1VSM2","M2VSS","H1VSS","R_PVSS_P")
for(DE_Method in DEG_comparisons){
  DE_Genes_table <-  fread(paste0("DEG/",DE_Method,"_Transcript_differential_expression.csv"))
  DE_Genes_table <- DE_Genes_table %>% dplyr::select(gene_name,`log2(fc)`)
  DE_Genes_table <- DE_Genes_table %>% group_by(gene_name) %>% summarise(enrichment=mean(`log2(fc)`,na.rm = T))
  colnames(DE_Genes_table) <- c("CandidateGenes_DE",paste0(DE_Method,"_","log2"))
  DE_SigTable_Unique <- left_join(DE_SigTable_Unique,DE_Genes_table)
}
DE_SigTable_Unique <- DE_SigTable_Unique %>% distinct()




#####






###KEGG Enrichment data
KEGG_HS <- read.csv("C:/Users/evan.long/USDA/REE-ARS-PWA-2054-SBGradedPool - General/BCTV_GPSeq_9-27/Analysis_result_new/Summary/summary/6_differential_expression/H1VSS/2_KEGG_Enrichment/H1VSS_KEGG_enrichment.csv")
KEGG_PP <- read.csv("C:/Users/evan.long/USDA/REE-ARS-PWA-2054-SBGradedPool - General/BCTV_GPSeq_9-27/Analysis_result_new/Summary/summary/6_differential_expression/R_PVSS_P/2_KEGG_Enrichment/R_PVSS_P_KEGG_enrichment.csv")

KEGG_HS <- KEGG_HS %>% mutate(RichFactor=S.gene.number/B.gene.number,Comparison="H vs S",GeneNumber=S.gene.number) %>% dplyr::select(pathway_name,RichFactor,pvalue,Comparison,GeneNumber) 
KEGG_PP <- KEGG_PP %>% mutate(RichFactor=S.gene.number/B.gene.number,Comparison="R_P vs S_P",GeneNumber=S.gene.number) %>% dplyr::select(pathway_name,RichFactor,pvalue,Comparison,GeneNumber) 

KEGG_input <- rbind.data.frame(KEGG_PP,KEGG_HS) %>% group_by(pathway_name) %>% filter(sum(pvalue<=0.05)>0)

KEGG_input
KEGG_levels <- KEGG_input %>% group_by(pathway_name) %>% summarize(pvalue=mean(pvalue),RichFactor=mean(RichFactor)) %>% arrange(-pvalue,RichFactor) %>% dplyr::select(pathway_name) 
KEGG_input$pathway_name <- factor(KEGG_input$pathway_name,levels = KEGG_levels$pathway_name)

levels(KEGG_input$pathway_name)
KEGG_input <- KEGG_input %>% filter(pathway_name %in% levels(pathway_name)[12:24])
# ggplot(KEGG_input, aes(y=pathway_name, x=RichFactor)) +
#   geom_point(aes(color=pvalue, size=GeneNumber + 3, shape=Comparison)) + # Increase base size here
#   scale_color_gradientn(colors = c("blue", "purple", "red"), 
#                         values = scales::rescale(c(0, 0.05, 0.2)), 
#                         limits = c(0, 0.2), 
#                         name="P-Value") +
#   scale_size_continuous(range=c(4, 8), name="Gene Number") + # Adjust range if needed
#   guides(shape = guide_legend(override.aes = list(size=5))) + # Adjust shape size here
#   labs(title="KEGG Pathway Enrichment",
#        x="Rich Factor",
#        y="Pathway Name") +
#   theme_minimal() +
#   theme(text=element_text(size=22))



KEGG_input$Comparison <- gsub("R_P","R(P)",KEGG_input$Comparison)
KEGG_input$Comparison <- gsub("S_P","S(P)",KEGG_input$Comparison)

KEGG_input <- KEGG_input %>% filter(!pathway_name %in% c("Photosynthesis - antenna proteins","Nitrogen metabolism","Tyrosine metabolism"))

ggplot(KEGG_input, aes(y=pathway_name, x=RichFactor)) +
  geom_point(aes(color=pvalue, size=GeneNumber, shape=Comparison)) +
  scale_color_gradientn(colors = c("blue", "purple", "red"), 
                        values = scales::rescale(c(0, 0.05, 0.2)), 
                        limits = c(0, 0.2), 
                        name="P-Value",breaks = c(0, 0.1, 0.2)) +
  scale_size_continuous(range=c(4, 20), name="Gene Number") + # Increased range for better gradation
  guides(shape = guide_legend(override.aes = list(size=5))) + # Adjust shape size here
  labs(title="KEGG Pathway Enrichment",
       x="Rich Factor",
       y="Pathway Name") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 26, face = "bold"),  # Title size
    axis.title.x = element_text(size = 28),  # X-axis title size
    axis.title.y = element_text(size = 28),  # Y-axis title size
    axis.text.x = element_text(size = 26),  # X-axis text size
    axis.text.y = element_text(size = 26),  # Y-axis text size
    legend.title = element_text(size = 26),  # Legend title size
    legend.text = element_text(size = 22)  # Legend text size
  )











### Viral Load Data

BCTV_Viral_Load <- read.table("Viral_readCount.txt")
colnames(BCTV_Viral_Load) <- c("Sample","Total_Reads","Viral_Reads","Viral_Ratio")

BCTV_Viral_Load <- BCTV_Viral_Load %>% mutate(Pool=gsub("_[^_]+$", "", Sample))
unique(BCTV_Viral_Load$Pool)
BCTV_Viral_Load$Pool <-gsub("I","1",BCTV_Viral_Load$Pool)

BCTV_Viral_Load$Pool <- gsub("H1","H",BCTV_Viral_Load$Pool)
BCTV_Viral_Load$Pool <- gsub("R_P","R(P)",BCTV_Viral_Load$Pool)
BCTV_Viral_Load$Pool <- gsub("S_P","S(P)",BCTV_Viral_Load$Pool)

BCTV_Viral_Load$Pool <- factor(BCTV_Viral_Load$Pool, levels = c("H","M1","M2","S","R(P)","S(P)"))
ggplot(BCTV_Viral_Load, aes(x=Pool, y=Viral_Ratio)) +
  geom_boxplot(aes(fill=Pool)) +
  labs(title="Viral Load in BCTV Samples",
       x="Pool",
       y="Viral Ratio") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 32, face = "bold"),  # Title size
    axis.title.x = element_text(size = 28),  # X-axis title size
    axis.title.y = element_text(size = 28),  # Y-axis title size
    axis.text.x = element_text(size = 16),  # X-axis text size
    axis.text.y = element_text(size = 26)  # Y-axis text size
  )
















# Load necessary libraries
# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(multcompView)
library(dplyr)

# Perform pairwise comparisons using t-tests
stat_test <- compare_means(Viral_Ratio ~ Pool, data = BCTV_Viral_Load, method = "t.test")

# Adjust p-values using Bonferroni correction

# Create a function to generate compact letter display (CLD)
generate_cld <- function(stat_test) {
  pvals <- stat_test$p
  names(pvals) <- paste(stat_test$group1, stat_test$group2, sep = "-")
  cld <- multcompLetters(pvals)$Letters
  return(cld)
}

# Generate CLD
cld <- generate_cld(stat_test)

# Create a data frame for CLD
cld_df <- data.frame(Pool = names(cld), group = cld)

BCTV_Viral_Load$Pool <- factor(BCTV_Viral_Load$Pool, levels = c("H","M1","M2","S","R(P)","S(P)"))

# Create the plot with statistical groupings
ggplot(BCTV_Viral_Load, aes(x=Pool, y=Viral_Ratio)) +
  geom_boxplot(aes(fill=Pool)) +
  geom_text(data = cld_df, aes(label = group, y = max(BCTV_Viral_Load$Viral_Ratio) + 0.005), size = 10) +
  labs(
       x="Pool",
       y="Viral Ratio") +
  theme_minimal() +
  theme(
    base_size = 28,  # Base font size for all text
    plot.title = element_text(size = 32, face = "bold"),  # Title size
    axis.title.x = element_text(size = 28),  # X-axis title size
    axis.title.y = element_text(size = 28),  # Y-axis title size
    axis.text.x = element_text(size = 26),  # X-axis text size
    axis.text.y = element_text(size = 26),  # Y-axis text size
    legend.title = element_text(size = 28),  # Legend title size
    legend.text = element_text(size = 26)  # Legend text size
  )
