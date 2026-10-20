setwd("C:/Users/evan.long/USDA/REE-ARS-PWA-2054-SBGradedPool - General/Assemblies/GENESPACE")
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)


Otholog_table <- fread("KDH4-9_ortholog.txt")
colnames(Otholog_table)

test<- Otholog_table %>% group_by(genome,id) %>% summarize(count=n())

length(unique(Otholog_table$id))
Otholog_table <- Otholog_table %>% filter(flag!="NSOrtho")
length(unique(Otholog_table$id))

Otholog_table <- Otholog_table %>% group_by(genome,id) %>% 
  mutate(samechrom=ifelse(chr==interpChr, T, F)) %>% 
  mutate(MultipleEntries=ifelse(n()>1,T,F)) %>%
  mutate(BestCopy=ifelse(samechrom==F & MultipleEntries==T,"no","yes"))

Otholog_table <- Otholog_table %>% filter(BestCopy=="yes") %>% select(-BestCopy,-samechrom,-MultipleEntries)
test<- Otholog_table %>% group_by(genome,id) %>% summarize(count=n())


annotate_pan_homology <- function(
    df,
    genome1,
    genome2,
    group_col = "pgID",            # set to "pgRepID" if that defines the ortholog group in your real data
    primary_flag = "PASS",
    array_flags = "array",
    allow_array_as_homolog = FALSE,  # keep FALSE per your requirement
    tie_break = c("by_ord", "by_start", "first"),
    normalize_flags = TRUE,
    classification_col = "classification",
    reason_col = "classification_reason"
) {
  tie_break <- match.arg(tie_break)
  
  # --- Validate columns ---
  req <- c(group_col, "genome", "flag", "id", "chr", "start", "end", "ord")
  miss <- setdiff(req, names(df))
  if (length(miss)) {
    stop("Input df missing required columns: ", paste(miss, collapse = ", "))
  }
  
  suppressPackageStartupMessages({
    library(dplyr)
  })
  
  # --- Filter to the two genomes of interest ---
  df2 <- df %>%
    filter(.data$genome %in% c(genome1, genome2)) %>%
    mutate(.row_id = dplyr::row_number())
  
  if (nrow(df2) == 0) {
    df2[[classification_col]] <- NA_character_
    df2[[reason_col]] <- NA_character_
    return(df2)
  }
  
  # --- Normalize flags (robust to case/whitespace) ---
  if (normalize_flags) {
    df2 <- df2 %>% mutate(flag_norm = tolower(trimws(.data$flag)))
    primary_flag_norm <- tolower(trimws(primary_flag))
    array_flags_norm  <- tolower(trimws(array_flags))
  } else {
    df2 <- df2 %>% mutate(flag_norm = .data$flag)
    primary_flag_norm <- primary_flag
    array_flags_norm  <- array_flags
  }
  
  grp <- dplyr::sym(group_col)
  
  # --- Presence metrics per group ---
  metrics <- df2 %>%
    group_by(!!grp) %>%
    summarise(
      n1 = sum(.data$genome == !!genome1),
      n2 = sum(.data$genome == !!genome2),
      .groups = "drop"
    )
  
  # --- Choose at most one homolog per (group, genome); arrays excluded unless allowed ---
  ordered <- df2 %>%
    mutate(
      is_array = .data$flag_norm %in% array_flags_norm,
      eligible = if (allow_array_as_homolog) TRUE else !is_array,  # arrays excluded by default
      priority = dplyr::case_when(
        eligible & (.data$flag_norm == primary_flag_norm) ~ 1L,   # PASS (non-array) first
        eligible                                           ~ 2L,   # other non-array flags
        TRUE                                               ~ 9L    # non-eligible (arrays)
      ),
      tb = dplyr::case_when(
        tie_break == "by_ord"   ~ as.numeric(.data$ord),
        tie_break == "by_start" ~ as.numeric(.data$start),
        TRUE                    ~ as.numeric(.data$ord)  # "first" handled later
      )
    ) %>%
    group_by(!!grp, .data$genome) %>%
    { if (tie_break == "first") . else arrange(., .data$priority, .data$tb, .data$id, .by_group = TRUE) } %>%
    summarise(
      pick_row_id = {
        cand <- ifelse(priority < 9, .row_id, NA_integer_)
        cand <- cand[!is.na(cand)]
        if (length(cand) >= 1) cand[1] else NA_integer_
      },
      .groups = "drop"
    ) %>%
    filter(!is.na(.data$pick_row_id))
  
  picks <- ordered$pick_row_id
  
  # --- Classify rows ---
  out <- df2 %>%
    left_join(metrics, by = group_col) %>%
    mutate(
      !!classification_col := dplyr::case_when(
        n1 > 0 & n2 == 0 & .data$genome == !!genome1 ~ "Unique",
        n2 > 0 & n1 == 0 & .data$genome == !!genome2 ~ "Unique",
        .data$.row_id %in% picks                      ~ "Homologous",
        TRUE                                          ~ "Duplicate"
      ),
      !!reason_col := dplyr::case_when(
        n1 > 0 & n2 == 0 & .data$genome == !!genome1 ~ paste0(
          group_col, " present only in ", genome1, " (n=", n1, "); absent from ", genome2
        ),
        n2 > 0 & n1 == 0 & .data$genome == !!genome2 ~ paste0(
          group_col, " present only in ", genome2, " (n=", n2, "); absent from ", genome1
        ),
        .data$.row_id %in% picks ~ paste0(
          "Chosen Homologous copy in ", .data$genome,
          if (!allow_array_as_homolog) " (arrays excluded)" else "",
          " by priority (PASS>other non-array), tie-break=", tie_break
        ),
        .data$flag_norm %in% array_flags_norm ~ "Tandem duplicate (flag in array_flags)",
        TRUE ~ "Duplicate context (more than one copy in genome for this group or no eligible non-array homolog)"
      )
    ) %>%
    select(-.data$flag_norm, -.data$.row_id) %>%
    arrange(!!grp, .data$genome, .data$ord)
  
  return(out)
}



# Annotated_Ortholog_Table <- annotate_pgid_genes_homology(Otholog_table,"EL102","KDH49")
# 
# Annotated_Ortholog_Table
# 

write.csv(Otholog_table,"test3.csv")



library(readr)
library(dplyr)

df2 <- readr::read_csv("test3.csv", show_col_types = FALSE)

# Use pgID as the group (matches your original spec)
annot_pgID <- annotate_pan_homology(
  df2,
  genome1 = "EL102",
  genome2 = "KDH49",
  group_col = "pgID",            # try "pgRepID" for your real dataset if duplicates span pgIDs
  primary_flag = "PASS",
  array_flags = "array",
  tie_break = "by_ord",
  allow_array_as_homolog = FALSE
)

annot_pgID %>%
  select(pgID, genome, id, flag, ord, classification, classification_reason)




annot_pgID <- annot_pgID %>% mutate(anno=ifelse(
  classification=="Homologous",
                         "mt:v;mc:black;ms:1;tp:0.02;ts:8;tf:Droid Sans;tc:black",
                         ifelse(classification=="Duplicate",
                                "mt:v;mc:yellow;ms:1;tp:0.02;ts:8;tf:Droid Sans;tc:yellow",
                                "mt:v;mc:red;ms:1;tp:0.02;ts:8;tf:Droid Sans;tc:red")))

EL102_Markers <- annot_pgID %>% filter(genome=="EL102") %>% select(pgID,chr,start,end,genome,id,classification)
KDH49_Markers <- annot_pgID %>% filter(genome=="KDH49") %>% select(pgID,chr,start,end,genome,id,classification)


syri_markers <- annot_pgID %>% select(chr,start,end,genome,anno)
write.table(syri_markers,"syri_markers.txt",sep="\t",quote=F,col.names=F,row.names=F)


####
#### Chromosome 3 Top QTL Region
####



Chrom3_QTL_pgID <- EL102_Markers[EL102_Markers$chr=="Chromosome3" & EL102_Markers$start>4000000 & EL102_Markers$end < 4400000,1] 
Chrom3_QTL_Otholog_table <- annot_pgID %>% filter(chr=="Chromosome3" & ((genome=="EL102" & start>4000000 & end < 4400000) | (genome=="KDH49" & start>4065998 & end < 4517556)))
Chrom3_QTL_count <- Chrom3_QTL_Otholog_table %>% 
  group_by(pgID,genome) %>% 
  summarise(count=n()) %>% 
  spread(genome, count)
Chrom3_QTL_count[is.na(Chrom3_QTL_count)] <- 0
sum(Chrom3_QTL_count$EL102==0 & Chrom3_QTL_count$KDH49 > 0) 


EL102_Markers_QTL <- EL102_Markers %>% filter(chr=="Chromosome3" & start>4000000 & end < 4400000)
KDH49_Markers_QTL <- KDH49_Markers %>% filter(chr=="Chromosome3" & start>4065998 & end < 4517556)
QTL_markers <- rbind(EL102_Markers_QTL,KDH49_Markers_QTL)

## Make plot with both genomes
QTL_markers$genome <- factor(QTL_markers$genome,levels=c("KDH49","EL102"))


ggplot(QTL_markers,aes(y=chr,x=start,color=classification)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~genome,scales="free",ncol=1) + 
  scale_color_manual(values=c("Unique"="red","Duplicate"="yellow","Homologous"="black")) + theme(legend.position = "top") + 
  labs(title="Gene Content Variation in Chromosome 3 QTL Region",x="Chromosome",y="Position (bp)",color="Gene Type") + 
  theme(plot.title = element_text(hjust = 0.5, size=16, face="bold"),
        axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        strip.text = element_text(size=14,face="bold"))
table(KDH49_Markers_QTL$classification)
table(EL102_Markers_QTL$classification)








ggplot(QTL_markers, aes(y = chr, x = start, color = classification)) +
  geom_point() +
  ggrepel::geom_text_repel(
    aes(label = pgID),
    size = 3,
    max.overlaps = 100,         # increase if you want more labels shown
    min.segment.length = 0,     # always draw a segment
    segment.size = 0.2,
    box.padding = 0.15,
    point.padding = 0.15,
    nudge_x = 0,                # usually not needed with repel
    show.legend = FALSE
  ) +
  theme_bw() +
  facet_wrap(~genome, scales = "free", ncol = 1) +
  scale_color_manual(values = c("Unique"="red","Duplicate"="yellow","Homologous"="black")) +
  theme(legend.position = "top") +
  labs(
    title = "Gene Content Variation in Chromosome 3 QTL Region",
    x = "Chromosome",
    y = "Position (bp)",
    color = "Gene Type"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold")
  ) +
  coord_cartesian(clip = "on")  # set to "off" if labels get clipped by facet boundaries


Unique_geneID <- QTL_markers %>% filter(classification=="Unique") %>% select(genome,id)
# write.table(Unique_geneID,"UniqueGeneID.txt",sep="\t",quote=F,col.names=F,row.names=F)





####
#### Chromosome 1 Top QTL Region
####

Chrom1_QTL_pgID <- EL102_Markers[EL102_Markers$chr=="Chromosome1" & EL102_Markers$start>18000000 & EL102_Markers$end < 19600000,1] 
AltGenomeRange <- annot_pgID %>% filter(genome== "KDH49" & chr=="Chromosome1") %>% filter(as.integer(pgID) %in% unlist(as.vector(Chrom1_QTL_pgID)))
Chrom1_QTL_Otholog_table <- annot_pgID %>% filter(chr=="Chromosome1" & ((genome=="EL102" & start>18000000 & end < 19600000) | (genome=="KDH49" & start>17865600 & end < 19325700)))
Chrom1_QTL_count <- Chrom1_QTL_Otholog_table %>% 
  group_by(pgID,genome) %>% 
  summarise(count=n()) %>% 
  spread(genome, count)
Chrom1_QTL_count[is.na(Chrom1_QTL_count)] <- 0
sum(Chrom1_QTL_count$EL102==0 & Chrom1_QTL_count$KDH49 > 0) 




EL102_Markers_QTL <- EL102_Markers %>% filter(chr=="Chromosome1" & start>18000000 & end < 19600000)
KDH49_Markers_QTL <- KDH49_Markers %>% filter(chr=="Chromosome1" & start>17865600 & end < 19325700)
QTL_markers <- rbind(EL102_Markers_QTL,KDH49_Markers_QTL)

## Make plot with both genomes
QTL_markers$genome <- factor(QTL_markers$genome,levels=c("KDH49","EL102"))

ggplot(QTL_markers, aes(y = chr, x = start, color = classification)) +
  geom_point() +
  ggrepel::geom_text_repel(
    aes(label = pgID),
    size = 3,
    max.overlaps = 100,         # increase if you want more labels shown
    min.segment.length = 0,     # always draw a segment
    segment.size = 0.2,
    box.padding = 0.15,
    point.padding = 0.15,
    nudge_x = 0,                # usually not needed with repel
    show.legend = FALSE
  ) +
  theme_bw() +
  facet_wrap(~genome, scales = "free", ncol = 1) +
  scale_color_manual(values = c("Unique"="red","Duplicate"="yellow","Homologous"="black")) +
  theme(legend.position = "top") +
  labs(
    title = "Gene Content Variation in Chromosome 1 QTL Region",
    x = "Chromosome",
    y = "Position (bp)",
    color = "Gene Type"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold")
  ) +
  coord_cartesian(clip = "on")  # set to "off" if labels get clipped by facet boundaries


Unique_geneID <- QTL_markers %>% filter(classification=="Unique") %>% select(genome,id)
write.table(Unique_geneID,"UniqueGeneID_QTLChrom1.txt",sep="\t",quote=F,col.names=F,row.names=F)



####
#### Chromosome 2 Top QTL Region
####

Chrom2_QTL_pgID <- EL102_Markers[EL102_Markers$chr=="Chromosome2" & EL102_Markers$start>44400000 & EL102_Markers$end < 44800000,1] 
AltGenomeRange <- annot_pgID %>% filter(genome== "KDH49" & chr=="Chromosome2") %>% filter(as.integer(pgID) %in% unlist(as.vector(Chrom2_QTL_pgID)))
Chrom2_QTL_Otholog_table <- annot_pgID %>% filter(chr=="Chromosome2" & ((genome=="EL102" & start>44400000 & end < 44800000) | (genome=="KDH49" & start>53775700 & end < 54163100)))
Chrom2_QTL_count <- Chrom2_QTL_Otholog_table %>% 
  group_by(pgID,genome) %>% 
  summarise(count=n()) %>% 
  spread(genome, count)
Chrom2_QTL_count[is.na(Chrom2_QTL_count)] <- 0
sum(Chrom2_QTL_count$EL102==0 & Chrom2_QTL_count$KDH49 > 0) 




EL102_Markers_QTL <- EL102_Markers %>% filter(chr=="Chromosome2" & start>44400000 & end < 44800000)
KDH49_Markers_QTL <- KDH49_Markers %>% filter(chr=="Chromosome2" & start>53775700 & end < 54163100)
QTL_markers <- rbind(EL102_Markers_QTL,KDH49_Markers_QTL)

## Make plot with both genomes
QTL_markers$genome <- factor(QTL_markers$genome,levels=c("KDH49","EL102"))

ggplot(QTL_markers, aes(y = chr, x = start, color = classification)) +
  geom_point() +
  ggrepel::geom_text_repel(
    aes(label = pgID),
    size = 3,
    max.overlaps = 100,         # increase if you want more labels shown
    min.segment.length = 0,     # always draw a segment
    segment.size = 0.2,
    box.padding = 0.15,
    point.padding = 0.15,
    nudge_x = 0,                # usually not needed with repel
    show.legend = FALSE
  ) +
  theme_bw() +
  facet_wrap(~genome, scales = "free", ncol = 1) +
  scale_color_manual(values = c("Unique"="red","Duplicate"="yellow","Homologous"="black")) +
  theme(legend.position = "top") +
  labs(
    title = "Gene Content Variation in Chromosome 2 QTL Region",
    x = "Chromosome",
    y = "Position (bp)",
    color = "Gene Type"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold")
  ) +
  coord_cartesian(clip = "on")  # set to "off" if labels get clipped by facet boundaries


Unique_geneID <- QTL_markers %>% filter(classification=="Unique") %>% select(genome,id)
write.table(Unique_geneID,"UniqueGeneID_QTLChrom2.txt",sep="\t",quote=F,col.names=F,row.names=F)


####
#### Chromosome 4 Top QTL Region
####

Chrom4_QTL_pgID <- EL102_Markers[EL102_Markers$chr=="Chromosome4" & EL102_Markers$start>62000000 & EL102_Markers$end < 62400000,1] 
AltGenomeRange <- annot_pgID %>% filter(genome== "KDH49" & chr=="Chromosome4") %>% filter(as.integer(pgID) %in% unlist(as.vector(Chrom4_QTL_pgID)))
Chrom4_QTL_Otholog_table <- annot_pgID %>% filter(chr=="Chromosome4" & ((genome=="EL102" & start>62000000 & end < 62400000) | (genome=="KDH49" & start>74733300 & end < 75107300)))
Chrom4_QTL_count <- Chrom4_QTL_Otholog_table %>% 
  group_by(pgID,genome) %>% 
  summarise(count=n()) %>% 
  spread(genome, count)
Chrom4_QTL_count[is.na(Chrom4_QTL_count)] <- 0
sum(Chrom4_QTL_count$EL102==0 & Chrom4_QTL_count$KDH49 > 0) 




EL102_Markers_QTL <- EL102_Markers %>% filter(chr=="Chromosome4" & start>62000000 & end < 62500000)
KDH49_Markers_QTL <- KDH49_Markers %>% filter(chr=="Chromosome4" & start>74733300 & end < 75207300)
QTL_markers <- rbind(EL102_Markers_QTL,KDH49_Markers_QTL)

## Make plot with both genomes
QTL_markers$genome <- factor(QTL_markers$genome,levels=c("KDH49","EL102"))

ggplot(QTL_markers, aes(y = chr, x = start, color = classification)) +
  geom_point() +
  ggrepel::geom_text_repel(
    aes(label = pgID),
    size = 3,
    max.overlaps = 100,         # increase if you want more labels shown
    min.segment.length = 0,     # always draw a segment
    segment.size = 0.2,
    box.padding = 0.15,
    point.padding = 0.15,
    nudge_x = 0,                # usually not needed with repel
    show.legend = FALSE
  ) +
  theme_bw() +
  facet_wrap(~genome, scales = "free", ncol = 1) +
  scale_color_manual(values = c("Unique"="red","Duplicate"="yellow","Homologous"="black")) +
  theme(legend.position = "top") +
  labs(
    title = "Gene Content Variation in Chromosome 4 QTL Region",
    x = "Chromosome",
    y = "Position (bp)",
    color = "Gene Type"
  ) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title   = element_text(size = 14),
    axis.text    = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold")
  ) +
  coord_cartesian(clip = "on")  # set to "off" if labels get clipped by facet boundaries


Unique_geneID <- QTL_markers %>% filter(classification=="Unique") %>% select(genome,id)
write.table(Unique_geneID,"UniqueGeneID_QTLChrom4.txt",sep="\t",quote=F,col.names=F,row.names=F)


