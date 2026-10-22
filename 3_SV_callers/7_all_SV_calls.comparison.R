library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(ggupset)

for (svcaller in c("sniffles","delly","cuteSV","HitSV")){

infile <- read.table(paste0("allinds.pangene_miniprot_summary.",svcaller,".csv"), header=FALSE, fill = TRUE,
  quote = "",
  comment.char = "", sep="\t",
  col.names = c("Gene", "Ind", "Pangene_count", "Miniprot_presence", "Contig", "Start", "Stop", "AS", "Query_coverage", "Frameshifts", "In_frame_stops", "TOGA_transcript", "Neighbor_contig_check", "NA_neighbor_check", "Genes_per_walk", "Lone_gene_check", "Same_contig_check", "Max_TOGA_transcript", "Transcript_check", "Deletion_genotype", "Species"))

infile$Hap <- infile$Ind
infile$Ind <- gsub("\\.[12]$", "", infile$Ind)

######### If the species is not read well because of an empty genotype, change the species accordingly and the genotype to 0/0

idx <- infile$Deletion_genotype %in% c("AI", "AC", "AW")
infile$Species[idx] <- infile$Deletion_genotype[idx]
infile$Deletion_genotype[idx] <- "0/0"

## Merge with CDS length and overlap

CDS_length <- read.table(paste0("allinds.CDS.",svcaller,".csv"), header=FALSE, fill = TRUE, quote = "", comment.char = "", sep="\t")
colnames(CDS_length)<- c("Gene", "Total_CDS_length", "Ind")

### Remember to split rows with multiple genes into different rows 

CDS_length <- CDS_length %>%
  separate_rows(Gene, sep = ",gene-") %>%
  mutate(
    Gene = ifelse(
      row_number() > 1 & !startsWith(Gene, "gene-"),
      paste0("gene-", Gene),

      Gene
    ),
    Gene = gsub("gene-", "", Gene)
  )

overlap_length <- read.table(paste0("allinds.CDS_del_overlap.",svcaller,".csv"), header=FALSE, fill = TRUE, quote = "", comment.char = "", sep="\t")
colnames(overlap_length)<- c("Gene", "Overlap_length", "Ind")

### Remember to split rows with multiple genes into different rows 

overlap_length <- overlap_length %>%
  separate_rows(Gene, sep = ",gene-") %>%
  mutate(
    Gene = ifelse(
      row_number() > 1 & !startsWith(Gene, "gene-"),
      paste0("gene-", Gene),
      Gene
    ),
    Gene = gsub("gene-", "", Gene)
  )

infile2 <- merge(merge(infile, CDS_length, by=c("Gene", "Ind"), all.x=TRUE), overlap_length, by=c("Gene", "Ind"), all.x=TRUE)

infile2$P_overlap <- infile2$Overlap_length/infile2$Total_CDS_length

### Retain only 1 gene per haplotype - the one with the best overlap length and, secondly, best deletion genotype ("1/1" taking priority)

infile2 <- infile2 %>%
  mutate(
    deletion_rank = case_when(
      Deletion_genotype == "1/1" ~ 4,
      Deletion_genotype == "0/1" ~ 3,
      Deletion_genotype == "0/0" ~ 2,
      is.na(Deletion_genotype) | Deletion_genotype == "" ~ 1,
      TRUE ~ 0
    )
  ) %>%
  group_by(Gene, Hap) %>%
  arrange(
    desc(Overlap_length),
    desc(deletion_rank),
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup() %>%
  select(-deletion_rank)

## Change deletion genotype to 0 in those genes where proportion of deletion >10% of CDSs' length

infile2[!is.na(infile2$P_overlap) & infile2$P_overlap < 0.1,]$Deletion_genotype <- "0/0"

if (any(infile2$Deletion_genotype =="", na.rm = TRUE)) {
  infile[infile2$Deletion_genotype == "", "Deletion_genotype"] <- "0/0"
}                                                                      ## Also empty deletion genotypes

### And finally rename dataframe so that we have a dataframe for each SVcaller, while also renaming the Deletion_genotype column

colnames(infile2) <- c("Gene", "Ind", "Pangene_count", "Miniprot_presence", "Contig", "Start", "Stop", "AS", "Query_coverage", "Frameshifts", "In_frame_stops", "TOGA_transcript", "Neighbor_contig_check", "NA_neighbor_check", "Genes_per_walk", "Lone_gene_check", "Same_contig_check", "Max_TOGA_transcript", "Transcript_check", paste0("Deletion_genotype_", svcaller), "Species", "Hap", "Total_CDS_length", paste0("Overlap_length_", svcaller), paste0("P_overlap_", svcaller))

assign(paste0("infile_", svcaller), infile2)

}

#####################################################################################

## Read the pangene CNV Rtable

CNV_file <- read.table("sj-94-a2.pg.CNV.Rtable", header=TRUE, sep="\t",fill = TRUE,
  quote = "",
  comment.char = "")

CNV_file <- pivot_longer(CNV_file, cols = -c(Gene), names_to="Hap", values_to="Gene_count_hap")
CNV_file$Ind <- gsub("\\.[12]$", "", CNV_file$Hap)
CNV_file <- CNV_file %>% group_by(Gene,Ind) %>% mutate(Gene_count_ind = sum(Gene_count_hap)) %>% ungroup()

## Get only single copy genes

CNV_file <- CNV_file[CNV_file$Gene_count_ind%in%c(0,1,2),]
CNV_file$Deletion_genotype_pangene <- ""
CNV_file[CNV_file$Gene_count_ind==0,]$Deletion_genotype_pangene <- "1/1"
CNV_file[CNV_file$Gene_count_ind==1,]$Deletion_genotype_pangene <- "0/1"
CNV_file[CNV_file$Gene_count_ind==2,]$Deletion_genotype_pangene <- "0/0"

####################################################################################

## First merge the SVcaller files by all common columns

#infile_sniffles <- subset(infile_sniffles, select=-c(Contig, Start,Stop,AS, TOGA_transcript, Max_TOGA_transcript, Total_CDS_length, Genes_per_walk))
#infile_delly <- subset(infile_delly, select=-c(Contig, Start,Stop,AS, TOGA_transcript, Max_TOGA_transcript, Total_CDS_length, Genes_per_walk))
#infile_cuteSV <- subset(infile_cuteSV, select=-c(Contig, Start,Stop,AS, TOGA_transcript, Max_TOGA_transcript, Total_CDS_length, Genes_per_walk))
#infile_HitSV <- subset(infile_HitSV, select=-c(Contig, Start,Stop,AS, TOGA_transcript, Max_TOGA_transcript, Total_CDS_length, Genes_per_walk))

infile_sniffles <- subset(infile_sniffles, select=c("Gene", "Ind", "Hap", "Species", "Deletion_genotype_sniffles"))
infile_delly <- subset(infile_delly, select=c("Gene", "Ind", "Hap", "Species", "Deletion_genotype_delly"))
infile_cuteSV <- subset(infile_cuteSV, select=c("Gene", "Ind", "Hap", "Species", "Deletion_genotype_cuteSV"))
infile_HitSV <- subset(infile_HitSV, select=c("Gene", "Ind", "Hap", "Species", "Deletion_genotype_HitSV"))


svcallers1 <- merge(infile_sniffles, infile_delly, by=intersect(names(infile_sniffles), names(infile_delly)), all=TRUE)
svcallers2 <- merge(infile_cuteSV, infile_HitSV, by=intersect(names(infile_cuteSV), names(infile_HitSV)), all=TRUE)

svcallers <- merge(svcallers1, svcallers2, by=intersect(names(svcallers1), names(svcallers2)), all=TRUE)

## Merge files

plotfile <- merge(CNV_file, svcallers, by=c("Gene", "Ind", "Hap"))

###############################################################################

plotfile <- read.table("Plotfile.1910.csv", header=TRUE, fill=TRUE,quote = "", comment.char = "", sep=" ")

#### Read raw_read BLAST results

BLAST_raw <- read.table("allinds.BLAST.retained_genes.txt", header=FALSE, fill = TRUE, quote = "\"", comment.char = "", sep=" ")
colnames(BLAST_raw)<- c("Gene", "Ind")

### Remember to split rows with multiple genes into different rows 

BLAST_raw <- BLAST_raw %>%
  separate_rows(Gene, sep = ",gene-") %>%
  mutate( Gene = ifelse(
      row_number() > 1 & !startsWith(Gene, "gene-"),
      paste0("gene-", Gene),Gene), 
      Gene = gsub("gene-", "", Gene)
  )

BLAST_raw <- unique(BLAST_raw)

### Add BLAST genotype column and then add the blast genotype after merging to those "NA" genes, which are 1/1 deleted

BLAST_raw$Deletion_genotype_blast_raw <- "0/0"
plotfile <- merge(plotfile, BLAST_raw, by=c("Gene", "Ind"), all=TRUE)
plotfile[is.na(plotfile$Deletion_genotype_blast_raw),]$Deletion_genotype_blast_raw <- "1/1"

### And remove genes not present in plotfile before (likely due to being multiple-copy genes)

plotfile <- plotfile[!is.na(plotfile$Hap),]

### For analyses, there are a few individuals whose reads we couldn't find:

no_read_inds <- c("AC_1873_10831", "AC_1873_20946", "AI_1363_74565", "AI_1603_79174", "AW_365337", "AW_365338", "AW_366497", "AW_366499", "CY_8788")
'%ni%' <- function(x,y)!('%in%'(x,y))
plotfile <- plotfile[plotfile$Ind %ni% no_read_inds,]

###############################################################################

### Comparison SV_callers_plot

#plotfile <- read.table("Plotfile.1910.csv", header=TRUE, fill=TRUE,quote = "", comment.char = "", sep=" ")

# Keep only whether each method calls 1/1

upset_data <- plotfile %>% mutate(pangene = ifelse(Deletion_genotype_pangene == "1/1", "Hom_del_pangene", NA), 
       sniffles = ifelse(Deletion_genotype_sniffles == "1/1", "Hom_del_sniffles", NA),
       delly = ifelse(Deletion_genotype_delly == "1/1", "Hom_del_delly", NA),
       cuteSV = ifelse(Deletion_genotype_cuteSV == "1/1", "Hom_del_cuteSV", NA),
       HitSV = ifelse(Deletion_genotype_HitSV == "1|1", "Hom_del_HitSV", NA),           # Remember for later: both haplotypes phases present in heterozygous
       BLAST_raw = ifelse(Deletion_genotype_blast_raw == "1/1", "Hom_del_blast_raw", NA))

upset_data <- upset_data %>% unite(col = "Hom_del", c(pangene, sniffles, delly, cuteSV, HitSV,BLAST_raw), sep = "; ",remove = TRUE,   na.rm = TRUE) %>%   # Unite in a column "Hom_del"
                    mutate(hom_del_list=as.list(strsplit(Hom_del, "; ")))

## Remove instances of non-1/1 in all callers

upset_data <- upset_data[upset_data$Hom_del!="",]
upset_data2 <- upset_data[upset_data$Deletion_genotype_pangene=="1/1",]

## Plot

ggplot(upset_data2, aes(x = hom_del_list)) + geom_bar() +
  scale_x_upset(sets = c("Hom_del_blast_raw", "Hom_del_pangene", "Hom_del_sniffles", "Hom_del_delly", "Hom_del_cuteSV", "Hom_del_HitSV")) + 
  labs(
    x = "Genotype callers reporting 1/1",
    y = "Number of lost genes across all individuals",
    title = "Concordance of homozygous deletion genotypes"
  ) + theme_minimal()


########################################################################################################################

## If we supose that SVcallers have problems calling homozygous vs heterozygous SVs

upset_data <- plotfile %>% mutate(pangene = ifelse(Deletion_genotype_pangene == "1/1", "Hom_del_pangene", NA), 
       sniffles = ifelse(Deletion_genotype_sniffles %in% c("0/1", "1/0", "1/1"), "Hom_del_sniffles", NA),
       delly = ifelse(Deletion_genotype_delly %in% c("0/1", "1/0", "1/1"), "Hom_del_delly", NA),
       cuteSV = ifelse(Deletion_genotype_cuteSV %in% c("0/1", "1/0", "1/1"), "Hom_del_cuteSV", NA),
       HitSV = ifelse(Deletion_genotype_HitSV %in% c("0|1", "1|0", "1|1"), "Hom_del_HitSV", NA),           # Remember for later: both haplotypes phases present in heterozygous
       BLAST_raw = ifelse(Deletion_genotype_blast_raw == "1/1", "Hom_del_blast_raw", NA))

upset_data <- upset_data %>% unite(col = "Hom_del", c(pangene, sniffles, delly, cuteSV, HitSV,BLAST_raw), sep = "; ",remove = TRUE,   na.rm = TRUE) %>%   # Unite in a column "Hom_del"
                    mutate(hom_del_list=as.list(strsplit(Hom_del, "; ")))

## Remove instances of non-1/1 in all callers

upset_data <- upset_data[upset_data$Hom_del!="",]
upset_data2 <- upset_data[upset_data$Deletion_genotype_pangene=="1/1",]

## Plot

ggplot(upset_data2, aes(x = hom_del_list)) + geom_bar() +
  scale_x_upset(sets = c("Hom_del_blast_raw", "Hom_del_pangene", "Hom_del_sniffles", "Hom_del_delly", "Hom_del_cuteSV", "Hom_del_HitSV")) + 
  labs(
    x = "Genotype callers reporting 1/1",
    y = "Number of lost genes across all individuals",
    title = "Concordance of homozygous deletion genotypes"
  ) + theme_minimal()


########################################################################################################################

### Heterozygous


plotfile <- read.table("Plotfile.1910.csv", header=TRUE, fill=TRUE,quote = "", comment.char = "", sep=" ")

# Keep only whether each method calls 1/1

upset_data <- plotfile %>% mutate(pangene = ifelse(Deletion_genotype_pangene == "0/1", "het_del_pangene", NA), 
       sniffles = ifelse(Deletion_genotype_sniffles == "0/1", "het_del_sniffles", NA),
       delly = ifelse(Deletion_genotype_delly == "0/1", "het_del_delly", NA),
       cuteSV = ifelse(Deletion_genotype_cuteSV == "0/1", "het_del_cuteSV", NA),
       HitSV = ifelse(Deletion_genotype_HitSV %in% c("0|1", "1|0"), "het_del_HitSV", NA))             # Both haplotypes phases present in heterozygous


upset_data <- upset_data %>% unite(col = "het_del", c(pangene, sniffles, delly, cuteSV, HitSV), sep = "; ",remove = TRUE,   na.rm = TRUE) %>%   # Unite in a column "het_del"
                    mutate(het_del_list=as.list(strsplit(het_del, "; ")))

## Remove instances of non-1/1 in all callers

#upset_data <- upset_data[upset_data$Deletion_genotype_pangene=="0/1",]
upset_data <- upset_data[upset_data$het_del!="",]

## Plot

ggplot(upset_data, aes(x = het_del_list)) + geom_bar() +
  scale_x_upset(sets = c("het_del_pangene", "het_del_sniffles", "het_del_delly", "het_del_cuteSV", "het_del_HitSV")) + 
  labs(
    x = "Genotype callers reporting 0/1",
    y = "Number of partially lost genes across all individuals",
    title = "Concordance of heterozygous deletion genotypes"
  ) + theme_minimal()


#########################################################################################

### Just number of SV callers supporting a call

upset_data <- plotfile %>% mutate(pangene = ifelse(Deletion_genotype_pangene == "1/1", "0", NA), 
       sniffles = ifelse(Deletion_genotype_sniffles == "1/1", "1", NA),
       delly = ifelse(Deletion_genotype_delly == "1/1", "1", NA),
       cuteSV = ifelse(Deletion_genotype_cuteSV == "1/1", "1", NA),
       HitSV = ifelse(Deletion_genotype_HitSV == "1|1", "1", NA))    

upset_data <- pivot_longer(upset_data, cols = c(pangene, sniffles, delly, cuteSV, HitSV), names_to="SVcaller", values_to="SVcall_value")

upset_data <- upset_data[upset_data$Deletion_genotype_pangene=="1/1",]
upset_data <- upset_data %>% subset(select=c(Gene, Hap, Species, SVcaller, SVcall_value)) 
upset_data <- unique(upset_data)

upset_data <- upset_data %>%  group_by(Gene, Hap) %>%  mutate(N_callers = sum(as.numeric(SVcall_value), na.rm = TRUE)) %>% ungroup()
upset_data <- upset_data %>% subset(select=-c(SVcaller, SVcall_value))         ### For the next step have to remove SVcaller
upset_data <- unique(upset_data)

upset_data <- upset_data %>% group_by(N_callers, Hap) %>% mutate(N_genes=n()) %>% ungroup()

boxplot_data <- subset(upset_data, select=c(Hap, Species, N_genes, N_callers))
boxplot_data <- unique(boxplot_data)
boxplot_data <- boxplot_data[boxplot_data$Species %in% c("AI", "AW", "AC"),]

ggplot(boxplot_data, aes(y=N_genes, x=as.factor(N_callers), fill=Species)) + geom_boxplot() +
  theme_minimal() + scale_fill_viridis_d()

### At least 1 call


upset_data <- plotfile %>% mutate(pangene = ifelse(Deletion_genotype_pangene == "1/1", "0", NA), 
       sniffles = ifelse(Deletion_genotype_sniffles == "1/1", "1", NA),
       delly = ifelse(Deletion_genotype_delly == "1/1", "1", NA),
       cuteSV = ifelse(Deletion_genotype_cuteSV == "1/1", "1", NA),
       HitSV = ifelse(Deletion_genotype_HitSV == "1|1", "1", NA))    

upset_data <- pivot_longer(upset_data, cols = c(pangene, sniffles, delly, cuteSV, HitSV), names_to="SVcaller", values_to="SVcall_value")

upset_data <- upset_data[upset_data$Deletion_genotype_pangene=="1/1",]
upset_data <- upset_data %>% subset(select=c(Gene, Hap, Species, SVcaller, SVcall_value)) 
upset_data <- unique(upset_data)

upset_data <- upset_data %>%  group_by(Gene, Hap) %>%  mutate(N_callers = sum(as.numeric(SVcall_value), na.rm = TRUE)) %>% ungroup()
upset_data <- upset_data %>% subset(select=-SVcaller)         ### For the next step have to remove SVcaller
upset_data <- unique(upset_data)

upset_data$N_callers_def <- "0"
upset_data[upset_data$N_callers > 0,]$N_callers_def <- ">0"

upset_data <- upset_data[!is.na(upset_data$SVcall_value),] #### Don't know why this error appeared

upset_data <- upset_data %>% group_by(N_callers_def, Hap) %>% mutate(N_genes=n()) %>% ungroup()

boxplot_data <- subset(upset_data, select=c(Hap, Species, N_genes, N_callers_def))
boxplot_data <- unique(boxplot_data)
boxplot_data <- boxplot_data[boxplot_data$Species %in% c("AI", "AW", "AC"),]

ggplot(boxplot_data, aes(y=N_genes, x=as.factor(N_callers_def), fill=Species)) + geom_boxplot() +
  theme_minimal() + scale_fill_viridis_d()







































###########################################################################

## Just for this one analysis
#plotfile <- subset(plotfile, select=c(Gene, Ind, Gene_count_ind, Deletion_genotype, Species))
#plotfile <- unique(plotfile)

plottable <- plotfile %>%
  count(Ind, Gene_count_ind, Deletion_genotype, Species)

e <-ggplot(plottable[plottable$Ind=="AI_1603_79203",], aes(x = Deletion_genotype, y = factor(Gene_count_ind), fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Sniffles inferred deletion genotype",
    y = "Pangene allele count",
    fill = "Number of genes"
  ) +
  theme_minimal()

ggplot(plottable[plottable$Gene_count_ind==1 & plottable$Deletion_genotype=="0/1",], aes(x=Species, y=n, fill=factor(Gene_count_ind))) + geom_boxplot() + theme_bw()

#################################
##### TYPE OF ERROR #############
#################################

pangene_miniprot <- plotfile

pangene_miniprot$Error_type <- "NA"

# Type A
idx <- pangene_miniprot$Gene_count_hap == 0 &
       pangene_miniprot$Same_contig_check == "True" &
       pangene_miniprot$Neighbor_contig_check == "False" &
       pangene_miniprot$NA_neighbor_check == "Neighbor_not_annotated"

pangene_miniprot$Error_type[idx] <- "Type_A"


# Type B
idx <- pangene_miniprot$Gene_count_hap == 0 &
       pangene_miniprot$Same_contig_check == "True" &
       pangene_miniprot$Neighbor_contig_check == "False" &
       pangene_miniprot$NA_neighbor_check == "Neighbor_annotated"

pangene_miniprot$Error_type[idx] <- "Type_B"


# Type C
idx <- pangene_miniprot$Gene_count_hap == 0 &
       pangene_miniprot$Same_contig_check == "True" &
       pangene_miniprot$Neighbor_contig_check == "True" &
       pangene_miniprot$NA_neighbor_check == "Neighbor_not_annotated"

pangene_miniprot$Error_type[idx] <- "Type_C"


# Type D
idx <- pangene_miniprot$Gene_count_hap == 0 &
       pangene_miniprot$Same_contig_check == "True" &
       pangene_miniprot$Neighbor_contig_check == "True" &
       pangene_miniprot$NA_neighbor_check == "Neighbor_annotated"

pangene_miniprot$Error_type[idx] <- "Type_D"


# Type E
idx <- pangene_miniprot$Gene_count_hap == 0 &
       pangene_miniprot$Same_contig_check == "False" &
       pangene_miniprot$Neighbor_contig_check == "False" &
       pangene_miniprot$NA_neighbor_check == "Neighbor_not_annotated"

pangene_miniprot$Error_type[idx] <- "Type_E"


# Type F
idx <- pangene_miniprot$Gene_count_hap == 0 &
       pangene_miniprot$Same_contig_check == "False" &
       pangene_miniprot$Neighbor_contig_check == "False" &
       pangene_miniprot$NA_neighbor_check == "Neighbor_annotated"

pangene_miniprot$Error_type[idx] <- "Type_F"


# Type G
idx <- pangene_miniprot$Gene_count_hap == 0 &
       pangene_miniprot$Same_contig_check == "False" &
       pangene_miniprot$Neighbor_contig_check == "True" &
       pangene_miniprot$NA_neighbor_check == "Neighbor_not_annotated"

pangene_miniprot$Error_type[idx] <- "Type_G"


# Type H
idx <- pangene_miniprot$Gene_count_hap == 0 &
       pangene_miniprot$Same_contig_check == "False" &
       pangene_miniprot$Neighbor_contig_check == "True" &
       pangene_miniprot$NA_neighbor_check == "Neighbor_annotated"

pangene_miniprot$Error_type[idx] <- "Type_H"

plotfile <- pangene_miniprot

#### And then also at Zebra finch annotation#####

#################################################
############# CONTIG LENGHT #####################
#################################################

ubi <- read.table("TOGA_locations.aphWoo1.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE,comment.char = "")
lengths <- read.table("aphWoo.v1.chr_lengths", header = TRUE, sep = ",", stringsAsFactors = FALSE,comment.char = "")

ubi_lengths <- merge(ubi, lengths)

ubi_lengths<-ubi_lengths[ubi_lengths$Chr!="aphWoo1#REF#SUPER_W_RagTag.trimmed",]                             ####HERE WE EXCLUDE CHRZ and CHRW because of funky stuff going on with their SFS
ubi_lengths<-ubi_lengths[ubi_lengths$Chr!="aphWoo1#REF#ScYP8k310HRSCAF43chZ_RagTag",]

plotfile <- merge(plotfile, ubi_lengths)

#plotfile <- plotfile %>%
#  group_by(Gene) %>%
#  mutate(presence_category = round(mean(allele_loss_count/30, na.rm = TRUE),digits = 1)) %>% ungroup()

plotfile <- plotfile %>%
  group_by(Gene) %>%
  mutate(
    Length_category = case_when(
      Length < 200000 ~ "<200Kb",
      Length >= 200000  & Length < 2000000 ~ "200Kb-2Mb",
      Length >= 2000000  & Length < 20000000 ~ "2Mb-20Mb",
      Length >= 20000000   ~ ">20Mb"
    )
  ) %>%
  ungroup()

plottable <- plotfile %>%
  count(Ind, Gene_count_ind, Deletion_genotype, Species, Length_category)

ggplot(plottable[plottable$Gene_count_ind==1 & plottable$Deletion_genotype=="0/1",], aes(x=Species, y=n, fill=factor(Length_category), colour=Deletion_genotype)) + geom_boxplot() + theme_bw()



### Mix in BLAST results

BLAST <- read.table("SJ.allinds.exoncount.prov.bed", header=TRUE, fill = TRUE,
  quote = "",
  comment.char = "", sep="\t")

BLAST$Ind <- BLAST$BLAST_ind

BLAST$Exon_prop <- as.numeric(BLAST$BLAST_exon_count)/as.numeric(BLAST$REF_exon_count)

BLAST$BLAST_filter <- ifelse(
  is.na(BLAST$Exon_prop),
  "NA",
  ifelse(BLAST$Exon_prop > 0.85,
         "Present_in_reads",
         "Absent_in_reads")
)

plotfile <- merge(merge(CNV_file, infile, by=c("Gene", "Ind", "Hap")), BLAST, by=c("TOGA_transcript","Ind"))

## Just for this one analysis
plotfile <- subset(plotfile, select=c(Gene, Ind, Gene_count_ind, Deletion_genotype, Species))
plotfile <- unique(plotfile)

## Get only single copies max

plotfile <- plotfile[plotfile$Gene_count_ind%in%c(0,1,2),]
plotfile[plotfile$Deletion_genotype=="",]$Deletion_genotype <- "0/0"

plottable <- plotfile %>%
  count(Ind, Gene_count_ind, Deletion_genotype, Species)

############################################################
### Accuracy and recall ####

plottable <- plotfile %>%
  count(Ind, Gene_count_ind, Deletion_genotype, Species, Length_category, Error_type)

plottable[plottable$Deletion_genotype=="0/0",]$Deletion_genotype <- "2"
plottable[plottable$Deletion_genotype=="0/1",]$Deletion_genotype <- "1"
plottable[plottable$Deletion_genotype=="1/1",]$Deletion_genotype <- "0"

plottable <- plottable[plottable$Species !="CY",]

## 1. PER SPECIES

## Let's calculate the false discovery proportion per species

FDR <- plottable %>% 
        mutate(concordant = Gene_count_ind == as.integer(Deletion_genotype)) %>%
        group_by(Ind, Species, Deletion_genotype) %>%
        summarise(concordance     = sum(n[concordant]) / sum(n),
                  false_discovery = sum(n[!concordant]) / sum(n))

ggplot(FDR, aes(x=Species, y=as.numeric(false_discovery), fill=as.factor(Deletion_genotype))) + geom_boxplot() + theme_bw() + scale_fill_viridis_d()

# Given that concordance might vary by category, maybe use Cohen's weighted kappa

## 2. PER LENGTH CATEGORY

FDR <- plottable %>% 
        mutate(concordant = Gene_count_ind == as.integer(Deletion_genotype)) %>%
        group_by(Ind, Length_category, Gene_count_ind) %>%
        summarise(concordance     = sum(n[concordant]) / sum(n),
                  false_discovery = sum(n[!concordant]) / sum(n))

FDR$Length_category <- factor(FDR$Length_category, levels = c(">20Mb", "2Mb-20Mb", "200Kb-2Mb", "<200Kb"))

ggplot(FDR, aes(x=Length_category, y=as.numeric(concordance), fill=as.factor(Gene_count_ind))) + geom_boxplot() + theme_bw() + scale_fill_viridis_d() + labs(x="Reference scaffold length", y="Concordance between Sniffles2 and pangene genotypes", fill="nº pangene copies per individual")

## 3. PER TYPE OF ERROR

FDR <- plottable %>% 
        mutate(concordant = Gene_count_ind == as.integer(Deletion_genotype)) %>%
        group_by(Ind, Error_type, Gene_count_ind) %>%
        summarise(concordance     = sum(n[concordant]) / sum(n),
                  false_discovery = sum(n[!concordant]) / sum(n))

ggplot(FDR, aes(x=Error_type, y=as.numeric(false_discovery), fill=as.factor(Gene_count_ind))) + geom_boxplot() + theme_bw() + scale_fill_viridis_d()

## 4. PER PANGENE GENOTYPE

FDR <- plottable %>% 
        mutate(concordant = Gene_count_ind == as.integer(Deletion_genotype)) %>%
        group_by(Ind, Deletion_genotype, Gene_count_ind) %>%
        summarise(concordance     = sum(n[concordant]) / sum(n),
                  false_discovery = sum(n[!concordant]) / sum(n))

ggplot(FDR, aes(x=as.factor(Gene_count_ind), y=as.numeric(false_discovery), fill=as.factor(Deletion_genotype))) + geom_boxplot() + theme_bw() + scale_fill_viridis_d()

### 5. Prova plot quadrat
plotfile$Block_unannotated_deletion <- "NA"
plotfile[plotfile$Error_type=="Type_A",]$Block_unannotated_deletion <- "No"
plotfile[plotfile$Error_type=="Type_B",]$Block_unannotated_deletion <- "No"
plotfile[plotfile$Error_type=="Type_C",]$Block_unannotated_deletion <- "No"
plotfile[plotfile$Error_type=="Type_D",]$Block_unannotated_deletion <- "No"
plotfile[plotfile$Error_type=="Type_E",]$Block_unannotated_deletion <- "Yes"
plotfile[plotfile$Error_type=="Type_F",]$Block_unannotated_deletion <- "Yes"
plotfile[plotfile$Error_type=="Type_G",]$Block_unannotated_deletion <- "Yes"
plotfile[plotfile$Error_type=="Type_H",]$Block_unannotated_deletion <- "Yes"

plottable <- plotfile %>%
  count(Ind, Gene_count_ind, Deletion_genotype, Species, Block_unannotated_deletion)

ggplot(plottable[plottable$Ind=="AI_1603_79203" & plottable$Block_unannotated_deletion=="No" & plottable$Gene_count_ind!=2,], aes(x = Deletion_genotype, y = factor(Gene_count_ind), fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(
    x = "Sniffles inferred deletion genotype",
    y = "Pangene allele count",
    fill = "Number of genes"
  ) +
  theme_minimal()

plottable[plottable$Deletion_genotype=="0/0",]$Deletion_genotype <- "2"
plottable[plottable$Deletion_genotype=="0/1",]$Deletion_genotype <- "1"
plottable[plottable$Deletion_genotype=="1/1",]$Deletion_genotype <- "0"

plottable <- plottable[plottable$Species !="CY",]

FDR <- plottable %>% 
        mutate(concordant = Gene_count_ind == as.integer(Deletion_genotype)) %>%
        group_by(Ind, Gene_count_ind, Block_unannotated_deletion) %>%
        summarise(concordance     = sum(n[concordant]) / sum(n),
                  false_discovery = sum(n[!concordant]) / sum(n))

ggplot(FDR, aes(x=as.factor(Block_unannotated_deletion), y=as.numeric(false_discovery), fill=as.factor(Gene_count_ind))) + geom_boxplot() + theme_bw() + scale_fill_viridis_d()
