rm(list=ls())
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(microbiome)
theme_set(theme_bw())

# ############################################################### (0) FUNCTIONS ##################################
PS.FILTER.F <- function(PS, DETECTION){
  for(n in 1:nsamples(PS)){otu_table(PS)[,n] <- otu_table(PS)[,n]/sample_data(PS)$reads.bacteria.bracken[n]} # divide each sample by the total number of bacterial reads
  PS <- transform_sample_counts(PS, function(OTU) ifelse(OTU < DETECTION/100, 0, OTU)) # for each sample separately, replace taxa with values with abundance below DETECTION in that sample by 0
  for(n in 1:nsamples(PS)){otu_table(PS)[,n] <- otu_table(PS)[,n]*sample_data(PS)$reads.bacteria.bracken[n]} # multiply each sample again by the total number of bacterial reads to return to raw abundance values
  PS <- prune_taxa(taxa_sums(PS) > 0, PS) # drop taxa that are now empty due to filtering
}
############################################################### (0) SETTINGS #######################################
setwd("C:/Users/cvandijck/OneDrive - ITG/ITM-PhD/STUDIES/PReGo/PReGoMetagenomics/Bioinformatics_Binning_AmrPlusPlus")

# load phyloseq object
BRACKEN <- readRDS("bracken.RDS"); BRACKEN
MEGARES <- readRDS("megares.RDS"); MEGARES


############################################################### (2) CUTOFF MINIMUM READ COUNT #############
ps <- BRACKEN
sample_data(ps)$reads.arg.cat <- ifelse(sample_data(ps)$reads.arg.all >0, "At least 1 ARG detected", "No ARG detected")
sample_data(ps)$reads.arg.cat <- factor(sample_data(ps)$reads.arg.cat, levels=c("At least 1 ARG detected", "No ARG detected"))

# Suppl_Fig1
ggplot(data=meta(ps), aes(x=log10(reads.nonhost+1), y=log10(1+1000000*reads.arg.all/reads.bacteria.bracken))) + 
  geom_point(aes(color=reads.arg.cat, size=reads.nonhost))+
  scale_shape_manual(values=c(1,19))+
  geom_vline(xintercept=log10(5500), color="grey")+
  annotate(geom="text", label=c("< 5,500 non-host reads", "\u2265 5,500 non-host reads"), x=c(log10(4500), log10(6500)), y=c(3.7,3.7), angle=90, col="grey", size=3) +
  labs(title="Supplementary Figure 1: \nNormalised read count of antimicrobial resistance genes \n(ARGs, in reads per million) versus number of non-host reads.", 
       subtitle="Every dot represents one sample. Dot size corresponds to number of non-host reads in the sample.",
       x="log10(Non-host read count + 1)", y="log10(Normalized ARG read count + 1)") +
  theme(legend.position = "bottom", legend.title=element_blank(), legend.spacing = unit(0,"pt"), 
        legend.box.background = element_rect(linetype = 2, size = 0.5, colour = "grey"),
        legend.box="vertical") +
  ylim(-0.1, 4.1)

sample_data(ps)$sufficientreads.cat <- ifelse(sample_data(ps)$reads.nonhost < 5500, "< 5,500 non-host reads", "\u2265 5,500 non-host reads")

# Suppl_Fig2
PS <- ps
PS <- aggregate_rare(PS, level="Genus", detection=0, prevalence=0)
PS <- transform(PS, transform="clr")
ord <- ordinate(PS, method="PCoA", distance = "euclidean")
plot_ordination(PS, ord, color = "sufficientreads.cat", axes=c(1,2)) + 
  geom_point(aes(size=reads.nonhost)) +
  scale_shape_manual(values=c(1,19))+
  scale_colour_manual(values=c(3,6))+
  stat_ellipse(aes(group=sufficientreads.cat)) +
  labs(title="Supplementary Figure 2: \nPrincipal component plot of all samples, coloured according \nto a cutoff of 5,500 non-host reads.", subtitle="PCA based on Genus-level Euclidean distances of centered log-ratio transformed taxonomic data. \nEvery dot represents one sample. Dot size corresponds to number of non-host reads in the sample.")+
  theme(legend.position = "bottom", legend.title=element_blank(), legend.spacing = unit(0,"pt"), 
        legend.box.background = element_rect(linetype = 2, size = 0.5, colour = "grey"),
        legend.box="vertical")

# Read statistics before discarding samples with low sequencing depth
df <- meta(ps)
df$sufficientreads.b <- ifelse(df$reads.nonhost >= 5500, 1, 0)
df <- df %>% mutate_at(.vars=vars(ends_with(".b")), .funs=funs(factor(., levels=c(0,1), labels=c("No", "Yes"), ordered=FALSE)))
summary(df$reads.raw)
summary(df$reads.trimmed)
summary(df$reads.nonhost)
summary(df$reads.arg.all)
df$arg.rpm <- 1000000*df$reads.arg.all/df$reads.bacteria.bracken
summary(df$arg.rpm)
cor.test(df$reads.arg.all, df$arg.rpm, method="spearman")

# Read statistics after discarding samples with low sequencing depth
table(df$sufficientreads.b)
df <- subset(df, sufficientreads.b=="Yes")
summary(df$reads.raw)
summary(df$reads.trimmed)
summary(df$reads.nonhost)
summary(df$reads.arg.all)
summary(df$arg.rpm)
cor.test(df$reads.arg.all, df$arg.rpm, method="spearman")


############################################################### (3) DECONTAMINATION #############
#---------------------------------------------------------------(3.1) BRACKEN ----------------
ps <- BRACKEN
table1::table1(~DNA.preparation.batch, data=meta(ps))
ps <- subset_samples(ps, reads.nonhost >= 5500) # drop samples with less than 5500 nonhost reads
ps.f <- PS.FILTER.F(PS=ps, DETECTION=0.1); ps.f # filter low-abundant taxa

PS <- ps.f
batchlist <- list()
for(batch in sort(unique(meta(PS)$DNA.preparation.batch))){
  batchlist[[paste0("batch", batch)]] <- subset_samples(PS, DNA.preparation.batch==batch) %>%
    aggregate_rare(detection=0, prevalence=0, level="Species") %>%
    # core(detection=1/100, prevalence=1/4, include.lowest=FALSE) %>%
    prevalence()
  batchlist[[paste0("not_batch", batch)]] <- subset_samples(PS, !DNA.preparation.batch==batch) %>%
    aggregate_rare(detection=0, prevalence=0, level="Species") %>%
    # core(detection=1/100, prevalence=1/4, include.lowest=FALSE) %>%
    prevalence()
}

df <- data.frame(unique(data.frame(tax_table(PS))[,"Species"]))
colnames(df) <- "tax"
for(batch in names(batchlist)){
  tmp.df <- data.frame(batchlist[[batch]])
  colnames(tmp.df) <- batch
  tmp.df$tax <- rownames(tmp.df)
  df <- merge(df, tmp.df, by="tax", all=TRUE)
}
df[is.na(df)] <- 0 # replace NA by 0
df <- df[-which(rowSums(df[,2:ncol(df)])==0),] # drop empty rows

library(ggrepel)
p02 <- ggplot(data=df, aes(x=batch02, y=not_batch02, label=tax)) +
  geom_point(color="dark grey") + 
  geom_point(data=subset(df, batch02 >= 0.25 & not_batch02 < 0.25), color="black") +
  geom_text_repel(data=subset(df, batch02 >= 0.25 & not_batch02 < 0.25), size=2, color="red", max.overlaps=100) +
  labs(title = "Likely contaminants in batch 02", x="Prevalence in batch 02", y="Prevalence in the remaining batches")

p03 <- ggplot(data=df, aes(x=batch03, y=not_batch03, label=tax)) +
  geom_point() + 
  geom_point(data=subset(df, batch03 >= 0.25 & not_batch03 < 0.25), color="red") + 
  geom_text_repel(data=subset(df, batch03 >= 0.25 & not_batch03 < 0.25), size=2, color="red", max.overlaps=100) +  
  labs(title = "Likely contaminants in batch 03")

p04 <- ggplot(data=df, aes(x=batch04, y=not_batch04, label=tax)) +
  geom_point() + 
  geom_point(data=subset(df, batch04 >= 0.25 & not_batch04 < 0.25), color="red") + 
  geom_text_repel(data=subset(df, batch04 >= 0.25 & not_batch04 < 0.25), size=2, color="red", max.overlaps=100) +  
  labs(title = "Likely contaminants in batch 04")

p05 <- ggplot(data=df, aes(x=batch05, y=not_batch05, label=tax)) +
  geom_point() + 
  geom_point(data=subset(df, batch05 >= 0.25 & not_batch05 < 0.25), color="red") + 
  geom_text_repel(data=subset(df, batch05 >= 0.25 & not_batch05 < 0.25), size=2, color="red", max.overlaps=100) +  
  labs(title = "Likely contaminants in batch 05")

p06 <- ggplot(data=df, aes(x=batch06, y=not_batch06, label=tax)) +
  geom_point() + 
  geom_point(data=subset(df, batch07 >= 0.25 & not_batch07 < 0.25), color="red") + 
  geom_text_repel(data=subset(df, batch07 >= 0.25 & not_batch07 < 0.25), size=2, color="red", max.overlaps=100) +  
  labs(title = "Likely contaminants in batch 06")

p07 <- ggplot(data=df, aes(x=batch07, y=not_batch07, label=tax)) +
  geom_point() + 
  geom_point(data=subset(df, batch07 >= 0.25 & not_batch07 < 0.25), color="red") + 
  geom_text_repel(data=subset(df, batch07 >= 0.25 & not_batch07 < 0.25), size=2, color="red", max.overlaps=100) +
  labs(title = "Likely contaminants in batch 07")

p11 <- ggplot(data=df, aes(x=batch11, y=not_batch11, label=tax)) +
  geom_point() + 
  geom_point(data=subset(df, batch11 >= 0.25 & not_batch11 < 0.25), color="red") + 
  geom_text_repel(data=subset(df, batch11 >= 0.25 & not_batch11 < 0.25), size=2, color="red", max.overlaps=100) +
  labs(title = "Likely contaminants in batch 11")

p12 <- ggplot(data=df, aes(x=batch12, y=not_batch12, label=tax)) +
  geom_point() + 
  geom_point(data=subset(df, batch12 >= 0.25 & not_batch12 < 0.25), color="red") + 
  geom_text_repel(data=subset(df, batch12 >= 0.25 & not_batch12 < 0.25), size=2, color="red", max.overlaps=100) +  
  labs(title = "Likely contaminants in batch 12")
  
plots.list = list(p02, p03, p04, p05, p06, p07, 
                  # p08, p09, p10, 
                  p11, p12)
li = structure(plots.list, class = c("gglist", "ggplot"))
print.gglist = function(x, ...) plyr::l_ply(x, print, ...)
# ggsave(li, file = "Suppl_Fig3.pdf", device="pdf", dpi=300, height=25, width=25, unit="cm")
  
PS <- ps.f
plots.list2 <- list()
for(batch in as.character(sort(unique(meta(PS)$DNA.preparation.batch)))){
  PSCOR <- PS %>% 
    subset_samples(DNA.preparation.batch==batch) %>%
    aggregate_rare(detection=0, prevalence=0.25, level="Species")
  PSCOR <- prune_taxa(taxa_sums(PSCOR) > 0, PSCOR)
  taxa_names(PSCOR) <- data.frame(tax_table(PSCOR))[,"Species"]
  OTU <- t(abundances(PSCOR))
  COR <- cor(OTU, method="spearman")
  CONF <- ggcorrplot::cor_pmat(OTU, sig.level=0.05)
  plots.list2[[batch]] <- ggcorrplot ::ggcorrplot(COR, hc.order = TRUE, hc.method="ward.D", type = "full", 
                                                  # p.mat = CONF, insig="blank",
                                                  lab=F, tl.cex=5, ggtheme=theme_bw) +
    labs(title=paste("Correlation plot batch", batch, "(n=", nsamples(PSCOR), "samples).")) + 
    scale_fill_gradient2(limit = c(-1,1), low = "blue", high =  "red", mid = "white", midpoint = 0)
}

li = structure(plots.list2, class = c("gglist", "ggplot"))
print.gglist = function(x, ...) plyr::l_ply(x, print, ...)
# ggsave(li, file = "Suppl_FigX.pdf", device="pdf", dpi=300, height=25, width=25, unit="cm")

# Suppl_Fig5
PS <- ps.f
PSCOR <- aggregate_rare(PS, detection=0, prevalence=0.1, level="Species")
PSCOR <- prune_taxa(taxa_sums(PSCOR) > 0, PSCOR)
taxa_names(PSCOR) <- data.frame(tax_table(PSCOR))[,"Species"]
OTU <- t(abundances(PSCOR))
COR <- cor(OTU, method="spearman")
CONF <- ggcorrplot::cor_pmat(OTU, sig.level=0.05)
ggcorrplot ::ggcorrplot(COR, hc.order = TRUE, hc.method="ward.D", type = "full", 
                                                # p.mat = CONF, insig="blank",
                                                lab=F, tl.cex=5, ggtheme=theme_bw) +
  labs(title=paste("Supplementary Figure 5: Spearman correlation of bacterial species (all batches combined).")) + 
  scale_fill_gradient2(limit = c(-1,1), low = "blue", high =  "red", mid = "white", midpoint = 0)
# ggsave(last_plot(), file = "Suppl_Fig5.jpg", device="jpg", dpi=300, height=25, width=25, unit="cm")

# which species to drop?
batch02 <-subset(df, select=c("tax", "batch02", "not_batch02"))
colnames(batch02) <- c("tax", "inbatch", "inotherbatches")
batch02$batch <- "batch02"
batch03 <-subset(df, select=c("tax", "batch03", "not_batch03"))
colnames(batch03) <- c("tax", "inbatch", "inotherbatches")
batch03$batch <- "batch03"
batch04 <-subset(df, select=c("tax", "batch04", "not_batch04"))
colnames(batch04) <- c("tax", "inbatch", "inotherbatches")
batch04$batch <- "batch04"
batch05 <-subset(df, select=c("tax", "batch05", "not_batch05"))
colnames(batch05) <- c("tax", "inbatch", "inotherbatches")
batch05$batch <- "batch05"
batch06 <-subset(df, select=c("tax", "batch06", "not_batch06"))
colnames(batch06) <- c("tax", "inbatch", "inotherbatches")
batch06$batch <- "batch06"
batch07 <-subset(df, select=c("tax", "batch07", "not_batch07"))
colnames(batch07) <- c("tax", "inbatch", "inotherbatches")
batch07$batch <- "batch07"
batch11 <-subset(df, select=c("tax", "batch11", "not_batch11"))
colnames(batch11) <- c("tax", "inbatch", "inotherbatches")
batch11$batch <- "batch11"
batch12 <-subset(df, select=c("tax", "batch12", "not_batch12"))
colnames(batch12) <- c("tax", "inbatch", "inotherbatches")
batch12$batch <- "batch12"

decontam <- rbind(batch02, batch03, batch04, batch05, batch06, batch07, batch11, batch12)
contam <- subset(decontam, inbatch>=0.25 & inotherbatches <0.25, select=c("batch", "tax"))
# write.csv2(contam, "contaminants.csv", row.names=F)  => manually add a column "contam" and indicate contaminants with 1, otherwise 0

contam <- read.csv2("Final/contaminants.csv")

a <- data.frame(pivot_wider(data=contam, names_from = batch, values_from = contam))
a[is.na(a)] <- 0
rownames(a) <- a$tax
a <- a[,-1]
a$total <- rowSums(a)
a <- subset(a, total>0)
a <- arrange(a, -total)

# Suppl_Table1
sjPlot::tab_df(a, show.rownames=T)

decontam <- decontam %>% full_join(contam, by=c("batch", "tax"))
decontam[is.na(decontam)] <- 0

# Suppl_Fig3
ggplot(data=decontam, aes(x=inbatch, y=inotherbatches, label=tax)) +
  geom_point(colour="dark grey") + 
  geom_point(data=subset(decontam, inbatch >= 0.25 & inotherbatches < 0.25), color="black") + 
  geom_point(data=subset(decontam, contam==1), color="red") + 
  geom_text_repel(data=subset(decontam, contam==1), size=2, color="red", max.overlaps=100) +  
  labs(title = "Supplementary Figure 3: Likely contaminant bacterial species, per analysis batch",
       subtitle="Dots represent bacterial species with a minimum abundance of 0.1%. Species with a prevalence \u2265 25% in one batch and < 25% in the remaining batches are \nindicated in black and red. Annotated species (red dots) were identified as likely contaminants based on prevalence and correlation analysis.",
       x="Prevalence in batch", y="Prevalence in remaining batches") +
  facet_wrap(~batch)
# ggsave(last_plot(), file = "2022-05-11_Suppl_Figure3_Bracken_contaminants_prevplots.jpg", device="jpg", dpi=300, height=30, width=30, unit="cm")


## exclude contaminants from their respective batches
PS <- ps.f
batchlist <- list()
for(batch in sort(unique(meta(PS)$DNA.preparation.batch))){
  bat <- paste0("batch", batch)
  dropTaxa <- subset(contam, batch==bat & contam==1)
  PS_sub <- subset_samples(PS, DNA.preparation.batch==batch) # subset batch of interest
  PS_sub <- prune_taxa(taxa_sums(PS_sub) > 0, PS_sub) # drop empty rows in count matrix
  allTaxa <- data.frame(tax_table(PS_sub))$Species
  keepTaxa <- allTaxa[!(allTaxa %in% dropTaxa$tax)]
  batchlist[[bat]] <- subset_taxa(PS_sub, Species %in% keepTaxa)
}

ps_decont <- merge_phyloseq(batchlist$batch02, batchlist$batch03, batchlist$batch04, batchlist$batch05, batchlist$batch06, batchlist$batch07, batchlist$batch11, batchlist$batch12)

# saveRDS(ps_decont, "bracken_decont.RDS")

#---------------------------------------------------------------(3.2) MEGARES ----------------
ps <- MEGARES
ps <- subset_samples(ps, reads.nonhost >= 5500) # exclude samples with less than 5500 nonhost reads
ps <- prune_taxa(taxa_sums(ps) > 0, ps); ps

PS <- ps
batchlist <- list()
for(batch in meta(PS)$DNA.preparation.batch){
  batchlist[[paste0("batch", batch)]] <- subset_samples(PS, DNA.preparation.batch==batch)%>%
    aggregate_rare(detection=0, prevalence=0, level="group") %>%
    prevalence()
  batchlist[[paste0("not_batch", batch)]] <- subset_samples(PS, !DNA.preparation.batch==batch) %>%
    aggregate_rare(detection=0, prevalence=0, level="group") %>%
    prevalence()
}

df <- data.frame(unique(data.frame(tax_table(PS))[,"group"]))
colnames(df) <- "tax"
for(batch in names(batchlist)){
  tmp.df <- data.frame(batchlist[[batch]])
  colnames(tmp.df) <- batch
  tmp.df$tax <- rownames(tmp.df)
  df <- merge(df, tmp.df, by="tax", all=TRUE)
}
df[is.na(df)] <- 0 # replace NA by 0
df <- df[-which(rowSums(df[,2:ncol(df)])==0),] # drop empty rows

library(ggrepel)
p02 <- ggplot(data=df, aes(x=batch02, y=not_batch02, label=tax)) +
  geom_point() + 
  geom_text_repel(data=subset(df, batch02 > 0.1 & not_batch02 == 0), size=2.5, max.overlaps=100) +
  labs(title = "Likely contaminants in batch 02")

p03 <- ggplot(data=df, aes(x=batch03, y=not_batch03, label=tax)) +
  geom_point() + 
  geom_text_repel(data=subset(df, batch03 > 0.1 & not_batch03 == 0), size=2.5, max.overlaps=100) +
  labs(title = "Likely contaminants in batch 03")

p04 <- ggplot(data=df, aes(x=batch04, y=not_batch04, label=tax)) +
  geom_point() + 
  geom_text_repel(data=subset(df, batch04 > 0.1 & not_batch04 == 0), size=2.5, max.overlaps=100) +
  labs(title = "Likely contaminants in batch 04")

p05 <- ggplot(data=df, aes(x=batch05, y=not_batch05, label=tax)) +
  geom_point() + 
  geom_text_repel(data=subset(df, batch05 > 0.1 & not_batch05 == 0), size=2.5, max.overlaps=100) +
  labs(title = "Likely contaminants in batch 05")

p06 <- ggplot(data=df, aes(x=batch06, y=not_batch06, label=tax)) +
  geom_point() + 
  geom_text_repel(data=subset(df, batch06 > 0.1 & not_batch06 == 0), size=2.5, max.overlaps=100) +
  labs(title = "Likely contaminants in batch 06")

p07 <- ggplot(data=df, aes(x=batch07, y=not_batch07, label=tax)) +
  geom_point() + 
  geom_text_repel(data=subset(df, batch07 > 0.1 & not_batch07 == 0), size=2.5, max.overlaps=100) +
  labs(title = "Likely contaminants in batch 07")

p11 <- ggplot(data=df, aes(x=batch11, y=not_batch11, label=tax)) +
  geom_point() + 
  geom_text_repel(data=subset(df, batch11 > 0.1 & not_batch11 == 0), size=2.5, max.overlaps=100) +
  labs(title = "Likely contaminants in batch 11")

p12 <- ggplot(data=df, aes(x=batch12, y=not_batch12, label=tax)) +
  geom_point() + 
  geom_text_repel(data=subset(df, batch12 > 0.1 & not_batch12 == 0), size=2.5, max.overlaps=100) +
  labs(title = "Likely contaminants in batch 12")

plots.list = list(p02, p03, p04, p05, p06, p07, 
                  # p08, p09, p10, 
                  p11, p12)
li = structure(plots.list, class = c("gglist", "ggplot"))
print.gglist = function(x, ...) plyr::l_ply(x, print, ...)
# ggsave(li, file = "Suppl_Fig.pdf", device="pdf", dpi=300, height=25, width=25, unit="cm")


# correlation plot
PS <- ps
OTU <- PS %>%
  aggregate_rare(detection=0, prevalence=0, level="group") %>%
  abundances() %>%
  t()
COR <- cor(OTU, method="spearman")
CONF <- ggcorrplot::cor_pmat(OTU, sig.level=0.05)
ggcorrplot ::ggcorrplot(COR, hc.order = TRUE, type = "full", 
                        # p.mat = CONF, insig="blank", 
                        tl.cex=8, ggtheme=theme_bw) +
  labs(title="Supplementary Figure 6: Spearman correlation of antimicrobial resistance genes \n(all batches combined).") + 
  scale_fill_gradient2(limit = c(-1,1), low = "blue", high =  "red", mid = "white", midpoint = 0)
# ggsave("Suppl_Fig6_Megares_contaminants_corrplot.jpg", plot=last_plot(), device="jpg",dpi=300, height=20, width=20, unit="cm")


# which ARG to drop?
batch02 <-subset(df, select=c("tax", "batch02", "not_batch02"))
colnames(batch02) <- c("tax", "inbatch", "inotherbatches")
batch02$batch <- "batch02"
batch03 <-subset(df, select=c("tax", "batch03", "not_batch03"))
colnames(batch03) <- c("tax", "inbatch", "inotherbatches")
batch03$batch <- "batch03"
batch04 <-subset(df, select=c("tax", "batch04", "not_batch04"))
colnames(batch04) <- c("tax", "inbatch", "inotherbatches")
batch04$batch <- "batch04"
batch05 <-subset(df, select=c("tax", "batch05", "not_batch05"))
colnames(batch05) <- c("tax", "inbatch", "inotherbatches")
batch05$batch <- "batch05"
batch06 <-subset(df, select=c("tax", "batch06", "not_batch06"))
colnames(batch06) <- c("tax", "inbatch", "inotherbatches")
batch06$batch <- "batch06"
batch07 <-subset(df, select=c("tax", "batch07", "not_batch07"))
colnames(batch07) <- c("tax", "inbatch", "inotherbatches")
batch07$batch <- "batch07"
batch11 <-subset(df, select=c("tax", "batch11", "not_batch11"))
colnames(batch11) <- c("tax", "inbatch", "inotherbatches")
batch11$batch <- "batch11"
batch12 <-subset(df, select=c("tax", "batch12", "not_batch12"))
colnames(batch12) <- c("tax", "inbatch", "inotherbatches")
batch12$batch <- "batch12"

decontam <- rbind(batch02, batch03, batch04, batch05, batch06, batch07, batch11, batch12)
dropARG <- c("APH3-DPRIME", "AAC6-PRIME", "L1", "IRI")
decontam$contam <- ifelse(decontam$batch=="batch11" & decontam$tax %in% dropARG,1,0)

# Suppl Fig 4
ggplot(data=decontam, aes(x=inbatch, y=inotherbatches, label=tax)) +
  geom_point(colour="dark grey") + 
  geom_point(data=subset(decontam, inbatch >= 0.1 & inotherbatches < 0.1), color="black") + 
  geom_point(data=subset(decontam, contam==1), color="red") + 
  geom_text_repel(data=subset(decontam, contam==1), size=2, color="red", max.overlaps=100) +  
  labs(title = "Supplementary Figure 4: Likely contaminant antimicrobial resistance genes (ARGs), per analysis batch",
       subtitle="Dots represent ARGs. ARGs with a prevalence \u2265 10% in one batch and < 10% in the remaining batches are \nindicated in black and red. Annotated ARGs (red dots) were identified as likely contaminants based on prevalence and correlation analysis.",
       x="Prevalence in batch", y="Prevalence in remaining batches") +
  facet_wrap(~batch)
# ggsave(last_plot(), file = "Suppl_Fig4_Megares_contaminants_prevplots.jpg", device="jpg", dpi=300, height=30, width=30, unit="cm")


## exclude contaminants from their respective batches
ps_decont <- subset_taxa(ps, !group %in% dropARG); ps_decont

# saveRDS(ps_decont, "megares_decont.RDS")

