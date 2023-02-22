rm(list=ls())
library(phyloseq)
library(ggplot2)
library(microbiome)
library(tidyverse)
theme_set(theme_bw())

ps <- readRDS("megares_decont.RDS")
sample_data(ps)$code <- factor(sample_data(ps)$code, levels = c("Employees [no AB]", "MSM [no AB]", "MSM [AB]"))

############################################################### (1) Calculation of RPM & merging by level of interest #############
ps.n <- ps; for(n in 1:nsamples(ps.n)){otu_table(ps.n)[,n] <- 1000000*otu_table(ps.n)[,n]/sample_data(ps.n)$reads.bacteria.bracken[n]}; ps.n # Normalization: divide otu column of each sample by total number of bacterial reads, multiplied by 1OE6 (reads/hits per million)
ps.group <- aggregate_rare(ps, level="group", detection = 0, prevalence = 0); ps.group
ps.group.n <- ps.group; for(n in 1:nsamples(ps.group.n)){otu_table(ps.group.n)[,n] <- 1000000*otu_table(ps.group.n)[,n]/sample_data(ps.group.n)$reads.bacteria.bracken[n]}; ps.group.n # Normalization: divide otu column of each sample by total number of bacterial reads, multiplied by 1OE6 (reads/hits per million)
ps.class <- aggregate_rare(ps, level="class", detection = 0, prevalence = 0); ps.class
ps.class.n <- ps.class; for(n in 1:nsamples(ps.class.n)){otu_table(ps.class.n)[,n] <- 1000000*otu_table(ps.class.n)[,n]/sample_data(ps.class.n)$reads.bacteria.bracken[n]}; ps.class.n # Normalization: divide otu column of each sample by total number of bacterial reads, multiplied by 1OE6 (reads/hits per million)

############################################################### (2) PREVALENCE #############
prevalence(ps.class.n, sort=T)
prevalence(ps.group.n, sort=T)

############################################################### (3) DIVERSITY #############
#####-----------------------------------------------------------(3.1) alpha -------------####
PS <- ps.n
# PS <- subset_samples(PS, sample_sums(PS)>0)

# Suppl_Fig7A
p <- plot_richness(PS, x="code", measures=c("Shannon", "InvSimpson"), color="code") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  labs(x="") +
  theme(legend.position="none")
p
# ggsave(plot=last_plot(), "Suppl_Fig7A_Megares_Alpha.jpg", device="jpg",dpi=300, height=10, width=10, unit="cm")

shannon <- subset(p$data, variable=="Shannon")
summary(shannon$value)
aggregate(value~code, FUN=quantile, data=shannon)
kruskal.test(value~code, data=shannon)
dunn.test::dunn.test(shannon$value, shannon$code, method="BH")

invsim <- subset(p$data, variable=="InvSimpson")
summary(invsim$value)
aggregate(value~code, FUN=quantile, data=invsim)
kruskal.test(value~code, data=invsim)
dunn.test::dunn.test(invsim$value, invsim$code, method="BH")

#####-----------------------------------------------------------(3.2) beta -------------####
#####------------------------------------------------------------------ pca -------------####
# BiocManager::install("MicrobiotaProcess") # https://github.com/YuLab-SMU/MicrobiotaProcess

PS <- ps.n
PS <- aggregate_rare(PS, level="group", prevalence=0, detection=0)
PS <- transform(PS, transform="clr")

# Suppl_Fig7B
pcares <- MicrobiotaProcess::get_pca(PS, method=NULL)
pcaplot1 <- MicrobiotaProcess::ggordpoint(obj=pcares, biplot=F, speciesannot=F, factorNames=c("code"), ellipse=TRUE, showsample=F, labelfactor="Sample.ID", max.overlaps=100) +
  guides(fill=guide_legend("Study group"))
pcaplot2 <- MicrobiotaProcess::ggordpoint(obj=pcares, pc=c(1, 3), biplot=F, speciesannot=TRUE, factorNames=c("code"), ellipse=TRUE) +
  theme(legend.position="none")
pcaplot1 | pcaplot2
# ggsave(plot=last_plot(), "Suppl_Fig7B_Megares_Beta.jpg", device="jpg",dpi=300, height=10, width=20, unit="cm")

#####------------------------------------------------------------------ permanova -------------####
library(vegan)
OTU <- t(abundances(PS))
DIST <- vegan::vegdist(OTU, method="euclidean")
permanova <- vegan::adonis(OTU ~ code, data=meta(PS), method="euclidean", permutations = 10000); permanova
permdisp <- vegan::betadisper(DIST, group=meta(PS)$code); permdisp
anova(permdisp)

############################################################### (4) ABUNDANCE #############
#####-----------------------------------------------------------(4.1) Heatmap -------------####
PS <- ps.n
PS <- tax_glom(PS, taxrank="group")
taxa_names(PS) <- data.frame(tax_table(PS))[,"group"]
PS <- microbiome::transform(PS, transform="clr")
otu <- abundances(PS)

df <- meta(PS)
df <- df %>% mutate_at(.vars=vars(ends_with(".b")), .funs=funs(factor(., levels=c(0,1), labels=c("No", "Yes"), ordered=FALSE)))

# Fig2
library(ComplexHeatmap)
col.heatmap <- colorRampPalette(c("white", "orange", "red"), space = "rgb")(99)
col_ha <- HeatmapAnnotation(df = df[,c("AB.183.ml.tsince.end", "AB.183.blcef.tsince.end", "AB.183.tc.tsince.end", "AB.183.fq.tsince.end")],
      col=list(
        AB.183.ml.tsince.end = circlize::colorRamp2(c(0, 183),c("blue", "white")),
        AB.183.blcef.tsince.end = circlize::colorRamp2(c(0, 183),c("green", "white")),
        AB.183.tc.tsince.end = circlize::colorRamp2(c(0, 183),c("red", "white")),
        AB.183.fq.tsince.end = circlize::colorRamp2(c(0, 183),c("yellow", "white"))
        ),
      na_col="white",
      annotation_label = c("Macrolide", "Betalactam", "Tetracycline", "Fluoroquinolone"),
      annotation_name_side = "left",
      border = TRUE,
      show_legend = FALSE
      )
row_ha = HeatmapAnnotation(Prevalence=anno_barplot(round(prevalence(PS),2), ylim=0:100, width=unit(2, "cm")), 
                           which="row",
                           annotation_name_side = "top",
                           annotation_name_rot = 0)
leg_ha = list(title="Relative abundance", 
              title_position = "lefttop-rot", 
              legend_height = unit(4, "cm"),
              border="black")

taxtab <- data.frame(tax_table(PS))
taxtab <- subset(taxtab, select=c("class", "group"))
taxtab <- taxtab[which(!duplicated(taxtab$group)),]
rowsplit= data.frame(group=rownames(otu))
rowsplit= plyr::join(rowsplit, taxtab, by="group")
rowsplit$class= ifelse(rowsplit$class=="Multi-drug_resistance", "MDR", rowsplit$class)
rowsplit$class= ifelse(rowsplit$class=="Aminoglycosides", "Ag", rowsplit$class)
rowsplit$class= ifelse(rowsplit$class=="Mupirocin", "Mp", rowsplit$class)
rowsplit$class= ifelse(rowsplit$class=="Phenicol", "Ph", rowsplit$class)
rowsplit$class= ifelse(rowsplit$class=="Fluoroquinolones", "Fq", rowsplit$class)
rowsplit$class= ifelse(rowsplit$class=="betalactams", "Bl", rowsplit$class)
rowsplit= factor(rowsplit$class)

colsplit= df$code

library(GetoptLong)
# pdf(qq("Fig2_ARG_group_class_heatmap.pdf"), width = 10, height = 6)
ht <- Heatmap(otu, name = "Relative abundance", col = col.heatmap, 
        bottom_annotation = col_ha,
        right_annotation = row_ha,
        
        show_column_names = FALSE, column_dend_reorder = F,
        show_row_names = TRUE, 
        
        cluster_column_slices=FALSE,
        
        show_column_dend = FALSE, column_dend_height = unit(2, "cm"), 
        show_row_dend=FALSE, row_dend_width = unit(2, "cm"),
        
        use_raster=FALSE, raster_quality=5, raster_resize_mat=max,
        column_names_gp = grid::gpar(fontsize = 8),
        row_names_gp = grid::gpar(fontsize = 8),
        
        column_split=colsplit,
        row_split = rowsplit,
        
        heatmap_legend_param = leg_ha
        )
draw(ht)
# dev.off()

#####-----------------------------------------------------------(4.2) Differential Abundance -------------####
# #------------------------------------------------------------------ ZINB & BOXPLOTS ####
library(pscl)

## ZINB regression
# overall
PS <- ps
sample_data(PS)$Abundance <- sample_sums(PS)
sample_data(PS)$code <- relevel(sample_data(PS)$code, "MSM [no AB]")
m.zinb <- zeroinfl(Abundance ~ code + offset(log(reads.bacteria.bracken)) | code,
               dist="negbin", link="logit", data=meta(PS))
sjPlot::tab_model(m.zinb)

# by group/class
LEVEL <- "class" # choose "class" or "group"

if(LEVEL=="group"){
  PS <- ps.group
  # PS <- subset_taxa(PS, class=="MLS")
  PS <- core(PS, detection=0, prevalence=0.1, include.lowest = F)
  sample_data(PS)$code <- relevel(sample_data(PS)$code, "MSM [no AB]")
} else if (LEVEL=="class"){
    PS <- ps.class
    PS <- core(PS, detection=0, prevalence=0.1, include.lowest = F)
    sample_data(PS)$code <- relevel(sample_data(PS)$code, "MSM [no AB]")
}
df <- psmelt(PS)
df[,LEVEL] <- gsub(" |-", "", df[,LEVEL])
df[,LEVEL] <- gsub(")", "", df[,LEVEL])
df[,LEVEL] <- gsub("\\(", "", df[,LEVEL])
model_list <- list()
for(LEV in df[,LEVEL]){
  DF <- df[which(df[,LEVEL]==LEV),]
  colnames(DF)[which(colnames(DF) == "Abundance")] <- LEV
  model_list[[LEV]] <- zeroinfl(as.formula(paste(LEV, "~ code + offset(log(reads.bacteria.bracken)) | code")), 
                                  dist="negbin", link="logit", data=DF)
}
sjPlot::tab_model(model_list)


## BOXPLOTS (Fig 1)
# overall
PS <- ps.n
sample_data(PS)$Abundance <- sample_sums(PS)
df <- meta(PS)

Fig1A <- ggplot(data = meta(PS), aes(x = code, y = Abundance)) +
  geom_boxplot(aes(color = code, fill=code, alpha=0.3), outlier.shape = NA) +
  geom_jitter(aes(color = code), height = 0, width = .2) +
  labs(x = NULL, y = "ARG abundance (RPM)", 
       # title = "Overall",
       tag="a")  +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.position = "none")+
  guides(alpha="none", color=guide_legend("Study group"), fill="none")
# Fig1A

# by class
PS <- ps.class.n
PS <- core(PS, detection=0, prevalence=0.1, include.lowest = F)
df <- psmelt(PS)
df$class <- ifelse(df$class=="Multi-drug_resistance", "MDR", df$class)
df$class <- ifelse(df$class=="betalactams", "Betalactams", df$class)

# comb <- split(t(combn(levels(sample_data(PS)$code), 2)), seq(nrow(t(combn(levels(sample_data(PS)$code), 2)))))
Fig1B <- ggplot(data = df, aes(x = code, y = Abundance)) +
  geom_boxplot(aes(color = code, fill=code, alpha=0.3), outlier.shape = NA) +
  geom_jitter(aes(color = code), height = 0, width = .2) +
  facet_wrap(~ class, nrow=1) + 
  labs(x = NULL, 
       y = "ARG abundance (RPM)",
       tag="b")+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())+
  guides(alpha="none", color="none", fill="none")
# Fig1B

# by ARG group
PS <- ps.group.n
PS <- core(PS, detection=0, prevalence=0.1, include.lowest = F)
df <- psmelt(PS)

# comb <- split(t(combn(levels(sample_data(PS)$code), 2)), seq(nrow(t(combn(levels(sample_data(PS)$code), 2)))))
Fig1C <- ggplot(data = df, aes(x = code, y = Abundance)) +
  geom_boxplot(aes(color = code, fill=code, alpha=0.3), outlier.shape = NA) +
  geom_jitter(aes(color = code), height = 0, width = .2) +
  facet_wrap(~ group, scales="free") + 
  labs(x = NULL, y = "ARG abundance (RPM)", 
       tag="c")+
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(), legend.position = "bottom")+
  guides(alpha="none", color=guide_legend("Study group"), fill="none")
# Fig1C

p <- gridExtra::grid.arrange(grobs=list(Fig1A, Fig1B, Fig1C), widths=c(1,2), layout_matrix=rbind(c(1,2), c(3,3), c(3,3)))
p
# ggsave("Fig1_ARG_boxplots.jpg", p, width=22, height=30, units="cm", dpi=300)


############################################################### (5) CORRELATION #############
####------------------------------------------------------------(5.1) SparCC #####
library(SpiecEasi)

# MEGARES correlation plot, according to SPARCC
MEGARES <- ps.n
OTU_MEG <- MEGARES %>%
  aggregate_rare(detection=0, prevalence=0.1, level="group") %>%
  abundances() %>%
  t()
net <- sparcc(OTU_MEG)
COR_MEG <- net$Cor
rownames(COR_MEG) <- colnames(COR_MEG) <- colnames(OTU_MEG)

# Suppl Fig 8A: SparCC correlation plot, taking into account bootstrapped p values, with MARco package (http://www.lifescipy.net/RcodeDB/MARco.html)
w.mb <- MARco::network.pipeline(otu=t(OTU_MEG),prevalence=0.1,alpha=0.05,bootstrap=99,cpu=2)
diag(w.mb) <- 1
ggcorrplot ::ggcorrplot(w.mb, hc.order = TRUE, type = "full", 
                        # p.mat = CONF, 
                        insig="blank", 
                        lab=TRUE,
                        tl.cex=8, ggtheme=theme_bw) +
  scale_fill_gradient2(limit = c(-1,1), low = "blue", high =  "red", mid = "white", midpoint = 0)
# ggsave("Suppl_Fig8A_Megares_SparCC_corrplot.jpg", plot=last_plot(), device="jpg",dpi=300, height=20, width=20, unit="cm")

# BRACKEN
BRACKEN <- readRDS("bracken_decont.RDS")

# BRACKEN correlation plot, according to SPARCC
OTU_BRA <- BRACKEN %>%
  aggregate_rare(detection=0, prevalence=0.25, level="Genus") %>%
  abundances() %>%
  t()
net <- sparcc(OTU_BRA)
COR_BRA <- net$Cor
rownames(COR_BRA) <- colnames(COR_BRA) <- colnames(OTU_BRA)

# Suppl Fig 8B: SparCC correlation plot, taking into account bootstrapped p values, with MARco package (http://www.lifescipy.net/RcodeDB/MARco.html)
w.mb <- MARco::network.pipeline(otu=t(OTU_BRA),prevalence=0,alpha=0.05,bootstrap=99,cpu=2)
diag(w.mb) <- 1
ggcorrplot ::ggcorrplot(w.mb, hc.order = TRUE, type = "full", 
                        # p.mat = CONF, insig="blank", 
                        lab=TRUE, lab_size=2,
                        tl.cex=8, ggtheme=theme_bw) +
  scale_fill_gradient2(limit = c(-1,1), low = "blue", high =  "red", mid = "white", midpoint = 0)
# ggsave("Suppl_Fig8B_Bracken_SparCC_corrplot.jpg", plot=last_plot(), device="jpg",dpi=300, height=20, width=20, unit="cm")

####------------------------------------------------------------(5.2) Spearman (Bracken vs Megares) #####
BRA <- data.frame(OTU_BRA)
BRA$Sample.ID <- rownames(BRA)
MEG <- data.frame(OTU_MEG)
MEG$Sample.ID <- rownames(MEG)
BRAMEG <- full_join(BRA, MEG, by="Sample.ID") %>%
  column_to_rownames("Sample.ID")

# Suppl Fig11
res <- psych::corr.test(x=OTU_MEG, y=OTU_BRA, method="spearman", adjust="none")
ggcorrplot ::ggcorrplot(res$r, hc.order = F, hc.method="ward.D", type = "full", 
                        p.mat = res$p.adj, insig="blank",
                        lab=TRUE, lab_size=2,
                        tl.cex=8,
                        ggtheme=theme_bw) +
  scale_fill_gradient2(limit = c(-1,1), low = "blue", high =  "red", mid = "white", midpoint = 0)
# ggsave(last_plot(), file = "Suppl_Fig11_Meg_Bra_Spearman_corrplot.jpg", device="jpg", dpi=300, height=20, width=20, unit="cm")
