rm(list=ls())
library(phyloseq)
library(ggplot2)
library(microbiome)
theme_set(theme_bw())

ps <- readRDS("bracken_decont.RDS")
sample_data(ps)$code <- factor(sample_data(ps)$code, levels = c("Employees [no AB]", "MSM [no AB]", "MSM [AB]"))

############################################################### (1) PREVALENCE #############
PS <- ps
PS <- aggregate_rare(PS, level="Genus", detection=0, prevalence=0.25)

prevalence(PS, detection=0, sort=TRUE)

# suppl fig 9A
Employees = data.frame(Genus = names(prevalence(PS)),
                       Prevalence = prevalence(subset_samples(PS, code=="Employees [no AB]")),
                       code="Employees [no AB]"
                       )
MSM_noAB = data.frame(Genus = names(prevalence(PS)),
                      Prevalence = prevalence(subset_samples(PS, code=="MSM [no AB]")),
                      code="MSM [no AB]"
                      )
MSM_AB = data.frame(Genus = names(prevalence(PS)),
                    Prevalence = prevalence(subset_samples(PS, code=="MSM [AB]")),
                    code="MSM [AB]"
                    )
df <- rbind(Employees, MSM_noAB, MSM_AB)
df$code <- factor(df$code, levels=c("Employees [no AB]", "MSM [no AB]", "MSM [AB]"))

ggplot(data=df, aes(x=reorder(Genus, Prevalence), y=Prevalence*100, fill=code)) +
  geom_bar(position="dodge", stat="identity") + coord_flip() +
  labs(y="Prevalence (%)", x="Genus") + guides(fill=guide_legend("Study group"))
# ggsave(plot=last_plot(), "Suppl_Fig9A_Bracken_Prevalence.jpg", device="jpg",dpi=300, height=12, width=15, unit="cm")

############################################################### (2) DIVERSITY #############
#####---------------------------------------------------------- (2.1) alpha -------------####
# library(patchwork); library(ggsignif)
PS <- ps
PS <- aggregate_rare(PS, level="Species", detection=0, prevalence=0)

# suppl fig 10A
p <- plot_richness(PS, x="code", measures=c("Shannon", "InvSimpson"), color="code") +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  labs(x="") +
  theme(legend.position="none")
p
# ggsave(plot=last_plot(), "Suppl_Fig10A_Bracken_Alpha.jpg", device="jpg",dpi=300, height=10, width=10, unit="cm")

shannon <- subset(p$data, variable=="Shannon")
summary(shannon$value)
aggregate(value~code, FUN=quantile, data=shannon)
kruskal.test(value~code, data=shannon)
dunn.test::dunn.test(shannon$value, shannon$code, method="BH", list=TRUE)

invsim <- subset(p$data, variable=="InvSimpson")
summary(invsim$value)
aggregate(value~code, FUN=quantile, data=invsim)
kruskal.test(value~code, data=invsim)
dunn.test::dunn.test(invsim$value, invsim$code, method="BH", list=TRUE)


#####-----------------------------------------------------------(2.2) beta -------------####
#####------------------------------------------------------------------ pca -------------####
# BiocManager::install("MicrobiotaProcess") # https://github.com/YuLab-SMU/MicrobiotaProcess

PS <- ps
PS <- aggregate_rare(PS, level="Species", detection=0, prevalence=0)
PS <- transform(PS, transform="clr")

# suppl fig 10B
pcares <- MicrobiotaProcess::get_pca(PS, method=NULL)
pcaplot1 <- MicrobiotaProcess::ggordpoint(obj=pcares, biplot=F, speciesannot=F, factorNames=c("code"), ellipse=TRUE, showsample=F, labelfactor="Sample.ID", max.overlaps=100) +
  guides(fill=guide_legend("Study group"))
pcaplot2 <- MicrobiotaProcess::ggordpoint(obj=pcares, pc=c(1, 3), biplot=F, speciesannot=TRUE, factorNames=c("code"), ellipse=TRUE) +
  theme(legend.position="none")
pcaplot1 | pcaplot2
# ggsave(plot=last_plot(), "Suppl_Fig10B_Bracken_Beta.jpg", device="jpg",dpi=300, height=10, width=20, unit="cm")

#####------------------------------------------------------------------ permanova -------------####
library(vegan)
OTU <- t(abundances(PS))
DIST <- vegan::vegdist(OTU, method="euclidean")
permanova <- vegan::adonis(OTU ~ code, data=meta(PS), method="euclidean", permutations = 10000); permanova
permdisp <- vegan::betadisper(DIST, group=meta(PS)$code); permdisp
anova(permdisp)

############################################################### (3) ABUNDANCE #############
#####-----------------------------------------------------------(3.1) Heatmap -------------####
PS <- ps
PS <- aggregate_rare(PS, level="Genus", detection=0, prevalence=0.25)
PS <- microbiome::transform(PS, transform="clr")
otu <- abundances(PS)

df <- meta(PS)
df <- df %>% mutate_at(.vars=vars(ends_with(".b")), .funs=funs(factor(., levels=c(0,1), labels=c("No", "Yes"), ordered=FALSE)))

# suppl Fig 9B
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
              # at = c(-10, 0, 10),
              legend_height = unit(4, "cm"),
              border="black")

colsplit= df$code

library(GetoptLong)
# pdf(qq("Suppl_Fig9B_Bracken_Abundance.pdf"), width = 10, height = 6)
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
              
              heatmap_legend_param = leg_ha
)
draw(ht)
# dev.off()

#####-----------------------------------------------------------(3.2) Differential Abundance -------------####
#####------------------------------------------------------------------ ANCOM-BC -------------####
# see also https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html
library(ANCOMBC)
PS <- ps
PS <- aggregate_rare(PS, level="Species", detection=0, prevalence=0.25)
sample_data(PS)$code <- factor(sample_data(PS)$code, labels = c("Employees", "MSM_noAB", "MSM_AB"), levels = c("Employees [no AB]", "MSM [no AB]", "MSM [AB]"))

## (1) ANCOM-BC with Employees as reference group
sample_data(PS)$code <- relevel(sample_data(PS)$code, "Employees")

ancom_da <- ANCOMBC::ancombc(PS, formula = "code",
                   p_adj_method = "BH", zero_cut = NULL, lib_cut = NULL, 
                   group="code", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                   max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)

employees <- data.frame(
  Species = row.names(ancom_da$res$beta),
  beta = unlist(ancom_da$res$beta),
  se = unlist(ancom_da$res$se),
  W = unlist(ancom_da$res$W),
  p_val = unlist(ancom_da$res$p_val),
  q_val = unlist(ancom_da$res$q_val),
  diff_abn = unlist(ancom_da$res$diff_abn),
  group="Reference = Employees [no AB]")

## (2) ANCOM-BC with MSM [no AB] as reference group
sample_data(PS)$code <- relevel(sample_data(PS)$code, "MSM_noAB")

ancom_da <- ANCOMBC::ancombc(PS, formula = "code",
                             p_adj_method = "BH", zero_cut = NULL, lib_cut = NULL, 
                             group="code", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                             max_iter = 100, conserve = TRUE, alpha = 0.05, global = FALSE)

msmnoab <- data.frame(
  Species = row.names(ancom_da$res$beta),
  beta = unlist(ancom_da$res$beta),
  se = unlist(ancom_da$res$se),
  W = unlist(ancom_da$res$W),
  p_val = unlist(ancom_da$res$p_val),
  q_val = unlist(ancom_da$res$q_val),
  diff_abn = unlist(ancom_da$res$diff_abn),
  group="Reference = MSM [no AB]")

df <- rbind(employees, msmnoab)

df$code <- ifelse(grepl("Employees", rownames(df)),
                  "Employees [no AB]", 
                  ifelse(grepl("MSM_noAB", rownames(df)),
                         "MSM [no AB]", 
                         "MSM [AB]")
                  )
df$code <- factor(df$code, levels=c("Employees [no AB]", "MSM [no AB]", "MSM [AB]"))
df$q_val_stars <- ifelse(df$q_val<0.001, "***", ifelse(df$q_val<0.01, "**", ifelse(df$q_val<0.05, "*", NA)))
df$q_val_stars_pos <- ifelse(df$beta<0, -6, 6)
keepTaxa <- unique(df$Species[which(df$diff_abn==TRUE, arr.ind=FALSE)]) # subset taxa indicated to be significant in at least one group
df <- df[df$Species %in% keepTaxa,]
df$struczero <- ifelse(df$q_val==0, 1, NA) # a q-value of 0 is assumed to be due to structural zeros
df$se <- ifelse(is.na(df$struczero), df$se, NA) # drop SE (error bars) for species with structural zeros

# Fig3
p <-  ggplot(data = df, aes(x = Species, y = beta, ymin = beta - 1.96*se, ymax = beta + 1.96*se, fill=code)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(width = 0.2, position = position_dodge(width=0.9),
                color = "black") + 
  labs(x = NULL, y = "Log fold change (95% confidence interval)", 
       title = "",
       subtitle="") + 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1)) + 
  geom_abline(intercept=0, slope=0) + 
  ylim(-6,6) +
  geom_text(aes(label=q_val_stars, y=q_val_stars_pos), position = position_dodge(width=0.9), color="black") +
  coord_flip() +
  guides(shape="none", fill=guide_legend("Study group"), color="none") + 
  facet_wrap(~group)
p
# ggsave(plot=p, "Fig3_Bracken_ANCOMBC.jpg", device="jpg",dpi=300, height=9, width=25, unit="cm")