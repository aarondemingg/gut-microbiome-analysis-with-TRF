---
title: "Gut microbiome analysis"
output: html_document
date: "2024-05-15"
---

# For paprica on linux. 
```{r}

#Full PAPRICA workflow

####Step 1: Paprica installation (Python)

git clone https://github.com/bowmanjeffs/paprica.git
cd paprica
chmod a+x *py
chmod a+x *sh


### the python package
conda create -n paprica_env python=3.10.12
pip install pandas
pip install biopython 
pip install joblib
pip install termcolor
pip install ete3
pip install seqmagick
conda install -c bioconda epa-ng -y
conda install -c bioconda gappa -y
conda install -c bioconda infernal -y
conda install -c bioconda hmmer -y 
conda install -c easel


## test
./paprica-run.sh test bacteria


####Step 2: individual excel files of seqtab for each sample  (R)
library(dplyr)
df<- readRDS(#your ps object) # input seqtab file, please change accordingly

a <- df %>% t %>% data.frame # convert to dataframe, transposed

for (x in 1:ncol(df)) { 
  b <- a[x]
  b
  c <- cbind(rownames(b), b)
  colnames(c) <- c('sequence', 'abundance')
  write.csv(c, paste0(colnames(a)[x], '.csv'), row.names = F)
}

# you should have individual csv files of all your samples from your seqtab file

####Step 3: Looping all CSV into Fasta files (R)
all_files <- list.files()
csv_files <- all_files[grepl("\\.csv$", all_files)]

for (i in csv_files){
  #print(i)
  
  df = read.csv(i, header=T, row.names=1)
  #head(df,2)
  
  m <- as.matrix(df)
  csum <- colSums(m)
  idx <- unlist(lapply(csum, seq_len), use.names=FALSE)
  res <- matrix(c(sprintf(">%s_%d", rep("seq", csum), idx), # id
                  rep(rownames(m)[row(m)], m)),                   # sequence
                nrow=2, byrow=TRUE)
  
  out <- paste0(substr(i, 1, nchar(i) - 4), ".fasta")
  
  writeLines(res, out)
  
}


# you should have individual fasta files after this step. 


####Step 4: PAPRICA

cd paprica
conda activate paprica_env
# I would suggest putting all your fasta files into a file and enter it

ls -1 | sed -e 's/\.fasta$//' > samples.txt #this would create a txt file of all your sample names

cd paprica #return to your paprica folder
nano my_analysis.sh #create a paprica loop command 
#!/bin/bash
while read f;do
   ./paprica-run.sh $f bacteria
done < samples.txt  #runs and loops all the samples in sample.txt into paprica


#### Step 5: Combining all paprica results into csv files

paprica-combine_results.py -domain bacteria -o 20240508out # change the name after -o to your preferred name for your files. 




```






# The list of libraries Ive used for this analysis (most if not all)

```{r}
library(BiocManager)
library(Rcpp)
library(dada2)
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggsignif)
library(knitr)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(microbiome)
library(dplyr)
library(ggpubr)
library(rstatix)
library(permute)
library(lattice)
library(vegan)
library(RColorBrewer)
library(reshape2)
library(extrafont)
library(Biostrings)
library(parallel)
library(BiocGenerics)
library(RSQLite)
library(DECIPHER)
library(nloptr)
library(magrittr)
library(qwraps2)
library(ANCOMBC)
library(DT)
library(biomformat)
library(microbiomeMarker)
library(data.table)
library(ALDEx2)
library(mixOmics)
library(viridis)
library(readr)
library(patchwork)
library(mixOmics)
library(phangorn)
library(mia)
library(TreeSummarizedExperiment)
library(bluster)
library(kableExtra)
library(scater)
library(miaViz)
library(ecodist)
library(microbiomeutilities)
library(ggpicrust2)
library(ggprism)
library(SpiecEasi)

```

# file load

# loading and phyloseq generation

```{r}

OTU = read.csv("newoutput.bacteria.edge_tally.csv", header = TRUE, row.names = 1)
TAX = as.matrix(read.csv("new.bacteria.taxon_map.csv", header = TRUE, row.names = 1))
meta <- read.csv('meta2.csv', header = TRUE, row.names = 1)


OTU[is.na(OTU)] <- 0



library(phyloseq)

tax <- tax_table(TAX)
otu <- otu_table(OTU, taxa_are_rows = FALSE)
mat = sample_data(meta)

#combined into a phyloseq object
combined <-phyloseq(otu, tax)
ps_obj = merge_phyloseq(combined, mat)


#filter taxa < 0.005% and present in >5% os the samples
filter <- phyloseq::genefilter_sample(ps_obj, filterfun_sample(function(x) x / sum(x) > 5e-3), A=0.05*nsamples(ps_obj))
ps <- prune_taxa(filter, ps_obj)




```

## Alpha diversity

```{r pressure, echo=FALSE}
ps_obj
ps.meta <- meta(ps_obj)
kable(head(ps.meta))
meta2 <- meta
tab <-microbiome::alpha(ps_obj, index = "all")
kable(head(tab))
meta2$Shannon <- tab$diversity_shannon
meta2$InverseSimpson <- tab$diversity_inverse_simpson #adding diversity table to metadata
meta2$Pielou <- tab$evenness_pielou
meta2$Chao1 <- tab$chao1

meta2[, c(1,3,4,5,6)] <- list(NULL)
meta2

df_diversity <- meta2 %>% 
  pivot_longer(cols = c(Shannon, InverseSimpson, Pielou, Chao1), 
               names_to = "Diversity", values_to = "value") 
df_diversity

# Plotting multiple barplots

df_diversity$Intervention <- as.factor(df_diversity$Intervention)



ggboxplot(df_diversity,  x= "Intervention", y = "value", fill="Diversity", palette = "jco",
          ylab = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  stat_compare_means(method = "anova") +
  facet_wrap(~Diversity, scales = "free") + theme_bw()


df_diversity$Intervention <- as.factor(df_diversity$Intervention)

stat.test2 <- df_diversity %>%
  group_by(Diversity) %>%
  tukey_hsd(value ~ Intervention) %>%
  add_xy_position()
stat.test2

stat.test2 <- apply(stat.test2,2,as.character)
as.data.frame(stat.test2)
write.csv(stat.test2, "alpha_stat.csv")
stat.test2 <- read.csv("alpha_stat.csv")

alpha <- ggboxplot(df_diversity,  x= "Intervention", y = "value", fill="Diversity", palette = "jco",
          ylab = "Diverstity") +
  stat_pvalue_manual(stat.test2, label = "p.adj.signif", 
                     hide.ns = TRUE, y.position = "y.position") +
  stat_compare_means(method = "anova")+ 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  facet_wrap(~Diversity, scales = "free") + theme_bw()
alpha #1200 x 670 probably rename to save space at labelling


alpha <- ggboxplot(df_diversity,  x= "Intervention", y = "value", fill="Diversity", palette = "jco",
          ylab = "Diverstity") +
  facet_wrap(~Diversity, scales = "free") + theme_bw()
alpha + stat_pvalue_manual(stat.test2, label = "p.adj.signif", 
                     hide.ns = TRUE, y.position = "y.position") #1200 x 670 probably rename to save space at labelling
alpha + 
  stat_pvalue_manual(stat.test2, label = "p.adj.signif", 
                     hide.ns = TRUE, y.position = "y.position") + theme(text = element_text(family = "serif"))

```

# beta diversity

```{r}

library(miaViz)
library(mia)
library(vegan)
tse <- makeTreeSummarizedExperimentFromPhyloseq(ps_obj)

tse <- transformAssay(tse,
                      method = "clr",
                      pseudocount = TRUE)

tse <- runMDS(tse,
              FUN = vegan::vegdist,
              method = "euclidean",
              assay.type = "clr",
              name = "MDS_aitchison")


plotReducedDim(tse, "MDS_aitchison", 
               colour_by = "Intervention") + theme_classic() 

#for supervised ordination
tse2 <- runRDA(tse, 
              assay.type = "clr",
              formula = assay ~Intervention,
              distance = "euclidean",
              name = "Aitchison")

rda_info <- attr(reducedDim(tse2, "Aitchison"), "significance")
rda_info

rda_info$permanova |>
  knitr::kable()

rda_info$homogeneity |>
  knitr::kable()

Aitchison_microbe <- plotRDA(tse2, "Aitchison", colour_by = "Intervention", label.size = 2) + theme_bw()
Aitchison_microbe

Aitchison_microbe <- plotRDA(tse2, "Aitchison", colour_by = "Intervention", label.size = 2,
                             add.significance = FALSE, add.expl.var = FALSE, 
                             parse.labels = FALSE, add.vectors = FALSE) + theme_bw()
Aitchison_microbe + theme(text = element_text(family = "serif"))

ps.meta <- sample_data(ps_obj) #global + pairwise adonis found on this line and below
perma.meta <- meta(ps_obj)

perma.otu <- abundances(ps_obj)
permanova <-adonis2(t(perma.otu) ~ Intervention, 
                   data = perma.meta, 
                   permutations = 999, method ="bray")

permanova <- adonis2(perma.otu ~ Intervention,
                     data = perma.meta) # i have no idea why removing the t works for this line of code

permanova
print(as.data.frame(permanova$aov.tab)["Intervention", "Pr(>F)"]) #general difference in all the groups

library(pairwiseAdonis)
adonis.pair <- pairwise.adonis(perma.otu, factors = ps.meta$Intervention, sim.function='vegdist',sim.method='euclidian',p.adjust.m='BH')
adonis.pair
write.csv(adonis.pair, "adonis.csv")


unique(tse$Intervention)
tse$Intervention %>% table()
tse_subset_by_sample <- tse[ , tse$Intervention %in% 
                               c("SOTT", "SO", "TRF", "TRFTT", "TT")]

dim(tse_subset_by_sample)


tse_subset_by_sample <- runRDA(tse_subset_by_sample, 
                               assay.type = "clr",
                               formula = assay ~Intervention,
                               distance = "euclidean",
                               name = "Aitchison")


rda_info_subset <- attr(reducedDim(tse_subset_by_sample, "Aitchison"), "significance")
rda_info_subset

rda_info_subset$permanova |>
  knitr::kable()

rda_info_subset$homogeneity |>
  knitr::kable()

Aitchison_microbe_subset <- plotRDA(tse_subset_by_sample, "Aitchison", colour_by = "Intervention", label.size =                               2, add.significance = FALSE, add.expl.var = FALSE, 
                             parse.labels = FALSE, add.vectors = FALSE) + theme_bw()
Aitchison_microbe_subset + theme(text = element_text(family = "serif"))



```


# drawing the heatmap to look at the community
```{r}

#ccommunity composition 
#tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
tse_phylum <- mergeFeaturesByRank(tse, rank ="Species", onRankOnly=TRUE)
top_taxa <- getTopFeatures(tse_phylum,top = 5, assay.type = "clr")
phylum_renamed <- lapply(rowData(tse)$Genus,
                         function(x){if (x %in% top_taxa) {x} else {"Other"}})

rowData(tse)$Phylum <- as.character(phylum_renamed)

plotAbundance(tse, assay.type="clr", rank = "Phylum",
              order_rank_by="abund")

# Add clr-transformation on samples
tse_phylum <- transformAssay(tse_phylum, assay.type = "counts",
                             method = "relabundance", pseudocount = 1)

tse_phylum <- transformAssay(tse_phylum, assay.type = "relabundance",
                             method = "clr")

# Add z-transformation on features (taxa)
tse_phylum <- transformAssay(tse_phylum, assay.type = "clr", 
                             MARGIN = "features",
                             method = "z", name = "clr_z")
df <- meltAssay(tse_phylum, assay.type = "clr_z")

# Determines the scaling of colours
maxval <- round(max(abs(df$clr_z)))
limits <- c(-maxval, maxval)
breaks <- seq(from = min(limits), to = max(limits), by = 0.5)
colours <- c("darkblue", "blue", "white", "red", "darkred")

# Creates a ggplot object
ggplot(df, aes(x = SampleID, y = FeatureID, fill = clr_z)) +
  geom_tile() +
  scale_fill_gradientn(name = "CLR + Z transform", 
                       breaks = breaks, limits = limits, colours = colours) + 
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.key.size = unit(1, "cm")) +
  labs(x = "Samples", y = "Taxa")







# Takes subset: only samples from Intervention
tse_phylum
tse_phylum_subset <- tse_phylum[ , tse_phylum$ColData$Intervention]

# Add clr-transformation
tse_phylum_subset <- transformAssay(tse_phylum_subset,
                                    method = "clr")

tse_phylum_subset <- transformAssay(tse_phylum_subset, assay.type = "clr",
                                    MARGIN = "features", 
                                    method = "z", name = "clr_z")

# Get n most abundant taxa, and subsets the data by them
top_taxa <- getTopFeatures(tse_phylum, top = 20)
tse_phylum_subset <- tse_phylum[top_taxa, ]

# Gets the assay table
mat <- assay(tse_phylum_subset, "clr_z")


# Creates the heatmap
library(pheatmap)
pheatmap(mat, fontfamily = "serif")


library(ape)
taxa_hclust <- hclust(dist(mat), method = "complete")

taxa_tree <- as.phylo(taxa_hclust)

library(ggtree)

# Plot taxa tree
taxa_tree <- ggtree(taxa_tree) + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins


# Get order of taxa in plot
taxa_ordered <- get_taxa_name(taxa_tree)

taxa_tree



taxa_clusters <- cutree(tree = taxa_hclust, k = 3)

# Converts into data frame
taxa_clusters <- data.frame(clusters = taxa_clusters)
taxa_clusters$clusters <- factor(taxa_clusters$clusters)

# Order data so that it's same as in phylo tree
taxa_clusters <- taxa_clusters[taxa_ordered, , drop = FALSE] 

# Prints taxa and their clusters
taxa_clusters


# Adds information to rowData
rowData(tse_phylum_subset)$clusters <- taxa_clusters[order(match(rownames(taxa_clusters), rownames(tse_phylum_subset))), ]

# Prints taxa and their clusters
rowData(tse_phylum_subset)$clusters

#this row is for sample data

# Hierarchical clustering
sample_hclust <- hclust(dist(t(mat)), method = "complete")

# Creates a phylogenetic tree
sample_tree <- as.phylo(sample_hclust)

# Plot sample tree
sample_tree <- ggtree(sample_tree) + layout_dendrogram() + 
  theme(plot.margin=margin(0,0,0,0)) # removes margins

# Get order of samples in plot
samples_ordered <- rev(get_taxa_name(sample_tree))

sample_tree

sample_clusters <- factor(cutree(tree = sample_hclust, k = 3))

# Converts into data frame
sample_data <- data.frame(clusters = sample_clusters)

# Order data so that it's same as in phylo tree
sample_data <- sample_data[samples_ordered, , drop = FALSE] 

# Order data based on 
tse_phylum_subset <- tse_phylum_subset[ , rownames(sample_data)]

# Add sample type data
sample_data$Intervention <- unfactor(colData(tse_phylum_subset$Intervention))

sample_data$Intervention <- tse_phylum_subset$Intervention

sample_data

sample_data <- subset(sample_data, select = Intervention)


anno.colors <- list(anno = c(Cont = "lightblue", 
                                          SO = "orange", 
                                          SOTT = "seagreen",
                                          TRF = "firebrick1",
                                          TRFTT = "purple",
                                          TT = "coral4"))

heatmap <- pheatmap(mat, annotation_col = sample_data,
                    fontfamily = "serif", 
                    fontsize = 11,
                    annotation_colors = anno.colors)

heatmap  



```


## Deseq2

```{r}
library(microbiomeMarker)
library(ggsci)

deseq <- run_deseq2(ps_obj,
                    group = "Intervention",
                    taxa_rank = "Species",
                    norm = "rarefy",
                    pvalue_cutoff = 0.05,
                    p_adjust = "BH")
deseq 
marker_table(deseq)
plot_abundance(deseq, group = "Intervention")
deseq_plot1 <- plot_ef_dot(deseq)
deseq_plot2 <- plot_ef_bar(deseq) + scale_color_npg()
deseq_plot2 + scale_color_jco()
deseq_plot2 + theme(text = element_text(family = "serif"))

plot_abundance(deseq, group = "Intervention") + theme(text = element_text(family = "serif"))

deseq_plot2 <- plot_ef_bar(deseq) 

deseq_plot2 +  scale_fill_manual(values = c("Cont" = "skyblue", "SO" = "tomato",
                                     "SOTT" = "darkseagreen3", "TRF" = "turquoise",
                                      "TRFTT" = "yellow3", "TT" = "palevioletred1"))
deseq_plot2 

    
lefse <- run_lefse(ps_obj,
                   group = "Intervention",
                   taxa_rank = "Species",
                   kw_cutoff = 0.05,
                   multigrp_strat = TRUE)
lefse

lefse_plot1 <- plot_ef_bar(lefse) + scale_fill_manual(values = c("Cont" = "lightblue", "SO" = "orange",
                                     "SOTT" = "seagreen", "TRF" = "firebrick1",
                                      "TRFTT" = "purple", "TT" = "coral4"))

ps_obj_species <- subset_samples(ps_obj,
                                 Intervention %in% c("Cont", "TRF"))

ancom_out <- run_ancombc(ps_obj,
                         group = "Intervention",
                         contrast = NULL,
                         p_adjust = "BH")
ancom_out

aldx2 <- run_aldex(ps_obj_species, #can only use for two conditions. Raw reads recommended for ALDEx2
                   group = "Intervention",
                   taxa_rank = "all",
                   transform = "log10p",
                   norm = "none", 
                   method = "wilcox.test",
                   p_adjust = "holm",
                   pvalue_cutoff = 0.05,
                   paired = FALSE)
aldx2
```


# ancombc bacteria DA

```{r}

#ancombc

output = ancombc2(ps_obj, assay_name = "counts", tax_level = "Species",
                  fix_formula = "Intervention", rand_formula = NULL,
                  p_adj_method = "holm", 
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "Intervention", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = TRUE, 
                  dunnet = TRUE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)
output

res <- output$res_pair
write.csv(res, "res.csv")
View(res)

head(res)



df_fig_pair1.1 = res %>%
  dplyr::filter(diff_InterventionVehicle_unvaccinated == 1 |
                  diff_InterventionTRF_unvaccinated == 1|
                  diff_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1 |
                  diff_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_InterventionVehicle_unvaccinated == 1,
                              round(lfc_InterventionVehicle_unvaccinated,2), 0),
                lfc2 = ifelse(diff_InterventionTRF_unvaccinated == 1,
                              round(lfc_InterventionTRF_unvaccinated, 2), 0),
                lfc3 = ifelse(diff_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1, 
                              round(lfc_InterventionVehicle_unvaccinated_InterventionControl_vaccinated, 2), 0),
                lfc4 = ifelse(diff_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1,
                              round(lfc_InterventionTRF_vaccinated_InterventionControl_vaccinated, 2), 0),
                lfc5 = ifelse(diff_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1,
                              round(lfc_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated, 2), 0),
                lfc6 = ifelse(diff_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1, 
                              round(lfc_InterventionVehicle_vaccinated_InterventionTRF_vaccinated, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6,
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_pair2.1 = res %>%
  dplyr::filter(diff_InterventionVehicle_unvaccinated == 1 |
                  diff_InterventionTRF_unvaccinated == 1|
                  diff_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1 |
                  diff_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1) %>%
  dplyr::mutate(lfc1 = ifelse(passed_ss_InterventionVehicle_unvaccinated ==1,
                              "green", "black"),
                lfc2 = ifelse(passed_ss_InterventionTRF_unvaccinated == 1,
                              "green", "black"),
                lfc3 = ifelse(passed_ss_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1, 
                              "green", "black"),
                lfc4 = ifelse(passed_ss_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1,
                              "green", "black"),
                lfc5 = ifelse(passed_ss_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1,
                              "green", "black"),
                lfc6 = ifelse(passed_ss_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1,
                              "green", "black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6,
                      names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_pair1.2 = df_fig_pair1.1 %>%
  dplyr::left_join(df_fig_pair2.1, by = c("taxon", "group"))

df_fig_pair1.2$group = recode(df_fig_pair1.2$group,
                           `lfc1` = "Cont vs SO",
                           `lfc2` = "Cont vs TRF",
                           `lfc3` = "TT vs SOTT",
                           `lfc4` = "TT vs TRFTT",
                           `lfc5` = "SO vs TRF",
                           `lfc6` = "SOTT vs TRFTT")

df_fig_pair1.2$group = factor(df_fig_pair1.2$group,
                           levels = c("Cont vs SO",
                                      "Cont vs TRF",
                                      "TT vs SOTT",
                                      "TT vs TRFTT",
                                      "SO vs TRF",
                                      "SOTT vs TRFTT"))
lo = floor(min(df_fig_pair1.2$value))
up = ceiling(max(df_fig_pair1.2$value))
mid = (lo1 + up1)/2

fig_pair2 = df_fig_pair1.2 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold change") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_pair2


library(readxl)
library(extrafont)
library(device = "win")
windowsFonts()

cont_trf <- read_xlsx("cont_trf.xlsx")
cont_so <- read_xlsx("cont_soy.xlsx")
trftt_tt <- read_xlsx("trftt_tt.xlsx")
sott_tt <- read_xlsx("sott_tt.xlsx")
so_trf <- read_xlsx("so_trf.xlsx")
sott_trftt <- read_xlsx("sott_trft.xlsx")


cont_trf1 <- ggplot(data = cont_trf, 
           aes(x = reorder(taxon, lfc_InterventionTRF_unvaccinated), y = lfc_InterventionTRF_unvaccinated, )) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black")  +
  geom_errorbar(aes(ymin = lfc_InterventionTRF_unvaccinated - se_InterventionTRF_unvaccinated, 
                    ymax = lfc_InterventionTRF_unvaccinated + se_InterventionTRF_unvaccinated), width = 0.2,
                position = position_dodge(0.5), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Cont-TRF") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y.left = element_text(vjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1,),
        family = "serif") + 
  scale_fill_manual(values=c("#FF8080","#87CEFA")) + coord_flip()
cont_trf1 + theme(text = element_text(family = "serif"))

cont_so1<- ggplot(data = cont_so, 
           aes(x = reorder(taxon, lfc_InterventionVehicle_unvaccinated), y = lfc_InterventionVehicle_unvaccinated)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black")  +
  geom_errorbar(aes(ymin = lfc_InterventionVehicle_unvaccinated - se_InterventionVehicle_unvaccinated, 
                    ymax = lfc_InterventionVehicle_unvaccinated + se_InterventionVehicle_unvaccinated), width = 0.2,
                position = position_dodge(0.5), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Cont-SO") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y.left = element_text(vjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1)) + 
  scale_fill_manual(values=c("#FF8080","#87CEFA")) + coord_flip()
cont_so1 + theme(text = element_text(family = "serif"))

trftt_tt1<- ggplot(data = trftt_tt, 
           aes(x = reorder(taxon, lfc_InterventionTRF_vaccinated_InterventionControl_vaccinated ), y = lfc_InterventionTRF_vaccinated_InterventionControl_vaccinated )) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black")  +
  geom_errorbar(aes(ymin = lfc_InterventionTRF_vaccinated_InterventionControl_vaccinated - se_InterventionTRF_vaccinated_InterventionControl_vaccinated, 
                    ymax = lfc_InterventionTRF_vaccinated_InterventionControl_vaccinated + se_InterventionTRF_vaccinated_InterventionControl_vaccinated), width = 0.2,
                position = position_dodge(0.5), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "TRFTT-TT") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y.left = element_text(vjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1)) + 
  scale_fill_manual(values=c("#FF8080","#87CEFA")) + coord_flip()
trftt_tt1 + theme(text = element_text(family = "serif"))

sott_tt1 <- ggplot(data = sott_tt, 
           aes(x = reorder(taxon, lfc_InterventionVehicle_vaccinated_InterventionControl_vaccinated  ), y = lfc_InterventionVehicle_vaccinated_InterventionControl_vaccinated  )) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black")  +
  geom_errorbar(aes(ymin = lfc_InterventionVehicle_vaccinated_InterventionControl_vaccinated  - se_InterventionVehicle_vaccinated_InterventionControl_vaccinated, 
                    ymax = lfc_InterventionVehicle_vaccinated_InterventionControl_vaccinated  + se_InterventionVehicle_vaccinated_InterventionControl_vaccinated), width = 0.2,
                position = position_dodge(0.5), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "SOTT-TT") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y.left = element_text(vjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1)) + 
  scale_fill_manual(values=c("#FF8080","#87CEFA")) + coord_flip()
sott_tt1 + theme(text = element_text(family = "serif"))


so_trf1 <- ggplot(data = so_trf, 
           aes(x = reorder(taxon, lfc_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated   ), y = lfc_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated   )) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black")  +
  geom_errorbar(aes(ymin = lfc_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated   - se_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated, 
                    ymax = lfc_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated   + se_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated), width = 0.2,
                position = position_dodge(0.5), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "SO-TRF") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y.left = element_text(vjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1)) + 
  scale_fill_manual(values=c("#FF8080","#87CEFA")) + coord_flip()
so_trf1 + theme(text = element_text(family = "serif"))

sott_trftt1 <- ggplot(data = sott_trftt, 
           aes(x = reorder(taxon, lfc_InterventionVehicle_vaccinated_InterventionTRF_vaccinated), y = lfc_InterventionVehicle_vaccinated_InterventionTRF_vaccinated)) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5, color = "black")  +
  geom_errorbar(aes(ymin = lfc_InterventionVehicle_vaccinated_InterventionTRF_vaccinated - se_InterventionVehicle_vaccinated_InterventionTRF_vaccinated, 
                    ymax = lfc_InterventionVehicle_vaccinated_InterventionTRF_vaccinated + se_InterventionVehicle_vaccinated_InterventionTRF_vaccinated), width = 0.2,
                position = position_dodge(0.5), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "SOTT-TRFTT") + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.y.left = element_text(vjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1)) + 
  scale_fill_manual(values=c("#FF8080","#87CEFA")) + coord_flip()
sott_trftt1 + theme(text = element_text(family = "serif"))
 

new_ancom <- read.csv("ancom_data_combined.csv")

new_ancom_plot <- ggplot(data = new_ancom, 
                         aes(x = reorder(taxon, lfc), y = lfc, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se), width = 0.2, position = position_dodge(0.7)) +
  theme_bw() + coord_flip()

new_ancom_plot


```


aldex2
```{r}

new_df <- as.data.frame(otu_table(ps_obj))
t_new_df <- transpose(new_df)
colnames(t_new_df) <- rownames(new_df)
rownames(t_new_df) <- colnames(new_df)

str(t_new_df)

tt_new_df<- t_new_df %>% 
  mutate_if(is.numeric, round)

meta

aldx2 <- aldex.clr(tt_new_df, meta$Intervention, mc.samples = 128,)


kw.aldx2 <- aldex.kw(aldx2)
kw.aldx2

glm1 <- aldex.glm(aldx2, verbose = FALSE, fdr.method = "BH") # I am not sure what to do with this 
glm1


ps_obj2 <- ps_obj

otu_table(ps_obj2) <- otu_table(tt_new_df, taxa_are_rows = TRUE)

tse3 <- makeTreeSummarizedExperimentFromPhyloseq(ps_obj2)

tse3
x <- aldex.clr(
  reads = assay(tse3),
  conds = colData(tse3)$Intervention, 
  # 128 recommened for ttest, 1000 for rigorous effect size calculation
  mc.samples = 128, 
  denom = "all",
  verbose = FALSE
)

# calculates expected values of the Welch's t-test and Wilcoxon rank test on
# the data returned by aldex.clr

kw.aldx <- aldex.kw(x)
kw.aldx

glm1 <- aldex.glm(x, verbose = FALSE, fdr.method = "BH")

```


chunk of code for something else?

```{r}
data(selex)
selex <- selex[1200:1600, ]

covariates <- data.frame("A" = sample(0:1, 14, replace = TRUE),
                         "B" = c(rep(0, 7), rep(1, 7))) #to subset the data into two variables
mm <- model.matrix(~ A + B, covariates)
x1 <- aldex(selex, mm, mc.samples = 128)
x <- aldex.clr(selex, mm, mc.sample = 4, denom = "all") 

glm.test <- aldex.glm(x)
glm.eff <- aldex.glm.effect(x)

aldex.glm.plot(glm.test, eff = glm.eff, contrast = "B", type = "MW", post.hoc = "holm")


```



For enzymes 
```{r}


OTU_pathway = read.csv("newoutput.bacteria.path_tally.csv", header = TRUE, row.names = 1)


meta <- read.csv('meta2.csv', header = TRUE, row.names = 1)
tax_pathway <- as.matrix(read.csv("path_tax.csv", header = TRUE, row.names = 1))

otu_pathway <- otu_table(OTU_pathway, taxa_are_rows = FALSE)
tax_pathway <- tax_table(tax_pathway)
mat = sample_data(meta)

#combined into a phyloseq object
ps_pathway <-phyloseq(otu_pathway, mat, tax_pathway)
ps_pathway

filter <- phyloseq::genefilter_sample(ps_pathway, filterfun_sample(function(x) x / sum(x) > 5e-3), A=0.05*nsamples(ps_obj))
ps_pathway_new <- prune_taxa(filter, ps_pathway)


```



beta div for enzymes

```{r}

ps_pathway
tse2 <- makeTreeSummarizedExperimentFromPhyloseq(ps_pathway)

tse2 <- transformAssay(tse2,
                      method = "clr",
                      pseudocount = TRUE)

tse2 <- runMDS(tse2,
              FUN = vegan::vegdist,
              method = "euclidean",
              assay.type = "clr",
              name = "MDS_aitchison")

plotReducedDim(tse2, "MDS_aitchison", 
               colour_by = "Intervention") + theme_classic() 

#for supervised ordination

tse2 <- runRDA(tse2, 
              assay.type = "clr",
              formula = assay ~Intervention,
              distance = "euclidean",
              name = "Aitchison")

rda_info1 <- attr(reducedDim(tse2, "Aitchison"), "significance")
rda_info1

rda_info1$permanova |>
  knitr::kable()

rda_info1$homogeneity |>
  knitr::kable()

aitchison_pathway <- plotRDA(tse2, "Aitchison", colour_by = "Intervention",
                             add.significance = FALSE, add.expl.var = FALSE,
                             parse.labels = FALSE, add.vectors = FALSE) + theme_bw()


aitchison_pathway + theme(text = element_text(family = "serif"))


tse2

clr_assay <- assays(tse)$clr
clr_assay <- t(clr_assay)

euclidean_dist <- vegan::vegdist(clr_assay, method = "euclidean")

euclidean_pcoa <- ecodist::pco(euclidean_dist)

euclidean_pcoa_df <- data.frame(pcoa1 = euclidean_pcoa$vectors[,1], 
                                pcoa2 = euclidean_pcoa$vectors[,2])

euclidean_Intervention_pcoa_df <- cbind(euclidean_pcoa_df,
                             Intervention = colData(tse2)$Intervention)


euclidean_plot <- ggplot(data = euclidean_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
  geom_point() +
  labs(x = "PC1",
       y = "PC2",
       title = "Euclidean PCoA with CLR transformation") +
  theme(title = element_text(size = 12)) # makes titles smaller

euclidean_plot







```

DA for enzymes

```{r}
library(microbiomeMarker)

deseq1 <- run_deseq2(ps_pathway,
                    group = "Intervention",
                    taxa_rank = "Species",
                    pvalue_cutoff = 0.01,
                    p_adjust = "BH")
deseq1

marker_table(deseq1)
plot_abundance(deseq1, group = "Intervention")
deseq_plot1_pathways <- plot_ef_dot(deseq1)
deseq_plot2_pathways <- plot_ef_bar(deseq1) + theme(text = element_text(family = "serif"))

deseq_plot2_pathways + scale_fill_manual(values = c("Cont" = "lightblue", "SO" = "orange",
                                                    "SOTT" = "seagreen", "TRF" = "firebrick1",
                                                    "TRFTT" = "purple", "TT" = "coral4"))

lefse1 <- run_lefse(ps_pathway,
                   group = "Intervention",
                   taxa_rank = "Species",
                   kw_cutoff = 0.01,
                   multigrp_strat = TRUE)
lefse1

lefse_plot_pathway <- plot_ef_bar(lefse1)
lefse_plot_pathway + theme(text = element_text(family = "serif")) +scale_fill_manual(values = c("Cont" = "lightblue", "SO" = "orange",
                                                    "SOTT" = "seagreen", "TRF" = "firebrick1",
                                                    "TRFTT" = "purple", "TT" = "coral4"))

limmavoom <- run_limma_voom(ps_pathway,
                   group = "Intervention",
                   taxa_rank = "Species",
                   transform = "log10p",
                   norm = "rarefy",
                   pvalue_cutoff = 0.05)
limmavoom 
limma_plot <- plot_ef_bar(limmavoom) #Error in .plot_ef(mm, label_level, max_label_len, markers, "bar") : The effect size must be one of lda, diff_mean, eta_squared, logFC, clr_diff_mean, clr_F_statistic, W, imp, LR or F
limma.dat <- as.data.frame(marker_table(limmavoom))
limma.dat <- as_tibble(limma.dat)

limma.dat %>%
  as.character(limma.dat$enrich_group) %>%
  mutate(enrich_group = fct(reorder(enrich_group, ef_F_statistic)))

limma_plot <- ggplot(limma.dat, aes(x = fct_reorder(feature, enrich_group), y = ef_F_statistic, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  coord_flip()
limma_plot + theme_bw() + theme(text = element_text(family = "serif")) + scale_fill_manual(values = c("Cont" = "lightblue", "SO" = "orange",
                                                    "SOTT" = "seagreen", "TRF" = "firebrick1",
                                                    "TRFTT" = "purple", "TT" = "coral4")) 

aldx2 <- run_aldex(ps_pathway, #can only use for two conditions
                   group = "Intervention",
                   taxa_rank = "all",
                   transform = "log10p",
                   norm = "none", 
                   method = "wilcox.test",
                   p_adjust = "holm",
                   pvalue_cutoff = 0.05,
                   paired = FALSE)
aldx2


output1 = ancombc2(ps_pathway, assay_name = "counts", rank = "Species",
                  fix_formula = "Intervention", rand_formula = NULL,
                  p_adj_method = "holm", 
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "Intervention", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = TRUE, 
                  dunnet = TRUE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)
output1

res2 <- output1$res_pair
write.csv(res2, "res_pathway.csv")
View(res2)

head(res2)

df_fig_pair1.1 = res2 %>%
  dplyr::filter(diff_InterventionVehicle_unvaccinated == 1 |
                  diff_InterventionTRF_unvaccinated == 1|
                  diff_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1 |
                  diff_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_InterventionVehicle_unvaccinated == 1,
                              round(lfc_InterventionVehicle_unvaccinated,2), 0),
                lfc2 = ifelse(diff_InterventionTRF_unvaccinated == 1,
                              round(lfc_InterventionTRF_unvaccinated, 2), 0),
                lfc3 = ifelse(diff_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1, 
                              round(lfc_InterventionVehicle_unvaccinated_InterventionControl_vaccinated, 2), 0),
                lfc4 = ifelse(diff_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1,
                              round(lfc_InterventionTRF_vaccinated_InterventionControl_vaccinated, 2), 0),
                lfc5 = ifelse(diff_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1,
                              round(lfc_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated, 2), 0),
                lfc6 = ifelse(diff_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1, 
                              round(lfc_InterventionVehicle_vaccinated_InterventionTRF_vaccinated, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6,
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_pair2.1 = res2 %>%
  dplyr::filter(diff_InterventionVehicle_unvaccinated == 1 |
                  diff_InterventionTRF_unvaccinated == 1|
                  diff_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1 |
                  diff_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1) %>%
  dplyr::mutate(lfc1 = ifelse(passed_ss_InterventionVehicle_unvaccinated ==1,
                              "aquamarine3", "black"),
                lfc2 = ifelse(passed_ss_InterventionTRF_unvaccinated == 1,
                              "aquamarine3", "black"),
                lfc3 = ifelse(passed_ss_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1, 
                              "aquamarine3", "black"),
                lfc4 = ifelse(passed_ss_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1,
                              "aquamarine3", "black"),
                lfc5 = ifelse(passed_ss_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1,
                              "aquamarine3", "black"),
                lfc6 = ifelse(passed_ss_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1,
                              "aquamarine3", "black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6,
                      names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_pair1.2 = df_fig_pair1.1 %>%
  dplyr::left_join(df_fig_pair2.1, by = c("taxon", "group"))

df_fig_pair1.2$group = recode(df_fig_pair1.2$group,
                           `lfc1` = "Cont-SO",
                           `lfc2` = "Cont-TRF",
                           `lfc3` = "TT-SOTT",
                           `lfc4` = "TT-TRFTT",
                           `lfc5` = "SO-TRF",
                           `lfc6` = "SOTT-TRFTT")

df_fig_pair1.2$group = factor(df_fig_pair1.2$group,
                           levels = c("Cont-SO",
                                      "Cont-TRF",
                                      "TT-SOTT",
                                      "TT-TRFTT",
                                      "SO-TRF",
                                      "SOTT-TRFTT"))
lo = floor(min(df_fig_pair1.2$value))
up = ceiling(max(df_fig_pair1.2$value))
mid = (lo1 + up1)/2

fig_pair2 = df_fig_pair1.2 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold change") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_pair2


```


for KEGG using EC identifiers (KO)

```{r}

library(clusterProfiler)
library(phyloseq)
library(microbiomeMarker)

ec <- read.csv("newoutput.bacteria.ec_tally.csv", header = TRUE, row.names = 1)
meta
ec_descrip <- as.matrix(read.csv("tax_table_ec.csv", header = TRUE, row.names = 1))


otu_ec <- otu_table(ec, taxa_are_rows = FALSE)
tax_ec <- tax_table(ec_descrip)
mat = sample_data(meta)

#combined into a phyloseq object
ps_ec <-phyloseq(otu_ec, mat, tax_ec)
ps_ec

#lefse
lefse_ec <- run_lefse(ps_ec,
                   group = "Intervention",
                   taxa_rank = "Species",
                   kw_cutoff = 0.01,
                   multigrp_strat = TRUE)
lefse_ec

lefse_plot_pathway <- plot_ef_bar(lefse1)

limma_ec <- run_limma_voom(ps_ec,
                           group = "Intervention",
                           norm = "RLE",
                           taxa_rank = "Species",
                           p_adjust = "BH",
                           pvalue_cutoff = 0.05)
limma_ec


limma.dat <- as.data.frame(marker_table(limma_ec))
limma.dat <- as_tibble(limma.dat)

dat <- as.matrix(marker_table(limma_ec))
write.csv(dat, "lim_ec.csv")

limma.dat2 <- read.csv("lim_ec_edited.csv", row.names = 1)

limma.dat %>%
  as.character(limma.dat$enrich_group) %>%
  mutate(enrich_group = fct(reorder(enrich_group, ef_F_statistic)))

limma_plot <- ggplot(limma.dat2, aes(x = fct_reorder(feature, enrich_group), y = ef_F_statistic, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  coord_flip()
limma_plot + theme(text = element_text(family = "serif")) + scale_fill_manual(values = c("Cont" = "lightblue", "SO" = "orange",
                                                    "SOTT" = "seagreen", "TRF" = "firebrick1",
                                                    "TRFTT" = "purple", "TT" = "coral4"))




#clusterprofiler with lim_ec? will it work
limma.dat2

library(clusterProfiler)
search_kegg_organism("ece", by = "kegg_code")

ecoli <- search_kegg_organism('Escherichia coli', by='scientific_name')
dim(ecoli)

kk <- enrichKEGG(gene = limma.dat$feature, 
                 organism = "hsa",
                 pvalueCutoff = 0.05)
kk









adlx1 <- subset_samples(ps_ec, Intervention %in% c("Cont", "TRF"))

adlx <- run_aldex(adlx1,
                  group = "Intervention",
                  transform = "log10p",
                  pvalue_cutoff = 0.05,
                  p_adjust = "none", 
                  paired = FALSE)

adlx





output1 = ancombc2(ps_ec, assay_name = "counts", rank = "Species",
                  fix_formula = "Intervention", rand_formula = NULL,
                  p_adj_method = "holm", 
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "Intervention", struc_zero = FALSE, neg_lb = FALSE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = FALSE, pairwise = TRUE, 
                  dunnet = TRUE, trend = FALSE,
                  iter_control = list(tol = 1e-5, max_iter = 20, 
                                      verbose = FALSE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = NULL, 
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), 
                  trend_control = NULL)
output1

res2 <- output1$res_pair
write.csv(res2, "res_pathway.csv")
View(res2)

head(res2)

df_fig_pair1.1 = res2 %>%
  dplyr::filter(diff_InterventionVehicle_unvaccinated == 1 |
                  diff_InterventionTRF_unvaccinated == 1|
                  diff_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1 |
                  diff_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_InterventionVehicle_unvaccinated == 1,
                              round(lfc_InterventionVehicle_unvaccinated,2), 0),
                lfc2 = ifelse(diff_InterventionTRF_unvaccinated == 1,
                              round(lfc_InterventionTRF_unvaccinated, 2), 0),
                lfc3 = ifelse(diff_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1, 
                              round(lfc_InterventionVehicle_unvaccinated_InterventionControl_vaccinated, 2), 0),
                lfc4 = ifelse(diff_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1,
                              round(lfc_InterventionTRF_vaccinated_InterventionControl_vaccinated, 2), 0),
                lfc5 = ifelse(diff_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1,
                              round(lfc_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated, 2), 0),
                lfc6 = ifelse(diff_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1, 
                              round(lfc_InterventionVehicle_vaccinated_InterventionTRF_vaccinated, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6,
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_pair2.1 = res2 %>%
  dplyr::filter(diff_InterventionVehicle_unvaccinated == 1 |
                  diff_InterventionTRF_unvaccinated == 1|
                  diff_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1 |
                  diff_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1 |
                  diff_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1) %>%
  dplyr::mutate(lfc1 = ifelse(passed_ss_InterventionVehicle_unvaccinated ==1,
                              "aquamarine3", "black"),
                lfc2 = ifelse(passed_ss_InterventionTRF_unvaccinated == 1,
                              "aquamarine3", "black"),
                lfc3 = ifelse(passed_ss_InterventionVehicle_vaccinated_InterventionControl_vaccinated == 1, 
                              "aquamarine3", "black"),
                lfc4 = ifelse(passed_ss_InterventionTRF_vaccinated_InterventionControl_vaccinated == 1,
                              "aquamarine3", "black"),
                lfc5 = ifelse(passed_ss_InterventionVehicle_unvaccinated_InterventionTRF_unvaccinated == 1,
                              "aquamarine3", "black"),
                lfc6 = ifelse(passed_ss_InterventionVehicle_vaccinated_InterventionTRF_vaccinated == 1,
                              "aquamarine3", "black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc6,
                      names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_pair1.2 = df_fig_pair1.1 %>%
  dplyr::left_join(df_fig_pair2.1, by = c("taxon", "group"))

df_fig_pair1.2$group = recode(df_fig_pair1.2$group,
                           `lfc1` = "Cont-SO",
                           `lfc2` = "Cont-TRF",
                           `lfc3` = "TT-SOTT",
                           `lfc4` = "TT-TRFTT",
                           `lfc5` = "SO-TRF",
                           `lfc6` = "SOTT-TRFTT")

df_fig_pair1.2$group = factor(df_fig_pair1.2$group,
                           levels = c("Cont-SO",
                                      "Cont-TRF",
                                      "TT-SOTT",
                                      "TT-TRFTT",
                                      "SO-TRF",
                                      "SOTT-TRFTT"))
lo = floor(min(df_fig_pair1.2$value))
up = ceiling(max(df_fig_pair1.2$value))
mid = (lo1 + up1)/2

fig_pair2 = df_fig_pair1.2 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = "none") +
  labs(x = NULL, y = NULL, title = "Log fold change") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_pair2

```

SCFA

```{r}

dfSCFA <- read.csv("scfa_copy.csv")
dfSCFA


library(ggplot2)
library(tidyr)

#dfSCFA <- dfSCFA[apply(df!=0, 1, all),] #we only need this if any of the vectors are non-numeric
#summarise(dfSCFA)
# Reshaping the data for plotting
df_long <- dfSCFA %>% 
  pivot_longer(cols = c(Acetic_acid, Butyric_acid, Propionic_acid, Isobutyric_acid), 
               names_to = "acid_type", values_to = "value") 

library(dplyr)

dfSCFA
# Calculate mean and standard error for each acid grouped by intervention
my_comparisons <- list( c("Vehicle", "Vehicle_vaccinated"), 
                        c("Vehicle", "TRF"), 
                        c("Vehicle", "TRF_vaccinated"),
                        c("Vehicle","Control_vaccinated"),
                        c("Vehicle","Control"),
                        c("Vehicle_vaccinated","TRF"),
                        c("Vehicle_vaccinated","TRF_vaccinated"),
                        c("vehicle_vaccinated","Control_vaccinated"),
                        c("Vehicle_vaccinated","Control_vaccinated"),
                        c("TRF","TRF_vaccinated"),
                        c("TRF","Control_vaccinated"),
                        c("TRF","Control"),
                        c("TRF_vaccinated","Control_vaccinated"),
                        c("TRF_vaccinated","Control"),
                        c("Control_vaccinated","Control"))


# Display the summary
#print(df_summary)

# Plotting multiple barplots

df_long$Intervention <- as.factor(df_long$Intervention)


ggboxplot(df_long, x="Intervention", y="value", fill="acid_type", shape = "acid_type") +
  stat_compare_means(comparisons = my_comparisons) +
  facet_wrap(~acid_type) 

df_long$Intervention <- as.factor(df_long$Intervention)


#using compact letter display
library(multcompView)
library(multcomp)

stat.test <- df_long %>%
  group_by(acid_type) %>%
  tukey_hsd(value ~ Intervention) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  add_xy_position() 
stat.test 

stat.test2 <- df_long %>% #doesnt work?
  group_by(acid_type) %>%
  nest %>%
  rowwise %>%
  mutate(aov_result = list(aov(value ~ Intervention, data = df_long)), 
         tukey = list(TukeyHSD(aov_result)),
         cld = list(multcompLetters4(aov_result, tukey)$acid_type$letters)) %>%
  unnest(cld) %>% 
  select(acid_type, LETTER = cld) %>%
  mutate(acid_type = names(LETTER))
   

stat.test2 <- df_long %>% #doesnt work?
  group_by(acid_type) %>%
  TukeyHSD(aov(formula = value ~ Intervention))


anova.rr <- aov(value ~  acid_type + Intervention, data = df_long)
summary(anova.rr)

tukey.rr <- TukeyHSD(anova.rr)
cld.rr <- multcompLetters4(anova.rr, tukey.rr)
View(cld.rr)

sc1 <- ggboxplot(df_long,  x= "Intervention", y = "value", fill="acid_type", palette = "jco",
          ylab = "Log transformed concentration (µM)") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+
  stat_compare_means(method = "anova") +
  facet_wrap(~acid_type, scales = "free") + theme_bw()

View(stat.test2)
sc2 <- ggboxplot(df_long,  x= "Intervention", y = "value", fill = "acid_type", palette = "jco",
          ylab = "Log transformed concentration (µM)") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", step.increase = 0.1, hide.ns = TRUE, y.position = 4) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  facet_wrap(~acid_type, scales = "free") +
  theme_bw()

sc2 + theme(text = element_text(family = "serif"))



sc2 <- ggboxplot(df_long,  x= "Intervention", y = "value", fill = "acid_type", palette = "jco",
          ylab = "Log transformed concentration (µM)") +
  facet_wrap(~acid_type, scales = "free") +
  theme_bw()

sc2 + stat_pvalue_manual(stat.test, label = "p.adj.signif", hide.ns = TRUE, tip.length = 0.01)

```



# Using the DIABLO framework

```{r}

immune2 <- read.csv("immune3.csv", row.names = 1)
immune2 <- data.matrix(immune2, rownames.force = 1)
asv_imm <- read.csv("asv_subset_copy2.csv", row.names = 1)
asv_imm <- data.matrix(asv_imm, rownames.force = 1)
scfa_imm <- read.csv("scfa_subset.csv", row.names = 1)
scfa_imm <- data.matrix(scfa_imm, rownames.force = 1)
meta_sub <- read.csv("meta_subset.csv", row.names = 1)

meta_sub$Intervention <- as.factor(meta_sub$Intervention)

diablo_imm <- list(asv = asv_imm,
                   imm = immune2,
                   scfa = scfa_imm,
                   meta = meta_sub$Intervention)


data <- list(asv = asv_imm,
             imm = immune2,
             scfa = scfa_imm)

Y = diablo_imm$meta # set the response variable as the Y df
summary(Y)

#selecting the arbitrary values and how do I do it? apparently the numbers represent top 25 features
list.keepX = c(25,25)
list.keepY = c(25,25)
list.keepY_scfa = c(4,4)

#generate three pairwise PLS models.Not sure if i have to logratio it
pls1 <- spls(data[["asv"]], data[["imm"]], 
             keepX = list.keepX, keepY = list.keepY) 
pls2 <- spls(data[["asv"]], data[["scfa"]], 
             keepX = list.keepX, keepY = list.keepY_scfa)
pls3 <- spls(data[["imm"]], data[["scfa"]], 
             keepX = list.keepX, keepY = list.keepY_scfa)


plotVar(pls1, cutoff = 0.5, title = "(a) Genes vs asv", 
        legend = c("genes", "Asv"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

plotVar(pls2, cutoff = 0.5, title = "(a) genes vs scfa", 
        legend = c("genes", "ASVs"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

plotVar(pls3, cutoff = 0.5, title = "(a) ASV vs scfa ", 
        legend = c("Asv", "Scfa"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))


cor(pls1$variates$X, pls1$variates$Y) 
cor(pls2$variates$X, pls2$variates$Y) 
cor(pls3$variates$X, pls3$variates$Y)


#designing the matrix unsure how


design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s

design

basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 4, design = design) 

perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 3, nrepeat = 10) 
plot(perf.diablo)

# set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 

test.keepX = list(asv = c(5:9, seq(10, 18, 2), seq(30,40,5)), 
                  imm = c(5:9, seq(10, 18, 2), seq(30,40,5)),
                  scfa = c(4:4, seq(10, 18, 2), seq(30,40,5)))

BPPARAM <- BiocParallel::MulticoreParam(workers = parallel::detectCores()-1) #not supported on windows apparently

snow <- BiocParallel::SnowParam(workers  = 8, type = "SOCK")

tune = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 3, nrepeat = 10,
                              dist = "centroids.dist",
                              BPPARAM = snow)
list.keepX = tune$choice.keepX # set the optimal values of features to retain
list.keepX


final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                          keepX = list.keepX, design = design)

final.diablo.model$design
#selecting the blocks
selectVar(final.diablo.model, block = 'asv', comp = 1)$asv$name 
selectVar(final.diablo.model, block = 'imm', comp = 1)$imm$name 

plotDiablo(final.diablo.model, ncomp = 2)

plotVar(final.diablo.model, var.names = FALSE, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 17, 15), cex = c(2,2,2), 
        col = c('darkorchid', 'brown1', 'lightgreen'))


circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1, size.variables = 0.6)


plotLoadings(final.diablo.model, comp = 1, contrib = 'max', method = 'median') 
 

cimDiablo(final.diablo.model, margins = c(10,13), legend.position = "right", size.legend = 0.75)  + theme(text = element_text(family = "serif"))


#checking the performance to see the error rate
perf.diablo = perf(final.diablo.model, validation = 'Mfold', 
                   M = 3, nrepeat = 10, 
                   dist = 'centroids.dist') 

perf.diablo$MajorityVote.error.rate

perf.diablo$WeightedVote.error.rate


```