setwd("~/q2xesto/454/seqs/phyloseq")
library("ggplot2")
library("phyloseq")
library("ape")
library("plyr")
library("dplyr")
library("vegan")
library("RColorBrewer")
library("svglite")

otu <- read.table(file = "otu_table.txt", header = TRUE, check.names=FALSE)
tax <- read.table(file = "taxonomy.tsv", sep = '\t', header = TRUE)
merged <- merge(otu, tax, by.x = c("OTUID"), by.y=c("OTUID"))
write.table(merged, file = "combined_otu_tax.tsv", sep = '\t', col.names = TRUE, row.names = FALSE)
otu_table = read.csv("otu_matrix.csv", sep=",", row.names=1,check.names=FALSE)
otu_table = as.matrix(otu_table)
taxonomy = read.csv("taxonomy.csv", sep=",", row.names=1)
taxonomy = as.matrix(taxonomy)
META<-import_qiime_sample_data("xesto454_metadata_for_R.tsv")
META
TREE=read_tree("tree.nwk")
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
taxa_names(TAX)
taxa_names(OTU)
sample_names(OTU)
sample_names(META)
xestops = phyloseq(OTU, TAX, META, TREE)

#First we'll Clean the unknowns then check how the mock communities compare to the Zymo standards
tax.clean <- data.frame(tax_table(xestops))
for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
  }
}
tax_table(xestops) <- as.matrix(tax.clean)

#Core microbiome
control <-subset_samples(xestops,Description == "healthy_no_disease")
diseased <-subset_samples(xestops,Description == "diseased_tissue")
HoD <-subset_samples(xestops,Description == "healthy_portion_diseased")
disease_band <-subset_samples(xestops,Description == "orange_band_tissue")

control_sub = subset_samples(control, Alias != "5_11_137H" )
TopNOTUs_control <- names(sort(taxa_sums(control_sub), TRUE)[1:20])
control20  <- prune_taxa(TopNOTUs_control, control_sub)

control_df<-psmelt(control_sub)
controlplot=ggplot(control_df, aes(fill=Phylum, y=Abundance,x=Sample)) + theme_classic() + geom_bar(stat="identity", position="fill")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
controlplot
controlplot<-controlplot + scale_fill_manual(values=cbb_pal_expand)
controlplot
png('xesto_control.png', units='in',width=7.5,height=5.5,res=600)
controlplot
dev.off()
ggsave(file="xesto_control.svg", plot=controlplot, dpi = 600)

#Set the palette
theme_set(theme_bw())
pal = "Dark2"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

all_rich<-plot_richness(xestops, x="Description", color="Location",measures=c("Observed","Simpson"))
all_rich
all_rich<-all_rich+geom_point(size=5)
ggsave(file="all_rich.svg", plot=all_rich, dpi = 600)
png('all_rich.png', units='in',width=5.68,height=3.62,res=600)
all_rich
dev.off()

#Now for the real data
#This alpha diversity needs to be done with the non-relative abundance data ie, the original phyloseq object
ALPHA<-estimate_richness(xestops, measures = c("Observed", "Shannon", "Simpson","Fisher"))

#Alpha Diversity for all samples
DATAFRAME<-data.frame(ALPHA,sample_data(xestops))
model<-lm(DATAFRAME$Simpson~DATAFRAME$Description)
summary(model)
s <- summary(model)
capture.output(s, file = "model1_output_xesto_simp__desc.txt")
model2<-lm(DATAFRAME$Simpson~1)
summary(model2)
s2 <- summary(model2)
capture.output(s2, file = "model2_output_xesto_simp_desc.txt")
anova(model,model2)
a<-summary(anova(model,model2))
a
capture.output(a, file = "anova_output_xesto_simp_desc.txt")

model_obs<-lm(DATAFRAME$Observed~DATAFRAME$Description)
summary(model_obs)
s_obs <- summary(model_obs)
capture.output(s_obs, file = "model1_output_xesto_obs.txt")
model2_obs<-lm(DATAFRAME$Observed~1)
summary(model2_obs)
s2_obs <- summary(model2_obs)
capture.output(s2_obs, file = "model2_output_xesto_obs.txt")
anova(model_obs,model2_obs)
a_obs<-summary(anova(model_obs,model2_obs))
capture.output(a_obs, file = "anova_output_xesto_obs.txt")

model_shannon<-lm(DATAFRAME$Shannon~DATAFRAME$Description)
summary(model_shannon)
s_shannon <- summary(model_shannon)
capture.output(s_shannon, file = "model1_output_xesto_shannon.txt")
model2_shannon<-lm(DATAFRAME$Shannon~1)
summary(model2_shannon)
s2_shannon <- summary(model2_shannon)
capture.output(s2_shannon, file = "model2_output_xesto_shannon.txt")
anova(model_shannon,model2_shannon)
a_shannon<-summary(anova(model_shannon,model2_shannon))
capture.output(a_shannon, file = "anova_output_xesto_shannon.txt")

model_fisher<-lm(DATAFRAME$Fisher~DATAFRAME$Description)
summary(model_fisher)
s_fisher <- summary(model_fisher)
capture.output(s_fisher, file = "model1_output_xesto_fisher.txt")
model2_fisher<-lm(DATAFRAME$Fisher~1)
summary(model2_fisher)
s2_fisher <- summary(model2_fisher)
capture.output(s2_fisher, file = "model2_output_xesto_fisher.txt")
anova(model_fisher,model2_fisher)
a_fisher<-summary(anova(model_fisher,model2_fisher))
capture.output(a_fisher, file = "anova_output_xesto_fisher.txt")


#Transforms into relative abundance
xestops_rel = transform_sample_counts(xestops, function(x) x / sum(x) )

#Beta Diversity and Ordination
##Ordination Plots
theme_set(theme_classic())
pal = "Dark2"
scale_colour_discrete <-  function(palname=pal, ...){
  scale_colour_brewer(palette=palname, ...)
}
scale_fill_discrete <-  function(palname=pal, ...){
  scale_fill_brewer(palette=palname, ...)
}

xesto.ord <- ordinate(xestops_rel, "NMDS", "bray")
g1 = plot_ordination(xestops_rel, xesto.ord, type="samples", color="Description")
g1
g2<-g1+ geom_point(size=3)
g2
png('xesto_ord_nonphylo_nmds_bray.png', units='in',width=5.68,height=3.62,res=600)
g2
dev.off()
ggsave(file="xesto_ord_nonphylo_nmds_bray.svg", plot=g2, dpi = 600)

#Phylogenetic-unweighted Unifrac

xesto.phylo.unw.ord <- ordinate(xestops_rel, "PCoA", "unifrac",weighted=FALSE)
g1 = plot_ordination(xestops_rel, xesto.phylo.unw.ord, type="samples", color="Description")
g1
g2<-g1+ geom_point(size=3)
g2
png('xesto_ord_unweighted-unifrac.png', units='in',width=5.68,height=3.62,res=600)
g2
dev.off()
ggsave(file="xesto_ord_unweighted-unifrac.svg", plot=g2, dpi = 600)

#Phylogenetic-weighted Unifrac

xesto.phylo.weight.ord <- ordinate(xestops_rel, "PCoA", "unifrac",weighted=TRUE)
g1 = plot_ordination(xestops_rel, xesto.phylo.weight.ord, type="samples", color="Description")
g1
g2<-g1+ geom_point(size=3)
g2
png('xesto_ord_weighted-unifrac.png', units='in',width=5.68,height=3.62,res=600)
g2
dev.off()
ggsave(file="xesto_ord_weighted-unifrac.svg", plot=g2, dpi = 600)

#Beta Diversity Analysis
##non-phylogenetic beta diversity metrics
set.seed(52)
all_bray <- phyloseq::distance(xestops_rel, method = "bray")
xestodf <- data.frame(sample_data(xestops_rel))
adonis_xesto_desc<-adonis(formula = all_bray ~ Description, data = xestodf)
adonis_xesto_desc
capture.output(adonis_xesto_desc, file="ADONIS_bray_desc.csv")

#Tests for homogeneity
beta_xesto_desc <- betadisper(all_bray, xestodf$Description)
beta_xesto_desc
capture.output(beta_xesto_desc, file="bdisper_bray_xesto_desc.csv")
pertest_desc<-permutest(beta_xesto_desc)
capture.output(pertest_desc, file="pertest_bray_xesto_desc.csv")

#phylogenetic beta diversity metrics-unweighted
all_unifrac_unw <- phyloseq::distance(xestops_rel, method = "unifrac")
xesto_phylodf <- data.frame(sample_data(xestops_rel))
adonis_xesto_desc_phylo_unw<-adonis(formula = all_unifrac_unw ~ Description, data = xesto_phylodf)
adonis_xesto_desc_phylo_unw
capture.output(adonis_xesto_desc_phylo_unw, file="ADONIS_unifrac_rel_desc_unw.csv")

#Tests for homogeneity-phylogentic
beta_xesto_uni_desc <- betadisper(all_unifrac_unw, xesto_phylodf$Description)
beta_xesto_uni_desc
capture.output(beta_xesto_uni_desc, file="bdisper_unifrac_xesto_desc_unw.csv")
pertest_uni_desc<-permutest(beta_xesto_uni_desc)
capture.output(pertest_uni_desc, file="pertest_unifrac_xesto_desc_unw.csv")

#phylogenetic beta diversity metrics-weighted
all_unifrac <- phyloseq::distance(xestops_rel, method = "wunifrac")
xesto_phylodf <- data.frame(sample_data(xestops_rel))
adonis_xesto_desc_phylo<-adonis(formula = all_unifrac ~ Description, data = xesto_phylodf)
adonis_xesto_desc_phylo
capture.output(adonis_xesto_desc_phylo, file="ADONIS_unifrac_rel_desc_weight.csv")

#Tests for homogeneity-phylogentic
beta_xesto_uni_desc <- betadisper(all_unifrac, xesto_phylodf$Description)
beta_xesto_uni_desc
capture.output(beta_xesto_uni_desc, file="bdisper_unifrac_xesto_desc_weight.csv")
pertest_uni_desc<-permutest(beta_xesto_uni_desc)
capture.output(pertest_uni_desc, file="pertest_unifrac_xesto_desc_unw.csv")

#Top taxonomy-top 50
TopNOTUs <- names(sort(taxa_sums(xestops_rel), TRUE)[1:50])
xesto50   <- prune_taxa(TopNOTUs, xestops_rel)


xesto50_df<-psmelt(xesto50)
level_order50 <- factor(xesto50_df$Description,level = c('healthy_no_disease','healthy_portion_diseased','orange_band_tissue','diseased_tissue','sediment','seawater'))
xestoplot50=ggplot(xesto50_df, aes(fill=Genus, y=Abundance, x=level_order50)) + theme_classic() + geom_bar(stat="identity", position="fill")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
xestoplot50
xestoplot<-xestoplot50 + scale_fill_manual(values=cbb_pal_expand)
xestoplot
png('xesto_top50.png', units='in',width=7.5,height=5.5,res=600)
xestoplot
dev.off()
ggsave(file="xesto_top50.svg", plot=xestoplot, dpi = 600)

#ANOVA on alpha diversity
ALL_DIV<-estimate_richness(merge_rare, measures = c("Observed","Simpson"))
GDNA_DIV<-estimate_richness(gDNA, measures = c("Observed","Simpson"))
CDNA_DIV<-estimate_richness(cDNA, measures = c("Observed","Simpson"))
all_na<-sample_data(merge_ns)$NA_type
all_tp<-sample_data(merge_ns)$Timepoint
gDNA_tp<-sample_data(gDNA)$Timepoint
cDNA_tp<-sample_data(cDNA)$Timepoint
summary(aov(ALL_DIV$Observed~all_na))
summary(aov(ALL_DIV$Simpson~all_na))
summary(aov(ALL_DIV$Observed~all_tp))
summary(aov(ALL_DIV$Simpson~all_tp))
summary(aov(GDNA_DIV$Observed~gDNA_tp))
summary(aov(GDNA_DIV$Simpson~gDNA_tp))
summary(aov(CDNA_DIV$Observed~cDNA_tp))
summary(aov(CDNA_DIV$Simpson~cDNA_tp))
rarecurve(t(otu_table(merge_ns)),step=100)
rarecurve(t(otu_table(merge_rare)),step=100)
rarecurve(t(otu_table(gDNA)),step=100)
rarecurve(t(otu_table(cDNA)),step=100)
