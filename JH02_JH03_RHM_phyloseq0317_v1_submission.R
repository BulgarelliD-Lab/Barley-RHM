#############################################################
#
# Ref to the ARTICLE
# 
# Senga Robertson & Davide Bulgarelli
# s.u.robertson@dundee.ac.uk
# d.bulgarelli@dundee.ac.uk
# 
# Revison March 2017
# 
# script to reproduce calculations and figures presented in the manuscript
# 
# Disclaimer: the manuscript is currently submitted, the script might be subjected to changes derived from the revision process
#
#############################################################

#############################################################
# Clean-up the memory and start a new session
#############################################################

rm(list=ls())
dev.off()

#############################################################
# Libraries required
#############################################################

#1st time installation
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
#biocLite("DESeq2")
#biocLite("PMCMR")

#required packages 
library("phyloseq")
library("DESeq2")
library("ggplot2")
library("vegan")
library ("ape")
library("PMCMR")
library("plyr")


#set the working directory


#############################################################
#import the count matrix and the desing file
#############################################################


#OTU table generated using QIIME 1.9.0. 
dat_info <- read.delim("JH02_JH03_RHM_otu_table_nc2.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#inspect the file 
dim(dat_info)
colnames(dat_info)

#extract the total number of reads clustered at OTU 97% identiy (the number underneath the Hv identifier represents the total number of reads clustered for that sample) 
OTU_97_reads <- sort(colSums(dat_info[, 1:55]))
OTU_97_reads

#total reads
OTU_97_reads_sum <- sum(colSums(dat_info[, 1:55]))
OTU_97_reads_sum

#design file
design <- read.delim("Map_JH02_JH03_RHM_phyloseq.txt", sep = "\t", header=TRUE, row.names=1)
design

#remove chloroplast and mitochondria OTUs from the original dataset
Chloroplast <- dat_info[grepl("Chloroplast", dat_info$ConsensusLineage), ]
dim(Chloroplast)

mitochondria <- dat_info[grepl("mitochondria", dat_info$ConsensusLineage), ]
dim(mitochondria)

#set a difference between the row names of the the three datasets: this information will be used to filter out Plant derived OTUs from the OTU table
noPlants <- setdiff(rownames(dat_info), c(rownames(Chloroplast), rownames(mitochondria)))

#inspect the results
length(rownames(dat_info))
length(noPlants)

#save the OTUids list generated at line 81. This will be used to usbset the OTU table and generate taxa-tables in QIIME
#write(noPlants, "JH02_JH03_RHM_noPlant_OTUs_id.txt")
#this corresponds to Worksheet ws3 in the Supplementary Database 1

#generate a new OTU table which will be devoid of Chloroplast and Mitochondria OTUs
dat_info_noPlants <- dat_info[noPlants, ]

#create a new count matrix without OTUs assigned to Choloplast and Mitochondria
dat_count <- dat_info[, rownames(design)]
dat_count_noplants <- dat_info[noPlants, rownames(design)]
dim(dat_count_noplants)

#and a new taxa table
dat_tax_noPlants <- as.data.frame(dat_info[rownames(dat_count_noplants), 56])
rownames(dat_tax_noPlants) <- rownames(dat_count_noplants)
#save the above file and in excel we will create a new tax table where each column represents a taxonomic rank
#write.table(dat_tax_noPlants, file="JH02_JH03_RHM_dat_tax_noPlants.txt", sep="\t")

#check the effect of mitochondria/chloroplast depletion on the new OTU table
#with plant sequences
dim(dat_count)

#w/o plant-derived sequences
dim(dat_count_noplants)

#total number of reads w/o plant sequences
OTU_97_reads_noPlants <- colSums(dat_count_noplants)
OTU_97_reads_noPlants

#now sort per sample
sort(OTU_97_reads_noPlants)

#total number of reads
OTU_97_reads_noPlants_sum <- sum(OTU_97_reads_noPlants)
OTU_97_reads_noPlants_sum 

#Define the proportion of non-plant reads in the original dataset 
useful_reads <- (OTU_97_reads_noPlants_sum/OTU_97_reads_sum)*100
useful_reads

#create a dataset to visualise the proportion of reads per sample before and after removing OTUs assigned to chloroplast and mitochondria
OTU_97_reads_noPlants <- as.data.frame(OTU_97_reads_noPlants)
OTU_97_microbial_reads_proportion <- as.data.frame(colSums(dat_count_noplants)/colSums(dat_count))*100

#rename the columns in the generated datasets
colnames(OTU_97_reads_noPlants) <- c("reads")
colnames(OTU_97_microbial_reads_proportion) <- c("Microbial_OTUs_reads")

#combine these datasets with the design file
design_info <- cbind(design, OTU_97_reads_noPlants)
design_info_2 <- cbind(design_info, OTU_97_microbial_reads_proportion)

#Create a new dataset with working samples. This is dataset Dataset_1
working_samples <- rownames(design_info_2)[which(design_info_2$Treatment!= "Control")]
Dataset_1 <- design_info_2[working_samples, ]
dim(Dataset_1)

#Use the working_samples id to define the proportion of useful reads in the experiment
dat_count_working_samples <- dat_count[, working_samples]
dat_count_noPlants_working_samples <- dat_count_noplants[, working_samples]
#total number of reads
OTU_97_working_samples_reads_sum <- sum(colSums(dat_count_working_samples))
OTU_97_working_samples_reads_sum 
OTU_97_working_samples_noPlants_reads_sum <- sum(colSums(dat_count_noPlants_working_samples))
OTU_97_working_samples_noPlants_reads_sum
useful_reads_working_samples <- OTU_97_working_samples_noPlants_reads_sum/OTU_97_working_samples_reads_sum *100
useful_reads_working_samples

#calculate the max, min and mean number of reads for the dataset
mean(Dataset_1$reads)
max(Dataset_1$reads)
min(Dataset_1$reads)

#Create a new dataset with rhizosphere samples. This is dataset Dataset_2
rhizo_samples <- rownames(Dataset_1)[which(Dataset_1$Microhabitat!= "Soil")]
Dataset_2 <- Dataset_1[rhizo_samples, ] 
dim(Dataset_2)

##############################################################
#Figure 1A: DNA concentration calculation, rhizosphere samples
#Data required: Dataset 2
##############################################################

#re-order the factors
Dataset_2$Description <- ordered(Dataset_2$Description, levels=c("Karat.quarry", "RHL1.quarry", "Dema.quarry", "RHP1.quarry", 
                                                                 "Karat.tayport", "RHL1.tayport", "Dema.tayport", "RHP1.tayport"))

#Figure S2
with(Dataset_2, boxplot(DNAConc ~ Description, xlab = "Genotype", ylab = "ng/ul",   main = "DNA concentration"))

#test the normality of the observed data
shapiro.test(Dataset_2$DNAConc)

#Test the soil effect on DNA concentration
wilcox.test(DNAConc ~ Soil, data = Dataset_2)

#Split the dataset and repeat the analysis independently for soil type
quarry_samples <- rownames(Dataset_2)[which(Dataset_2$Soil== "Quarryfield")]
tayport_samples <- rownames(Dataset_2)[which(Dataset_2$Soil== "Tayport")]

Dataset_2_quarryfield <- Dataset_2[quarry_samples, ]
Dataset_2_tayport <- Dataset_2[tayport_samples, ]

#Stats soil corrected
#Quarryfield
kruskal.test(DNAConc ~ Genotype, data = Dataset_2_quarryfield)
posthoc.kruskal.dunn.test (x=Dataset_2_quarryfield$DNAConc, g=Dataset_2_quarryfield$Genotype, p.adjust.method="BH")

#Tayport
kruskal.test(DNAConc ~ Genotype, data = Dataset_2_tayport)
posthoc.kruskal.dunn.test (x=Dataset_2_tayport$DNAConc, g=Dataset_2_tayport$Genotype, p.adjust.method="BH")

#############################################################
#Figure 1B: reads distribution averall and whitin rhizosphre samples
#Data required: Dataset 2
#############################################################

#data visualisation: rhizosphere samples only
with(Dataset_2, boxplot(Microbial_OTUs_reads ~ Description, xlab = "Samples", ylab = " sequencing reads proportion",   main = "Reads assigned to Microbial OTUs"))

#test the normality of the observed data
shapiro.test(Dataset_2$Microbial_OTUs_reads)

#use two-way anova to asses significant differences between samples
Microbial_OTUs_stats <- aov(Dataset_2$Microbial_OTUs_reads ~ Soil * Genotype, data = Dataset_2)
summary(Microbial_OTUs_stats)

#############################################################
#Figure 1C: Dry weight calculation
#Data required: Dataset 2
#############################################################

#Data visualisation
with(Dataset_2, boxplot(Dryweight ~ Description, xlab = "Genotype", ylab = "dry weight g",   main = "Aboveground biomass"))

#test the normality of the observed data
shapiro.test(Dataset_2$Dryweight)

#Test the soil effect on aboveground biomass
wilcox.test(Dryweight ~ Soil, data = Dataset_2)

#Split the dataset and repeat the analysis independently for soil type
quarry_samples <- rownames(Dataset_2)[which(Dataset_2$Soil== "Quarryfield")]
tayport_samples <- rownames(Dataset_2)[which(Dataset_2$Soil== "Tayport")]

Dataset_2_quarryfield <- Dataset_2[quarry_samples, ]
Dataset_2_tayport <- Dataset_2[tayport_samples, ]

#Stats
#Quarryfield
kruskal.test(Dryweight ~ Genotype, data = Dataset_2_quarryfield)
posthoc.kruskal.dunn.test (x=Dataset_2_quarryfield$Dryweight, g=Dataset_2_quarryfield$Genotype, p.adjust.method="BH")

#Tayport
kruskal.test(Dryweight ~ Genotype, data = Dataset_2_tayport)
posthoc.kruskal.dunn.test (x=Dataset_2_tayport$Dryweight, g=Dataset_2_tayport$Genotype, p.adjust.method="BH")

#############################################################
#Figure 2: High taxonomic ranks distribution
#Data required: Description_otu_table_L2 
#############################################################

#Import average value phylum level
dat_info_taxa_Phylum <- read.delim("Description_otu_table_L2.txt", skip=1, sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)
colnames(dat_info_taxa_Phylum)
dim(dat_info_taxa_Phylum)

#identify the samples for the analysis
RHM_10_samples <- c("Bulk.quarry", "Karat.quarry", "RHL1.quarry", "Dema.quarry", "RHP1.quarry",
                "Bulk.tayport", "Karat.tayport", "RHL1.tayport", "Dema.tayport", "RHP1.tayport")

#create a new dataset containing only the 10 experimental samples
dat_info_taxa_Phylum_2 <- dat_info_taxa_Phylum[ ,RHM_10_samples]
colnames(dat_info_taxa_Phylum_2)

# transform the data in % for visualisation
dat_count_taxa_Phylum <- dat_info_taxa_Phylum_2
dat_norm_Phylum <- ((dat_count_taxa_Phylum)/colSums(dat_count_taxa_Phylum,na=T)) * 100 
dat_norm_Phylum[1:5, ]
colSums(dat_norm_Phylum)

# determine the average % of reads for each phylum
Phylum_mean_sorted <- dat_norm_Phylum[(order(-rowSums(dat_norm_Phylum))), ] 

#Identifythe phyla whose average abundance in above 1% and their aggregated relative abundance
Phylum_mean_topRank <- Phylum_mean_sorted[rownames(Phylum_mean_sorted)[which(rowMeans(Phylum_mean_sorted) > 1)], ]
dim(Phylum_mean_topRank)
colSums(Phylum_mean_topRank)

#overall
mean(colSums(Phylum_mean_topRank))

#individual phyla aggregated mean, corrected for microhabitat
#Bulk samples
Phylum_mean_topRank_Bulk <- Phylum_mean_topRank[, (c("Bulk.quarry", "Bulk.tayport"))]
rowMeans(Phylum_mean_topRank_Bulk)

#Rhizosphere samples
Phylum_mean_topRank_Rhizosphere <- Phylum_mean_topRank[, (setdiff(colnames(Phylum_mean_topRank), colnames(Phylum_mean_topRank_Bulk)))]
rowMeans(Phylum_mean_topRank_Rhizosphere)

#transform in matrix for plotting purposes
Phylum_mean_topRank <- as.matrix(Phylum_mean_topRank) 
colnames(Phylum_mean_topRank) 

# Stacked Bar Plot with Colors and Legend
barplot(Phylum_mean_topRank, main="Phylum Distribution",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkblue","red", "green", "cyan", "orange", "brown", "black", "magenta", "grey", "yellow", "purple"), beside=FALSE,   legend = rownames(Phylum_mean_topRank))

#Due to size limits the legend covers part of the graph, save the graph as .eps file and in illustrator uou can adjust this (and other)  graphical issues
#barplot_no legend: you can save this as image and then reconstruct the leged using the previous figure. Note: the order of samples can be inferred from the command on line 166
barplot(Phylum_mean_topRank, main="Phylum Distribution",
        xlab="Samples", ylab = "Reads %", ylim = c(0,100), col=c("darkblue","red", "green", "cyan", "orange", "brown", "black", "magenta", "grey", "yellow", "purple"), beside=FALSE)

#############################################################
#Genererate the phyloseq object
#Data required: dat_count_noplants; design, JH02_JH03_dat_tax_noPlants_ordered.txt, and 97_otus.tree.gz
#############################################################

#The OTU Table counts
JH02_JH03_RHM_OTU <- otu_table(dat_count_noplants, taxa_are_rows=TRUE)

#The taxonomy information
#Note that the file JH02_JH03_dat_tax_noPlants_ordered.txt has been generated from the output of lines 96-100  
#it is a tab-delimited file with 8 columns, the column headers are: OTU id; "Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species",
JH02_JH03_RHM_taxa_ordered <- read.delim ("JH02_JH03_RHM_dat_tax_noPlants_ordered.txt", sep = "\t", row.names=1, header=T, blank.lines.skip = FALSE)
JH02_JH03_RHM_taxa <- tax_table(as.matrix(JH02_JH03_RHM_taxa_ordered))
dim(JH02_JH03_RHM_taxa)

#The mapping file 
JH02_JH03_RHM_map <- sample_data(design)

#The phylogenetic tree: the OTU table has been generated using a closed reference approach agains the greengenes 13_05 database use the corresponding phylogenetic tree 
JH02_JH03_RHM_tree <- read_tree_greengenes("97_otus.tree.gz")

#check whether the tree is rooted
is.rooted(JH02_JH03_RHM_tree)

#merge the files and create the phyloseq object
JH02_JH03_RHM_data_phyloseq <- merge_phyloseq(JH02_JH03_RHM_OTU, JH02_JH03_RHM_taxa, JH02_JH03_RHM_map,  JH02_JH03_RHM_tree)

#inspect the generated data
JH02_JH03_RHM_data_phyloseq
sum(colSums(otu_table(JH02_JH03_RHM_data_phyloseq)))
dim(dat_count_noplants)
sum(colSums(dat_count_noplants))

#remove the control samples from the dataset
JH02_JH03_RHM_data_phyloseq_2 <- subset_samples(JH02_JH03_RHM_data_phyloseq, Treatment=="Sample")
design_2 <- design[colnames(otu_table(JH02_JH03_RHM_data_phyloseq_2)), ]
JH02_JH03_RHM_data_phyloseq_2 

#############################################################
#Figure 3: Alphadiversity calculations
#Data required: design_2; JH02_JH03_RHM_data_phyloseq_rare_table_counts_2.txt
#############################################################

#rarefy the dataset
#JH02_JH03_RHM_data_phyloseq_rare <- rarefy_even_depth(JH02_JH03_RHM_data_phyloseq_2, rngseed=TRUE)

#extract and save the OTU table for reproducibility of the code
#JH02_JH03_RHM_data_phyloseq_rare_table <- as.data.frame(otu_table(JH02_JH03_RHM_data_phyloseq_rare))

#inspect the generated file
#class(JH02_JH03_RHM_data_phyloseq_rare_table)
#dim(JH02_JH03_RHM_data_phyloseq_rare_table)

#save the file for the reproducibility of the code
#write.table(JH02_JH03_RHM_data_phyloseq_rare_table, file="JH02_JH03_RHM_data_phyloseq_rare_table_counts.txt", sep="\t")

#import the rarefied OTU counts (note file name counts2.txt this file has been generated in excel and includes the #OTU ID as header of the first column)
dat_count_rare <- read.delim("JH02_JH03_RHM_data_phyloseq_rare_table_counts_2.txt", sep = "\t", header=T, row.names= 1, blank.lines.skip = FALSE)

#inspect the generated file
dim(dat_count_rare)
colSums(dat_count_rare)

#generate a new phyloseq object wich will contain only the rarefied counts and the design file (only these two pieces of information are required for alphadiversity calculation)
JH02_JH03_RHM_OTU_rare <- otu_table(dat_count_rare, taxa_are_rows=TRUE)
JH02_JH03_RHM_map_rare <- sample_data(design_2)
JH02_JH03_RHM_data_rare_phyloseq <- merge_phyloseq(JH02_JH03_RHM_OTU_rare, JH02_JH03_RHM_map_rare)

#Inspect the generated file
JH02_JH03_RHM_data_rare_phyloseq
sample_sums(JH02_JH03_RHM_data_rare_phyloseq)

#Index calculations
JH02_JH03_RHM_alpha_rare <-  estimate_richness(JH02_JH03_RHM_data_rare_phyloseq, measures = c("Observed", "Shannon"))

#generate a new dataframes for data visualisation

#Sample information

#Genotype
design_genotype <- as.data.frame(design_2[, 3])
rownames(design_genotype) <- rownames(design_2)
colnames(design_genotype) <- c("Genotype")

#Sample
design_description <- as.data.frame(design_2[, 8])
rownames(design_description) <- rownames(design_2)
colnames(design_description) <- c("Description")

#data frame Genotype_Description
design_GD <- cbind(design_genotype, design_description)

#Observed OTUs
JH02_JH03_RHM_alpha_rare_Observed <- as.data.frame(JH02_JH03_RHM_alpha_rare[ ,1])
rownames(JH02_JH03_RHM_alpha_rare_Observed) <- rownames(JH02_JH03_RHM_alpha_rare)
colnames(JH02_JH03_RHM_alpha_rare_Observed) <- c("Observed")

#Combine the dataset sample description and Observed OTUs
JH02_JH03_RHM_alpha_rare_Observed_GD <- cbind(design_GD, JH02_JH03_RHM_alpha_rare_Observed)

#Order the levels according to a defined order
JH02_JH03_RHM_alpha_rare_Observed_GD$Description <- ordered(JH02_JH03_RHM_alpha_rare_Observed_GD$Description, levels=c("Bulk.quarry", "Karat.quarry", "RHL1.quarry", "Dema.quarry", "RHP1.quarry",
                                                                                                                                         "Bulk.tayport", "Karat.tayport", "RHL1.tayport", "Dema.tayport", "RHP1.tayport"))
#plotting
p <- ggplot(JH02_JH03_RHM_alpha_rare_Observed_GD, aes(x=Description, y=Observed)) + geom_point(aes(colour = factor(Genotype)))
p + scale_colour_manual(values = c("brown","orange", "cyan", "darkblue", "magenta"))

#Shannon
JH02_JH03_RHM_alpha_rare_Shannon <- as.data.frame(JH02_JH03_RHM_alpha_rare[ ,2])
rownames(JH02_JH03_RHM_alpha_rare_Shannon) <- rownames(JH02_JH03_RHM_alpha_rare)
colnames(JH02_JH03_RHM_alpha_rare_Shannon) <- c("Shannon")

#Combine the dataset sample description and Shannon OTUs
JH02_JH03_RHM_alpha_rare_Shannon_GD <- cbind(design_GD, JH02_JH03_RHM_alpha_rare_Shannon)

#Order the levels according to a defined order
JH02_JH03_RHM_alpha_rare_Shannon_GD$Description <- ordered(JH02_JH03_RHM_alpha_rare_Shannon_GD$Description, levels=c("Bulk.quarry", "Karat.quarry", "RHL1.quarry", "Dema.quarry", "RHP1.quarry",
                                                                                                                 "Bulk.tayport", "Karat.tayport", "RHL1.tayport", "Dema.tayport", "RHP1.tayport"))
#plotting
p <- ggplot(JH02_JH03_RHM_alpha_rare_Shannon_GD, aes(x=Description, y=Shannon)) + geom_point(aes(colour = factor(Genotype)))
p + scale_colour_manual(values = c("brown","orange", "cyan", "darkblue", "magenta"))

#generate a new dataframe for statistical analysis
JH02_JH03_RHM_alpha_rare_info <- cbind(design_2, JH02_JH03_RHM_alpha_rare)

#check the new dataset: it contains both the description of the samples and alpha 
JH02_JH03_RHM_alpha_rare_info

#Test the soil effect on alphadiversity indices

#Soil effect
#Observed
wilcox.test(Observed ~ Soil, data = JH02_JH03_RHM_alpha_rare_info)
#Shannon
wilcox.test(Shannon ~ Soil, data = JH02_JH03_RHM_alpha_rare_info)

#Genotype not corrected for soil (use the parameter 'Description' in the dataset)
#Observed
kruskal.test(Observed ~ Description, data = JH02_JH03_RHM_alpha_rare_info)
posthoc.kruskal.dunn.test (x=JH02_JH03_RHM_alpha_rare_info$Observed, g=JH02_JH03_RHM_alpha_rare_info$Description, p.adjust.method="BH")

#Shannon
kruskal.test(Shannon ~ Description, data = JH02_JH03_RHM_alpha_rare_info)
posthoc.kruskal.dunn.test (x=JH02_JH03_RHM_alpha_rare_info$Shannon, g=JH02_JH03_RHM_alpha_rare_info$Description, p.adjust.method="BH")

#############################################################
#Figure 4: Betadiversity calculations
#Data required: design_2; JH02_JH03_RHM_data_phyloseq_2
#############################################################

#abundance filtering
#Remove OTUs not seen more than 5 times in at least 20% of the samples
JH02_JH03_RHM_data_phyloseq_3 = filter_taxa(JH02_JH03_RHM_data_phyloseq_2, function(x) sum(x > 5) > (0.2*length(x)), TRUE)
JH02_JH03_RHM_data_phyloseq_3
sum(sample_sums(JH02_JH03_RHM_data_phyloseq_2))
sum(sample_sums(JH02_JH03_RHM_data_phyloseq_3))

#ratio filtered reads/total reads
ratio <- sum(sample_sums(JH02_JH03_RHM_data_phyloseq_3))/sum(sample_sums(JH02_JH03_RHM_data_phyloseq_2))*100

#Transform the count in relative abundance cpm
JH02_JH03_RHM_data_phyloseq_prop <- transform_sample_counts(JH02_JH03_RHM_data_phyloseq_3,  function(x) 1e+06 * x/sum(x))

#PCoA weighted unifrac distance
JH02_JH03_RHM_data_phyloseq_prop_wunifrac <- ordinate(JH02_JH03_RHM_data_phyloseq_prop, "PCoA", "unifrac", weighted = TRUE)
plot_ordination(JH02_JH03_RHM_data_phyloseq_prop, JH02_JH03_RHM_data_phyloseq_prop_wunifrac , color = "Description")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH02_JH03_RHM_data_phyloseq_prop, JH02_JH03_RHM_data_phyloseq_prop_wunifrac , shape ="Soil", color = "Genotype")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("brown","orange", "cyan", "darkblue", "magenta"))
p + ggtitle("PCoA 16S data, weighted Unifrac distance")

#PCoA bray distance
JH02_JH03_RHM_data_phyloseq_prop_bray <- ordinate(JH02_JH03_RHM_data_phyloseq_prop, "PCoA", "bray")
plot_ordination(JH02_JH03_RHM_data_phyloseq, JH02_JH03_RHM_data_phyloseq_prop_bray , color = "Description")

#assign shapes to Soil type and color to Sample type
p=plot_ordination(JH02_JH03_RHM_data_phyloseq_prop, JH02_JH03_RHM_data_phyloseq_prop_bray , shape ="Soil", color = "Genotype")
p = p + geom_point(size = 4, alpha = 0.75)
p = p + scale_colour_manual(values = c("brown","orange", "cyan", "darkblue", "magenta"))
p + ggtitle("PCoA 16S data, Bray distance")

#adonis calculations

#rhizosphere effect

#WU distance
WU <- phyloseq::distance(JH02_JH03_RHM_data_phyloseq_prop, "unifrac", weighted= TRUE)
adonis(WU ~ Soil * Microhabitat, data= design_2, permutations = 5000)

#BC distance
BC <- phyloseq::distance(JH02_JH03_RHM_data_phyloseq_prop, "bray")
adonis(BC ~ Soil * Microhabitat, data= design_2, permutations = 5000)

#genotype effect (rhizosphere only)
#Subsetting
JH02_JH03_RHM_data_phyloseq_prop_rhizo <- subset_samples(JH02_JH03_RHM_data_phyloseq_prop, Microhabitat == "Rhizosphere")
design_rhizosphere <- design_2[colnames(otu_table(JH02_JH03_RHM_data_phyloseq_prop_rhizo)), ]

#WU distance
WU <- phyloseq::distance(JH02_JH03_RHM_data_phyloseq_prop_rhizo, "unifrac", weighted= TRUE)
adonis(WU ~ Soil * Genotype, data= design_rhizosphere, permutations = 5000)

#BC distance
BC <- phyloseq::distance(JH02_JH03_RHM_data_phyloseq_prop_rhizo, "bray")
adonis(BC ~ Soil * Genotype, data= design_rhizosphere, permutations = 5000)

#############################################################
#Figure 5: OTUs differentially recruited among samples
#Data required: design_2; JH02_JH03_RHM_data_phyloseq_2; JH02_JH03_RHM_taxa_ordered
#############################################################

#data analysis correct for Soil type
JH02_JH03_RHM_data_phyloseq_Quarry <- subset_samples(JH02_JH03_RHM_data_phyloseq_2, Soil=="Quarryfield")
JH02_JH03_RHM_data_phyloseq_Tayport <- subset_samples(JH02_JH03_RHM_data_phyloseq_2, Soil=="Tayport")

#independent abundance filtering
JH02_JH03_RHM_data_phyloseq_Quarry = filter_taxa(JH02_JH03_RHM_data_phyloseq_Quarry, function(x) sum(x > 5) > (0.2*length(x)), TRUE)
JH02_JH03_RHM_data_phyloseq_Tayport = filter_taxa(JH02_JH03_RHM_data_phyloseq_Tayport, function(x) sum(x > 5) > (0.2*length(x)), TRUE)

############
#quarryfield
############

#extract count data 
JH02_JH03_RHM_OTU_counts_integer <- otu_table(JH02_JH03_RHM_data_phyloseq_Quarry)
countData = as.data.frame(JH02_JH03_RHM_OTU_counts_integer)
colnames(JH02_JH03_RHM_OTU_counts_integer)

#the design file containing sample information
colData = design_2[colnames(JH02_JH03_RHM_OTU_counts_integer), ]
rownames(colData)

#construct a DESeq dataset combining count data and sample information
JH02_JH03_RHM_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Description)

#execute the differential count analysis with the function DESeq 
JH02_JH03_RHM_cds_test <- DESeq(JH02_JH03_RHM_cds, fitType="local", betaPrior = FALSE) 

#define the OTUs significantly enriched in the rhizosphere samples
Soil_Rhizosphere_Karat <- results(JH02_JH03_RHM_cds_test, contrast = c("Description",  "Bulk.quarry", "Karat.quarry")) 
Soil_Rhizosphere_Dema <- results(JH02_JH03_RHM_cds_test , contrast = c("Description",  "Bulk.quarry", "Dema.quarry")) 
Soil_Rhizosphere_RHL1 <- results(JH02_JH03_RHM_cds_test , contrast = c("Description",  "Bulk.quarry", "RHL1.quarry"))
Soil_Rhizosphere_RHP1 <- results(JH02_JH03_RHM_cds_test , contrast = c("Description",  "Bulk.quarry", "RHP1.quarry"))

#inspect a result file
Soil_Rhizosphere_Karat  
mcols(Soil_Rhizosphere_Karat  , use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.01. 
Soil_Rhizosphere_Karat_FDR_001 <- Soil_Rhizosphere_Karat[(rownames(Soil_Rhizosphere_Karat)[which(Soil_Rhizosphere_Karat$padj <0.01)]), ]
Soil_Rhizosphere_Dema_FDR_001 <- Soil_Rhizosphere_Dema[(rownames(Soil_Rhizosphere_Dema)[which(Soil_Rhizosphere_Dema$padj <0.01)]), ]
Soil_Rhizosphere_RHL1_FDR_001 <- Soil_Rhizosphere_RHL1[(rownames(Soil_Rhizosphere_RHL1)[which(Soil_Rhizosphere_RHL1$padj <0.01)]), ]
Soil_Rhizosphere_RHP1_FDR_001 <- Soil_Rhizosphere_RHP1[(rownames(Soil_Rhizosphere_RHP1)[which(Soil_Rhizosphere_RHP1$padj <0.01)]), ]

#Identify OTUs enriched in the rhizosphere (second term of the comparison, negative fold change)
Soil_Rhizosphere_Karat_enriched <-  Soil_Rhizosphere_Karat[(rownames(Soil_Rhizosphere_Karat)[which(Soil_Rhizosphere_Karat$log2FoldChange < 0)]), ]
Soil_Rhizosphere_Dema_enriched <-  Soil_Rhizosphere_Dema[(rownames(Soil_Rhizosphere_Dema)[which(Soil_Rhizosphere_Dema$log2FoldChange < 0)]), ]
Soil_Rhizosphere_RHL1_enriched <-  Soil_Rhizosphere_RHL1[(rownames(Soil_Rhizosphere_RHL1)[which(Soil_Rhizosphere_RHL1$log2FoldChange < 0)]), ]
Soil_Rhizosphere_RHP1_enriched <-  Soil_Rhizosphere_RHP1[(rownames(Soil_Rhizosphere_RHP1)[which(Soil_Rhizosphere_RHP1$log2FoldChange < 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Soil_Rhizosphere_Karat_enriched_FDR001 <- intersect(rownames(Soil_Rhizosphere_Karat_FDR_001), rownames(Soil_Rhizosphere_Karat_enriched))
Soil_Rhizosphere_Dema_enriched_FDR001 <- intersect(rownames(Soil_Rhizosphere_Dema_FDR_001), rownames(Soil_Rhizosphere_Dema_enriched))
Soil_Rhizosphere_RHL1_enriched_FDR001 <- intersect(rownames(Soil_Rhizosphere_RHL1_FDR_001), rownames(Soil_Rhizosphere_RHL1_enriched))
Soil_Rhizosphere_RHP1_enriched_FDR001 <- intersect(rownames(Soil_Rhizosphere_RHP1_FDR_001), rownames(Soil_Rhizosphere_RHP1_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(Soil_Rhizosphere_Karat_enriched_FDR001)
length(Soil_Rhizosphere_Dema_enriched_FDR001)
length(Soil_Rhizosphere_RHL1_enriched_FDR001)
length(Soil_Rhizosphere_RHP1_enriched_FDR001)

#rename the dataset for plotting purposes
Soil_Rhizosphere_Karat_enriched_FDR001_Quarry <- Soil_Rhizosphere_Karat_enriched_FDR001
Soil_Rhizosphere_Dema_enriched_FDR001_Quarry <- Soil_Rhizosphere_Dema_enriched_FDR001
Soil_Rhizosphere_RHL1_enriched_FDR001_Quarry <- Soil_Rhizosphere_RHL1_enriched_FDR001
Soil_Rhizosphere_RHP1_enriched_FDR001_Quarry <- Soil_Rhizosphere_RHP1_enriched_FDR001

##################
#Karat vs. RHL1
Karat_RHL1 <- results(JH02_JH03_RHM_cds_test, contrast = c("Description",  "Karat.quarry", "RHL1.quarry"))

#inspect a result file
Karat_RHL1 
mcols(Karat_RHL1, use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.01. 
Karat_RHL1_FDR_001 <- Karat_RHL1[(rownames(Karat_RHL1)[which(Karat_RHL1$padj <0.01)]), ]

#Identify OTUs enriched in Karat (first term of the comparison)
Karat_enriched_vs_RHL1 <-  Karat_RHL1[(rownames(Karat_RHL1)[which(Karat_RHL1$log2FoldChange > 0)]), ]

#Identify OTUs enriched in RHL1 (second term of the comparison)
RHL1_enriched_vs_Karat <-  Karat_RHL1[(rownames(Karat_RHL1)[which(Karat_RHL1$log2FoldChange < 0)]), ]

#intersect the datasets
Karat_enriched_vs_RHL1_FDR_001 <- intersect(rownames(Karat_RHL1_FDR_001), rownames(Karat_enriched_vs_RHL1))
RHL1_enriched_vs_Karat_FDR_001 <- intersect(rownames(Karat_RHL1_FDR_001), rownames(RHL1_enriched_vs_Karat))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
Karat_enriched_vs_RHL1_FDR_001_2_Quarry <- intersect(Karat_enriched_vs_RHL1_FDR_001, Soil_Rhizosphere_Karat_enriched_FDR001_Quarry)
RHL1_enriched_vs_Karat_FDR_001_2_Quarry <- intersect(RHL1_enriched_vs_Karat_FDR_001, Soil_Rhizosphere_RHL1_enriched_FDR001_Quarry)

#inspect the files
length(Karat_enriched_vs_RHL1_FDR_001_2_Quarry)
length(RHL1_enriched_vs_Karat_FDR_001_2_Quarry)


##################
#Dema vs. RHP1
Dema_RHP1 <- results(JH02_JH03_RHM_cds_test, contrast = c("Description",  "Dema.quarry", "RHP1.quarry"))

#inspect a result file
Dema_RHP1 
mcols(Dema_RHP1 , use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.01 
Dema_RHP1_FDR_001 <- Dema_RHP1[(rownames(Dema_RHP1)[which(Dema_RHP1$padj <0.01)]), ]

#Identify OTUs enriched in Dema (first term of the comparison)
Dema_enriched_vs_RHP1 <-  Dema_RHP1[(rownames(Dema_RHP1)[which(Dema_RHP1$log2FoldChange > 0)]), ]

#Identify OTUs enriched in RHP1 (second term of the comparison)
RHP1_enriched_vs_Dema <-  Dema_RHP1[(rownames(Dema_RHP1)[which(Dema_RHP1$log2FoldChange < 0)]), ]

#intersect the datasets
Dema_enriched_vs_RHP1_FDR_001 <- intersect(rownames(Dema_RHP1_FDR_001), rownames(Dema_enriched_vs_RHP1))
RHP1_enriched_vs_Dema_FDR_001 <- intersect(rownames(Dema_RHP1_FDR_001), rownames(RHP1_enriched_vs_Dema))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
Dema_enriched_vs_RHP1_FDR_001_2_Quarry <- intersect(Dema_enriched_vs_RHP1_FDR_001, Soil_Rhizosphere_Dema_enriched_FDR001_Quarry)
RHP1_enriched_vs_Dema_FDR_001_2_Quarry <- intersect(RHP1_enriched_vs_Dema_FDR_001, Soil_Rhizosphere_RHP1_enriched_FDR001_Quarry)

#inspect the files
length(Dema_enriched_vs_RHP1_FDR_001_2_Quarry)
length(RHP1_enriched_vs_Dema_FDR_001_2_Quarry)

##################
#Karat vs. Dema
Karat_Dema <- results(JH02_JH03_RHM_cds_test, contrast = c("Description",  "Karat.quarry", "Dema.quarry"))

#inspect a result file
Karat_Dema 
mcols(Karat_Dema , use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.01 
Karat_Dema_FDR_001 <- Karat_Dema[(rownames(Karat_Dema)[which(Karat_Dema$padj <0.01)]), ]

#Identify OTUs enriched in Karat (first term of the comparison)
Karat_enriched_vs_Dema <-  Karat_Dema[(rownames(Karat_Dema)[which(Karat_Dema$log2FoldChange > 0)]), ]

#Identify OTUs enriched in Dema (second term of the comparison)
Dema_enriched_vs_Karat <-  Karat_Dema[(rownames(Karat_Dema)[which(Karat_Dema$log2FoldChange < 0)]), ]

#intersect the datasets
Karat_enriched_vs_Dema_FDR_001 <- intersect(rownames(Karat_Dema_FDR_001), rownames(Karat_enriched_vs_Dema))
Dema_enriched_vs_Karat_FDR_001 <- intersect(rownames(Karat_Dema_FDR_001), rownames(Dema_enriched_vs_Karat))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
Karat_enriched_vs_Dema_FDR_001_2_Quarry <- intersect(Karat_enriched_vs_Dema_FDR_001, Soil_Rhizosphere_Karat_enriched_FDR001_Quarry)
Dema_enriched_vs_Karat_FDR_001_2_Quarry <- intersect(Dema_enriched_vs_Karat_FDR_001, Soil_Rhizosphere_Dema_enriched_FDR001_Quarry)

#inspect the files
length(Karat_enriched_vs_Dema_FDR_001_2_Quarry)
length(Dema_enriched_vs_Karat_FDR_001_2_Quarry)

#Pair-wise wt-mutants comparisons data visualisation
#########################################

#create new datasets combining statisical and tanonomy information
#Karat
Karat_enriched_RHL1_taxa_Quarry <- JH02_JH03_RHM_taxa_ordered[Karat_enriched_vs_RHL1_FDR_001_2_Quarry, ]
Karat_enriched_vs_RHL1_FDR_001_2_Quarry_taxa <- cbind(as.data.frame(Karat_RHL1[Karat_enriched_vs_RHL1_FDR_001_2_Quarry, ]), Karat_enriched_RHL1_taxa_Quarry)
#RHL1
RHL1_enriched_Karat_taxa_Quarry <- JH02_JH03_RHM_taxa_ordered[RHL1_enriched_vs_Karat_FDR_001_2_Quarry, ]
RHL1_enriched_vs_Karat_FDR_001_2_Quarry_taxa <- cbind(as.data.frame(Karat_RHL1[RHL1_enriched_vs_Karat_FDR_001_2_Quarry, ]), RHL1_enriched_Karat_taxa_Quarry)
#Dema
Dema_enriched_RHP1_taxa_Quarry <- JH02_JH03_RHM_taxa_ordered[Dema_enriched_vs_RHP1_FDR_001_2_Quarry, ]
Dema_enriched_vs_RHP1_FDR_001_2_Quarry_taxa <- cbind(as.data.frame(Dema_RHP1[Dema_enriched_vs_RHP1_FDR_001_2_Quarry, ]), Dema_enriched_RHP1_taxa_Quarry)
#RHP1
RHP1_enriched_Dema_taxa_Quarry <- JH02_JH03_RHM_taxa_ordered[RHP1_enriched_vs_Dema_FDR_001_2_Quarry, ]
RHP1_enriched_vs_Dema_FDR_001_2_Quarry_taxa <- cbind(as.data.frame(Dema_RHP1[RHP1_enriched_vs_Dema_FDR_001_2_Quarry, ]), RHP1_enriched_Dema_taxa_Quarry)

#save the files and combine them into Supplementary Database 1
#worksheet ws 7
#write.table(Karat_enriched_vs_RHL1_FDR_001_2_Quarry_taxa, file="Karat_enriched_RHL1_taxa_Quarry_FDR001.txt ", sep="\t")
#worksheet ws 8
#write.table(RHL1_enriched_vs_Karat_FDR_001_2_Quarry_taxa, file="RHL1_enriched_Karat_taxa_Quarry_FDR001.txt ", sep="\t")
#worksheet ws 9
#write.table(Dema_enriched_vs_RHP1_FDR_001_2_Quarry_taxa, file="Dema_enriched_RHP1_taxa_Quarry_FDR001.txt ", sep="\t")
#worksheet ws 10
#write.table(RHP1_enriched_vs_Dema_FDR_001_2_Quarry_taxa, file="RHP1_enriched_Dema_taxa_Quarry_FDR001.txt ", sep="\t")

#identify unique OTUs and their taxonomic assignments
Karat_RHL1_Quarry <- union(Karat_enriched_vs_RHL1_FDR_001_2_Quarry, RHL1_enriched_vs_Karat_FDR_001_2_Quarry)
Dema_RHP1_Quarry <- union(Dema_enriched_vs_RHP1_FDR_001_2_Quarry, RHP1_enriched_vs_Dema_FDR_001_2_Quarry)
Quarry_differentially_regulated <- unique(union(Karat_RHL1_Quarry, Dema_RHP1_Quarry))
Quarry_differentially_regulated_taxa <- JH02_JH03_RHM_taxa_ordered[Quarry_differentially_regulated , ]
Abundant_Orders_Quarry <- count(Quarry_differentially_regulated_taxa, "Order")
Abundant_Orders_Quarry

#data visualisation: pair-wise comparison wild type/mutant Quarryfield

#color coding
#Generate a dataset that comprosises the taxonomy information of the OTUs enriched in both terms of the comparisons
#Karat_RHL1
Karat_RHL1_taxa <- rbind(Karat_enriched_vs_RHL1_FDR_001_2_Quarry_taxa, RHL1_enriched_vs_Karat_FDR_001_2_Quarry_taxa)

#create distinct datasets comprising ony the major taxonomic groups at order level with at least 10 representative OTUs
Karat_RHL1_taxa_Burkholderiales <- rownames(Karat_RHL1_taxa)[which(Karat_RHL1_taxa$Order == " o__Burkholderiales")]
Karat_RHL1_taxa_Rhizobiales <- rownames(Karat_RHL1_taxa)[which(Karat_RHL1_taxa$Order == " o__Rhizobiales")]                                                   
Karat_RHL1_taxa_Xanthomonadales <- rownames(Karat_RHL1_taxa)[which(Karat_RHL1_taxa$Order == " o__Xanthomonadales")]
Karat_RHL1_taxa_Sphingomonadales <- rownames(Karat_RHL1_taxa)[which(Karat_RHL1_taxa$Order == " o__Sphingomonadales")]
Karat_RHL1_taxa_Actinomycetales <- rownames(Karat_RHL1_taxa)[which(Karat_RHL1_taxa$Order == " o__Actinomycetales")]
Karat_RHL1_taxa_other_taxa <- setdiff(rownames(Karat_RHL1_taxa), c(Karat_RHL1_taxa_Burkholderiales,Karat_RHL1_taxa_Rhizobiales,Karat_RHL1_taxa_Xanthomonadales,Karat_RHL1_taxa_Sphingomonadales,Karat_RHL1_taxa_Actinomycetales))

fig5a_Karat_RHL1_colors <- ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_Burkholderiales , "darkorange",
                                  ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_Rhizobiales , "deepskyblue",
                                         ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_Xanthomonadales , "forestgreen",
                                                ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_Sphingomonadales , "magenta",
                                                       ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_Actinomycetales , "red",                                           
                                                              ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_other_taxa, "blueviolet", "grey"))))))                                    


#Dema_RHP1
Dema_RHP1_taxa <- rbind(Dema_enriched_vs_RHP1_FDR_001_2_Quarry_taxa, RHP1_enriched_vs_Dema_FDR_001_2_Quarry_taxa)

Dema_RHP1_taxa_Burkholderiales <- rownames(Dema_RHP1_taxa)[which(Dema_RHP1_taxa$Order == " o__Burkholderiales")]
Dema_RHP1_taxa_Rhizobiales <- rownames(Dema_RHP1_taxa)[which(Dema_RHP1_taxa$Order == " o__Rhizobiales")]                                                   
Dema_RHP1_taxa_Xanthomonadales <- rownames(Dema_RHP1_taxa)[which(Dema_RHP1_taxa$Order == " o__Xanthomonadales")]
Dema_RHP1_taxa_Sphingomonadales <- rownames(Dema_RHP1_taxa)[which(Dema_RHP1_taxa$Order == " o__Sphingomonadales")]
Dema_RHP1_taxa_Actinomycetales <- rownames(Dema_RHP1_taxa)[which(Dema_RHP1_taxa$Order == " o__Actinomycetales")]
Dema_RHP1_taxa_other_taxa <- setdiff(rownames(Dema_RHP1_taxa), c(Dema_RHP1_taxa_Burkholderiales,Dema_RHP1_taxa_Rhizobiales,Dema_RHP1_taxa_Xanthomonadales,Dema_RHP1_taxa_Sphingomonadales,Dema_RHP1_taxa_Actinomycetales))

fig5a_Dema_RHP1_colors <- ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_Burkholderiales , "darkorange",
                                 ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_Rhizobiales , "deepskyblue",
                                        ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_Xanthomonadales , "forestgreen",
                                               ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_Sphingomonadales , "magenta",
                                                      ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_Actinomycetales , "red",                                           
                                                             ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_other_taxa, "blueviolet", "grey"))))))

#plotting                                      
dev.off()
par(mfrow=c(2,1))            
plotMA(Karat_RHL1, alpha = 0.01, colSig = fig5a_Karat_RHL1_colors, colNonSig = "grey", ylim = c(-6, 6),  main="Karat vs. RHL1, Quarryfield")
plotMA(Dema_RHP1, alpha = 0.01, colSig = fig5a_Dema_RHP1_colors, colNonSig = "grey", ylim = c(-6, 6),   main="Dema vs. RHP1, Quarryfield")

#save the plots as .eps and generate Figure 5a

#################
#Tayport
################

#extract count data 
JH02_JH03_RHM_OTU_counts_integer <- otu_table(JH02_JH03_RHM_data_phyloseq_Tayport)
countData = as.data.frame(JH02_JH03_RHM_OTU_counts_integer)
colnames(JH02_JH03_RHM_OTU_counts_integer)

#the design file containing sample information
colData = design_2[colnames(JH02_JH03_RHM_OTU_counts_integer), ]

#construct a DESeq dataset combining count data and sample information
JH02_JH03_RHM_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Description)

#execute the differential count analysis with the function DESeq 
JH02_JH03_RHM_cds_test <- DESeq(JH02_JH03_RHM_cds, fitType="local", betaPrior = FALSE) 

#define the OTUs significantly enriched in the rhizosphere samples
Soil_Rhizosphere_Karat <- results(JH02_JH03_RHM_cds_test, contrast = c("Description",  "Bulk.tayport", "Karat.tayport")) 
Soil_Rhizosphere_Dema <- results(JH02_JH03_RHM_cds_test , contrast = c("Description",  "Bulk.tayport", "Dema.tayport")) 
Soil_Rhizosphere_RHL1 <- results(JH02_JH03_RHM_cds_test , contrast = c("Description",  "Bulk.tayport", "RHL1.tayport"))
Soil_Rhizosphere_RHP1 <- results(JH02_JH03_RHM_cds_test , contrast = c("Description",  "Bulk.tayport", "RHP1.tayport"))

#inspect a result file
Soil_Rhizosphere_Karat  
mcols(Soil_Rhizosphere_Karat  , use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.01 
Soil_Rhizosphere_Karat_FDR_001 <- Soil_Rhizosphere_Karat[(rownames(Soil_Rhizosphere_Karat)[which(Soil_Rhizosphere_Karat$padj <0.01)]), ]
Soil_Rhizosphere_Dema_FDR_001 <- Soil_Rhizosphere_Dema[(rownames(Soil_Rhizosphere_Dema)[which(Soil_Rhizosphere_Dema$padj <0.01)]), ]
Soil_Rhizosphere_RHL1_FDR_001 <- Soil_Rhizosphere_RHL1[(rownames(Soil_Rhizosphere_RHL1)[which(Soil_Rhizosphere_RHL1$padj <0.01)]), ]
Soil_Rhizosphere_RHP1_FDR_001 <- Soil_Rhizosphere_RHP1[(rownames(Soil_Rhizosphere_RHP1)[which(Soil_Rhizosphere_RHP1$padj <0.01)]), ]

#Identify OTUs enriched in the rhizosphere (second term of the comparison, negative fold change)
Soil_Rhizosphere_Karat_enriched <-  Soil_Rhizosphere_Karat[(rownames(Soil_Rhizosphere_Karat)[which(Soil_Rhizosphere_Karat$log2FoldChange < 0)]), ]
Soil_Rhizosphere_Dema_enriched <-  Soil_Rhizosphere_Dema[(rownames(Soil_Rhizosphere_Dema)[which(Soil_Rhizosphere_Dema$log2FoldChange < 0)]), ]
Soil_Rhizosphere_RHL1_enriched <-  Soil_Rhizosphere_RHL1[(rownames(Soil_Rhizosphere_RHL1)[which(Soil_Rhizosphere_RHL1$log2FoldChange < 0)]), ]
Soil_Rhizosphere_RHP1_enriched <-  Soil_Rhizosphere_RHP1[(rownames(Soil_Rhizosphere_RHP1)[which(Soil_Rhizosphere_RHP1$log2FoldChange < 0)]), ]

#To identify OTUs significantly enriched in the compartment of interest an intersection of these lists is required
Soil_Rhizosphere_Karat_enriched_FDR001 <- intersect(rownames(Soil_Rhizosphere_Karat_FDR_001), rownames(Soil_Rhizosphere_Karat_enriched))
Soil_Rhizosphere_Dema_enriched_FDR001 <- intersect(rownames(Soil_Rhizosphere_Dema_FDR_001), rownames(Soil_Rhizosphere_Dema_enriched))
Soil_Rhizosphere_RHL1_enriched_FDR001 <- intersect(rownames(Soil_Rhizosphere_RHL1_FDR_001), rownames(Soil_Rhizosphere_RHL1_enriched))
Soil_Rhizosphere_RHP1_enriched_FDR001 <- intersect(rownames(Soil_Rhizosphere_RHP1_FDR_001), rownames(Soil_Rhizosphere_RHP1_enriched))

#Define the number of OTUs significantly enriched in and differentiating between each plant-associated compartment and unplanted soil
length(Soil_Rhizosphere_Karat_enriched_FDR001)
length(Soil_Rhizosphere_Dema_enriched_FDR001)
length(Soil_Rhizosphere_RHL1_enriched_FDR001)
length(Soil_Rhizosphere_RHP1_enriched_FDR001)

#rename the dataset for plotting purposes
Soil_Rhizosphere_Karat_enriched_FDR001_Tayport <- Soil_Rhizosphere_Karat_enriched_FDR001
Soil_Rhizosphere_Dema_enriched_FDR001_Tayport <- Soil_Rhizosphere_Dema_enriched_FDR001
Soil_Rhizosphere_RHL1_enriched_FDR001_Tayport <- Soil_Rhizosphere_RHL1_enriched_FDR001
Soil_Rhizosphere_RHP1_enriched_FDR001_Tayport <- Soil_Rhizosphere_RHP1_enriched_FDR001

##################
#Karat vs. RHL1
Karat_RHL1 <- results(JH02_JH03_RHM_cds_test, contrast = c("Description",  "Karat.tayport", "RHL1.tayport"))

#inspect a result file
Karat_RHL1   
mcols(Karat_RHL1, use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.01. 
Karat_RHL1_FDR_001 <- Karat_RHL1[(rownames(Karat_RHL1)[which(Karat_RHL1$padj <0.01)]), ]

#Identify OTUs enriched in Karat (first term of the comparison)
Karat_enriched_vs_RHL1 <-  Karat_RHL1[(rownames(Karat_RHL1)[which(Karat_RHL1$log2FoldChange > 0)]), ]

#Identify OTUs enriched in RHL1 (second term of the comparison)
RHL1_enriched_vs_Karat <-  Karat_RHL1[(rownames(Karat_RHL1)[which(Karat_RHL1$log2FoldChange < 0)]), ]

#intersect the datasets
Karat_enriched_vs_RHL1_FDR_001 <- intersect(rownames(Karat_RHL1_FDR_001), rownames(Karat_enriched_vs_RHL1))
RHL1_enriched_vs_Karat_FDR_001 <- intersect(rownames(Karat_RHL1_FDR_001), rownames(RHL1_enriched_vs_Karat))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
Karat_enriched_vs_RHL1_FDR_001_2_tayport <- intersect(Karat_enriched_vs_RHL1_FDR_001, Soil_Rhizosphere_Karat_enriched_FDR001_Tayport)
RHL1_enriched_vs_Karat_FDR_001_2_tayport <- intersect(RHL1_enriched_vs_Karat_FDR_001, Soil_Rhizosphere_RHL1_enriched_FDR001_Tayport)

#inspect the generated files
length(Karat_enriched_vs_RHL1_FDR_001_2_tayport)
length(RHL1_enriched_vs_Karat_FDR_001_2_tayport)

##################
#Dema vs. RHP1
Dema_RHP1 <- results(JH02_JH03_RHM_cds_test, contrast = c("Description",  "Dema.tayport", "RHP1.tayport"))

#inspect the result file
Dema_RHP1 
mcols(Dema_RHP1, use.names=TRUE)

# extract  OTUs whose adjusted p.value in a given comparison is below 0.01. 
Dema_RHP1_FDR_001 <- Dema_RHP1[(rownames(Dema_RHP1)[which(Dema_RHP1$padj <0.01)]), ]

#Identify OTUs enriched in Dema (first term of the comparison)
Dema_enriched_vs_RHP1 <-  Dema_RHP1[(rownames(Dema_RHP1)[which(Dema_RHP1$log2FoldChange > 0)]), ]

#Identify OTUs enriched in RHP1 (second term of the comparison)
RHP1_enriched_vs_Dema <-  Dema_RHP1[(rownames(Dema_RHP1)[which(Dema_RHP1$log2FoldChange < 0)]), ]

#intersect the datasets
Dema_enriched_vs_RHP1_FDR_001 <- intersect(rownames(Dema_RHP1_FDR_001), rownames(Dema_enriched_vs_RHP1))
RHP1_enriched_vs_Dema_FDR_001 <- intersect(rownames(Dema_RHP1_FDR_001), rownames(RHP1_enriched_vs_Dema))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
Dema_enriched_vs_RHP1_FDR_001_2_tayport <- intersect(Dema_enriched_vs_RHP1_FDR_001, Soil_Rhizosphere_Dema_enriched_FDR001_Tayport)
RHP1_enriched_vs_Dema_FDR_001_2_tayport <- intersect(RHP1_enriched_vs_Dema_FDR_001, Soil_Rhizosphere_RHP1_enriched_FDR001_Tayport)

#inspect the files
length(Dema_enriched_vs_RHP1_FDR_001_2_tayport)
length(RHP1_enriched_vs_Dema_FDR_001_2_tayport)

##################
#Karat vs. Dema
Karat_Dema <- results(JH02_JH03_RHM_cds_test, contrast = c("Description",  "Karat.tayport", "Dema.tayport"))

#inspect a result file
Karat_Dema 
mcols(Karat_Dema, use.names=TRUE)

#extract  OTUs whose adjusted p.value in a given comparison is below 0.01 
Karat_Dema_FDR_001 <- Karat_Dema[(rownames(Karat_Dema)[which(Karat_Dema$padj <0.01)]), ]

#Identify OTUs enriched in Karat (first term of the comparison)
Karat_enriched_vs_Dema <-  Karat_Dema[(rownames(Karat_Dema)[which(Karat_Dema$log2FoldChange > 0)]), ]

#Identify OTUs enriched in Dema (second term of the comparison)
Dema_enriched_vs_Karat <-  Karat_Dema[(rownames(Karat_Dema)[which(Karat_Dema$log2FoldChange < 0)]), ]

#intersect the datasets
Karat_enriched_vs_Dema_FDR_001 <- intersect(rownames(Karat_Dema_FDR_001), rownames(Karat_enriched_vs_Dema))
Dema_enriched_vs_Karat_FDR_001 <- intersect(rownames(Karat_Dema_FDR_001), rownames(Dema_enriched_vs_Karat))

#further filtering for OTUs significantly enriched vs. soil in the respective comparision
Karat_enriched_vs_Dema_FDR_001_2_tayport <- intersect(Karat_enriched_vs_Dema_FDR_001, Soil_Rhizosphere_Karat_enriched_FDR001_Tayport)
Dema_enriched_vs_Karat_FDR_001_2_tayport <- intersect(Dema_enriched_vs_Karat_FDR_001, Soil_Rhizosphere_Dema_enriched_FDR001_Tayport)

#inspect the files
length(Karat_enriched_vs_Dema_FDR_001_2_tayport)
length(Dema_enriched_vs_Karat_FDR_001_2_tayport)

#Pair-wise wt-mutants comparisons data visualisation
#########################################

#create new datasets combining statisical and tanonomy information
Karat_enriched_RHL1_taxa_tayport <- JH02_JH03_RHM_taxa_ordered[Karat_enriched_vs_RHL1_FDR_001_2_tayport, ]
Karat_enriched_vs_RHL1_FDR_001_2_tayport_taxa <- cbind(as.data.frame(Karat_RHL1[Karat_enriched_vs_RHL1_FDR_001_2_tayport, ]), Karat_enriched_RHL1_taxa_tayport)
RHL1_enriched_Karat_taxa_tayport <- JH02_JH03_RHM_taxa_ordered[RHL1_enriched_vs_Karat_FDR_001_2_tayport, ]
RHL1_enriched_vs_Karat_FDR_001_2_tayport_taxa <- cbind(as.data.frame(Karat_RHL1[RHL1_enriched_vs_Karat_FDR_001_2_tayport, ]), RHL1_enriched_Karat_taxa_tayport)
Dema_enriched_RHP1_taxa_tayport <- JH02_JH03_RHM_taxa_ordered[Dema_enriched_vs_RHP1_FDR_001_2_tayport, ]
Dema_enriched_vs_RHP1_FDR_001_2_tayport_taxa <- cbind(as.data.frame(Dema_RHP1[Dema_enriched_vs_RHP1_FDR_001_2_tayport, ]), Dema_enriched_RHP1_taxa_tayport)
RHP1_enriched_Dema_taxa_tayport <- JH02_JH03_RHM_taxa_ordered[RHP1_enriched_vs_Dema_FDR_001_2_tayport, ]
RHP1_enriched_vs_Dema_FDR_001_2_tayport_taxa <- cbind(as.data.frame(Dema_RHP1[RHP1_enriched_vs_Dema_FDR_001_2_tayport, ]), RHP1_enriched_Dema_taxa_tayport)

Karat_enriched_vs_Dema_taxa_tayport <- JH02_JH03_RHM_taxa_ordered[Karat_enriched_vs_Dema_FDR_001_2_tayport, ]
Karat_enriched_vs_Dema_FDR_001_2_tayport_taxa <- cbind(as.data.frame(Karat_Dema[Karat_enriched_vs_Dema_FDR_001_2_tayport, ]), Karat_enriched_vs_Dema_taxa_tayport)
Dema_enriched_vs_Karat_taxa_tayport <- JH02_JH03_RHM_taxa_ordered[Dema_enriched_vs_Karat_FDR_001_2_tayport, ]
Dema_enriched_vs_Karat_FDR_001_2_tayport_taxa <- cbind(as.data.frame(Karat_Dema[Dema_enriched_vs_Karat_FDR_001_2_tayport, ]), Dema_enriched_vs_Karat_taxa_tayport)

#save the files and combine them into Supplementary Database 1
#worksheet ws 11
#write.table(Karat_enriched_vs_RHL1_FDR_001_2_tayport_taxa, file="Karat_enriched_RHL1_taxa_tayport_FDR001.txt ", sep="\t")
#worksheet ws 12
#write.table(RHL1_enriched_vs_Karat_FDR_001_2_tayport_taxa, file="RHL1_enriched_Karat_taxa_tayport_FDR001.txt ", sep="\t")
#worksheet ws 13
#write.table(Dema_enriched_vs_RHP1_FDR_001_2_tayport_taxa, file="Dema_enriched_RHP1_taxa_tayport_FDR001.txt ", sep="\t")
#worksheet ws 14
#write.table(RHP1_enriched_vs_Dema_FDR_001_2_tayport_taxa, file="RHP1_enriched_Dema_taxa_tayport_FDR001.txt ", sep="\t")
#worksheet ws 15
#write.table(Karat_enriched_vs_Dema_FDR_001_2_tayport_taxa, file="Karat_enriched_Dema_taxa_tayport_FDR001.txt ", sep="\t")
#worksheet ws 16
#write.table(Dema_enriched_vs_Karat_FDR_001_2_tayport_taxa, file="Dema_enriched_Karat_taxa_tayport_FDR001.txt ", sep="\t")

#identify unique OTUs and their taxonomic assignments
Karat_RHL1_tayport <- union(Karat_enriched_vs_RHL1_FDR_001_2_tayport, RHL1_enriched_vs_Karat_FDR_001_2_tayport)
Dema_RHP1_tayport <- union(Dema_enriched_vs_RHP1_FDR_001_2_tayport, RHP1_enriched_vs_Dema_FDR_001_2_tayport)
tayport_differentially_regulated <- unique(union(Karat_RHL1_tayport, Dema_RHP1_tayport))
tayport_differentially_regulated_taxa <- JH02_JH03_RHM_taxa_ordered[tayport_differentially_regulated , ]
Abundant_Orders_tayport <- count(tayport_differentially_regulated_taxa, "Order")
Abundant_Orders_tayport 

#data visualisation: pair-wise comparison wild type/mutant Tayport

#color coding
#Generate a dataset that comprosises the taxonomy information of the OTUs enriched in both terms of the comparisons
#Karat_RHL1
Karat_RHL1_taxa <- rbind(Karat_enriched_vs_RHL1_FDR_001_2_tayport_taxa, RHL1_enriched_vs_Karat_FDR_001_2_tayport_taxa)

#create distinct datasets comprising ony the major taxonomic groups at order level with at least 10 representative OTUs
Karat_RHL1_taxa_Burkholderiales <- rownames(Karat_RHL1_taxa)[which(Karat_RHL1_taxa$Order == " o__Burkholderiales")]
Karat_RHL1_taxa_Rhizobiales <- rownames(Karat_RHL1_taxa)[which(Karat_RHL1_taxa$Order == " o__Rhizobiales")]                                                   
Karat_RHL1_taxa_Xanthomonadales <- rownames(Karat_RHL1_taxa)[which(Karat_RHL1_taxa$Order == " o__Xanthomonadales")]
Karat_RHL1_taxa_Sphingomonadales <- rownames(Karat_RHL1_taxa)[which(Karat_RHL1_taxa$Order == " o__Sphingomonadales")]
Karat_RHL1_taxa_Actinomycetales <- rownames(Karat_RHL1_taxa)[which(Karat_RHL1_taxa$Order == " o__Actinomycetales")]
Karat_RHL1_taxa_other_taxa <- setdiff(rownames(Karat_RHL1_taxa), c(Karat_RHL1_taxa_Burkholderiales,Karat_RHL1_taxa_Rhizobiales,Karat_RHL1_taxa_Xanthomonadales,Karat_RHL1_taxa_Sphingomonadales,Karat_RHL1_taxa_Actinomycetales))
  
fig5b_Karat_RHL1_colors <- ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_Burkholderiales , "darkorange",
                       ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_Rhizobiales , "deepskyblue",
                              ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_Xanthomonadales , "forestgreen",
                                     ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_Sphingomonadales , "magenta",
                                          ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_Actinomycetales , "red",                                           
                                            ifelse(rownames(Karat_RHL1) %in% Karat_RHL1_taxa_other_taxa, "blueviolet", "grey"))))))                                    


#Dema_RHP1
Dema_RHP1_taxa <- rbind(Dema_enriched_vs_RHP1_FDR_001_2_tayport_taxa, RHP1_enriched_vs_Dema_FDR_001_2_tayport_taxa)

Dema_RHP1_taxa_Burkholderiales <- rownames(Dema_RHP1_taxa)[which(Dema_RHP1_taxa$Order == " o__Burkholderiales")]
Dema_RHP1_taxa_Rhizobiales <- rownames(Dema_RHP1_taxa)[which(Dema_RHP1_taxa$Order == " o__Rhizobiales")]                                                   
Dema_RHP1_taxa_Xanthomonadales <- rownames(Dema_RHP1_taxa)[which(Dema_RHP1_taxa$Order == " o__Xanthomonadales")]
Dema_RHP1_taxa_Sphingomonadales <- rownames(Dema_RHP1_taxa)[which(Dema_RHP1_taxa$Order == " o__Sphingomonadales")]
Dema_RHP1_taxa_Actinomycetales <- rownames(Dema_RHP1_taxa)[which(Dema_RHP1_taxa$Order == " o__Actinomycetales")]
Dema_RHP1_taxa_other_taxa <- setdiff(rownames(Dema_RHP1_taxa), c(Dema_RHP1_taxa_Burkholderiales,Dema_RHP1_taxa_Rhizobiales,Dema_RHP1_taxa_Xanthomonadales,Dema_RHP1_taxa_Sphingomonadales,Dema_RHP1_taxa_Actinomycetales))

fig5b_Dema_RHP1_colors <- ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_Burkholderiales , "darkorange",
                                  ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_Rhizobiales , "deepskyblue",
                                         ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_Xanthomonadales , "forestgreen",
                                                ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_Sphingomonadales , "magenta",
                                                       ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_Actinomycetales , "red",                                           
                                                              ifelse(rownames(Dema_RHP1) %in% Dema_RHP1_taxa_other_taxa, "blueviolet", "grey"))))))

#plotting                                      
dev.off()
par(mfrow=c(2,1))            
plotMA(Karat_RHL1, alpha = 0.01, colSig = fig5b_Karat_RHL1_colors, colNonSig = "grey", ylim = c(-6, 6),  main="Karat vs. RHL1, Tayport")
plotMA(Dema_RHP1, alpha = 0.01, colSig = fig5b_Dema_RHP1_colors, colNonSig = "grey", ylim = c(-6, 6),   main="Dema vs. RHP1, Tayport")

#save the plots as .eps and combine them with plots generated at lines 739-743 to generate Figure 5
#############################################################

#############################################################
#Pair-wise comparison Karat vs. Karat Control
#############################################################

#data analysis corrected for Soil type: use the original phyloseq object to include Karat controls
JH02_JH03_RHM_data_phyloseq_Tayport <- subset_samples(JH02_JH03_RHM_data_phyloseq, Soil=="Tayport")

#independent abundance filtering
JH02_JH03_RHM_data_phyloseq_Tayport = filter_taxa(JH02_JH03_RHM_data_phyloseq_Tayport, function(x) sum(x > 5) > (0.2*length(x)), TRUE)

#extract count data 
JH02_JH03_RHM_OTU_counts_integer <- otu_table(JH02_JH03_RHM_data_phyloseq_Tayport)
countData = as.data.frame(JH02_JH03_RHM_OTU_counts_integer)
colnames(JH02_JH03_RHM_OTU_counts_integer)

#the design file containing sample information
colData = design[colnames(JH02_JH03_RHM_OTU_counts_integer), ]

#construct a DESeq dataset combining count data and sample information
JH02_JH03_RHM_cds <- DESeqDataSetFromMatrix(countData =countData, colData=colData , design= ~Description)

#execute the differential count analysis with the function DESeq 
JH02_JH03_RHM_cds_test <- DESeq(JH02_JH03_RHM_cds, fitType="local", betaPrior = FALSE)

#set up the contrasts
Karat_Karat_control <- results(JH02_JH03_RHM_cds_test, contrast = c("Description",  "Karat.control", "Karat.tayport")) 

#inspect the result file
Karat_Karat_control 
mcols(Karat_Karat_control, use.names=TRUE)

#identify OTUs differentially recruited
Karat_Karat_control_FDR_001 <- Karat_Karat_control[(rownames(Karat_Karat_control)[which(Karat_Karat_control$padj <0.01)]), ]

#define how do the differentially recruited OTU relate with OTU enriched in muntant/wild type plants
KC_RHL1 <- intersect(RHL1_enriched_vs_Karat_FDR_001_2_tayport, rownames(Karat_Karat_control_FDR_001))
KC_Karat <- intersect(Karat_enriched_vs_RHL1_FDR_001_2_tayport, rownames(Karat_Karat_control_FDR_001))

#define the proportion of reads assigned to these OTUs in the dataset
KC_RHL1_reads <- (colSums(dat_count_noplants[KC_RHL1, ])/colSums(dat_count_noplants))*100
KC_Karat_reads <- (colSums(dat_count_noplants[KC_Karat, ])/colSums(dat_count_noplants))*100

#############################################################
