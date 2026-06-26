# Load packages
library(tidyverse)
library(vegan)
library(viridis)
library(RColorBrewer)
library(factoextra)
library(treemap)
library(ggpubr)
library(reshape2)
library(qiime2R)
library(phyloseq)
library(fantaxtic)
library(microbiome)
library(biomehorizon)  
library(palmerpenguins)
library(ggbreak) 
library(patchwork)

## Load 18S count and taxonomy files
# Load 18S ASV count table
table <- read_qza(file="table_18S.qza")
count_tab <- table$data %>% as.data.frame() # Convert to data frame 

# Load 18S taxonomy file
taxonomy <- read_qza(file="taxonomy_18S.qza")
tax_tab <- taxonomy$data %>% # Convert to data frame, tab separate and rename taxa levels, and remove row with confidence values
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Kingdom","Supergroup","Division","Class","Order","Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  dplyr::select(-Confidence)

# Load metadata file
sample_info_tab <- read.table(file="metadata_18S.tsv", header=TRUE, row.names=1, check.names=F, fill = TRUE, sep="\t")

# Merge into phyloseq object
ps <- phyloseq(tax_table(as.matrix(tax_tab)), otu_table(count_tab, taxa_are_rows = T), sample_data(sample_info_tab))

# Remove unwanted groups (vertebrates, land plants, and macroalgae)
ps_new = subset_taxa(ps, Class != "Craniata" |is.na(Class)) 
ps_new = subset_taxa(ps_new, Division !="Streptophyta" |is.na(Division))
ps_new = subset_taxa(ps_new, Division !="Rhodophyta" |is.na(Division))

# Add an unassigned label to NA values in the taxonomy table
ps_new = name_na_taxa(ps_new) 

# Remove the controls for now (samples with < 100,000 reads)
ps_sub = prune_samples(sample_sums(ps_new )>=100000, ps_new) 

# Remove other blanks and samples that appear to be outliers (3C, 2B, and 11B - check out the stacked bar plots in qiime2 view)
ps_sub = subset_samples(ps_sub, sample_names(ps_sub) != "GMT1_18S_3C" & sample_names(ps_sub) != "GMT1_18S_2B" & sample_names(ps_sub) != "GMT1_18S_11B" & sample_names(ps_sub) != "GMT1_18S_B1" & sample_names(ps_sub) != "GMT1_18S_B2")

# Remove singletons from the dataset
ps_filt = filter_taxa(ps_sub, function (x) {sum(x) > 1}, prune=TRUE)

# Rarefy to even sampling depth
ps_rare <- rarefy_even_depth(ps_filt, sample.size = min(sample_sums(ps_filt)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# PCoA ordination of sampling intervals
ordu = ordinate(ps_rare, "PCoA", "bray")
nb.cols <- 13 # Set up a color palette
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols) # Set up a color palette
p<-plot_ordination(ps_rare, ordu, color="date_range")+theme_bw() + theme(text = element_text(size=16)) +
  geom_point(aes(fill=date_range),size = 5, shape = 21, colour = "black") + scale_fill_manual(values=mycolors)
p$data$date_range <- factor(p$data$date_range, levels = c("12/7-12/14","12/14-12/21","12/21-01/4","01/4-01/18","01/18-02/01","02/01-02/15","02/15-03/01","03/01-03/15", "03/15-03/29","03/29-04/12","04/12-04/26","04/26-05/10","05/10-05/24")) # Set the order of the sampling range 
p
ggsave(filename = "18S_PCoA.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 8, height = 5, dpi = 150) # Save figure as pdf

# Repeat PCoA without May and March dates
ps_sub2 = subset_samples(ps_rare, sample_names(ps_rare) != "GMT1_18S_13A" & sample_names(ps_rare) != "GMT1_18S_13B" & sample_names(ps_rare) != "GMT1_18S_13C" & sample_names(ps_rare) != "GMT1_18S_10A" & sample_names(ps_rare) != "GMT1_18S_10B"& sample_names(ps_rare) != "GMT1_18S_10C") # Subset out May/March samples     
ordu = ordinate(ps_sub2, "PCoA", "bray")
p<-plot_ordination(ps_sub2, ordu, color="date_range")+theme_bw() +  theme(text = element_text(size=16)) +
  geom_point(aes(fill=date_range),size = 5, shape = 21, colour = "black") + scale_fill_manual(values=mycolors)
p$data$date_range <- factor(p$data$date_range, levels = c("12/7-12/14","12/14-12/21","12/21-01/4","01/4-01/18","01/18-02/01","02/01-02/15","02/15-03/01","03/01-03/15", "03/15-03/29","03/29-04/12","04/12-04/26","04/26-05/10","05/10-05/24"))
p
ggsave(filename = "18S_PCoA_zoomed.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 8, height = 5, dpi = 150)

# Estimate ASV richness and Shannon diversity index for each sampling interval
rich_18S <- estimate_richness(ps_sub, measures=c("Observed", "Shannon"))
date = sample_data(ps_sub)$date_range # Add sampling range as a variable
rich_all <- data.frame(rich_18S,date)# Combine diversity and sampling range
df2 = melt(rich_all) # Melt data

# Plot richness and diversity
p <- ggplot(df2, aes(x=factor(date), y=value, fill=date))
p$data$date <- factor(p$data$date, levels = c("12/7-12/14","12/14-12/21","12/21-01/4","01/4-01/18","01/18-02/01","02/01-02/15","02/15-03/01","03/01-03/15", "03/15-03/29","03/29-04/12","04/12-04/26","04/26-05/10","05/10-05/24"))
p + geom_boxplot() + theme_classic() + 
  theme(text = element_text(size=16)) + ylab("18S diversity values") + theme(legend.position="right")+ scale_fill_manual(values=mycolors) + #scale_fill_manual(values=ap)+
  geom_point(aes(fill=date), size =2, shape = 21, colour = "black", position=position_jitterdodge(0.4))+
  theme(axis.text.x=element_text(angle=45, hjust=1, color="black"))+
  theme(axis.title.x =element_blank()) + facet_wrap(~variable,scales="free_y")

ggsave(filename = "Diversity_18S.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height =5, dpi = 150) 

# Set up for taxonomy tree plots
ps <- tax_glom(ps_rare, "Genus", NArm=FALSE)
ps0 <- transform_sample_counts(ps, function(x) x / sum(x)) # Convert to relative abundance
x1 = psmelt(ps0) # Melt data
focus <- c("Opisthokonta", "Alveolata", "Archaeplastida", "Rhizaria", "Stramenopiles", "Excavata","Unknown Eukaryota (Kingdom)") # Focus on major groups at Supergroup level
x1$Supergroup <- ifelse(x1$Supergroup %in% focus, x1$Supergroup, "Others") # Others category for plotting
x1<- x1%>%  # Relabel the Unassigned group
  mutate(Supergroup = as.character(Supergroup)) %>%
  mutate(Supergroup  = replace(Supergroup , Supergroup  == 'Unknown Eukaryota (Kingdom)', 'Unassigned'))%>%
  mutate(Class = as.character(Class)) %>%
  mutate(Order = replace(Order , Order  == 'Unknown Eukaryota (Kingdom)', 'Unassigned'))%>%
  mutate(Order = as.character(Order)) %>%
  mutate(Order = replace(Order , Order == 'Unknown Eukaryota (Kingdom)', 'Unassigned'))

# Plot 18S tree map
pdf("treemap_foram.pdf", width = 10, height = 3)
treemap(dtf = x1,
        title = "18S all samples", 
        algorithm = "pivotSize", border.lwds = c(3,0.5,0.1),
        border.col = c("black", "black", "black"),
        lowerbound.cex.label=0,
        mapping = c(0,0,0),
        index = c("Family", "Species"),
        vSize = "Abundance",
        vColor = "Family",
        palette = "Dark2",
        type="categorical",
        fontsize.labels=c(14,10),
        fontface.labels=c(2,2),
        overlap.labels=0.5,
        bg.labels=220, 
        position.legend = "right",
        align.labels=list(
          c("left", "top"), 
          c("left", "bottom"),
          c("center", "center")),
        inflate.labels=F,
        force.print.labels = F) 
dev.off()

# 18S taxa abundance over time for top 6 groups - class level
class_18S <- tax_glom(ps_foram, taxrank = "Class")
class_18S <- transform_sample_counts(class_18S, function(x)100* x / sum(x)) # Transform to relative abundance
OTU <- otu_table(class_18S)
TAX <- tax_table(class_18S)[,"Class"]
Average <- as.data.frame(rowMeans(OTU))# Estimate average abundance at the order level
names(Average) <- c("Mean") # Rename to mean
Table <- merge(TAX, Average, by=0, all=TRUE)
Table$Row.names = NULL # View abundances at class level - we are interested in the top 6 most RELATIVELY abundant groups

# Temporal abundance plots for the major 18S groups
abund_18S <- psmelt(class_18S) # Melt the data
abund_18S$Class <- as.character(abund_18S$Class) # Convert to character
abund_18S <- abund_18S%>%  # Rename some groups
  mutate(Class = as.character(Class)) %>%
  mutate(Class = replace(Class, Class == 'Unknown Chlorophyta (Division)', 'Unknown_Chlorophyta'))
abund_18S_new <- subset(abund_18S, grepl("Syndiniales|Carpediemonadea|Unknown_Chlorophyta|Bacillariophyta|Dinophyceae|Filosa-Granofilosea", abund_18S$Class)) # Subset to major groups
abund_18S_new$Class <- factor(abund_18S_new$Class,levels=c("Syndiniales", "Carpediemonadea","Unknown_Chlorophyta", "Bacillariophyta", "Dinophyceae","Filosa-Granofilosea")) # Order groups in the plot
p <- ggplot(abund_18S_new, aes(factor(date_range), Abundance, group=Class))
p$data$date_range <- factor(p$data$date_range, levels = c("12/7-12/14","12/14-12/21","12/21-01/4","01/4-01/18","01/18-02/01","02/01-02/15","02/15-03/01","03/01-03/15", "03/15-03/29","03/29-04/12","04/12-04/26","04/26-05/10","05/10-05/24"))
p + geom_point(aes(fill = Class)) + theme(text = element_text(size=16))+
  geom_smooth(se=TRUE, linetype="solid", size=1, level=0.95, fill = "gray45", color = "black", alpha = 0.4, method = "loess")+ 
  facet_wrap(~ Class, ncol =2, scales = "free_y") + theme_bw() + ylab("Relative Abundance (%)") + 
  theme(legend.position = "none")+ guides(fill=guide_legend(ncol=2))+ theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(filename = "18S_abund.pdf", plot = last_plot(),  device = "pdf", path = NULL, scale = 1, width = 7, height = 8, dpi = 150) # Save and export ggplot figure as .eps file 

#### Repeat same analysis with 16S samples (filtered)
# Load ASV count table
#table_16S <- read_qza(file="table_16S.qza")
table_16S <- read_qza(file.choose())
count_tab_16S <- table_16S$data %>% as.data.frame() # Convert to data frame 

# Load taxonomy
#taxonomy_16S <- read_qza(file="taxonomy_16S.qza")
taxonomy_16S <- read_qza(file.choose())
tax_tab_16S <- taxonomy_16S$data %>% 
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
  column_to_rownames("Feature.ID") %>%
  dplyr::select(-Confidence)
tax_tab_16S$Kingdom <- gsub("^.{0,3}", "", tax_tab_16S$Kingdom) # Clean up the taxonomy names for each level
tax_tab_16S$Phylum <- gsub("^.{0,4}", "", tax_tab_16S$Phylum)
tax_tab_16S$Class <- gsub("^.{0,4}", "", tax_tab_16S$Class)
tax_tab_16S$Order <- gsub("^.{0,4}", "", tax_tab_16S$Order)
tax_tab_16S$Family <- gsub("^.{0,4}", "", tax_tab_16S$Family)
tax_tab_16S$Genus <- gsub("^.{0,4}", "", tax_tab_16S$Genus)
tax_tab_16S$Species<- gsub("^.{0,4}", "", tax_tab_16S$Species)

# Load metadata file
#sample_info_tab_16S <- read.table(file="metadata_16S.txt", header=TRUE, row.names=1, check.names=F, fill = TRUE, sep="\t")
sample_info_tab_16S <- read.table(file.choose(), header=TRUE, row.names=1, check.names=F, fill = TRUE, sep="\t")

# Merge into phyloseq object
ps_16S <- phyloseq(tax_table(as.matrix(tax_tab_16S)), otu_table(count_tab_16S, taxa_are_rows = T), sample_data(sample_info_tab_16S))

# Remove chloroplast, mitochondria, and eukaryote reads 
ps_new = subset_taxa(ps_16S, Order != "Chloroplast" |is.na(Order))
ps_new = subset_taxa(ps_new, Family !="Mitochondria" |is.na(Family))
ps_new = subset_taxa(ps_new, Kingdom !="Eukaryota" |is.na(Kingdom))

# Add an unassigned label to NA values in the taxonomy table
ps_new = name_na_taxa(ps_new)

# Remove the controls for now 
ps_sub = subset_samples(ps_new, sample_names(ps_new) != "GMT1_16S_B1" & sample_names(ps_new) != "GMT1_16S_B2" & sample_names(ps_new) != "GMT1_16S_NTC" & sample_names(ps_new) != "RTSF_NTC_1585")

# Remove singletons
ps_filt = filter_taxa(ps_sub, function (x) {sum(x) > 1}, prune=TRUE)

# Rarefy to even sampling depth
ps_rare <- rarefy_even_depth(ps_filt, sample.size = min(sample_sums(ps_filt)), rngseed = 714, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# PCoA ordination for 16S samples
ordu = ordinate(ps_rare, "PCoA", "bray")
p<-plot_ordination(ps_rare, ordu, color="date_range")+theme_bw() + theme(text = element_text(size=16))+
  geom_point(aes(fill=date_range),size = 5, shape = 21, colour = "black") + scale_fill_manual(values=mycolors)
p$data$date_range <- factor(p$data$date_range, levels = c("12/7-12/14","12/14-12/21","12/21-01/4","01/4-01/18","01/18-02/01","02/01-02/15","02/15-03/01","03/01-03/15", "03/15-03/29","03/29-04/12","04/12-04/26","04/26-05/10","05/10-05/24"))
p
ggsave(filename = "16S_PCoA_by_dates.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 8, height = 5, dpi = 150)
pcup<-plot_ordination(ps_rare, ordu, color="cup")+theme_bw() + theme(text = element_text(size=16))+
  geom_point(aes(fill=cup),size = 5, shape = 21, colour = "black") + scale_fill_manual(values=mycolors)
pcup
ggsave(filename = "16S_PCoA_by_cups.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 8, height = 5, dpi = 150)

# Estimate ASV richness and Shannon diversity index for each sampling interval
rich_16S <- estimate_richness(ps_sub, measures=c("Observed", "Shannon"))
cup = sample_data(ps_sub)$cup
date = sample_data(ps_sub)$date_range
rich_all <- data.frame(rich_16S,date,cup)
df2 = melt(rich_all) # Melt data
p <- ggplot(rich_all, aes(x=factor(cup), y=value, fill=cup))
p
p$data$date<- factor(p$data$date, levels = c("12/7-12/14","12/14-12/21","12/21-01/4","01/4-01/18","01/18-02/01","02/01-02/15","02/15-03/01","03/01-03/15", "03/15-03/29","03/29-04/12","04/12-04/26","04/26-05/10","05/10-05/24"))
p$data$cup<- factor(p$data$cup, levels = c(1,2,3,4,5,6,7,8,9,10,11,12,13))
mycolors <- c("#17385A","#427475", "#0E9594","#1F5041","#C2CFB2","#F9DFBE","#ECCBD3","#827191","#92416C","#52113C","#962C2C","#F2673D","#FF9B71")

#Box and whisker for shannon diversity

cup <- c('1', '2', '3', '4','5','6','7','8','9','10','11','12','13')
mass_flux<- c(".261",".868",".497",".32",".165",".078",".105",".016",".058","1.03",".149","0.093","0.016")
geochem<-data.frame(cup,mass_flux)
p <- ggplot(rich_all,aes(x=factor(cup, level=c('1', '2', '3', '4','5','6','7','8','9','10','11','12','13')), y=Shannon, width=0.5)) +
  geom_boxplot(width=0.9,aes(x=factor(cup, level=c('1', '2', '3', '4','5','6','7','8','9','10','11','12','13')), y=Shannon,fill=cup)) +
  theme(text = element_text(size=30)) + ylab("16S Shannon Diversity Values") + xlab("Cup") + theme(legend.position="right")+ scale_fill_manual(values=mycolors,breaks=c('1', '2', '3', '4','5','6','7','8','9','10','11','12','13')) +
  geom_point(aes(fill=cup), size =5, shape = 21, colour = "black", position=position_jitterdodge(0.4)) +
  stat_boxplot(geom = "errorbar", width=0, size=1.2) +
  theme(axis.text.x=element_text(angle=0, color="black",size=20)) 
p

#Boxplot for OTU abundance
p1 <- ggplot(rich_all,aes(x=factor(cup, level=c('1', '2', '3', '4','5','6','7','8','9','10','11','12','13')), y=Observed, width=0.5)) +
  geom_boxplot(width=0.9,aes(x=factor(cup, level=c('1', '2', '3', '4','5','6','7','8','9','10','11','12','13')), y=Observed,fill=cup)) +
  theme(text = element_text(size=30)) + ylab("# of Observed ASVs") + xlab("Cup") + theme(legend.position="right")+ scale_fill_manual(values=mycolors,breaks=c('1', '2', '3', '4','5','6','7','8','9','10','11','12','13')) +
  geom_point(aes(fill=cup), size =5, shape = 21, colour = "black", position=position_jitterdodge(0.4)) +
  stat_boxplot(geom = "errorbar", width=0, size=1.2) +
  theme(axis.text.x=element_text(angle=0, size=20, color="black")) 
p1
ggsave(filename = "Diversity_16S.pdf", plot = last_plot(), device = "pdf", path = NULL, scale = 1, width = 10, height =5, dpi = 150) 

# 16S taxa tree plots
ps <- tax_glom(ps_rare, "Genus", NArm=FALSE)
ps0 <- transform_sample_counts(ps, function(x) x / sum(x))
x1 = psmelt(ps0)
focus <- c("Proteobacteria", "Bacteroidota", "Desulfobacterota", "Firmicutes", "Campilobacterota", "Planctomycetota", "Unknown Bacteria (Kingdom)") # Focus on the top 16S groups
x1$Phylum <- ifelse(x1$Phylum %in% focus, x1$Phylum, "Others") # Others category for plotting

# Plot 16S tree map
pdf("treemap_16S.pdf", width = 10, height = 3)
tree<-treemap(dtf = x1,
              title = "16S all samples", 
              algorithm = "pivotSize", border.lwds = c(3,0.5,0.1),
              border.col = c("black", "black", "black"),
              lowerbound.cex.label=0,
              mapping = c(0,0,0),
              index = c("Class", "Genus"),
              vSize = "Abundance",
              vColor = "Phylum",
              palette = "Set2",
              type="categorical",
              fontsize.labels=c(14,10),
              fontface.labels=c(2,2),
              overlap.labels=0.8,
              bg.labels=220, 
              position.legend = "right",
              align.labels=list(
                c("left", "top"), 
                c("left", "bottom"),
                c("center", "center")),
              inflate.labels=F,
              force.print.labels = F) 
dev.off()

# Taxa abundance over time for top 6 groups
class_16S <- tax_glom(ps_rare, taxrank = "Class")
class_16S <- transform_sample_counts(class_16S, function(x)100* x / sum(x)) # Transform to relative abundance
OTU <- otu_table(class_16S)
TAX <- tax_table(class_16S)[,"Class"]
Average <- as.data.frame(rowMeans(OTU))# Estimate average abundance at the class level
names(Average) <- c("Mean") # Rename to mean
Table <- merge(TAX, Average, by=0, all=TRUE)
Table$Row.names = NULL # View abundances at class level - we are interested in the top 6 most relatively abundant groups

# Temporal abundance plots for the major 16S groups
abund_16S <- psmelt(class_16S) # Melt the data
abund_16S$Class <- as.character(abund_16S$Class) # Convert to character
abund_16S_new <- subset(abund_16S, grepl("Gammaproteobacteria|Bacteroidia|Alphaproteobacteria|Campylobacteria|Clostridia|Desulfobacteria", abund_16S$Class)) # Subset to major groups
abund_16S_new$Class <- factor(abund_16S_new$Class,levels=c("Gammaproteobacteria", "Bacteroidia","Alphaproteobacteria", "Campylobacteria", "Clostridia","Desulfobacteria")) # Order groups in the plot
p <- ggplot(abund_16S_new, aes(factor(date_range), Abundance, group=Class))
p$data$date_range <- factor(p$data$date_range, levels = c("12/7-12/14","12/14-12/21","12/21-01/4","01/4-01/18","01/18-02/01","02/01-02/15","02/15-03/01","03/01-03/15", "03/15-03/29","03/29-04/12","04/12-04/26","04/26-05/10","05/10-05/24"))
p + geom_point(aes(fill = Class)) + theme(text = element_text(size=16))+
  geom_smooth(se=TRUE, linetype="solid", size=1, level=0.95, fill = "gray45", color = "black", alpha = 0.4, method = "loess")+ 
  facet_wrap(~ Class, ncol =2, scales = "free_y") + theme_bw() + ylab("Relative Abundance (%)") + 
  theme(legend.position = "none")+ guides(fill=guide_legend(ncol=2))+ theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(filename = "16S_abund.pdf", plot = last_plot(),  path = NULL, scale = 1, width = 7, height = 8, dpi = 150) 

########Additional analysis to those that Sean Anderson performed
# First aglomerate the ASVs at the order level using the phyloseq function, tax_glom
phylumGlommed = tax_glom(ps_rare, "Phylum")
                                     
# and plot
gpt <- subset_taxa(phylumGlommed, Kingdom=="Bacteria")
gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:300]), gpt)
plot_heatmap(gpt, sample.label="cup")
p1 <- ggplot(rich_all,aes(x=factor(cup, level=c('1', '2', '3', '4','5','6','7','8','9','10','11','12','13')), y=Observed, width=0.5)) +
  geom_boxplot(width=0.9,aes(x=factor(cup, level=c('1', '2', '3', '4','5','6','7','8','9','10','11','12','13')), y=Observed,fill=cup)) +
  theme(text = element_text(size=30)) + ylab("# of Observed ASVs") + xlab("Cup") + theme(legend.position="right")+ scale_fill_manual(values=mycolors,breaks=c('1', '2', '3', '4','5','6','7','8','9','10','11','12','13')) +
  geom_point(aes(fill=cup), size =5, shape = 21, colour = "black", position=position_jitterdodge(0.4)) +
  stat_boxplot(geom = "errorbar", width=0, size=1.2) +
  theme(axis.text.x=element_text(angle=0, size=20, color="black")) 
plot_bar(phylumGlommed, fill = "Phylum")+
  facet_grid(rows="cup", scales="free_y",space="free")+coord_flip()+
  ylab("Relative abundance") + theme(axis.title = element_text(size=20, face="bold")) +
  xlab("Phylum") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=15)) + 
  theme(legend.text = element_text(size = 15)) + theme(legend.title = element_text(size = 10, face="bold")) + 
  theme(strip.text.x = element_text(size=15, face="bold")) 
  #theme(legend.key.size = unit(1, "cm"))+
  #theme(legend.key = element_rect(size=10)) + 
  #theme(legend.direction = "vertical",legend.box = "vertical")
bubcolors<- c("#277da1","#604798","#f67894","#13aa8b","#707787", "#f1801a","#70be6d","#f00149")
xx<-ggplot(x1, aes(x = cup , y = x1$Abundance)) + 
  geom_point(aes(size = x1$Abundance , fill = Phylum ), alpha =0.8, shape = 21) +
  guides(fill = guide_legend(override.aes = list(size = 10))) +
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_rect(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") +
  scale_fill_manual(values = mycolors) +
  scale_size_continuous(range = c(0,25),name = "Abundance") 
xx

#Average the triplicate cup values for each OTU for a simplified bubble plot
#subset data
newdata <- x1[c(1,3,4
                )]
#for each unique cup value in newdata, make a new DF with the mean for each OTU
test1<-(aggregate(Abundance ~ OTU + cup,  data=newdata,FUN=mean,na.rm=TRUE))
colnames(test1)[3] ="Average_Abundance"
Avgab<-test1
                                     
#Add taxonomic information based on OTU 
testing<-merge(x1,Avgab)
total <- merge(x1,Avgab,by=c("OTU","cup"))
subset<-total[total$Average_Abundance != 0, ]
avgabund<-subset[,-4]
avgabund<-subset[,-3]

#keep only unique rows
avgabund %>% distinct()
Avgab<-avgabund

#Replot a simplified bubble plot
mycolors <- c("#17385A","#427475", "#0E9594","#1F5041","#C2CFB2","#F9DFBE","#ECCBD3","#827191","#92416C","#52113C","#962C2C","#F2673D","#FF9B71")
mycolors1 <- c("#C2CFB2","#1F5041","#427475","#17385A","#52113C","#92416C","#827191","#F9DFBE")            
bubb1<-ggplot(Avgab, aes(x=factor(cup, level=c('1', '2', '3', '4','5','6','7','8','9','10','11','12','13')) , y = Average_Abundance)) + 
  geom_point(aes(size = Average_Abundance , fill = Phylum ), alpha =0.8, shape = 21) +
  ylim(0.05,.46) +
  scale_y_break(c(0.25,0.30), scales = 0.20) + 
  guides(fill = guide_legend(override.aes = list(size = 15))) +
  labs( x= "", y = "", size = "Average Relative Abundance (%)", fill = "")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 25, face = "bold"), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 25), 
        legend.text = element_text(size = 25, face ="bold", colour ="black"), 
        legend.title = element_text(size = 25, face = "bold"), 
        panel.background = element_rect(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") +
  scale_fill_manual(values = mycolors1) +
  scale_size_continuous(range = c(0,25),name = "Average Abundance") 
bubb1
resultList <- list()
for(x in 1:nrow(x1)){
  key = x1[x, ]$OTU
  value = x1[x, ]$Abundance
  cup = x1[x, ]$cup
  if(key == '8af50705664f5e1b1ed664535738d942'){
    if(!(is.list(resultList[[key]]))){
      resultList[[key]] <- list()
    }
    if(!(is.list(resultList[[key]][[cup]]))){
      resultList[[key]][[cup]] = list()
    }
    resultList[[key]][[cup]] = c(resultList[[key]][[cup]], value)
  }
}
finalResultList <- data.frame()
for(x in 1:length(resultList)){
  for(y in 1:length(resultList[[x]])){
    counter = 0
    total = 0.0
    for(z in 1:length(resultList[[x]][[y]])){
      total = (resultList[[x]][[y]][[z]] + total)
      counter = (counter + 1)
    }
    finalResultList[nrow(finalResultList) + 1,] <- c(resultList[[x]], y, (total / counter))
  }
}
