rm(list = ls())
# Load the necessary packages
library(tidyverse)
library(stringr)
################################################################################
# Information from dbSNPs
################################################################################
# Load the data
snp_data<- read.delim("C:/Users/ksanders5/Desktop/BMI5330/Final Project/snp_result_tab.txt", header = T)
View(snp_data)
## Cleaning the dataset ##
clean_snp_df<- snp_data %>% filter(variation != "variation")
# Now I will add another column to segregate the mutation type
clean_snp_df$mutation_type <- ifelse(clean_snp_df$variant_type %in% c("snv", "mnv"), "SNP", "indel")
# Next I will create three datasets.
# This is for SNP
snp_only <- clean_snp_df %>%
  filter(mutation_type == "SNP")
indel_only <- clean_snp_df %>%
  filter(mutation_type == "indel")
### Summarize counts and findings
summary_table_mutation <- clean_snp_df %>%
  group_by(mutation_type) %>%
  summarize(
    count = n(),
    percentage = n() / nrow(clean_snp_df) * 100
  ) %>%
  arrange(desc(count))

print(summary_table_mutation)
###############################################################
# SNP summary
summary_table_snp <- snp_only %>%
  group_by(variant_type,variation) %>%
  summarize(
    count = n(),
    percentage = n() / nrow(snp_only) * 100
  ) %>%
  arrange(desc(count))

print(summary_table_snp)
view(summary_table_snp)
# Visualize the SNPs

# Visualization
# Create a new column that extracts the letter before the ">" symbol
summary_table_snp <- summary_table_snp %>% mutate(base = substring(variation, 1, 2))
head(summary_table_snp)

ggplot(data = summary_table_snp, aes(x = reorder(variation, count), y = count, fill = base))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Dark2")+
  coord_flip()+
  labs(title = "Distribution of Nucleotide Variants", x = "Polymorphs", y = "Occurrence", 
       fill = "Base Nucleotide(s)", caption = "Data source: https://www.ncbi.nlm.nih.gov/snp/ ")+
  theme_classic()
# Mutation types
snp_func <- snp_only %>%
  filter(grepl("(coding_sequence_variant|missense_variant|stop_gained|stop_loss)", function_class))

snp_func_count <- snp_func %>%
  mutate(function_class = str_extract_all(function_class, "(missense_variant|stop_gained|synonymous_variant|stop_loss)")) %>%
  filter(length(function_class) > 0) %>%
  unnest(function_class) %>%
  group_by(variation, function_class) %>%
  summarize(count = n())
head(snp_func_count)

# Create a new column that extracts the letter before the ">" symbol
snp_func_count <- snp_func_count %>% mutate(base = substring(variation, 1, 2))
# Visualization

ggplot(data = snp_func_count, aes(x = reorder(function_class, count), y = count, fill = base))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_brewer(palette = "Dark2")+
  labs(title = "Distribution of Variations across Function Classes", x = "SNP Functional Classes", y = "Occurrence", 
       fill = "Base Nucleotide(s)", caption = "Data source: https://www.ncbi.nlm.nih.gov/snp/ ")+
  theme_classic()

snp_func_count <- snp_func %>%
  filter(grepl("(missense_variant|stop_gained|synonymous_variant|stop_loss)", function_class)) %>%
  group_by(variation, function_class) %>%
  summarize(count = n())
head(snp_func_count)

# indels
head(indel_only)
### Summarize counts and findings
summary_table_indel <- indel_only %>%
  group_by(variant_type) %>%
  summarize(
    count = n(),
    percentage = n() / nrow(indel_only) * 100
  ) %>%
  arrange(desc(count))
head(summary_table_indel)

# Visualization
ggplot(data = summary_table_indel, aes(x = reorder(variant_type, count), y = count, fill = variant_type))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Set2")+
  coord_flip()+
  labs(title = "Distribution of Indel Variants", x = "Variation Type", y = "Occurrence", 
       fill = "Indels Type", caption = "Data source: https://www.ncbi.nlm.nih.gov/snp/ ")+
  theme_classic()

indel_func <- indel_only %>%
  filter(grepl("(coding_sequence_variant|inframe_insertion|inframe_deletion|frameshift_variant)", function_class))

indel_func_count <- indel_func %>%
  mutate(function_class = str_extract_all(function_class, "(inframe_insertion|inframe_deletion|frameshift_variant)")) %>%
  filter(length(function_class) > 0) %>%
  unnest(function_class) %>%
  group_by(variant_type, function_class) %>%
  summarize(count = n())
head(indel_func_count)

# Visualization
ggplot(data = indel_func_count, aes(x = reorder(function_class, count), y = count, fill = variant_type))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_fill_brewer(palette = "Set2")+
  labs(title = "Distribution of Indel Variants Across Function Classes", x = "Functional Class Type", y = "Occurrence", 
       fill = "Indels Type", caption = "Data source: https://www.ncbi.nlm.nih.gov/snp/ ")+
  theme_classic()
################################################################################
# dbVar
################################################################################
# Read the dbVar Data
dbVar_data<- read_csv("dbVar.csv")
head(dbVar_data)
dbVar_Origin <- read_csv("dbVar_Origins.csv")
head(dbVar_Origin)
# Visualization
# Create vector of colors
colors <- c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#00FFFF", "#FF00FF", "#800000", "#008000", "#000080", "#808000", "#800080")
# dbVar
ggplot(data = dbVar_data, aes(x = reorder(Variation_Type, Occurance), y = Occurance, fill = Variation_Type))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colors) +
  coord_flip()+
  labs(title = "Distribution of Structural Variants", x = "Variation Type", y = "Occurrence", 
       fill = "Variation Type", caption = "Data source: https://www.ncbi.nlm.nih.gov/dbvar/ ")+
  theme_classic()
# db Origins
colors2 <- c("#0072B2", "#009E73", "#D55E00", "#CC79A7", "#F0E442", "#56B4E9", "#800000", "#999999", "#000000", "#FF7F00", "#8B008B", "#008000","#FF0000")
ggplot(data = dbVar_Origin, aes(x = reorder(Origin, Count), y = Count, fill = Origin))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = colors2) +
  coord_flip()+
  labs(title = "Distribution of the Origins of Structural Variants", x = "Origin", y = "Occurrence", 
       fill = "Origins", caption = "Data source: https://www.ncbi.nlm.nih.gov/dbvar/ ")+
  theme_classic()
################################################################################
# ClinVar
################################################################################
## Read in ClinVar data
clinvar_data<- read.delim("C:/Users/ksanders5/Desktop/BMI5330/Final Project/clinvar_result.txt", header = T)
head(clinvar_data)
clean_clinvar <- clinvar_data %>%
  select(Gene.s.,Protein.change,Condition.s., Clinical.significance..Last.reviewed.) %>%
  filter(str_detect(Clinical.significance..Last.reviewed., regex("pathogen", ignore_case = TRUE)))
write_csv(clean_clinvar, "clean_clinvar.csv")
clean_clinvar<-read_csv("clean_clinvar.csv")
summary_clean_clinvar <- clean_clinvar %>%
  group_by(Condition.s.) %>%
  summarize(count = n())
  

# Visualize
colors3 <- c("#F0E442","#00FF00","#008000","#00FFFF","#56B4E9","#0072B2","#FF0000","#800000",
            "#800080",'#A93226',"#76448A","#E74C3C","#D55E00","#000080","#6D4C41","#BCAAA4",
            "#3E2723","#808000","#FBC02D","#F57F17","#999999","#C2185B","#000000")
ggplot(summary_clean_clinvar, aes(x = reorder(Condition.s., -count), y = count, fill = Condition.s.))+
  geom_bar(stat = "identity", color = "black")+
  labs(title = "Frequency of Pathogenic Conditions Related to APOE Variants", x = "Condition", y = "Occurrence", 
       fill = "Pathogenic Conditions", caption = "Data source: https://www.ncbi.nlm.nih.gov/clinvar/ ")+
  scale_fill_manual(values = colors3) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
#################################################################################
# Manhattan Plot
################################################################################
AD_GWAS_sample<- read_delim("C:/Users/ksanders5/Desktop/RdmGWAS/GWAS_AD_Wightman_exclude_UKBB_23andme.txt")
#Read in the AD_GWAS data frame and filter out X and Y chromosomes
data <- AD_GWAS_sample %>% filter(chr %in% 1:22)

#Trying to remove 0 pvalues
AD_GWAS_sample2 <- AD_GWAS_sample
AD_GWAS_sample2 <- subset(AD_GWAS_sample2, pvalue != 0)
# Group SNPs by chromosome and position
data <- data %>%
  group_by(chr) %>%
  mutate(position = cumsum(snp_pos) - snp_pos[1] + 1)

# Calculate -log10(p-value)
data$pval <- -log10(data$pvalue)

# Set significance threshold (p < 0.05)
threshold <- -log10(0.05)

# Define order of chromosomes
chr_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22")

# Convert "chr" column to factor with desired levels
data$chr <- factor(data$chr, levels = chr_order)

# Create Manhattan plot
ggplot(data, aes(x = position, y = pval, color = factor(chr))) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("#CC0000", "#0066CC", "#FFCC00", "#33CC33", "#CC33CC", "#FF6600", "#999999", "#000000", "#FF33CC", "#00CC99", "#9999CC", "#CCCC99", "#FF9933", "#66CCCC", "#669933", "#660066", "#FFCC99", "#00FF00", "#993333", "#33CCFF", "#9900CC", "#CC9999", "#FF6666", "#FF6632"))+
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "top") +
  xlab("Chromosome") +
  ylab("-log10(p-value)") +
  ggtitle("GWAS Manhattan plot") +
  geom_hline(yintercept = threshold, color = "red")+
  scale_x_discrete(labels = 1:22)
###############################################################################
# Sample Manhattanplot 
library(tidyverse)
library(ggtext)
library(normentR)
install.packages("janitor")
library(janitor)
devtools::install_github("norment/normentR")
set.seed(2404)

gwas_data_load <- simulateGWAS(nSNPs = 1e6, nSigCols = 3) %>% 
  janitor::clean_names()
#
sig_data <- gwas_data_load %>% 
  subset(p < 0.05)
notsig_data <- gwas_data_load %>% 
  subset(p >= 0.05) %>%
  group_by(chr) %>% 
  sample_frac(0.1)
gwas_data <- bind_rows(sig_data, notsig_data)
#
data_cum <- gwas_data %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(bp)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(chr, bp_add)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)
#
axis_set <- gwas_data %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(p == min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)

sig <- 5e-8
## Plot the data
manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(p), 
                                  color = as_factor(chr), size = -log10(p))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )
manhplot
print(manhplot)
# Lets try with my data
#
sig_data <- AD_GWAS_sample2 %>% 
  subset(pvalue < 0.05)
notsig_data <- AD_GWAS_sample2 %>% 
  subset(pvalue >= 0.05) %>%
  group_by(chr) %>% 
  sample_frac(0.1)
gwas_data <- bind_rows(sig_data, notsig_data)
#
#
data_cum <- gwas_data %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(snp_pos)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(chr, bp_add)

gwas_data <- gwas_data %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = snp_pos + bp_add)
#
axis_set <- gwas_data %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

ylim <- gwas_data %>% 
  filter(pvalue == min(pvalue)) %>% 
  mutate(ylim = abs(floor(log10(pvalue))) + 2) %>% 
  pull(ylim)
#
sig <- 5e-8
manhplot <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(pvalue), 
                                  color = as_factor(chr), size = -log10(pvalue))) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$chr)))) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 75, size = 8, vjust = 0.5)
  )
manhplot
print(manhplot)