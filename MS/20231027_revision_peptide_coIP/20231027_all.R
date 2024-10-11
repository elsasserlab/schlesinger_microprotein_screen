#trying out reading in data



library("DEP")
library("dplyr")
library("shiny")
library("tidyr")
 

setwd("~/Desktop/Mass_Spec/20231027_revision_peptide_coIP")
run_app("LFQ") #supposedly for shiny but never really worked...

#20231027----
# load data
my_data <- read.delim("proteinGroups.txt")
data <-my_data
colnames(data) # use this to make experimental design file
experimental_design <- read.delim("20231027_experiment_all.txt")

# We filter for contaminant proteins and decoy database hits, which are indicated by "+" in the columns "Potential.contaminants" and "Reverse", respectively. 
data <- filter(data, Reverse != "+", Potential.contaminant != "+")


dim(data)
## [1] 455 418


# Are there any duplicated gene names?
data$Gene.names %>% duplicated() %>% any()
## [1] TRUE
# Make a table of duplicated gene names
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
#A tibble: 2 × 2
#Gene.names frequency
#<chr>          <int>
#  1 ""                 9
#2 "HNRNPR"           2
#For further analysis these proteins must get unique names. Additionally,
#some proteins do not have an annotated gene name and for those we will use the Uniprot ID.

# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Are there any duplicated names?
data$name %>% duplicated() %>% any()
## [1] FALSE

# Generate a SummarizedExperiment object using an experimental design
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se <- make_se(data_unique, LFQ_columns, experimental_design)

#Generate a SummarizedExperiment object by parsing condition information from the column names
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
data_se_parsed <- make_se_parse(data_unique, LFQ_columns)

# Let's have a look at the SummarizedExperiment object
data_se
#class: SummarizedExperiment 
#dim: 455 48 
#metadata(0):
#  assays(1): ''
#rownames(455): EFCAB7 C17orf100 ... CFAP20 IGF2BP2
#rowData names(372): Protein.IDs Majority.protein.IDs ... name ID
#colnames(48): parental_1 parental_2 ... cand19_5 cand19_6
#colData names(4): label ID condition replicate
# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
#data_filt <- filter_missval(data_se, thr = 0)
#filter_missval
# OR FOR Less stringent filtering:
#Filter for proteins that are identified in 3 out of 6 replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 3)


#filter_missval filters a proteomics dataset based on missing values.
#The dataset is filtered for proteins that have a maximum of ’thr’ missing values in at least one condition.


# Less stringent filtering:
# Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
#data_filt2 <- filter_missval(data_se, thr = 1) here have only 2 replicates so shouldn't need this
#After filtering, the number of identified proteins per sample can be plotted as well as the overlap in identifications between samples.

# Plot a barplot of the number of identified proteins per samples
plot_numbers(data_filt)

# Plot a barplot of the protein identification overlap between samples
plot_coverage(data_filt)

#The data is background corrected and normalized by variance stabilizing transformation (vsn).

# Normalize the data
data_norm <- normalize_vsn(data_filt)
#The normalization can be inspected by checking the distributions of the samples before and after normalization.

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm)
#-> looks worse than without I think, so don't normalize

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(data_filt, fun = "")


# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
#data_imp <- impute(data_norm, fun = "MinProb", q = 0.01) #with norm
#data_imp <- impute(data_norm, fun = "minDet") #with norm
#more explanation: https://rdrr.io/bioc/MSnbase/man/impute-methods.html


#OR 

#without
#data_imp <- impute(data_filt, fun = "MinProb", q = 0.01) #without norm
data_imp <- impute(data_filt, fun = "min") #without norm
#data_imp <- impute(data_filt, fun = "MinDet") #without norm
#data_imp <- impute(data_filt, fun = "zero") #without norm

#data_imp <- impute(data_filt, fun = "MinProb", q = 0.01) #without norm
# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
#data_imp_man <- impute(data_filt, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
#data_imp_knn <- impute(data_filt, fun = "knn", rowmax = 0.9)

# Plot intensity distributions before and after imputation
#plot_imputation(data_norm, data_imp) #with norm
#OR 
plot_imputation(data_filt, data_imp) #without norm

# Plot intensity distributions before and after imputation
#plot_imputation(data_norm, data_imp_man)

# Plot intensity distributions before and after imputation
#plot_imputation(data_norm, data_imp_knn)


# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "parental")
##Tested contrasts: cand16_vs_parental
# Test all possible comparisons of samples
#data_diff_all_contrasts <- test_diff(data_imp, type = "all")
# Test manually defined comparisons
#data_diff_manual <- test_diff(data_imp, type = "manual", 
                         #     test = c("RPL11_vs_cand1", "cand1_vs_RPL11"))
## Tested contrasts: Ubi4_vs_Ctrl, Ubi6_vs_Ctrl
##Finally, significant proteins are defined by user-defined cutoffs using add_rejections.


# Test every sample versus control without imputation--> didn't work, gave error
#data_diff <- test_diff(data_norm, type = "control", control = "parental")
##Tested contrasts: cand16_vs_parental, RPL11_vs_parental
# Test all possible comparisons of samples
#data_diff_all_contrasts <- test_diff(data_norm, type = "all")

# Denote significant proteins based on user defined cutoffs, vs control
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))
#OR
# Denote significant proteins based on user defined cutoffs, vs all
#dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = log2(1.5)) 

# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 54, point_size = 4)
## Warning: Use of `pca_df[[indicate[1]]]` is discouraged. Use
## `.data[[indicate[1]]]` instead.
## Warning: Use of `pca_df[[indicate[2]]]` is discouraged. Use
## `.data[[indicate[2]]]` instead.

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 3, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 3, col_limit = 10, show_row_names = FALSE)

# Plot a volcano plot for the contrast " vs Ctrl""
pdf('20231027_allvsparental_filt3out6reps_NOnorm_imp_min_lfc_1_5_v2.pdf')
plot_volcano(dep, contrast = "GFP_link_vs_parental", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "cand1_vs_parental", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "cand2_vs_parental", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "cand11_vs_parental", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "cand12_vs_parental", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "cand16_vs_parental", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "cand19_vs_parental", label_size = 2, add_names = TRUE)
dev.off()


# Generate a results table
data_results <- get_results(dep)

#save the results table
write.table(data_results, "20231027_allvsparental_filt3out6reps_NOnorm_imp_min_lfc_1_5.txt", sep="\t")
