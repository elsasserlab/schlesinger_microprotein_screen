#trying out reading in data


library("DEP")
library("dplyr")
library("shiny")

run_app("LFQ") #supposedly for shiny but never really worked...
setwd("~/Desktop/Mass_Spec/20210702:20210730/Analysis")
#20210730  whole script normalized----
# load data
my_data <- read.delim("DS_proteinGroups.txt")
data <-my_data
colnames(data) # use this to make experimental design file
experimental_design <- read.delim("20210730_experiment_hexuplicates.txt")

# We filter for contaminant proteins and decoy database hits, which are indicated by "+" in the columns "Potential contaminants" and "Reverse", respectively. 
data <- filter(data, Reverse != "+", Potential.contaminant != "+")

dim(data)
[1] 1773  322


# Are there any duplicated gene names?
data$Gene.names %>% duplicated() %>% any()
## [1] TRUE
# Make a table of duplicated gene names
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
## A tibble: 16 × 2
##Gene.names frequency
##<chr>          <int>
##  1 ""                21
##2 "ATAD3A"           2
##3 "HLA-A"            2
##4 "HMGA1"            2
##5 "HNRNPA3"          2
##6 "HNRNPD"           2
##7 "HNRNPK"           2
##8 "HNRNPR"           2
##9 "ILF3"             2
##10 "IMMT"             2
##11 "LMNA"             2
##12 "PCBP2"            2
##13 "RPL11"            2
##14 "SERPINE2"         2
##15 "TMPO"             2
##16 "YBX3"             2
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
##class: SummarizedExperiment 
##dim: 1773 36 
##metadata(0):
##  assays(1): ''
##rownames(1773): >MMP24OS UBA6 ... IGF2BP2 SEC23IP
##rowData names(288): Protein.IDs Majority.protein.IDs ... name ID
##colnames(36): Ctrl_1 Ctrl_2 ... MMP24OS_5 MMP24OS_6
##colData names(4): label ID condition replicate

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)

# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)

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
#does not look super good after normalization 

# Visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm) #normalization does not look very good :/

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt)

# All possible imputation methods are printed in an error, if an invalid function name is given.
impute(data_norm, fun = "")

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
#data_imp <- impute(data_norm, fun = "MinProb", q = 0.01) #with norm
#data_imp <- impute(data_norm, fun = "min") #with norm


data_imp <- impute(data_filt, fun = "min") #without norm
# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
#data_imp_man <- impute(data_filt, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
#data_imp_knn <- impute(data_filt, fun = "knn", rowmax = 0.9)

# Plot intensity distributions before and after imputation
plot_imputation(data_filt, data_imp)

# Plot intensity distributions before and after imputation
#plot_imputation(data_norm, data_imp_man)

# Plot intensity distributions before and after imputation
#plot_imputation(data_norm, data_imp_knn)


# Differential enrichment analysis  based on linear models and empherical Bayes statistics

# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Ctrl")
##Tested contrasts:cand1_vs_Ctrl, cand11_vs_Ctrl, cand12_vs_Ctrl, cand19_vs_Ctrl, MMP24OS_vs_Ctrl
# Test all possible comparisons of samples
#data_diff_all_contrasts <- test_diff(data_imp, type = "all")
## Tested contrasts: cand1_vs_cand16, cand1_vs_Ctrl, cand1_vs_RPL11, cand16_vs_Ctrl, cand16_vs_RPL11, Ctrl_vs_RPL11
# Test manually defined comparisons
#data_diff_manual <- test_diff(data_imp, type = "manual", 
                              #test = c("RPL11_vs_cand1", "cand1_vs_RPL11"))
## Tested contrasts: Ubi4_vs_Ctrl, Ubi6_vs_Ctrl
##Finally, significant proteins are defined by user-defined cutoffs using add_rejections.

# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.5))


# Plot the first and second principal components
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4) #actually looks quite good in comparison to the other data one's
## Warning: Use of `pca_df[[indicate[1]]]` is discouraged. Use
## `.data[[indicate[1]]]` instead.
## Warning: Use of `pca_df[[indicate[2]]]` is discouraged. Use
## `.data[[indicate[2]]]` instead.

# Plot the Pearson correlation matrix
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = TRUE,
             indicate = c("condition", "replicate"))

# Plot a heatmap of all significant proteins (rows) and the tested contrasts (columns)
plot_heatmap(dep, type = "contrast", kmeans = TRUE, 
             k = 6, col_limit = 10, show_row_names = TRUE)

# Plot a volcano plot"
pdf('20210730__plotted230301_hexuplicates_MSMSfalse_NOnorm_imp_min_lfc_1.5.pdf')
par(mfrow=c(1,5))
plot_volcano(dep, contrast = "MMP24OS_vs_Ctrl", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "cand1_vs_Ctrl", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "cand11_vs_Ctrl", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "cand12_vs_Ctrl", label_size = 2, add_names = TRUE)
plot_volcano(dep, contrast = "cand19_vs_Ctrl", label_size = 2, add_names = TRUE)
dev.off()


# Generate a results table
data_results <- get_results(dep)

#save the results table
write.table(data_results, "20210730__plotted230301_hexuplicates_MSMSfalse_NOnorm_imp_min_lfc_1.5.txt", sep="\t")

