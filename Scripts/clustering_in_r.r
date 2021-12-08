library(cluster)
library(factoextra)
library(igraph)
library(readxl)
library(ontologyIndex)
library(ontologySimilarity)
library("writexl")
library(ggplot2)
if (!requireNamespace("BiocManager", quietly = TRUE))
{
  install.packages("BiocManager")
  BiocManager::install("AnnotationDbi")
}

get_non_zero_columns <- function(df)
{
  col_sums <- colSums(df[,1:ncol(df)])
  df <- df[,col_sums>0]
  return(df)
}

df = read.csv("hpo_cols_with_labels.csv",  header = TRUE)

df_external_id <- df$Variant.type
df_pnumber <- df$P.number
df = df[,0:(ncol(df)-2)]

hpo_terms <- colnames(df)
hpo_terms <- sub("\\.", ":",hpo_terms)
hpo_terms <- hpo_terms[grepl("HP*", hpo_terms)]

groepen_hpo <- read_excel("hpo_groups.xlsx")


df_chd3 <- df

for (i in 1:length(groepen_hpo$groep))
{
  if (is.na(groepen_hpo$groep[i]) == FALSE)
  {
    df_chd3 <- fncols(df_chd3,groepen_hpo$groep[i]) 
  }
}

for (i in 1:length(groepen_hpo$groep))
{
  hpo_wrong_format <- groepen_hpo[[i,"hpo_terms"]]
  hpo_right_format <-  sub("\\:", ".",hpo_wrong_format)
  groep_naam <- groepen_hpo[[i,"groep"]]
  if (is.na(groep_naam) == FALSE)
  {
    df_chd3[,groep_naam] <- df_chd3[,groep_naam] + df_chd3[,hpo_right_format]
  }
}

groepen_hpo[is.na(groepen_hpo[,"groep"]), "groep"] <- groepen_hpo[is.na(groepen_hpo[,"groep"]), "hpo_terms"]

selected_columns = gsub("\\:", ".",groepen_hpo$groep)

names.use <- names(df_chd3)[(names(df_chd3) %in% selected_columns)]
df_chd3_groups <- df_chd3[,names.use]

set.seed(123) # for reproducibility
# Compute PAM

df_dn_inher <- get_non_zero_columns(df_chd3_groups)

row.names(df_dn_inher) <- paste0(df_external_id, '_', row.names(df_dn_inher))
rownames(df_dn_inher) <- df_pnumber
pam.res <- pam(df_dn_inher, 2)
Observed <- df_external_id
fviz_cluster(pam.res, xlab='PC1', ylab='PC2',labelsize=14, repel=TRUE,geom = NULL) + theme(axis.text.x = element_text(size = 8)) + geom_point(aes(shape = Observed,colour = factor(pam.res$cluster)))
df_dn_inher$cluster= pam.res$cluster-1

df_dn_inher$labels <- df_external_id
df_dn_inher$labels[df_dn_inher$labels == 'Missense'] = 1
df_dn_inher$labels[df_dn_inher$labels == 'PTV'] = 0


ggsave(file="cluster_plot.eps", dpi=300, width=20, height=10)
ggsave(file="cluster_plot.png", dpi=300, width=20, height=10)

df_dn_inher$labels_reversed <- 1-as.numeric(as.character(df_dn_inher$labels))

obs_score = sum(df_dn_inher$labels == df_dn_inher$cluster)
if (sum(df_dn_inher$labels_reversed == df_dn_inher$cluster) > obs_score)
{
  obs_score = sum(df_dn_inher$labels_reversed == df_dn_inher$cluster)
}

get_random_draw <- function(n1, n_patients) {
  draw <- rbinom(n_patients, 1,(n1/n_patients))
  while (sum(draw) != n1){
    draw <- rbinom(n_patients, 1,(n1/n_patients))
  }
  return(draw)
}

num_trials = 100000
scores_all <- vector("list", num_trials)
pb = txtProgressBar(min = 0, max = num_trials, initial = 0) 
for (i in 1:num_trials)
{ 
  n_patients = nrow(df_dn_inher)
  prior_missense =sum(df_dn_inher[,"labels"] == 1) 
  draw <- get_random_draw(prior_missense, n_patients)
  score = sum(draw == df_dn_inher$labels)
  scores_all[i] = score
  setTxtProgressBar(pb,i)
}

df_to_write <- df_dn_inher[c("cluster","labels","labels_reversed")]
df_to_write["index"] <- df_external_id[df$de_novo_probands_inherited > -1]

write_xlsx(df_to_write,"ptv_missense_with_clusters.xlsx")

scores_all <- as.numeric(scores_all)

hist(as.numeric(scores_all))
max(as.numeric(scores_all))
min(as.numeric(scores_all))
print(mean(as.numeric(scores_all)))
print(obs_score)

print("P value ")
mean(as.numeric(scores_all)>=obs_score)


