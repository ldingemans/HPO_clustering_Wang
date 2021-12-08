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

fncols <- function(data, cname) {
  add <-cname[!cname%in%names(data)]
  
  if(length(add)!=0) data[add] <- 0
  data
}


library(HPOSim)
customWang <- function(x,y){ 
  z <- getSimWang(x,y)
  return(z)
}

df = read.csv("\\\\umcfs020\\ANTRGdata$\\Genetica Projecten\\Facial Recognition\\Studenten en Onderzoekers\\Lex\\Projecten\\ANKRD11_Elke\\Data\\hpo_cols_with_labels.csv",  header = TRUE)

df_external_id <- df$Variant.type
df_pnumber <- df$P.number
df = df[,0:(ncol(df)-2)]

hpo_terms <- colnames(df)
hpo_terms <- sub("\\.", ":",hpo_terms)
hpo_terms <- hpo_terms[grepl("HP*", hpo_terms)]

sim_mat <- outer(hpo_terms, hpo_terms, Vectorize(customWang))

list_of_vars_similar_hpo <- list()
list_of_vars_similar <- list()

for (i in 1:length(hpo_terms)){
  list_of_vars_similar_hpo[[i]] <- vector()
  list_of_vars_similar[[i]] <- vector()
  for (y in 1:length(hpo_terms))
  {
    if (sim_mat[y, i] > 0.5 & y != i)
    {
      list_of_vars_similar_hpo[[i]] <- c(list_of_vars_similar_hpo[[i]],hpo_terms[y])
      list_of_vars_similar[[i]] <- c(list_of_vars_similar[[i]],y)
    }
  }
}

df_hpo_groups <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("index", "hpo_terms"))

for (y in 1:length(list_of_vars_similar_hpo))
{
  df_hpo_groups[nrow(df_hpo_groups) + 1,] = y
  df_hpo_groups[nrow(df_hpo_groups),2] = toString(list_of_vars_similar_hpo[y][1])
}

write.csv(df_hpo_groups,"similar_hpo.csv")

dfr1 <- data.frame(source=numeric(), target=numeric())

for (i in 1:length(hpo_terms)){
  for (y in 1:length(list_of_vars_similar[[i]]))
  {
    dfr1[nrow(dfr1) + 1,] = c(i,list_of_vars_similar[[i]][y])
  }
}

dfr1 <- na.omit(dfr1)

net <- graph_from_data_frame(d=dfr1, directed=T) 

facial_congenital_hpo <- read_excel("list_facial_hpo.xlsx")

hpo_names <- vector("list", 1)

for (i in 1:length(hpo_terms))
{
  if (grepl("HP", hpo_terms[i], fixed=TRUE))
  {
    name <- facial_congenital_hpo[facial_congenital_hpo$hpo_id == hpo_terms[i],]$hpoName
    hpo_names[i] <- name
  }
}

plot(simplify(net), 
     vertex.label = hpo_names,
     edge.width = 1,
     edge.arrow.width = 0.3,
     vertex.size = 5,
     edge.arrow.size = 0.5,
     vertex.size2 = 3,
     vertex.label.cex = 0.8,
     asp = 0.35,
     margin = 0.001)

ggsave("hpo_terms_cluster.png", width = 20, height = 20, units = "cm")

df_hpo <- cbind(hpo_names,tail(hpo_terms, length(hpo_names)))
write.table(df_hpo,"hpo_names.csv", sep=';', row.names = FALSE)
