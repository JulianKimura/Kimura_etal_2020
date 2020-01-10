#Making a heatmap of TPM values
library(gplots)

#import TPM matrix
mean_TPM<-read.table("/Users/JulianKimura/Documents/Lab/Broseq/Broseq_2019_analysis/heatmap/Pairwise_atleast1_sig/wardD_clust/qval_filt_sig_pair_gene_TPM.txt",header = TRUE, stringsAsFactors=FALSE, row.names = 1)
#convert to matrix and create distance metric
mat_mean_TPM <- as.matrix(mean_TPM)
getVar <- apply(mat_mean_TPM[, -1], 1, var)
param <- 0
mean_TPM_variance_trimmed <- mat_mean_TPM[getVar > param & !is.na(getVar), ]
hc <- hclust(as.dist(1 - cor(t(mean_TPM))), method="ward.D")

#determine the height of the tree to cut by to create clusters
mycl <- cutree(hc, h=max(hc$height/30))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]

#picking colors
my_palette <- colorRampPalette(c("royalblue","white", "chocolate"))(n=64)

#plotting heatmap
heatmap.2(as.matrix(mean_TPM), dendrogram = c("row"), Rowv = as.dendrogram(hc), Colv = FALSE, labRow = "", col = my_palette, margins = c(10,12), scale = "row", trace = "none", RowSideColors = myClusterSideBar)
