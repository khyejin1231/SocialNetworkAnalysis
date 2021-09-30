#https://igraph.org/r/doc/igraph.pdf

install.packages("readxl")
install.packages("linkcomm")
install.packages("igraph")
install.packages("NetData")
install.packages(c("sna","triads","psych",'nFactors','GPArotation','NetCluster'))
install.packages("leiden")
install.packages("RColorBrewer")
install.packages("ExPosition")
install.packages("cppRouting")
install.packages("ggplot2")
library(ExPosition)
library(RColorBrewer)
library(lenden)
library(readxl)
library(linkcomm)
library(igraph)
library(sna)
library(triads) 
library(psych)
library(nFactors)
library(GPArotation)
library(NetCluster)
library(cppRouting)
library(dplyr)
library(ggplot2)


data <- read_excel("Department_Data.xlsx")
data <- data[,2:30]

isSymmetric(as.matrix(data)) #This matrix is not symmetric

############################################################################Data
A <- data.frame("V3" = NA)
colnames(data)[1:29] <- "V3"

for(i in 1:29){
  A <- rbind(data[,i], A, make.row.names = TRUE)
  A <- na.omit(A)

}

A$V1 <- rep(1:29,each=29)
A$V2 <- rep(1:29, 29)
A <- subset(A, A$V3>0)
A <- A[,c("V1", "V2", "V3")]

######################################################################Question1.
data <- as.matrix(data)
d <- as.dist(data)
?as.dist
g1 <- graph_from_adjacency_matrix(as.matrix(data), mode = c("directed"), 
                                  weighted = TRUE, diag = FALSE)
cc <- hclust(d, method = "ward.D") #Which method?
plot(cc)
clusters.list <- rect.hclust(cc, k= 2, border = "blue")
clusters <- cutree(cc, k= 2)
node.cols <- brewer.pal(max(c(3,clusters)), "Pastel1")[clusters]
plot(g1,  vertex.color = node.cols,layout = coords)


https://igraph.org/r/doc/sample_pa.html

fc <- cluster_fast_greedy(g1)
plot_dendrogram(fc)

######################################################################node-based
rownames(data) <- c(1:29)
colnames(data) <- c(1:29)

g1 <- graph_from_adjacency_matrix(as.matrix(data), mode = c("directed"), 
                                  weighted = TRUE, diag = FALSE)
coords <- layout_with_fr(g1)
plot(g1, layout=coords, vertex.label = NA, vertex.size = 5)

#Which clusterings to use?
#Greedy community detection, Spectral community detection, betewenness community detection
#Hierarchical clustering, Optimal community
B <- get.adjacency(g1, sparse = FALSE)

#What is modularity?

#question 1. 
#a. 2,3,4,5,6,7 clusters
# 2 cluster
?cluster_optimal

partition <- leiden(B)
table(B)
node.cols <- brewer.pal(max(c(3,partition)), "Pastel1")[partition]
plot(g1, vertex.color = node.cols)
#Who is the leader?

#b. Matrix normalization?????
data_norm <- data
row_sum <- apply(data, 1, function(x)(sum(x)))
data_norm <- (data[1,]/row_sum[1])
for (i in 1:29){
  data_norm[i,] <- data[i,]/row_sum[i]
}


g2 <- graph_from_adjacency_matrix(as.matrix(data_norm), mode = c("directed"), 
                                  weighted = TRUE, diag = FALSE)

C <- get.adjacency(g2, sparse = FALSE)

partition <- leiden(C)

node.cols <- brewer.pal(max(c(3,partition)), "Pastel1")[partition]
plot(g2, vertex.color = node.cols)

#question 2.
#a.
cc <- hclust(data, method = "ward") #Which method?
plot(cc)
clusters.list <- rect.hclust(cc, k= 5, border = "blue")
clusters <- cutree(cc, k= 5)
node.cols <- brewer.pal(max(c(3,clusters)), "Pastel1")[clusters]
plot(g1,  vertex.color = node.cols,layout = coords)

#b. link
A$V1 <- as.numeric(A$V1)
A$V2 <- as.numeric(A$V2)
A$V3 <- as.numeric(A$V3)
A <- as.matrix(A)

?getLinkCommunities
lc <- getLinkCommunities(A, hcmethod = "ward" ,directed = TRUE)

#question3.

##############################################################Dijkstra algorithm
A[,3] <- 1/A[,3] #First, we inverse the weights
head(A)
colnames(A) <- c("from", "to", 'weight')
graph <- makegraph(A, directed = T, coords = NULL)
graph$nbnode
graph$dict$ref

get_distance_matrix(graph, 1, 2, algorithm = "mch", allcores = TRUE)
get_path_pair(graph, 22,1,algorithm = "Dijkstra")

######################################################Edge.betweenness.community
edge_cluster <- cluster_edge_betweenness(g1, weights = E(g1)$weight, directed = TRUE,
                         edge.betweenness = TRUE, merges = TRUE, bridges = TRUE,
                         modularity = FALSE, membership = TRUE)

membership(edge_cluster)
cut <- cutat(edge_cluster, 2)
colors <- rainbow(10)
plot(g1, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=8)


cut <- cutat(edge_cluster, 3)
plot(g1, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=8)

cut <- cutat(edge_cluster, 4)
plot(g1, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=8)

cut <- cutat(edge_cluster, 5)
plot(g1, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=8)

cut <- cutat(edge_cluster, 6)
plot(g1, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=8)

cut <- cutat(edge_cluster, 7)
plot(g1, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=8)

plot_dendrogram(edge_cluster)

cut <- cutat(edge_cluster, 7)
plot(g1,edge.width = 0, layout = layout_in_circle, 
     vertex.color = colors[cut], xlim = c(-3.5, 3.5),
     ylim = c(-3.5, 3.5),vertex.size=degree(g1, mode="all")*0.5, edge.color= "white", 
     vertex.label=NA)
par(new=TRUE)



##norm
edge_cluster <- cluster_edge_betweenness(g2, weights = E(g2)$weight, directed = TRUE,
                                         edge.betweenness = TRUE, merges = TRUE, bridges = TRUE,
                                         modularity = FALSE, membership = TRUE)



membership(edge_cluster)==1
cut <- cutat(edge_cluster, 2)
colors <- rainbow(10)
plot(g2, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=8, edge.width = log(edge.betweenness(g2)))

cut <- cutat(edge_cluster, 3)
plot(g2, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=8)

cut <- cutat(edge_cluster, 4)
plot(g2, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=8)

cut <- cutat(edge_cluster, 5)
plot(g2, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=8)

cut <- cutat(edge_cluster, 6)
plot(g2, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=8)

cut <- cutat(edge_cluster, 7)
plot(g2, vertex.color = colors[cut],  edge.arrow.size=0.5, vertex.size=degree(g2, mode="all")/2,
     xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))

plot_dendrogram(edge_cluster)
#image size: 1128,825
(membership(edge_cluster)==1)
#visualization_norm
cut <- cutat(edge_cluster, 2)
plot(g2,edge.width = 0, layout = layout_in_circle, 
     vertex.color = colors[cut], xlim = c(-0.75, 0.75),
     ylim = c(-0.75, 0.75),vertex.size=degree(g2, mode="all")*0.5, edge.color= "white", 
     vertex.label=NA)
par(new=TRUE)
cut <- cutat(edge_cluster, 3)
plot(g2,edge.width = 0, layout = layout_in_circle, 
     vertex.color = colors[cut], xlim = c(-1, 1),
     ylim = c(-1, 1),vertex.size=degree(g2, mode="all")*0.5, edge.color= "white", 
     vertex.label=NA)
par(new=TRUE)
cut <- cutat(edge_cluster, 4)
plot(g2,edge.width = 0, layout = layout_in_circle, 
     vertex.color = colors[cut], xlim = c(-1.4, 1.4),
     ylim = c(-1.4, 1.4),vertex.size=degree(g2, mode="all")*0.5, edge.color= "white", 
     vertex.label=NA)
par(new=TRUE)
cut <- cutat(edge_cluster, 5)
plot(g2,edge.width = 0, layout = layout_in_circle, 
     vertex.color = colors[cut], xlim = c(-2, 2),
     ylim = c(-2, 2),vertex.size=degree(g2, mode="all")*0.5, edge.color= "white", 
     vertex.label=NA)
par(new=TRUE)
cut <- cutat(edge_cluster, 6)
plot(g2,edge.width = 0, layout = layout_in_circle, 
     vertex.color = colors[cut], xlim = c(-2.7, 2.7),
     ylim = c(-2.7, 2.7),vertex.size=degree(g2, mode="all")*0.5, edge.color= "white", 
     vertex.label=NA)
par(new=TRUE)
cut <- cutat(edge_cluster, 7)
plot(g2,edge.width = 0, layout = layout_in_circle, 
     vertex.color = colors[cut], xlim = c(-3.5, 3.5),
     ylim = c(-3.5, 3.5),vertex.size=degree(g2, mode="all")*0.5, edge.color= "white", 
     vertex.label=NA)



#subgraph
g3 <- induced_subgraph(g2, which(cutat(edge_cluster, 7)==1))
plot(g3,edge.width = log(edge.betweenness(g3)), layout = layout_in_circle, vertex.color = colors[1], xlim = c(-1, 1),
     ylim = c(-1, 1),vertex.size=degree(g3, mode="all")/2, edge.color= "white")
par(new=TRUE)

centr_degree(g3, mode = c("all"), normalized = TRUE)$res #node level


########################################################################walktrap
# Mean recurrence time of a markov chain
# Graph --> P--> M
# How to P --> M: Formula
# will I get block in M?
# Do we take average and make M into a symmetric matrix?
#


source("Mean_recurrence.R")
install.packages(c("expm","matlib"))
library(expm)
library(matlib)
X <- as.matrix(data_norm)
Y <- MRT_distance(X, 10)
Y <- as.dist(Y, diag=TRUE)
hc <- hclust(Y)
plot(hc)

#####################################################################Question 2.
g1 <- graph_from_adjacency_matrix(as.matrix(data), mode = c("directed"), 
                                  weighted = TRUE, diag = FALSE)
coords <- layout_with_fr(g1)
plot(g1, layout=coords, vertex.label = NA, vertex.size = 5)

#Which clusterings to use?
#Greedy community detection, Spectral community detection, betewenness community detection
#Hierarchical clustering, Optimal community
B <- get.adjacency(g1, sparse = FALSE)

#What is modularity?

#question 1. 
#a. 2,3,4,5,6,7 clusters
# 2 cluster
?cluster_optimal

partition <- leiden(B)
table(B)
node.cols <- brewer.pal(max(c(3,partition)), "Pastel1")[partition]
plot(g1, vertex.color = node.cols)
#Who is the leader?

#b. Matrix normalization?????
data_norm <- data
row_sum <- apply(data, 1, function(x)(sum(x)))
data_norm <- (data[1,]/row_sum[1])
for (i in 1:29){
  data_norm[i,] <- data[i,]/row_sum[i]
}


g2 <- graph_from_adjacency_matrix(as.matrix(data_norm), mode = c("directed"), 
                                  weighted = TRUE, diag = FALSE)

C <- get.adjacency(g2, sparse = FALSE)

partition <- leiden(C)

