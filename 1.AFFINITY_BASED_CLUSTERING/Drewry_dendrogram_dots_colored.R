rawData <- t(read.csv("Drewry-data.csv"))

library('ape')

kinaseDist <- dist(rawData)
kinaseClust <- hclust(kinaseDist)
my_tree = as.phylo(kinaseClust)

write.tree(my_tree,file='hclust_tree.tre')

inhData<- read.csv("data-Drewry-90inh.csv")

mycols <- grDevices::rainbow(406)

tipcol <- rep('black', length(my_tree$tip.label))
for (i in 1:length(my_tree$tip.label)){tipcol[grep(my_tree$tip.label[i],my_tree$tip.label)] <- mycols[inhData[ ,my_tree$tip.label[i]] ]}

setEPS()
postscript("Drewry-dendrogram-rainbow-real-labels.eps")
plot(my_tree, "u", 
     use.edge.length = FALSE,
	 tip.color = tipcol,
	 edge.width = 0.1,
	 cex=0.05)
dev.off()

my_tree = as.phylo(kinaseClust)

tipcol <- rep('black', length(my_tree$tip.label))
for (i in 1:length(my_tree$tip.label)){tipcol[grep(my_tree$tip.label[i],my_tree$tip.label)] <- mycols[inhData[ ,my_tree$tip.label[i]] ]}

dotLabels <- rep('.', length(my_tree$tip.label))
my_tree_dotLabels<- my_tree
my_tree_dotLabels$tip.label <- dotLabels

setEPS()
postscript("Drewry-dendrogram-rainbow-dot-labels.eps")
plot(my_tree_dotLabels, "u", 
     use.edge.length = FALSE,
	 tip.color = tipcol,
	 edge.width = 1.5,
	 cex=3)
dev.off()