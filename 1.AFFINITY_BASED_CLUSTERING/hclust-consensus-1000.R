rawData <- t(read.csv("Drewry-data.csv",header=TRUE))

colnames(rawData) <- 1:dim(rawData)[2]

rawDataBoot <- lapply(1:1000,function(dx) rawData[,sample(colnames(rawData),length(colnames(rawData)),replace=TRUE)])

kinaseBoot_hclust <- lapply(rawDataBoot,function(dx) hclust(dist(dx)))

kinaseBoot_hclust <- lapply(kinaseBoot_hclust, as.phylo)
class(kinaseBoot_hclust) <- 'multiPhylo'
write.tree(kinaseBoot_hclust,file='hclust-consensus-1000.tre')

#the program 'consense' was then used to create a consensus tree from these 1000 trees

hcd <- ReadDendrogram('bootstrap/hclust-1000-rooted/outtree')

pdf(file="Drewry-consensus-three-colors.pdf",width=20, height=20)
fviz_dend(hcd,k = 3, k_colors = "jco",
          type = "circular", repel = TRUE)
dev.off()
