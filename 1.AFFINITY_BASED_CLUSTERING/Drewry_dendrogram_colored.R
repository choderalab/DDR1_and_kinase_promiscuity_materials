rawData <- t(read.csv("Drewry-data.csv"))
kinaseDist <- dist(rawData)
kinaseClust <- hclust(kinaseDist)

inhData<- read.csv("data-Drewry-75inh.csv")
pdf(file="Drewry-dendrogram-colors.pdf",width=40, height=10)
hcd <-as.dendrogram(kinaseClust)
mycols <- grDevices::rainbow(length(inhData))
colLab <- function(n) {
    if(is.leaf(n)) {
      a <- attributes(n)
	  label <- a$label
      attr(n, "nodePar") <-
        c(a$nodePar, list(lab.col = mycols[inhData[ ,label]]))
      attr(n, "frame.plot") <- TRUE
    }
    n
  }
clusDendro = dendrapply(hcd, colLab)
par(cex=0.6, font=2,mar=c(12, 4, 4, 0))
plot(clusDendro, xlab="", ylab="", main="", sub="", pin=c(8,1),mai=rep(0.1,4),axes=FALSE)
par(cex=0.9)
title(xlab="Kinases", ylab="", main="Kinases Clustered by Inhibitor Activity and then colored by # of drugs with >75% Inhibition")
axis(2)
dev.off()