library(enrichR)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library('enrichR')
library(clusterProfiler)
library(org.Mm.eg.db)



#BiocManager::install("enrichplot")
#install.packages("ggpubr")
#aenrichmentDotPlot(enriched)
#enrichR::plotEnrich(enriched)
#enrichR::printEnrich(enriched)
#asd = enrichR::printEnrich(enriched)
#
#dbs <- c("KEGG_2019_Mouse")
#enriched <- enrichr(clus0, dbs)
#plotEnrich(enriched$KEGG_2019_Mouse)


mapping = select(org.Mm.eg.db, clus0, columns = "ENTREZID", keytype = "SYMBOL")
mapping = mapping['ENTREZID']
de <- mapping$ENTREZID
genlist = c()
for (i in de) {
  genlist = c(genlist, as.integer(i))
}

ego <- enrichGO(na.omit(genlist), OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE)
barplot(ego)

b1 = barplot(mapp(clus0), showCategory = 3) 
b2 = barplot(mapp(clus1), showCategory = 3)
b3 = barplot(mapp(clus2), showCategory = 3)
b4 = barplot(mapp(clus3), showCategory = 3)
b5 = barplot(mapp(clus4), showCategory = 3)
b6 = barplot(mapp(clus5), showCategory = 3)
b7 = barplot(mapp(clus6), showCategory = 3)
b8 = barplot(mapp(clus7), showCategory = 3)
b9 = barplot(mapp(clus8), showCategory = 3)

barplot(ego) + theme(axis.text.y =element_text(size=15,face="bold"),
                           axis.title=element_text(size=20,face="bold"),
                           legend.text = element_text(size=20))

ego
ggarrange(b1,b2,b3,b4,b6,b7,b8,b9,
          ncol = 3, nrow = 3)

mapp <- function(cluster)
{
  
  mapping = select(org.Mm.eg.db, cluster, columns = "ENTREZID", keytype = "SYMBOL")
  mapping = mapping['ENTREZID']
  de <- mapping$ENTREZID
  genlist = c()
  for (i in de) {
    genlist = c(genlist, as.integer(i))
  }
  
  return(enrichGO(genlist, OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE))
}
