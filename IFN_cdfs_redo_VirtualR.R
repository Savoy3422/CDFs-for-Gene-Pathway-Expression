### CDF of IFN signal calculations 11-19-19
## ref: https://www.nature.com/articles/s41590-019-0386-1#data-availability
library(plyr)
library(dplyr)
library(limma) # for DE analysis and "avereps" function
library(openxlsx)
library(ComplexHeatmap)
require("Seurat") # used for reading single cell data
library(tidyverse)
require("ggplot2")
library(reshape)
options(bitmapType='cairo') # for X11 to work so PNG is recognized



## Load data
data<-readRDS("/Volumes/cgm/sarc/scrna/analyses/20190919_analysis_dont_use/SCT_normalized_data_63K.Rds")

### Cell types
cell.types=unique(data$celltype)
cell.types=cell.types[-c(4,13,17,18,20,21)]

### Gene list of IFN genes
#IFN=read.table("/Volumes/cgm/sarc/scrna/ifn_inducible.txt")
IFN=read.table("/s/cgm/Savoy/IFN-response-genes.txt")
IFN=as.vector(IFN$V1)
IFN=IFN[which(IFN %in% rownames(data))]


# Remove subjects
keep=read.table("/s/cgm/Savoy/keep_samples.txt")
data2=data # create copy since it takes awhile to import




### Locating Ubiquitous Gene set
all.genes=data.frame()
## Find genes that are expressed a high percent of all cell types (using a higher count sum to indicate its exsistence in many patient cells)
for (cell.type in cell.types){
  genes.1=as.data.frame(Matrix::rowSums(GetAssayData(data[,which(data$celltype==cell.type & data$library_id %in% keep$V1)],assay="SCT",slot="data")!=0))# count number nonzero
  genes=(genes.1/ncol(data[["SCT"]]@counts[,data$celltype==cell.type]))*100
  genes$name=rownames(genes)
  colnames(genes)=c("percent","genes")
  genes2=genes[which(genes$percent >=30),2]
  if(cell.type=="CD4 Memory T"){
    all.genes=genes2
  }else{
    all.genes=intersect(all.genes,genes2)
  }
}
## looking for low variance genes
for (cell.type in cell.types){
  averages=as.data.frame(Matrix::rowMeans(GetAssayData(data[all.genes,which(data$celltype==cell.type)],assay="SCT",slot="data")))
  averages$genes=rownames(averages)
  if(cell.type == "CD4 Memory T"){
    all.avg=averages
  }else{
    all.avg=rbind(all.avg,averages)
  }
}
all.avg2=all.avg[order(all.avg$genes),]
colnames(all.avg2)=c("mean","genes")
all.var=ddply(all.avg2, .(genes), summarise, variance=var(mean))


all.genes.0=all.genes[!all.genes %in% IFN] # removing IFN genes
# save both sets of genes for future reference
sets=cbind(all.genes.0,IFN) 
write.xlsx(sets,"/s/cgm/Savoy/CDF_gene_set_lists.xlsx")


# Use caution when calling GetAssayData. supply assay and slot arguments to make sure we're pulling the values we want. 
# If I use SCT/data, then I shouldn't have to worry about normalizing from cell to cell - it's already taken care of

# The average expression from case cells and controls cells for the two gene sets
cdf.ifn<-data.frame()
cdf.ubiquitous=data.frame()
for (cell.type in cell.types){
  # get the average of the intereferon genes and unbiquitous genes across all cells of this type, for these patients.
  ifn.case<- Matrix::rowMeans(GetAssayData(data[IFN,which(data$celltype==cell.type & data$stim=="Case")], assay="SCT",slot="data"))
  ifn.control=Matrix::rowMeans(GetAssayData(data[IFN,which(data$celltype==cell.type & data$stim=="Control")], assay="SCT",slot="data"))
  ubiquitous.case=Matrix::rowMeans(GetAssayData(data[all.genes.0,which(data$celltype==cell.type & data$stim=="Case")], assay="SCT",slot="data"))
  ubiquitous.control=Matrix::rowMeans(GetAssayData(data[all.genes.0,which(data$celltype==cell.type & data$stim=="Control")], assay="SCT",slot="data"))
  # store the result by averaging all the average ifn scores
  cdf.ifn<-rbind(cdf.ifn, data.frame(Case.IFN.avg= ifn.case, Control.IFN.avg = ifn.control,cell.type=cell.type))
  cdf.ubiquitous=rbind(cdf.ubiquitous,data.frame(Case.Ubiquitous.avg=ubiquitous.case, Control.Ubiquitous.avg=ubiquitous.control, cell.type=cell.type))
}



##log2((normalized.counts(case)+1)/(normalized.counts(control)+1))
cdf.ifn2=data.frame(Gene=row.names(cdf.ifn),Log.Ratio=log((cdf.ifn$Case.IFN.avg+1)/(cdf.ifn$Control.IFN.avg+1),2),cell.type=cdf.ifn$cell.type)
cdf.ubiquitous2=data.frame(Gene=row.names(cdf.ubiquitous),Log.Ratio=log((cdf.ubiquitous$Case.Ubiquitous.avg+1)/(cdf.ubiquitous$Control.Ubiquitous.avg+1),2),cell.type=cdf.ubiquitous$cell.type)

## Plot cdfs
png("/s/cgm/Savoy/IFN_cdfs.png",width=4000,height=4000,res=300)
# start plot
par(mfrow=c(5,3))
for(cell.type in cell.types){
  plot(ecdf(cdf.ifn2[which(cdf.ifn2$cell.type==cell.type),2]), verticals=TRUE,do.points=FALSE,col="red",xlab="log2{[normalized.counts(Sarc)+1]/[normalized.counts(Control)+1]}",
       ylab="Cumulative Distribution",main=cell.type)
  plot(ecdf(cdf.ubiquitous2[which(cdf.ubiquitous2$cell.type==cell.type),2]), verticals=TRUE,do.points=FALSE, add=TRUE)
}
dev.off()

##Obtain cumulative probabilities table
# define fuction
cumprob=function(y){
  fun=function(y,x) length(y[y<x])/length(y)
  prob=sapply(y,fun,y=y)
  data=data.frame(value=unique(y[order(y)]),prob=unique(prob[order(prob)]))
}
# Obtain tables
OUT = createWorkbook()
for(cell.type in cell.types){
  cp.ifn=cumprob(cdf.ifn2[which(cdf.ifn2$cell.type==cell.type),2])
  addWorksheet(OUT, cell.type)
  writeData(OUT, sheet=cell.type, x=cp.ifn, rowNames=FALSE)
}# end cell.type
saveWorkbook(OUT,"/s/cgm/Savoy/IFN_cdf_values.xlsx")

## Obtain pvalues for curve differences
ks.p=data.frame()
for (cell.type in cell.types){
  ks=ks.test(cdf.ifn2[which(cdf.ifn2$cell.type==cell.type),2], cdf.ubiquitous2[which(cdf.ubiquitous2$cell.type==cell.type),2], alternative = "two.sided")
  ks.p<-rbind(ks.p, data.frame(P.value= ks$p.value, Cell.Type=cell.type))
} # end cell.type
write.xlsx(ks.p,"/s/cgm/Savoy/CDF_comparison_pvalues.xlsx")
