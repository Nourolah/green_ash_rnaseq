##-----------------------------------------------------------------
## dendro.R is an Rscript used to 
## normalize htseq counts and generate clustered output in
## "De novo assembly of the green ash transcriptome and 
## identification of genes responding to abiotic and biotic stress"
## Lane et. al. 2015
##
## dendro.R assumes a directory in home called "dendro" exists
## dendro.R also assumes a sub-directory in home called counts-sorted
## containing the output of htseq runs accross 55 green ash libraries.
##
## The rlog step takes a long time to complete.
##
## dendro.R will create a single output Tiff "dendro.tiff"
##

# set up R environment
library("DESeq2")
library("pheatmap")
setwd("~/dendro/counts-sorted")
directory<-"~/dendro/counts-sorted"

# load in files as a list
sampleFiles <- list.files()

# create list of branch names
sampleName<-c("Leaf, cold then recovery ",
              "Petiole, cold then recovery",
              "Root, cold then recovery",
              "Leaflets, control",
              "Leaflets, mechanical wounding",
              "Twigs, control",
              "Twigs, mechanical wounding",
              "Xylem",
              "Seeds with wings",
              "Axial buds",
              "Terminal buds",
              "Leaf, ozone + wounding - 80ppb, 29 days",
              "Leaf, ozone + wounding - 125ppb, 29 days",
              "Leaf, ozone + wounding - 225ppb, 29 days",
              "Leaf, ozone + wounding - control, 29 days",
              "Leaf, cold ",
              "Leaf, control",
              "Petiole, control",
              "Root, control",
              "Petiole, cold",
              "Root, cold",
              "Leaf, drought",
              "Petiole, drought",
              "Root, drought",
              "Leaf, heat",
              "Petiole, heat",
              "Root, heat",
              "Leaf, mechanical wounding, 24hrs",
              "Leaf, mechanical wounding, 5hrs",
              "Petiole, mechanical wounding, 24hrs",
              "Petiole, mechanical wounding, 5hrs",
              "Leaf, ozone - control, 7 hours",
              "Leaf, ozone - 80ppb, 28 days",
              "Leaf, ozone - 125ppb, 28 days",
              "Leaf, ozone - 225ppb, 28 days",
              "Leaf, ozone - 80ppb, 7 hours",
              "Leaf, ozone - 125ppb, 7 hours",
              "Leaf, ozone - 225ppb, 7 hours",
              "Leaf, ozone - control, 14 days",
              "Leaf, ozone - 80ppb, 14 days",
              "Leaf, ozone - 125ppb, 14 days",
              "Leaf, ozone - 225ppb, 14 days",
              "Leaf, ozone - control, 28 days",
              "Bark and phloem Tree 19 control",
              "Bark and phloem Tree 19 after EAB feeding",
              "Bark and phloem Tree 21 control",
              "Bark and phloem Tree 21 after EAB feeding",
              "Bark and phloem Tree 22 control",
              "Bark and phloem Tree 22 after EAB feeding",
              "Bark and phloem Tree 24 control",
              "Bark and phloem Tree 24 after EAB feeding",
              "Bark and phloem Tree 36 control",
              "Bark and phloem Tree 36 after EAB feeding",
              "Bark and phloem Tree Summit control",
              "Bark and phloem Tree Summit after EAB feeding")

# create data frame for input into DESeq2
sampleTable<-data.frame(
sampleName=sampleName,
filename=sampleFiles)

# load htseq count files into deseq2 object
ddsNOFACTOR <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                        directory=directory,
                                        design= ~ 1)

# create deseq2 object
ddsNOFACTOR <- DESeq(ddsNOFACTOR)

# run rlog normalization, takes a long time to run on this data
rldNOFACTOR <- rlog(ddsNOFACTOR)

# create distance array of normalized values
sampleDistsNOFACTOR <- dist(t(assay(rldNOFACTOR)))

# transform array into matrix
sampleDistsMatrixNOFACTOR <- as.matrix(sampleDistsNOFACTOR)

# draw dendrogram using pheatmap
setwd("~/dendro")
tiff("dendro.tiff", width = 10, height = 10, units = 'in', res = 300)
pheatmap(sampleDistsMatrixNOFACTOR,
         clustering_distance_rows=sampleDistsNOFACTOR,
         clustering_distance_cols=sampleDistsNOFACTOR,
         treeheight_col=0,
         cellwidth=0,
         legend=FALSE,
         show_colnames=FALSE,
         treeheight_row=300)
dev.off()
