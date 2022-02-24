#=====================================================================================
#
#  Code chunk 0
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "D:/Work/Diabetes_Beta_cell/Analysis";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read data set
library("readxl")
Data = read_excel("Data_coexpress.xlsx",1);


datExpr = as.data.frame(t(Data[, -c(1:2)]));
names(datExpr) = Data$GeneIDs;
rownames(datExpr) = names(Data)[-c(1:2)];

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# traitdata

traitData = read_excel("Traits.xlsx");

# remove columns that hold information we do not need.
allTraits = traitData;
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

Samples = rownames(datExpr);
traitRows = match(Samples, allTraits$Samples);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits$Samples;

collectGarbage();

save(datExpr, datTraits, file = "dataInput.RData")




#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================



# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


net = blockwiseModules(datExpr, power = 10,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM", 
                       verbose = 3)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "networkConstruction-auto.RData")
table(moduleColors)
table(net$colors)
table(moduleLabels)
#write.csv(moduleLabels,"module_numbers.csv")
#write.csv(moduleColors,"module_colors.csv")



#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "networkConstruction-auto.RData");
lnames


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);




sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



col1=c("Navy","Navy","Navy","Navy","Navy","Navy","Navy","Navy","Navy","Navy","#B22222","#B22222","#B22222","#B22222","#B22222","#B22222","#B22222","#B22222","#B22222","#B22222")
which.module1="turquoise"
ME1=MEs[, paste("ME",which.module1, sep="")]
which.module2="blue"
ME2=MEs[, paste("ME",which.module2, sep="")]
which.module3="brown"
ME3=MEs[, paste("ME",which.module3, sep="")]
which.module4="yellow"
ME4=MEs[, paste("ME",which.module4, sep="")]


barplot(ME2, col=col1, ylim=c(-0.4,0.6), main="Module 2", cex.main=1.35,
        ylab="ME expression", font.axis=1.25, cex.lab=1.5, font.lab=2,
        xpd=TRUE, las=1, lwd=1, cex.names=0.5, cex.axis=1.25)



par(mfrow=c(2,2))
barplot(ME1, col=col1, ylim=c(-0.4,0.6), main="Module 1", cex.main=1.5,
        ylab="ME expression")

barplot(ME2, col=col1, ylim=c(-0.4,0.6), main="Module 2", cex.main=1.5,
        ylab="ME expression")

barplot(ME3, col=col1, ylim=c(-0.4,0.9), main="Module 3", cex.main=1.5,
        ylab="ME expression")

barplot(ME4, col=col1, ylim=c(-0.4,0.6), main="Module 4", cex.main=1.5,
        ylab="ME expression")



tiff(file="ME1.tiff",
     width=6, height=5, units="in", res=600)
barplot(ME1, col=col1, ylim=c(-0.4,0.6),
        xpd=TRUE, las=1, lwd=1, cex.names=0.5, cex.axis=1.3)
dev.off()

tiff(file="ME2.tiff",
     width=6, height=5, units="in", res=600)
barplot(ME2, col=col1, ylim=c(-0.4,0.6),
        xpd=TRUE, las=1, lwd=1, cex.names=0.5, cex.axis=1.3)
dev.off()

tiff(file="ME3.tiff",
     width=6, height=5, units="in", res=600)
barplot(ME3, col=col1, ylim=c(-0.4,0.8),
        xpd=TRUE, las=1, lwd=1, cex.names=0.5, cex.axis=1.3)
dev.off()

tiff(file="ME4.tiff",
     width=6, height=5, units="in", res=600)
barplot(ME4, col=col1, ylim=c(-0.4,0.6),
        xpd=TRUE, las=1, lwd=1, cex.names=0.5, cex.axis=1.3)
dev.off()








