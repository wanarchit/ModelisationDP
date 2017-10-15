library(ggplot2)
library(stringr)

args = commandArgs(trailingOnly = TRUE)
inputfile = args[1]
outputPath = str_match(inputfile,"(.*\\\\)Results_(.*)\\.")[2]
outputname = str_match(inputfile,"(.*\\\\)Results_(.*)\\.")[3]
outputname = paste("Plot_",outputname,sep="")
outputname = paste(outputPath,outputname,sep="")
outputname = paste(outputname,".png",sep="")

tableTest = read.csv(inputfile, header=TRUE, sep=";")

plotDist <- ggplot(tableTest, aes(x = Cluster, y = Effective, group = Frame_t, color = Frame_t)) 
plotDist <- plotDist + scale_colour_gradientn(colours=rainbow(6))
plotDist <- plotDist + geom_point()
plotDist <- plotDist + geom_smooth(fill = NA, span = 0.1)
plotDist <- plotDist + theme_bw()

png(file = outputname, width = 800, height = 700)
plotDist
dev.off()