library(ggplot2)
library(stringr)

# Retrieving of the PATH of input CSV file : counting table
args = commandArgs(trailingOnly = TRUE)
inputfile = args[1]

# Building of output PNG file name with the execution program date and hour
outputPath = str_match(inputfile,"(.*\\\\)Results_(.*)\\.")[2]
outputname = str_match(inputfile,"(.*\\\\)Results_(.*)\\.")[3]
outputname = paste("Plot_",outputname,sep="")
outputname = paste(outputPath,outputname,sep="")
outputname = paste(outputname,".png",sep="")

# Loading of data in data.frame : 3 columns : time frame, cluster and effective.
tableTest = read.csv(inputfile, header=TRUE, sep=";")

# Building of plot with cluster distance on the x-axis and effective on the y-axis
plotDist <- ggplot(tableTest, aes(x = Cluster, y = Effective, group = Frame_t, color = Frame_t))
# Coloring of the different groups
plotDist <- plotDist + scale_colour_gradientn(colours=rainbow(6))
plotDist <- plotDist + geom_point()
# Adding a trend line for each time frame
plotDist <- plotDist + geom_smooth(fill = NA, span = 0.1)
# Deleting of backgroung grey color and the color around the trand line
plotDist <- plotDist + theme_bw()
# Adding legend on the plot
plotDist <- plotDist + ggtitle("Distances traveled by reference atoms according to the range of time")
plotDist <- plotDist + xlab("Distance (A)")
plotDist <- plotDist + ylab("Effective (distance, time)")

# Saving of output (plot) in PNG format
png(file = outputname, width = 800, height = 700)
plotDist
dev.off()