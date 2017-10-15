library(ggplot2)
args = commandArgs(trailingOnly = TRUE)
#filePath = "D:/Cours/Master 2/Modelisation Bioinformatique/p3_p4_p5_p8/Results.csv"
tableTest = read.csv(args[1], header=TRUE, sep=";")
#View(tableTest)
#ggplot(tableTest, aes(x = Cluster, y = Effective, color = Frame_t)) + geom_point() + theme_classic() 

p <- ggplot(tableTest, aes(x = Cluster, y = Effective, group = Frame_t, color = Frame_t)) 
p <- p + scale_colour_gradientn(colours=rainbow(6))
p <- p + geom_point()
p <- p + geom_smooth(fill = NA, span = 0.1)
p <- p + theme_bw()
#p
png(file = "toto.png", width = 800, height = 700)
p
dev.off()


#frameT <- unique(tableTest$Frame_t)
