library(ggplot2)
library(RColorBrewer)
#filePath = "C:/Users/Asus/Documents/Master 2 GPhy/Semestre 1/Modelisation/Results.csv"
filePath = "D:/Cours/Master 2/Modelisation Bioinformatique/p3_p4_p5_p8/Results.csv"
tableTest = read.csv(filePath, header=TRUE, sep=";")
View(tableTest)
#ggplot(tableTest, aes(x = Cluster, y = Effectif, color = Frame_t)) + geom_point() + theme_classic() 

p <- ggplot(tableTest, aes(x = Cluster, y = Effectif, group = Frame_t, color = Frame_t)) 
p <- p + scale_colour_gradientn(colours=rainbow(6))
p <- p + geom_point()
p <- p + geom_smooth(fill = NA, span = 0.1)
p <- p + theme_bw()
p

#frameT <- unique(tableTest$Frame_t)
