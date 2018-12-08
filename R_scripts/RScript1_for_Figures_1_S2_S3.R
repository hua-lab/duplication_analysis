
#####################

#To run the script, save Data S1 "Dataset1_for_Figures_1_S2_S3.tab" and the file of this script,"RScript1_for_Figures_1_S2_S3.R", in the same directory and type "Rscript RScript1_for_Figures_1_S2_S3.R" in terminal under the same directory to generate all Figures.

#####################

#required packages

library(gplots)
library(RColorBrewer)


#Figures 1A. Size distribution of the ubiquiton family in 50 plant genomes
#Figure 1B. Sizes of top 8 ubiquiton subfamilies
#Figure S2A. Correlation of the size of the ubiquiton family with the genome size of a plant

	data<-read.table("Dataset1_for_Figures_1_S2_S3.tab",header=T)

	d<-data
	d<-d[,-1]
	d<-d[,-dim(d)[2]]

#subfamily_size

	colsum<-colSums(d)
	o_colsum<-colsum[order(colsum)]

#ubl_size_each_species

	rowsum<-rowSums(d)
	o_rowsum<-rowsum[order(rowsum)]

	data<-data[order(data$Genome_Size),]

#correlation between genomesize and ubl size
	trend<-lm( data$Sum ~ data$Genome_Size)
	summary(trend)

######

# Call:
# lm(formula = data$Sum ~ data$Genome_Size)

# Residuals:
#      Min       1Q   Median       3Q      Max 
# -133.963  -32.477   -5.389   18.763  220.463 

# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      8.058e+01  1.137e+01   7.086 5.42e-09 ***
# data$Genome_Size 8.820e-08  1.990e-08   4.433 5.39e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 53.74 on 48 degrees of freedom
# Multiple R-squared:  0.2905,	Adjusted R-squared:  0.2757 
# F-statistic: 19.65 on 1 and 48 DF,  p-value: 5.393e-05



#####

	cor.test(data$Genome_Size,data$Sum, method = "spearman",
 	exact = NULL, co.level = 0.95,)


######
# 	Spearman's rank correlation rho

# data:  data$Genome_Size and data$Sum
# S = 9400.8, p-value = 3.714e-05
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5485808

# Warning message:
# In cor.test.default(data$Genome_Size, data$Sum, method = "spearman",  :
#   Cannot compute exact p-value with ties
# > 

######


	pdf ("Figure 1A_1B_S2A.pdf", family="Times", height=10, width=10)

	par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
	layout(matrix(c(1,2,2,1,2,2,1,3,3),nrow=3,ncol=3))

	barplot(o_rowsum, cex.names = 0.2)
	barplot(o_colsum, cex.names = 0.2)

	plot(data$Genome_Size,data$Sum,,cex=1, pch=16, col=c("blue"))
	text(data$Genome_Size, data$Sum,labels = rownames(data), pos = 1, cex=0.7)
	trend<-lm( data$Sum ~ data$Genome_Size)
	abline(trend, col="blue",lwd=2)
	text(1.5e+9, 350, "rho=0.55, p-value=3.7e-05")

	dev.off()



#############

# Figures 1C. Differential expansion of ubiquiton subfamilies resulted in 4 clusters of plant species 
# Figure S3. Clustering analysis of ubiquiton subfamilies in 50 plant genomes


	d<-read.table("Dataset1_for_Figures_1_S2_S3.tab",header=T)

### wo_ub heatmap

	d<-d[,-1]
	d<-d[,-dim(d)[2]]
	d<-d[,-dim(d)[2]]
	d<-as.matrix(d)

	rowDistance=dist(d,method="manhattan")
	rowCluster = hclust(rowDistance,method="ward.D2")
	rowDend = as.dendrogram(rowCluster)
	rowDend = reorder(rowDend, rowSums(d))

	colDistance=dist(t(d),method="manhattan")
	colCluster = hclust(colDistance,method="ward.D2")
	colDend = as.dendrogram(colCluster)
	colDend = reorder(colDend, colSums(d))

	my_palette<-colorRampPalette(c("light blue","yellow","red"))(n=39)
	col_breaks=c(seq(0,5,length=20),seq(5.5,30,length=15),seq(30.1,max(d),length=5))

	pdf ("Figures 1C and S3_Version1.pdf", family="Times", height=10, width=10)
	par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
	heatmap.2(d,Rowv=rowDend, Colv=colDend,col=my_palette,breaks=col_breaks,keysize=2,margins=c(20,20),trace=c("none"),density.info=c("none"),cexRow=0.5,cexCol=0.5)
	dev.off()

	h<-heatmap.2(d,Rowv=rowDend, Colv=colDend,col=my_palette,breaks=col_breaks,keysize=2,margins=c(20,20),trace=c("none"),density.info=c("none"),cexRow=0.5,cexCol=0.5)

### w_ub heatmap

	d_ub<-read.table("Dataset1_for_Figures_1_S2_S3.tab",header=T)
	d_ub<-d_ub[,-1]
	d_ub<-d_ub[,-dim(d_ub)[2]]
	d_ub<-as.matrix(d_ub)
	d_ub_order<-d_ub[rev(colnames(h$carpet)),]

	colDistance=dist(t(d_ub_order),method="manhattan")
	colCluster = hclust(colDistance,method="ward.D2")
	colDend = as.dendrogram(colCluster)
	colDend = reorder(colDend, colSums(d_ub_order))

	my_palette<-colorRampPalette(c("light blue","yellow","red"))(n=39)
	col_breaks=c(seq(0,25,length=20),seq(25.1,200,length=15),seq(200.1,max(d_ub_order),length=5))

	pdf ("Figures 1C and S3_Version2.pdf", family="Times", height=10, width=10)
	par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
	heatmap.2(d_ub_order,Rowv=FALSE, Colv=colDend, 	dendrogram=c("column"),col=my_palette,breaks=col_breaks,keysize=2,margins=c(20,20),trace=c("none"),density.info=c("none"),cexRow=0.5,cexCol=0.5)
	dev.off()



#########


# Figure 1D. Size correlation analysis of 17 ubiquiton subfamilies

	cor<-cor(d_ub_order, method = c("spearman"))

	rowDistance=dist(cor,method="manhattan")
	rowCluster = hclust(rowDistance,method="ward.D2")
	rowDend = as.dendrogram(rowCluster)
	rowDend = reorder(rowDend, rowSums(cor))

	colDistance=dist(t(cor),method="manhattan")
	colCluster = hclust(colDistance,method="ward.D2")
	colDend = as.dendrogram(colCluster)
	colDend = reorder(colDend, colSums(cor))

	mycol <- colorpanel(n=9,low="blue",mid="light yellow",high="red")

	pdf ("Figure 1D.pdf", family="Times", height=10, width=10)
	par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
	heatmap.2(cor,Rowv=rowDend, Colv=colDend, col=mycol,  keysize=2, margins=c(20,20),trace=c("none"),density.info=c("none"),cexRow=0.5,cexCol=0.5)
	dev.off()



###########
#Figure S2B. Spearman’s correlation relationships between 17 ubiquiton subfamilies and the scales of 50 plant genomes

#get cor data

	d<-read.table("Dataset1_for_Figures_1_S2_S3.tab",header=T)

	cor_test<-c();

	for(i in 1:dim(d)[2]){
		tst<-cor.test(d[,i],d$Genome_Size, method = "spearman",exact = FALSE, co.level = 0.95,)		
		tst_value<-c(colnames(d)[i],tst$p.value,tst$estimate)		
		cor_test<-rbind(cor_test,tst_value)				
		}
		
	write.table(cor_test, file="Figure S2B.tab")







