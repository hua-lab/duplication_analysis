
#####################

#To run the script, save Data S3 "Dataset3_for_Figures_4_S7_S8.tab" and the file of this script,"RScript3_for_Figures_4_S7_S8.R", in the same directory and type "Rscript RScript3_for_Figures_4_S7_S8.R" in terminal under the same directory to generate all Figures.

#####################

d<-read.table("Dataset3_for_Figures_4_S7_S8.tab",header=T);

wgd<-d[(d[,5]!="tandem"),]

wgd<-wgd[(wgd[,5]!="unknown"),]

total<-d

total_species<-substring(total[,1],1,3)
total_species_distribution<-table(total_species)

total_wgd_species<-substring(wgd[,1],1,3)
total_wgd_species_distribution<-table(total_wgd_species)

#Figure 4A. Spearman’s correlation test between the number of WGD-derived ubiquiton genes and the total number of ubiquiton genes in 23 plant genomes.

	distribution<-rbind(total_species_distribution,total_wgd_species_distribution)

	cor.test(distribution[1,],distribution[2,],method = "spearman",exact = FALSE, co.level = 0.95)

		# 	Spearman's rank correlation rho

		# data:  distribution[1, ] and distribution[2, ]
		# S = 84.846, p-value = 7.116e-13
		# alternative hypothesis: true rho is not equal to 0
		# sample estimates:
		#       rho 
		# 0.9580799 
		

#Figure 4B. The differential retention of ancient WGD-derived ubiquiton genes played an opposite role in the expansion of the plant ubiquiuton family.  

###Cretaceous_wgd

	Cretaceous_wgd<-wgd[(wgd[,5]=="Cretaceous"),]
	Cretaceous_wgd_species<-substring(Cretaceous_wgd[,1],1,3)
	Cretaceous_wgd_species_distribution<-table(Cretaceous_wgd_species)
	
	distribution2<-rbind(total_species_distribution,Cretaceous_wgd_species_distribution)

	cor.test(distribution2[1,],distribution2[2,],method = "spearman",exact = FALSE, co.level = 0.95)
		
		# 	Spearman's rank correlation rho

		# data:  distribution2[1, ] and distribution2[2, ]
		# S = 2955.4, p-value = 0.02714
		# alternative hypothesis: true rho is not equal to 0
		# sample estimates:
		#        rho 
		# -0.4601879 

#Figure 4C. K-Pg WGDs played a random role in the expansion of the plant ubiquiuton family.  

##### Paleo
	
	Paleo_wgd<-wgd[(wgd[,5]=="paleo"),]
	Paleo_wgd_species<-substring(Paleo_wgd[,1],1,3)
	Paleo_wgd_species_distribution<-table(Paleo_wgd_species)
	
	total_species_distribution_paleo<-total_species_distribution[names(total_species_distribution) %in% names(Paleo_wgd_species_distribution)]
	total_species_distribution_paleo<-total_species_distribution_paleo[order(names(Paleo_wgd_species_distribution))]
	Paleo_wgd_species_distribution<-Paleo_wgd_species_distribution[order(names(Paleo_wgd_species_distribution))]

	distribution3<-rbind(total_species_distribution_paleo,Paleo_wgd_species_distribution)

	cor.test(distribution3[1,],distribution3[2,],method = "spearman",exact = FALSE, co.level = 0.95)

		
		# 	Spearman's rank correlation rho

		# data:  distribution3[1, ] and distribution3[2, ]
		# S = 484.93, p-value = 0.6339
		# alternative hypothesis: true rho is not equal to 0
		# sample estimates:
		#       rho 
		# 0.1340483 

#Figure 4D. Recent WGDs positively increased the size of the ubiquiton family in 3 plant genomes.

##neo distribution
	
	neo_wgd<-wgd[(wgd[,5]=="neo"),]
	neo_wgd_species<-substring(neo_wgd[,1],1,3)
	neo_wgd_species_distribution<-table(neo_wgd_species)
	
	total_species_distribution_neo<-total_species_distribution[names(total_species_distribution) %in% names(neo_wgd_species_distribution)]
	total_species_distribution_neo<-total_species_distribution_neo[order(names(neo_wgd_species_distribution))]
	neo_wgd_species_distribution<-neo_wgd_species_distribution[order(names(neo_wgd_species_distribution))]

	distribution4<-rbind(total_species_distribution_neo,neo_wgd_species_distribution)

	cor.test(distribution4[1,],distribution4[2,],method = "spearman",exact = FALSE, co.level = 0.95)
		
		# 	Spearman's rank correlation rho

		# data:  distribution4[1, ] and distribution4[2, ]
		# S = 0, p-value < 2.2e-16
		# alternative hypothesis: true rho is not equal to 0
		# sample estimates:
		# rho 
		#   1 
		

#Figure 4. WGDs played a predominant role in the expansion of the ubiquiton family in plants

	pdf("Figure 4.pdf",family="Times", width=10,height=10)

	par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)

	layout(matrix(c(1,2,3,4), nrow=2,ncol=2))

	plot(distribution[2,],distribution[1,],cex=2, pch=19, col="blue", ylim=c(40,160), xlim=c(0,160), xlab="total_wgd_species_distribution", ylab="total_species_distribution")
	text(distribution[2,],distribution[1,],labels = names(total_wgd_species_distribution), pos = 1, cex=0.5)
	trend<-lm( distribution[1,] ~ distribution[2,])
	abline(trend, col="blue",lwd=2)
	text(100,140,"rho=0.96,p-value=7.116e-13")

	plot(distribution2[2,],distribution2[1,],cex=2, pch=19, col="blue", ylim=c(40,160), xlim=c(0,100), xlab="Cretaceous_wgd_species_distribution ", ylab="total_species_distribution")	
	text(distribution2[2,],distribution2[1,],labels = names(total_wgd_species_distribution), pos = 1, cex=0.5)
	trend<-lm( distribution2[1,] ~ distribution2[2,])
	abline(trend, col="blue",lwd=2)
	text(60,140,"rho=-0.4601879 ,p-value=0.02714")

	plot(distribution3[2,],distribution3[1,],cex=2, pch=19, col="blue", ylim=c(40,160), xlim=c(0,100), xlab="Paleo_wgd_species_distribution ", ylab="total_species_distribution")
	text(distribution3[2,],distribution3[1,],labels = names(total_species_distribution_paleo), pos = 1, cex=0.5)
	trend<-lm( distribution3[1,] ~ distribution3[2,])
	abline(trend, col="blue",lwd=2)
	text(60,140,"rho=0.1340483,p-value=0.6339")

	plot(distribution4[2,],distribution4[1,],cex=2, pch=19, col="blue", ylim=c(40,160), xlim=c(0,150), xlab="neo_wgd_species_distribution ", ylab="total_species_distribution")
	text(distribution4[2,],distribution4[1,],labels = names(total_species_distribution_neo), pos = 1, cex=0.5)
	trend<-lm( distribution4[1,] ~ distribution4[2,])
	abline(trend, col="blue",lwd=2)
	text(60,140,"rho=1,p-value< 2.2e-16")

	dev.off()


	#####
		
		
#Figure S7. Differential retention rates of ancient WGD-derived ubiquiton genes in 23 plant genomes.


	pdf("Figure S7.pdf",family="Times", width=10,height=10)

	barplot(distribution2, col=c("gray","light blue"),main="number_of_ubls",
		xlab="Species",ylim=c(0,150), beside=TRUE,cex.axis=0.3, cex.names=0.4,)

	dev.off()

	retention_rate<-distribution2[2,]/distribution2[1,]*100
	retention_rate<-as.matrix(retention_rate)

	species_only_w_cretaceous<-c("Atr","Csa","Csi","Fve","Rco","Tca","Vvi")
	species_only_w_cretaceous_retention_rate<-retention_rate[rownames(retention_rate) %in% species_only_w_cretaceous, ]
	species_only_w_cretaceous_retention_rate<-as.matrix(species_only_w_cretaceous_retention_rate)

	species_w_other_wgds_retention_rate<-retention_rate[!(rownames(retention_rate) %in% species_only_w_cretaceous), ]
	species_w_other_wgds_retention_rate<-as.matrix(species_w_other_wgds_retention_rate)

	wilcox.test(species_only_w_cretaceous_retention_rate[,1],species_w_other_wgds_retention_rate[,1], alternative="greater")

		#	Wilcoxon rank sum test

		#	data:  species_only_w_cretaceous_retention_rate[, 1] and species_w_other_wgds_retention_rate[, 1]
		#	W = 112, p-value = 4.079e-06
		#	alternative hypothesis: true location shift is greater than 0

		mean(species_only_w_cretaceous_retention_rate)
		#[1] 75.90478

		mean(species_w_other_wgds_retention_rate)
		#[1] 27.93418


#####		
		

#Figure S8. Spearman’s correlation test between the number of tandem-duplicated ubiquiton genes and the size of the ubiquiton family in 23 plant genomes.

##tandem
	
	tandem<-d[(d[,5]=="tandem"),]
	tandem_species<-substring(tandem[,1],1,3)
	tandem_species_distribution<-table(tandem_species)
	
	total_species_distribution_neo<-total_species_distribution[names(total_species_distribution) %in% names(tandem_species_distribution)]
	total_species_distribution_neo<-total_species_distribution_neo[order(names(tandem_species_distribution))]
	tandem_species_distribution<-tandem_species_distribution[order(names(tandem_species_distribution))]

	distribution5<-rbind(total_species_distribution_neo,tandem_species_distribution)

	cor.test(distribution5[1,],distribution5[2,],method = "spearman",exact = FALSE, co.level = 0.95)


		# 	Spearman's rank correlation rho

		# data:  distribution5[1, ] and distribution5[2, ]
		# S = 654.93, p-value = 0.006426
		# alternative hypothesis: true rho is not equal to 0
		# sample estimates:
		#      rho 
		# 0.5747229 
		

	pdf("Figure S8.pdf",family="Times", width=10,height=10)

	plot(distribution5[2,],distribution5[1,],cex=2, pch=19, col="blue", ylim=c(40,160), xlim=c(0,30), xlab="tandem_dup_species_distribution ", ylab="total_species_distribution")
	text(distribution5[2,],distribution5[1,],labels = names(total_species_distribution_neo), pos = 1, cex=0.5)
	trend<-lm( distribution5[1,] ~ distribution5[2,])
	abline(trend, col="blue",lwd=2)
	text(20,60,"rho=0.57,p-value=0.006426")

	dev.off()








