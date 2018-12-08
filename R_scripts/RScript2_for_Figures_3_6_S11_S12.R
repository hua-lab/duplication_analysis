
#####################

#To run the script, save Data S2 "Dataset2_for_Figures_3_6_S11_S12.tab" and the file of this script,"RScript2_for_Figures_3_6_S11_S12.R", in the same directory and type "Rscript RScript2_for_Figures_3_6_S11_S12.R" in terminal under the same directory to generate all Figures.

#####################

d<-read.table("Dataset2_for_Figures_3_6_S11_S12.tab",header=T) #OrthoMCL numbers 4912 UBLS

dim(d)[1] 

		#4912 total UBLs with Ks and Ka/Ks values

D<-d[(substring(d$ID,4,4)=="_"),]

dim(D)[1] 

		#4550 originally annotated UBLs with clear intron information

new<-d[(substring(d$ID,4,4)!="_"),]

dim(new)[1] 

	#362 new annotated with no clear intron information

#required packages

library(gplots)
library(RColorBrewer)
library(plotly)


#Nomenclature of 6 OrthoMCL clusters

	#cluster1="Green Plants"
	#cluster2="Land Plant I"
	#cluster3="Land Plant II"
	#Brassicaceae="Brasicaceae Plants"
	#monocluster="Poaceae Plants"
	#Else="Plants with Stochastic retaining/duplicated UBLs"
	
# Functions

# Distribution of a specific UBL in 6 OrthoMCL clusters

	ubl_cluster_distribution<-function(ubl,d){

		clusters<-c("cluster1","cluster2","cluster3","Brassicaceae","monocluster","Else")

		ubl_cluster_orthoMCL_count<-c()
		
		for (i in 1:6){

			name<-clusters[i]
			clusteri<-d[(d$Cluster==name),]

			clusteri_ubl<-clusteri[(clusteri$UBL_DM==ubl),]

			orthoi_ubl<-table(clusteri_ubl$OrthoMCL)

			orthoi_ubl<-orthoi_ubl[orthoi_ubl>0]
	
			ubl_ortho_count<-length(orthoi_ubl)
		
			ubl_cluster_orthoMCL_count<-cbind(ubl_cluster_orthoMCL_count,ubl_ortho_count)
		
				}

		
		ubl_cluster_orthoMCL_count<-as.matrix(ubl_cluster_orthoMCL_count)

		colnames(ubl_cluster_orthoMCL_count)<-clusters
		
		ubl_cluster_orthoMCL_count
		
		}
		
#Calculate the OrthoMCL groups that contain a specific UBL domain

	OrthoMCL_Clusters_w_UBL<-function(ubl,d){

		ubls<-d[(d$UBL_DM==ubl),]
		ubl_orthos<-table(ubls$OrthoMCL)
		ubl_orthos<-ubl_orthos[ubl_orthos>0]
		ubl_orthos<-ubl_orthos[order(names(ubl_orthos))]

		ubl_orthos

			}

#The proportion of 4 UBLs in each cluster

	ubl_proportions<-function(x){

		ubl<-d[(d$UBL_DM==x),]	
		ubl_proportion<-table(ubl$Cluster)/dim(ubl)[1]	
		ubl_proportion

			}

# ubl_count_by_intron
	
	ubl_count_by_introns<-function(x){

		x_count_by_introns<-c()
		x<-x[(substring(x$ID,4,4)=="_"),]
	
		for(i in 0:3){

			if(i<3){
				xi<-x[(x$Introns==i),]			
					}
			else{
				xi<-x[(x$Introns>=i),]
				}
		
				xi_count<-dim(xi)[1]				
				x_count_by_introns<-cbind(x_count_by_introns,xi_count)
	 		
				}
			
		x_count_by_introns	

		} 
			

#fisher exact test

	fisher_exact_test<-function(x){

		fisher_exact_test_pvalues<-c()

		for(i in 2:dim(x)[1]){

			ubl<-rownames(x)[i]

			for(j in 1:dim(x)[2]){

				intron<-colnames(x)[j]

				all_count_w_assigned_intron<-x[1,j]
				all_count_wo_assigned_intron<-rowSums(x)[1]-all_count_w_assigned_intron

				ubl_count_w_assigned_intron<-x[i,j]
				ubl_count_wo_assigned_intron<-rowSums(x)[i]-ubl_count_w_assigned_intron

				test<-matrix(c(ubl_count_w_assigned_intron,all_count_w_assigned_intron,
								ubl_count_wo_assigned_intron,all_count_wo_assigned_intron), 
								nr=2, dimnames=list( c("a","b"),c("w","wo") ) )

				fe_greater_test<-fisher.test(test,alternative="greater")
				fe_less_test<-fisher.test(test,alternative="less")					
				fe_greater_test_pvalue<-c(paste(ubl,intron,"greater"),fe_greater_test$p.value)
				fe_less_test_pvalue<-c(paste(ubl,intron,"less"),fe_less_test$p.value)

				fisher_exact_test_pvalues<-rbind(fisher_exact_test_pvalues,fe_greater_test_pvalue)	
				fisher_exact_test_pvalues<-rbind(fisher_exact_test_pvalues,fe_less_test_pvalue)

								}

					}
			rownames(fisher_exact_test_pvalues)<-fisher_exact_test_pvalues[,1]

			pvalue<-p.adjust(fisher_exact_test_pvalues[,2], method = "bonferroni")			
			pvalue<-pvalue[(pvalue<0.05)]
			pvalue<-as.matrix(pvalue)

			pvalue			

		}


##########			


# Figure 3A: Distribution of the OrthoMCL groups from each subfamily in 6 major plant groups 

	
	atg8_orthos_in_clusters<-ubl_cluster_distribution("Atg8",d)
	mub_orthos_in_clusters<-ubl_cluster_distribution("Rad60-SLD2",d)
	sumo_orthos_in_clusters<-ubl_cluster_distribution("Rad60-SLD",d)
	ub_orthos_in_clusters<-ubl_cluster_distribution("ubiquitin",d)

	ubl_orthos_in_clusters<-rbind(atg8_orthos_in_clusters,mub_orthos_in_clusters,sumo_orthos_in_clusters,ub_orthos_in_clusters)

	rownames(ubl_orthos_in_clusters)<-c("Atg8","MUB","SUMO","Ub")

	x<-ubl_orthos_in_clusters

	rowDistance = dist(x, method = "manhattan")
	rowCluster = hclust(rowDistance, method = "ward.D")
	rowDend = as.dendrogram(rowCluster)
	rowDend = reorder(rowDend, rowMeans(x))

	colDistance = dist(t(x), method = "manhattan")
	colCluster = hclust(colDistance, method = "ward.D")
	colDend = as.dendrogram(colCluster)
	colDend = reorder(colDend, colMeans(x))

	my_palette<-colorRampPalette(c("light blue","blue","red"))(n=39)

	col_breaks=c(seq(1,5,length=10),
			seq(5.1,20,length=10),
			seq(20.1,30,,length=20))


	pdf("Figure 3A.pdf",family="Times",height=10,width=10)

	par(mar=c(2.75,2.75,0.5,.5),mgp=c(1.7,.7,0))

	heatmap.2(x, Rowv=rowDend , Colv=FALSE,
 		col=my_palette, 
		breaks=col_breaks,
		keysize=2, 
		#RowSideColors=ROWcolors, 
		#ColSideColors=colcolors,
		trace=("none"), margins=c(3,5), 
		dendrogram=c("row"),
		#labRow="",
		#labCol="",
		density.info=c("none")

			)

	dev.off()


#Figure 3B. Number distribution of each subfamily in 6 major plant groups 


#1 Atg8

	OrthoMCL_Clusters_w_UBL("Atg8",d)

#2 MUB

	OrthoMCL_Clusters_w_UBL("Rad60-SLD2",d)

#3 SUMO

	OrthoMCL_Clusters_w_UBL("Rad60-SLD",d)

#4 Ubiquitin

	ub_orthos<-OrthoMCL_Clusters_w_UBL("ubiquitin",d)

	length(ub_orthos) # in total 264 OrthoMCLs containing Ubiquitin domains.

	total_OrthoMCLs<-table(d[,1])

	length(total_OrthoMCLs) # total OrthoMCL groups: 375


	
	Atg8_proportion<-ubl_proportions("Atg8")
	MUB_proportion<-ubl_proportions("Rad60-SLD2")
	SUMO_proportion<-ubl_proportions("Rad60-SLD")
	Ub_proportion<-ubl_proportions("ubiquitin")


	ub<-d[(d$UBL_DM=="ubiquitin"),]

	ubi_proportions<-function(x){


		if(x<5){	
			ubi<-ub[(ub$Moiety==x),]		
			}		
		else{	
			ubi<-ub[(ub$Moiety>=x),]
			}	
		
		ubi_proportion<-table(ubi$Cluster)/dim(ubi)[1]	
		ubi_proportion
	
		}
	ub1_proportion<-ubi_proportions(1)
	ub2_proportion<-ubi_proportions(2)
	ub3_proportion<-ubi_proportions(3)
	ub4_proportion<-ubi_proportions(4)
	ub5_proportion<-ubi_proportions(5)


	ubl_ubi_clusters<-rbind(Atg8_proportion,MUB_proportion,SUMO_proportion,ub1_proportion,
							ub2_proportion,ub3_proportion,ub4_proportion,ub5_proportion)							
		
	ubl_ub_clusters<-cbind(ubl_ubi_clusters[,2],ubl_ubi_clusters[,3],
							ubl_ubi_clusters[,4],ubl_ubi_clusters[,1],
							ubl_ubi_clusters[,6],ubl_ubi_clusters[,5])

	colnames(ubl_ub_clusters)<-c("Cluster1","Cluster2","Cluster3","Brassicaceae","Poaceae","Stochastic")


	x<-ubl_ub_clusters

	my_palette<-colorRampPalette(c("light blue","yellow","red"))(n=29)
	col_breaks=c(seq(0,0.1,length=10),
			seq(0.101,0.6,length=10),
			seq(0.601,1,length=10))

	rowDistance = dist(x, method = "manhattan")
	rowCluster = hclust(rowDistance, method = "ward.D")
	rowDend = as.dendrogram(rowCluster)
	rowDend = reorder(rowDend, rowMeans(x))
	
pdf("Figure 3B.pdf",family="Times",height=10,width=10)

			heatmap.2(x, Rowv=rowDend, Colv=FALSE,
			 		col=my_palette, 
					breaks=col_breaks,
					keysize=2, 
					#RowSideColors=ROWcolors, 
					#ColSideColors=colcolors,
					trace=("none"), margins=c(3,5), 
					dendrogram=c("row"),
					#labRow="",
					#labCol="",
					density.info=c("none")

						)

	dev.off()


#Figure 3C. Size expansion correlation analysis 

	x<-cor(t(ubl_ub_clusters))

	my_palette<-colorRampPalette(c("blue","yellow","red"))(n=29)
	col_breaks=c(seq(-1,0,length=10),
			seq(0.01,0.8,length=10),
			seq(0.81,1,length=10))

	rowDistance = dist(x, method = "manhattan")
	rowCluster = hclust(rowDistance, method = "ward.D")
	rowDend = as.dendrogram(rowCluster)
	rowDend = reorder(rowDend, rowMeans(x))

	colDistance = dist(t(x), method = "manhattan")
	colCluster = hclust(colDistance, method = "ward.D")
	colDend = as.dendrogram(colCluster)
	colDend = reorder(colDend, colMeans(x))
	


pdf("Figure 3C.pdf",family="Times",height=10,width=10)

	
	par(mar=c(2.75,2.75,0.5,0.5),mgp=c(1.7,.7,0))
	#layout( matrix(c(1,2,3,4,5,6),nrow=2,ncol=3),widths=c(1,1,1))

	heatmap.2(x, Rowv=rowDend, Colv=colDend,

 		col=my_palette, 
		breaks=col_breaks,
		keysize=2, 
		#RowSideColors=ROWcolors, 
		#ColSideColors=colcolors,
		trace=("none"), margins=c(3,5), 
		#dendrogram=c("row"),
		#labRow="",
		#labCol="",
		density.info=c("none")

			)

	dev.off()




################

#Figure S12. Fraction of poly-ubiquiton genes 

#Atg8

	atg8<-d[(d$UBL_DM=="Atg8"),]
	multi_atg8<-atg8[(atg8$Moiety>1),]
	multi_atg8_ratio<-dim(multi_atg8)[1]/dim(atg8)[1]*100
	
	# [1] 4.632153
	

#MUB
	mub<-d[(d$UBL_DM=="Rad60-SLD2"),]
	multi_mub<-mub[(mub$Moiety>1),]
	multi_mub_ratio<-dim(multi_mub)[1]/dim(mub)[1]*100	
	
	# [1] 0

#SUMO
	sumo<-d[(d$UBL_DM=="Rad60-SLD"),]
	multi_sumo<-sumo[(sumo$Moiety>1),]
	multi_sumo_ratio<-dim(multi_sumo)[1]/dim(sumo)[1]*100	
	
	# [1] 7.541899

	
	ub<-d[(d$UBL_DM=="ubiquitin"),]
	multi_ub<-ub[(ub$Moiety>1),]
	multi_ub_ratio<-dim(multi_ub)[1]/dim(ub)[1]*100

	# [1] 25.23478

	
	multi_ubl_ratios<-cbind(multi_atg8_ratio,multi_mub_ratio,multi_sumo_ratio,multi_ub_ratio)

pdf("Figure S12.pdf",family="Times",height=10,width=10)
	
	barplot(multi_ubl_ratios, main="multi_ubl_ratios", ylab="Proportion of Group", ylim=c(0,30))
	
	dev.off()



########################################

#Figure S11 An anti-correlation relationship between the number of introns and the moieties of Ub coding regions in Ub genes.

	ub<-d[(d$UBL_DM=="ubiquitin"),]
	
	ub<-ub[(substring(ub[,3],4,4)=="_"),]

	dim(ub)[1]
	
		# [1] 3035

	cor.test(ub$Introns,ub$Moiety,method = "spearman",exact = FALSE, co.level = 0.95)
	
		# Spearman's rank correlation rho

		# data:  ub$Introns and ub$Moiety
		# S = 6635100000, p-value < 2.2e-16
		# alternative hypothesis: true rho is not equal to 0
		# sample estimates:
		#        rho 
		# -0.4240432 

 	plot(ub$Introns,ub$Moiety,cex=2, pch=19, col="blue", ylim=c(0,15), xlim=c(0,60), xlab="Number of Introns", ylab="Number of Ub Moieties")

	text(30,14,"rho=-0.42,p-value < 2.2e-16")

pdf("Figure S11.pdf",family="Times",height=5,width=5)

	plot(ub$Introns,ub$Moiety,cex=2, pch=19, col="blue", ylim=c(0,15), xlim=c(0,60), xlab="Number of Introns", ylab="Number of Ub Moieties")
	text(30,14,"rho=-0.42,p-value < 2.2e-16")

	dev.off()



################



####################

#Figure 6A, The fraction of ubiquiton genes from the ATG8, MUB, SUMO, and Ub subfamilies with 0, 1, 2, or ≥3 introns

	atg8<-d[(d$UBL_DM=="Atg8"),]
	mub<-d[(d$UBL_DM=="Rad60-SLD2"),]
	sumo<-d[(d$UBL_DM=="Rad60-SLD"),]
	ub<-d[(d$UBL_DM=="ubiquitin"),]

	atg8_count_by_introns<-ubl_count_by_introns(atg8)
	mub_count_by_introns<-ubl_count_by_introns(mub)
	sumo_count_by_introns<-ubl_count_by_introns(sumo)
	ub_count_by_introns<-ubl_count_by_introns(ub)
	all_count_by_intron<-ubl_count_by_introns(d)
	

	ubls_count_by_introns<-rbind(all_count_by_intron,atg8_count_by_introns,
								mub_count_by_introns,sumo_count_by_introns,ub_count_by_introns)

	colnames(ubls_count_by_introns)<-c("Intron0","Intron1","Intron2","Intron3_or_greater")

	rownames(ubls_count_by_introns)<-c("All","Atg8","MUB","SUMO","Ubiquitin")

	ubls_count_by_introns
		
			#          Intron0 Intron1 Intron2 Intron3_or_greater
			# All           784     610     693               2463
			# Atg8            1       8       7                322
			# MUB             6      10     217                 35
			# SUMO           48      19     153                106
			# Ubiquitin     705     564     292               1474
		 	

	ubls_ratio_by_introns<-ubls_count_by_introns/rowSums(ubls_count_by_introns)*100		
	ubls_ratio_by_introns
	
			# 	      Intron0   Intron1   Intron2 Intron3_or_greater
			# All       17.230769 13.406593 15.230769           54.13187
			# Atg8       0.295858  2.366864  2.071006           95.26627
			# MUB        2.238806  3.731343 80.970149           13.05970
			# SUMO      14.723926  5.828221 46.932515           32.51534
			# Ubiquitin 23.228995 18.583196  9.621087           48.56672

	ratios<-c()
	for(i in 1:4){

		axis<-cbind(0,colnames(ubls_ratio_by_introns)[i])
		ratios<-rbind(ratios,axis)

		for(j in 1:5){

			ubl_ratio<-cbind(ubls_ratio_by_introns[j,i],colnames(ubls_ratio_by_introns)[i])

			ratios<-rbind(ratios,ubl_ratio)

			}

		}

				
	t<-c("0","15","30","45","60","75","90",
			"105","120","135","150","165","180",
			"195","210","225","240","255","270",
			"285","300","315","330","345")

	r<-ratios[,1]
	nms<-ratios[,2]		
	ratios<-cbind(r,t,nms)	
	ratios<-data.frame(ratios)

	p <- plot_ly(ratios, r = ~r, t = ~t, sizes=c(0,100)) %>% add_area
	layout(p, radialaxis = list(ticksuffix = "%"), orientation = 270)


	fisher_exact_test(ubls_count_by_introns)
			
			#                                            [,1]
			# Atg8 Intron0 less                  4.377439e-24
			# Atg8 Intron1 less                  1.262927e-10
			# Atg8 Intron2 less                  5.671267e-14
			# Atg8 Intron3_or_greater greater    5.714132e-60
			# MUB Intron0 less                   4.587523e-13
			# MUB Intron1 less                   4.851581e-06
			# MUB Intron2 greater               1.586929e-114
			# MUB Intron3_or_greater less        3.441762e-41
			# SUMO Intron1 less                  4.179138e-04
			# SUMO Intron2 greater               1.779892e-36
			# SUMO Intron3_or_greater less       7.571714e-13
			# Ubiquitin Intron0 greater          2.902429e-09
			# Ubiquitin Intron1 greater          2.575939e-08
			# Ubiquitin Intron2 less             9.728555e-12
			# Ubiquitin Intron3_or_greater less  3.607330e-05


##########

##Figure 6b, Differential intron-richness among the Ub subfamilies encoding 1, 2, 3, 4, or ≥5 Ub domains

	
	ub<-d[(d$UBL_DM=="ubiquitin"),]
	ub1<-ub[(ub$Moiety==1),]
	ub2<-ub[(ub$Moiety==2),]
	ub3<-ub[(ub$Moiety==3),]
	ub4<-ub[(ub$Moiety==4),]
	ub5<-ub[(ub$Moiety>=5),]

	
	ub1_count_by_introns<-ubl_count_by_introns(ub1)
	ub2_count_by_introns<-ubl_count_by_introns(ub2)
	ub3_count_by_introns<-ubl_count_by_introns(ub3)
	ub4_count_by_introns<-ubl_count_by_introns(ub4)
	ub5_count_by_introns<-ubl_count_by_introns(ub5)
	ub_count_by_introns<-ubl_count_by_introns(ub)

	ubis_count_by_introns<-rbind(ub_count_by_introns,
								ub1_count_by_introns,
								ub2_count_by_introns,
								ub3_count_by_introns,
								ub4_count_by_introns,
								ub5_count_by_introns
								)

	colnames(ubis_count_by_introns)<-c("Intron0","Intron1","Intron2","Intron3_or_greater")

	rownames(ubis_count_by_introns)<-c("ub","ub1","ub2","ub3","ub4","ub5")

	ubis_count_by_introns

			
			#     Intron0 Intron1 Intron2 Intron3_or_greater
			# ub      705     564     292               1474
			# ub1     423     276     177               1427
			# ub2      77     140      89                 30
			# ub3      91      34       8                  3
			# ub4      34      40       2                  4
			# ub5      80      74      16                 10

	
	ubis_ratio_by_introns<-ubis_count_by_introns/rowSums(ubis_count_by_introns)*100

	ubis_ratio_by_introns


			#      Intron0  Intron1   Intron2 Intron3_or_greater
			# ub  23.22900 18.58320  9.621087          48.566722
			# ub1 18.36735 11.98437  7.685627          61.962657
			# ub2 22.91667 41.66667 26.488095           8.928571
			# ub3 66.91176 25.00000  5.882353           2.205882
			# ub4 42.50000 50.00000  2.500000           5.000000
			# ub5 44.44444 41.11111  8.888889           5.555556

	ubi_ratios<-c()
	for(i in 1:4){

			axis<-cbind(0,colnames(ubis_ratio_by_introns)[i])
			ubi_ratios<-rbind(ubi_ratios,axis)

			for(j in 1:6){

					ubi_ratio<-cbind(ubis_ratio_by_introns[j,i],colnames(ubis_ratio_by_introns)[i])

					ubi_ratios<-rbind(ubi_ratios,ubi_ratio)

				}

				}

	t<-c("0","13","26","39","52","65","78","90",
			"103","116","129","142","155","168","180",
			"193","206","219","232","245","258","270",
			"283","296","309","322","335","348")

	r<-ubi_ratios[,1]
	nms<-ubi_ratios[,2]		
	ubi_ratios<-cbind(r,t,nms)	
	ubi_ratios<-data.frame(ubi_ratios)

	p <- plot_ly(ubi_ratios, r = ~r, t = ~t, sizes=c(0,100)) %>% add_area
	layout(p, radialaxis = list(ticksuffix = "%"), orientation = 270)


	fisher_exact_test(ubis_count_by_introns)

			#                                        [,1]
			# ub1 Intron0 less               3.498630e-04
			# ub1 Intron1 less               8.743252e-10
			# ub1 Intron3_or_greater greater 4.483719e-21
			# ub2 Intron1 greater            1.675837e-18
			# ub2 Intron2 greater            5.971474e-15
			# ub2 Intron3_or_greater less    1.877552e-49
			# ub3 Intron0 greater            3.313844e-24
			# ub3 Intron3_or_greater less    9.535146e-32
			# ub4 Intron0 greater            5.285553e-03
			# ub4 Intron1 greater            1.718324e-08
			# ub4 Intron3_or_greater less    9.325882e-16
			# ub5 Intron0 greater            4.445118e-08
			# ub5 Intron1 greater            4.882023e-10
			# ub5 Intron3_or_greater less    8.339623e-34
			


###################

##Figure 6C_D. Comparison of Ks values of Ub genes with different number of Ub moieties and introns

#boxplot

#ubiquitin #ks kaks

		ks_df<-c()
		kaks_df<-c()
		ks_pvalues<-c()
		kaks_pvalues<-c()		

		ub<-d[(d$UBL_DM=="ubiquitin"),]
		ub<-ub[(substring(ub$ID,4,4)=="_"),]		
		ub<-ub[(ub$Ks<=5),]
		
		ub_ks_df<-data.frame(Measure1=ub$Ks,Group1="ub")
		
		ub_kaks_df<-data.frame(Measure2=ub$Ka_Ks,Group2="ub")
		
		ks_df<-rbind(ks_df,ub_ks_df)
		kaks_df<-rbind(kaks_df,ub_kaks_df)
		
		for(j in 1:5){

			if(j<5){
				ubj<-ub[(ub$Moiety==j),]
							}
			else{
				ubj<-ub[(ub$Moiety>=j),]	
						}
		
			ubij_count<-c()
					
				introns<-c("intron0","intron1","intron2","intron3");			
		
			for(i in 0:3){
					
				name<-introns[i+1]
				if(i<3){
					ubij<-ubj[(ubj$Introns==i),]
					}
				else{
					ubij<-ubj[(ubj$Introns>=i),]		
					}
					
				size<-dim(ubij)[1]
				
				ubij_ks_df<-data.frame(Measure1=ubij$Ks, Group1=paste(j,name,size,sep="_"))				
				ubij_kaks_df<-data.frame(Measure2=ubij$Ka_Ks, Group2=paste(j,name,size,sep="_"))
								
				ks_df<-rbind(ks_df,ubij_ks_df)
				kaks_df<-rbind(kaks_df,ubij_kaks_df)	
				
				wilcox_ks_greater_test<-wilcox.test(ubij$Ks,ub$Ks,alternative="greater")			
				wilcox_ks_less_test<-wilcox.test(ubij$Ks,ub$Ks,alternative="less")
			
				ks_greater_pvalue<-c(paste("Ub",j,name,"greater",sep="_"),wilcox_ks_greater_test$p.value)
				ks_less_pvalue<-c(paste("Ub",j,name,"less",sep="_"),wilcox_ks_less_test$p.value)								
								
				wilcox_kaks_greater_test<-wilcox.test(ubij$Ka_Ks,ub$Ka_Ks,alternative="greater")								
				wilcox_kaks_less_test<-wilcox.test(ubij$Ka_Ks,ub$Ka_Ks,alternative="less")
			
				kaks_greater_pvalue<-c(paste("Ub",j,name,"greater",sep="_"),wilcox_kaks_greater_test$p.value)
				kaks_less_pvalue<-c(paste("Ub",j,name,"less",sep="_"),wilcox_kaks_less_test$p.value)

				ks_pvalues<-rbind(ks_pvalues,ks_greater_pvalue,ks_less_pvalue)
				kaks_pvalues<-rbind(kaks_pvalues,kaks_greater_pvalue,kaks_less_pvalue)
							
					}
					
				}


	rownames(ks_pvalues)<-ks_pvalues[,1]

	ks_pvalue<-p.adjust(ks_pvalues[,2], method = "bonferroni")			
	ks_pvalue<-ks_pvalue[(ks_pvalue<0.05)]
	as.matrix(ks_pvalue)

			#                              [,1]
			# Ub_4_intron0_greater 5.465363e-04
			# Ub_5_intron0_greater 4.932267e-15
			# Ub_5_intron1_greater 8.982147e-11
			# Ub_5_intron2_greater 2.006791e-02


	rownames(kaks_pvalues)<-kaks_pvalues[,1]

	kaks_pvalue<-p.adjust(kaks_pvalues[,2], method = "bonferroni")			
	kaks_pvalue<-kaks_pvalue[(kaks_pvalue<0.05)]
	as.matrix(kaks_pvalue)

			#                              [,1]
			# Ub_1_intron2_greater 3.545532e-02
			# Ub_2_intron2_less    2.274611e-18

	attach(ks_df)
	attach(kaks_df)

pdf("Figure 6C_D.pdf", family="Times",height=10, width=10)

	par(mar=c(5.1,4.1,4.1,2.1), mgp=c(3,1,0), las=0)
	layout(matrix(c(1,2), nrow=2,ncol=1))

	boxplot(Measure1 ~ Group1,ylim=c(0,5), cex.names=0.3,main="ubi_intron_ks_comparison",col=c("lightblue","gray","red","yellow","green"),outline=FALSE)
	boxplot(Measure2 ~ Group2,ylim=c(0,1.5),cex.names=0.3,main="ubi_intron_kaks_comparison",col=c("lightblue","gray","red","yellow","green"),outline=FALSE)

	dev.off()
   



####
