
#####################

#To run the script, save Data S4 "Dataset4_for_Figures_5_S6_S9_S10.tab" and the file of this script,"RScript4_for_Figures_5_S6_S9_S10.R", in the same directory and type "Rscript RScript4_for_Figures_5_S6_S9_S10.R" in terminal under the same directory to generate all Figures.

#####################

d<-read.table("Dataset4_for_Figures_5_S6_S9_S10.tab",header=T);

d<-d[(d$Ks<=5),] #only genes with Ks <=5 are analyzed to avoid sequences that are too divergent.

#required packages

	library(gplots)

#count the proportion of WGDs in a certain group of UBLs

# Nomenclature 

		# x: group names
		# y: groups
		# z: WGD events
		# wgd: duplication dataset

#Functions

#Calculate the distribution of UBLs in each group
		
	percent_wgd_events<-function(wgd,x,y,z){

		clusteri_wgd_types<-c()
	
		for (i in 1:length(x)){

			clusteri_wgd<-wgd[(wgd[,as.numeric(y)] == x[i]), ]

			clusteri_wgd_type<-(table(clusteri_wgd[,as.numeric(z)])/dim(clusteri_wgd)[1])*100

			clusteri_wgd_types<-rbind(clusteri_wgd_types,clusteri_wgd_type)

			}

		cluster_total_wgd_type<-(table(wgd[,as.numeric(z)])/dim(wgd)[1])*100

		clusteri_wgd_types<-rbind(cluster_total_wgd_type,clusteri_wgd_types)

		clusters<-c("All",x)
		rownames(clusteri_wgd_types)<-clusters
		clusteri_wgd_types
	
	}
	
#Count the number of UBLs in each group

	count_wgd_events<-function(wgd,x,y,z){

		clusteri_wgd_types<-c()

		for (i in 1:length(x)){

				clusteri_wgd<-wgd[(wgd[,as.numeric(y)] == x[i]), ]

				clusteri_wgd_type<-table(clusteri_wgd[,as.numeric(z)])

				clusteri_wgd_types<-rbind(clusteri_wgd_types,clusteri_wgd_type)

				}

		cluster_total_wgd_type<-table(wgd[,as.numeric(z)])

		clusteri_wgd_types<-rbind(cluster_total_wgd_type,clusteri_wgd_types)

		clusters<-c("All",x)
		rownames(clusteri_wgd_types)<-clusters
		clusteri_wgd_types

		}

#Fisher exact test with bonferroni pvalue 0.05 cutoff 

	fisher_test<-function(wgd,x,y,z){

		fisher_exact_test<-c()
	
		wgd_events<-c("Cretaceous","paleo","neo","tandem","unknown","na")
	
		for(i in 1:length(wgd_events)){
	
			wgd_event<-wgd_events[i]
		
			all_wgd_event<-wgd[wgd$WGD==wgd_event,]
			all_wgd_event_count<-dim(all_wgd_event)[1]
			non_all_wgd_event_count<-dim(wgd)[1]-dim(all_wgd_event)[1]
		
			for(j in 1:length(x)){
		
				name<-x[j]
				clusteri<-wgd[(wgd[,as.numeric(y)] == name), ]			
				clusteri_wgd_event<-clusteri[clusteri$WGD==wgd_event,]
				clusteri_wgd_event_count<-dim(clusteri_wgd_event)[1]
				non_clusteri_wgd_event_count<-dim(clusteri)[1]-dim(clusteri_wgd_event)[1]
			
				test<-matrix(c(clusteri_wgd_event_count,all_wgd_event_count,
							non_clusteri_wgd_event_count,non_all_wgd_event_count), 
							nr=2, dimnames=list( c("a","b"),c("w","wo") ) )
			
				fe_greater_test<-fisher.test(test,alternative="greater")
				fe_less_test<-fisher.test(test,alternative="less")
			
				fe_greater_test<-c(paste(name,wgd_event,"greater"),fe_greater_test$p.value)
				fe_less_test<-c(paste(name,wgd_event,"less"),fe_less_test$p.value)
			
				fisher_exact_test<-rbind(fisher_exact_test,fe_greater_test)	
				fisher_exact_test<-rbind(fisher_exact_test,fe_less_test)
			
							}
					}
	
			
			rownames(fisher_exact_test)<-fisher_exact_test[,1]

			pvalue<-p.adjust(fisher_exact_test[,2], method = "bonferroni")			
			pvalue<-pvalue[(pvalue<0.05)]
			pvalue<-as.matrix(pvalue)

			pvalue	
		
			}


#make a ks or ka/ks data frame for boxplot
	
	k_df<-function(wgd,x,y,z){
	
		clusteri_wgd_k<-c()
	
		for (i in 1:length(x)){
		
			name<-x[i]
			clusteri_wgd<-wgd[(wgd[,as.numeric(y)] == name), ]
			clusteri_wgd_k_df<-data.frame(Measure=clusteri_wgd[,as.numeric(z)],Group=name)		
			clusteri_wgd_k<-rbind(clusteri_wgd_k,clusteri_wgd_k_df)

			}	
	
	
		cluster_total_wgd_k_df<-data.frame(Measure=wgd[,as.numeric(z)],Group="All")	

		clusteri_wgd_k<-rbind(cluster_total_wgd_k_df,clusteri_wgd_k)
		
		clusteri_wgd_k
	
		}


#wilcox test to compare the Ks or Ka/Ks values with bonferroni pvalue 0.01 cutoff 

	wilcoxtest<-function(wgd,x,y,z){

		wilcox_test<-c()
		
		total_wgd_k<-wgd[,as.numeric(z)]	

		for (i in 1:length(x)){
	
			name<-x[i]
			clusteri_wgd<-wgd[(wgd[,as.numeric(y)] == name), ]
			clusteri_wgd_k<-clusteri_wgd[,as.numeric(z)]
				
			w_great_k_test<-wilcox.test(clusteri_wgd_k,total_wgd_k, alternative="greater" )
			w_greater_test<-c(paste(name,"greater"),w_great_k_test$p.value)
			wilcox_test<-rbind(wilcox_test,w_greater_test)	
							
			w_less_k_test<-wilcox.test(clusteri_wgd_k,total_wgd_k,alternative="less" )
			w_less_test<-c(paste(name,"less"),w_less_k_test$p.value)
			wilcox_test<-rbind(wilcox_test,w_less_test)	
		
			}
			
		rownames(wilcox_test)<-wilcox_test[,1]
		
		pvalue<-p.adjust(wilcox_test[,2], method = "bonferroni")			
		pvalue<-pvalue[(pvalue<0.05)]
		pvalue<-as.matrix(pvalue)

		pvalue
			
		}


#Figure 5A. Contribution of differential duplication events to the expansion of the ubiquiton genes at 6 different conservational levels as defined in Figure 2.

	#proportion of WGDs in each OrthoMCL cluster

	x<-c("cluster1","cluster2","cluster3","Brassicaceae","monocluster","Else")	
	y<-2
	z<-5
	wgd<-d
	
	percent_wgd_events(wgd,x,y,z)
	
		##percent_wgd_events result
		#        	     Cretaceous        na       neo    paleo    tandem  unknown
		# All            37.77030 0.8649688 11.821240 28.44786  7.640557 13.45507
		# cluster1       38.24885 0.6144393 13.671275 28.72504  5.529954 13.21045
		# cluster2       41.74041 0.7374631 13.421829 30.53097  2.212389 11.35693
		# cluster3       46.68508 1.9337017  8.839779 27.62431  1.657459 13.25967
		# Brassicaceae   18.64407 0.8474576  7.627119 42.37288 20.338983 10.16949
		# monocluster    15.71429 0.0000000  7.142857 34.28571 25.714286 17.14286
		# Else           25.74257 0.4950495  9.900990 11.88119 29.702970 22.27723
		##

#Figure S9A. Co-evolution analysis of ubiquiton genes at different conservational levels 

	percent<-percent_wgd_events(wgd,x,y,z)	
	cor<-cor(t(percent))
	
	my_palette<-colorRampPalette(c("blue","yellow","red"))(n=29)
	col_breaks=c(seq(-1,0,length=10),
				seq(0.01,0.8,length=10),
				seq(0.81,1,length=10))

	rowDistance = dist(cor, method = "manhattan")
	rowCluster = hclust(rowDistance, method = "ward.D")
	rowDend = as.dendrogram(rowCluster)
	rowDend = reorder(rowDend, rowMeans(cor))


	colDistance = dist(t(cor), method = "manhattan")
	colCluster = hclust(colDistance, method = "ward.D")
	colDend = as.dendrogram(colCluster)
	colDend = reorder(colDend, colMeans(cor))

	pdf("Figure S9A.pdf",family="Times",height=10,width=10)

		par(mar=c(2.75,2.75,0.5,0.5),mgp=c(1.7,.7,0))
		
	
		heatmap.2(cor, Rowv=rowDend, Colv=colDend,
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

	count_wgd_events(wgd,x,y,z)
		
		#             Cretaceous na neo paleo tandem unknown
		# All                 786 18 246   592    159     280
		# cluster1            249  4  89   187     36      86
		# cluster2            283  5  91   207     15      77
		# cluster3            169  7  32   100      6      48
		# Brassicaceae         22  1   9    50     24      12
		# monocluster          11  0   5    24     18      12
		# Else                 52  1  20    24     60      45


	fisher_test(wgd,x,y,z) #with Bonferroni correction

		#                                   [,1]
		# Brassicaceae Cretaceous less 6.626763e-04
		# monocluster Cretaceous less  4.400828e-03
		# Else Cretaceous less         2.593657e-02
		# Else paleo less              3.238501e-06
		# cluster2 tandem less         1.810314e-06
		# cluster3 tandem less         1.151059e-04
		# Brassicaceae tandem greater  1.384592e-03
		# monocluster tandem greater   4.770397e-04
		# Else tandem greater          5.572548e-16
		# 
		
#Figure 5B. Comparison of Ks and Ka/Ks values for the groups in Figure 5A.


	#Ks comparison

	x<-c("cluster1","cluster2","cluster3","Brassicaceae","monocluster","Else")	
	y<-2	
	z<-6
	wgd<-d
	
	wilcoxtest(wgd,x,y,z)

		##wilcoxtest result
	    #                     [,1]
		# cluster1 greater    3.948991e-06
		# Brassicaceae less   5.256429e-09
		# monocluster greater 4.630509e-04
		# Else less           2.317387e-03
		##	
		
	cluster_ks_df<-c()
	cluster_ks_df<-k_df(wgd,x,y,z)
	
	attach(cluster_ks_df)
		
	pdf("Figure 5B_Ks.pdf", family="Times",height=7, width=2.5)

	boxplot(cluster_ks_df$Measure ~ cluster_ks_df$Group, ylim=c(0,1.5), 	 
	cex.names=0.3,main="clusteri_wgd_ks_comparison",col=c("red","yellow"),outline=FALSE)
	
	dev.off()
	
	
#Ka/Ks comparison

	x<-c("cluster1","cluster2","cluster3","Brassicaceae","monocluster","Else")	
	y<-2	
	z<-7
	wgd<-d

	wilcoxtest(wgd,x,y,z)

		##wilcoxtest result
		#                      [,1]
		# cluster1 less        6.244881e-18
		# Brassicaceae greater 1.825416e-10
		# Else greater         6.045595e-26
		##

	cluster_kaks_df<-c()
	cluster_kaks_df<-k_df(wgd,x,y,z)	
	attach(cluster_kaks_df)
	
	pdf("Figure 5B_KaKs.pdf", family="Times",height=7, width=2.5)

	boxplot(cluster_kaks_df$Measure ~ cluster_kaks_df$Group, ylim=c(0,1.5), 
	cex.names=0.3,main="clusteri_wgd_kaks_comparison",col=c("red","yellow"),outline=FALSE)
	
	dev.off()
	
###########

#proportion of WGDs in ATG8(Atg8), MUB(Rad60-SLD2),SUMO(Rad60-SLD), and Ub(ubiquitin) subfamilies

#Figure 5C. Enrichment comparison of 4 ubiquiton subfamily members with the total members that were derived from different duplication mechanisms 	
	
	x<-c("Atg8","Rad60-SLD2","Rad60-SLD","ubiquitin")
	y<-4
	z<-5
	wgd<-d
	
	percent_wgd_events(wgd,x,y,z)
	
		##percent_wgd_events result

		#        	  Cretaceous        na       neo    paleo    tandem   unknown
		# All          37.77030 0.8649688 11.821240 28.44786  7.640557 13.455070
		# Atg8         36.12903 1.9354839 14.838710 36.77419  1.290323  9.032258
		# Rad60-SLD2   43.30709 0.0000000 19.685039 29.13386  1.574803  6.299213
		# Rad60-SLD    38.85350 0.6369427  8.917197 20.38217 19.108280 12.101911
		# ubiquitin    36.75400 0.9461426 11.280932 28.09316  9.097525 13.828239	
		##
		
	percent_ubls<-percent_wgd_events(wgd,x,y,z)
	
	count_wgd_events(wgd,x,y,z)
	
		#            Cretaceous na neo paleo tandem unknown
		# All               786 18 246   592    159     280
		# Atg8               56  3  23    57      2      14
		# Rad60-SLD2         55  0  25    37      2       8
		# Rad60-SLD          61  1  14    32     30      19
		# ubiquitin         505 13 155   386    125     190
		
	
	fisher_test(wgd,x,y,z)
		
		#                                 [,1]
		# Atg8 tandem less         0.0281768375
		# Rad60-SLD tandem greater 0.0003984467
		#	
		
#Figure 5D. Comparison of Ks and Ka/Ks values for the groups in Figure 5C.

#Ks comparison	

	x<-c("Atg8","Rad60-SLD2","Rad60-SLD","ubiquitin")
	y<-4
	z<-6
	wgd<-d
	
	wilcoxtest(wgd,x,y,z)

		##wilcoxtest result
		#                      [,1]
		# Rad60-SLD2 less 0.03453814
		##
		
	cluster_ks_df<-c()
	cluster_ks_df<-k_df(wgd,x,y,z)	
	attach(cluster_ks_df)	
	
	pdf("Figure 5D_Ks.pdf", family="Times",height=7, width=2.5)

	boxplot(cluster_ks_df$Measure ~ cluster_ks_df$Group, ylim=c(0,1.5), 
	cex.names=0.3,main="clusteri_wgd_ks_comparison",col=c("red","yellow"),outline=FALSE)
	
	dev.off()

#Ka/Ks comparison	

	x<-c("Atg8","Rad60-SLD2","Rad60-SLD","ubiquitin")	
	y<-4	
	z<-7
	wgd<-d
	
	wilcoxtest(wgd,x,y,z)
	
		##wilcoxtest result
		#                         [,1]
		# Atg8 less         5.856962e-19
		# ubiquitin greater 1.231931e-02
		#
		
	cluster_kaks_df<-c()
	cluster_kaks_df<-k_df(wgd,x,y,z)	
	attach(cluster_kaks_df)
	
	pdf("Figure 5D_KaKs.pdf", family="Times",height=7, width=2.5)

	boxplot(cluster_kaks_df$Measure ~ cluster_kaks_df$Group, ylim=c(0,1), 
	cex.names=0.3,main="clusteri_wgd_kaks_comparison",col=c("red","yellow"),outline=FALSE)
	
	dev.off()


#####
#proportion of WGDs in ub subfamilies containing different number of ubiquitin moieties

#Figure 5E. Enrichment comparison of 5 Ub subfamily members with the total members as described in (a) and (c) that were derived from different duplication mechanisms.  
	
	x<-c("ubiquitin_1","ubiquitin_2","ubiquitin_3","ubiquitin_4","ubiquitin_5")
	y<-3
	z<-5
	wgd<-d

	percent_wgd_events(wgd,x,y,z)
	
		#	percent_wgd_events result

		#             Cretaceous        na       neo     paleo    tandem  unknown
		# All           37.77030 0.8649688 11.821240 28.447862  7.640557 13.45507
		# ubiquitin_1   37.08609 0.8514664 12.109745 28.760643  8.420057 12.77200
		# ubiquitin_2   37.03704 2.4691358  9.259259 32.716049  5.555556 12.96296
		# ubiquitin_3   43.63636 0.0000000  7.272727 16.363636 14.545455 18.18182
		# ubiquitin_4   50.00000 0.0000000 13.333333  6.666667  6.666667 23.33333
		# ubiquitin_5   18.75000 0.0000000  3.125000 18.750000 21.875000 37.50000
	
	#distribution correlation test among ATG8, MUB, SUMO and 5 ubiquitin subfamilies

#Figure S9B. The ATG8, MUB, and SUMO subfamilies along with 5 Ub subfamilies that encode 1, 2, 3, 4, or 5 or more Ub-domains are grouped into 3 clusters that have experienced different duplication mechanisms.
		
	percent_ubi<-percent_wgd_events(wgd,x,y,z)

	percent_ubls<-percent_ubls[-5,]
	percent_ubls<-percent_ubls[-1,]
	percent_ubi<-percent_ubi[-1,]
			
	percent<-rbind(percent_ubls,percent_ubi)
			
	cor<-cor(t(percent))

	
	my_palette<-colorRampPalette(c("blue","yellow","red"))(n=29)
	col_breaks=c(seq(-1,0,length=10),
				seq(0.01,0.8,length=10),
				seq(0.81,1,length=10))



	rowDistance = dist(cor, method = "manhattan")
	rowCluster = hclust(rowDistance, method = "ward.D")
	rowDend = as.dendrogram(rowCluster)
	rowDend = reorder(rowDend, rowMeans(cor))


	colDistance = dist(t(cor), method = "manhattan")
	colCluster = hclust(colDistance, method = "ward.D")
	colDend = as.dendrogram(colCluster)
	colDend = reorder(colDend, colMeans(cor))


	pdf("Figure S9B.pdf",family="Times",height=10,width=10)

		par(mar=c(2.75,2.75,0.5,0.5),mgp=c(1.7,.7,0))
		
		heatmap.2(cor, Rowv=rowDend, Colv=colDend,
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
		
		
	cound_ubi<-count_wgd_events(wgd,x,y,z)
	
	cound_ubi

		
		#            Cretaceous na neo paleo tandem unknown
		# All                786 18 246   592    159     280
		# ubiquitin_1        392  9 128   304     89     135
		# ubiquitin_2         60  4  15    53      9      21
		# ubiquitin_3         24  0   4     9      8      10
		# ubiquitin_4         15  0   4     2      2       7
		# ubiquitin_5          6  0   1     6      7      12
		
	
	percent

		#             Cretaceous        na       neo    paleo    tandem   unknown
		# Atg8          36.53846 1.9230769 14.743590 36.53846  1.282051  8.974359
		# Rad60-SLD2    43.30709 0.0000000 19.685039 29.13386  1.574803  6.299213
		# Rad60-SLD     38.75000 0.6250000  8.750000 21.25000 18.750000 11.875000
		# ubiquitin_1   37.08791 0.8241758 11.904762 29.02930  8.241758 12.912088
		# ubiquitin_2   36.80982 2.4539877  9.202454 32.51534  5.521472 13.496933
		# ubiquitin_3   40.32258 0.0000000  6.451613 19.35484 14.516129 19.354839
		# ubiquitin_4   43.58974 0.0000000 10.256410 12.82051  7.692308 25.641026
		# ubiquitin_5   19.60784 0.0000000  7.843137 23.52941 17.647059 31.372549

#Figure S10. Fraction of genes with unknown duplication mechanisms in the ATG8, MUB, and SUMO subfamilies along with 5 Ub subfamilies that encode 1, 2, 3, 4, or 5 or more Ub-domains.

		pdf("Figure S10.pdf",family="Times",height=10,width=10)
	
		barplot(percent[,6],ylim=c(0,40),main="all_ubiquitin1_2_3_4_5_unknown_dup")

		dev.off()
	
	
	fisher_test(wgd,x,y,z)	
	
		#                                   [,1]
		# ubiquitin_5 unknown greater 0.03957457

#Figure 5F. Comparison of Ks and Ka/Ks values for the groups in Figure 5E.		

#Ks comparison	

	x<-c("ubiquitin_1","ubiquitin_2","ubiquitin_3","ubiquitin_4","ubiquitin_5")
	y<-3
	z<-6
	wgd<-d
	
	wilcoxtest(wgd,x,y,z)
	
		##wilcoxtest result
		#                            [,1]
		# ubiquitin_4 greater 6.783895e-05
		# ubiquitin_5 greater 1.114604e-07	
		##
		
	cluster_ks_df<-c()
	cluster_ks_df<-k_df(wgd,x,y,z)	
	attach(cluster_ks_df)	
	
	pdf("Figure 5F_Ks.pdf", family="Times",height=7, width=2.5)

	boxplot(cluster_ks_df$Measure ~ cluster_ks_df$Group, ylim=c(0,5), 
	cex.names=0.3,main="clusteri_wgd_ks_comparison",col=c("red","yellow"),outline=FALSE)
	
	dev.off()

#Ka/Ks comparison	
	
	x<-c("ubiquitin_1","ubiquitin_2","ubiquitin_3","ubiquitin_4","ubiquitin_5")
	y<-3	
	z<-7
	wgd<-d
	
	wilcoxtest(wgd,x,y,z)

		##wilcoxtest result	
		#                            [,1]
		# ubiquitin_1 greater 0.0040950285
		# ubiquitin_5 greater 0.0002849618
		##
		
	cluster_kaks_df<-c()
	cluster_kaks_df<-k_df(wgd,x,y,z)	
	attach(cluster_kaks_df)
	
pdf("Figure 5F_KaKs.pdf", family="Times",height=7, width=2.5)

	boxplot(cluster_kaks_df$Measure ~ cluster_kaks_df$Group, ylim=c(0,1.5), 
	cex.names=0.3,main="clusteri_wgd_kaks_comparison",col=c("red","yellow"),outline=FALSE)
	
	dev.off()



####

#Figure S6A. Fractions of ubiquiton genes that were originated from 6 indicated gene duplication mechanisms.  

#proportion of each WGD group


	wgd<-d

	dups<-table(wgd[5])/dim(wgd)[1]*100


	pdf("Figure S6A.pdf", family="Times",height=10, width=5)

	barplot(dups, ylim=c(0,50),main="DUP Types",col=c("red","gray"))
	
	dev.off()	

#Figure S6B. Ks comparison among different WGD groups

	x<-c("Cretaceous","paleo","neo","na","tandem","unknown")
	y<-5
	z<-6
	wgd<-d

	wilcoxtest(wgd,x,y,z)
	

		#                           [,1]
		# Cretaceous greater 4.819951e-05
		# paleo less         2.359616e-04
		# neo less           1.046235e-06
		# unknown greater    1.128757e-02
		
		
	cluster_ks_df<-c()
	cluster_ks_df<-k_df(wgd,x,y,z)	
	attach(cluster_ks_df)	


	pdf("Figure S6B.pdf", family="Times",height=7, width=2.5)

	boxplot(cluster_ks_df$Measure ~ cluster_ks_df$Group, ylim=c(0,3), 	
	cex.names=0.3,main="clusteri_wgd_ks_comparison",col=c("red","yellow"),outline=FALSE)

	dev.off()

#Figure S6C. Ka/Ks comparison among different WGD groups

	x<-c("Cretaceous","paleo","neo","na","tandem","unknown")
	y<-5	
	z<-7
	wgd<-d

	wilcoxtest(wgd,x,y,z)

		
		##wilcoxtest result		
		#                        [,1]
		# tandem greater 8.517936e-14

	cluster_kaks_df<-c()
	cluster_kaks_df<-k_df(wgd,x,y,z)	
	attach(cluster_kaks_df)

	pdf("Figure S6C.pdf", family="Times",height=7, width=2.5)

	boxplot(cluster_kaks_df$Measure ~ cluster_kaks_df$Group, ylim=c(0,1.5), 
	cex.names=0.3,main="clusteri_wgd_kaks_comparison",col=c("red","yellow"),outline=FALSE)


	dev.off()

###






