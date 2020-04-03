### FASTQC TAB ###

listsamples = function(sDir) {
				aGroupDirs = list.dirs(sDir, full.names = T, recursive = F)
				aGroupDirs = mixedsort(aGroupDirs)
				
				aSampleNames = vector()
				for (i in 1:length(aGroupDirs)) {
					sG = aGroupDirs[i]
					sFile = Sys.glob(paste(sG,"*_fastqc/Images","*.per_base_quality.png", sep="/"))[1]
					sFile = gsub("_1_1_sequence.per_base_quality.png", "", basename(sFile))
					aSampleNames = c(aSampleNames, sFile)
				}
				
				names(aGroupDirs) = aSampleNames
				aGroupDirs = aGroupDirs[mixedorder(aSampleNames)]
				
				return(aGroupDirs)
			  }

defaultsample = function(sDir) {
					aGroupDirs = listsamples(sDir)
					sSampleID = names(aGroupDirs)[1]
					return(sSampleID)
				}

listimages = function(sDir) {
				sImage = Sys.glob(paste(sDir,"*_fastqc/Images","*.per_base_quality.png", sep="/"))[1]
				aImages = Sys.glob(paste(dirname(sImage),"*.png", sep="/"))
				
				names(aImages) = basename(aImages)
				aImages = mixedsort(aImages)
				
				return(aImages)
			 }

defaultimage = function(sDir) {
					sImage = Sys.glob(paste(sDir,"*_fastqc/Images","*.per_base_quality.png", sep="/"))[1]
					sImage = basename(sImage)
					return(sImage)
			   }

### ALIGNMENT TAB ###

create_barplot1 = function(oDF) {
					oP = ggplot(data = oDF, aes(x = Sample.ID, y = (Total.Mapped.Reads/1000000), fill = Sample.ID))
					oP = oP + geom_bar(stat = "identity", position = "dodge")
					oP = oP + ylab("Total Mapped Reads (in millions)")
					oP = oP + theme_bw()
					oP = oP + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
					return(oP)
				 }

create_barplot2 = function(oDF) {
					oDAT = melt(oDF, id.vars = c("Sample.ID"))
					oP = ggplot(data = oDAT, aes(x = Sample.ID, y = (value*100), colour = variable, fill = variable))
					oP = oP + geom_bar(stat = "identity", position = "fill")
					oP = oP + ylab("Percent Mapped Reads")
					oP = oP + theme_bw()
					oP = oP + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "top", legend.title=element_blank())
					return(oP)
				 }

### GENE EXPRESSION TAB ###

readCountFiles = function(aFiles,nColNum=2){
					## read all files
					oRawCounts = lapply(1:length(aFiles),
										function(x) {
											read.delim(aFiles[x], header=FALSE, stringsAsFactor=FALSE)
										}
								 )
					
					## check to make sure the same genes are in the files
					sGenes = oRawCounts[[1]][,1]
					for (i in 2:length(oRawCounts)) {
						if( !all(oRawCounts[[i]][,1] == sGenes) ) {
							stop("The list of genes in the first column of file ", i, "does not exactly match that of file 1.\n")
						}
					}
					
					## make a nice table
					oCounts = as.data.frame(sapply(oRawCounts, function(x) x[,nColNum]))
					rownames(oCounts) = sGenes
					colnames(oCounts) = unlist(lapply(basename(aFiles), FUN=function(x) {strsplit(x, "[.]")[[1]][1]}))
					oCounts = oCounts[,mixedorder(colnames(oCounts))]
					return(oCounts)
				 }

generate_count_matrix = function(sDir) {
							aFiles = Sys.glob(paste(sDir,"g*","*.counts", sep="/"))
							oCountTable = readCountFiles(aFiles,nColNum=2)
							
							#idx = which(rownames(oCountTable)=="no_feature")
							idx = grep("no_feature", rownames(oCountTable))
							oCountTable = oCountTable[1:(idx-1),]
							
							idx = which((rowSums(oCountTable)/ncol(oCountTable)) >= 1)
							oCountTable = oCountTable[idx,]
							
							return(oCountTable)
						}

generate_cpm_matrix = function(oCount) {
						oCPM = (as.matrix(oCount + 1) %*% diag(1/colSums(oCount))) * 1000000
						oCPM = log10(oCPM)
						oCPM = apply(oCPM, 2, function(x) {round(x, digits = 4)})
						oCPM = as.data.frame(oCPM)
						rownames(oCPM) = rownames(oCount)
						colnames(oCPM) = colnames(oCount)
						return(oCPM)
					  }

filter_cpm_matrix = function(oCPM, nCutOff) {
						idx = which(rowSums(oCPM >= nCutOff) > 0)
						oCPM.filter = oCPM[idx,]
						oCPM.filter = apply(oCPM.filter, 2, function(x) {round(x, digits = 4)})
						oCPM.filter = as.data.frame(oCPM.filter)
						rownames(oCPM.filter) = rownames(oCPM)[idx]
						colnames(oCPM.filter) = colnames(oCPM)
						return(oCPM.filter)
					}

plot_density = function(oDF, nCutOff) {
					oMLT = cbind(rownames(oDF), oDF)
					colnames(oMLT) = c("Gene.ID", colnames(oDF))
					oMLT = melt(oMLT, id=c("Gene.ID"))
					colnames(oMLT) = c("Gene.ID", "Sample.ID", "Log10.CPM")
					
					oP = ggplot(data = oMLT, aes(x = Log10.CPM, color = Sample.ID))
					oP = oP + geom_density(show.legend = FALSE)
					#oP = oP + scale_fill_hue()
					oP = oP + theme_bw()
					oP = oP + theme(legend.key = element_blank())
					oP = oP + geom_vline(xintercept = nCutOff)
					oP = oP + stat_density(aes(x = Log10.CPM, color = Sample.ID), geom = "line", position = "identity")
					return(oP)
			   }

plot_boxplot = function(oDF) {
					oMLT = cbind(rownames(oDF), oDF)
					colnames(oMLT) = c("Gene.ID", colnames(oDF))
					oMLT = melt(oMLT, id=c("Gene.ID"))
					colnames(oMLT) = c("Gene.ID", "Sample.ID", "Expression")
					
					return(
						boxplot(Expression ~ Sample.ID,
							data = oMLT,
							las = 2,
							par(mar = c(8, 4, 4, 2) + 0.1),
							outline = F,
							col = rainbow(ncol(oDF)),
							main = "Gene Expression distribution across samples"
						)
					)
			   }

plot_pcaplot = function(oDF) {
					oEXPR = t(as.matrix(oDF))
					oPCA = prcomp(oEXPR)
					
					oEXPR = as.data.frame(cbind(colnames(oDF), oEXPR))
					colnames(oEXPR)[1] = "Sample.ID"
					oEXPR$Sample.ID = factor(oEXPR$Sample.ID, levels=mixedsort(oEXPR$Sample.ID))
					
					oP = autoplot(oPCA, data = oEXPR, colour = 'Sample.ID')
					oP = oP + theme_bw()
					oP = oP + xlab("Principle Component 1")
					oP = oP + ylab("Principle Component 2")
					return(oP)
			   }

### DIFFERENTIAL GENE EXPRESSION TAB ###

listcomparisons = function(sDir) {
						aFiles = Sys.glob(paste(sDir,"g*","*.de_genes.txt", sep="/"))
						
						names(aFiles) = gsub(".de_genes.txt", "", basename(aFiles))
						aFiles = aFiles[mixedorder(names(aFiles))]
						
						return(aFiles)
				  }

defaultcomparison = function(sDir) {
						aComparisons = listcomparisons(sDir)
						sComparison = names(aComparisons)[1]
						return(sComparison)
					}

filter_de_genes = function(oDEG, nFDR, nPVAL, nRCP, nLFC) {
					nRC = quantile(oDEG[,2], probs = nRCP)[1]
					idx = which((oDEG[,8] <= nFDR) & (oDEG[,7] <= nPVAL) & ((oDEG[,3] >= nRC) | (oDEG[,4] >= nRC)) & (abs(oDEG[,6]) >= nLFC))
					oDEG.filter = oDEG[idx,]
					rownames(oDEG.filter) = rownames(oDEG)[idx]
					colnames(oDEG.filter) = colnames(oDEG)
					return(oDEG.filter)
				  }

filter_norm_count = function(oDEG, oNRM) {
						oNRM.filter = oNRM[oDEG$ID,]
						return(oNRM.filter)
					}

plot_maplot = function(oDEG, nFDR, nPVAL, nRCP, nLFC) {
					nRC = quantile(oDEG[,2], probs = nRCP)[1]
					bDEG = ((oDEG[,8] <= nFDR) & (oDEG[,7] <= nPVAL) & ((oDEG[,3] >= nRC) | (oDEG[,4] >= nRC)) & (abs(oDEG[,6]) >= nLFC))
					aCOL = rep("grey", nrow(oDEG))
					aCOL[oDEG[,6] <= (-1*nLFC)] = "deeppink3"
					aCOL[oDEG[,6] >= nLFC] = "chartreuse3"
					
					colnames(oDEG)[c(2,6)] = c("Avg.Norm.Expr","Log2FC")
					oDEG$COL = aCOL
					
					oP = ggplot(data = oDEG, aes(x = Avg.Norm.Expr, y = Log2FC, color = COL))
					oP = oP + geom_point(size = 0.5, show.legend = FALSE)
					oP = oP + theme_bw()
					oP = oP + xlab("Average Normalized Expression")
					oP = oP + ylab("Log2(Fold-change)")
					oP = oP + scale_x_log10()
					oP = oP + scale_color_manual(values = levels(factor(oDEG$COL)))
					oP = oP + theme(legend.key = element_blank())
					oP = oP + geom_hline(yintercept = c(-1*nLFC, 0, nLFC))
					return(oP)
			  }

plot_heatmap = function(oDF) {
					oMAT = as.matrix(oDF[,c(2:ncol(oDF))])
					aNames = colnames(oDF)[c(2:ncol(oDF))]
					
					aSD = apply(oMAT, 1, sd)
					oMAT = oMAT[order(-aSD),]
					if (nrow(oMAT) > 2000) {
						oMAT = oMAT[1:2000,]
					}
					
					my_palette = colorRampPalette(c("red", "yellow", "green"))(n = 299)
					return(
						heatmap.2(
							oMAT,
							main = paste("Heatmap of DE Genes", paste("(top ", nrow(oMAT), " genes)", sep=""), sep="\n"), 
							trace="none", 
							margins=c(10,2), 
							col=my_palette, 
							dendrogram="both", 
							scale="row", 
							keysize = 1.5, 
							density.info="none", 
							cexRow=0.8, cexCol=0.8, 
							labRow=NA, 
							labCol=aNames, 
						)
					)
			   }

overlap_de_genes = function(aLST, nFDR, nPVAL, nRCP, nLFC) {
						sFile = aLST[1]
						sID = gsub(".de_genes.txt", "", basename(sFile))
						oDAT = read.delim(sFile, header=T, sep="\t", stringsAsFactor=F)
						colnames(oDAT)[1] = "ID"
						rownames(oDAT) = oDAT$ID
						
						### Eliminate Infinity values
						oDAT = eliminate_infinity_values(oDAT)
						oDEG = filter_de_genes(oDAT, nFDR, nPVAL, nRCP, nLFC)
						colnames(oDEG)[6] = sID
						oOVLP = oDEG[,c(1,6)]
						
						for (nI in c(2:length(aLST))) {
							sFile = aLST[nI]
							sID = gsub(".de_genes.txt", "", basename(sFile))
							oDAT = read.delim(sFile, header=T, sep="\t", stringsAsFactor=F)
							colnames(oDAT)[1] = "ID"
							rownames(oDAT) = oDAT$ID
							
							### Eliminate Infinity values
							oDAT = eliminate_infinity_values(oDAT)
							oDEG = filter_de_genes(oDAT, nFDR, nPVAL, nRCP, nLFC)
							colnames(oDEG)[6] = sID
							
							oOVLP = merge(oOVLP, oDEG[,c(1,6)], by="ID", all=T)
						}
						
						return(oOVLP)
				   }

plot_venn = function(oDF) {
				oID = list()
				for (nI in c(2:ncol(oDF))) {
					sCID = colnames(oDF)[nI]
					aGID = oDF[!(is.na(oDF[,nI])), "ID"]
					oID[[sCID]] = unique(aGID)
				}
				
				oVenn = venn(oID, small=nCS, simplify=FALSE, show.plot=FALSE)
				oPC = format(oVenn[, "num"] / sum(oVenn[, "num"]) * 100, digits=2, nsmall=2)
				oPC = gsub(" ", "", oPC)
				oVenn[, "num"] = paste(oVenn[, "num"], "\n(", oPC, "%)", sep="")
				
				return(
					oVenn
				)
			}

paircomparisons = function(oDF) {
						aComparisons = colnames(oDF)[2:ncol(oDF)]
						return(aComparisons)
				  }

plot_quadrant = function(oDF, sCMP1, sCMP2) {
					oDEG = oDF[,c("ID",sCMP1,sCMP2)]
					oDEG[is.na(oDEG)] = 0
					oDEG = oDEG[rowSums(oDEG[,c(sCMP1,sCMP2)] == 0) < 2,]
					
					if(nrow(oDEG) > 1) {
						colnames(oDEG) = c("ID","CMP1","CMP2")
						aBIN = paste(ifelse(oDEG$CMP1<0,"1",ifelse(oDEG$CMP1>0,"2","0")),
									 ifelse(oDEG$CMP2<0,"1",ifelse(oDEG$CMP2>0,"2","0")),
									 sep="")
						oDEG$BIN = aBIN
						aLevels = c("22","12","11","21","20","02","10","01")
						oDEG$BIN = factor(oDEG$BIN, levels = aLevels)
						
						aColors = c("darkgreen","darkblue","darkred","goldenrod1","tan3","cyan","magenta","orange")
						names(aColors) = aLevels
						
						aLabels = c("UP|UP","DOWN|UP","DOWN|DOWN","UP|DOWN","UP|NONE","NONE|UP","DOWN|NONE","NONE|DOWN")
						aLabels = paste(aLabels, unlist(table(oDEG$BIN)), sep=":")
						names(aColors) = aLevels
						
						oP = ggplot(oDEG, aes(CMP1, CMP2, colour=BIN))
						oP = oP + geom_point(size=0.5, alpha=0.5)
						oP = oP + scale_colour_manual(values = aColors, labels = aLabels)
						oP = oP + theme_bw()
						oP = oP + geom_vline(xintercept = 0, color="black", alpha=0.2)
						oP = oP + geom_hline(yintercept = 0, color="black", alpha=0.2)
						oP = oP + ggtitle("")
						oP = oP + labs(x=paste(sCMP1, ": Log2(Fold-change)", sep=""), y=paste(sCMP2, ": Log2(Fold-change)", sep=""))
						oP = oP + theme(legend.position="right", legend.title=element_blank())
						return(oP)
					}
				}
							
### MULTIPLE TABS ###

eliminate_infinity_values = function(oDF) {
								bINF = is.infinite(oDF[,6])
								aLFC = oDF[!(bINF),6]
								oDF[(bINF & oDF[,6] < 0),6] = floor(min(aLFC))
								oDF[(bINF & oDF[,6] > 0),6] = ceiling(max(aLFC))
								return(oDF)
							}
								

add_ensembl_links = function(oDF) {
						oDF[,1] = paste(
									"<a href='http://useast.ensembl.org/Multi/Search/Results?q=", 
									oDF[,1], 
									";site=ensembl' ", 
									"target='_blank'>", 
									oDF[,1], 
									"</a>", 
									sep = ""
								  )
						return(oDF)
					}

