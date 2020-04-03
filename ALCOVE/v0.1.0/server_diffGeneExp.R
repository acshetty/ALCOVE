#diffGeneExpTable_data <- reactive({
### rm(list=ls())

### Invoking requisite R libraries
library(shinydashboard, warn.conflicts = F)
library(shiny, warn.conflicts = F)
library(DT,warn.conflicts = F)

### Define path to Alignment Summary
output$selectuser2 = renderText({paste("Username :", rUser())})
output$selectproject2 = renderText({paste("Project :", rProject())})
output$selectpipeline2 = renderText({paste("Pipeline ID :", rPipelineID())})
output$repositoryroot2 = renderText({paste("Repository :", rRepoRoot())})

rDEGeneDir = reactive({paste(rRepoRoot(), "/output_repository/deseq/", rPipelineID(), "_*/i1", sep="")})

output$degeneset = renderText({ifelse(rRepoRoot() == '', "", rGeneExpDir())})

### Generate 'Comparison' Drop Down list
output$Select_Comparison = renderUI({
	selectInput(
		"comparison",
		label = h4("Select Comparison"),
		choices = listcomparisons(rDEGeneDir()),
		selected = defaultcomparison(rDEGeneDir()),
		width = "100%"
	)
})

output$degenedir = renderText({ifelse(rRepoRoot() == '', "", dirname(input$comparison))})

rDEGeneFile = reactive({
	input$comparison
})

output$degenefile = renderText({ifelse(rRepoRoot() == '', "", rDEGeneFile())})

### Generate DE Genes Matrix
rDEG = reactive({
	oDAT = read.delim(rDEGeneFile(), header=T, sep="\t", stringsAsFactor=F)
	colnames(oDAT)[1] = "ID"
	rownames(oDAT) = oDAT$ID
	
	### Eliminate Infinity values
	bINF = is.infinite(oDAT[,6])
	aLFC = oDAT[!(bINF),6]
	oDAT[(bINF & oDAT[,6] < 0),6] = floor(min(aLFC))
	oDAT[(bINF & oDAT[,6] > 0),6] = ceiling(max(aLFC))
	
	oDAT = oDAT[order(oDAT[,6]),]
	return(oDAT)
})

output$degenesummary = renderText({
	ifelse(rRepoRoot() == '', "", 
		paste(
			paste("Number of columns: ", ncol(rDEG()), sep="\t"),
			paste("Number of genes: ", nrow(rDEG()), sep="\t"),
			sep="\n"
		)
	)
})

### Generate Normalized Counts Matrix
rNormCount = reactive({
	sNormFile = paste(dirname(rDEGeneFile()), "all_counts_noZero_normalized", sep="/")
	oDAT = read.delim(sNormFile, header=T, sep="\t", stringsAsFactor=F)
	colnames(oDAT)[1] = "ID"
	rownames(oDAT) = oDAT$ID
	
	return(oDAT)
})

output$normcount = renderText({
	ifelse(rRepoRoot() == '', "", 
		paste(
			paste("Number of columns: ", ncol(rNormCount()), sep="\t"),
			paste("Number of genes: ", nrow(rNormCount()), sep="\t"),
			sep="\n"
		)
	)
})

### Render LFC Cut-off Slider
output$LFC_Slider <- renderUI({
	sliderInput(
		"LFC", 
		"Select Log2(Fold-Change) cut-off", 
		value = 1 , 
		min = 0, 
		max = max(abs(rDEG()[,6])), 
		animate = FALSE,
		round = TRUE
	)
})

rFDR_cutoff = eventReactive(input$filter, {
	input$FDR
})

rPVAL_cutoff = eventReactive(input$filter, {
	input$PVAL
})

rRCP_cutoff = eventReactive(input$filter, {
	input$RCP
})

rLFC_cutoff = eventReactive(input$filter, {
	input$LFC
})

### Generate Filtered DE Genes Matrix
rFilteredDEG = reactive({
	filter_de_genes(rDEG(), rFDR_cutoff(), rPVAL_cutoff(), rRCP_cutoff(), rLFC_cutoff())
})

output$degenefiltered = renderText({
	ifelse(rRepoRoot() == '', "", 
		paste(
			paste("Read count cut-off: ", quantile(rDEG()[,2], probs = rRCP_cutoff())[1], sep="\t"),
			paste("Number of columns: ", ncol(rFilteredDEG()), sep="\t"),
			paste("Number of genes: ", nrow(rFilteredDEG()), sep="\t"),
			sep="\n"
		)
	)
})

rFilteredNormCount = reactive({
	filter_norm_count(rFilteredDEG(), rNormCount())
})

output$normcountfiltered = renderText({
	ifelse(rRepoRoot() == '', "", 
		paste(
			paste("Number of columns: ", ncol(rFilteredNormCount()), sep="\t"),
			paste("Number of genes: ", nrow(rFilteredNormCount()), sep="\t"),
			sep="\n"
		)
	)
})

### MA Plot
rMAPlot = reactive({
	oGraph = plot_maplot(rDEG(), rFDR_cutoff(), rPVAL_cutoff(), rRCP_cutoff(), rLFC_cutoff())
	return(oGraph)
})

output$MAPlot1 = renderPlot({
	print(rMAPlot())
})

### Download MA Plot
output$downloadPlot8 = downloadHandler(
							filename = function() {
								paste("MA_Plot", ".png", sep="")
							},
							content = function(file) {
								ggsave(file, rMAPlot())
							}
					   )

### Heat Map
rHeatMap = reactive({
	oGraph = plot_heatmap(rFilteredNormCount())
	return(oGraph)
})

output$HeatMap2 = renderPlot({
	print(rHeatMap())
})

### Download Heat Map
output$downloadPlot9 = downloadHandler(
							filename = function() {
								paste("Heatmap", ".png", sep="")
							},
							content = function(file) {
								ggsave(file, rHeatMap())
							}
					   )

### Display Table
output$diff_gene_expression = DT::renderDataTable({
	DT::datatable(
		rFilteredDEG(),
		options = list(scrollX = TRUE)
	)
})

### Download Table
output$downloadTable3 = downloadHandler(
							filename = function() {
								paste("Differential_Gene_Expression", 'txt', sep='.')
							},
							content = function(file) {
								write.table(rFilteredDEG(),
											file, append=F, quote=F, sep="\t", eol="\n", na="NA", row.names=F, col.names=T
								)
							}
						)

