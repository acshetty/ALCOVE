### rm(list=ls())

### Invoking requisite R libraries
library(shinydashboard, warn.conflicts = F)
library(shiny, warn.conflicts = F)
library(DT,warn.conflicts = F)

### Define path to Alignment Summary
rGeneExpDir = reactive({paste(rRepoRoot(), "/output_repository/htseq/", rPipelineID(), "_*/i1", sep="")})

output$expression = renderText({ifelse(rRepoRoot() == '', "", rGeneExpDir())})

### Generate Counts Matrix
rRaw.Count = reactive({
	generate_count_matrix(rGeneExpDir())
})

output$dimensions = renderText({
	ifelse(rRepoRoot() == '', "", 
		paste(
			paste("Number of samples: ", ncol(rRaw.Count()), sep="\t"),
			paste("Number of genes: ", nrow(rRaw.Count()), sep="\t"),
			sep="\n"
		)
	)
})

### Generate CPM Matrix
rCPM = reactive({
	generate_cpm_matrix(rRaw.Count())
})

output$normalized = renderText({
	ifelse(rRepoRoot() == '', "", 
		paste(
			paste("Number of samples: ", ncol(rCPM()), sep="\t"),
			paste("Number of genes: ", nrow(rCPM()), sep="\t"),
			paste("Minimum log10(CPM): ", round(min(rCPM()), digits = 2), sep="\t"),
			paste("Maximum log10(CPM): ", round(max(rCPM()), digits = 2), sep="\t"),
			sep="\n"
		)
	)
})

### Render CPM Cut-off Slider
output$CPM_Slider = renderUI({
	sliderInput(
		"CPM_cutoff",
		"Select log10(CPM) cut-off",
		value = 0,
		min = round(min(rCPM()), digits = 2),
		max = round(max(rCPM()), digits = 2),
		animate = FALSE,
		round = TRUE
	)
})

rCPM_cutoff = eventReactive(input$submit, {
	input$CPM_cutoff
})

### Filter CPM Matrix
rCPM.filtered = reactive({
	filter_cpm_matrix(rCPM(), rCPM_cutoff())
})

output$filtered = renderText({
	ifelse(rRepoRoot() == '', "", 
		paste(
			paste("log10(CPM) filter: ", rCPM_cutoff(), sep="\t"),
			paste("Number of samples: ", ncol(rCPM.filtered()), sep="\t"),
			paste("Number of genes: ", nrow(rCPM.filtered()), sep="\t"),
			paste("Minimum log10(CPM): ", round(min(rCPM.filtered()), digits = 2), sep="\t"),
			paste("Maximum log10(CPM): ", round(max(rCPM.filtered()), digits = 2), sep="\t"),
			sep="\n"
		)
	)
})

### CPM Density Plot
rDensityPlot = reactive({
	oGraph = plot_density(rCPM(), rCPM_cutoff())
	return(oGraph)
})

output$DensityPlot1 = renderPlot({
	print(rDensityPlot())
})

### Download Density Plot
output$downloadPlot5 = downloadHandler(
							filename = function() {
								paste("Density_Plot.CPM", ".png", sep="")
							},
							content = function(file) {
								ggsave(file, rDensityPlot())
							}
					   )

### Select Dataset for Plots
rPlotData = reactive({
	switch(
		input$normalization,
		"CPM Data" = rCPM.filtered(),
		"Raw Data" = rRaw.Count()
	)
})

rNormalization = reactive({
	switch(
		input$normalization,
		"CPM Data" = "Filtered CPM Values",
		"Raw Data" = "Raw Count Values"
	)
})

### Box Plot
rBoxPlot = reactive({
	oGraph = plot_boxplot(rPlotData())
	return(oGraph)
})

output$BoxPlot2 = renderPlot({
	print(rBoxPlot())
})

### Download Box Plot
output$downloadPlot6 = downloadHandler(
							filename = function() {
								paste("Box_Plot", ".png", sep="")
							},
							content = function(file) {
								ggsave(file, rBoxPlot())
							}
					   )

### PCA Plot
rPCAPlot = reactive({
	oGraph = plot_pcaplot(rPlotData())
	return(oGraph)
})

output$PCAPlot3 = renderPlot({
	print(rPCAPlot())
})

### Download PCA Plot
output$downloadPlot7 = downloadHandler(
							filename = function() {
								paste("PCA_Plot", ".png", sep="")
							},
							content = function(file) {
								ggsave(file, rPCAPlot())
							}
					   )

### Display Table
output$gene_expression = DT::renderDataTable({
	DT::datatable(
		rPlotData(),
		options = list(scrollX = TRUE)
	)
})

### Download Table
output$downloadTable2 = downloadHandler(
							filename = function() {
								paste("Gene_Expression", 'txt', sep='.')
							},
							content = function(file) {
								write.table(rPlotData(),
											file, append=F, quote=F, sep="\t", eol="\n", na="NA", row.names=F, col.names=T
								)
							}
						)
