### rm(list=ls())

### Invoking requisite R libraries
library(shinydashboard, warn.conflicts = F)
library(shiny, warn.conflicts = F)
library(DT,warn.conflicts = F)

### Define path to Alignment Summary
rAlignSummary = reactive({paste(rRepoRoot(), "/output_repository/wrapper_align/", rPipelineID(), "_wrap/Summary.txt", sep="")})

output$alignsummary = renderText({ifelse(rRepoRoot() == '', "", rAlignSummary())})

### Read Alignment Summary
oAlignment = reactive({
	oDAT = read.delim(rAlignSummary(), header=T, sep="\t", stringsAsFactor=F)
	colnames(oDAT)[1] = "Sample.ID"
	rownames(oDAT) = oDAT$Sample.ID
	oDAT = oDAT[mixedorder(oDAT$Sample.ID),]
	return(oDAT)
})

### Alignment Graph 1
rAlignGraph1 = reactive({
	oDAT = subset(oAlignment(), select = c("Sample.ID","Total.Reads","Total.Mapped.Reads","Percent.Mapped.Reads"))
	oGraph = create_barplot1(oDAT)
	return(oGraph)
})

output$AlignmentGraph1 = renderPlot({
	print(rAlignGraph1())
})

### Download Graph 1
output$downloadPlot3 = downloadHandler(
							filename = function() {
								paste("Alignment_Map_Stats.1", ".png", sep="")
							},
							content = function(file) {
								ggsave(file, rAlignGraph1())
							}
					   )

### Alignment Graph 2
rAlignGraph2 = reactive({
	oDAT = oAlignment()[,c(1,7:(ncol(oAlignment()) - 3))]
	oGraph = create_barplot2(oDAT)
	return(oGraph)
})

output$AlignmentGraph2 = renderPlot({
	print(rAlignGraph2())
})

### Download Graph 2
output$downloadPlot4 = downloadHandler(
							filename = function() {
								paste("Alignment_Map_Stats.2", ".png", sep="")
							},
							content = function(file) {
								ggsave(file, rAlignGraph2())
							}
					   )

### Select Table Columns
output$column_select = renderUI({
	checkboxGroupInput(
		'show_columns',
		label = NULL,
		names(oAlignment()),
		selected = names(oAlignment()),
		inline = FALSE
	)
})

### Display Table
rColumns = reactive({input$show_columns})
output$alignment_table = DT::renderDataTable({
	DT::datatable(
		subset(oAlignment(), select = rColumns()),
		options = list(scrollX = TRUE)
	)
})

### Download Table
output$downloadTable1 = downloadHandler(
							filename = function() {
								paste("Alignment_Summary", 'txt', sep='.')
							},
							content = function(file) {
								write.table(subset(oAlignment(), select = rColumns()),
											file, append=F, quote=F, sep="\t", eol="\n", na="NA", row.names=F, col.names=T
								)
							}
						)

