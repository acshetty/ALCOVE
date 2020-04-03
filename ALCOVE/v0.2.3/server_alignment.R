### rm(list=ls())

### Invoking requisite R libraries
library(shinydashboard, warn.conflicts = F)
library(shiny, warn.conflicts = F)
library(DT,warn.conflicts = F)
library(gridExtra)
### Upload Info File 



### Define path to Alignment Summary

#observeEvent(input$useergatis,{
rAlignSummary = renderText({paste(rRepoRoot(), "/output_repository/wrapper_align/", rPipelineID(), "_wrap/Summary.txt", sep="")})
output$summary = renderText({ rAlignSummary() })
output$alignsummary = renderText({ifelse(rRepoRoot() == '', "", rAlignSummary())})

oAlignment = reactive({
    oDAT = read.delim(rAlignSummary(), header=T, sep="\t", stringsAsFactor=F)
    colnames(oDAT)[1] = "Sample.ID"
    rownames(oDAT) = oDAT$Sample.ID
    oDAT = oDAT[mixedorder(oDAT$Sample.ID),]
    return(oDAT)
})

### Read Info File
oInfo = reactive({
    oDAT = read.delim(rInfoFile(), header=F, sep="\t", stringsAsFactor=F)
    colnames(oDAT)[1] = "Sample.ID"
    colnames(oDAT)[2] = "Condition"
    return(oDAT)
    })
#### Some variable to check if everything works
oPullInfo = reactive({
    #oDAT1 = subset(oAlignment(), select = c("Sample.ID","Percent.Mapped.Reads"))
    oDAT1 = subset(oAlignment(), select = c("Sample.ID","Total.Reads","Total.Mapped.Reads", "Percent.Mapped.Reads"))
    oDAT2 = subset(oInfo(), select = c("Sample.ID", "Condition"))
    oDAT = merge(oDAT1, oDAT2, by = "Sample.ID")
    #colnames(oDAT)<-c("SampleID", "Total.Reads", "Total.Mapped", "Percent.Mapped")
    oDAT$Unmapped = (oDAT$Total.Reads-oDAT$Total.Mapped)
    #oDAT = melt(oDAT, id.vars = c("Sample.ID", "Condition"))
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

###Alignment with info file

rAlignInfo1 = reactive({
    oDAT1 = subset(oAlignment(), select = c("Sample.ID","Total.Reads","Total.Mapped.Reads"))
    oDAT2 = subset(oInfo(), select = c("Sample.ID", "Condition"))
    oDAT = merge(oDAT1, oDAT2, by = "Sample.ID")
    #oDAT = melt(oDAT, id.vars = c("Sample.ID", "Condition"))
    oGraph = create_info_barplot1(oDAT)
    return(oGraph)
    })

output$AlignInfoGraph1 = renderPlot({
        print(rAlignInfo1())
})

rAlignInfo2 = reactive({
        oDAT1 = subset(oAlignment(), select = c("Sample.ID","Percent.Mapped.Reads"))
        oDAT2 = subset(oInfo(), select = c("Sample.ID", "Condition"))
        oDAT = merge(oDAT1, oDAT2, by = "Sample.ID")
        oGraph = create_info_barplot2(oDAT)
        return(oGraph)
        })

output$AlignInfoGraph2 = renderPlot({
        print(rAlignInfo2())
})

rAlignCombo1 = reactive({
    oDAT = subset(oAlignment(), select = c("Sample.ID","Total.Reads","Total.Mapped.Reads", "Percent.Mapped.Reads"))
    colnames(oDAT) = c("SampleID", "Total.Reads", "Total.Mapped", "Percent.Mapped")
    oDAT$Unmapped = (oDAT$Total.Reads-oDAT$Total.Mapped)
    oGraph = create_combo_barplot(oDAT)
    return(oGraph)
})

output$AlignComboGraph1 = renderPlot({
        print(rAlignCombo1())
})

## Download Graph 1
output$downloadPlot3 = downloadHandler(
                            filename = function() {
                                paste("Alignment_Map_Stats.1", ".png", sep="")
                            },
                            content = function(file) {
                                ggsave(file, rAlignGraph1())
                            }
                       )
#Download Combo Plot1
output$downloadPlot5_1 = downloadHandler(
                            filename = function() {
                                paste("Alignment_Reads_Grp.1", ".png", sep="")
                            },
                            content = function(file) {
                                ggsave(file, rAlignInfo1())
                            }
                       )
#Downlad Combo plot 2
output$downloadPlot6_1 = downloadHandler(
                            filename = function() {
                                paste("Alignment_combo.1", ".png", sep="")
                            },
                            content = function(file) {
                                ggsave(file, rAlignCombo1())
                            }
                       )

output$downloadPlot7_1 = downloadHandler(
                            filename = function() {
                                paste("Alignment_Pmapped_Grp.1", ".png", sep="")
                            },
                            content = function(file) {
                                ggsave(file, rAlignInfo2())
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

### Info Table Display

output$info_table = DT::renderDataTable({
        DT::datatable(oPullInfo(),
                options = list(scrollX = TRUE)
        )
})

output$info_table1 = DT::renderDataTable({
        DT::datatable(oInfo(),
                options = list(scrollX = TRUE)
        )
})

#})
#############################
###Alignment for BDBAG#######
#############################

observeEvent(input$unzip,{
rAlignSummary <- renderText({paste(rReporoot(), "/data/Summary.txt", sep="")})
output$summary = renderText({ rAlignSummary() })
output$alignsummary = renderText({ifelse(rReporoot() == '', "", rAlignSummary())})

### Read Alignment Summary
oAlignment = reactive({
	oDAT = read.delim(rAlignSummary(), header=T, sep="\t", stringsAsFactor=F)
	colnames(oDAT)[1] = "Sample.ID"
	rownames(oDAT) = oDAT$Sample.ID
	oDAT = oDAT[mixedorder(oDAT$Sample.ID),]
	return(oDAT)
})

#observeEvent(input$infoF,{
#                    output$infoF <- renderText({ input$infoF$datapath })
#                    rInfofile <- renderText({ input$infoF$datapath })

### Read Info File
oInfo = reactive({
	#oDAT = read.delim(rInfofile(), header=F, sep="\t", stringsAsFactor=F)
	oDAT = read.delim(input$infoF$datapath, header=F, sep="\t", stringsAsFactor=F)
    colnames(oDAT)[1] = "Sample.ID"
	colnames(oDAT)[2] = "Condition"
	return(oDAT)
	})

#### Some variable to check if everything works
oPullInfo = reactive({
	#oDAT1 = subset(oAlignment(), select = c("Sample.ID","Percent.Mapped.Reads"))
	oDAT1 = subset(oAlignment(), select = c("Sample.ID","Total.Reads","Total.Mapped.Reads", "Percent.Mapped.Reads"))
	oDAT2 = subset(oInfo(), select = c("Sample.ID", "Condition"))
	oDAT = merge(oDAT1, oDAT2, by = "Sample.ID")
	#colnames(oDAT)<-c("SampleID", "Total.Reads", "Total.Mapped", "Percent.Mapped")
	oDAT$Unmapped = (oDAT$Total.Reads-oDAT$Total.Mapped)
	#oDAT = melt(oDAT, id.vars = c("Sample.ID", "Condition"))
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

###Alignment with info file

rAlignInfo1 = reactive({
	oDAT1 = subset(oAlignment(), select = c("Sample.ID","Total.Reads","Total.Mapped.Reads"))
	oDAT2 = subset(oInfo(), select = c("Sample.ID", "Condition"))
	oDAT = merge(oDAT1, oDAT2, by = "Sample.ID")
	#oDAT = melt(oDAT, id.vars = c("Sample.ID", "Condition"))
	oGraph = create_info_barplot1(oDAT)
	return(oGraph)
	})

output$AlignInfoGraph1 = renderPlot({
        print(rAlignInfo1())
})

rAlignInfo2 = reactive({
        oDAT1 = subset(oAlignment(), select = c("Sample.ID","Percent.Mapped.Reads"))
        oDAT2 = subset(oInfo(), select = c("Sample.ID", "Condition"))
        oDAT = merge(oDAT1, oDAT2, by = "Sample.ID")
        oGraph = create_info_barplot2(oDAT)
        return(oGraph)
        })

output$AlignInfoGraph2 = renderPlot({
        print(rAlignInfo2())
})

rAlignCombo1 = reactive({
	oDAT = subset(oAlignment(), select = c("Sample.ID","Total.Reads","Total.Mapped.Reads", "Percent.Mapped.Reads"))
	colnames(oDAT) = c("SampleID", "Total.Reads", "Total.Mapped", "Percent.Mapped")
	oDAT$Unmapped = (oDAT$Total.Reads-oDAT$Total.Mapped)	
	oGraph = create_combo_barplot(oDAT)
	return(oGraph)
})

output$AlignComboGraph1 = renderPlot({
        print(rAlignCombo1())
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



#Download Combo Plot5
output$downloadPlot5_1 = downloadHandler(
                            filename = function() {
                                paste("Alignment_Reads_Grp.1", ".png", sep="")
                            },
                            content = function(file) {
                                ggsave(file, rAlignInfo1())
                            }
                       )
#Downlad Combo plot 7
output$downloadPlot6_1 = downloadHandler(
                            filename = function() {
                                paste("Alignment_combo.1", ".png", sep="")
                            },
                            content = function(file) {
                                ggsave(file, rAlignCombo1())
                            }
                       )

#Download Plot 6
output$downloadPlot7_1 = downloadHandler(
                            filename = function() {
                                paste("Alignment_Pmapped_Grp.1", ".png", sep="")
                            },
                            content = function(file) {
                                ggsave(file, rAlignInfo2())
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

### Info Table Display

output$info_table = DT::renderDataTable({
        DT::datatable(oPullInfo(),
                options = list(scrollX = TRUE)
        )
})

output$info_table1 = DT::renderDataTable({
        DT::datatable(oInfo(),
                options = list(scrollX = TRUE)
        )
})
})
