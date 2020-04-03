### rm(list=ls())

### Invoking requisite R libraries
library(shinydashboard, warn.conflicts = F)
library(shiny, warn.conflicts = F)

### Define path to FASTQC directory
rFastQCPath = reactive({paste(rRepoRoot(), "/output_repository/fastqc_stats/", rPipelineID(), "_fastqc/i1", sep="")})

output$fastqcdir = renderText({ifelse(rRepoRoot() == '', "", rFastQCPath())})

### Generate 'Samples' Drop Down list
output$sample_select = renderUI({
	selectInput(
		"sample",
		label = h4("Select Sample to be displayed"),
		choices = listsamples(rFastQCPath()),
		selected = defaultsample(rFastQCPath()),
		width = "100%"
	)
})

### Generate 'Images' Drop Down list
rSample = reactive({input$sample})
output$image_select = renderUI({
	selectInput(
		"image",
		label = h4("Select the image to be displayed"),
		choices = listimages(rSample()),
		selected = defaultimage(rSample()),
		width = "100%"
	)
})

### Generate Mate 1 Image
rImage1 = reactive({input$image})

output$set1Img = renderImage({
	list(src = rImage1(), width = "100%")
}, deleteFile = F)

### Generate Mate 2 Image
rImage2 = reactive({gsub("1_1_sequence","2_1_sequence",input$image)})

output$set2Img = renderImage({
	list(src = rImage2(), width = "100%")
}, deleteFile = F)

### Download Mate 1 Image
output$downloadPlot1 = downloadHandler(
							filename = function() {
								paste("Set1_", basename(rImage1()), sep="")
							},
							content = function(file) {
								src = rImage1()
								
								# temporarily switch to the temp dir, in case you do not have write
								# permission to the current working directory
								owd = setwd(tempdir())
								on.exit(setwd(owd))
								file.copy(src, file)
							},
							contentType = 'image/png'
						)

### Download Mate 2 Image
output$downloadPlot2 = downloadHandler(
							filename = function() {
								paste("Set2_", basename(rImage2()), sep="")
							},
							content = function(file) {
								src = rImage2()
								
								# temporarily switch to the temp dir, in case you do not have write
								# permission to the current working directory
								owd = setwd(tempdir())
								on.exit(setwd(owd))
								file.copy(src, file)
							},
							contentType = 'image/png'
						)

