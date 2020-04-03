### rm(list=ls())

### Invoking requisite R libraries
library(shinydashboard, warn.conflicts = F)
library(shiny, warn.conflicts = F)

### Define path to FASTQC directory
###Determine if BDBag or not
#By default assume ergatis, if bdbag click button and set flag to True. Process accordingly.

output$fqc<-renderText({ if(input$usebdbag){return(TRUE)} else {return(FALSE)}})

#observeEvent(input$useergatis,{
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
#})

##### If BDBAG#####
observeEvent(input$unzip,{
rFastQCPath = reactive({paste(rReporoot(),"/data/FastQC_Files/" , sep="")})
#rFastQCPath = reactive(gsub("//","/",{paste(rReporoot(),"/data/Files/" , sep="")}))
output$fastqcdir = renderText({rFastQCPath()})
### Generate 'ImageType' Drop Down list
output$sample_select = renderUI({
	selectInput(
		"type",
		label = h4("Select Type to be displayed"),
		choices = listImageTypes(rFastQCPath()),
		width = "100%"
	)
})
### Generate 'Samples' Drop Down list
rSample = reactive({input$type})
rSample_s = reactive({paste0(rFastQCPath(), "/",rSample())})
output$image_select = renderUI({
	selectInput(
		"image",
		label = h4("Select the sample to be displayed"),
        choices = listbdbagSamples(rSample()),
		width = "100%"
	)
})

### Generate Mate 1 Image
rImage = reactive({input$image})
get_images = reactive({ paste0(rSample(), "/", rImage(), "*") })
rImage1 = renderText({getImages1(get_images())})
output$holdimg = renderText({rImage1()})
output$set1Img = renderImage({
	list(src = rImage1(), width = "100%")
}, deleteFile = F)


### Generate Mate 2 Image
rImage2 = reactive({gsub("1_1_sequence","2_1_sequence",rImage1())})

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
})
### Table with explanations of images
StaticTable = reactive({
			Image_Type = c("Adapter Content", "Base Quality", "Duplication Levels", "Kmer Profiles", "Per Base N Content", "Per Base Sequence Content", "Per Sequence GC Content", "Per Sequence Quality Scores", "Per Tile Sequence Quality", "Sequence Length Distribution")
			Description = c("Illustrates the total proportion of your library that contains Adapter Kmers",
                         "Illustrates the distribution of per base quality across all reads for the samples",
                         "Illustrates the degree of duplication for every sequence in a library",
                         "Illustrate the occurrence of the top 6 kmers across all reads for the sample",
                         "Illustrates the percentage of base calls at each position for which an N was called",
                         "Illustrates the proportion of each base position for which each of the four normal DNA bases has been called.",
                         "Illustrates GC content across the whole length of each sequence",
                         "Illustrates if a subset of your sequences have universally low quality values",
                         "Illustrates the quality scores from each tile across all the bases",
                         "Illustrates the distribution of fragment sizes")

			STable = data.frame(Image_Type, Description)
			return(STable)
})
  
output$static_table = DT::renderDataTable({
        DT::datatable(StaticTable(),
                options = list(scrollX = TRUE)
        )
})
