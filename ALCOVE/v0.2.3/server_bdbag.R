### rm(list=ls())

### Invoking requisite R libraries
library(shinydashboard, warn.conflicts = F)
library(shiny, warn.conflicts = F)

### Reading tab-delimited file containing the following:
### Column 1: User ID
### Column 2: Project ID
### Column 3: Repository Root
### Column 4: Pipeline ID
#oProjects = read.delim("/local/projects/RNASEQ/pipeline_info.txt", sep="\t", header = TRUE, stringsAsFactor = FALSE)
#oProjects = oProjects[1:48,]
#colnames(oProjects) = c("User.ID","Project.ID","Repository.Root","Pipeline.ID","Sample.Info","Pipeline.Config")

### Upload BdBag
unzipPath <-"/local/scratch/"
bdbag_path<-NULL
makeReactiveBinding("bdbag_path")

observeEvent(input$bdbag,{
        bdbag_path <<- renderText({ input$bdbag$datapath })

observeEvent(input$load,{
        output$zipped <- renderTable({
        unzip(bdbag_path(), list= T)})
        output$bdbagPath = renderText({bdbag_path()})
})

observeEvent(input$unzip, {
        ldb<-unzip(bdbag_path(), exdir= unzipPath, unzip="unzip")
})

### Read Info File
oInfo = reactive({
    #oDAT = read.delim(rInfofile(), header=F, sep="\t", stringsAsFactor=F)
    oDAT = read.delim(input$infoF$datapath, header=F, sep="\t", stringsAsFactor=F)
    #parseFilePaths("/home", input$infoF)
    colnames(oDAT)[1] = "Sample.ID"
    colnames(oDAT)[2] = "Condition"
    return(oDAT)
    })

#### Some variable to check if everything works
#oPullInfo = reactive({
    #oDAT1 = subset(oAlignment(), select = c("Sample.ID","Percent.Mapped.Reads"))
 #   oDAT1 = subset(oAlignment(), select = c("Sample.ID","Total.Reads","Total.Mapped.Reads", "Percent.Mapped.Reads"))
  #  oDAT2 = subset(oInfo(), select = c("Sample.ID", "Condition"))
   # oDAT = merge(oDAT1, oDAT2, by = "Sample.ID")
    #colnames(oDAT)<-c("SampleID", "Total.Reads", "Total.Mapped", "Percent.Mapped")
   # oDAT$Unmapped = (oDAT$Total.Reads-oDAT$Total.Mapped)
    #oDAT = melt(oDAT, id.vars = c("Sample.ID", "Condition"))
    #return(oDAT)
   # })

#output$pname<-renderText({ if(input$load){ return(tools::file_path_sans_ext(input$bdbag$name))} else {return (tools::file_path_sans_ext(basename(input$lbdbag)))}})
#prjname<-renderText({if(input$load){ return(tools::file_path_sans_ext(input$bdbag$name))} else {return (tools::file_path_sans_ext(basename(input$lbdbag)))}})

output$pname<-renderText({tools::file_path_sans_ext(input$bdbag$name)})
prjname<-renderText({tools::file_path_sans_ext(basename(input$bdbag$name))})
})
##rReporoot<-renderText({ paste0(unzipPath,"/",prjname())})
#output$repositoryroot <- renderText({paste("Repository :", rRepoRoot())})
##repoRoot<- renderText({ paste0(unzipPath,"/",prjname())})
##output$reporoot <- renderText({ rReporoot() })
#output$rmode <- renderText({ if(input$load || input$lload){return(TRUE)} else {return(FALSE)}})


#})


observeEvent(input$lbdbag, {
        bdbag_path <<- renderText({ input$lbdbag })


observeEvent(input$load,{
        output$zipped <- renderTable({
        unzip(bdbag_path(), list= T)})
        output$bdbagPath = renderText({bdbag_path()})
})

observeEvent(input$unzip, {
        ldb<-unzip(bdbag_path(), exdir= unzipPath, unzip="unzip")
})
oInfo = reactive({
    #oDAT = read.delim(rInfofile(), header=F, sep="\t", stringsAsFactor=F)
    oDAT = read.delim(input$infoF$datapath, header=F, sep="\t", stringsAsFactor=F)
    #parseFilePaths("/home", input$infoF)
    colnames(oDAT)[1] = "Sample.ID"
    colnames(oDAT)[2] = "Condition"
    return(oDAT)
    })

observeEvent(input$infoF,{
        output$infoV <- renderTable({
	read.delim(input$infoF$datapath, header=F, sep="\t", stringsAsFactor=F)
        #unzip(bdbag_path(), list= T)

})
})

#output$pname<-renderText({ if(input$load){ return(tools::file_path_sans_ext(input$bdbag$name))} else {return (tools::file_path_sans_ext(basename(input$lbdbag)))}})
##output$pname<-renderText({tools::file_path_sans_ext(input$lbdbag)})
#prjname<-renderText({if(input$load){ return(tools::file_path_sans_ext(input$bdbag$name))} else {return (tools::file_path_sans_ext(basename(input$lbdbag)))}})
##prjname<-renderText({tools::file_path_sans_ext(basename(input$lbdbag))})
})
#output$pname<-renderText({ if(input$bdbag){ return(tools::file_path_sans_ext(input$bdbag$name))} else {return (tools::file_path_sans_ext(basename(input$lbdbag)))}})
#prjname<-renderText({ if(input$bdbag){ return(tools::file_path_sans_ext(input$bdbag$name))} else {return (tools::file_path_sans_ext(basename(input$lbdbag)))}})
output$pname<-renderText({tools::file_path_sans_ext(basename(bdbag_path()))})
prjname<-renderText({tools::file_path_sans_ext(basename(bdbag_path()))})
rReporoot<-renderText({ paste0(unzipPath,"/",prjname())})
output$reporoot <- renderText({ rReporoot() })
#output$rmode <- renderText({ if(input$load || input$lload){return(TRUE)} else {return(FALSE)}})
output$rmode <- renderText({ if(input$load){return(TRUE)} else {return(FALSE)}})

#})
