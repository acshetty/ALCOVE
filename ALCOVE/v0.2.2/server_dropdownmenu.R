### rm(list=ls())

### Invoking requisite R libraries
library(shinydashboard, warn.conflicts = F)
library(shiny, warn.conflicts = F)

### Reading tab-delimited file containing the following:
### Column 1: User ID
### Column 2: Project ID
### Column 3: Repository Root
### Column 4: Pipeline ID
oProjects = read.delim("./user_projects/RNASEQ.pipelines.txt", sep="\t", header = TRUE, stringsAsFactor = FALSE)

### Generate 'Users' Drop Down list
output$Users = renderUI({
	selectInput(
		"username",
		label = h4("1. Select an User"),
		choices = c(unique(oProjects$User.ID),"-"),
		selected = "-"
	)
})

### Generate 'Projects' Drop Down list
rUser = reactive({input$username})
output$Projects = renderUI({
	aProjects = oProjects$Project.ID[oProjects$User.ID == rUser()]
	selectInput(
		"projectname",
		label = h4("2. Select a Project"),
		choices = c(aProjects,"-"),
		selected = "-"
	)
})

### Generate 'Pipelines' Drop Down list
rProject = reactive({input$projectname})
output$Pipelines = renderUI({
	aPipelines = oProjects$Pipeline.ID[oProjects$User.ID == rUser() & oProjects$Project.ID == rProject()]
	selectInput(
		"pipelineID",
		label = h4("3. Select a Pipeline"),
		choices = c(aPipelines,"-"),
		selected = "-"
	)
})

### Identify Repository Root
rPipelineID = reactive({input$pipelineID})
output$selectuser = renderText({paste("Username :", rUser())})
output$selectproject = renderText({paste("Project :", rProject())})
output$selectpipeline = renderText({paste("Pipeline ID :", rPipelineID())})

rRepoRoot = reactive({oProjects$Repository.Root[oProjects$User.ID == rUser() & oProjects$Project.ID == rProject() & oProjects$Pipeline.ID == rPipelineID()]})
output$repositoryroot = renderText({paste("Repository :", rRepoRoot())})

