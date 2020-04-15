rm(list=ls())

### Invoking requisite R libraries
library(shinydashboard, warn.conflicts = F)
library(shiny, warn.conflicts = F)

### Reading tab-delimited file containing the following:
### Column 1: User ID
### Column 2: Project ID
### Column 3: Repository Root
### Column 4: Pipeline ID
#oProjects = read.delim("./user_projects/RNASEQ.pipelines.txt", sep="\t", header = TRUE, stringsAsFactor = FALSE)
oProjects = read.delim("/local/projects/RNASEQ/pipeline_info.txt", sep="\t", header = FALSE, stringsAsFactor = FALSE)
#colnames(oProjects)<- c("User.ID", "Project.ID", "Repository.Root", "Pipeline.ID")

### Setting HTML CSS
dbCSS = list(
					tags$head(
						tags$style(
							HTML("
								.multicol {
									height: 10%;
									-webkit-column-count: 4; /* Chrome, Safari, Opera */
									-moz-column-count: 4;    /* Firefox */
									column-count: 4;
									-webkit-column-fill: balance;
									-moz-column-fill: balance;
									column-fill: balance;
								}
								div.checkbox {margin-top: 0px;}
							")
						)
					)
				)

### Dashboard Header
dbHeader = dashboardHeader(title = "Alcove")

### Dashboard Sidebar
dbSidebar = dashboardSidebar(
					uiOutput("Users"),
					uiOutput("Projects"),
					uiOutput("Pipelines"),
					
                    sidebarMenu(
						menuItem("RNASeq Pipeline", tabName = "rnaseq", icon = icon("th-list")),
                        menuItem("Upload BDBAG", tabName = "BDBag", icon = icon("cog", lib = "glyphicon")),
						menuItem("FastQC", tabName = "fastqc", icon = icon("cog", lib = "glyphicon")),
						menuItem("Alignment", tabName = "alignment", icon = icon("cog", lib = "glyphicon")),
						menuItem("GeneExpression", tabName = "geneexp", icon = icon("cog", lib = "glyphicon")),
						menuItem("Differential Gene Expression", tabName = "diffgeneexp", icon = icon("cog", lib = "glyphicon"))
					)
				)

### dashboard Body RNASEQ Tab
#ergatisTab = tabItem(tabName = "ergatis",
#                    uiOutput("Users"),
#                    uiOutput("Projects"),
#                    uiOutput("Pipelines"),
                    #sidebarMenu(
                    #    menuItem("RNASeq Pipeline", tabName = "rnaseq", icon = icon("th-list")),
                    #    menuItem("FastQC", tabName = "fastqc", icon = icon("cog", lib = "glyphicon")),
                    #    menuItem("Alignment", tabName = "alignment", icon = icon("cog", lib = "glyphicon")),
                    #    menuItem("GeneExpression", tabName = "geneexp", icon = icon("cog", lib = "glyphicon")),
                    #    menuItem("Differential Gene Expression", tabName = "diffgeneexp", icon = icon("cog", lib = "glyphicon"))
                    #)
#                )

rnaseqTab = tabItem(tabName = "rnaseq",
						h2("Alcove"),
						h3("Transcriptomics Pipeline Results Exploration Tool"),
						helpText("BETA VERSION STATEMENT:"),
						helpText("This data exploration tool is intended for use by the Institute for Genome Sciences. ",
								"The data presented here are outputs from the in-house transcriptomics pipeline. ",
								"You may select your username, project and pipeline ID from the dropdown menu in the side bar to explore your results. ",
								"Please send questions, comments, suggestions for improvements, and error reports via email to ",
								"Amol Carl Shetty (ashetty@som.umaryland.edu). "
								),
						br(),
						br(),
						helpText("Developers:"),
						helpText("Amol Carl Shetty, Karthik Jallawaram, Ankita Parihar, Apaala Chatterjee"),
						br(),
						br(),
						helpText("Selected Repository:"),
						verbatimTextOutput("selectuser"),
						verbatimTextOutput("selectproject"),
						verbatimTextOutput("selectpipeline"),
						verbatimTextOutput("repositoryroot"),
						verbatimTextOutput("sampleinfo")
					)

uploadbdbagTab = tabItem(tabName = "BDBag",
                        textInput("lbdbag", "path/to/sample.zip", ""),
                        verbatimTextOutput("bdbagPath"),
                        fileInput("bdbag", "OR upload Zip file (<5mb)", accept = ".zip"),
                        actionButton("load", "View Uploaded Archive"),
                        actionButton("unzip", "Unzip Uploaded Archive"),
                        fileInput("infoF", "Upload Info File"),
                        #shinyFileButton("infoF", "Path/To/InfoFile", "InfoFile", F),
			tableOutput("filedf"),
                        tableOutput("zipped"),
                        verbatimTextOutput("pname"),
                        verbatimTextOutput("rmode"),
                        verbatimTextOutput("reporoot")
            )

fastqcTab =  tabItem(tabName = "fastqc",
                        helpText("FastQC Directory:"),
                        verbatimTextOutput("fastqcdir"),
                        helpText("These graphs were created using the quality control software FastQC ",
							"(Documentation: https://www.bioinformatics.babraham.ac.uk/projects/ fastqc/)."),	
                        fluidRow(
							column(6,uiOutput("sample_select")),
							column(6,uiOutput("image_select"))
						),
						fluidRow(
						   column(6,imageOutput("set1Img", width = "100%", height = "60%")),
						   column(6,imageOutput("set2Img", width = "100%", height = "60%"))
						),
						fluidRow(
						   column(6,downloadButton('downloadPlot1', 'Download Plot1')),
						   column(6,downloadButton('downloadPlot2', 'Download Plot2'))
						),
						fluidRow(
                                                                column(12, DT::dataTableOutput("static_table", width = "100%"))
                                                        )
                    )

alignmentTab = tabItem(tabName = "alignment",
                            helpText("Alignment Summary:"),
							verbatimTextOutput("alignsummary"),
						    verbatimTextOutput("summary"),
							fluidRow(
								column(12, plotOutput("AlignmentGraph1", width = "100%"))
							),
							fluidRow(
								column(6, downloadButton('downloadPlot3', 'Download Plot3'))
							),
							
							HTML("<br><br>"),
							
							fluidRow(
								column(12, plotOutput("AlignmentGraph2", width = "100%"))
							),
							fluidRow(
								column(6, downloadButton('downloadPlot4', 'Download Plot4'))
							),
                            fluidRow(
                                column(12, plotOutput("AlignInfoGraph1", width = "100%"))
                            ),
                            fluidRow(
                                column(6, downloadButton('downloadPlot5_1', 'Download Plot5'))
                            ),
							
							HTML("<br><br>"),
							fluidRow(
                            column(12, plotOutput("AlignInfoGraph2", width = "100%"))
                            ),
                            fluidRow(
                                column(6, downloadButton('downloadPlot7_1', 'Download Plot6'))
                            ),

							HTML("<br><br>"),
							fluidRow(
                                column(12, plotOutput("AlignComboGraph1", width = "100%"))
                            ),
                            fluidRow(
                                column(6, downloadButton('downloadPlot6_1', 'Download Plot7'))
                            ),

							HTML("<br><br>"),
							fluidRow(
								column(12, tags$div(align = 'left', class = 'multicol', uiOutput("column_select")))
							),
							
							HTML("<br><br>"),
							
							fluidRow(
								column(12, DT::dataTableOutput("alignment_table", width = "100%"))
							),
							fluidRow(
								column(6, downloadButton('downloadTable1', 'Download Table1'))
							)
							#fluidRow(
                                                        #        column(12, DT::dataTableOutput("info_table", width = "100%"))
                                                        #),
							#fluidRow(
                                                         #       column(12, DT::dataTableOutput("info_table1", width = "100%"))
                                                        #),

							#fluidRow(
                                                        #        column(12, plotOutput("AlignInfoGraph2", width = "100%"))
                                                        #)
						)

geneexpTab = tabItem(tabName = "geneexp",
                            helpText("Gene Expression Set:"),
							verbatimTextOutput("expression"),
							verbatimTextOutput("dimensions"),
							
							helpText("Normalized Expression Set:"),
							verbatimTextOutput("normalized"),
							fluidRow(
								column(12,
									tags$div(align = 'left', 
										radioButtons(
											"normalization",
											h4("Gene Expression Values"),
											choices = c("Raw Data","CPM Data"),
											selected = "CPM Data",
											inline = TRUE
										)
									)
								)
							),
							
							conditionalPanel("input.normalization == 'CPM Data'",
								fluidRow(
									column(12, tags$div(align = 'left', uiOutput("CPM_Slider")))
								),
								fluidRow(
									column(12, tags$div(align = 'left', actionButton("submit", "Apply Filter")))
								),
								
								helpText("Filtered Expression Set:"),
								verbatimTextOutput("filtered"),
								
                                fluidRow(
                                    column(12, plotOutput("DensityPlotInfo"))
                                ), 
								fluidRow(
									column(6, downloadButton('downloadPlot5', 'Download Plot8'))
								),
								HTML("<br><br>")
							),
							
							fluidRow(
								column(12, plotOutput("BoxPlot2"))
							),
							fluidRow(
							   column(6, downloadButton('downloadPlot6', 'Download Plot9'))
							),
							
							HTML("<br><br>"),
                            fluidRow(
                                column(12, plotOutput("PCAPlotInfo"))
                            ),
							fluidRow(
							   column(6,downloadButton('downloadPlot7', 'Download Plot10'))
							),
							
							HTML("<br><br>"),
	
							fluidRow(
								column(12, DT::dataTableOutput("gene_expression", width = "100%"))
							),
							fluidRow(
                                column(12, DT::dataTableOutput("PCA_Input", width = "100%"))
                            ),
							HTML("<br><br>"),
							fluidRow(
								column(6, downloadButton('downloadTable2', 'Download Table2'))
							),
                                                       fluidRow(
                                                                column(12, DT::dataTableOutput("cpm_table", width = "100%"))
                                                        )
						)

diffgeneexpTab = tabItem(tabName = "diffgeneexp",
                            helpText("Differential Gene Expression:"),
							verbatimTextOutput("degeneset"),
							
							tabsetPanel(
								type = "tabs",
								tabPanel(
									title = "Within Comparisons",
									
									column(12, tags$div(align = 'left', uiOutput("Select_Comparison"))),
									
									helpText("Differential Gene Expression Set:"),
									verbatimTextOutput("degenedir"),
									verbatimTextOutput("degenefile"),
									verbatimTextOutput("degenesummary"),
									
									helpText("Normalized Expression Set:"),
									verbatimTextOutput("countf"),
                                    verbatimTextOutput("normcount"), 
									
									helpText("Differential Gene Expression Filters:"),
									fluidRow(
										column(6, sliderInput("FDR", 'Select FDR cut-off', value = 0.05, min = 0, max = 1)),
										column(6, sliderInput("PVAL", 'Select P-value cut-off', value = 0.05, min = 0, max = 1))
									),
									fluidRow(
										column(6, sliderInput("RCP", 'Select Read Count Percentile cut-off',value = 0.1 , min = 0, max = 1)),
										column(6, uiOutput("LFC_Slider"))
									),
									fluidRow(
										column(12, tags$div(align = 'left', actionButton("filter", "Apply Filters")))
									),
									
									helpText("Differential Gene Expression Filtered Set:"),
									verbatimTextOutput("degenefiltered"), 
									
									helpText("Normalized Expression Filtered Set:"),
									verbatimTextOutput("normcountfiltered"), 
									
									fluidRow(
										column(12, plotOutput("MAPlot1"))
									),
									fluidRow(
									   column(6, downloadButton('downloadPlot8', 'Download Plot11'))
									), 
									
									HTML("<br><br>"),
									
									fluidRow(
										column(12, plotOutput("HeatMap2"))
									),
									fluidRow(
									   column(6, downloadButton('downloadPlot9', 'Download Plot12'))
									), 
									
									HTML("<br><br>"),
									
									fluidRow(
										column(12, DT::dataTableOutput("diff_gene_expression", width = "100%"))
									),
									fluidRow(
										column(6, downloadButton('downloadTable3', 'Download Table3'))
									)
								),
								tabPanel(
									title = "Between Comparisons",
									
									column(12, tags$div(align = 'left', uiOutput("Select_Multiple_Comparisons"))),
									
									helpText("Differential Gene Expression Sets:"),
									verbatimTextOutput("numsets"),
									verbatimTextOutput("degenelist"),
									verbatimTextOutput("overlapset"),
									
									helpText("Venn Diagram:"),
									
									fluidRow(
										column(12, plotOutput("Venn", height = "600px"))
									),
									fluidRow(
									   column(6, downloadButton('downloadPlot10', 'Download Plot13'))
									), 
									
									HTML("<br><br>"),
									
									helpText("Quadrant Plot:"),
									
									fluidRow(
										column(6, uiOutput("Comparison1")),
										column(6, uiOutput("Comparison2"))
									),
									fluidRow(
										column(12, plotOutput("Quadrant"))
									),
									fluidRow(
									   column(6, downloadButton('downloadPlot11', 'Download Plot14'))
									), 
									
									HTML("<br><br>")
								)
							)
						)

### Dashboard Body
dbBody = dashboardBody(
					dbCSS,
					tabItems(
                        uploadbdbagTab,
						rnaseqTab,
						fastqcTab,
						alignmentTab,
						geneexpTab,
						diffgeneexpTab
					)
				)

dashboardPage(dbHeader, dbSidebar, dbBody)

#dashboardPage(
#
#   dashboardSidebar(
#      # selectInput("user", label = h3("Select User"),
#      #             choices = c(unique(project_assign$userName),"-"),
#      #             selected = "-"),
#      #
#      # uiOutput("projectControl"),
#      # uiOutput("data"),
#
#      sidebarMenu(
#         menuItem("transcriptomics", tabName = "transcriptomics", icon = icon("dashboard")),
#
#         menuItem("Login", tabName = "login", icon = icon("th"),
#                  sidebarMenuOutput("user_select"),
#                  sidebarMenuOutput("project_select"),
#                  sidebarMenuOutput("data_select")
#         ),
#
#
#         menuItem("RNASEQ", tabName = "rnaseq", icon = icon("th"),
#                  menuSubItem("FastQC", tabName = "fastqc"),
#                  menuSubItem("Alignment", tabName = "alignment"),
#                  menuSubItem("GeneExpression", tabName = "geneexp"),
#                  menuSubItem("Differential Gene Expression", tabName = "diffgeneexp")
#                  #menuSubItem("Isoform Expression", tabName = "isoformexp"),
#                  #menuSubItem("Differential Isoform Expression", tabName = "diffisoformexp")
#         )
#      )
#   ),
#
#   dashboardBody(tweaks,
#                 tabItems(
#                    tabItem(tabName = "transcriptomics",
#                            tabItem(tabName = "dashboard",
#                                    selectInput("user", label = h3("Select User"),
#                                                choices = c(unique(project_assign$User.ID),"-"),
#                                                selected = "-"),
#
#                                    uiOutput("projectControl"),
#                                    uiOutput("data")
#                            )),
#
#                    tabItem(tabName = "fastqc",
#                            fluidRow(
#                               column(6,uiOutput("project_option_select")
#                               ),
#                               column(6,uiOutput("imageSelect"))
#                            ),
#
#
#
#                            fluidRow(
#                               column(6,imageOutput("set1Img", width = "100%", height = "60%")),
#                               column(6,imageOutput("set2Img", width = "100%", height = "60%"))),
#                            fluidRow(
#                               column(6,downloadButton('downloadPlot1', 'Download Plot1')),
#                               column(6,downloadButton('downloadPlot2', 'Download Plot2')))
#                    ),
#
#
#                    tabItem(tabName = "alignment",
#                            fluidRow(column(12,tags$div(align = 'left', class = 'multicol', uiOutput("selection")))),
#                            fluidRow(textOutput("Alignment_path")),
#
#                            HTML("<br><br>"),
#                            fluidRow(column(12, plotOutput("readGraph", width = "100%")),downloadButton('downloadPlot3', 'Download Plot')),
#                            HTML("<br><br>"),
#                            #fluidRow(column(2, numericInput("groups", "Number of groups",value = 1 ,min = 1, max = nrow(alignment), step = 1)),
#                            #column(6,uiOutput("group_selector"))),
#
#                            #fluidRow(column(3,uiOutput("row_selector")),
#                            #column(1,actionButton("select1", "Select Rows")),
#                            #column(1,actionButton("clear1", "Clear Rows"))),
#
#                            #HTML("<br><br>"),
#                            fluidRow(column(width = 12, DT::dataTableOutput("mytable", width = "100%"))),
#                            fluidRow(column(6,downloadButton('downloadTable_alignment', 'Download Table')))
#                    ),
#
#                    tabItem(tabName = "geneexp",
#                            radioButtons("normalize", h3("Dataset to be used for analysis"), choices = c("Raw Data","CPM Data"),selected = "CPM Data", inline = TRUE),
#                            uiOutput("CPM_Slider"),
#                            conditionalPanel("input.normalize == 'CPM Data'",fluidRow(column(12, plotOutput("densityGene"))), fluidRow(
#                               column(6,downloadButton('downloadPlot_densityPlot', 'Download Plot'))),
#                               HTML("<br><br>")),
#                            fluidRow(column(12, plotOutput("boxGene"))),
#                            fluidRow(
#                               column(6,downloadButton('downloadPlot_boxPlot', 'Download Plot'))),
#                            HTML("<br><br>"),
#                            fluidRow(column(12, plotOutput("pcaPlot_geneExp"))),
#                            fluidRow(
#                               column(6,downloadButton('downloadPlot_pcaPlot', 'Download Plot'))),
#                            HTML("<br><br>"),
#
#                            fluidRow(DT::dataTableOutput("geneExp"), width = "100%"),
#                            fluidRow(column(6,downloadButton('downloadTable_geneExp', 'Download Table')))
#
#
#                    ),
#
#                    tabItem(tabName = "diffgeneexp",
#                            tabsetPanel(type = "tabs",
#                                        tabPanel(title = "Heatmap",
#                                                 uiOutput("Select_Comparison"),
#                                                 fluidRow(
#                                                    column(6 , sliderInput("FDR_slider", 'FDR Slider',value = 0.05, min = 0, max = 1)),
#                                                    column(6, sliderInput("Prob_slider", 'P-value Slider',value = 0.05 , min = 0, max = 1))
#                                                 ),
#                                                 fluidRow(
#                                                    column(6, tags$style(type="text/css",
#                                                                         ".shiny-output-error { visibility: hidden; }",
#                                                                         ".shiny-output-error:before { visibility: hidden; }"), uiOutput("LFC")),
#                                                    column(6, tags$style(type="text/css",
#                                                                         ".shiny-output-error { visibility: hidden; }",
#                                                                         ".shiny-output-error:before { visibility: hidden; }"), uiOutput("RC_Slider"))
#                                                 ),
#                                                 fluidRow(
#
#                                                    column(8,textOutput("Heatmap_dim"))),
#
#                                                 HTML("<br><br>"),
#
#                                                 fluidRow(column(12,plotOutput("geneExpHeatmap1"))),
#                                                    column(8,textOutput("Heatmap_dim"))),
#
#                                                 HTML("<br><br>"),
#
#                                                 fluidRow(column(12,plotOutput("geneExpHeatmap1"))),
#                                                 fluidRow(column(6,downloadButton('downloadPlot_geneExpHeatMap1', 'Download Plot'))
#                                                 )
#                                                 ,
#                                                 fluidRow(column(12,plotOutput("geneExpHeatmap2"))),
#                                                 br(),
#                                                 fluidRow(column(12, plotOutput("diffGeneExp_scatterPlot"))),
#
#                                                br(),
#
#                                                 fluidRow(
#                                                    DT::dataTableOutput("diffgeneExp"), width = "100%", height = "100%"),
#                                                 fluidRow(column(6,downloadButton('downloadTable_diffgeneExp', 'Download Table')))
#                                        ),
#                                        tabPanel("Quad Plot",
#                                                 uiOutput("Quad_plot_comparison")
#                                                 )
#                            )
#                    )
#                 )
#   )
#
#)
