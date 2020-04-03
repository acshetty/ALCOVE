userinp <- reactive({input$user})

output$projectControl <- renderUI({

   projects <- project_assign$Project.ID[project_assign$User.ID == userinp()]
   memo <- project_assign$Pipeline.ID[project_assign$Pipeline.ID == projects]
   selectInput("project", label = h3("Choose Project"),
               choices = c(projects,"-"),
               selected = "-")
})

project.select <- reactive({input$project})

output$data <- renderUI({

   memo <- project_assign$Pipeline.ID[project_assign$User.ID == userinp() & project_assign$Project.ID == project.select()]

   selectInput("pipeline", label = h3("Choose Pipeline"), choices = c(memo,"-"), selected = "-")
})

data.select <- reactive({input$pipeline})

output$user_select <- renderMenu({
   selectInput("user_input", label = h5("Select User"),
               choices = c(unique(project_assign$User.ID),"-"),
               selected = "-")
})

user_inp_select <- reactive({input$user_input})

output$project_select <- renderMenu({
   projects <- project_assign$Project.ID[project_assign$User.ID == user_inp_select()]
   memo <- project_assign$Pipeline.ID[project_assign$Project.ID == projects]

   selectInput("project_input", label = h5("Select Project"),
               choices = c(projects,"-"),
               selected = "-")
})

project_input_select <- reactive({input$project_input})

output$data_select <- renderMenu({
   memo <- project_assign$Pipeline.ID[project_assign$User.ID == user_inp_select() & project_assign$Project.ID == project_input_select()]
   selectInput("data_select", label = h5("Choose Dataset"), choices = c(memo,"-"), selected = "-")
})

output$project_option_select <- renderUI({selectInput("select", label = h3("Select Project"),
                                                      choices = mixedsort(fastqc_dir_names))
})

base_path = reactive({normalizePath(file.path(paste0(project_assign$Repository.Root[project_assign$User.ID == user_inp_select() & project_assign$Project.ID == project_input_select() & project_assign$Pipeline.ID == data.select()]), paste0("output_repository")))})