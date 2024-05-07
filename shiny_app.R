library(shiny)
library(shinydashboard)
library(shinybusy)
library(plotly)
library(rmfanova)
library(data.table)

mean_fun_point <- function(x, values = FALSE, type = "l", lty = 1, 
                           main = "Sample mean functions", ...) {
  p <- ncol(x[[1]])
  l <- length(x)
  means <- matrix(NA, nrow = p, ncol = l)
  for (i in seq_len(l)) {
    means[, i] <- colMeans(x[[i]])
  }
  
  p <- plot_ly()
  for (i in 1:ncol(means)) {
    p <- add_trace(p, x = 1:nrow(means), y = means[,i],
                   type = 'scatter', mode = 'lines',
                   name = paste("Group", i))
  }
  
  p <- layout(p, xaxis = list(title = "t"), yaxis = list(title = "FA"),
              title = "Sample mean functions by group")
  return(p)
}

ssa_point <- function(x, values = FALSE, 
                      type = "l", ylab = "", main = "SSA(t)", ...) {
  n <- nrow(x[[1]])
  means_gr <- sapply(x, colMeans)
  means_all <- rowMeans(means_gr)
  ssa <- n * rowSums((means_gr - means_all)^2)
  p <- plot_ly()
  p <- add_trace(p, x = 1:length(ssa), y = ssa,
                 type = 'scatter', mode = 'lines')
  p <- layout(p, xaxis = list(title = "t"),
              title = "SSA(t)")
  return(p)
}

f_point <- function(x, values = FALSE, 
                    type = "l", ylab = "", main = "F(t)", ...) {
  n <- nrow(x[[1]])
  p <- ncol(x[[1]])
  k <- length(x)
  means_gr <- sapply(x, colMeans)
  means_all <- rowMeans(means_gr)
  means_sub <- matrix(0, nrow = p, ncol = n)
  for (ii_s in seq_len(n)) {
    means_sub_temp <- 0
    for (jj_s in seq_len(k)) {
      means_sub_temp <- means_sub_temp + x[[jj_s]][ii_s, ]
    }
    means_sub[, ii_s] <- means_sub_temp / k
  }
  SSA <- rowSums(n * (means_gr - means_all)^2)
  SSS <- rowSums(k * (means_sub - means_all)^2)
  SST <- rowSums(sapply(lapply(x, function(x) (t(x) - means_all)^2), rowSums))
  SSE <- SST - SSA - SSS
  f_point <- (SSA / (k - 1)) / (SSE / ((n - 1) * (k - 1)))
  f_point <- f_point[is.finite(f_point)]
  f_point <- ifelse(f_point < .Machine$double.eps, 0, f_point)
  p <- plot_ly()
  p <- add_trace(p, x = 1:length(f_point), y = f_point,
                 type = 'scatter', mode = 'lines')
  p <- layout(p, xaxis = list(title = "t"),
              title = "F(t)")
  return(p)
}

header_img <- div(class = "my-title",
                  h1('Functional repeated measures analysis of variance'),
                  tags$style(".my-title :is(h1){color: white; text-align: center; margin-top: 10px; font-size: 24px;}")
                  )
# dashboard_header <- dashboardHeader(
#   title = h4(span("Loading Files", style = "margin-top: 50px;"))
# )
header <-  htmltools::tagQuery(dashboardHeader(title="Loading Files"))

# header$title <- div(
#   header$title,
#   tags$style("margin-top: 15px;")
# )
header <- header$
  addAttrs(style = "position: relative")$ # add some styles to the header 
  find(".navbar.navbar-static-top")$ # find the header right side
  append(header_img)$ # inject our img
  allTags()

ui <- dashboardPage(
  # dashboardHeader(titleWidth='100%',
  #                 title = div(
  #                   column(12, class="title-box", 
  #                          tags$h1(style='margin-top:10px;', 'Functional repeated measures analysis of variance')
  #                   )
  #                 ),
  #                 dropdownMenuOutput("helpMenu")),
  header,
  dashboardSidebar(
    fileInput("file1", "Choose CSV File",
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
    tags$hr(),
    checkboxInput("header", "Header", TRUE),
    radioButtons("sep", "Separator",
                 choices = c(Comma = ",",
                             Semicolon = ";",
                             Tab = "\t"),
                 selected = ","),
    radioButtons("quote", "Quote",
                 choices = c(None = "",
                             "Double Quote" = '"',
                             "Single Quote" = "'"),
                 selected = '"'),
    tags$hr()
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
      .main-header .logo {
        font-size: 20px;
      }
    '))),
    add_busy_spinner(spin = "fading-circle", 
                     margins = c(5, 5),
                     height = "40px",
                     width = "40px",
                     color = "white"
    ),
    tabsetPanel(
      tabPanel(
        title = "Info",
        uiOutput("informations")
      ),
      tabPanel(
        title = "Data set",
        br(), #instead of br() I can use fluidRow(...,style = "padding-top:20px")
        fluidRow(
          valueBoxOutput("dataset_l"), #number of samples
          valueBoxOutput("dataset_n"), #number of observations
          valueBoxOutput("dataset_p") #number of time points
        ),
        div(
          style = "overflow-x: auto;",
          DT::dataTableOutput("data_table")
        )
        #DT::dataTableOutput("data_table")
      ),
      tabPanel(
        title = "Data visualisation",
        br(),
        uiOutput("input_df_plots")
      ),
      tabPanel(
        title = "Summary plots",
        br(),
        plotlyOutput("mean_functions"),
        plotlyOutput("ssa_statistics"),
        plotlyOutput("f_statistics")
      ),
      # tabPanel(
      #   title = "rmfanova Summary",
      #   verbatimTextOutput("rmfanova")
      # )
      tabPanel(
        title = "Hypothesis testing",
        br(),
        # fluidRow(
        #   valueBoxOutput("rmfanova_summary_l"), #number of samples
        #   valueBoxOutput("rmfanova_summary_n"), #number of observations
        #   valueBoxOutput("rmfanova_summary_p") #number of time points
        # ),
        fluidRow(
          #valueBoxOutput("rmfanova_summary_method"),
          uiOutput("test_stat_table")
        ),
        #fluidRow(uiOutput("test_stat_table"))
        # fluidRow(
        #   DT::dataTableOutput("test_stat")
        # ),
        fluidRow(
          uiOutput("p_values_table")
        ),
        fluidRow(
          uiOutput("p_values_pc_table")
        )
      )
    )
  )
)

server <- function(input, output) {
  output$informations<- renderUI({
    tagList(
      p("Here are some useful links related to rmfanova:"),
      div(
        #img(src = "Rlogo.png"), # Add your image path
        p("Feel free to explore these resources for more information:"),
        a("CRAN rmfanova package", href = "https://cran.r-project.org/web/packages/rmfanova/index.html"),
        p("\n"),
        a("Functional repeated measures analysis of variance and its application article", href = "https://arxiv.org/abs/2306.03883"),
      )
    )
  })
  data <- reactive({
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    df[, -1] <- sapply(df[, -1], as.numeric)
    
    splited_df <- split(df[, -1], df[, 1]) #split by first column
    matrix <- lapply(splited_df, as.matrix)
    
    return(list(df = df, matrix = matrix))
  })
  
  output$data_table <- DT::renderDataTable({
    req(data())
    return(data()$df)
  })
  
  output$input_df_plots <- renderUI({
    req(data())
    df <- data()$df
    group_names <- unique(df[,1])
    
    plots <- lapply(1:length(group_names), function(i) {
      group_data <- df[df[,1] == group_names[i], -1]
      p <- plot_ly()
      for (j in 1:nrow(group_data)) {
        p <- add_trace(p, x = seq(1, ncol(group_data)), y = as.numeric(group_data[j,]), 
                       type = 'scatter', mode = 'lines', name = paste("Observation", j))
      }
      p <- layout(p, xaxis = list(title = "X-axis"), yaxis = list(title = "Y-axis"), 
                  title = paste("Group", group_names[i]))
      return(p)
    })
    return(plots)
  })
  
  output$mean_functions <- renderPlotly({
    req(data())
    yy <- data()$matrix
    p<-mean_fun_point(yy, values = FALSE, 
                      col = 1:4, xlab = "t", ylab = "FA", xaxt = "n")
    return(p)
  })
  
  output$ssa_statistics <- renderPlotly({
    req(data())
    yy <- data()$matrix
    ssa <- ssa_point(yy, xlab = "t", xaxt = "n")
    return(ssa)
  })
  
  output$f_statistics <- renderPlotly({
    req(data())
    yy <- data()$matrix
    f <- f_point(yy, xlab = "t", xaxt = "n")
    return(f)
  })
  
  rmfanova_result <- reactive({
    req(data())
    df <- data()$df
    yy <- split(df[, -1], df[, 1]) #split by first column
    yy <- lapply(yy, as.matrix)
    res <- rmfanova(yy)
    return(res)
  })
  
  output$dataset_n <- renderValueBox({
    req(data())
    x <- data()$matrix
    valueBox(
      value = nrow(x[[1]]),
      subtitle = "Number of Observations",
      color = "blue"
    )
  })
  
  output$dataset_p <- renderValueBox({
    req(data())
    x <- data()$matrix
    valueBox(
      value = ncol(x[[1]]),
      subtitle = "Number of Design Time Points",
      color = "blue"
    )
  })
  
  output$dataset_l <- renderValueBox({
    req(data())
    x <- data()$matrix
    valueBox(
      value = length(x),
      subtitle = "Number of Samples",
      color = "blue"
    )
  })
  
  # output$rmfanova_summary_n <- renderValueBox({
  #   res <- rmfanova_result()
  #   valueBox(
  #     value = res$n,
  #     subtitle = "Number of Observations",
  #     color = "blue"
  #   )
  # })
  # 
  # output$rmfanova_summary_p <- renderValueBox({
  #   res <- rmfanova_result()
  #   valueBox(
  #     value = res$p,
  #     subtitle = "Number of Design Time Points",
  #     color = "blue"
  #   )
  # })
  # 
  # output$rmfanova_summary_l <- renderValueBox({
  #   res <- rmfanova_result()
  #   valueBox(
  #     value = res$l,
  #     subtitle = "Number of Samples",
  #     color = "blue"
  #   )
  # })
  
  # output$rmfanova_summary_method <- renderValueBox({
  #   res <- rmfanova_result()
  #   valueBox(
  #     value = res$method,
  #     subtitle = "Adjustment Method",
  #     color = "blue"
  #   )
  # })
  
  output$test_stat <- DT::renderDataTable({
    res <- rmfanova_result()
    df <- as.data.frame(res$test_stat)
    return(DT::datatable(df, options = list(dom = 't', pageLength = 5), rownames = FALSE))
  })
  
  output$test_stat_table <- renderUI({
    res<- rmfanova_result()
    box(
      title = "Overall test statistics",
      status = "primary",
      solidHeader = TRUE,
      DT::dataTableOutput("test_stat")
    )
  })
  
  output$p_values <- DT::renderDataTable({
    res <- rmfanova_result()
    df <- as.data.frame(res$p_values)
    return(DT::datatable(df, options = list(dom = 't', pageLength = 5), rownames = FALSE))
  })
  
  output$p_values_table <- renderUI({
    res<- rmfanova_result()
    box(
      title = "Overall p-values",
      status = "primary",
      solidHeader = TRUE,
      width = "auto",
      DT::dataTableOutput("p_values")
    )
  })
  
  output$p_values_pc <- DT::renderDataTable({
    res <- rmfanova_result()
    df <- as.data.frame(res$p_values_pc)
    return(DT::datatable(df, options = list(dom = 't', pageLength = 5)))
  })
  
  output$p_values_pc_table <- renderUI({
    res<- rmfanova_result()
    box(
      title = "Pairwise comparison p-values",
      status = "primary",
      solidHeader = TRUE,
      width = "auto",
      DT::dataTableOutput("p_values_pc")
    )
  })
}

shinyApp(ui, server)
