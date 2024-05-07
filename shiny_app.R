library(shiny)
library(shinydashboard)
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
              title = "Mean functions by group")
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

ui <- dashboardPage(
  dashboardHeader(title = "Uploading Files"),
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
    tabsetPanel(
      tabPanel(
        title = "Data set",
        div(
          style = "overflow-x: auto;",
          DT::dataTableOutput("data_table")
        )
        #DT::dataTableOutput("data_table")
      ),
      tabPanel(
        title = "Data visualisation",
        uiOutput("input_df_plots")
      ),
      tabPanel(
        title = "Statistical plots",
        plotlyOutput("mean_functions"),
        plotlyOutput("ssa_statistics"),
        plotlyOutput("f_statistics")
      ),
      # tabPanel(
      #   title = "rmfanova Summary",
      #   verbatimTextOutput("rmfanova")
      # )
      tabPanel(
        title = "rmfanova Summary",
        fluidRow(
          valueBoxOutput("rmfanova_summary_l"), #number of samples
          valueBoxOutput("rmfanova_summary_n"), #number of observations
          valueBoxOutput("rmfanova_summary_p") #number of time points
        ),
        fluidRow(
          valueBoxOutput("rmfanova_summary_method"),
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
                       type = 'scatter', mode = 'lines', name = paste("Trajectory", j))
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

  output$rmfanova_summary_n <- renderValueBox({
    res <- rmfanova_result()
    valueBox(
      value = res$n,
      subtitle = "Number of Observations",
      color = "blue"
    )
  })
  
  output$rmfanova_summary_p <- renderValueBox({
    res <- rmfanova_result()
    valueBox(
      value = res$p,
      subtitle = "Number of Design Time Points",
      color = "blue"
    )
  })
  
  output$rmfanova_summary_l <- renderValueBox({
    res <- rmfanova_result()
    valueBox(
      value = res$l,
      subtitle = "Number of Samples",
      color = "blue"
    )
  })
  
  output$rmfanova_summary_method <- renderValueBox({
    res <- rmfanova_result()
    valueBox(
      value = res$method,
      subtitle = "Adjustment Method",
      color = "blue"
    )
  })
  
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
