library(shiny)
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
  #matplot(means, type = type, lty = lty, main = main, ...)
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
  # plot(ssa, type = type, ylab = ylab, main = main, ...)
  # if (values) {
  #   return(ssa)
  # }
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
  # plot(f_point, type = type, ylab = ylab, main = main, ...)
  # if (values) {
  #   return(f_point)
  # }
}


ui <- fluidPage(
  
  titlePanel("Uploading Files"),
  sidebarLayout(
    sidebarPanel(
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
      tags$hr(),
      radioButtons("disp", "Display",
                   choices = c(Head = "head",
                               All = "all"),
                   selected = "head")
      
    ),
    mainPanel(
      tabsetPanel(id = "tabs",
                  tabPanel(
                    title = "Data set",
                    DT::dataTableOutput("data_table")
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
                  tabPanel(
                    title = "rmfanova Summary",
                    # conditionalPanel(
                    #   condition = "!output.loaded",
                    #  # condition = "output.loaded == false",
                    #   #HTML('<img src="loading.gif" width="50" height="50" alt="loading...">')
                    #   tags$img(src = "loading.png", width = 50, height = 50, alt = "loading...")
                    # ),
                    verbatimTextOutput("rmfanova")
                  )
      )
    )
  )
)

server <- function(input, output) {
  #output$loaded <- reactiveVal(FALSE)
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
  output$data_table<-DT::renderDataTable({
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
  output$ssa_statistics<-renderPlotly({
    req(data())
    yy <- data()$matrix
    ssa<-ssa_point(yy, xlab = "t", xaxt = "n")
    return(ssa)
  })
  output$f_statistics<-renderPlotly({
    req(data())
    yy <- data()$matrix
    f<-f_point(yy, xlab = "t", xaxt = "n")
    return(f)
  })
  # output$rmfanova<- renderText({
  #   req(data())
  #   df <- data()$df
  #   yy <- split(df[, -1], df[, 1]) #split by first column
  #   yy <- lapply(yy,as.matrix)
  #   res <- rmfanova(yy)
  #   return(summary(res, digits = 3))
  # })
  output$rmfanova <- renderPrint({
      req(data())
      df <- data()$df
      yy <- split(df[, -1], df[, 1]) #split by first column
      yy <- lapply(yy,as.matrix)
      res <- rmfanova(yy)
      
    cat("Number of observations:", res$n, "\n")
    cat("Number of design time points:", res$p, "\n")
    cat("Number of samples:", res$l, "\n")
    cat("Adjustment method for pairwise comparison tests:", res$method, "\n")
    print(res$test_stat)
    print(res$p_values)
    print(res$p_values_pc)
    #output$loaded(TRUE)
  })
  
}

shinyApp(ui, server)
