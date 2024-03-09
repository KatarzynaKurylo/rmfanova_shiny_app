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
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select a file ----
      fileInput("file1", "Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Checkbox if file has header ----
      checkboxInput("header", "Header", TRUE),
      
      # Input: Select separator ----
      radioButtons("sep", "Separator",
                   choices = c(Comma = ",",
                               Semicolon = ";",
                               Tab = "\t"),
                   selected = ","),
      
      # Input: Select quotes ----
      radioButtons("quote", "Quote",
                   choices = c(None = "",
                               "Double Quote" = '"',
                               "Single Quote" = "'"),
                   selected = '"'),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Select number of rows to display ----
      radioButtons("disp", "Display",
                   choices = c(Head = "head",
                               All = "all"),
                   selected = "head")
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Data file ----
      plotlyOutput("mean_functions"),
      plotlyOutput("ssa_statistics"),
      plotlyOutput("f_statistics")#,
      #textOutput("rmfanova")
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
    groups <- unique(df[, 1])

    return(list(df = df, groups = groups))
  })
  
  means <-reactive({
    req(data())
    df <- data()$df
    groups<-data()$groups
    
    mean_gr <- data.frame(matrix(nrow = length(groups), ncol = ncol(df) - 1))
    
    for (i in seq_along(groups)) {
      group_data <- df[df[, 1] == groups[i], -1]  # Exclude group column
      mean_values <- colMeans(group_data)
      mean_gr[i, ] <- mean_values
      n<-nrow(df[df[, 1] == groups[i],])
    }
    #print(mean_gr)
    mean_all<-colMeans(mean_gr)
    
    return(list(mean_gr = mean_gr, mean_all = mean_all,groups=groups,n=n))
  })
  
  output$mean_functions <- renderPlotly({
    # req(means())
    # groups <- means()$groups
    # mean_gr<-means()$mean_gr
    # 
    # p <- plot_ly()
    # 
    # for (i in seq_along(groups)) {
    #   mean_values<-unlist(mean_gr[i, ])
    #   p <- add_trace(p, x = 1:length(mean_values), y = mean_values, 
    #                  type = 'scatter', mode = 'lines', 
    #                  name = paste("Group", groups[i]))
    # }
    # 
    # p <- layout(p, xaxis = list(title = "t"), yaxis = list(title = "FA"),
    #             title = "Mean functions by group")
    req(data())
    df <- data()$df
    yy <- split(df[, -1], df[, 1]) #split by first column
    yy <- lapply(yy,as.matrix)
    p<-mean_fun_point(yy, values = FALSE, 
                   col = 1:4, xlab = "t", ylab = "FA", xaxt = "n")
    return(p)
  })
  output$ssa_statistics<-renderPlotly({
    # req(means())
    # means_gr<-t(means()$mean_gr)
    # means_all<-means()$mean_all
    # n<-means()$n
    # ssa <- n * rowSums((means_gr - means_all)^2)
    # p <- plot_ly()
    # p <- add_trace(p, x = 1:length(ssa), y = ssa, 
    #                type = 'scatter', mode = 'lines')
    # p <- layout(p, xaxis = list(title = "t"), #yaxis = list(title = "FA"),
    #             title = "SSA(t)")
    # return(p)
    req(data())
    df <- data()$df
    yy <- split(df[, -1], df[, 1]) #split by first column
    yy <- lapply(yy,as.matrix)
    ssa_point(yy, xlab = "t", xaxt = "n")
    ssa<-ssa_point(yy, xlab = "t", xaxt = "n")
    return(ssa)
  })
  output$f_statistics<-renderPlotly({
    req(data())
    df <- data()$df
    yy <- split(df[, -1], df[, 1]) #split by first column
    yy <- lapply(yy,as.matrix)
    f<-f_point(yy, xlab = "t", xaxt = "n")
    return(f)
  })
  output$rmfanova<- renderText({
    req(data())
    df <- data()$df
    yy <- split(df[, -1], df[, 1]) #split by first column
    yy <- lapply(yy,as.matrix)
    res <- rmfanova(yy)
    return(summary(res, digits = 3))
  })
  
}

shinyApp(ui, server)
