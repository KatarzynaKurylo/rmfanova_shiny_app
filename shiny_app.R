library(shiny)
library(shinydashboard)
library(shinybusy) # for loading spinner
library(bsplus) #for info popover
library(plotly)
library(rmfanova)
library(data.table)

mean_fun_point <- function(x, is_legend, is_color, values = FALSE, type = "l", lty = 1, ...) {
  p <- ncol(x[[1]])
  l <- length(x)
  means <- matrix(NA, nrow = p, ncol = l)
  for (i in seq_len(l)) {
    means[, i] <- colMeans(x[[i]])
  }
  
  p <- plot_ly()
  for (i in 1:ncol(means)) {
    if (is_color) {
      p <- add_trace(p, x = 1:nrow(means), y = means[,i],
                     type = 'scatter', mode = 'lines',
                     name = paste("Group", i),
                     showlegend = is_legend)
    } else {
      p <- add_trace(p, x = 1:nrow(means), y = means[,i],
                     type = 'scatter', mode = 'lines',
                     name = paste("Group", i),
                     showlegend = is_legend,
                     line = list(color = "black"))  # Set line color to black
    }
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
                 type = 'scatter', mode = 'lines',
                 line = list(color = "black"))
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
                 type = 'scatter', mode = 'lines',
                 line = list(color = "black"))
  p <- layout(p, xaxis = list(title = "t"),
              title = "F(t)")
  return(p)
}

main_header <- div(class = "title",
                  h1('Functional repeated measures analysis of variance'),
                  tags$style(".title :is(h1){color: white; text-align: center; margin-top: 10px; font-size: 24px;}")
                  )

header <-  htmltools::tagQuery(dashboardHeader(title="Uploading Files"))

header <- header$
  addAttrs(style = "position: relative")$ 
  find(".navbar.navbar-static-top")$
  append(main_header)$ # inject our main header
  allTags()

ui <- dashboardPage(
  header,
  dashboardSidebar(
    # fileInput("file1", "Choose CSV File",
    #           multiple = FALSE,
    #           accept = c("text/csv",
    #                      "text/comma-separated-values,text/plain",
    #                      ".csv")),
    use_bs_popover(),
    fileInput("file1", "Choose CSV File",
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")) %>%
      shinyInput_label_embed(
        shiny_iconlink() %>%
          bs_embed_popover(
            title = "File Format", 
            content = "CSV File should have the sample number in the first column. Each of the other columns should indicate a discrete time point. The number of rows should be the number of samples multiplied by the number of observations.", 
            placement ="right"
          )
      ),
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
    tags$style(HTML('.popover-title {color:black;}
                               .popover-content {color:black;max-width: 400px;min-width: 200px; text-align: justify;}
                               .main-sidebar {z-index:auto;}
                                }')),
    tags$head(tags$style(HTML('
      .main-header .logo {
        font-size: 20px;
      }
      .author-text {
        text-align: right;
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
      ),
      tabPanel(
        title = "Data visualisation",
        br(),
        # fluidRow(
        #   column(4, textInput("x_axis", "Provide x axis name:")),
        #   column(4, textInput("y_axis", "Provide y axis name:")),
        #   column(2, checkboxInput("color", "Color", TRUE),style = "margin-top: 20px;"),
        #   column(2, checkboxInput("legend", "Legend", TRUE),style = "margin-top: 20px;")
        # ),
        textInput("x_axis", "Provide x axis name:",value="Time"),
        textInput("y_axis", "Provide y axis name:",value="Value"),
        checkboxInput("legend", "Legend", TRUE),
        checkboxInput("color", "Color", TRUE),
        uiOutput("input_df_plots")
      ),
      tabPanel(
        title = "Summary plots",
        br(),
        checkboxInput("mean_functions_legend", "Legend", TRUE),
        checkboxInput("mean_functions_color", "Color", TRUE),
        plotlyOutput("mean_functions"),
        plotlyOutput("ssa_statistics"),
        plotlyOutput("f_statistics")
      ),
      tabPanel(
        title = "Hypothesis testing",
        br(),
        numericInput("n_perm", "Number of permutation replicates:", value = 1000, min = 1, step = 1),
        numericInput("n_boot", "Number of bootstrap replicates:", value = 1000, min = 1, step = 1),
        checkboxInput("parallel", "Parallel computing", FALSE),
        checkboxInput("multi_gen", "Multiple generations", FALSE),
        fluidRow(
          uiOutput("test_stat_table")
        ),
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
      tags$style(
        HTML("
        .text-justified {
          text-align: justify;
        }
        .info-box {
          border: 1px solid #ccc;
          padding: 10px;
          margin-bottom: 10px;
        }
      ")
      ),
      br(),
      div(
        class = "info-box",
        p("Functional data analysis (FDA) is a branch of statistics which analyzes observations treated as functions, curves, or surfaces. To represent the data in such a way, one needs only to measure some variable over time or space, which is a scenario encountered in many fields. Then the discrete data observed at so-called design time points can be transformed into functional data. Such a representation allows us to avoid many problems of classical multivariate statistical methods, for example, the curse of dimensionality and missing data. ", class = "text-justified"),
        p("To compare the results for different samples, we thus consider functional repeated measures analysis of variance. For this purpose, a pointwise test statistic is constructed by adapting the classical test statistic for one-way repeated measures analysis of variance to the functional data framework. By integrating and taking the supremum of the pointwise test statistic, we create two global test statistics. Apart from verifying the general null hypothesis on the equality of mean functions corresponding to different objects, we also propose a simple method for post hoc analysis.", class = "text-justified"),
        p("Feel free to explore these resources for more information:"),
        a("Functional repeated measures analysis of variance and its application article", href = "https://arxiv.org/abs/2306.03883"),
        p("\n"),
        a("CRAN rmfanova package", href = "https://cran.r-project.org/web/packages/rmfanova/index.html"),
        p("\n"),
        p("Authors: Katarzyna Kurylo & Lukasz Smaga",class="author-text")
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
  
  output$data_table <- DT::renderDataTable({
    req(data())
    return(data()$df)
  })
  
  output$input_df_plots <- renderUI({
    req(data())
    df <- data()$df
    group_names <- unique(df[,1])
    is_legend<-input$legend
    is_color<-input$color
    x_axis<-input$x_axis
    y_axis<-input$y_axis
    plots <- lapply(1:length(group_names), function(i) {
      group_data <- df[df[,1] == group_names[i], -1]
      p <- plot_ly()
      for (j in 1:nrow(group_data)) {
        if (is_color) {
          p <- add_trace(p, x = seq(1, ncol(group_data)), y = as.numeric(group_data[j,]), 
                         type = 'scatter', mode = 'lines', name = paste("Observation", j), 
                         showlegend = is_legend)  #Use default colors
        } else {
          p <- add_trace(p, x = seq(1, ncol(group_data)), y = as.numeric(group_data[j,]), 
                         type = 'scatter', mode = 'lines', name = paste("Observation", j), 
                         showlegend = is_legend,
                         line = list(color = "black"))  # Set line color to black
        }
      }
      p <- layout(p, xaxis = list(title = x_axis), yaxis = list(title = y_axis), 
                  title = paste("Group", group_names[i]))
      return(p)
    })
    return(plots)
  })
  
  output$mean_functions <- renderPlotly({
    req(data())
    yy <- data()$matrix
    is_legend<-input$mean_functions_legend
    is_color<-input$mean_functions_color
    p<-mean_fun_point(yy, is_legend, is_color, values = FALSE, 
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
    n_perm<-input$n_perm
    n_boot<-input$n_boot
    parallel<-input$parallel
    multi_gen<-input$multi_gen
    df <- data()$df
    yy <- split(df[, -1], df[, 1]) #split by first column
    yy <- lapply(yy, as.matrix)
    res <- rmfanova(yy,
                    n_perm = n_perm,
                    n_boot = n_boot,
                    parallel = parallel,
                    #n_cores = NULL,
                    multi_gen = multi_gen)
    return(res)
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
