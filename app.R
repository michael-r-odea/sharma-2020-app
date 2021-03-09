#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(DT)

sharma <- readRDS("sharma_slim.rds")
gene.choices <- rownames(sharma)

all.groups <- readRDS("all_markers.rds")
p5m.p0m <- readRDS("p5m_p0m.rds")
p5m.p0u <- readRDS("p5m_p0u.rds")
p5m.p5u <- readRDS("p5m_p5u.rds")
p5u.p0m <- readRDS("p5u_p0m.rds")
p5u.p0u <- readRDS("p5u_p0u.rds")
p0m.p0u <-readRDS("p0m_p0u.rds")


# Define UI for application
ui <- fluidPage(
    titlePanel("Sharma et al. 2020: P0 & P5 data"),
    sidebarLayout(
    sidebarPanel(
        selectInput("task", "What would you like to do?", 
                    c("Explore a candidate gene", "Explore DEGs between groups")),
            
        conditionalPanel(condition = 'input.task == "Explore a candidate gene"',
                             selectInput("gene", "Search gene of interest:", gene.choices, multiple = TRUE)),
        
        conditionalPanel(condition = 'input.task == "Explore a candidate gene"',
                         checkboxGroupInput("comp1", "Choose groups to compare expression between (select 2 or 4):", 
                          c("P5 - Myelinated", "P5 - Unmyelinated", "P0 - Myelinated", "P0 - Unmyelinated"))),
            
        conditionalPanel(condition = 'input.task == "Explore DEGs between groups"', 
                           checkboxGroupInput("comp2", "Choose the groups to compare (select 2 or 4):", 
                            c("P5 - Myelinated", "P5 - Unmyelinated", "P0 - Myelinated", "P0 - Unmyelinated"))),

    
        conditionalPanel(condition = 'input.task == "Explore DEGs between groups"',
              sliderInput("logfc", "Log2 Fold-Change Threshold:", min = 0, max = 1, step = 0.05, value = 0.10, round = -2))),

        mainPanel(
          plotOutput("general_plots"),
          plotOutput("expression_plots"),
          DT::dataTableOutput("expression_table")
    )
))

server <- function(input, output){

    output$general_plots <- renderPlot({ 
        p1 <- DimPlot(sharma) + ggtitle("Types")
        p2 <- DimPlot(sharma, group.by = "ginty_ids") + ggtitle("Subtypes")
        grid.arrange(p1, p2, ncol = 2)
    })
    
    output$expression_plots <- renderPlot({
        if (!is.null(input$gene)) {
            p3 <- FeaturePlot(sharma, features = input$gene)
            p4 <- VlnPlot(sharma, features = input$gene, idents= input$comp1, pt.size = 0) + xlab("") + theme(legend.position = "none")
            grid.arrange(p3, p4, ncol = 2)
        } else {
        }
        })
    
    output$expression_table = DT::renderDataTable({
      if (input$task == "Explore a candidate gene" & !is.null(input$gene)) {
        if (!is.null(input$comp1) & length(input$comp1) == 2 & all(input$comp1 %in% c("P5 - Myelinated", "P5 - Unmyelinated"))){
          p5m.p5u %>% select(-1) %>%  filter(rownames(p5m.p5u) == input$gene)
          } else if ((!is.null(input$comp1) & length(input$comp1) == 2 & all(input$comp1 %in% c("P5 - Myelinated", "P0 - Myelinated")))){
          p5m.p0m %>% select(-1) %>% filter(rownames(p5m.p0m) == input$gene)
        } else if ((!is.null(input$comp1) & length(input$comp1) == 2 & all(input$comp1 %in% c("P5 - Myelinated", "P0 - Unmyelinated")))){
          p5m.p0u %>% select(-1) %>% filter(rownames(p5m.p0u) == input$gene)
        } else if ((!is.null(input$comp1) & length(input$comp1) == 2 & all(input$comp1 %in% c("P5 - Unmyelinated", "P0 - Myelinated")))){
          p5u.p0m %>% select(-1) %>% filter(rownames(p5u.p0m) == input$gene)
        } else if ((!is.null(input$comp1) & length(input$comp1) == 2 & all(input$comp1 %in% c("P5 - Unmyelinated", "P0 - Unmyelinated")))){
          p5u.p0u %>% select(-1)%>% filter(rownames(p5u.p0u) == input$gene)
        } else if ((!is.null(input$comp1) & length(input$comp1) == 2 & all(input$comp1 %in% c("P0 - Myelinated", "P0 - Unmyelinated")))){
          p0m.p0u %>% select(-1) %>% filter(rownames(p0m.p0u) == input$gene)
        } else if ((!is.null(input$comp1) & length(input$comp1) == 4 & all(input$comp1 %in% c("P5 - Myelinated", "P5 - Unmyelinated", "P0 - Myelinated", "P0 - Unmyelinated")))){
          all.groups %>% select(-1) %>% filter(rownames(all.groups) == input$gene)
        }  
      } else if (input$task == "Explore DEGs between groups"){
        if (!is.null(input$comp2) & length(input$comp2) == 2 & all(input$comp2 %in% c("P5 - Myelinated", "P5 - Unmyelinated"))){
          p5m.p5u %>% select(-1) %>% filter(p_val_adj < 0.001) %>%  filter(abs(avg_logFC) >= input$logfc)
        } else if ((!is.null(input$comp2) & length(input$comp2) == 2 & all(input$comp2 %in% c("P5 - Myelinated", "P0 - Myelinated")))){
          p5m.p0m %>% select(-1) %>% filter(p_val_adj < 0.001) %>%  filter(abs(avg_logFC) >= input$logfc)
        } else if ((!is.null(input$comp2) & length(input$comp2) == 2 & all(input$comp2 %in% c("P5 - Myelinated", "P0 - Unmyelinated")))){
          p5m.p0u %>% select(-1) %>% filter(p_val_adj < 0.001) %>%  filter(abs(avg_logFC) >= input$logfc)
        } else if ((!is.null(input$comp2) & length(input$comp2) == 2 & all(input$comp2 %in% c("P5 - Unmyelinated", "P0 - Myelinated")))){
          p5u.p0m %>% select(-1) %>% filter(p_val_adj < 0.001) %>%  filter(abs(avg_logFC) >= input$logfc)
        } else if ((!is.null(input$comp2) & length(input$comp2) == 2 & all(input$comp2 %in% c("P5 - Unmyelinated", "P0 - Unmyelinated")))){
          p5u.p0u %>% select(-1) %>% filter(p_val_adj < 0.001) %>%  filter(abs(avg_logFC) >= input$logfc)
        } else if ((!is.null(input$comp2) & length(input$comp2) == 2 & all(input$comp2 %in% c("P0 - Myelinated", "P0 - Unmyelinated")))){
          p0m.p0u %>% select(-1) %>% filter(p_val_adj < 0.001) %>%  filter(abs(avg_logFC) >= input$logfc)
        } else if ((!is.null(input$comp2) & length(input$comp2) == 4 & all(input$comp2 %in% c("P5 - Myelinated", "P5 - Unmyelinated", "P0 - Myelinated", "P0 - Unmyelinated")))){
          all.groups %>% select(-1) %>% filter(p_val_adj < 0.001) %>%  filter(abs(avg_logFC) >= input$logfc)
        } 
        } else {
      }
    }, callback = htmlwidgets::JS("
        var tips = ['Gene', 'Average log2 fold-change. Positive value indicates expression is higher in the first group compared to the second; negative indicates value is lower in the first group compared to the second. First/second order is determined by order in the input checkbox.', 
        'Percent of cells in the first group which express at least one read of the gene. First/second order is determined by order in the input checkbox.', 
        'Percent of cells in the second group which express at least one read of the gene. First/second order is determined by order in the input checkbox.', 
        'P-value of Wilcoxon rank-sum test, adjusted for multiple comparisons using Bonferroni correction.'],
            header = table.columns().header();
        for (var i = 0; i < tips.length; i++) {
          $(header[i]).attr('title', tips[i]);
        }
"), 
    extensions = 'Buttons', 
    options = list(dom = 'Blfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
    pageLength = 10,
    lengthMenu = c(10, 25, 50, 100, 100000))
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
