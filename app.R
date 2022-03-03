# profvis({
# 
# credentials <- data.frame(
#   user = c("nagellab", "shinymanager"), # mandatory
#   password = c("Circadian_2", "BonnotGillard"),
#   stringsAsFactors = FALSE
# )


# Packages ----
library(shiny)
library(shinymanager)
library(shinydashboard)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(scales)
library(readr)
library(leaflet)
library(gridExtra)
library(svglite)

library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggplotify)
library(patchwork)
library(gt)

library(visNetwork)
library(geomnet)
library(igraph)

# Loading data ----
#setwd("C:/Users/titou/Documents/Post doc/Shiny app/CASTR/")


#========= Data Tab "Single genes" ========
#==========================================#

exp <- read_csv("Time_of_day.csv")
exp$Temp <- as.factor(exp$Temp)

rec<-read_csv("Recovery.csv")
#rec$Time = as.factor(rec$Time)
rec$Temperature = as.factor(rec$Temperature)


group.labs <- c("Total" = "Transcriptome",
                "TRAP"= "Translatome")

expB<-exp
expB$Time = as.factor(expB$Time)

phase_metacycle<-read_csv("Phases_Metacycle.csv")
phase_metacycle$meta2d_BH.Q<-round(phase_metacycle$meta2d_BH.Q, digits = 4)

#========= Data Tab "Multiple genes" ======
#==========================================#

# import the list of transcription factors (users can visualize a particular TF family)
TF <- read.table("Arabidopsis TFs.txt", header = T)

# Import expression data (data at 22C)
rlog.mean.TOT <- read.csv("rlog.mean.TOT.csv", header = T)
rlog.mean.TRAP <- read.csv("rlog.mean.TRAP.csv", header = T)
rlog.mean.TOT$Phase2 <- as.character(rlog.mean.TOT$Phase2)
rlog.mean.TRAP$Phase2 <- as.character(rlog.mean.TRAP$Phase2)

# Import the temperature response data (37C vs 22C)
data.LFC.TOT <- read.csv("data.LFC.TOT.csv")
data.LFC.TRAP <- read.csv("data.LFC.TRAP.csv")
data.LFC.TOT$Phase2 <- as.character(data.LFC.TOT$Phase2)
data.LFC.TRAP$Phase2 <- as.character(data.LFC.TRAP$Phase2)

# Import the recovery data
data.Recovery.TOT <- read.csv("Recovery_LFC_TOT.csv", header = T)
data.Recovery.TRAP <- read.csv("Recovery_LFC_TRAP.csv", header = T)
data.Recovery.TOT$Phase2 <- as.character(data.Recovery.TOT$Phase2)
data.Recovery.TRAP$Phase2 <- as.character(data.Recovery.TRAP$Phase2)

paletteLength <- 50
myColor <- colorRampPalette(c("#3300FF","yellow"))(paletteLength)
myBreaks <- c(seq((-3), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(3/paletteLength, 3, length.out=floor(paletteLength/2)))
myColor2 <- colorRampPalette(c("Blue","White", "Red"))(paletteLength)
myBreaks2 <- c(seq((-3), 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(3/paletteLength, 3, length.out=floor(paletteLength/2)))

annotation_colors = list(
  Phase = c('0'='#FFCC00', '1.5'='#FFCC00',
            '3'='#FF9900', '4.5'='#FF9900',
            '6'='#FF6600', '7.5'='#FF6600',
            '9'='#FF0000', '10.5'='#FF0000',
            '12'='#CC0099', '13.5'='#CC0099',
            '15'='#9900FF', '16.5'='#9900FF',
            '18'='#0000FF', '19.5'='#0000FF',
            '21'='#000066', '22.5'='#000066'),
  Cycling = c(Yes="Black", No="White"))

#========= Data Tab "Phases" =============
#=========================================#

Phase_enrichment <- read_csv("Phases.csv")
Ref_list<-as.vector(unique(Phase_enrichment$Dataset))
Ref_list<-sort(Ref_list)
For_PhasEnr_Legend<-read.table("Upload_test_phases.txt", header=T)

#========= Data Tab "Network" =============
#=========================================#

TFfamilylist <- read.table("Arabidopsis TFs.txt", header = T)
TF_list<-as.vector(unique(TFfamilylist$Family))
TF_list<-sort(TF_list)

Nodes_original <- read.csv("Nodes2.csv", header = T)
# Nodes_original$Phases_transcriptome<-as.factor(Nodes_original$Phases_transcriptome)
# Nodes_original$Phases_transcriptome <- factor(Nodes_original$Phases_transcriptome, levels = c("0-1.5","3-4.5","6-7.5","9-10.5","12-13.5","15-16.5","18-19.5","21-22.5","Not circadian"))
# Nodes_original$Phases_translatome <- factor(Nodes_original$Phases_translatome, levels = c("0-1.5","3-4.5","6-7.5","9-10.5","12-13.5","15-16.5","18-19.5","21-22.5","Not circadian"))
# Nodes_original$Heat_stress_transcriptome <- factor(Nodes_original$Heat_stress_transcriptome, levels = c("> 3", "2-3", "1-2", "0-1", "(-1) - 0", "(-1) - (-2)", "(-2) - (-3)", "< (-3)"))
# Nodes_original$Heat_stress_translatome <- factor(Nodes_original$Heat_stress_translatome, levels = c("> 3", "2-3", "1-2", "0-1", "(-1) - 0", "(-1) - (-2)", "(-2) - (-3)", "< (-3)"))

Edges_original <- read.csv("edges.csv", header = T)
clock <- read.csv("clock_nodes.csv", header = T)

#========= Data Tab "Multispecies circadian oscillations" ========
#=================================================================#

Arabidopsis_data <- read.csv("arabidopsis_data.csv", header = T)
Rice_data <- read.csv("rice_data.csv", header = T)
Brassica_data <- read.csv("brassica_data.csv", header = T)
Barley_data <- read.csv("barley_data.csv", header = T)
Grey_areas <- read.csv("Grey_areas.csv", header = T)

Ara_dataset <- as.vector(unique(Arabidopsis_data$Condition))
Rice_dataset <- as.vector(unique(Rice_data$Condition))
Bra_dataset <- as.vector(unique(Brassica_data$Condition))
Barley_dataset <- as.vector(unique(Barley_data$Condition))


#========= Other =============

#https://stackoverflow.com/questions/14452465/how-to-create-textarea-as-input-in-a-shiny-webapp-in-r
textareaInput <- function(id, label, value, rows=20, cols=35, class="form-control"){
  tags$div(
    class="form-group shiny-input-container",
    tags$label('for'=id,label),
    tags$textarea(id=id,class=class,rows=rows,cols=cols,value))
}

#-----------------------------------------------------------------------------------------------#
# ui.R ---- contains information about the layout of the app as it appears in web browser ------
#-----------------------------------------------------------------------------------------------#

ui <- fluidPage(
  # tags$h2("Secured access"),
  # verbatimTextOutput("auth_output"),
  tags$head(tags$style(HTML('
         #sidebar {background-color: #bcd3e7;font-size: 14px}
        .tabbable > .nav > li > a {background-color: #E4EDF5; color:#6b8eb7;font-size: 15px}
        .tabbable > .nav > li[class=active]    > a {background-color:  #3c8dbc; color:white;font-weight: bold;font-size: 15px}
        .tabbable > .nav > li > a[data-value="t1"] {background-color: #ECF0F5; color:#585858;font-size: 15px}
        .tabbable > .nav > li[class=active]    > a[data-value="t1"] {background-color:  #3c8dbc; color:white;font-weight: bold;font-size: 15px}
        .shiny-notification {position: fixed; top: 20% ;left: 30%; right: 40%}
        # .form-group {
        #     margin-bottom: 0 !important;
        #   }
        
        body, label, input, button, select { 
          font-family: "Arial";
        }')
  )),
  dashboardPage(
    dashboardHeader(title= "CAST-R: A shiny application to visualize and identify Circadian And heat STress-Responsive genes", titleWidth = 1000,
                    tags$li(class = "dropdown",
                            tags$a(href="http://snap.uaf.edu", target="_blank", 
                                   tags$img(height= "40px", alt="Logo", src="logo6.png", style="vertical-align:middle;margin:-10px -10px")
                            ))
    ),
    dashboardSidebar(disable=T),
    dashboardBody(
      tabsetPanel(type="tabs",
                  #### UI Tab "Introduction" ####
                  tabPanel(value="t1", title = "Introduction", 
                           fluidRow(box(width=6, title = span("Welcome", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                        htmlOutput("Welcome"), status="primary"),
                                    box(width=6, title = span("Quick start", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                        htmlOutput("Summary"), status="primary"))
                  ),
                  #### UI Tab 1 "Single genes" ####
                  tabPanel(title = "Single genes",
                           sidebarLayout(
                             sidebarPanel(id="sidebar",
                                          textInput(inputId = "AGI",label = "1. Enter gene ID", value = "", placeholder = ""),
                                          div(HTML("ex: AT1G01060"),style = "margin-bottom:6px"),
                                          div(HTML("<i>Make sure there is no extra space before characters</i>"),style = "margin-bottom:15px"),
                                          checkboxGroupInput(inputId = "Plot_to_show", label = "2. Select data to show",
                                                             choices =  c("Time course at 22\u00B0C"="TC",
                                                                          "Heat stress response"="HSR",
                                                                          "Recovery following heat stress"="REC"), selected = "TC"),
                                          actionButton(inputId = "Submit_Tab1", label = "Submit"),
                                          width=2),
                             mainPanel(
                               fluidRow(box(width=10,title = span("Instructions and methodological details", style = "color: #000000; font-weight: bold; font-size: 20px"),
                                            htmlOutput("Tab1_description"), status="warning", solidHeader = F, collapsible = T, collapsed = T,
                                            div(style = "margin-top:20px"))),
                               conditionalPanel(condition = "input.Plot_to_show.includes('TC')",
                               fluidRow(box(width=10,title = span("Time course at 22\u00B0C", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                            textOutput("beginTab1"),
                                            plotOutput("One", height="330px"), status="primary", solidHeader = F, collapsible = T,
                                            div(style = "margin-top:20px"),
                                            htmlOutput("Legend_single_timecourse"),
                                            downloadButton("Sing_Time_course_22PNG","Plot.png"),
                                            downloadButton("Sing_Time_course_22SVG","Plot.svg"),
                                            downloadButton("Sing_Time_course_22PDF","Plot.pdf"),
                                            downloadButton("Sing_Data_Time_course_22", "Data.txt")))
                               ),
                               conditionalPanel(condition = "input.Plot_to_show.includes('HSR')",
                               fluidRow(box(width=10,title = span("Heat stress response", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                            plotOutput("Two"), status="primary", collapsible = T,
                                            htmlOutput("Legend_single_heat"),
                                            downloadButton("Sing_Heat_stressPNG","Plot.png"),
                                            downloadButton("Sing_Heat_stressSVG","Plot.svg"),
                                            downloadButton("Sing_Heat_stressPDF","Plot.pdf"),
                                            downloadButton("Sing_Data_Heat_stress", "Data.txt")))
                               ),
                               conditionalPanel(condition = "input.Plot_to_show.includes('REC')",
                               fluidRow(box(width=8,id="Recov_box",
                                            title = span("Recovery following heat stress", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                            plotOutput("Three"), status="primary", collapsible = T,
                                            div(style = "margin-top:20px"),
                                            htmlOutput("Legend_single_recovery"),
                                            downloadButton("Sing_RecoveryPNG","Plot.png"),
                                            downloadButton("Sing_RecoverySVG","Plot.svg"),
                                            downloadButton("Sing_RecoveryPDF","Plot.pdf"),
                                            downloadButton("Sing_Data_Recovery", "Data.txt")))
                               )
                             ))),
                  #### UI Tab 2 "Multiple genes" ####
                  tabPanel(title = "Multiple genes",
                           sidebarLayout(
                             sidebarPanel(id="sidebar",
                                          radioButtons(inputId="List_gen",label="1. Choose or paste a list of genes", selected=character(0),
                                                       choices=c("Choose a transcription factor family" = "fam",
                                                                 "Paste a list of genes" = "paste")),
                                          conditionalPanel('input.List_gen === "fam"', selectizeInput(inputId = "TFfamil", label = "2. Pick a TF family", 
                                                                                                      choices = TF_list, selected = character(0), multiple =T, options=list(placeholder = '', maxItems=1))),
                                          conditionalPanel('input.List_gen === "paste"', textareaInput(id="Genes_list_Tab2","2. Paste a list of genes","",rows=20),
                                                           div(HTML("<i>One AGI per line, no separator</i>"),style = "margin-bottom:15px")),
                                          conditionalPanel('input.List_gen === "fam"', radioButtons(inputId = "Sort_by", label = "3. Sort by", selected="order_by_transcriptome_ph",
                                                                                                    choices = c("Phase at the transcriptome level" = "order_by_transcriptome_ph",
                                                                                                                "Phase at the translatome level" = "order_by_translatome_ph",
                                                                                                                "Heat stress response at the transcriptome level" = "Heat_transcriptome",
                                                                                                                "Heat stress response at the translatome level" = "Heat_translatome",
                                                                                                                "Heat stress recovery at the transcriptome level" = "Recovery_transcriptome",
                                                                                                                "Heat stress recovery at the translatome level" = "Recovery_tranlatome"))),
                                          conditionalPanel('input.List_gen === "paste"', radioButtons(inputId = "Sort_by1", label = "3. Sort by", selected="order_by_transcriptome_ph",
                                                                                                      choices = c("Phase at the transcriptome level" = "order_by_transcriptome_ph",
                                                                                                                  "Phase at the translatome level" = "order_by_translatome_ph",
                                                                                                                  "Heat stress response at the transcriptome level" = "Heat_transcriptome",
                                                                                                                  "Heat stress response at the translatome level" = "Heat_translatome",
                                                                                                                  "Heat stress recovery at the transcriptome level" = "Recovery_transcriptome",
                                                                                                                  "Heat stress recovery at the translatome level" = "Recovery_tranlatome",
                                                                                                                  "Input order" = "order_by_input"))),
                                          actionButton(inputId = "Submit_Tab2", label = "Submit"),
                                          width=2),
                             mainPanel(
                               fluidRow(box(width=10,title = span("Instructions and methodological details", style = "color: #000000; font-weight: bold; font-size: 20px"),
                                            htmlOutput("Tab2_description"), status="warning", solidHeader = F, collapsible = T, collapsed = T,
                                            div(style = "margin-top:20px"))),
                               fluidRow(box(width=10,title = span("Time course at 22\u00B0C", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                            textOutput("beginTab2"),
                                            plotOutput("Four"),
                                            #dataTableOutput("mytable"),
                                            status="primary", solidHeader = F,collapsible = T,
                                            htmlOutput("Legend_heatmap_timecourse"),
                                            downloadButton("Mult_Time_course_22PNG","Plot.png"),
                                            downloadButton("Mult_Time_course_22SVG","Plot.svg"),
                                            downloadButton("Mult_Time_course_22PDF","Plot.pdf"),
                                            downloadButton("Mult_DataTOT_Time_course_22", "Data_Transcriptome.txt"),
                                            downloadButton("Mult_DataTRAP_Time_course_22", "Data_Translatome.txt"))),
                               
                               #div(style = "margin-top:20px"),
                               fluidRow(box(width=10,title = span("Heat stress response", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                            plotOutput("Five"), status="primary",collapsible = T,
                                            htmlOutput("Legend_heatmap_heat"),
                                            downloadButton("Mult_Heat_stressPNG","Plot.png"),
                                            downloadButton("Mult_Heat_stressSVG","Plot.svg"),
                                            downloadButton("Mult_Heat_stressPDF","Plot.pdf"),
                                            downloadButton("Mult_DataTOT_Heat_stress", "Data_Transcriptome.txt"),
                                            downloadButton("Mult_DataTRAP_Heat_stress", "Data_Translatome.txt"))),
                               
                               #div(style = "margin-top:20px"),
                               fluidRow(box(width=10,id="Recov_box",
                                            title = span("Recovery following heat stress", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                            plotOutput("Six"), status="primary",collapsible = T,
                                            htmlOutput("Legend_heatmap_recovery"),
                                            downloadButton("Mult_RecoveryPNG","Plot.png"),
                                            downloadButton("Mult_RecoverySVG","Plot.svg"),
                                            downloadButton("Mult_RecoveryPDF","Plot.pdf"),
                                            downloadButton("Mult_DataTOT_Recovery", "Data_Transcriptome.txt"),
                                            downloadButton("Mult_DataTRAP_Recovery", "Data_Translatome.txt")))
                             ))),
                  #### UI Tab 4 "Network" ####
                  tabPanel(title = "Network",
                           sidebarLayout(
                             sidebarPanel(id="sidebar",
                                          checkboxGroupInput(inputId = "Clock_to_show", label = "1. Select clock proteins",
                                                             choices =  c("CCA1","LHY","LUX","PRR5","PRR7","PRR9","TOC1"), selected = c("CCA1","LHY","LUX","PRR5","PRR7","PRR9","TOC1")),
                                          radioButtons(inputId="Chose_list",label="2. Select input genes", selected=character(0),
                                                       choices=c("Transcription factor family" = "TFfamily",
                                                                 "Other list of genes" = "paste_list_network")),
                                          conditionalPanel('input.Chose_list === "TFfamily"', selectizeInput(inputId = "TFfam", label = "Pick a TF family", 
                                                                                                             choices = TF_list, selected = character(0), multiple =T, options=list(placeholder = '', maxItems=1))),
                                          conditionalPanel('input.Chose_list === "paste_list_network"', textareaInput("Genes_list_Tab4","Paste a list of genes","",rows=20),
                                                           div(HTML("<i>One AGI per line, no separator</i>"),style = "margin-bottom:15px")),
                                          actionButton(inputId = "Submit_Tab4", label = "Submit"),
                                          div(style = "margin-top:20px"),
                                          radioButtons(inputId = "Color_choice", label = "3. Criteria for node color",selected="color_Transcriptome_phases",
                                                       choices = c("Phases at the transcriptome level" = "color_Transcriptome_phases",
                                                                   "Phases at the translatome level" = "color_Translatome_phases",
                                                                   "Heat stress response at the transcriptome level" = "color_Transcriptome_heatstress",
                                                                   "Heat stress response at the translatome level" = "color_Translatome_heatstress",
                                                                   "No color" = "No_col")),
                                          width=2),
                             mainPanel(column(width=12,
                                              fluidRow(box(width=9,title = span("Instructions and methodological details", style = "color: #000000; font-weight: bold; font-size: 20px"),
                                                           htmlOutput("Tab4_description"), status="warning", solidHeader = F, collapsible = T, collapsed = T,
                                                           div(style = "margin-top:20px"),
                                                           
                                                           fluidRow(column(width=6,
                                                                           gt_output("Clock_gene_AGI")
                                                                           #tableOutput("mytable")
                                                           ),
                                                           column(width=6,
                                                                  uiOutput("Clock_genes_image"))),
                                                           div(style = "margin-top:20px"),
                                                           fluidRow(column(width=12,
                                                                           htmlOutput("tab4_description2")
                                                           )),),
                                              fluidRow(column(width=9,box(width=12, title = span("Network of genes targeted by clock proteins", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                                                          textOutput("beginTab4"),
                                                                          visNetworkOutput("Eleven"), status="primary",
                                                                          htmlOutput("Legend_network")),
                                                              downloadButton("Nodes_file", "Data_nodes.txt"),
                                                              downloadButton("Edges_file", "Data_edges.txt")),
                                                       column(width=3, align="center",
                                                              htmlOutput("Selection_info_Network"),
                                                              div(style = "margin-top:20px"),
                                                              plotOutput("Tab4_Legend"),
                                                              div(style = "margin-top:20px"),
                                                              textOutput("node_info"),
                                                              uiOutput("click_ui"), status="primary"))
                                                           )
                                              )
                                              )
                             )
                           ),
                  #### UI Tab 3 "Phase enrichment" ####
                  tabPanel(title = "Phase enrichment",
                           sidebarLayout(
                             sidebarPanel(id="sidebar",
                                          radioButtons(inputId="Chose_ref",label="1. Select phase reference dataset type", selected=character(0),
                                                       choices=c("Existing reference" = "exist",
                                                                 "Use your own reference" = "brow")),
                                          conditionalPanel('input.Chose_ref === "exist"', selectizeInput(inputId = "Existing_ref", label = "2. Pick a dataset", 
                                                                                                         choices = Ref_list, selected = character(0), multiple =T, options=list(placeholder = '', maxItems=1))),
                                          conditionalPanel('input.Chose_ref === "brow"', fileInput("file1", "2. Upload CSV File", accept = c("text/csv",  "text/comma-separated-values,text/plain", ".csv")),
                                                           div(style = "margin-top:-30px"),
                                                           htmlOutput("Instructions_uploadCSV"),
                                                           div(style = "margin-top:20px")),
                                          textareaInput("Genes_list_Tab3","3. Paste a list of genes","",rows=20),
                                          div(HTML("<i>One AGI per line, no separator</i>"),style = "margin-bottom:15px"),
                                          div(HTML("We suggest providing a minimal list size of 100 genes for a meaningful enrichment calculation"),style = "margin-bottom:15px"),
                                          actionButton(inputId = "Submit_Tab3", label = "Submit"),
                                          width=2),
                             mainPanel(
                               fluidRow(box(width=12,title = span("Instructions and methodological details", style = "color: #000000; font-weight: bold; font-size: 20px"),
                                            htmlOutput("Tab3_description"), status="warning", solidHeader = F, collapsible = T, collapsed = T,
                                            div(style = "margin-top:20px"))),
                               column(width=6,
                                      fluidRow(box(width=12,title = span("Phase distribution in the defined phase reference dataset", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                                   textOutput("beginTab3"),
                                                   plotOutput("Seven"), status="primary",
                                                   div(style = "margin-top:20px"),
                                                   htmlOutput("Legend_phase_ref"),
                                                   downloadButton("Phase_distrib_refPNG","Plot.png"),
                                                   downloadButton("Phase_distrib_refSVG","Plot.svg"),
                                                   downloadButton("Phase_distrib_refPDF","Plot.pdf"),
                                                   downloadButton("Genes_phase_ref", "Genes_phase_ref.txt"))),
                                      fluidRow(box(width=12, title = span("Phase enrichment", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                                   plotOutput("Nine", height="450px"), status="primary",
                                                   htmlOutput("Legend_phase_enrichment"),
                                                   downloadButton("Phase_enrichmentPNG","Plot.png"),
                                                   downloadButton("Phase_enrichmentSVG","Plot.svg"),
                                                   downloadButton("Phase_enrichmentPDF","Plot.pdf")))
                               ),
                               column(width=6,
                                      fluidRow(box(width=12,title = span("Phase distribution in the user subset of genes", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                                   plotOutput("Eight"), status="primary",
                                                   div(style = "margin-top:20px"),
                                                   htmlOutput("Legend_phase_subset"),
                                                   downloadButton("Phase_distrib_subsetPNG","Plot.png"),
                                                   downloadButton("Phase_distrib_subsetSVG","Plot.svg"),
                                                   downloadButton("Phase_distrib_subsetPDF","Plot.pdf"),
                                                   downloadButton("Genes_phase_subset", "Genes_phase_subset.txt"))),
                                      htmlOutput("Selection_info"),
                                      div(style = "margin-top:20px"),
                                      fluidRow(box(width=12,title = span("Phase enrichment summary", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                                   gt_output("Ten"), status="primary",
                                                   # dataTableOutput("mytable"), status="primary",
                                                   div(style = "margin-top:20px"),
                                                   htmlOutput("Legend_summary_table"),
                                                   downloadButton("Phase_enrich_summary", "Data.txt")))
                               )
                             )
                           )
                  ),
                  #### UI Tab 5 "Multispecies circadian oscillations" ####
                  tabPanel(title = "Multispecies circadian oscillations",
                           sidebarLayout(
                             sidebarPanel(id="sidebar", width=2,
                                          selectizeInput(inputId = "Existing_ref1", label = "1. Pick an Arabidopsis time course dataset", 
                                                         choices = Ara_dataset, selected = character(0), multiple =T, options=list(placeholder = '', maxItems=1)),
                                          textInput(inputId = "AGI_oscillation",label = "Enter Arabidopsis gene ID", value = "", placeholder = ""),
                                          div(HTML("ex: AT1G01060"),style = "margin-bottom:6px"),
                                          div(HTML("<i>Make sure there is no extra space before characters</i>"),style = "margin-bottom:15px"),
                                          actionButton(inputId = "Submit_Tab5", label = "Submit"),
                                          div(style = "margin-top:30px"),
                                          
                                          selectizeInput(inputId = "Existing_ref3", label = "2. Pick a Brassica time course dataset", 
                                                         choices = Bra_dataset, selected = character(0), multiple =T, options=list(placeholder = '', maxItems=1)),
                                          textInput(inputId = "AGI_oscillation3",label = "Enter Brassica gene ID", value = "", placeholder = ""),
                                          div(HTML("ex: BraA10g01800R"),style = "margin-bottom:6px"),
                                          div(HTML("<i>Make sure there is no extra space before characters</i>"),style = "margin-bottom:15px"),
                                          actionButton(inputId = "Submit_Tab5.2", label = "Submit"),
                                          div(style = "margin-top:30px"),
                                          
                                          selectizeInput(inputId = "Existing_ref4", label = "3. Pick a Barley time course dataset", 
                                                         choices = Barley_dataset, selected = character(0), multiple =T, options=list(placeholder = '', maxItems=1)),
                                          textInput(inputId = "AGI_oscillation4",label = "Enter Barley gene ID", value = "", placeholder = ""),
                                          div(HTML("ex: HORVU7Hr1G070870"),style = "margin-bottom:6px"),
                                          div(HTML("<i>Make sure there is no extra space before characters</i>"),style = "margin-bottom:15px"),
                                          actionButton(inputId = "Submit_Tab5.3", label = "Submit"),
                                          div(style = "margin-top:30px"),
                                          
                                          selectizeInput(inputId = "Existing_ref2", label = "4. Pick a Rice time course dataset", 
                                                         choices = Rice_dataset, selected = character(0), multiple =T, options=list(placeholder = '', maxItems=1)),
                                          textInput(inputId = "AGI_oscillation2",label = "Enter Rice gene ID", value = "", placeholder = ""),
                                          div(HTML("ex: LOC_Os08g06110"),style = "margin-bottom:6px"),
                                          div(HTML("<i>Make sure there is no extra space before characters</i>"),style = "margin-bottom:15px"),
                                          actionButton(inputId = "Submit_Tab5.1", label = "Submit"),
                                          div(style = "margin-top:30px")),
                             mainPanel(
                               fluidRow(box(width=12,title = span("Instructions and methodological details", style = "color: #000000; font-weight: bold; font-size: 20px"),
                                            htmlOutput("Tab5_description"), status="warning", solidHeader = F, collapsible = T, collapsed = T,
                                            div(style = "margin-top:20px"))),
                               column(width=12,
                                      fluidRow(
                                        column(width=6,
                                               box(width=12,title = span(HTML("Time course in <em>Arabidopsis thaliana</em>"), style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                                   textOutput("beginTab5"),
                                                   plotOutput("Twelve"), status="primary",
                                                   htmlOutput("Legend_Arabido"),
                                                   downloadButton("Arabidopsis_PNG","Plot.png"),
                                                   downloadButton("Arabidopsis_SVG","Plot.svg"),
                                                   downloadButton("Arabidopsis_PDF","Plot.pdf"),
                                                   downloadButton("Arabidopsis_timecourse", "Data.txt"))),
                                        column(width=6,
                                               box(width=12,title = span(HTML("Time course in <em>Brassica rapa</em>"), style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                                   textOutput("beginTab5.2"),
                                                   plotOutput("Fourteen"), status="primary",
                                                   htmlOutput("Legend_Brassica"),
                                                   downloadButton("Brassica_PNG","Plot.png"),
                                                   downloadButton("Brassica_SVG","Plot.svg"),
                                                   downloadButton("Brassica_PDF","Plot.pdf"),
                                                   downloadButton("Brassica_timecourse", "Data.txt")))),
                                      
                                      fluidRow(
                                        column(width=6, 
                                               box(width=12,title = span(HTML("Time course in <em>Hordeum vulgare</em>"), style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                                   textOutput("beginTab5.3"),
                                                   plotOutput("Fifteen"), status="primary",
                                                   htmlOutput("Legend_Barley"),
                                                   downloadButton("Barley_PNG","Plot.png"),
                                                   downloadButton("Barley_SVG","Plot.svg"),
                                                   downloadButton("Barley_PDF","Plot.pdf"),
                                                   downloadButton("Barley_timecourse", "Data.txt"))), 
                                        column(width=6,
                                               box(width=12,title = span(HTML("Time course in <em>Oryza sativa</em>"), style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                                   textOutput("beginTab5.1"),
                                                   plotOutput("Thirteen"), status="primary",
                                                   htmlOutput("Legend_Rice"),
                                                   downloadButton("Rice_PNG","Plot.png"),
                                                   downloadButton("Rice_SVG","Plot.svg"),
                                                   downloadButton("Rice_PDF","Plot.pdf"),
                                                   downloadButton("Rice_timecourse", "Data.txt"))))
                               )
                               
                               
                             )
                           )
                  ),
                  #### UI Tab 6 "How to cite" ####
                  tabPanel(value="t1", title = "How to cite",
                           fluidRow(box(width=6, title = span("References of datasets and tools", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                        htmlOutput("cite_refs"), status="primary"),
                                    box(width=6, title = span("How to cite CAST-R", style = "color: #3c8dbc; font-weight: bold; font-size: 20px"),
                                        htmlOutput("cite_CASTR"), status="primary"))
                  )
                  
      ))
  )
)

#ui <- secure_app(ui)


#---------------------------------------------------------------------------------------------------------------------------------------------------#
# server.R ---- contains information about the computation of the app, creating plots, tables, maps etc. using information provided by the user ----
#---------------------------------------------------------------------------------------------------------------------------------------------------#

server <- function(input, output,session) {
  # res_auth <- secure_server(
  #   check_credentials = check_credentials(credentials)
  #)
  
  # output$auth_output <- renderPrint({
  #   reactiveValuesToList(res_auth)
  #})
  
  ####========= Server Tab "Introduction ========
  #=======================================#
  
  output$Welcome <-renderText({# Welcome
    paste("<p style='text-align:justify'>","Here we present a Shiny application, CAST-R, that allows the user to identify and visualize circadian and heat stress-responsive transcripts.
          It notably uses our recently published datasets obtained at the transcriptome and translatome levels in Arabidopsis over a 24 h time course, in response to heat stress and during the plant recovery following the stress (<a href='https://doi.org/10.1093/plcell/koab113'>Bonnot and Nagel, 2021</a>).
          More specifically, this application allows the user to generate and export profiles and heatmaps representing transcript abundance of a single or of multiple genes in these datasets. In addition, the application takes advantage of published Arabidopsis ChIP-Seq datasets to visualize in an interactive network the connections between clock proteins and their targets.
          Nodes can be sorted either by timing of expression or by heat stress responsiveness. We also implemented a tool to perform phase (i.e. timing of expression) enrichment analyses.
          This functionality combines statistical analyses and graphical representations to identify significantly over- and under-represented phases within a subset of genes.
          Lastly, this Shiny app allows the user to visualize circadian oscillations in other published datasets in <i>Arabidopsis thaliana</i>, <i>Brassica rapa</i>, <i>Hordeum vulgare</i> and <i>Oryza sativa</i>.",
          "</p>","<p style='text-align:justify'>",
          'At the top of each tab, a box entitled "Instructions and methodological details" provides detailed information about the functionalities available in the tab and how to use it, as well as methodological details. 
          The main points to have in mind to properly navigate the tabs are summarized on the right, in the "Quick start" section.',
          "</p>","<p style='text-align:justify'>",
          '<b>The user session will timeout after 30 minutes when not in active use, to prevent unnecessary consumption of server resources. In addition, server disconnection or slow app can also come from user internet connection issues or browser performances.</b>',
          "</p"
    )
  })
  
  output$Summary <-renderText({# Summary/Quick start
    paste("<p style='text-align:justify'>",
          "The following short description summarizes the main functionalities of the Shiny application. For more details, please read the full description for each individual tab in the corresponding box 'Instructions and methodological details'.",
          "</p>","<p style='text-align:justify'>",
          '- <b>Tab "Single genes"</b>: From an Arabidopsis AGI number, i) identify if the gene exhibits circadian oscillations, ii) identify the gene phase (timing of peak expression), iii) visualize the transcript abundance during the day, in response to heat stress and during the plant recovery following heat stress, at the transcriptome and translatome levels.',
          "</p>","<p style='text-align:justify'>",
          '- <b>Tab "Multiple genes"</b>: From a list of Arabidopsis AGI numbers, i) identify if the genes exhibit circadian oscillations, ii) identify their phase (timing of peak expression), iii) visualize as heatmaps the transcript abundance during the day, in response to heat stress and during the plant recovery following heat stress, at the transcriptome and translatome levels.',
          "</p>","<p style='text-align:justify'>",
          '- <b>Tab "Network"</b>: From a list of Arabidopsis AGI numbers, i) identify what genes are targeted by clock proteins based on published ChIP-Seq data and ii) visualize these interactions in an interactive network, highlighting genes based on their phase or their heat stress response.',
          "</p>","<p style='text-align:justify'>",
          '- <b>Tab "Phase enrichment"</b>: From a list of Arabidopsis AGI numbers, identify how many genes exhibit significant circadian oscillations and if any phases (timing of peak expression) are significantly over-represented in this list of genes.',
          "</p>","<p style='text-align:justify'>",
          '- <b>Tab "Multispecies circadian oscillations"</b>: Visualize transcript profiles of <i>Arabidopsis thaliana</i>, <i>Brassica rapa</i>, <i>Hordeum vulgare</i>, and <i>Oryza sativa</i> genes in multiple time course datasets.',
          "<br>",
          "If you would like your plant circadian dataset incorporated into CAST-R, please contact us at dawnn@ucr.edu.",
          "</p>"
    )
  })
  
  ####========= Server Tab 1 "Single genes" ========
  #=======================================#
  output$Tab1_description <-renderText({# Single genes
    paste("<p style='text-align:justify'>",
          "This tab generates plots from the Arabidopsis transcriptome and translatome datasets published in <a href='https://doi.org/10.1093/plcell/koab113'>Bonnot and Nagel (2021)</a>.",
          "<br>",
          "You need to enter an Arabidopsis AGI number (ATXGXXXXX) in the designated area and click on submit. For example: AT1G01060 (<i>LHY</i> gene).",
          "<br>",
          "All data were obtained from 12 days old seedlings (whole seedlings) that were grown in light (12 h) and dark cycles (12 h) at 22\u00B0C for 10 days and then transferred to free-running conditions (continuous light and temperature) for two days before sampling. 
          Time of day is referred to as Zeitgeber Time (ZT) and corresponds to the hours after moving the seedlings in constant conditions (light and temperature). 
          At each time point, total mRNAs and ribosome-associated mRNAs (using Translating Ribosome Affinity Purification, TRAP) were isolated, and mRNA-Seq (transcriptome) and TRAP-Seq (translatome) was performed.",
          "</p>","<p style='text-align:justify'>",
          "- <b>Time course at 22\u00B0C:</b> A 24 h time course (with samples every 3 h) was analyzed under control conditions at 22\u00B0C. A detection of circadian oscillations was performed using the R package <a href='https://doi.org/10.1093/bioinformatics/btw405'>'Metacycle'</a>. 
          This analysis calculated a 'Phase' and an 'Adjusted p-value' for each transcript. In <a href='https://doi.org/10.1093/plcell/koab113'>Bonnot and Nagel (2021)</a>, oscillations were considered significant for genes: with an Adjusted p-value < 0.01; with an 0.01 < Adjusted p-value < 0.05 and that overlapped with previously published lists of circadian transcripts. For details, please see the Methods section in <a href='https://doi.org/10.1093/plcell/koab113'>Bonnot and Nagel (2021)</a>.",
          "</p>","<p style='text-align:justify'>",
          "- <b>Heat stress response:</b> Heat stress treatments were applied at different times of day (from ZT48, subjective dawn, to ZT69, end of subjective night) on separate sets of plants, and correspond to treatments at 37\u00B0C for 1 h.
          To compare the stress vs control conditions, a differential expression analysis was performed using the R package <a href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8'>'DESeq2'</a>, and FDR and Log2 Fold Change (log2 FC) values (37\u00B0C vs 22\u00B0C) were calculated.
          In <a href='https://doi.org/10.1093/plcell/koab113'>Bonnot and Nagel (2021)</a>, genes were considered as significantly differentially regulated when the FDR < 0.05 and the log2 FC > |1|.",
          "</p>","<p style='text-align:justify'>",
          "- <b>Recovery following heat stress:</b> Heat stress treatment at 37\u00B0C for 1 h was applied from ZT53 to ZT54, corresponding to the middle of the light period. 
          Following heat stress, seedlings were either sampled (0 h of recovery), or moved back to 22\u00B0C for recovery. To compare the stress (or recovery) vs control conditions, a differential expression analysis was performed using the R package <a href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8'>'DESeq2'</a>, and FDR and Log2 Fold Change (log2 FC) values were calculated.
          In <a href='https://doi.org/10.1093/plcell/koab113'>Bonnot and Nagel (2021)</a>, genes were considered as significantly differentially regulated when the FDR < 0.05 and the log2 FC > |1|.",
          "<br>", 
          "</p>","<p style='text-align:justify'>"
    )
  })
  
  output$beginTab1 <- renderText({
    validate(
      need(input$AGI, "Enter gene ID and click on the Submit button")
    )
  })
  
  Single_gene_data_ctrl<-eventReactive(input$Submit_Tab1, {
    req(input$AGI)
    temp<-c("22")
    list<-toupper(input$AGI)
    ctrl<-filter(exp, Temp %in% temp)
    Single_gene_data_ctrl<-filter(ctrl, AGI %in% list)
  })
  
  Single_gene_data_tot<-eventReactive(input$Submit_Tab1, {
    req(input$AGI)
    tot<-phase_metacycle[phase_metacycle$Type == "Total",]
    list<-toupper(input$AGI)
    Single_gene_data_tot<-filter(tot, AGI %in% list)
  })
  
  Single_gene_data_trap<-eventReactive(input$Submit_Tab1, {
    req(input$AGI)
    trap<-phase_metacycle[phase_metacycle$Type == "TRAP",]
    list<-toupper(input$AGI)
    Single_gene_data_trap<-filter(trap, AGI %in% list)
  })
  
  # Single_gene_data_tot<-Single_gene_data_tot()
  # Single_gene_data_trap<-Single_gene_data_trap()
  # 
  # 
  # g<-ggplot(data=Single_gene_data_ctrl()
  
  #=========== Tab 1 Box 1 - Time course at 22 ("One") =====
  
  Nb_One <- eventReactive(input$Submit_Tab1,{
    
    validate(
      need(input$AGI != "", "Provide a gene ID before clicking on the Submit button")
    )
    
    Single_gene_data_tot<-Single_gene_data_tot()
    Single_gene_data_trap<-Single_gene_data_trap()
    
    validate(
      need(nrow(Single_gene_data_trap)>0, "This gene is miswritten (check instructions) or below the threshold of detection based on our criteria")
    )
    
    withProgress(message = 'Making plot 1', value = 1, {
      
      g<-ggplot(data=Single_gene_data_ctrl(), aes(x = Time, y = rlog, group= Temp)) +
        geom_rect(aes(xmin=60, xmax=72, ymin=-Inf, ymax=Inf), 
                  fill="grey", alpha=0.05, color = NA)+
        geom_point(aes(shape = Temp, colour = Temp),
                   position=position_dodge(width=0.1), size = 4)+
        scale_shape_manual(values=c(16),labels=c("22\u00B0C"))+
        geom_line(aes(linetype = Temp, colour = Temp), size = 1,
                  position=position_dodge(width=0.1))+
        scale_linetype_manual(values=c("solid", "blank"),labels=c("22\u00B0C","37\u00B0C"))+
        geom_errorbar(aes(ymin=rlog-SD, ymax=rlog+SD, colour = Temp),width=0.5, size=0.8,
                      position=position_dodge(width=0.1),show.legend=F)+
        scale_color_manual(values=c("#666666","#FF6600"),labels=c("22\u00B0C","37\u00B0C"))+
        scale_x_continuous(breaks = seq(48,72,3))+
        scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
        theme_bw()+
        facet_wrap(~Type,labeller=as_labeller(group.labs))+
        xlab("Time of day (ZT)")+
        ylab("Transcript abundance\n(rlog normalized counts)")+
        theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
              axis.text=element_text(size=13, colour = "black"),
              axis.text.x=element_text(margin = unit(c(3, 0, 0, 0), "mm")),
              axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
              strip.text = element_text(size = 15),
              axis.title = element_text(size = 14),
              legend.title=element_text(size=13), 
              legend.text=element_text(size=13),
              plot.title = element_text(face="bold", size=18,hjust=-0.14,vjust=0),
              strip.background =element_rect(colour=NA, fill="white"),
              strip.text.x = element_text(margin = margin(0, 0, 0, 0, "cm"),color="white",size=4))+
        guides(colour=guide_legend(""),
               shape = guide_legend(""),
               linetype = guide_legend(""))
      
      lg<- cowplot::get_legend(g)
      lg<-as_ggplot(lg)
      g<-g+theme(legend.position = "none")
      
      na<-ggplot()+
        geom_rect()+
        theme_void()
      na
      
      pg<-plot_grid(g,na,lg,na,ncol=4,rel_widths=c(5.5,-0.1,1.1,0.1))
      
      #======
      
      #d <- head(tot[tot$AGI == input$AGI,3:5])
      d <- head(Single_gene_data_tot[,3:5])
      colnames(d) <- c("Phase","Adjusted \np-value","Cycling")
      j<-ggtexttable(d[1,], rows = NULL, 
                     theme = ttheme("classic", base_size = 13))
      title_j <- ggplot() + 
        labs(title = "Transcriptome")+
        theme_minimal()+
        theme(plot.title=element_text(hjust=0.5,vjust=0,margin = margin(-0.3, 0, 0, 0, "cm"),size=16,face="bold"))
      
      #d1 <- head(trap[trap$AGI == input$AGI,3:5])
      d1 <- head(Single_gene_data_trap[,3:5])
      colnames(d1) <- c("Phase","Adjusted \np-value","Cycling")
      j1<-ggtexttable(d1[1,], rows = NULL,
                      theme = ttheme("classic", base_size = 13))
      title_j1 <- ggplot() + 
        labs(title = "Translatome")+
        theme_minimal()+
        theme(plot.title=element_text(hjust=0.5,vjust=0,margin = margin(-0.3, 0, 0, 0, "cm"),size=16,face="bold"))
      
      pg2<-plot_grid(na,title_j,na,title_j1,na,
                     na,j,na,j1,na,nrow=2,ncol=5,rel_widths=c(0.5,2,0.05,2,0.95),rel_heights = c(0.5,1),align="v")
      

      plot_grid(pg2,na,pg,nrow=3, rel_heights=c(1.2,0.2,3.2), align="h",axis = "r")
      
    })
    
  })
  
  output$One <- renderPlot({
    print(Nb_One())
  })
  
  output$Legend_single_timecourse <-renderText({
    paste("<p style='text-align:justify'>","The phase is defined as the timing of peak abundance (a phase of 0 and 12 indicates a peak abundance at subjective dawn and the beginning of subjective night, respectively). 
          In the column 'Cycling', 'Yes' or 'No' indicate if the circadian oscillation was detected as significant or not in <a href='https://doi.org/10.1093/plcell/koab113'>Bonnot and Nagel (2021)</a>, respectively. 
          Of note, because the phase is calculated as the time point of highest transcript abundance, a gene can have a phase but not be significantly cycling.
          Adjusted p-values < 5e-5 are noted as '0'. The grey areas represent the subjective night. Data are means +/- sd for n = 3 biological replicates.",
          "</p>"
    )
  })
  
  output$Sing_Time_course_22PNG <- downloadHandler(
    filename = "Time_course_22.png.png",
    content = function(file){
      ggsave(file, plot = Nb_One(), width = 7, height = 5, device = "png")
    })
  
  output$Sing_Time_course_22SVG <- downloadHandler(
    filename = "Time_course_22.svg.svg",
    content = function(file){
      ggsave(file, plot = Nb_One(), width = 8, height = 5, device = "svg")
    })
  
  output$Sing_Time_course_22PDF <- downloadHandler(
    filename = "Time_course_22.pdf.pdf",
    content = function(file){
      pdf(file,width = 8, height = 5)
      print(Nb_One())
      dev.off()
    })
  
  
  # Downloadable txt of selected dataset
  output$Sing_Data_Time_course_22 <- downloadHandler(
    filename = "Data_Time_course_22.txt.txt",
    content = function(file) {
      #write.table(datasetInputOne(), file, row.names = FALSE)
      write.table(Single_gene_data_ctrl(), file, row.names = FALSE)
    }
  )
  
  #=========== Tab 1 Box 2 - Heat stress response ("Two") ======
  
  Nb_Two <- eventReactive(input$Submit_Tab1,{
    
    validate(
      need(input$AGI != "", "Provide a gene ID before clicking on the Submit button")
    )
    
    Single_gene_data_ctrl<-Single_gene_data_ctrl()
    
    validate(
      need(nrow(Single_gene_data_ctrl)>0, "This gene is miswritten (check instructions) or below the threshold of detection based on our criteria")
    )
    
    withProgress(message = 'Making plot 2', value = 1, {
      
      req(input$AGI)
      g0<-ggplot(data=exp[exp$AGI == toupper(input$AGI),], aes(x = Time, y = rlog, group= Temp)) +
        geom_rect(aes(xmin=60, xmax=72, ymin=-Inf, ymax=Inf), 
                  fill="grey", alpha=0.02, color = NA)+
        geom_point(aes(shape = Temp, colour = Temp),
                   position=position_dodge(width=0.1), size = 4)+
        scale_shape_manual(values=c(16,17),labels=c("22\u00B0C","37\u00B0C"))+
        geom_line(aes(linetype = Temp, colour = Temp), size = 1,
                  position=position_dodge(width=0.1))+
        scale_linetype_manual(values=c("solid", "blank"),labels=c("22\u00B0C","37\u00B0C"))+
        geom_errorbar(aes(ymin=rlog-SD, ymax=rlog+SD, colour = Temp),width=.5, size=0.8,
                      position=position_dodge(width=0.1),show.legend=F)+
        scale_color_manual(values=c("#666666","#FF6600"),labels=c("22\u00B0C","37\u00B0C"))+
        scale_x_continuous(breaks = seq(48,72,3))+
        scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
        theme_bw()+
        facet_wrap(~Type,labeller=as_labeller(group.labs))+
        xlab("Time of day (ZT)")+
        ylab("Transcript abundance\n(rlog normalized counts)")+
        #ggtitle("Heat response")+
        theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
              axis.text=element_text(size=13, colour = "black"),
              axis.text.x=element_text(margin = unit(c(3, 0, 0, 0), "mm")),
              axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
              strip.text = element_text(size = 15),
              axis.title = element_text(size = 14),
              legend.title=element_text(size=13), 
              legend.text=element_text(size=13),
              plot.title = element_text(face="bold", size=18,hjust=-0.14,vjust=0),
              strip.background =element_rect(colour=NA, fill="white"),
              strip.text.x = element_text(margin = margin(0, 0, 0.5, 0, "cm"),face="bold", color="black",size=16))+
        guides(colour=guide_legend(""),
               shape = guide_legend(""),
               linetype = guide_legend(""))
      
      lg<- cowplot::get_legend(g0)
      lg<-as_ggplot(lg)
      g0<-g0+theme(legend.position = "none")
      
      na<-ggplot()+
        geom_rect()+
        theme_void()
      na
      
      #pg0<-plot_grid(g0,na,lg,na,ncol=4,rel_widths=c(5,-0.1,1,0.3))
      
      g1<-ggplot(data=expB[expB$AGI == toupper(input$AGI)&expB$Type == "Total"&expB$Time%in% c(48,51,54,57,60,63,66,69),], aes(x = Time, y = Type)) +
        geom_tile(aes(fill = logFC), color="grey") + 
        geom_text(aes(label = round(logFC, 1)),size=4.5)+
        scale_fill_gradient2(low = "blue", mid="white", high = "red", midpoint=0,limits=c(-4,4),oob=squish)+
        xlab("Time of day (ZT)")+
        ylab(expression("log"[2]*" FC"))+
        scale_x_discrete(expand = c(0, 0))+
        scale_y_discrete(expand = expansion(mult = c(0.65, 0.65),add = -1))+
        #ggtitle("Total mRNAs")+
        theme_minimal()+
        theme(plot.title=element_text(hjust=0.5),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=13, color="black"),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_line(colour = "black", size = 1),
              axis.title.y = element_text(angle = 0, vjust=0.5),
              axis.title.x = element_text(angle = 0, vjust=0.5, size=13),
              panel.grid.major = element_blank(),
              axis.title = element_text(size = 14),
              legend.position = "none",
              plot.margin = unit(c(-1,-0.2,0.5,0.1), "cm"))
      
      g2<-ggplot(data=expB[expB$AGI == toupper(input$AGI)&expB$Type == "TRAP"&expB$Time%in% c(48,51,54,57,60,63,66,69),], aes(x = Time, y = Type)) +
        geom_tile(aes(fill = logFC), color="grey") + 
        geom_text(aes(label = round(logFC, 1)),size=4.5)+
        scale_fill_gradient2(low = "blue", mid="white", high = "red", midpoint=0,limits=c(-4,4),oob=squish)+
        xlab("Time of day (ZT)")+
        ylab("")+
        labs(fill = expression("log"[2]*" FC"))+
        scale_x_discrete(expand = c(0, 0))+
        scale_y_discrete(expand = expansion(mult = c(0.65, 0.65),add = -1))+
        #ggtitle("TRAP mRNAs")+
        theme_minimal()+
        theme(plot.title=element_text(hjust=0.5),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=13, color="black"),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_line(colour = "black", size = 1),
              axis.title.y = element_text(colour = "transparent", angle = 0, vjust=0.5),
              axis.title.x = element_text(angle = 0, vjust=0.5),
              panel.grid.major = element_blank(),
              #legend.position = "none",
              axis.title = element_text(size = 14),
              legend.title=element_text(size=13), 
              legend.text=element_text(size=13),
              legend.text.align = 1,
              plot.margin = unit(c(-1,0.9,0.5,0.9), "cm"))
      
      lg1<-cowplot::get_legend(g2)
      lg1<-as_ggplot(lg1)
      
      g2<-g2+theme(legend.position = "none")
      
      g3<-ggplot(data=expB[expB$AGI == toupper(input$AGI)&expB$Type == "Total"&expB$Time%in% c(48,51,54,57,60,63,66,69),], aes(x = Time, y = Type)) +
        geom_tile(aes(fill = FDR <0.05), color="transparent") + 
        scale_fill_manual(values = setNames(c("transparent", "transparent"),c(F, T)))+
        geom_point(shape=16,size=5,show.legend =T,aes(colour=FDR<0.05))+
        scale_color_manual(values = setNames(c("lightgrey", "limegreen"),c(F, T)))+
        #xlab("Time of day (ZT)")+
        ylab("FDR")+
        scale_x_discrete(expand = c(0, 0))+
        scale_y_discrete(expand = expansion(mult = c(1, 1),add = -1))+
        #ggtitle("TRAP mRNAs")+
        theme_minimal()+
        theme(plot.title=element_text(hjust=0.5),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_text(color="black",angle = 0, vjust=0.5),
              axis.title.x = element_blank(),
              panel.grid.major = element_blank(),
              axis.title = element_text(size = 14),
              plot.margin = unit(c(-0.7,-0.2,0,0.8), "cm"),
              legend.position = "none")
      
      g4<-ggplot(data=expB[expB$AGI == toupper(input$AGI)&expB$Type == "TRAP"&expB$Time%in% c(48,51,54,57,60,63,66,69),], aes(x = Time, y = Type)) +
        geom_tile(aes(fill = FDR <0.05), color="transparent") + 
        scale_fill_manual(values = setNames(c("transparent", "transparent"),c(F, T)))+
        geom_point(shape=16,size=5,show.legend =T,aes(colour=FDR<0.05))+
        scale_color_manual(values = setNames(c("lightgrey", "limegreen"),c(F, T)))+
        #xlab("Time of day (ZT)")+
        ylab("")+
        scale_x_discrete(expand = c(0, 0))+
        scale_y_discrete(expand = expansion(mult = c(1, 1),add = -1))+
        #ggtitle("TRAP mRNAs")+
        theme_minimal()+
        theme(plot.title=element_text(hjust=0.5),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_text(color="white",angle = 0, vjust=0.5),
              axis.title.x = element_blank(),
              panel.grid.major = element_blank(),
              plot.margin = unit(c(-0.7,1,0,0.9), "cm"))
      
      #lg2<-cowplot::get_legend(g2)
      #lg2<-as_ggplot(lg2)
      
      g4<-g4+theme(legend.position = "none")
      
      gx<-ggplot(data=rec[rec$AGI == "AT1G01010"&rec$Type == "TRAP",], aes(x = Time, y = Type)) +
        geom_tile(aes(fill = FDR <0.05), color="white") + 
        scale_fill_manual(values = setNames(c("white", "white"),c(F, T)))+
        geom_point(shape=16,size=5,show.legend =T,aes(colour=FDR<0.05))+
        scale_color_manual(values = setNames(c("lightgrey", "limegreen"),c(F, T)))+
        theme(legend.title=element_text(size=13),
              legend.text=element_text(size=13),
              legend.key = element_rect(colour = "white", fill = NA))
      
      lg2<-cowplot::get_legend(gx)
      lg2<-as_ggplot(lg2)
      lg2
      
      lay <- rbind(c(1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                   c(1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  6,6,6,6,6,6, 6,6,6,6),
                   c(1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
                   c(1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  1,1,1,1,1,1,1,1,1,1,  8,8,8,8,8,8,8,8,8,8),
                   c(NA,2,2,2,2,2,2,2,2,2, 2,2,2,2,2,2,2,2,2,2,  2,2,2,2,NA,  3,3,3,3,3,  3,3,3,3,3,3,3,3,3,3,  3,3,3,3,3,3,3,3,3,NA,  7,7,7,7,7,7,7,7,7,7),
                   c(NA,4,4,4,4,4,4,4,4,4,  4,4,4,4,4,4,4,4,4,4, 4,4,4,4,NA,  5,5,5,5,5,  5,5,5,5,5,5,5,5,5,5,  5,5,5,5,5,5,5,5,5,NA,  7,7,7,7,7,7,7,7,7,7))
      
      ptlist <- list(g0,g3,g4,g1,g2,lg,lg1,lg2)
      grid.arrange(grobs=ptlist,layout_matrix = lay)
      
    })
    
  })
  
  
  output$Two <- renderPlot({print(Nb_Two())})
  
  output$Legend_single_heat <-renderText({
    paste("<p style='text-align:justify'>","The grey areas represent the subjective night. Data are means +/- sd for n = 3 biological replicates. 
    FDR and Log2 Fold Change (log2 FC) values correspond to the result of a differential expression analysis, comparing 37\u00B0C (orange triangles) to 22\u00B0C (grey circles). 
    Colored (green) dots indicate if the FDR < 0.05, and the color scale (blue to red) represents log2 FC values.
          Upregulation and downregulation are represented in red and blue, respectively.
          ",
          "</p>"
    )
  })
  
  output$Sing_Heat_stressPNG <- downloadHandler(
    filename = "Heat_stress.png.png",
    content = function(file){
      ggsave(file, plot = Nb_Two(), width = 9, height = 5, device = "png")
    })
  
  output$Sing_Heat_stressSVG <- downloadHandler(
    filename = "Heat_stress.svg.svg",
    content = function(file){
      ggsave(file, plot = Nb_Two(), width = 9, height = 5, device = "svg")
    })
  
  output$Sing_Heat_stressPDF <- downloadHandler(
    filename = "Heat_stress.pdf.pdf",
    content = function(file){
      ggsave(file, plot = Nb_Two(), width = 9, height = 5, device = "pdf")
    })
  
  # Reactive value for selected dataset
  datasetInputTwo <- reactive({
    expB[expB$AGI == input$AGI,]
  })
  
  # Downloadable txt of selected dataset
  output$Sing_Data_Heat_stress <- downloadHandler(
    filename = "Data_Heat_stress.txt.txt",
    content = function(file) {
      write.table(datasetInputTwo(), file, row.names = FALSE)
    }
  )
  
  #=========== Tab 1 Box 3 - Recovery following heat stress ("Three") ======
  
  Nb_Three <- eventReactive(input$Submit_Tab1,{
    
    validate(
      need(input$AGI != "", "Provide a gene ID before clicking on the Submit button")
    )
    req(input$AGI)
    
    Single_gene_data_ctrl<-Single_gene_data_ctrl()
    
    validate(
      need(nrow(Single_gene_data_ctrl)>0, "This gene is miswritten (check instructions) or below the threshold of detection based on our criteria")
    )
    
    withProgress(message = 'Making plot 3', value = 1, {
      
      
      
      g5<-ggplot(data=rec[rec$AGI == toupper(input$AGI),], aes(x = Time, y = Mean, group= Temperature)) +
        geom_rect(aes(xmin=6, xmax=7, ymin=-Inf, ymax=Inf), 
                  fill="grey", alpha=0.05, color = NA)+
        geom_point(aes(shape = Temperature, colour = Temperature),
                   position=position_dodge(width=0.1), size = 4)+
        geom_line(position=position_dodge(width=0.1),
                  aes(linetype = Temperature, colour = Temperature))+
        scale_linetype_manual(values=c("solid", "blank"),labels=c("22\u00B0C","37\u00B0C"))+
        scale_color_manual(values=c("#666666","#FF6600"),labels=c("22\u00B0C","37\u00B0C"))+
        scale_shape_manual(values=c(16,17),labels=c("22\u00B0C","37\u00B0C"))+
        geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd, colour = Temperature),width=.5, size=0.8,
                      position=position_dodge(width=0.1),show.legend=F)+
        scale_x_continuous(breaks = seq(0,7,1))+
        facet_wrap(~Type,labeller=as_labeller(group.labs))+
        theme_bw()+
        xlab("Time of recovery (h)\n0h = ZT54")+
        ylab("Transcript abundance\n(rlog normalized counts)")+
        theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
              axis.text=element_text(size=13, colour = "black"),
              axis.text.x=element_text(margin = unit(c(3, 0, 0, 0), "mm")),
              axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
              strip.text = element_text(size = 15),
              axis.title = element_text(size = 14),
              legend.title=element_text(size=13), 
              legend.text=element_text(size=13),
              plot.title = element_text(face="bold", size=18,hjust=-0.14,vjust=0),
              strip.background =element_rect(colour=NA, fill="white"),
              strip.text.x = element_text(margin = margin(0, 0, 0.5, 0, "cm"),face="bold", color="black",size=16))+
        guides(colour=guide_legend(""),
               shape = guide_legend(""),
               linetype = guide_legend(""))
      
      lg3<- cowplot::get_legend(g5)
      lg3<-as_ggplot(lg3)
      g5<-g5+theme(legend.position = "none")
      
      
      g6<-ggplot(data=rec[rec$AGI == toupper(input$AGI)&rec$Type == "Total"&rec$Time%in% c(0,1,3,6),], aes(x = factor(Time), y = Type)) +
        geom_tile(aes(fill = logFC), color="grey") + 
        geom_text(aes(label = round(logFC, 1)),size=4.5)+
        scale_fill_gradient2(low = "blue", mid="white", high = "red", midpoint=0,limits=c(-4,4),oob=squish)+
        xlab("Time of recovery (h)\n0h = ZT54")+
        ylab(expression("log"[2]*" FC"))+
        scale_x_discrete(expand = c(0, 0))+
        scale_y_discrete(expand = expansion(mult = c(0.65, 0.65),add = -1))+
        theme_minimal()+
        theme(plot.title=element_text(hjust=0.5),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=13, color="black"),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_line(colour = "black", size = 1),
              axis.title.y = element_text(angle = 0, vjust=0.5),
              axis.title.x = element_text(angle = 0, vjust=0.5),
              axis.title = element_text(size = 14),
              panel.grid.major = element_blank(),
              legend.position = "none",
              #axis.ticks = element_line(colour = "black", size = 1),
              plot.margin = unit(c(-1,0.5,0,0.9), "cm"))
      
      
      
      
      g7<-ggplot(data=rec[rec$AGI == toupper(input$AGI)&rec$Type == "TRAP"&rec$Time%in% c(0,1,3,6),], aes(x = factor(Time), y = Type)) +
        geom_tile(aes(fill = logFC), color="grey") + 
        geom_text(aes(label = round(logFC, 1)),size=4.5)+
        scale_fill_gradient2(low = "blue", mid="white", high = "red", midpoint=0,limits=c(-4,4),oob=squish)+
        xlab("Time of recovery (h)\n0h = ZT54")+
        ylab(expression("log"[2]*" FC"))+
        labs(fill = expression("log"[2]*" FC"))+
        scale_x_discrete(expand = c(0, 0))+
        scale_y_discrete(expand = expansion(mult = c(0.65, 0.65),add = -1))+
        theme_minimal()+
        theme(plot.title=element_text(hjust=0.5),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=13, color="black"),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_line(colour = "black", size = 1),
              axis.title.y = element_text(angle = 0, vjust=0.5),
              axis.title.x = element_text(angle = 0, vjust=0.5),
              panel.grid.major = element_blank(),
              axis.title = element_text(size = 14),
              #legend.position = "none",
              legend.title=element_text(size=13), 
              legend.text=element_text(size=13),
              legend.text.align = 1,
              plot.margin = unit(c(-1,0.5,0,0.9), "cm"))
      
      lg4<-cowplot::get_legend(g7)
      lg4<-as_ggplot(lg4)
      
      g7<-g7+theme(legend.position = "none")
      
      g8<-ggplot(data=rec[rec$AGI == toupper(input$AGI)&rec$Type == "Total"&rec$Time%in% c(0,1,3,6),], aes(x = factor(Time), y = Type)) +
        geom_tile(aes(fill = FDR <0.05), color="transparent") + 
        scale_fill_manual(values = setNames(c("transparent", "transparent"),c(F, T)))+
        geom_point(shape=16,size=5,show.legend =T,aes(colour=FDR<0.05))+
        scale_color_manual(values = setNames(c("lightgrey", "limegreen"),c(F, T)))+
        ylab("FDR")+
        scale_x_discrete(expand = c(0, 0))+
        scale_y_discrete(expand = expansion(mult = c(1, 1),add = -1))+
        theme_minimal()+
        theme(plot.title=element_text(hjust=0.5),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_text(color="black",angle = 0, vjust=0.5),
              axis.title.x = element_blank(),
              panel.grid.major = element_blank(),
              axis.title = element_text(size = 14),
              plot.margin = unit(c(-1.2,0.5,-0.5,1.6), "cm"),
              legend.position = "none")
      
      g9<-ggplot(data=rec[rec$AGI == toupper(input$AGI)&rec$Type == "TRAP"&rec$Time%in% c(0,1,3,6),], aes(x = factor(Time), y = Type)) +
        geom_tile(aes(fill = FDR <0.05), color="transparent") + 
        scale_fill_manual(values = setNames(c("transparent", "transparent"),c(F, T)))+
        geom_point(shape=16,size=5,show.legend =T,aes(colour=FDR<0.05))+
        scale_color_manual(values = setNames(c("lightgrey", "limegreen"),c(F, T)))+
        ylab("FDR")+
        scale_x_discrete(expand = c(0, 0))+
        scale_y_discrete(expand = expansion(mult = c(1, 1),add = -1))+
        theme_minimal()+
        theme(plot.title=element_text(hjust=0.5),
              axis.text.y = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_text(color="black",angle = 0, vjust=0.5),
              axis.title.x = element_blank(),
              panel.grid.major = element_blank(),
              axis.title = element_text(size = 14),
              plot.margin = unit(c(-1.2,0.5,-0.5,1.6), "cm"),
              legend.position = "none")
      
      gxy<-ggplot(data=rec[rec$AGI == "AT1G01010"&rec$Type == "TRAP",], aes(x = Time, y = Type)) +
        geom_tile(aes(fill = FDR <0.05), color="white") + 
        scale_fill_manual(values = setNames(c("white", "white"),c(F, T)))+
        geom_point(shape=16,size=5,show.legend =T,aes(colour=FDR<0.05))+
        scale_color_manual(values = setNames(c("lightgrey", "limegreen"),c(F, T)))+
        theme(legend.title=element_text(size=13),
              legend.text=element_text(size=13),
              legend.key = element_rect(colour = "white", fill = NA))
      
      lg5<-cowplot::get_legend(gxy)
      lg5<-as_ggplot(lg5)
      lg5 
      
      
      
      
      lay <- rbind(c(1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  NA,NA,NA,NA,NA),
                   c(1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  6,6,6,6,6),
                   c(1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  NA,NA,NA,NA,NA),
                   c(1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  1,1,1,1,1,  8,8,8,8,8),
                   c(NA,2,2,2,2, 2,2,2,2,2, 2,2,3,3, 3,3,3,3,3,  3,3,3,3,NA,NA, 7,7,7,7,7),
                   c(NA,4,4,4,4, 4,4,4,4,4, 4,4,5,5, 5,5,5,5,5,  5,5,5,5,NA,NA, 7,7,7,7,7))
      
      
      ptlist <- list(g5,g8,g9,g6,g7,lg3,lg4,lg5)
      grid.arrange(grobs=ptlist,layout_matrix = lay)
      
      
    })
    
  })
  
  output$Three <- renderPlot({print(Nb_Three())})
  
  output$Legend_single_recovery <-renderText({
    paste("<p style='text-align:justify'>","The grey areas correspond to the subjective night. Data are means +/- sd for n = 3 biological replicates. 
    FDR and Log2 Fold Change (log2 FC) values correspond to the result of a differential expression analysis, comparing the stress (orange triangles) vs control (grey circles) conditions. 
    Colored (green) dots indicate if the FDR < 0.05, and the color scale (blue to red) represents log2 FC values.
          Upregulation and downregulation are represented in red and blue, respectively.",
          "</p>"
    )
  })
  
  output$Sing_RecoveryPNG <- downloadHandler(
    filename = "Recovery.png.png",
    content = function(file){
      ggsave(file, plot = Nb_Three(), width = 7, height = 5, device = "png")
    })
  
  output$Sing_RecoverySVG <- downloadHandler(
    filename = "Recovery.svg.svg",
    content = function(file){
      ggsave(file, plot = Nb_Three(), width = 7, height = 5, device = "svg")
    })
  
  output$Sing_RecoveryPDF <- downloadHandler(
    filename = "Recovery.pdf.pdf",
    content = function(file){
      ggsave(file, plot = Nb_Three(), width = 8, height = 5, device = "pdf")
    })
  
  # Reactive value for selected dataset
  datasetInputThree <- reactive({
    rec[rec$AGI == input$AGI,]
  })
  
  # Downloadable csv of selected dataset
  output$Sing_Data_Recovery <- downloadHandler(
    filename = "Data_Recovery.txt.txt",
    content = function(file) {
      write.table(datasetInputThree(), file, row.names = FALSE)
    }
  )
  
  
  
  ###========= Server Tab 2 "Multiple genes" ============
  #=============================================#
  output$Tab2_description <-renderText({# Multiple genes
    paste("<p style='text-align:justify'>", 'Similar to the "Single genes" tab, this tab generates plots from the Arabidopsis transcriptome and translatome datasets published in <a href="https://doi.org/10.1093/plcell/koab113">Bonnot and Nagel (2021)</a>.
          You need to either click on:',"<br>",
          '- "choose a transcription factor family" (1.) then pick a TF family (2.). The TF families were defined in <a href="https://doi.org/10.1016/j.celrep.2014.06.033">Pruneda-Paz <i>et al.,</i> (2014)</a>.',
          "<br>", '- "paste a list of genes" (1.) then paste a list of Arabidopsis AGI numbers (2.) in the indicated area.',
          "<br>", 'You can choose to sort the genes by (3.) their phase (timing of peak abundance), their average heat stress response over the day, their average heat stress recovery (between 1 h and 6 h of recovery following heat stress), either in the transcriptome or translatome data, or by input order (if a list of genes has been pasted above).',
          "</p>","<p style='text-align:justify'>",
          'All data were obtained from 12 days old seedlings (whole seedlings) that were grown in light (12 h) and dark cycles (12 h) at 22\u00B0C for 10 days and then transferred to free-running conditions (continuous light and temperature) for two days before sampling. 
          Time of day is referred to as Zeitgeber Time (ZT) and corresponds to the hours after moving the seedlings in constant conditions (light and temperature). At each time point, total mRNAs and ribosome-associated mRNAs (using Translating Ribosome Affinity Purification, TRAP) were isolated, and mRNA-Seq (transcriptome) and TRAP-Seq (translatome) was performed.',
          "</p>","<p style='text-align:justify'>",
          '- <b>Time course at 22\u00B0C</b>: A 24 h time course (with samples every 3 h) was analyzed under control conditions at 22\u00B0C. A detection of circadian oscillations was performed using the R package <a href="https://doi.org/10.1093/bioinformatics/btw405">"Metacycle"</a>. 
          This analysis calculated a "Phase" and an "Adjusted p-value" for each transcript. In <a href="https://doi.org/10.1093/plcell/koab113">Bonnot and Nagel (2021)</a>, oscillations were considered significant for genes: with an Adjusted p-value < 0.01; with an 0.01 < Adjusted p-value < 0.05 and that overlapped with previously published lists of circadian transcripts. 
          For details, please see the Methods section in <a href="https://doi.org/10.1093/plcell/koab113">Bonnot and Nagel (2021).</a>',
          "<br>", '- <b>Heat stress response</b>: Heat stress treatments were applied at different times of day (from ZT48, subjective dawn, to ZT69, end of subjective night) on separate sets of plants, and correspond to treatments at 37\u00B0C for 1 h.
          To compare the stress vs control conditions, a differential expression analysis was performed using the R package <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8">"DESeq2"</a>, and FDR and Log2 Fold Change (log2 FC) values (37\u00B0C vs 22\u00B0C) were calculated.
          In <a href="https://doi.org/10.1093/plcell/koab113">Bonnot and Nagel (2021)</a>, genes were considered as significantly differentially regulated when the FDR < 0.05 and the log2 FC > |1|.',
          "<br>", '- <b>Recovery following heat stress</b>: Heat stress treatment at 37\u00B0C for 1 h was applied from ZT53 to ZT54, corresponding to the middle of the light period. 
          Following heat stress, seedlings were either sampled (0 h of recovery), or moved back to 22\u00B0C for recovery. To compare the stress (or recovery) vs control conditions, a differential expression analysis was performed using the R package <a href="https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8">"DESeq2"</a>, and FDR and Log2 Fold Change (log2 FC) values were calculated.
          In <a href="https://doi.org/10.1093/plcell/koab113">Bonnot and Nagel (2021)</a>, genes were considered as significantly differentially regulated when the FDR < 0.05 and the log2 FC > |1|.',
          "</p>","<p style='text-align:justify'>"
        
    )
  })
  
  output$beginTab2 <- renderText({
    validate(
      need(!is.null(input$List_gen), "Set parameters and click on the Submit button")
    )
  })
  
  
  
  #=========== Tab 2 Box 1 - Time course at 22 ("Four")=====
  
  
  min <- c(1,34,44,54,64,74,84,94,104,114)
  max <- c(35,45,55,65,75,85,95,105,115,26000)
  font <- c(10,9,8,7,6,5,4,3,2,1)
  font_range <- data.frame(min, max, font)
  
  
  rlog.TOT <- eventReactive(input$Submit_Tab2, {
    req(input$List_gen)
    if (input$List_gen == "fam") {
      
      validate(
        need(input$TFfamil != "", "Select a TF family before clicking on the Submit button")
      )
      
      list<-input$TFfamil
      #nb_genes<-length(list)
      ID<-filter(TFfamilylist, Family %in% list)
      data.rlog.TOT<-filter(rlog.mean.TOT, AGI %in% ID$AGI)
      
      nb_genes <- nrow(data.rlog.TOT)
      fontsize <- font_range[font_range$min < nb_genes & font_range$max > nb_genes,3]
      data.rlog.TOT$Font <- fontsize
      
      if (input$Sort_by == "order_by_transcriptome_ph") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$Phase),]
      }
      else if (input$Sort_by == "order_by_translatome_ph") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$Phase_TRAP),]
      }
      else if (input$Sort_by == "Heat_transcriptome") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$HS_TOT),]
      }  
      else if (input$Sort_by == "Heat_translatome") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$HS_TRAP),]
      } 
      else if (input$Sort_by == "Recovery_transcriptome") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$Recovery_TOT),]
      }  
      else if (input$Sort_by == "Recovery_translatome") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$Recovery_TRAP),]
      } 
      
    } else if (input$List_gen == "paste") {
      
      validate(
        need(input$Genes_list_Tab2!= "", "Provide a list of genes before clicking on the Submit button")
      )
      
      list<-unlist(strsplit(toupper(input$Genes_list_Tab2), "\n"))
      #nb_genes<-length(list)
      data.rlog.TOT<-filter(rlog.mean.TOT, AGI %in% list)
      
      nb_genes <- nrow(data.rlog.TOT)
      fontsize <- font_range[font_range$min < nb_genes & font_range$max > nb_genes,3]
      data.rlog.TOT$Font <- fontsize
      
      if (input$Sort_by1 == "order_by_transcriptome_ph") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$Phase),]
      }
      else if (input$Sort_by1 == "order_by_translatome_ph") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$Phase_TRAP),]
      }
      else if (input$Sort_by1 == "Heat_transcriptome") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$HS_TOT),]
      }  
      else if (input$Sort_by1 == "Heat_translatome") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$HS_TRAP),]
      }  
      else if (input$Sort_by1 == "Recovery_transcriptome") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$Recovery_TOT),]
      }  
      else if (input$Sort_by1 == "Recovery_translatome") {
        data.rlog.TOT <- data.rlog.TOT[order(data.rlog.TOT$Recovery_TRAP),]
      } 
      else if (input$Sort_by1 == "order_by_input") {
        list<-as.data.frame(unlist(strsplit(toupper(input$Genes_list_Tab2), "\n")))
        names(list)<-c("AGI")
        data.rlog.TOT<-merge(list, data.rlog.TOT, by="AGI", sort=F)
      }
    }
  })
  
  rlog.TRAP <- eventReactive(input$Submit_Tab2, {
    req(input$List_gen)
    if (input$List_gen == "fam") {
      
      list<-input$TFfamil
      #nb_genes<-length(list)
      ID<-filter(TFfamilylist, Family %in% list)
      data.rlog.TRAP<-filter(rlog.mean.TRAP, AGI %in% ID$AGI)
      
      if (input$Sort_by == "order_by_transcriptome_ph") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$Phase_TOT),]
      }
      else if (input$Sort_by == "order_by_translatome_ph") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$Phase),]
      }
      else if (input$Sort_by == "Heat_transcriptome") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$HS_TOT),]
      }
      else if (input$Sort_by == "Heat_translatome") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$HS_TRAP),]
      }
      else if (input$Sort_by == "Recovery_transcriptome") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$Recovery_TOT),]
      }  
      else if (input$Sort_by == "Recovery_translatome") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$Recovery_TRAP),]
      } 
      
      
    } else if (input$List_gen == "paste") {
      
      list<-unlist(strsplit(toupper(input$Genes_list_Tab2), "\n"))
      
      #validate(
      #  need(nrow(list)<=2, "Please provide more than one gene")
      #)
      
      #nb_genes<-length(list)
      data.rlog.TRAP<-filter(rlog.mean.TRAP, AGI %in% list)
      
      if (input$Sort_by1 == "order_by_transcriptome_ph") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$Phase_TOT),]
      }
      else if (input$Sort_by1 == "order_by_translatome_ph") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$Phase),]
      }
      else if (input$Sort_by1 == "Heat_transcriptome") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$HS_TOT),]
      }
      else if (input$Sort_by1 == "Heat_translatome") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$HS_TRAP),]
      }
      else if (input$Sort_by1 == "Recovery_transcriptome") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$Recovery_TOT),]
      }  
      else if (input$Sort_by1 == "Recovery_translatome") {
        data.rlog.TRAP <- data.rlog.TRAP[order(data.rlog.TRAP$Recovery_TRAP),]
      } 
      else if (input$Sort_by1 == "order_by_input") {
        list<-as.data.frame(unlist(strsplit(toupper(input$Genes_list_Tab2), "\n")))
        names(list)<-c("AGI")
        data.rlog.TRAP<-merge(list, data.rlog.TRAP, by="AGI", sort=F)
      }
      
    }
  })
  
  Nb_Four <- eventReactive(input$Submit_Tab2, {
    validate(
      need(input$List_gen != "", "Provide a list of genes before clicking on the Submit button")
    )
    req(input$List_gen)
    
    Gene.list.rlog.TOT<-rlog.TOT() 
    
    validate(
      need(nrow(Gene.list.rlog.TOT)>0, "These genes are miswritten (check instructions), or below the threshold of detection based on our criteria")
    )
    
    row.names(Gene.list.rlog.TOT)<-Gene.list.rlog.TOT[,1]
    annotation_row.TOT <- Gene.list.rlog.TOT[,c(12,18)]
    names(annotation_row.TOT)[2] <- "Phase"
    col.index1 <- which(colnames(Gene.list.rlog.TOT) == "Font") 
    fontsize <- Gene.list.rlog.TOT[1,col.index1]
    
    withProgress(message = "Making plot 1", value = 1, {
    
    g9 <- as.ggplot(pheatmap(Gene.list.rlog.TOT[,2:10], color=myColor,
                             cluster_cols=F, cluster_rows = F, show_rownames = T, breaks = myBreaks, scale = 'row',
                             annotation_row = annotation_row.TOT, annotation_colors = annotation_colors,
                             labels_col  = seq(48,72,3),angle_col = 45, main = "Transcriptome", fontsize_row = fontsize))
    
    
    Gene.list.rlog.TRAP<-rlog.TRAP() 
    row.names(Gene.list.rlog.TRAP)<-Gene.list.rlog.TRAP[,1]
    annotation_row.TRAP <- Gene.list.rlog.TRAP[,c(12,18)]
    names(annotation_row.TRAP)[2] <- "Phase"
    
    g10 <- as.ggplot(pheatmap(Gene.list.rlog.TRAP[,2:10], color=myColor, 
                              cluster_cols=F, cluster_rows = F, show_rownames = T, breaks = myBreaks, scale = 'row',
                              annotation_row = annotation_row.TRAP, annotation_colors = annotation_colors,
                              labels_col  = seq(48,72,3),angle_col = 45, main = "Translatome", fontsize_row = fontsize))
    
    lay <- rbind(c(NA,1,1,1,1,1,1,1,1,1,1,1,1,
                   2,2,2,2,2,2,2,2,2,2,2,2))
    
    ptlist <- list(g9,g10)
    grid.arrange(grobs=ptlist,layout_matrix = lay)
    })
  })
  
  
  output$Four <- renderPlot({
    print(Nb_Four())
  })
  
  # output$mytable = renderDataTable({
  #   # Sel_ref()
  #   rlog.TOT()
  # })
  
  output$Mult_Time_course_22PNG <- downloadHandler(
    filename = "Mult_Time_course_22.png.png",
    content = function(file){
      ggsave(file, plot = Nb_Four(), width = 10, height = 5, device = "png")
    })
  
  output$Mult_Time_course_22SVG <- downloadHandler(
    filename = "Mult_Time_course_22.svg.svg",
    content = function(file){
      ggsave(file, plot = Nb_Four(), width = 10, height = 5, device = "svg")
    })
  
  output$Mult_Time_course_22PDF <- downloadHandler(
    filename = "Mult_Time_course_22.pdf.pdf",
    content = function(file){
      ggsave(file, plot = Nb_Four(), width = 10, height = 5, device = "pdf")
    })
  
  output$Legend_heatmap_timecourse <-renderText({
    paste("<p style='text-align:justify'>",
    "The column 'Phase' defines the timing of peak abundance (a phase of 0 and 12 indicates a peak abundance at subjective dawn and the beginning of subjective night, respectively). 
    In the column 'Cycling', 'Yes' or 'No' indicate if the circadian oscillation was detected as significant (black) or not (white) in <a href='https://doi.org/10.1093/plcell/koab113'>Bonnot and Nagel (2021)</a>, respectively. 
    Of note, because the phase is calculated as the time point of highest transcript abundance, a gene can have a phase but not be significantly cycling.
    In columns '48' to '72', the color scale from blue to yellow reflects the normalized transcript abundance (rlog normalized counts). High and low transcript abundance are represented in yellow and purple, respectively.
    Times 48 - 60 and 60 - 72 correspond to the day and the subjective night, respectively.",
          "</p>"
          
    )
  })
  
  # Downloadable txt of selected dataset
  output$Mult_DataTOT_Time_course_22 <- downloadHandler(
    filename = "DataTOT_Time_course_22.txt.txt",
    content = function(file) {
      write.table(rlog.TOT(), file, row.names = FALSE)
    }
  )
  
  output$Mult_DataTRAP_Time_course_22 <- downloadHandler(
    filename = "DataTRAP_Time_course_22.txt.txt",
    content = function(file) {
      write.table(rlog.TRAP(), file, row.names = FALSE)
    }
  )
  
  
  #=========== Tab 2 Box 2 - Heat stress response ("Five") =====
  
  LFC.TOT <- eventReactive(input$Submit_Tab2, {
    req(input$List_gen)
    if (input$List_gen == "fam") {
      
      validate(
        need(input$TFfamil != "", "Select a TF family before clicking on the Submit button")
      )
      
      list<-input$TFfamil
      ID<-filter(TFfamilylist, Family %in% list)
      data.LFC.TOT<-filter(data.LFC.TOT, AGI %in% ID$AGI)
      
      nb_genes <- nrow(data.LFC.TOT)
      fontsize <- font_range[font_range$min < nb_genes & font_range$max > nb_genes,3]
      data.LFC.TOT$Font <- fontsize
      
      if (input$Sort_by == "order_by_transcriptome_ph") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$Phase),]
      }
      else if (input$Sort_by == "order_by_translatome_ph") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$Phase_TRAP),]
      }
      else if (input$Sort_by == "Heat_transcriptome") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$HS_TOT),]
      }  
      else if (input$Sort_by == "Heat_translatome") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$HS_TRAP),]
      } 
      else if (input$Sort_by == "Recovery_transcriptome") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$Recovery_TOT),]
      }  
      else if (input$Sort_by == "Recovery_translatome") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$Recovery_TRAP),]
      } 
      
    } else if (input$List_gen == "paste") {
      
      validate(
        need(input$Genes_list_Tab2!= "", "Provide a list of genes before clicking on the Submit button")
      )
      
      list<-unlist(strsplit(toupper(input$Genes_list_Tab2), "\n"))
      data.LFC.TOT<-filter(data.LFC.TOT, AGI %in% list)
      
      nb_genes <- nrow(data.LFC.TOT)
      fontsize <- font_range[font_range$min < nb_genes & font_range$max > nb_genes,3]
      data.LFC.TOT$Font <- fontsize
      
      if (input$Sort_by1 == "order_by_transcriptome_ph") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$Phase),]
      }
      else if (input$Sort_by1 == "order_by_translatome_ph") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$Phase_TRAP),]
      }
      else if (input$Sort_by1 == "Heat_transcriptome") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$HS_TOT),]
      }  
      else if (input$Sort_by1 == "Heat_translatome") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$HS_TRAP),]
      } 
      else if (input$Sort_by1 == "Recovery_transcriptome") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$Recovery_TOT),]
      }  
      else if (input$Sort_by1 == "Recovery_translatome") {
        data.LFC.TOT <- data.LFC.TOT[order(data.LFC.TOT$Recovery_TRAP),]
      } 
      else if (input$Sort_by1 == "order_by_input") {
        list<-as.data.frame(unlist(strsplit(toupper(input$Genes_list_Tab2), "\n")))
        names(list)<-c("AGI")
        data.LFC.TOT<-merge(list, data.LFC.TOT, by="AGI", sort=F)
      }
    }
  })
  
  LFC.TRAP <- eventReactive(input$Submit_Tab2, {
    req(input$List_gen)
    if (input$List_gen == "fam") {
      
      list<-input$TFfamil
      ID<-filter(TFfamilylist, Family %in% list)
      data.LFC.TRAP<-filter(data.LFC.TRAP, AGI %in% ID$AGI)
      
      if (input$Sort_by == "order_by_transcriptome_ph") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$Phase_TOT),]
      }
      else if (input$Sort_by == "order_by_translatome_ph") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$Phase),]
      }
      else if (input$Sort_by == "Heat_transcriptome") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$HS_TOT),]
      }  
      else if (input$Sort_by == "Heat_translatome") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$HS_TRAP),]
      } 
      else if (input$Sort_by == "Recovery_transcriptome") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$Recovery_TOT),]
      }  
      else if (input$Sort_by == "Recovery_translatome") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$Recovery_TRAP),]
      } 
      
    } else if (input$List_gen == "paste") {
      
      list<-unlist(strsplit(toupper(input$Genes_list_Tab2), "\n"))
      data.LFC.TRAP<-filter(data.LFC.TRAP, AGI %in% list)
      
      if (input$Sort_by1 == "order_by_transcriptome_ph") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$Phase_TOT),]
      }
      else if (input$Sort_by1 == "order_by_translatome_ph") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$Phase),]
      }
      else if (input$Sort_by1 == "Heat_transcriptome") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$HS_TOT),]
      }  
      else if (input$Sort_by1 == "Heat_translatome") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$HS_TRAP),]
      } 
      else if (input$Sort_by1 == "Recovery_transcriptome") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$Recovery_TOT),]
      }  
      else if (input$Sort_by1 == "Recovery_translatome") {
        data.LFC.TRAP <- data.LFC.TRAP[order(data.LFC.TRAP$Recovery_TRAP),]
      } 
      else if (input$Sort_by1 == "order_by_input") {
        list<-as.data.frame(unlist(strsplit(toupper(input$Genes_list_Tab2), "\n")))
        names(list)<-c("AGI")
        data.LFC.TRAP<-merge(list, data.LFC.TRAP, by="AGI", sort=F)
      }
    }
  })
  
  Nb_Five <- eventReactive(input$Submit_Tab2, {
    
    validate(
      need(input$List_gen != "", "Provide a list of genes before clicking on the Submit button")
    )
    req(input$List_gen)
    
    Gene.list.LFC.TOT<-LFC.TOT() 
    
    validate(
      need(nrow(Gene.list.LFC.TOT)>0, "These genes are miswritten (check instructions) or below the threshold of detection based on our criteria")
    )
    
    row.names(Gene.list.LFC.TOT)<-Gene.list.LFC.TOT[,1]
    annotation_row.TOT <- Gene.list.LFC.TOT[,c(11,17)]
    names(annotation_row.TOT)[2] <- "Phase"
    col.index1 <- which(colnames(Gene.list.LFC.TOT) == "Font") 
    fontsize <- Gene.list.LFC.TOT[1,col.index1]
    
    withProgress(message = "Making plot 2", value = 1, {
    
    g11 <- as.ggplot(pheatmap(Gene.list.LFC.TOT[,2:9], color=myColor2,
                              cluster_cols=F, cluster_rows = F, show_rownames = T, breaks = myBreaks2,
                              labels_col  = seq(48,69,3),angle_col = 45,
                              annotation_row = annotation_row.TOT, annotation_colors = annotation_colors,
                              annotation_legend = T, main = "Transcriptome", fontsize_row = fontsize))
    
    
    Gene.list.LFC.TRAP<-LFC.TRAP()
    row.names(Gene.list.LFC.TRAP)<-Gene.list.LFC.TRAP[,1]
    annotation_row.TRAP <- Gene.list.LFC.TRAP[,c(11,17)]
    names(annotation_row.TRAP)[2] <- "Phase"
    
    g12 <- as.ggplot(pheatmap(Gene.list.LFC.TRAP[,2:9], color=myColor2, 
                              cluster_cols=F, cluster_rows = F, show_rownames = T, breaks = myBreaks2,
                              labels_col  = seq(48,69,3),angle_col = 45,
                              annotation_row = annotation_row.TRAP, annotation_colors = annotation_colors, 
                              annotation_legend = T, main = "Translatome", fontsize_row = fontsize))
    
    lay <- rbind(c(NA,1,1,1,1,1,1,1,1,1,1,1,1,
                   2,2,2,2,2,2,2,2,2,2,2,2))
    
    ptlist <- list(g11,g12)
    grid.arrange(grobs=ptlist,layout_matrix = lay)
    })
    
  })
  
  output$Five <- renderPlot({
    print(Nb_Five())
  })
  
  output$Mult_Heat_stressPNG <- downloadHandler(
    filename = "Mult_Heat_stress.png.png",
    content = function(file){
      ggsave(file, plot = Nb_Five(), width = 10, height = 5, device = "png")
    })
  
  output$Mult_Heat_stressSVG <- downloadHandler(
    filename = "Mult_Heat_stress.svg.svg",
    content = function(file){
      ggsave(file, plot = Nb_Five(), width = 10, height = 5, device = "svg")
    })
  
  output$Mult_Heat_stressPDF <- downloadHandler(
    filename = "Mult_Heat_stress.pdf.pdf",
    content = function(file){
      ggsave(file, plot = Nb_Five(), width = 10, height = 5, device = "pdf")
    })
  
  output$Legend_heatmap_heat <-renderText({
    paste("<p style='text-align:justify'>",
    "The column 'Phase' defines the timing of peak abundance (a phase of 0 and 12 indicates a peak abundance at subjective dawn and the beginning of subjective night, respectively). 
    In the column 'Cycling', 'Yes' or 'No' indicate if the circadian oscillation was detected as significant (black) or not (white) in <a href='https://doi.org/10.1093/plcell/koab113'>Bonnot and Nagel (2021)</a>, respectively.
    Of note, because the phase is calculated as the time point of highest transcript abundance, a gene can have a phase but not be significantly cycling.
    In columns '48' to '69', the color scale from blue to red represents the transcript response to heat stress (Log2 Fold Change values, 37\u00B0C vs 22\u00B0C). Blue and red indicate a downregulation and an upregulation, respectively. 
    Times 48 - 60 and 60 - 69 correspond to the day and the subjective night, respectively.",
          "</p>"
          
    )
  })
  
  # Downloadable txt of selected dataset
  output$Mult_DataTOT_Heat_stress <- downloadHandler(
    filename = "DataTOT_Heat_stress.txt.txt",
    content = function(file) {
      write.table(LFC.TOT(), file, row.names = FALSE)
    }
  )
  
  output$Mult_DataTRAP_Heat_stress <- downloadHandler(
    filename = "DataTRAP_Heat_stress.txt.txt",
    content = function(file) {
      write.table(LFC.TRAP(), file, row.names = FALSE)
    }
  )
  
  #=========== Tab 2 Box 3 - Recovery following heat stress ("Six") =====
  
  Recovery.TOT <- eventReactive(input$Submit_Tab2, {
    req(input$List_gen)
    if (input$List_gen == "fam") {
      
      validate(
        need(input$TFfamil != "", "Select a TF family before clicking on the Submit button")
      )
      
      list<-input$TFfamil
      ID<-filter(TFfamilylist, Family %in% list)
      data.Recovery.TOT<-filter(data.Recovery.TOT, AGI %in% ID$AGI)
      
      nb_genes <- nrow(data.Recovery.TOT)
      fontsize <- font_range[font_range$min < nb_genes & font_range$max > nb_genes,3]
      data.Recovery.TOT$Font <- fontsize
      
      if (input$Sort_by == "order_by_transcriptome_ph") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$Phase),]
      }
      else if (input$Sort_by == "order_by_translatome_ph") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$Phase_TRAP),]
      }
      else if (input$Sort_by == "Heat_transcriptome") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$HS_TOT),]
      }  
      else if (input$Sort_by == "Heat_translatome") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$HS_TRAP),]
      } 
      else if (input$Sort_by == "Recovery_transcriptome") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$Recovery_TOT),]
      }  
      else if (input$Sort_by == "Recovery_translatome") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$Recovery_TRAP),]
      } 
      
    } else if (input$List_gen == "paste") {
      
      validate(
        need(input$Genes_list_Tab2!= "", "Provide a list of genes before clicking on the Submit button")
      )
      list<-unlist(strsplit(toupper(input$Genes_list_Tab2), "\n"))
      data.Recovery.TOT<-filter(data.Recovery.TOT, AGI %in% list)
      
      nb_genes <- nrow(data.Recovery.TOT)
      fontsize <- font_range[font_range$min < nb_genes & font_range$max > nb_genes,3]
      data.Recovery.TOT$Font <- fontsize
      
      if (input$Sort_by1 == "order_by_transcriptome_ph") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$Phase),]
      }
      else if (input$Sort_by1 == "order_by_translatome_ph") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$Phase_TRAP),]
      }
      else if (input$Sort_by1 == "Heat_transcriptome") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$HS_TOT),]
      }  
      else if (input$Sort_by1 == "Heat_translatome") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$HS_TRAP),]
      } 
      else if (input$Sort_by1 == "Recovery_transcriptome") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$Recovery_TOT),]
      }  
      else if (input$Sort_by1 == "Recovery_translatome") {
        data.Recovery.TOT <- data.Recovery.TOT[order(data.Recovery.TOT$Recovery_TRAP),]
      } 
      else if (input$Sort_by1 == "order_by_input") {
        list<-as.data.frame(unlist(strsplit(toupper(input$Genes_list_Tab2), "\n")))
        names(list)<-c("AGI")
        data.Recovery.TOT<-merge(list, data.Recovery.TOT, by="AGI", sort=F)
      }
    }
  })
  
  Recovery.TRAP <- eventReactive(input$Submit_Tab2, {
    req(input$List_gen)
    if (input$List_gen == "fam") {
      
      list<-input$TFfamil
      ID<-filter(TFfamilylist, Family %in% list)
      data.Recovery.TRAP<-filter(data.Recovery.TRAP, AGI %in% ID$AGI)
      
      if (input$Sort_by == "order_by_transcriptome_ph") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$Phase_TOT),]
      }
      else if (input$Sort_by == "order_by_translatome_ph") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$Phase),]
      }
      else if (input$Sort_by == "Heat_transcriptome") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$HS_TOT),]
      }  
      else if (input$Sort_by == "Heat_translatome") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$HS_TRAP),]
      } 
      else if (input$Sort_by == "Recovery_transcriptome") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$Recovery_TOT),]
      }  
      else if (input$Sort_by == "Recovery_translatome") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$Recovery_TRAP),]
      } 
      
    } else if (input$List_gen == "paste") {
      
      list<-unlist(strsplit(toupper(input$Genes_list_Tab2), "\n"))
      data.Recovery.TRAP<-filter(data.Recovery.TRAP, AGI %in% list)
      
      if (input$Sort_by1 == "order_by_transcriptome_ph") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$Phase_TOT),]
      }
      else if (input$Sort_by1 == "order_by_translatome_ph") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$Phase),]
      }
      else if (input$Sort_by1 == "Heat_transcriptome") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$HS_TOT),]
      }  
      else if (input$Sort_by1 == "Heat_translatome") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$HS_TRAP),]
      } 
      else if (input$Sort_by1 == "Recovery_transcriptome") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$Recovery_TOT),]
      }  
      else if (input$Sort_by1 == "Recovery_translatome") {
        data.Recovery.TRAP <- data.Recovery.TRAP[order(data.Recovery.TRAP$Recovery_TRAP),]
      } 
      else if (input$Sort_by1 == "order_by_input") {
        list<-as.data.frame(unlist(strsplit(toupper(input$Genes_list_Tab2), "\n")))
        names(list)<-c("AGI")
        data.Recovery.TRAP<-merge(list, data.Recovery.TRAP, by="AGI", sort=F)
      }
    }
  })
  
  Nb_Six <- eventReactive(input$Submit_Tab2, {
    
    validate(
      need(input$List_gen != "", "Provide a list of genes before clicking on the Submit button")
    )
    req(input$List_gen)
    
    Gene.list.Recovery.TOT<-Recovery.TOT()
    validate(
      need(nrow(Gene.list.Recovery.TOT)>0, "These genes are miswritten (check instructions) or below the threshold of detection based on our criteria")
    )
    
    row.names(Gene.list.Recovery.TOT)<-Gene.list.Recovery.TOT[,1]
    annotation_row.TOT <- Gene.list.Recovery.TOT[,c(7,13)]
    names(annotation_row.TOT)[2] <- "Phase"
    col.index1 <- which(colnames(Gene.list.Recovery.TOT) == "Font") 
    fontsize <- Gene.list.Recovery.TOT[1,col.index1]
    
    withProgress(message = "Making plot 3", value = 1, {
    
    g13 <- as.ggplot(pheatmap(Gene.list.Recovery.TOT[,2:5], color=myColor2, 
                              cluster_cols=F, cluster_rows = F, show_rownames = T, breaks = myBreaks2,
                              labels_col  = c(0,1,3,6),angle_col = 45,
                              annotation_row = annotation_row.TOT, annotation_colors = annotation_colors, 
                              annotation_legend = T, main = "Transcriptome", fontsize_row = fontsize))
    
    Gene.list.Recovery.TRAP<-Recovery.TRAP() 
    row.names(Gene.list.Recovery.TRAP)<-Gene.list.Recovery.TRAP[,1]
    annotation_row.TRAP <- Gene.list.Recovery.TRAP[,c(7,13)]
    names(annotation_row.TRAP)[2] <- "Phase"
    
    g14 <- as.ggplot(pheatmap(Gene.list.Recovery.TRAP[,2:5], color=myColor2, 
                              cluster_cols=F, cluster_rows = F, show_rownames = T, breaks = myBreaks2,
                              labels_col  = c(0,1,3,6),angle_col = 45,
                              annotation_row = annotation_row.TRAP, annotation_colors = annotation_colors, 
                              annotation_legend = T, main = "Translatome", fontsize_row = fontsize))
    
    lay <- rbind(c(NA,1,1,1,1,1,1,1,1,1,1,1,1,
                   2,2,2,2,2,2,2,2,2,2,2,2))
    
    ptlist <- list(g13,g14)
    grid.arrange(grobs=ptlist,layout_matrix = lay)
    })
    
  })
  
  output$Six <- renderPlot({
    print(Nb_Six())
  })
  
  output$Mult_RecoveryPNG <- downloadHandler(
    filename = "Recovery.png.png",
    content = function(file){
      ggsave(file, plot = Nb_Six(), width = 8, height = 5, device = "png")
    })
  
  output$Mult_RecoverySVG <- downloadHandler(
    filename = "Recovery.svg.svg",
    content = function(file){
      ggsave(file, plot = Nb_Six(), width = 8, height = 5, device = "svg")
    })
  
  output$Mult_RecoveryPDF <- downloadHandler(
    filename = "Mult_Recovery.pdf.pdf",
    content = function(file){
      ggsave(file, plot = Nb_Six(), width = 10, height = 5, device = "pdf")
    })
  
  output$Legend_heatmap_recovery <-renderText({
    paste("<p style='text-align:justify'>",
          "The column 'Phase' defines the timing of peak abundance (a phase of 0 and 12 indicates a peak abundance at subjective dawn and the beginning of subjective night, respectively). 
    In the column 'Cycling', 'Yes' or 'No' indicate if the circadian oscillation was detected as significant (black) or not (white) in <a href='https://doi.org/10.1093/plcell/koab113'>Bonnot and Nagel (2021)</a>, respectively.
    Of note, because the phase is calculated as the time point of highest transcript abundance, a gene can have a phase but not be significantly cycling.
          In columns '0' to '6', the color scale from blue to red represents the transcript response during the plant recovery following heat stress (Log2 Fold Change values, stress vs control conditions). Blue and red indicate a downregulation and an upregulation, respectively.",
          "</p>"
    )
  })
  
  # Downloadable txt of selected dataset
  output$Mult_DataTOT_Recovery <- downloadHandler(
    filename = "DataTOT_Recovery.txt.txt",
    content = function(file) {
      write.table(Recovery.TOT(), file, row.names = FALSE)
    }
  )
  
  output$Mult_DataTRAP_Recovery <- downloadHandler(
    filename = "DataTRAP_Recovery.txt.txt",
    content = function(file) {
      write.table(Recovery.TRAP(), file, row.names = FALSE)
    }
  )
  
  ###========= Server Tab 3 "Phases" ============
  #============================================#
  output$Instructions_uploadCSV <-renderText({# Phase enrichment
    paste("To use your own reference, provide a .csv file with two columns - Gene names in the first one, the phase when they peak in the second one - with column headers, and with comma as column separator, i.e.:","<br>",
          "AGI, Phase","<br>",
          "AT1G01060, 0","<br>",
          "AT2G31010,	9","<br>",
          "AT2G31410,	21","<br>",
          "AT4G24920, 16.5","<br>",
          "AT3G57090, 12"
    )
  })
  
  output$Tab3_description <-renderText({# Phase enrichment
    paste(
      "<p style='text-align:justify'>",
      "This tab allows you to identify over-represented of under-represented phases (timing of peak abundance) from a list of genes of interest.",
      "<p style='text-align:justify'>",
      "1. First, select a dataset that will be used as the reference. You can either select an existing dataset or your own reference.",
      "</p>","<p style='text-align:justify'>",
      '2. If "Existing reference" has been selected, pick a dataset in the list provided. 
      These datasets correspond to <b>all transcripts that exhibited significant circadian oscillations</b> in the corresponding published datasets:',
      "<br>",
      '- <b>Bonnot and Nagel_Transcriptome_LL_LDHH</b>: 8028 circadian transcripts identified at the transcriptome level in <a href="https://doi.org/10.1093/plcell/koab113"> Bonnot and Nagel (2021)</a>.
      Data were obtained from 12 days old seedlings (whole seedlings) that were grown in light (12 h) and dark (12 h) cycles at 22\u00B0C for 10 days and then transferred to free-running conditions (continuous light and temperature) for two days before sampling.
      Data correspond to mRNA-Seq, performed on a 24 h time course (with samples every 3 h) under control conditions at 22\u00B0C. The detection of circadian oscillations was performed using the R package <a href="https://doi.org/10.1093/bioinformatics/btw405">"Metacycle"</a>.',
      "<br>",
      '- <b>Bonnot and Nagel_Translatome_LL_LDHH</b>: 10657 circadian transcripts identified at the translatome level in <a href="https://doi.org/10.1093/plcell/koab113"> Bonnot and Nagel (2021)</a>.
      Data were obtained from 12 days old seedlings (whole seedlings) that were grown in light (12 h) and dark (12 h) cycles at 22\u00B0C for 10 days and then transferred to free-running conditions (continuous light and temperature) for two days before sampling.
      Data correspond to TRAP-Seq, performed on a 24 h time course (with samples every 3 h) under control conditions at 22\u00B0C. The detection of circadian oscillations was performed using the R package <a href="https://doi.org/10.1093/bioinformatics/btw405">"Metacycle"</a>.',
      "<br>",
      '- <b>Covington and Harmer_LL_LDHH</b>: 7858 circadian transcripts identified at the transcriptome level in <a href="https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050222"> Covington and Harmer (2007)</a>.
      Data were obtained from whole seedlings that were grown in light (12 h) and dark (12 h) cycles at 22\u00B0C for 7 days and then transferred to free-running conditions (continuous light and temperature).
      Data correspond to a microarray analysis, performed from ZT24 (subjective dawn) on day 9 and during 44 h (with samples every 4 h). For this dataset, the phase information was collected from the tool <a href="http://phaser.mocklerlab.org/">"Phaser"</a>, using a correlation cutoff of 0.7.',
      "<br>",
      '- <b>Edwards et al_LL_LDHH</b>: 9940 circadian transcripts identified at the transcriptome level in <a href="https://doi.org/10.1105/tpc.105.038315">Edwards et al (2006)</a>.
      Data were obtained from whole seedlings that were grown in light (12 h) and dark cycles (12 h) at 22\u00B0C for 7 days and then transferred to free-running conditions (continuous light and temperature).
      Data correspond to a microarray analysis, performed from ZT26 (2 h after subjective dawn) on day 9 and during two days (with samples every 4 h). For this dataset, the phase information was collected from the tool <a href="http://phaser.mocklerlab.org/">"Phaser"</a>, using a correlation cutoff of 0.7.',
      "<br>",
      '- <b>Hsu and Harmer_LL_LDHH</b>: 7124 circadian transcripts identified at the transcriptome level in <a href="https://doi.org/10.1371/journal.pone.0049853">Hsu and Harmer (2012)</a>.
      Data were obtained from whole seedlings that were grown in light (12 h) and dark cycles (12 h) at 22\u00B0C for 7 days and then transferred to free-running conditions (continuous light and temperature).
      Data correspond to a microarray analysis, performed on days 11 and 12 (fourth and fifth days of free run) with samples every 4 h. The detection of circadian oscillations was performed using the tools <a href="https://doi.org/10.1177/0748730410379711">"JTK_CYCLE"</a> and <a href="https://doi.org/10.1016/S0076-6879(04)83007-6">"COSOPT"</a>',
      "<br>",
      '- <b>Michael et al_LL_LDHC</b>: 8909 circadian transcripts identified at the transcriptome level in <a href="https://doi.org/10.1371/journal.pgen.0040014">Michael et al. (2008)</a>.
      Data were obtained from whole seedlings that were grown in light/dark cycles and thermocycles (12 h light at 22\u00B0C/12 h dark at 12\u00B0C) for 7 days and then transferred to free-running conditions (continuous light and temperature).
      Data correspond to a microarray analysis, performed from ZT0 (subjective dawn) and during 44 h (with samples every 4 h). For this dataset, the phase information was collected from the tool <a href="http://phaser.mocklerlab.org/">"Phaser"</a>, using a correlation cutoff of 0.7.',
      "<br>",
      '- <b>Michael et al_LL_LLHC</b>: 7955 circadian transcripts identified at the transcriptome level in <a href="https://doi.org/10.1371/journal.pgen.0040014">Michael et al. (2008)</a>.
      Data were obtained from whole seedlings that were grown in continuous light and thermocycles (12 h at 22\u00B0C/12 h at 12\u00B0C) for 7 days and then transferred to free-running conditions (continuous light and temperature).
      Data correspond to a microarray analysis, performed from ZT0 (subjective dawn) and during 44 h (with samples every 4 h). For this dataset, the phase information was collected from the tool <a href="http://phaser.mocklerlab.org/">"Phaser"</a>, using a correlation cutoff of 0.7.',
      "<br>",
      '- <b>Romanowski et al_LL_LDHH</b>: 9128 circadian transcripts identified at the transcriptome level in <a href="https://doi.org/10.1111/tpj.14776">Romanowski et al. (2020)</a>.
      Data were obtained from whole seedlings that were grown in light (12 h) and dark (12 h) cycles at 22\u00B0C for 12 days, and then transferred to free-running conditions (continuous light and temperature).
      Data correspond to mRNA-seq, performed from ZT24 (subjective dawn) on day 14, and during 44 h (with samples every 4 h). The detection of circadian oscillations was performed using the tools <a href="https://doi.org/10.1177/0748730410379711">"JTK_CYCLE"</a>.',
      "<br>",
      "<u>Note</u>: For this reference, phases 24 and 26 from the initial dataset have been changed in 0 and 2, respectively.",
      "</p>","<p style='text-align:justify'>",
      "3. Paste the list of genes you are interested in.",
      "<br>",
      "When clicking on Submit, this will generate three different plots and a table. 
      For the fold enrichment calculation, proportions of phases within the user subset of genes are compared with those of the selected reference.
      Chi-Square tests are then performed and significance is judged at P-value < 0.05.",
      "</p>","<p style='text-align:justify'>"
      
    )
  })
  
  output$beginTab3 <- renderText({
    validate(
      need(!is.null(input$Chose_ref), "Set parameters and click on the Submit button")
    )
  })
  
  Sel_ref<-eventReactive(input$Submit_Tab3, {
    # req(input$Chose_ref)
    validate(
      need(input$Chose_ref != "", "Select a type of phase reference dataset before clicking on the Submit button")
    )
    
    if (input$Chose_ref=="exist")
    {
      validate(
        need(input$Existing_ref != "", "Pick a phase reference dataset before clicking on the Submit button")
      )
      list<-input$Existing_ref
      Sel_ref<-filter(Phase_enrichment, Dataset %in% list)
    }
    
    else {
      inFile <- input$file1
      
      validate(
        need(input$file1 != "", "Upload a phase reference dataset before clicking on the Submit button")
      )
      
      Sel_ref<- read.csv(inFile$datapath, header=T, sep=",")
      
      validate(
        need(ncol(Sel_ref) == 2, "Verify the instructions related to the format of the dataset to upload in the introduction tab. The phase reference dataset provided should be composed of two columns, with comma as column separator.")
      )
      names(Sel_ref) <- c("AGI","Phase")
      
    }
    
    Sel_ref
  })
  
  Sel_genes <- eventReactive(input$Submit_Tab3, {
    Phase_enrichment1<-Sel_ref()
    
    if (input$Chose_ref=="exist")
    {
      list<-unlist(strsplit(toupper(input$Genes_list_Tab3), "\n"))
      validate(
        need(length(list) > 0, "Provide a list of genes before clicking on the Submit button")
      )
      Sel_genes<-filter(Phase_enrichment1, AGI %in% list)
      validate(
        need(nrow(Sel_genes) > 0, "None of the genes provided at Step 3 are circadian, or they are not present in the reference dataset")
      )
    }
    
    else if (input$Chose_ref=="brow")
      {
      list<-unlist(strsplit(input$Genes_list_Tab3, "\n"))
      validate(
        need(length(list) > 0, "Provide a list of genes before clicking on the Submit button")
      )
      Sel_genes<-filter(Phase_enrichment1, AGI %in% list)
      validate(
        need(nrow(Sel_genes) > 0, "None of the genes provided at Step WTF are circadian, or they are not present in the reference dataset")
      )
    }
    Sel_genes
  })
  
  output$Genes_phase_subset <- downloadHandler(
    filename = "Genes_phase_subset.txt.txt",
    content = function(file) {
      write.table(Sel_genes(), file, row.names = FALSE)
    }
  )
  
  output$Genes_phase_ref <- downloadHandler(
    filename = "Genes_phase_ref.txt.txt",
    content = function(file) {
      write.table(Sel_ref(), file, row.names = FALSE)
    }
  )
  
  nb_gene_X<-eventReactive(input$Submit_Tab3,{
    list<-unlist(strsplit(toupper(input$Genes_list_Tab3), "\n"))
    nb<-nrow(as.data.frame(list))
    nb
  })
  
  nb_gene_Y<-eventReactive(input$Submit_Tab3,{
    nb<-nrow(as.data.frame(Sel_genes()))
    nb
  })
  
  
  output$Selection_info <-renderText({
    Upload<-Sel_genes()
    
    # validate(
    #   need(input$Existing_ref != "", ""),
    #   need(nrow(Upload) > 0, "")
    # )
    nb_gene_X<-as.character(nb_gene_X())
    nb_gene_Y<-as.character(nb_gene_Y())
    
    paste(
      "Number of genes in the user subset:", nb_gene_X,"<br>",
      "Number of genes from the user subset that are circadian and therefore included in the analysis:", nb_gene_Y)
  })
  
  #=========== Tab 3 Box 1 - Phase distrib ref ("Seven")=====
  
  Nb_Seven <- reactive ({
    
    Reference<-Sel_ref()
    
    Ref.phases <- data.frame(table(Reference$Phase))
    names(Ref.phases) <- c("Phase","Freq")
    
    withProgress(message = 'Making plot 1', value = 1, {
      
      Plot.ref <- ggplot(Ref.phases, aes(x = as.numeric(as.character(Phase)), y = Freq))+
        geom_bar(stat = "identity")+
        theme_bw()+
        coord_polar( start= -pi/24.5)+
        scale_x_continuous(breaks = seq(0, 23, 1), limits=c(-1,24),expand=c(0,-0.5)) +
        theme(axis.text = element_text(color = "black", size = 13),
              # axis.text.x = element_text(margin=unit(c(0, 0, 0, 0), "cm")),
              axis.title = element_text(size = 14),
              # panel.grid = element_line(color = "grey"),
              panel.grid.major.x = element_line(size = 0.5, linetype = 'solid', colour = "grey70"),
              panel.grid.major.y = element_line(size= 0.5, linetype = 'dotted', colour = "grey70"),
              panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90"), 
              title = element_text(size = 13))+
        xlab("Phase")+
        ylab("Counts")
    })
    
  })
  
  output$Seven <- renderPlot({
    print(Nb_Seven())
  })
  
  output$Check <-renderText({
    paste("Check instructions for format"
    )
  })
  
  output$Legend_phase_ref <-renderText({
    paste("<p style='text-align:justify'>","The phase is defined as the timing of peak abundance (a phase of 0 and 12 indicates a peak abundance at subjective dawn and the beginning of subjective night, respectively).
    The circular barplot represents the counts of the different phases identified in the selected reference.",
          "</p>"
    )
  })
  
  output$Phase_distrib_refPNG <- downloadHandler(
    filename = "Phase_distrib_ref.png.png",
    content = function(file){
      ggsave(file, plot = Nb_Seven(), width = 5, height = 5, device = "png")
    })
  
  output$Phase_distrib_refSVG <- downloadHandler(
    filename = "Phase_distrib_ref.svg.svg",
    content = function(file){
      ggsave(file, plot = Nb_Seven(), width = 5, height = 5, device = "svg")
    })
  
  output$Phase_distrib_refPDF <- downloadHandler(
    filename = "Phase_distrib_ref.pdf.pdf",
    content = function(file){
      pdf(file,width = 5, height = 5)
      print(Nb_Seven())
      dev.off()
    })
  
  #=========== Tab 3 Box 2 - Phase distrib subset ("Eight")=====
  
  Nb_Eight <- eventReactive(input$Submit_Tab3, {
    
    Upload<-Sel_genes()
    # validate(
    #   need(nrow(Upload) > 0, "Provide a list of genes before clicking on the Submit button")
    # )
    
    Upload <- data.frame(table(Upload$Phase))
    names(Upload) <- c("Phase","Freq")
    
    withProgress(message = 'Making plot 2', value = 1, {
      
      Plot.Upload <- ggplot(Upload, aes(x = as.numeric(as.character(Phase)), y = Freq))+
        geom_bar(stat = "identity")+
        theme_bw()+
        coord_polar(start = -pi/24.5)+
        scale_x_continuous(breaks = seq(0, 23, 1), limits=c(-1,24),expand=c(0,-0.5)) +
        theme(axis.text = element_text(color = "black", size = 13),
              axis.title = element_text(size = 14),
              # axis.line.y=element_line(color="white"),
              panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey70"),
              panel.grid.major.y = element_line(size= 0.5, linetype = 'dotted', colour = "grey70"),
              panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90"), 
              title = element_text(size = 13))+
        xlab("Phase")+
        ylab("Counts")
    })
    
  })
  
  output$Eight <- renderPlot({
    print(Nb_Eight())
  })
  
  output$Legend_phase_subset <-renderText({
    paste("<p style='text-align:justify'>","The phase is defined as the timing of peak abundance (a phase of 0 and 12 indicates a peak abundance at subjective dawn and the beginning of subjective night, respectively).
    The circular barplot represents the counts of the different phases identified in the user subset of genes, from the information found in the defined phase reference dataset. 
          Only circadian genes are represented on this plot.",
          "</p>"
    )
  })
  
  output$Phase_distrib_subsetPNG <- downloadHandler(
    filename = "Phase_distrib_subset.png.png",
    content = function(file){
      ggsave(file, plot = Nb_Eight(), width = 5, height = 5, device = "png")
    })
  
  output$Phase_distrib_subsetSVG <- downloadHandler(
    filename = "Phase_distrib_subset.svg.svg",
    content = function(file){
      ggsave(file, plot = Nb_Eight(), width = 5, height = 5, device = "svg")
    })
  
  output$Phase_distrib_subsetPDF <- downloadHandler(
    filename = "Phase_distrib_subset.pdf.pdf",
    content = function(file){
      pdf(file,width = 5, height = 5)
      print(Nb_Eight())
      dev.off()
    })
  
  #=========== Tab 3 Box 3 - Phase enrichment ("Nine")======
  
  Calculated_data <- eventReactive(input$Submit_Tab3, {
    Reference<-Sel_ref()
    Ref.phases <- data.frame(table(Reference$Phase))
    names(Ref.phases) <- c("Phase","Freq")
    
    Upload<-Sel_genes()
    # validate(
    #   need(nrow(Upload) > 0, "Provide a list of genes before clicking on the Submit button")
    # )
    Upload <- data.frame(table(Upload$Phase))
    names(Upload) <- c("Phase","Freq")
    
    # merge results from the Upload with the reference
    names(Upload)[2] <- "Upload"
    Results.phase <- left_join(Ref.phases, Upload, by = "Phase")
    Results.phase<-na.omit(Results.phase)
    
    # Calculate enrichment
    #---------------------
    
    # genes that don't have the phase
    Results.phase$Ref_N <- sum(Results.phase$Freq) - Results.phase$Freq
    Results.phase$Upload_N <- sum(Results.phase$Upload) - Results.phase$Upload
    
    # Calculate proportions of phases
    Results.phase$Ref_proportion <- Results.phase$Freq/sum(Results.phase$Freq)
    Results.phase$Upload_proportion <- Results.phase$Upload/sum(Results.phase$Upload)
    
    # Calculate enrichment
    Results.phase$Enrichment <- Results.phase$Upload_proportion/Results.phase$Ref_proportion
    
    # Calculate significance
    #-----------------------
    list.y.var <- unique(Results.phase$Phase)
    
    Sub = list()
    result = list()
    
    for(i in list.y.var)
    {
      sub = Results.phase[Results.phase$Phase == i, c(3,5,2,4)]
      names(sub) <- rep(c("N","non_N"),2)
      sub = rbind(sub[,1:2],sub[,3:4])
      sub = as.matrix(sub)
      khi2 <- chisq.test(sub)
      result[[i]] <- chisq.test(sub)
      Sub[[i]] <- khi2$p.value
    }
    
    Chi_square <- do.call("rbind", lapply(Sub, as.data.frame))
    Chi_square <- data.frame(row.names(Chi_square),Chi_square)
    names(Chi_square) <- c("Phase", "Chi2")
    
    Results.phase <- merge.data.frame(Results.phase, Chi_square, by = "Phase")
    
    # define the threshold for significance
    Results.phase$Color <- ifelse(Results.phase$Enrichment < 1, "Under-represented","Over-represented")
    Results.phase$Signif <- ifelse(Results.phase$Chi2 < 0.05, "Significant","Non-significant")
    Results.phase$ColSign <- paste(Results.phase$Color, ", ", Results.phase$Signif)
    Results.phase
    
  })
  
  Nb_Nine <- reactive({
    
    withProgress(message = 'Making plot 3', value = 1, {
      
      Legend<-Plot.enrichment <- ggplot(For_PhasEnr_Legend, aes(x = as.numeric(as.character(Phase)), y = Enrichment,  color = ColSign, alpha = ColSign, pch=ColSign))+
        scale_alpha_manual(values = c(1,1,1,1))+
        geom_point(size = 4, stroke=1.5)+
        theme_bw()+
        theme(legend.text = element_text(size = 13, color = "black", margin = margin(t = 4, b=4)),
              legend.title = element_blank(),
              legend.position = "bottom",
              legend.box = "horizontal",
              legend.direction = "vertical",
              legend.margin = margin(0,1,0.5,1, unit="cm"),
              title = element_text(size = 13))+
        scale_color_manual(values = c("purple","purple","limegreen","limegreen"))+
        scale_shape_manual(values = c(1,19,1,19))
      
      lg<- cowplot::get_legend(Legend)
      #lg<-as_ggplot(Legend)
      
      
      Plot.enrichment <- ggplot(Calculated_data(), aes(x = as.numeric(as.character(Phase)), y = Enrichment,  color = Color, alpha = Signif, pch=Signif))+
        geom_hline(aes(yintercept=1), lwd=1, lty=1) + 
        scale_alpha_manual(values = c(1,1,1,1))+
        coord_polar(theta = "x", start=0)+
        geom_point(size = 4, stroke=1.5)+
        scale_x_continuous(breaks = seq(0, 23, 1), limits=c(0,24)) +
        theme_bw()+
        theme(axis.text = element_text(color = "black", size = 13),
              panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "grey70"),
              panel.grid.major.y = element_line(size= 0.5, linetype = 'dotted', colour = "grey70"),
              panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90"), 
              axis.title = element_text(size = 14),
              # legend.text = element_text(size = 13, color = "black", margin = margin(t = 4, b=4)),
              # legend.title = element_blank(),
              # legend.position = "bottom",
              # legend.box = "horizontal",
              # legend.direction = "vertical",
              # legend.margin = margin(0,1,0.5,1, unit="cm"),
              legend.position = "none",
              title = element_text(size = 13))+
        scale_color_manual(values = c("purple","limegreen"))+
        scale_shape_manual(values = c(1,19))+
        xlab("Phase")+
        ylab("Fold enrichment")
      
      # plot_grid(Plot.enrichment,lg,nrow=2, rel_heights=c(11,2.5), align="vh")
      
      plot_grid(Plot.enrichment,lg,nrow=2, rel_heights=c(11,2.5), align="vh")
      
      #rel_heights=c(1.2,0.2,3.2), axis = "r"
      
    })
    
  })
  
  output$Nine <- renderPlot({
    print(Nb_Nine())
  })
  
  output$Legend_phase_enrichment <-renderText({
    paste("<p style='text-align:justify'>","The phase is defined as the timing of peak abundance (a phase of 0 and 12 indicates a peak abundance at subjective dawn and the beginning of subjective night, respectively). 
    Proportions of the different phases are compared between the user subset of genes and the defined phase reference data, and the circular bubble plot represents the over- and under-represented phases in the user subset of genes as compared to the defined phase reference dataset. 
          Only genes identified as circadian are considered for this analysis. 
          Chi-Square tests are performed and significance is judged at P-value < 0.05.",
          "</p>"
    )
  })
  
  output$Phase_enrichmentPNG <- downloadHandler(
    filename = "Phase_enrichment.png.png",
    content = function(file){
      ggsave(file, plot = Nb_Nine(), width = 5, height = 6, device = "png")
    })
  
  output$Phase_enrichmentSVG <- downloadHandler(
    filename = "Phase_enrichment.svg.svg",
    content = function(file){
      ggsave(file, plot = Nb_Nine(), width = 5, height = 6, device = "svg")
    })
  
  output$Phase_enrichmentPDF <- downloadHandler(
    filename = "Phase_enrichment.pdf.pdf",
    content = function(file){
      pdf(file,width = 5, height = 6)
      print(Nb_Nine())
      dev.off()
    })
  
  #=========== Tab 3 Box 4 - Summary table ("Ten")=====
  
  Nb_Ten <- reactive({
    
    withProgress(message = 'Computing', value = 1, {
      
      Results.phase<-Calculated_data()
      Khi2.table <- Results.phase[,c(1:3,8,9)]
      names(Khi2.table) <- c("Phase","Ref.List","User.List","Fold.enrichment","P_value")
      
      Khi2.table <- Khi2.table[order(Khi2.table$Phase),]
      
      Khi2.table.plot <- gt(Khi2.table) %>% 
        tab_options(table.border.top.color = "Grey") %>%
        tab_style(
          style = list(
            cell_fill(color = "lightgreen"),
            cell_text(weight = "bold")),
          locations = cells_body(
            columns = vars("P_value"),
            rows = P_value< 0.05))%>%
        fmt_number(columns=vars("Fold.enrichment","P_value"),
                   rows= "P_value" > 0.001, decimals=3)%>%
        cols_align(align = "center")%>%
        cols_label(
          Phase = html("Phase<br> "),
          Ref.List = html("Reference list"),
          User.List = html("User list<br> "),
          Fold.enrichment= html("Fold enrichment"),
          P_value= html("P value<br> ")
        )
      
      Khi2.table.plot 
    })
  })
  
  
  output$Ten = render_gt({
    Nb_Ten()
  })
  
  # mydf<-reactive({
  #   Gene.list.rlog.TOT<-rlog.TOT() 
  #   
  #   if (input$Sort_by == "order_by_transcriptome_ph") {
  #     Gene.list.rlog.TOT <- Gene.list.rlog.TOT[order(Gene.list.rlog.TOT$Phase),]
  #   }
  #   else if (input$Sort_by == "order_by_translatome_ph") {
  #     Gene.list.rlog.TOT <- Gene.list.rlog.TOT[order(Gene.list.rlog.TOT$Phase_TRAP),]
  #   }
  #   else {NULL}
  # })
  
  # output$mytable = renderDataTable({
  # Sel_genes()
  # })
  
  for_txt_table <- reactive({
    Results.phase<-Calculated_data()
    Khi2.table <- Results.phase[,c(1:3,8,9)]
    names(Khi2.table) <- c("Phase","Ref_List","User_List","Fold_enrichment","P_value")
    
    Khi2.table <- Khi2.table[order(Khi2.table$Phase),]
  })
  
  output$Phase_enrich_summary <- downloadHandler(
    filename = "Summary_phase_enrichment.txt.txt",
    content = function(file) {
      write.table(for_txt_table(), file, row.names = FALSE)
    }
  )
  
  
  output$Legend_summary_table <-renderText({
    paste("<p style='text-align:justify'>","Phase: timing of peak abundance (a phase of 0 and 12 indicates a peak abundance at subjective dawn and the beginning of subjective night, respectively).",
          "<br>",
          "Reference list: counts of the different phases identified in the defined phase reference dataset.",
          "<br>",
          "User list: counts of the different phases identified in the user subset of genes, from the information found in the defined phase reference dataset.",
          "<br>",
          "Fold enrichment: Ratio between the proportion of the phase in the user subset of genes vs in the defined phase reference dataset.",
          "<br>",
          "P value: result of a Chi-square test. P values indicated in bold and highlighted in green are < 0.05 and are considered significantly over- (if fold enrichment > 1) or under-represented (if fold enrichment < 1) in the user subset of genes.",
          "</p>"
    )
  })
  
  
  
  ####========= Server Tab 4 "Network" ========
  #=======================================#
  output$Tab4_description <-renderText({# Network
    paste("<p style='text-align:justify'>",
          'This tab allows you to visualize if your gene(s) or gene family are targets of clock proteins in Arabidopsis. 
          Clock proteins accumulate at different times of day and regulate the expression of clock genes, as shown in the simplified model below:'   
    )
  })
  
  Cl_genes_id<-clock[,c(1,19)]
  
  Cl_genes_table <- gt(Cl_genes_id) %>% 
    cols_align(align = "center")%>%
    tab_style(
      style = cell_text(size = px(13)),
      locations = list(
        cells_column_labels(columns = everything()),
        cells_body(columns = everything()))
    )%>%
    # tab_style(
    #   style = cell_text(size = px(13)),
    #   locations = cells_body(
    #     columns = vars(id, AGI)))%>%
    cols_label(
      id = html("Protein name"))
  
  
  output$Clock_gene_AGI <- render_gt({
    Cl_genes_table
  })
  
  output$Clock_genes_image <- renderUI({img(src="Clock_genes.png", height = "150px")})
  
  output$tab4_description2 <-renderText({# Network
    paste("<p style='text-align:justify'>",
          'Connections between clock proteins and their targets were identified in published ChIP-Seq data, for CCA1 (<a href="https://doi.org/10.1073/pnas.1513609112">Nagel <i>et al.,</i> 2015</a>; 
                                                                                                                     <a href="https://doi.org/10.1105/tpc.15.00737">Kamioka <i>et al.,</i> 2016</a>), 
          LHY (<a href="https://doi.org/10.1111/nph.15415">Adams <i>et al.,</i> 2018</a>), 
          PRR5 (<a href="https://doi.org/10.1073/pnas.1205156109">Nakamichi <i>et al.,</i> 2012</a>), 
          PRR7 (<a href="https://doi.org/10.1111/tpj.12276">Liu <i>et al.,</i> 2013</a>), 
          PRR9 (<a href="https://doi.org/10.1104/pp.15.01562">Liu <i>et al.,</i> 2016</a>), 
          TOC1 (<a href="https://science.sciencemag.org/content/336/6077/75/">Huang <i>et al.,</i> 2012</a>), 
          and LUX (<a href="https://www.nature.com/articles/nplants201787">Ezer <i>et al.,</i> 2017</a>).',
          "</p>","<p style='text-align:justify'>",
          'To build a network, you need to:',
          "<br>",
          '1. Select the clock proteins for which you want to see the targets.',
          "</p>","<p style='text-align:justify'>",
          '2. Select the input genes you want to visualize in the network, by either selecting:',
          "<br>", '- "transcription factor family", then pick a TF family. The TF families were defined in (<a href="https://doi.org/10.1016/j.celrep.2014.06.033">Pruneda-Paz <i>et al.,</i> 2014</a>).',
          "<br>", '- "other list of genes", which allows you to paste a list of Arabidopsis AGI numbers in the indicated area.',
          "<br>",
          'Click on Submit to generate the network. This network is generated using the R package <a href="https://cran.r-project.org/web/packages/visNetwork/index.html">"visNetwork"</a>.',
          "</p>","<p style='text-align:justify'>",
          '3. By default, nodes are colored based on their phase (timing of peak abundance) in the transcriptome dataset published in <a href="https://doi.org/10.1093/plcell/koab113">Bonnot and Nagel (2021)</a>. 
          The phase was determined from RNA-Seq data obtained in a 24 h time course (with samples every 3 h) under control conditions at 22\u00B0C. The R package <a href="https://doi.org/10.1093/bioinformatics/btw405">"Metacycle"</a> was used to detect circadian oscillations, and calculated a "Phase" and an "Adjusted p-value" for each transcript. 
          In <a href="https://doi.org/10.1093/plcell/koab113">Bonnot and Nagel (2021)</a>, oscillations were considered significant for genes: with an Adjusted p-value < 0.01; with an 0.01 < Adjusted p-value < 0.05 and that overlapped with previously published lists of circadian transcripts.
          For details, please see the Methods section in <a href="https://doi.org/10.1093/plcell/koab113">Bonnot and Nagel (2021)</a>.',
          "<br>",
          'Nodes can also be colored based on their phase in the translatome data or based on their heat stress response in the transcriptome and translatome datasets (<a href="https://doi.org/10.1093/plcell/koab113">Bonnot and Nagel, 2021</a>). 
          The heat stress response corresponds to data obtained in response to heat stress treatments that were applied at different times of day (from ZT48, subjective dawn, to ZT69, end of subjective night) on separate sets of 12 days old Arabidopsis seedlings (whole seedlings), and that correspond to treatments at 37\u00B0C for 1 h. 
          For this network analysis, Log2 Fold Change (log2 FC) values (37\u00B0C vs 22\u00B0C) calculated at the different times of day were averaged to determine the average transcript response to heat stress.',
          "<br>",
          'In addition, the last option in the "criteria for node color" allows you to not color the nodes based on a specific criteria.',
          "</p>","<p style='text-align:justify'>",
          '<u>Note 1</u>: The network is interactive. Nodes can be moved, and new nodes and connections can be manually added. We are not responsible for the accuracy of the network if new nodes and connections have been added by the user.',
          "<p style='text-align:justify'>",
          "<u>Note 2</u>: Currently, the R package <a href='https://cran.r-project.org/web/packages/visNetwork/index.html'>'visNetwork'</a> does not allow exporting the network image as a .pdf. 
          If this option is integrated in future versions of the package, we will implement this functionality and update the Shiny app. 
          Right now, the best option seems to be taking a screenshot that saves as a .jpeg or .tiff, then exporting the image as a .pdf.",
          "</p>","<p style='text-align:justify'>"
          
    )
  })
  
  output$beginTab4 <- renderText({
    validate(
      need(!is.null(input$Chose_list), "Set parameters and click on the Submit button")
    )
  })
  
 
  New_clock<-eventReactive(input$Submit_Tab4, {
    validate(
      need(input$Clock_to_show, "Select at least one clock protein")
    )
    
    clock<-clock[,1:17]

    list<-input$Clock_to_show
    new_clock<-filter(clock, id %in% list)
    
  })
  
  New_edges<-eventReactive(input$Submit_Tab4, {
    
    list_clock_genes<-input$Clock_to_show
    
    req(input$Chose_list)
    if (input$Chose_list == "TFfamily") {
      list<-input$TFfam
      ID<-filter(TFfamilylist, Family %in% list)
      Edges<-filter(Edges_original, to %in% ID$AGI)
      Edges<-filter(Edges, from %in% list_clock_genes)

    } else if (input$Chose_list == "paste_list_network") {
      
      list<-unlist(strsplit(toupper(input$Genes_list_Tab4), "\n"))
      Edges<-filter(Edges_original, to %in% list)
      Edges<-filter(Edges, from %in% list_clock_genes)
    }
  })
  
  
  New_nodes<-eventReactive(input$Submit_Tab4, {
    
    req(input$Chose_list)
    clock<-New_clock()
    Nodes_original<-Nodes_original[,1:17]
    edges<-New_edges()
    
    if (input$Chose_list == "TFfamily") {
      list<-input$TFfam
      ID<-filter(TFfamilylist, Family %in% list)
      Nodes<-filter(Nodes_original, id %in% ID$AGI)
      validate(
        need(nrow(Nodes)>0, "No identified connections between the selected genes and the selected clock proteins in the ChIP-Seq datasets used to build the network")
      )
      Nodes<-filter(Nodes, id %in% edges$to)
      validate(
        need(nrow(Nodes)>0, "No identified connections between the selected genes and the selected clock proteins in the ChIP-Seq datasets used to build the network")
      )
      Nodes <- Nodes[order(Nodes$id),]
      Nodes <- rbind(clock, Nodes)
    } else if (input$Chose_list == "paste_list_network") {
      
      list<-unlist(strsplit(toupper(input$Genes_list_Tab4), "\n"))
      Nodes<-filter(Nodes_original, id %in% list)
      Nodes<-filter(Nodes, id %in% edges$to)
      Nodes <- Nodes[order(Nodes$id),]
      Nodes <- rbind(clock, Nodes)
    }
  })
  
    

  
  Col_Nodes<-reactive({
    #req(input$Chose_list)
    Nodes<-New_nodes()
    Nodes1 <- data.frame(Nodes)
    
    if (input$Color_choice == "color_Transcriptome_phases") {
      
      No <- Nodes[,c(1:3,4:5,17)]
      names(No)[4] <- "group"
      names(No)[5] <- "color"
      # No[,4] <- factor(No[,4], levels = c("0-1.5","3-4.5","6-7.5","9-10.5","12-13.5","15-16.5","18-19.5","21-22.5","Not circadian"))
      # return(No)
    } else if (input$Color_choice == "color_Translatome_phases") {
      
      No <- Nodes[,c(1:3,6:7,17)]
      names(No)[4] <- "group"
      names(No)[5] <- "color"
    } else if (input$Color_choice == "color_Transcriptome_heatstress") {
      
      No <- Nodes[,c(1:3,8:9,17)]
      names(No)[4] <- "group"
      names(No)[5] <- "color"
    } else if (input$Color_choice == "color_Translatome_heatstress") {
      
      No <- Nodes[,c(1:3,10:11,17)]
      names(No)[4] <- "group"
      names(No)[5] <- "color"
    } else if (input$Color_choice == "No_col") {
      
      No <- Nodes[,c(1:3,12:13,17)]
      names(No)[4] <- "group"
      names(No)[5] <- "color"
    }
    No
  })
  
  #=========== Tab 4 Box 1 - visNetwork ("Eleven")=====
  
  #https://datastorm-open.github.io/visNetwork/
  #https://visjs.github.io/vis-network/docs/network/
  
 
  output$Eleven <- renderVisNetwork({
    
    withProgress(message = 'Computing network', value = 1, {
    
    Nodes<-Col_Nodes()
    Edges<-New_edges()
      
      visNetwork(Nodes, Edges, width = "300px", height = 700) %>%
        visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T),
                   selectedBy = "group",
                   nodesIdSelection = TRUE,
                   manipulation = TRUE)%>%
        visEdges(arrows = "to")%>%
        # visPhysics(enabled = FALSE)%>%
        visEvents(stabilizationIterationsDone="function () {this.setOptions( { physics: false } );}",
                  # click="function(ng_nodes){Shiny.onInputChange('got_network_current_node_id',ng_nodes);}",
                  selectNode = "function(nodes){Shiny.onInputChange('click', nodes.nodes[0]);}")
      #%>%
        #visExport(type = "jpeg", name = "export-network", float = "left", label = "Save network", background = "white", style= "") 
    })
    
  })
  
  # observeEvent(input$store_position, {
  #   visNetworkProxy("Eleven") %>% visGetPositions()
  # })
  
  output$Legend_network <-renderText({
    paste("<p style='text-align:justify'>","Nodes correspond to clock proteins (squares) and their known targets from the list of input genes that has been selected or provided (circles). 
          Edges correspond to connections from the clock proteins to their targets, based on published ChIP-Seq data (for details, see 'Instructions and Methodological details).
          If 'Phases at the transcriptome level' or 'Phases at the translatome level' have been selected for the criteria for node selection, the 'Phase' defines the timing of peak abundance (a phase of 0 and 12 indicates a peak abundance at subjective dawn and the beginning of subjective night, respectively).
          Genes identified as 'Not circadian' did not exhibit significant circadian oscillations in <a href='https://doi.org/10.1093/plcell/koab113'>Bonnot and Nagel (2021)</a>.
          If 'Heat stress response at the transcriptome level' or 'Heat stress response at the translatome level' have been selected for the criteria for node selection, log2 FC corresponds to the average of different Log2 Fold Change (log2 FC) values (37\u00B0C vs 22\u00B0C) calculated at eight different times of day (from ZT48, beginning of day, to ZT69, end of night).
          Blue and red indicate an average downregulation and upregulation in response to heat stress, respectively.",
          "</p>"
    )
  })
  
  output$node_info <- renderText({
    paste("ID of the node selected by click: ", input$click)
  })
  
  output$click_ui <- renderUI({
    if (is.null(input$got_network_current_node_id) )
    {
      ""
    }
    else if (length(input$got_network_current_node_id$node) == 0)
    {
      ""
    }
    else
    {
      nodeid <- input$got_network_current_node_id$nodes
      tempnodeid <-unlist(nodeid)
      nodedata = subset(ng_nodes, (ng_nodes$id %in% tempnodeid))
      x<-as.character(nodedata$label)
      
      print(paste0("Selected Nodes:",paste(unlist(x), collapse = ", ")))
    }
    
  })
  
  
  Dat_nodes <- reactive({
    Noeuds<-Col_Nodes()
    Noeuds <- Noeuds[,c(1,4)]
    
    if (input$Color_choice == "color_Transcriptome_phases") {
      
      names(Noeuds) <- c("id", "Phases_transcriptome")
    } else if (input$Color_choice == "color_Translatome_phases") {
      
      names(Noeuds) <- c("id","Phases_translatome")
    } else if (input$Color_choice == "color_Transcriptome_heatstress") {
      
      names(Noeuds) <- c("id", "Heat_stress_transcriptome")
    } else if (input$Color_choice == "color_Translatome_heatstress") {
      
      names(Noeuds) <- c("id","Heat_stress_translatome")
    }
    Noeuds
  })
  
  Dat_edges <- reactive({
    Ed<-New_edges()
    Ed <- Ed[,c(1,2)]
  })
  
  output$Nodes_file <- downloadHandler(
    filename = "Data_nodes.txt.txt",
    content = function(file) {
      write.table(Dat_nodes(), file, row.names = FALSE)
    })
  
  output$Edges_file <- downloadHandler(
    filename = "Data_edges.txt.txt",
    content = function(file) {
      write.table(Dat_edges(), file, row.names = FALSE)
    })
  
  
  nb_gene_X1<-eventReactive(input$Submit_Tab4,{
    req(input$Chose_list)
    if (input$Chose_list == "TFfamily") {
      list<-input$TFfam
      ID<-filter(TFfamilylist, Family %in% list)
      nb<-nrow(as.data.frame(ID))
      nb
    } else if (input$Chose_list == "paste_list_network") {
      list<-unlist(strsplit(toupper(input$Genes_list_Tab4), "\n"))
      list<-unique(list)
      nb<-nrow(as.data.frame(list))
      nb
    }
  })
  
  nb_gene_Y1<-eventReactive(input$Submit_Tab4,{
    data<-New_edges()
    counts<-data %>% count(to)
    nb<-nrow(counts)
  })
  
  nb_gene_Y2<-eventReactive(input$Submit_Tab4,{
    data<-New_edges()
    nb<-nrow(data)
  })
  
  nb_gene_CCA1<-eventReactive(input$Submit_Tab4,{
    data<-New_edges()
    if ( "CCA1" %in% input$Clock_to_show) {
      counts<- length(which(data$from=="CCA1"))
    }
    else {data1=NA}
  })
  
  # output$mytable = renderDataTable({
  #   New_edges()
  # })

  nb_gene_LHY<-eventReactive(input$Submit_Tab4,{
    data<-New_edges()
    if ( "LHY" %in% input$Clock_to_show) {
      counts<- length(which(data$from=="LHY"))
    }
    else {data1=NA}
  })
  
  nb_gene_LUX<-eventReactive(input$Submit_Tab4,{
    data<-New_edges()
    if ( "LUX" %in% input$Clock_to_show) {
      counts<- length(which(data$from=="LUX"))
    }
    else {data1=NA}
  })
  
  nb_gene_PRR5<-eventReactive(input$Submit_Tab4,{
    data<-New_edges()
    if ( "PRR5" %in% input$Clock_to_show) {
      counts<- length(which(data$from=="PRR5"))
    }
    else {data1=NA}
  })
  
  nb_gene_PRR7<-eventReactive(input$Submit_Tab4,{
    data<-New_edges()
    if ( "PRR7" %in% input$Clock_to_show) {
      counts<- length(which(data$from=="PRR7"))
    }
    else {data1=NA}
  })
  
  nb_gene_PRR9<-eventReactive(input$Submit_Tab4,{
    data<-New_edges()
    if ( "PRR9" %in% input$Clock_to_show) {
      counts<- length(which(data$from=="PRR9"))
    }
    else {data1=NA}
  })
  
  nb_gene_TOC1<-eventReactive(input$Submit_Tab4,{
    data<-New_edges()
    if ( "TOC1" %in% input$Clock_to_show) {
      counts<- length(which(data$from=="TOC1"))
    }
    else {data1=NA}
  })
  
  output$Selection_info_Network <-renderText({

    nb_gene_X1<-as.character(nb_gene_X1())
    nb_gene_Y1<-as.character(nb_gene_Y1())
    nb_gene_Y2<-as.character(nb_gene_Y2())
    
    nb_gene_CCA1<-as.character(nb_gene_CCA1())
    nb_gene_LHY<-as.character(nb_gene_LHY())
    nb_gene_LUX<-as.character(nb_gene_LUX())
    nb_gene_PRR5<-as.character(nb_gene_PRR5())
    nb_gene_PRR7<-as.character(nb_gene_PRR7())
    nb_gene_PRR9<-as.character(nb_gene_PRR9())
    nb_gene_TOC1<-as.character(nb_gene_TOC1())
    
    paste(
      "Number of input genes: ", nb_gene_X1,"<br>",
      "Number of edges: ", nb_gene_Y2, "<br>",
      "Number of input genes targeted by selected clock proteins: ", nb_gene_Y1,"<br>","<br>",
      "Number of input genes targeted by:","<br>",
      "CCA1: ", nb_gene_CCA1,"<br>",
      "LHY: ", nb_gene_LHY,"<br>",
      "LUX: ", nb_gene_LUX,"<br>",
      "PRR5: ", nb_gene_PRR5,"<br>",
      "PRR7: ", nb_gene_PRR7,"<br>",
      "PRR9: ", nb_gene_PRR9,"<br>",
      "TOC1: ", nb_gene_TOC1)
  })
  
  #=========== Tab 4 Box 2 - text ("Tab4_Legend")=====
  
  output$Tab4_Legend <-renderPlot({
    DF<-New_nodes()
    if (input$Color_choice == "color_Transcriptome_phases") {
      DF$Phases_transcriptome <- factor(DF$Phases_transcriptome, levels = c("0-1.5","3-4.5","6-7.5","9-10.5","12-13.5","15-16.5","18-19.5","21-22.5","Not circadian"))
      names(DF)[names(DF) == "Phases_transcriptome"] <- "Phase"
      g<-ggplot(DF, aes(x=id,y=label,shape=Type,color=Phase))+
        geom_point(size=5,stroke=1.25) + 
        theme_classic(base_size=15)+
        theme(legend.key = element_rect(size = 5, color=NA),
              legend.key.size = unit(1.5, 'lines'))+
        scale_color_manual(values=c("0-1.5"="#FFCC00","3-4.5"="#FF9900","6-7.5"="#FF6600","9-10.5"="#FF0000",
                                    "12-13.5"="#CC0099","15-16.5"="#9900FF","18-19.5"="#0000FF","21-22.5"="#000066","Not circadian"="#CCCCCC"))+
        scale_shape_manual(values=c(0,1))+
        guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))
      
      legend <- cowplot::get_legend(g)
      leg<-as.ggplot(legend)
    }
    
    else if (input$Color_choice == "color_Translatome_phases") {
      DF$Phases_translatome <- factor(DF$Phases_translatome, levels = c("0-1.5","3-4.5","6-7.5","9-10.5","12-13.5","15-16.5","18-19.5","21-22.5","Not circadian"))
      names(DF)[names(DF) == "Phases_translatome"] <- "Phase"
      g<-ggplot(DF, aes(x=id,y=label,shape=Type,color=Phase))+
        geom_point(size=5,stroke=1.25) + 
        theme_classic(base_size=15)+
        theme(legend.key = element_rect(size = 5, color=NA),
              legend.key.size = unit(1.5, 'lines'))+
        scale_color_manual(values=c("0-1.5"="#FFCC00","3-4.5"="#FF9900","6-7.5"="#FF6600","9-10.5"="#FF0000",
                                    "12-13.5"="#CC0099","15-16.5"="#9900FF","18-19.5"="#0000FF","21-22.5"="#000066","Not circadian"="#CCCCCC"))+
        scale_shape_manual(values=c(0,1))+
        guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))
      
      legend <- cowplot::get_legend(g)
      leg<-as.ggplot(legend)
    }
    
    else if (input$Color_choice == "color_Transcriptome_heatstress") {
      DF$Heat_stress_transcriptome <- factor(DF$Heat_stress_transcriptome, levels = c("> 3", "[2 - 3]", "[1 - 2]", "[0 - 1]", "[(-1) - 0]", "[(-2) - (-1)]", "[(-3) - (-2)]", "< (-3)"))
      g<-ggplot(DF, aes(x=id,y=label,shape=Type,color=Heat_stress_transcriptome))+
        geom_point(size=5,stroke=1.25) + 
        theme_classic(base_size=15)+
        theme(legend.key = element_rect(size = 5, color=NA),
              legend.key.size = unit(1.5, 'lines'))+
        scale_color_manual(values=c("> 3"="#F40000","[2 - 3]"="#FF6A6A","[1 - 2]"="#FFB2B2","[0 - 1]"="#FFEBEB",
                                    "[(-1) - 0]"="#E5EAFF","[(-2) - (-1)]"="#C0CCFF","[(-3) - (-2)]"="#92A7FF","< (-3)"="#204AFF"))+
        scale_shape_manual(values=c(0,1))+
        labs(shape="Type", color="log2 FC \n(37\u00B0C vs 22\u00B0C)")+
        guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))
      
      legend <- cowplot::get_legend(g)
      leg<-as.ggplot(legend)
    }
    else if (input$Color_choice == "color_Translatome_heatstress") {
      DF$Heat_stress_translatome <- factor(DF$Heat_stress_translatome, levels = c("> 3", "[2 - 3]", "[1 - 2]", "[0 - 1]", "[(-1) - 0]", "[(-2) - (-1)]", "[(-3) - (-2)]", "< (-3)"))
      g<-ggplot(DF, aes(x=id,y=label,shape=Type,color=Heat_stress_translatome))+
        geom_point(size=5,stroke=1.25) + 
        theme_classic(base_size=15)+
        theme(legend.key = element_rect(size = 5, color=NA),
              legend.key.size = unit(1.5, 'lines'))+
        scale_color_manual(values=c("> 3"="#F40000","[2 - 3]"="#FF6A6A","[1 - 2]"="#FFB2B2","[0 - 1]"="#FFEBEB",
                                    "[(-1) - 0]"="#E5EAFF","[(-2) - (-1)]"="#C0CCFF","[(-3) - (-2)]"="#92A7FF","< (-3)"="#204AFF"))+
        scale_shape_manual(values=c(0,1))+
        labs(shape="Type", color="log2 FC \n(37\u00B0C vs 22\u00B0C)")+
        guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))
      
      legend <- cowplot::get_legend(g)
      leg<-as.ggplot(legend)
    }
    else if (input$Color_choice == "No_col") {
      g<-ggplot(DF, aes(x=id,y=label,shape=Type, color=Type))+
        geom_point(size=5,stroke=1.25) + 
        theme_classic(base_size=15)+
        theme(legend.key = element_rect(size = 5, color=NA),
              legend.key.size = unit(1.5, 'lines'))+
        scale_shape_manual(values=c(15,16))+
        scale_color_manual(values=c("#7F7F7F","#D0CECE"))
      
      legend <- cowplot::get_legend(g)
      leg<-as.ggplot(legend)
    }
    
    leg
  })
  

  
  
  ####========= Server Tab 5 "Multispecies circadian oscillations" ========
  #=======================================#
  output$Tab5_description <-renderText({# Multispecies circadian oscillations
    paste(
      "<p style='text-align:justify'>",
      'This tab allows you to represent gene expression patterns during the day, in different published datasets and in different species: 
      <i>Arabidopsis thaliana</i>, <i>Brassica rapa</i>, <i>Hordeum vulgare</i>, and <i>Oryza sativa</i>.',
      "</p>","<p style='text-align:justify'>",
      "1. Pick an Arabidopsis time course dataset: these datasets correspond to all expressed genes in <i>Arabidopsis thaliana</i> published datasets from:",
      "<br>",
      '- <b>Bonnot and Nagel_Transcriptome_LL_LDHH</b>: Twelve days old seedlings (whole seedlings) were grown in light (12 h) and dark (12 h) cycles at 22\u00B0C for 10 days and then transferred to free-running conditions (continuous light and temperature) for two days before sampling (<a href="https://doi.org/10.1093/plcell/koab113">Bonnot and Nagel, 2021</a>).
      Data correspond to mRNA-Seq, performed on a 24 h time course (with samples every 3 h) under control conditions at 22\u00B0C.',
      "<br>",
      '- <b>Bonnot and Nagel_Translatome_LL_LDHH</b>: Twelve days old seedlings (whole seedlings) were grown in light (12 h) and dark (12 h) cycles at 22\u00B0C for 10 days and then transferred to free-running conditions (continuous light and temperature) for two days before sampling (<a href="https://doi.org/10.1093/plcell/koab113">Bonnot and Nagel, 2021</a>).
      Data correspond to TRAP-Seq, performed on a 24 h time course (with samples every 3 h) under control conditions at 22\u00B0C.',
      "<br>",
      '- <b>Covington and Harmer_LL_LDHH</b>: Data were obtained from whole seedlings that were grown in light (12 h) and dark (12 h) cycles at 22\u00B0C for 7 days and then transferred to free-running conditions (continuous light and temperature, <a href="https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050222"> Covington and Harmer, 2007</a>).
      Data correspond to a microarray analysis, performed from ZT24 (subjective dawn) on day 9 and during 44 h (with samples every 4 h).',
      "<br>",
      '- <b>Edwards et al_LL_LDHH</b>: Data were obtained from whole seedlings that were grown in light (12 h) and dark cycles (12 h) at 22\u00B0C for 7 days and then transferred to free-running conditions (continuous light and temperature, <a href="https://doi.org/10.1105/tpc.105.038315">Edwards et al., 2006</a>).
      Data correspond to a microarray analysis, performed from ZT26 (2 h after subjective dawn) on day 9 and during two days (with samples every 4 h).',
      "<br>",
      '- <b>Michael et al_LL_LDHC</b>: Data were obtained from whole seedlings that were grown in light/dark cycles and thermocycles (12 h light at 22\u00B0C/12 h dark at 12\u00B0C) for 7 days and then transferred to free-running conditions (continuous light and temperature, <a href="https://doi.org/10.1371/journal.pgen.0040014">Michael et al., 2008</a>).
      Data correspond to a microarray analysis, performed from ZT0 (subjective dawn) and during 44 h (with samples every 4 h).',
      "<br>",
      '- <b>Michael et al_LL_LLHC</b>: Data were obtained from whole seedlings that were grown in continuous light and thermocycles (12 h at 22\u00B0C/12 h at 12\u00B0C) for 7 days and then transferred to free-running conditions (continuous light and temperature, <a href="https://doi.org/10.1371/journal.pgen.0040014">Michael et al., 2008</a>).
      Data correspond to a microarray analysis, performed from ZT0 (subjective dawn) and during 44 h (with samples every 4 h).',
      "<br>",
      '- <b>Romanowski et al_LL_LDHH</b>: Data were obtained from whole seedlings that were grown in light (12 h) and dark (12 h) cycles at 22\u00B0C for 12 days, and then transferred to free-running conditions (continuous light and temperature, <a href="https://doi.org/10.1111/tpj.14776">Romanowski et al., 2020</a>).
      Data correspond to mRNA-seq, performed from ZT24 (subjective dawn) on day 14, and during 44 h (with samples every 4 h).',
      "<br>",
      "Enter an AGI number in the designated area and click on submit to visualize the transcript profile in the selected dataset, in the window 'Time course in <i>Arabidopsis thaliana</i>'.",
      "<p style='text-align:justify'>",
      "2. Pick a Brassica time course dataset: these datasets correspond to all expressed genes in <i>Brassica rapa</i>, from <a href='https://elifesciences.org/articles/58993'>Greenham <em>et al.</em> (2020)</a>.",
      "<br>",
      "Plants were grown in light (12 h) and dark (12 h) cycles and constant 20\u00B0C (<b>Greenham et al_LL_LDHH</b>), or in continuous light and thermocycle (12 h at 20\u00B0C/12 h at 10\u00B0C) conditions (<b>Greenham et al_LL_LLHC</b>) for 15 days after sowing, before transfer to free-running conditions (continuous light and temperature at 20\u00B0C).
      Following 24 h in constant conditions, leaf tissues were collected from ZT24 (subjective dawn), then every 2 h for 48 h, and mRNA-Seq was performed.",
      "<br>",
      "Enter a locus number in the designated area and click on submit to visualize the transcript profile in the selected dataset, in the window 'Time course in <i>Brassica rapa</i>'.",
      "</p>","<p style='text-align:justify'>",
      "3. Pick a Barley time course dataset: the provided dataset (<b>Muller et al_LL_LDHH</b>) corresponds to all expressed genes in <i>Hordeum vulgare</i>, from <a href='https://academic.oup.com/plphys/article/183/2/765/6116414'>Müller <em>et al.</em> (2020)</a>.",
      "<br>",
      "Plants were grown for three weeks in light (12 h) and dark (12 h) cycles, then transferred to free-running conditions (continuous light and temperature at 20\u00B0C).
      Leaf tissues were collected every 4 h for 36 h starting from the first subjective night, and mRNA-Seq was performed.",
      "<br>",
      "Enter a locus number in the designated area and click on submit to visualize the transcript profile in the window 'Time course in <i>Hordeum vulgare</i>'.",
      "</p>","<p style='text-align:justify'>",
      "4. Pick a Rice time course dataset: these datasets correspond to all expressed genes in <i>Oryza sativa</i>, from <a href='https://doi.org/10.1371/journal.pone.0016907'>Filichkin <em>et al.</em> (2011)</a>.",
      "<br>",
      "Plants were grown for three months: in light (12 h) and dark (12 h) cycles at a constant daytime temperature (31\u00B0C, <b>Filichkin et al_LLHH_LDHH</b>); 
      in light (12 h) and dark (12 h) cycles and thermocycles (12 h at 31\u00B0C/12 h at 20\u00B0C, <b>Filichkin et al_LLHH_LDHC</b>); 
      or in thermocycles (12 h at 31\u00B0C/12 h at 20\u00B0C) and continuous light (<b>Filichkin et al_LLHH_LLHC</b>). 
      Plants were then transferred to free-running conditions (continuous light and temperature). 
      Following 48 h in constant conditions, leaf tissues were collected from ZT0 (subjective dawn), then every 4 h for 48 h, and a microarray analysis was performed.",
      "<br>",
      "Enter a locus number in the designated area and click on submit to visualize the transcript profile in the window 'Time course in <i>Oryza sativa</i>'.",
      "<p style='text-align:justify'>"
    )
  })
  
  output$beginTab5 <- renderText({
    validate(
      need(input$AGI_oscillation, "Set parameters and click on the Submit button")
    )
  })
  
  output$beginTab5.1 <- renderText({
    validate(
      need(input$AGI_oscillation2, "Set parameters and click on the Submit button")
    )
  })
  
  output$beginTab5.2 <- renderText({
    validate(
      need(input$AGI_oscillation3, "Set parameters and click on the Submit button")
    )
  })
  
  output$beginTab5.3 <- renderText({
    validate(
      need(input$AGI_oscillation4, "Set parameters and click on the Submit button")
    )
  })
  
  Arabi_ref<-eventReactive(input$Submit_Tab5, {
    # req(input$Existing_ref1)
    validate(
      need(input$Existing_ref1 != "", "Pick an Arabidopsis time course dataset before clicking on the Submit button")
    )
    list<-input$Existing_ref1
    Arabi_ref<-filter(Arabidopsis_data, Condition %in% list)
  })
  
  Arabi_gene<-eventReactive(input$Submit_Tab5, {
    # req(input$AGI_oscillation)
    Arabi_ref<-Arabi_ref()
    validate(
      need(input$AGI_oscillation != "", "Provide an Arabidopsis gene ID before clicking on the Submit button")
    )
    list<-input$AGI_oscillation
    Arabi_gene<-filter(Arabi_ref, AGI %in% list)
  })
  
  Rice_ref<-eventReactive(input$Submit_Tab5.1, {
    # req(input$Existing_ref2)
    validate(
      need(input$Existing_ref2 != "", "Pick a Rice time course dataset before clicking on the Submit button")
    )
    list<-input$Existing_ref2
    Rice_ref<-filter(Rice_data, Condition %in% list)
  })
  
  Rice_gene<-eventReactive(input$Submit_Tab5.1, {
    # req(input$AGI_oscillation2)
    Rice_ref<-Rice_ref()
    validate(
      need(input$AGI_oscillation2 != "", "Provide a Rice gene ID before clicking on the Submit button")
    )
    list<-input$AGI_oscillation2
    Rice_gene<-filter(Rice_ref, AGI %in% list)
  })
  
  Bra_ref<-eventReactive(input$Submit_Tab5.2, {
    # req(input$Existing_ref2)
    validate(
      need(input$Existing_ref3 != "", "Pick a Brassica time course dataset before clicking on the Submit button")
    )
    list<-input$Existing_ref3
    Bra_ref<-filter(Brassica_data, Condition %in% list)
  })
  
  Bra_gene<-eventReactive(input$Submit_Tab5.2, {
    # req(input$AGI_oscillation2)
    Bra_ref<-Bra_ref()
    validate(
      need(input$AGI_oscillation3 != "", "Provide a Brassica rapa gene ID before clicking on the Submit button")
    )
    list<-input$AGI_oscillation3
    Bra_gene<-filter(Bra_ref, AGI %in% list)
  })  
  
  Barley_ref<-eventReactive(input$Submit_Tab5.3, {
    # req(input$Existing_ref2)
    validate(
      need(input$Existing_ref4 != "", "Pick a Barley time course dataset before clicking on the Submit button")
    )
    list<-input$Existing_ref4
    Barley_ref<-filter(Barley_data, Condition %in% list)
  })
  
  Barley_gene<-eventReactive(input$Submit_Tab5.3, {
    # req(input$AGI_oscillation2)
    Barley_ref<-Barley_ref()
    validate(
      need(input$AGI_oscillation4 != "", "Provide a Barley gene ID before clicking on the Submit button")
    )
    list<-input$AGI_oscillation4
    Barley_gene<-filter(Barley_ref, AGI %in% list)
  })  
  
  #=========== Tab 5 Box 1 - Arabidopsis oscillation ("Twelve")=====
  
  Nb_Twelve <- eventReactive(input$Submit_Tab5, {
    Gene1 <- Arabi_gene()
    validate(
      need(nrow(Gene1)>0, "This gene is miswritten (check instructions) or not included in the dataset")
    )
    Gene1 <- merge.data.frame(Grey_areas, Gene1, by = "Condition")
    #xmin <- Gene1[1,5]
    #xmax <- Gene1[1,6]
    #Start2 <- Gene1[1,7]
    #End2 <- Gene1[1,8]
    
    withProgress(message = "Making plot", value = 1, {
    
    Plot.1 <- ggplot(Gene1, aes(x = Time, y = value)) +
      geom_rect(aes(xmin=Start1, xmax=End1, ymin=-Inf, ymax=Inf), 
                fill="grey", alpha=0.06, color = NA)+
      geom_rect(aes(xmin=Start2, xmax=End2, ymin=-Inf, ymax=Inf), 
                fill="grey", alpha=0.06, color = NA)+
      #annotate("rect", xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, alpha = 0.5)+
      #annotate("rect", xmin = Start2, xmax = End2, ymin = -Inf, ymax = Inf, alpha = 0.5)+
      geom_point(size = 2)+
      geom_line()+
      theme_bw()+
      xlab("Time of day (ZT)")+
      ylab("Transcript abundance")+
      theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
            axis.text=element_text(size=10, colour = "black"),axis.text.x=element_text(margin = unit(c(3, 0, 3, 0), "mm")),
            axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
            axis.title = element_text(size = 12))+
      scale_x_continuous(breaks = Gene1$Time)
    })
    
  })
  
  output$Twelve <- renderPlot({
    print(Nb_Twelve())
  })
  
  output$Legend_Arabido <-renderText({
    paste("<p style='text-align:justify'>","The grey area represents the subjective night. The black line represents the mean of different biological replicates",
          "</p>"
    )
  })
  
  output$Arabidopsis_PNG <- downloadHandler(
    filename = "Arabidopsis_oscillation.png",
    content = function(file){
      ggsave(file, plot = Nb_Twelve(), width = 5, height = 5, device = "png")
    })
  
  output$Arabidopsis_SVG <- downloadHandler(
    filename = "Arabidopsis_oscillation.svg",
    content = function(file){
      ggsave(file, plot = Nb_Twelve(), width = 5, height = 5, device = "svg")
    })
  
  output$Arabidopsis_PDF <- downloadHandler(
    filename = "Arabidopsis_oscillation.pdf",
    content = function(file){
      pdf(file,width = 5, height = 5)
      print(Nb_Twelve())
      dev.off()
    })
  
  output$Arabidopsis_timecourse <- downloadHandler(
    filename = "Arabidopsis_timecourse.txt",
    content = function(file) {
      #write.table(datasetInputOne(), file, row.names = FALSE)
      write.table(Arabi_gene(), file, row.names = FALSE)
    }
  )
  
  
  #=========== Tab 5 Box 2 - Rice oscillation ("Thirteen")=====
  
  Nb_Thirteen <- eventReactive(input$Submit_Tab5.1, {
    Gene2 <- Rice_gene()
    validate(
      need(nrow(Gene2)>0, "This gene is miswritten (check instructions) or not included in the dataset")
    )
    
    Gene2 <- merge.data.frame(Grey_areas, Gene2, by = "Condition")
    
    withProgress(message = "Making plot", value = 1, {
    
    Plot.2 <- ggplot(Gene2, aes(x = Time, y = value, group = Probe)) +
      geom_rect(aes(xmin=Start1, xmax=End1, ymin=-Inf, ymax=Inf), 
                fill="grey", alpha=0.02, color = NA)+
      geom_rect(aes(xmin=Start2, xmax=End2, ymin=-Inf, ymax=Inf), 
                fill="grey", alpha=0.02, color = NA)+
      geom_point(size = 2)+
      geom_line()+
      theme_bw()+
      xlab("Time of day (ZT)")+
      ylab("Transcript abundance")+
      theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
            axis.text=element_text(size=10, colour = "black"),axis.text.x=element_text(margin = unit(c(3, 0, 3, 0), "mm")),
            axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
            axis.title = element_text(size = 12))+
      scale_x_continuous(breaks = Gene2$Time)
    })
    
  })
  
  output$Thirteen <- renderPlot({
    print(Nb_Thirteen())
  })
  
  output$Legend_Rice <-renderText({
    paste("<p style='text-align:justify'>","The grey area represents the subjective night. The black line represents the mean of different biological replicates. If multiple profiles are represented on the plot, they correspond to different probes in the microarray analysis.",
          "</p>"
    )
  })
  
  output$Rice_PNG <- downloadHandler(
    filename = "Rice_oscillation.png",
    content = function(file){
      ggsave(file, plot = Nb_Thirteen(), width = 5, height = 5, device = "png")
    })
  
  output$Rice_SVG <- downloadHandler(
    filename = "Rice_oscillation.svg",
    content = function(file){
      ggsave(file, plot = Nb_Thirteen(), width = 5, height = 5, device = "svg")
    })
  
  output$Rice_PDF <- downloadHandler(
    filename = "Rice_oscillation.pdf",
    content = function(file){
      pdf(file,width = 5, height = 5)
      print(Nb_Thirteen())
      dev.off()
    })
  
  output$Rice_timecourse <- downloadHandler(
    filename = "Rice_timecourse.txt",
    content = function(file) {
      #write.table(datasetInputOne(), file, row.names = FALSE)
      write.table(Rice_gene(), file, row.names = FALSE)
    }
  )
  
  
  #=========== Tab 5 Box 3 - Brassica oscillation ("Fourteen")=====
  
  Nb_Fourteen <- eventReactive(input$Submit_Tab5.2, {
    Gene3 <- Bra_gene()
    validate(
      need(nrow(Gene3)>0, "This gene is miswritten (check instructions) or not included in the dataset")
    )
    
    Gene3 <- merge.data.frame(Grey_areas, Gene3, by = "Condition")
    
    withProgress(message = "Making plot", value = 1, {
    
    Plot.3 <- ggplot(Gene3, aes(x = Time, y = value)) +
      geom_rect(aes(xmin=Start1, xmax=End1, ymin=-Inf, ymax=Inf), 
                fill="grey", alpha=0.02, color = NA)+
      geom_rect(aes(xmin=Start2, xmax=End2, ymin=-Inf, ymax=Inf), 
                fill="grey", alpha=0.02, color = NA)+
      geom_point(size = 2)+
      geom_line()+
      theme_bw()+
      xlab("Time of day (ZT)")+
      ylab("Transcript abundance")+
      theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
            axis.text=element_text(size=10, colour = "black"),axis.text.x=element_text(margin = unit(c(3, 0, 3, 0), "mm")),
            axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
            axis.title = element_text(size = 12))+
      scale_x_continuous(breaks = Gene3$Time)
    })
  })
  
  output$Fourteen <- renderPlot({
    print(Nb_Fourteen())
  })
  
  output$Legend_Brassica <-renderText({
    paste("<p style='text-align:justify'>","The grey area represents the subjective night. The black line represents the mean of different biological replicates.",
          "</p>"
    )
  })
  
  output$Brassica_PNG <- downloadHandler(
    filename = "Brassica_oscillation.png",
    content = function(file){
      ggsave(file, plot = Nb_Fourteen(), width = 5, height = 5, device = "png")
    })
  
  output$Brassica_SVG <- downloadHandler(
    filename = "Brassica_oscillation.svg",
    content = function(file){
      ggsave(file, plot = Nb_Fourteen(), width = 5, height = 5, device = "svg")
    })
  
  output$Brassica_PDF <- downloadHandler(
    filename = "Brassica_oscillation.pdf",
    content = function(file){
      pdf(file,width = 5, height = 5)
      print(Nb_Fourteen())
      dev.off()
    })
  
  output$Brassica_timecourse <- downloadHandler(
    filename = "Brassica_timecourse.txt",
    content = function(file) {
      #write.table(datasetInputOne(), file, row.names = FALSE)
      write.table(Bra_gene(), file, row.names = FALSE)
    }
  )
  
  
  #=========== Tab 5 Box 4 - Barley oscillation ("Fifteen")=====
  
  Nb_Fifteen <- eventReactive(input$Submit_Tab5.3, {
    Gene4 <- Barley_gene()
    validate(
      need(nrow(Gene4)>0, "This gene is miswritten (check instructions) or not included in the dataset")
    )
    
    Gene4 <- merge.data.frame(Grey_areas, Gene4, by = "Condition")
    
    withProgress(message = "Making plot", value = 1, {
    
    Plot.4 <- ggplot(Gene4, aes(x = Time, y = value)) +
      geom_rect(aes(xmin=Start1, xmax=End1, ymin=-Inf, ymax=Inf), 
                fill="grey", alpha=0.06, color = NA)+
      geom_rect(aes(xmin=Start2, xmax=End2, ymin=-Inf, ymax=Inf), 
                fill="grey", alpha=0.06, color = NA)+
      geom_point(size = 2)+
      geom_line()+
      theme_bw()+
      xlab("Time of day (ZT)")+
      ylab("Transcript abundance")+
      theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
            axis.text=element_text(size=10, colour = "black"),axis.text.x=element_text(margin = unit(c(3, 0, 3, 0), "mm")),
            axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
            axis.title = element_text(size = 12))+
      scale_x_continuous(breaks = Gene4$Time)
    })
  })
  
  output$Fifteen <- renderPlot({
    print(Nb_Fifteen())
  })
  
  output$Legend_Barley <-renderText({
    paste("<p style='text-align:justify'>","The grey area represents the subjective night. The black line represents the mean of different biological replicates.",
          "</p>"
    )
  })
  
  
  output$Barley_PNG <- downloadHandler(
    filename = "Barley_oscillation.png",
    content = function(file){
      ggsave(file, plot = Nb_Fifteen(), width = 5, height = 5, device = "png")
    })
  
  output$Barley_SVG <- downloadHandler(
    filename = "Barley_oscillation.svg",
    content = function(file){
      ggsave(file, plot = Nb_Fifteen(), width = 5, height = 5, device = "svg")
    })
  
  output$Barley_PDF <- downloadHandler(
    filename = "Barley_oscillation.pdf",
    content = function(file){
      pdf(file,width = 5, height = 5)
      print(Nb_Fifteen())
      dev.off()
    })
  
  output$Barley_timecourse <- downloadHandler(
    filename = "Barley_timecourse.txt",
    content = function(file) {
      #write.table(datasetInputOne(), file, row.names = FALSE)
      write.table(Barley_gene(), file, row.names = FALSE)
    }
  )
  
  
  ####========= Server Tab 6 "How to cite" ========
  #=======================================#
  output$cite_refs <-renderText({# Welcome
    paste("<p>","<b>R Core Team</b> (2020). R: A Language and Environment for Statistical Computing.",
          "<br>",
          "<b>Chang W, Cheng J, Allaire J, Xie Y, McPherson J </b> (2020). shiny: Web Application Framework for R.",
          "<br>",
          "<b>Wickham H</b> (2016) ggplot2: Elegant Graphics for Data Analysis.",
          
          "</p>","<p style='text-align:justify'>",
          "<p>", "<u>Tab 'Single genes'</u>",
          "</p>","<p style='text-align:justify'>",
          "<p>", 
          "<b>Bonnot T, Nagel DH</b> (2021). Time of day prioritizes the pool of translating mRNAs in response to heat stress. <i>Plant Cell.</i> doi: <a href='https://doi.org/10.1093/plcell/koab113'>10.1093/plcell/koab113 </a>",
          "<br>",
          "<b>Love MI, Huber W, Anders S</b> (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. <i> Genome Biol.</i> doi: <a href='https://doi.org/10.1186/s13059-014-0550-8'>10.1186/s13059-014-0550-8</a>",
          "<br>",
          "<b>Wu G, Anafi RC, Hughes ME, Kornacker K, Hogenesch JB</b> (2016). MetaCycle: an integrated R package to evaluate periodicity in large scale data. <i>Bioinformatics.</i> doi: <a href='https://doi.org/10.1093/bioinformatics/btw405'>10.1093/bioinformatics/btw405</a>",
          "</p>","<p style='text-align:justify'>",
          
          "<p>", "<u>Tab 'Multiple genes'</u>",
          "</p>","<p style='text-align:justify'>",
          "<p>", 
          "<b>Bonnot T, Nagel DH</b> (2021). Time of day prioritizes the pool of translating mRNAs in response to heat stress. <i>Plant Cell.</i> doi: <a href='https://doi.org/10.1093/plcell/koab113'>10.1093/plcell/koab113 </a>",
          "<br>",
          "<b>Love MI, Huber W, Anders S</b> (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. <i> Genome Biol.</i> doi: <a href='https://doi.org/10.1186/s13059-014-0550-8'>10.1186/s13059-014-0550-8</a>",
          "<br>",
          "<b>Pruneda-Paz JL, Breton G, Nagel DH, Kang SE, Bonaldi K, Doherty CJ, Ravelo S, Galli M, Ecker JR, Kay SA</b> (2014). A Genome-Scale Resource for the Functional Characterization of Arabidopsis Transcription Factors. <i>Cell Rep</i> 8: 622-632, doi: <a href='https://doi.org/10.1016/j.celrep.2014.06.033'>10.1016/j.celrep.2014.06.033</a>",
          "<br>",
          "<b>Wu G, Anafi RC, Hughes ME, Kornacker K, Hogenesch JB</b> (2016). MetaCycle: an integrated R package to evaluate periodicity in large scale data. <i>Bioinformatics.</i> doi: <a href='https://doi.org/10.1093/bioinformatics/btw405'>10.1093/bioinformatics/btw405</a>",
          "</p>", "<p style='text-align:justify'>",
          
          "<p>", "<u>Tab 'Network'</u>",
          "</p>","<p style='text-align:justify'>",
          "<p>", "Data:","<br>",
          "<b>Adams S, Grundy J, Veflingstad SR, Dyer NP, Hannah MA, Ott S, Carré IA</b> (2018). Circadian control of abscisic acid biosynthesis and signalling pathways revealed by genome-wide analysis of LHY binding targets. <i>New Phytol</i> 220: 893-907. doi: <a href='https://doi.org/10.1111/nph.15415'>10.1111/nph.15415</a>",
          "<br>",
          "<b>Bonnot T, Nagel DH</b> (2021). Time of day prioritizes the pool of translating mRNAs in response to heat stress. <i>Plant Cell.</i> doi: <a href='https://doi.org/10.1093/plcell/koab113'>10.1093/plcell/koab113 </a>",
          "<br>",
          "<b>Ezer D, Jung J-H, Lan H, Biswas S, Gregoire L, Box MS, Charoensawan V, Cortijo S, Lai X, Stöckle D, et al</b> (2017). The evening complex coordinates environmental and endogenous signals in Arabidopsis. <i>Nat Plants</i> 3: 17087. doi: <a href='https://www.nature.com/articles/nplants201787'>10.1038/nplants.2017.87</a>",
          "<br>",
          "<b>Huang W, Pérez-García P, Pokhilko A, Millar AJ, Antoshechkin I, Riechmann JL, Mas P</b> (2012). Mapping the Core of the Arabidopsis Circadian Clock Defines the Network Structure of the Oscillator. <i>Science</i> (80- ) 336: 75-79. doi: <a href='https://science.sciencemag.org/content/336/6077/75/'>10.1126/science.1219075</a>",
          "<br>",
          "<b>Kamioka M, Takao S, Suzuki T, Taki K, Higashiyama T, Kinoshita T, Nakamichi N</b> (2016). Direct Repression of Evening Genes by CIRCADIAN CLOCK-ASSOCIATED1 in the Arabidopsis Circadian Clock. <i>Plant Cell</i> 28: 696-711. doi: <a href='https://doi.org/10.1105/tpc.15.00737'>10.1105/tpc.15.00737</a>",
          "<br>",
          "<b>Liu T, Carlsson J, Takeuchi T, Newton L, Farré EM</b> (2013). Direct regulation of abiotic responses by the Arabidopsis circadian clock component PRR7. <i>Plant</i> J 76: 101-14. doi: <a href='https://doi.org/10.1111/tpj.12276'>10.1111/tpj.12276</a>",
          "<br>",
          "<b>Liu TL, Newton L, Liu M-J, Shiu S-H, Farré EM</b> (2016). A G-Box-Like Motif Is Necessary for Transcriptional Regulation by Circadian Pseudo-Response Regulators in Arabidopsis. <i>Plant Physiol</i> 170: 528-539. doi: <a href='https://doi.org/10.1104/pp.15.01562'>10.1104/pp.15.01562</a>",
          "<br>",
          "<b>Nagel DH, Doherty CJ, Pruneda-Paz JL, Schmitz RJ, Ecker JR, Kay SA</b> (2015). Genome-wide identification of CCA1 targets uncovers an expanded clock network in Arabidopsis. <i>Proc Natl Acad Sci</i> 112: E4802-E4810. doi: <a href='https://doi.org/10.1073/pnas.1513609112'>10.1073/pnas.1513609112</a>",
          "<br>",
          "<b>Nakamichi N, Kiba T, Kamioka M, Suzuki T, Yamashino T, Higashiyama T, Sakakibara H, Mizuno T</b> (2012). Transcriptional repressor PRR5 directly regulates clock-output pathways. <i>Proc Natl Acad Sci</i> 109: 17123-17128. doi: <a href='https://doi.org/10.1073/pnas.1205156109'>10.1073/pnas.1205156109</a>",
          "<br>",
          "<b>Pruneda-Paz JL, Breton G, Nagel DH, Kang SE, Bonaldi K, Doherty CJ, Ravelo S, Galli M, Ecker JR, Kay SA</b> (2014). A Genome-Scale Resource for the Functional Characterization of Arabidopsis Transcription Factors. <i>Cell Rep</i> 8: 622-632, doi: <a href='https://doi.org/10.1016/j.celrep.2014.06.033'>10.1016/j.celrep.2014.06.033</a>",
          "<br>",
          "<b>Wu G, Anafi RC, Hughes ME, Kornacker K, Hogenesch JB</b> (2016). MetaCycle: an integrated R package to evaluate periodicity in large scale data. <i>Bioinformatics.</i> doi: <a href='https://doi.org/10.1093/bioinformatics/btw405'>10.1093/bioinformatics/btw405</a>",
          "</p>", "<p style='text-align:justify'>",
          
          "</p>","<p style='text-align:justify'>",
          "<p>", "Plots:",
          "<br>",
          "<b>Almende BV, Thieurmel B, Robert T</b> (2019). visNetwork: Network Visualization using 'vis.js' Library. ",
          "</p>","<p style='text-align:justify'>",
          
          "<p>", "<u>Tab 'Phase enrichment'</u>",
          "</p>","<p style='text-align:justify'>",
          "<p>",
          "<b>Bonnot T, Nagel DH</b> (2021) Time of day prioritizes the pool of translating mRNAs in response to heat stress. <i>Plant Cell.</i> doi: <a href='https://doi.org/10.1093/plcell/koab113'>10.1093/plcell/koab113 </a>",
          "<br>",
          "<b>Covington MF, Harmer SL</b> (2007). The Circadian Clock Regulates Auxin Signaling and Responses in Arabidopsis. <i> PLoS Biology.</i> doi: <a href='https://doi.org/10.1371/journal.pbio.0050222'>10.1371/journal.pbio.0050222</a>",
          "<br>",
          "<b>Edwards KD, Anderson PE, Hall A, Salathia NS, Locke JCW, Lynn JR, Straume M, Smith JQ, Millar AJ</b> (2006). FLOWERING LOCUS C Mediates Natural Variation in the High-Temperature Response of the Arabidopsis Circadian Clock. <i> Plant Cell.</i> doi: <a href='https://doi.org/10.1105/tpc.105.038315'>10.1105/tpc.105.038315</a>",
          "<br>",
          "<b>Hsu PY, Harmer SL</b> (2012) Circadian Phase Has Profound Effects on Differential Expression Analysis. <i>PLoS One</i> 7: e49853. doi: <a href=' https://doi.org/10.1371/journal.pone.0049853'>10.1371/journal.pone.0049853</a>",
          "<br>",
          "<b>Hughes ME, Hogenesch JB, Kornacker K</b> (2010). JTK_CYCLE: an efficient nonparametric algorithm for detecting rhythmic components in genome-scale data sets. <i> J Biol Rhythms.</i> doi: <a href='https://doi.org/10.1177/0748730410379711'>10.1177/0748730410379711</a>",
          "<br>",
          "<b>Michael TP, Mockler TC, Breton G, McEntee C, Byer A, Trout JD, Hazen SP, Shen R, Priest HD, Sullivan CM, Givan SA, Yanovsky M, Hong F, Kay SA, Chory J</b> (2008). Network Discovery Pipeline Elucidates Conserved Time-of-Day-Specific cis-Regulatory Modules. <i> PLoS Genet.</i> doi: <a href='https://doi.org/10.1371/journal.pgen.0040014'>10.1371/journal.pgen.0040014</a>",
          "<br>",
          "<b>Mockler TC, Michael TP, Priest HD, Shen R, Sullivan CM, Givan SA, McEntee C, Kay SA, Chory J</b> (2007) The Diurnal Project: Diurnal and Circadian Expression Profiling, Model-based Pattern Matching, and Promoter Analysis. Cold Spring Harb Symp Quant Biol 72: 353-363. doi: <a href=' http://symposium.cshlp.org/content/72/353'>10.1101/sqb.2007.72.006</a>",
          "<br>",
          "<b>Romanowski A, Schlaen RG, Perez-Santangelo S, Mancini E, Yanovsky MJ</b> (2020) Global transcriptome analysis reveals circadian control of splicing events in Arabidopsis thaliana. Plant J 103: 889-902. doi: <a href='https://doi.org/10.1111/tpj.14776'>10.1111/tpj.14776</a>",
          "<br>",
          "<b>Straume M</b> (2004). DNA microarray time series analysis: automated statistical assessment of circadian rhythms in gene expression patterning. <i>Methods Enzymol.</i> doi: <a href='https://doi.org/10.1016/S0076-6879(04)83007-6'>10.1016/S0076-6879(04)83007-6</a>",
          "<br>",
          "<b>Wu G, Anafi RC, Hughes ME, Kornacker K, Hogenesch JB</b> (2016). MetaCycle: an integrated R package to evaluate periodicity in large scale data. <i>Bioinformatics.</i> doi: <a href='https://doi.org/10.1093/bioinformatics/btw405'>10.1093/bioinformatics/btw405</a>",
          "</p>","<p style='text-align:justify'>",
          
          "<p>", "<u>Tab 'Multispecies circadian oscillations'</u>",
          "</p>","<p style='text-align:justify'>",
          "<p>", 
          "<b>Bonnot T, Nagel DH</b> (2021) Time of day prioritizes the pool of translating mRNAs in response to heat stress. <i>Plant Cell.</i> doi: <a href='https://doi.org/10.1093/plcell/koab113'>10.1093/plcell/koab113 </a>",
          "<br>",
          "<b>Covington MF, Harmer SL</b> (2007). The Circadian Clock Regulates Auxin Signaling and Responses in Arabidopsis. <i> PLoS Biology.</i> doi: <a href='https://doi.org/10.1371/journal.pbio.0050222'>10.1371/journal.pbio.0050222</a>",
          "<br>",
          "<b>Edwards KD, Anderson PE, Hall A, Salathia NS, Locke JCW, Lynn JR, Straume M, Smith JQ, Millar AJ</b> (2006). FLOWERING LOCUS C Mediates Natural Variation in the High-Temperature Response of the Arabidopsis Circadian Clock. <i> Plant Cell.</i> doi: <a href='https://doi.org/10.1105/tpc.105.038315'>10.1105/tpc.105.038315</a>",
          "<br>",
          "<b>Filichkin SA, Breton G, Priest HD, Dharmawardhana P, Jaiswal P, Fox SE, Michael TP, Chory J, Kay SA, Mockler TC</b> (2011). Global Profiling of Rice and Poplar Transcriptomes Highlights Key Conserved Circadian-Controlled Pathways and cis-Regulatory Modules. <i>PLoS One.</i> doi: <a href='https://doi.org/10.1371/journal.pone.0016907'>10.1371/journal.pone.0016907</a>",
          "<br>",
          "<b>Greenham K, Sartor RC, Zorich S, Lou P, Mockler TC, McClung CR</b> (2020) Expansion of the circadian transcriptome in <i>Brassica rapa</i> and genome-wide diversification of paralog expression patterns. eLife 2020;9:e58993. doi: <a href='https://elifesciences.org/articles/58993'>10.7554/eLife.58993</a>",
          "<br>",
          "<b>Michael TP, Mockler TC, Breton G, McEntee C, Byer A, Trout JD, Hazen SP, Shen R, Priest HD, Sullivan CM, Givan SA, Yanovsky M, Hong F, Kay SA, Chory J</b> (2008). Network Discovery Pipeline Elucidates Conserved Time-of-Day-Specific cis-Regulatory Modules. <i> PLoS Genet.</i> doi: <a href='https://doi.org/10.1371/journal.pgen.0040014'>10.1371/journal.pgen.0040014</a>",
          "<br>",
          "<b>Mockler TC, Michael TP, Priest HD, Shen R, Sullivan CM, Givan SA, McEntee C, Kay SA, Chory J</b> (2007) The Diurnal Project: Diurnal and Circadian Expression Profiling, Model-based Pattern Matching, and Promoter Analysis. Cold Spring Harb Symp Quant Biol 72: 353-363. doi: <a href=' http://symposium.cshlp.org/content/72/353'>10.1101/sqb.2007.72.006</a>",
          "<br>",
          "<b>Müller LM, Mombaerts L, Pankin A, Davis SJ, Webb AAR, Goncalves J, von Korff M</b> (2020) Differential Effects of Day/Night Cues and the Circadian Clock on the Barley Transcriptome. Plant Physiology 183: 765-779. doi: <a href='https://academic.oup.com/plphys/article/183/2/765/6116414'>10.1104/pp.19.01411</a>",
          "<br>",
          "<b>Romanowski A, Schlaen RG, Perez-Santangelo S, Mancini E, Yanovsky MJ</b> (2020) Global transcriptome analysis reveals circadian control of splicing events in Arabidopsis thaliana. Plant J 103: 889-902. doi: <a href='https://doi.org/10.1111/tpj.14776'>10.1111/tpj.14776</a>",
          "</p>"

    )
  })
  
  output$cite_CASTR <-renderText({# Welcome
    paste("<b>Bonnot T, Gillard MB, Nagel DH</b> (20XX) CAST-R: A shiny application to visualize Circadian And heat STress-Responsive genes in plants. doi:xxxxx"
    )
  })
  
  
}#end of server

#})#for profvis; remove

# Run the app ----
shinyApp(ui = ui, server = server)

