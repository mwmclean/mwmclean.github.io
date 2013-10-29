library(shiny)

# Define UI for dataset viewer application
shinyUI(pageWithSidebar(
  #includeCSS('McLean.css'),
  # Application title.
  headerPanel("Comparison of FLM and FGAM"),
  # tags$head(tags$link(rel="stylesheet", type="text/css", href="mysidepanel.css")),
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view. The helpText function is also used to 
  # include clarifying text. Most notably, the inclusion of a 
  # submitButton defers the rendering of output until the user 
  # explicitly clicks the button (rather than doing it immediately
  # when inputs change). This is useful if the computations required
  # to render output are inordinately time-consuming.
  sidebarPanel(
    tags$head(
      tags$style(type="text/css", "select { max-width: 200px; }"),
      tags$style(type="text/css", "textarea { max-width: 185px; }"),
      tags$style(type="text/css", ".jslider { max-width: 200px; }"),
      tags$style(type='text/css', ".well { max-width: 310px; }"),
      tags$style(type='text/css', ".span4 { max-width: 310px; }")
    ),
    helpText('For introduction to FLM and FGAM, see the Research tab'),
    helpText(strong('Report questions/bugs to:')),
    includeHTML('myemail.html'),
    h4(''),
    selectInput("dataset", "Choose a dataset:", 
                choices = c("Emissions", "DTI", "Tecator", "Temperature", "Precipitation")),
    
#    numericInput("obs", "Number of obs. to highlight:", 8, min = 3, max = 12, step = 1),
     div(class='row',
       div(class="shiny-bound-input offset1", numericInput("obs", "Number of obs. to highlight:", 8, min = 3, max = 12, step = 1))
     ),
     tags$head(tags$style(type="text/css", "#obs { max-width: 185px; }")),
  #  tags$style(type='text/css', "#obs label { margin-left: 20px;}"),
#     helpText("Note: while the data view will show only the specified",
#              "number of observations, the fits will still be based",
#              "on the full dataset."),
    helpText(strong('FLM Paramaters')),
    sliderInput("FLMnbf", "Number of basis functions to use for functional coefficient:", value = 15, min = 5, max = 25, step = 1),
    selectInput("sptypeFLM", "Method for Estimating Smoothing Parameters:", 
                choices = c("REML", "GCV", "Fixed")),
    conditionalPanel(condition = "input.sptypeFLM == 'Fixed'",
      sliderInput("FLMsp", "base-10 log smoothing parameter for coefficient function", value = 1, min = -5, max = 10, step = 1)),
    helpText(strong('FGAM Paramaters')),
    uiOutput('xbfs'),
    uiOutput('tbfs'),
#     sliderInput("FGAMxnbf", "Number of basis functions to use for x-axis basis for FGAM:", value = 10, min = 5, max = 10, step = 1),
#     sliderInput("FGAMtnbf", "Number of basis functions to use for t-axis basis for FGAM:", value = 10, min = 5, max = 10, step = 1),
    selectInput("sptypeFGAM", "Method for Estimating Smoothing Parameters:", 
                choices = c("REML", "GCV", "Fixed")),
    conditionalPanel(condition = "input.sptypeFGAM == 'Fixed'",
    sliderInput("FGAMxsp", "base-10 log smoothing parameter x-axis", value = 1, min = -5, max = 10, step = 1),
    sliderInput("FGAMtsp", "base-10 log smoothing parameter t-axis", value = 1, min = -5, max = 10, step = 1)),
    radioButtons("plottype", "Type of plot for FGAM:",
                 c('persp', 'contour'), selected='persp'),
#     radioButtons("pred", "Perform out of sample prediction?",
#                  c('Yes', 'No'), selected='No'),
    checkboxInput("pred", "Perform out of sample prediction?", FALSE),
    
    uiOutput('n.test'),
    checkboxInput("htest", "Conduct hypothesis test?", FALSE),
    conditionalPanel(condition = "input.htest==1",
                     checkboxGroupInput('tests', 'Tests to use:', c('RLRT1', 'RLRT2', 'Bootstrap'),
                                        selected = 'RLRT1')),
#    conditionalPanel(condition='input.pred',
#                     numericInput("test.prop", "Proportion of samples to use for test set:", .2, min = 0, max = .5, step = .05),#),
    #conditionalPanel(condition='input.pred',
#    helpText("Note: Selecting Zero for the Proportion of test samples results in one sample being used for testing.",
#             "on the full dataset.")),
    actionButton('fitaction', "Fit Models")
  ),
  
  # Show a summary of the dataset and an HTML table with the requested
  # number of observations. Note the use of the h4 function to provide
  # an additional header above each output section.
  mainPanel(
 #   includeCSS('McLean.css'),
#    includeScript('toggle.js'),
    #includeHTML('index.html'),
    h4("Observed Functional Covariates"),
    verbatimTextOutput("summary"),
    plotOutput("matPlot"),
    h4("FLM fit"),
    plotOutput('FLMplot'),
    h4("FGAM fit"),
    plotOutput('FGAMplot'),
#     conditionalPanel(condition = "input.plottype=='persp'",
#                      numericInput('theta','theta:', 0),
#                      numericInput('phi','phi:', 30),
#                      helpPopup('About', 'These define the viewing angles for the perspective plot (azimuth and coaltitude, respectively)',
#                                            placement=c('right', 'top', 'left', 'bottom'),
#                                            trigger='hover')), 
    conditionalPanel(condition = "input.plottype=='persp'&&input.fitaction",
                     div(class='row',
                     div(class="span2 offset1", numericInput('theta','theta:', value=10, step=5)),
                     div(class="span2", numericInput('phi','phi:', value=30, step=5),
                         helpPopup('About', 'These define the viewing angles for the perspective plot (azimuth and coaltitude, respectively)',
                               placement='right', trigger='hover'))
                         ),
                     tags$style(type="text/css", '#theta {width: 50px;}'),
                     tags$style(type="text/css", '#phi {width: 50px;}')
                    ), 
    conditionalPanel(condition = "input.pred==1",
    h4('Out-of-Sample Prediction Root Mean Square Error')),
    conditionalPanel(condition = "input.pred==1",
    verbatimTextOutput("prederror")),
    conditionalPanel(condition = "input.htest==1",
    h4('Hypothesis Test of H_0: FLM vs H_1: FGAM')),
    conditionalPanel(condition = "input.htest==1",
    verbatimTextOutput("testres"))
  )
))