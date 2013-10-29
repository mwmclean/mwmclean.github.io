library(shiny)
library(fda)
library(RColorBrewer)
library(mgcv)
library(refund)
library(lattice)
library(RLRsim)
library(lme4)
library(fda.usc)
source('C:\\Users\\mmclean\\Documents\\testshiny\\fgamm4v2.R')
helpPopup <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
  tagList(
    singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })")
      )
    ),
    tags$a(
      href = "#", class = "btn btn-mini", `data-toggle` = "popover",
      title = title, `data-content` = content, `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok=TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],
      
      tags$i(class="icon-question-sign")
    )
  )
}

data(DTI)
data(tecator)
DTI.ind <- match(unique(DTI$ID), DTI$ID)
DTI.X <- DTI$cca[DTI.ind, ]
DTI.y <- DTI$pasat[DTI.ind]
DTI.X <- DTI.X[-59, ]
DTI.y <- DTI.y[-59]
DTI.X <- DTI.X[!is.na(DTI.y), ]
DTI.y <- DTI.y[!is.na(DTI.y)]

# Define server logic required to summarize and view the selected dataset
shinyServer(function(input, output, session) {
#  includeCSS('McLean.css')
  # Return the requested dataset
  datasetInput <- reactive({
    switch(input$dataset,
           "Emissions" = as.matrix(read.table('PM_flm.dat',sep=',',header=TRUE)[,2:42]),
           "DTI" = DTI.X,
           "Tecator" = tecator$absorp.fdata[[1]],
           "Temperature" = t(CanadianWeather$dailyAv[, , 1]),
           "Precipitation" = t(CanadianWeather$dailyAv[, , 2]))
  })
  
  response <- reactive({
    switch(input$dataset,
           "Emissions" = drop(read.table('PM_flm.dat',sep=',',header=TRUE)[,1]),
           "DTI" = DTI.y,
           "Tecator" = tecator$y[, 1],
           "Temperature" = rowMeans(t(CanadianWeather$dailyAv[, , 2])),
           "Precipitation" = rowMeans(t(CanadianWeather$dailyAv[, , 1])))
  })
  
#   N.all <- reactive({
#     switch(input$dataset,
#            "Emissions" = 157,
#            "Temperature" = 35,
#            "Precipitation" = 35)
#     })
#  seed <- reactive({sample.int(.Machine$integer.max, 1)})
  
  
#   n.max <- switch(get(input$dataset),
#                   'Emissions'=80, 'Temperature'=5, 'Precipitation'=5)
  
#   datasetInput <- reactiveValues()
#   datasetInput$X <-      switch(input$dataset,
#                 "Emissions" = read.table('PM_flm.dat',sep=',',header=TRUE)[,2:42],
#                 "Temperature" = t(CanadianWeather$dailyAv[, , 1]),
#                 "Precipitation" = t(CanadianWeather$dailyAv[, , 2]))
#   datasetInput$seed = sample.int(.Machine$integer.max, 1)
#     datasetInput$response <- reactive({
#       switch(input$dataset,
#              "Emissions" = drop(read.table('PM_flm.dat',sep=',',header=TRUE)[,1]),
#              "Temperature" = rowMeans(t(CanadianWeather$dailyAv[, , 2])),
#              "Precipitation" = rowMeans(t(CanadianWeather$dailyAv[, , 1])))
#     })
  
  output$xbfs <- renderUI({
    N <- length(response())
    xbfmax <- floor(sqrt(N))
    val <- floor((xbfmax+5)/2)
    if (xbfmax==5){
      numericInput("FGAMxnbf", "Number of basis functions to use for x-axis basis for FGAM:", value = val, min = 5, max = xbfmax, step = 0)
    }else{
      sliderInput("FGAMxnbf", "Number of basis functions to use for x-axis basis for FGAM:", value = val, min = 5, max = xbfmax, step = 1)
    }
    
  })
  
  output$tbfs <- renderUI({
    input$FGAMxnbf
    N <- length(response())
    tbfmax <- floor(sqrt(N))
    tbfmax <- tbfmax+(N-tbfmax^2)%/%tbfmax
    val <- floor((tbfmax+5)/2)
    sliderInput("FGAMtnbf", "Number of basis functions to use for t-axis basis for FGAM:", value = val, min = 5, max = tbfmax, step = 1)
  })
  
  output$n.test <- renderUI({
   # predbool=input$pred
    #save(predbool, 'temp.RData')
    
    if(!input$pred)
      return()

    N <- length(response())   # N.all()
    Nmax <- N-input$FGAMxnbf*input$FGAMtnbf  # switch(input$dataset,
              #        'Emissions'=80, 'Temperature'=5, 'Precipitation'=5)
    Nstart <- floor(Nmax/2)
    if(Nmax <= 0)
      return()
   # observe({updateNumericInput(session, 'dataset', input$dataset)})
#    numericInput("n.testing", "Number samples to use for test set:", 1, min = 1, max = Nmax, step = 1)
    sliderInput("n.testing", "Number samples to use for test set:", Nstart, min = 1, max = Nmax, step = 1)

#     print(input$n.testing)
#     helpText("Note: Selecting Zero for the Proportion of test samples results in one sample being used for testing.",
#              "on the full dataset.")
  })

  # Generate a summary of the dataset
  output$summary <- renderPrint({
    #dataset <- datasetInput()
    cat(paste(
    switch(input$dataset,
           "Emissions" = "X(t): Velocities for 157 trucks over 40 seconds of driving",
           "DTI" = "X(t): One-dimensional brain image from diffusion tensor imaging",
           "Tecator" = "X(t): Absorbances of ground pork samples for different wavelengths",
           "Temperature" = "X(t): Ave. daily temperature at 35 weather stations for one year",
           "Precipitation" = "X(t): Ave. daily rainfall at 35 weather stations for one year"),
    switch(input$dataset,
           "Emissions" = "Y: Particulate Matter exhaust emissions at end of 40 seconds",
           "DTI" = "Y: Score on cognitive test",
           "Tecator" = "Y: Fat content of ground pork sample",
           "Temperature" = "Y: Average yearly rainfall",
           "Precipitation" = "Y: Average yearly temperature"),
    sep='\n'))
  })
  
  # Show the first "n" observations
  output$matPlot <- renderPlot({
    data <- datasetInput()
    y <- response()
    nc <- input$obs #number of curves to plot
    #sset.seed(seed())
    use <- sample(nrow(data), nc)
    cols <- brewer.pal(nc, 'Set3')
    
    xlab <- switch(input$dataset,
           "Emissions" = "Time (seconds)",
           "Temperature" = "Day",
           "Precipitation" = "Day")
    
    matplot(t(data), col=rgb(0,0,0,.1), lty=1, type='l', xlab = xlab, ylab = 'X(t)')
    matlines(t(data[use,]), col=cols, lty=1, lwd=2)
    
#     matplot(t(resMC$data$Xsparse), col=rgb(0,0,0,.1), pch=19, cex=.5, 
#             ylim = ylim, bty="n",
#             xlab="Hour [t]", ylab="Log-Bid [X(t)]")
#     matlines(t(resMC$means$X), col=rgb(0,0,0,.1), lty=1)
#     
#     matplot(x=seq(0, 10, l=50), X10$X, type="l", lwd=2.5,  lty=1, ylim=range(X10$Xsparse, na.rm=TRUE),
#             xlab=expression(italic(t)), ylab=expression(italic(X(t))), main="a)", col = brewer.pal(8, 'Set3'),
#             bty="n", cex.lab = 1.5, cex.main = 2.1)
  })
  
  output$FLMplot <- renderPlot({
    if(!input$fitaction){
      return('')
    }
    input$spFLMsp
    input$FLMnbf
    isolate({
    dat <- list()
    dat$X <- datasetInput()
    dat$y <- response()
    N <- nrow(dat$X)
    J <- ncol(dat$X)
    if(isolate(input$pred)){
      # observe({print(input$test.prop)})
      # n.samps <- switch(as.character(input$test.prop), '0'=1, ceiling(N*input$test.prop))
#       if(is.null(input$n.testing)){
#         n.samps <- 1
#       }else{
#         n.samps <- input$n.testing  # switch(as.character(input$n.testing), '0'=1, input$n.testing) 
#       }
      test.samps <- sample(N, isolate(input$n.testing))
      assign('test.samps', test.samps, globalenv())
      dat$X <- dat$X[-test.samps, ]
      dat$y <- dat$y[-test.samps]
     # cat('hi', sum(dat$y))
      dat$L.X <- dat$X/J
    #  print(class(dat$L.X))
      dat$tmat <- matrix(seq(0, 1, l=J), N-isolate(input$n.testing), J, byrow=TRUE)
      #cat('flm', length(dat$y), dim(dat$L.X), dim(dat$tmat))
    }else{
      dat$L.X=matrix(1/J, N, J)*dat$X
      dat$tmat <- matrix(seq(0, 1, l=J), N, J, byrow=TRUE)
    }
    if (input$sptypeFLM == 'Fixed')
      assign('sp', input$FLMsp, globalenv())
    dat$k <- input$FLMnbf
    #fit <- with(dat, fgam(y~lf(X, splinepars=list(k=k), presmooth=FALSE), sp=exp(get('sp', env=globalenv()))))


    fit <- with(dat, switch(input$sptypeFLM,
                            'REML'=gam(y~s(tmat, by=L.X, k=k), method = 'REML'),
                            'GCV'=gam(y~s(tmat, by=L.X, k=k), method = 'GCV.Cp'),
                            'Fixed'=gam(y~s(tmat, by=L.X, k=k), method = 'REML', sp=exp(get('sp', env=globalenv())))))
    assign('fitFLM', fit, globalenv())
    #observe({print('FLMfit')})
    plot.gam(fit, rug=FALSE, main=expression(hat(beta)(t)), xlab='t', ylab='')
    })
    #plot(k,sp,main=class(sp))

  })
  
  output$FGAMplot <- renderPlot({
    if(!input$fitaction){
      return('')
    }
    isolate({
    dat <- list()
    dat$X <- datasetInput()
    dat$y <- response()
    N <- nrow(dat$X)
    J <- ncol(dat$X)

   # temp <- input$test.prop*3
    if(isolate(input$pred)){
     # observe({print(input$n.testing)})
    #  test.samps <- sample(N, input$n.testing)
     # invisible(input$n.testing)
      get('test.samps', globalenv())
    #  n.samps <- length(test.samps)
      dat$X <- dat$X[-test.samps, ]
      dat$y <- dat$y[-test.samps]
      dat$L <- matrix(1/J, N-isolate(input$n.testing), J)
      dat$tmat <- matrix(seq(0, 1, l=J), N-isolate(input$n.testing), J, byrow=TRUE)
      #print(sum(dat$y[test.samps]))
    # cat('fgam', length(dat$y), dim(dat$X), dim(dat$L), dim(dat$tmat))
    }else{
      dat$L <- matrix(1/J, N, J)
      dat$tmat <- matrix(seq(0, 1, l=J), N, J, byrow=TRUE)
    }
    #print(dim(dat$L))
    #print(dim(dat$tmat))
    if (input$sptypeFGAM == 'Fixed'){
      assign('xsp', input$FGAMxsp, globalenv())
      assign('tsp', input$FGAMtsp, globalenv())
    }
    dat$k <- c(input$FGAMxnbf, input$FGAMtnbf)

#     fit <- with(dat, fgam(y~af(X, splinepars=list(k=k), presmooth=FALSE), 
#                           sp=c(exp(get('xsp', env=globalenv())), exp(get('tsp', env=globalenv())))))
      fit <- with(dat, switch(input$sptypeFGAM,
                              'REML'=gam(y~te(X, tmat, by=L, k=k), method = 'REML'),
                              'GCV'=gam(y~te(X, tmat, by=L, k=k), method = 'GCV.Cp'),
                              'Fixed'=gam(y~te(X, tmat, by=L, k=k), method = 'REML',
                                         sp=c(exp(get('xsp', env=globalenv())), exp(get('tsp', env=globalenv()))))))
    assign('fitFGAM', fit, globalenv())
#     switch(input$plottype,
#       'persp'=refund::vis.fgam(fit, plot.type=input$plottype),
#                                'contour'=vis.gam(fit, plot.type=input$plottype, main = expression(hat(F)(x,t),
#                                                   xlab='x', ylab='t')))
    })
    phi <- ifelse(is.null(input$phi), 30, input$phi)
    theta <- ifelse(is.null(input$phi), 10, input$theta)

    switch(input$plottype,
      'persp'=vis.gam(fit, plot.type=input$plottype, phi=phi, theta=theta, ticktype='detailed', 
                                main = expression(hat(F)(x,t)), zlab='', xlab='x', ylab='t'),
                               'contour'=vis.gam(fit, plot.type=input$plottype, main = expression(hat(F)(x,t)),
                                                  xlab='x', ylab='t'))

#vis.gam(fit, plot.type=input$plottype, main = expression(hat(F)(x,t), xlab='x', ylab='t'))    
    #plot(k,sp,main=class(sp))
  })
  
  output$prederror <- renderPrint({
    input$fitaction
    isolate({
    if(!input$pred || is.null(input$n.testing))
      return(cat('Select number of test samples and click the \'Fit Models\' button'))
    #cat(input$n.testing)
 
    dat <- list()
    dat$X <- datasetInput()
    min.X <- min(dat$X)
    max.X <- max(dat$X)
    dat$y <- response()
    N <- nrow(dat$X)
    J <- ncol(dat$X)
    #temp <- input$test.prop*3
    # invisible(input$pred)
    # invisible(input$n.testing)
    
    get('fitFGAM', globalenv())
    get('fitFLM', globalenv())
    get('test.samps', globalenv())

     testX <- dat$X[test.samps, , drop=FALSE]
     testX[testX < min.X] <- min.X
     testX[testX > max.X] <- max.X

     L <- matrix(1/J, input$n.testing, J)
     tmat <- matrix(seq(0, 1, l=J), input$n.testing, J, byrow=TRUE)
    L.X <- L*testX
    newd <- list(tmat=tmat, L.X=L.X)
     predFLM <- predict(fitFLM, newd)
    if (input$n.testing==1){ # bug in mgcv here. something to do with qr decomp
      tmat <- rbind(tmat, tmat)
      L <- rbind(L, L)
      testX <- rbind(testX, testX)
    }
    newd <- list(tmat=tmat, L=L, X=testX)
    predFGAM <- predict(fitFGAM, newd)[1:input$n.testing]
   # cat(dim(testX), dim(L), dim(L.X), dim(tmat), length(input$n.testing))
   # cat(predFLM)
    
#     cat('Number of test samples = ',length(test.samps),'.\nRMSE for FGAM = ', 
#         signif(sqrt(mean((dat$y[test.samps]-predFGAM))), 4),
#         '.\nRMSE for FLM = ', signif(sqrt(mean((dat$y[test.samps]-predFLM))), 4))
    cat('RMSE for FGAM = ', 
        signif(sqrt(mean((dat$y[test.samps]-predFGAM)^2)), 4),
        '\nRMSE for FLM = ', signif(sqrt(mean((dat$y[test.samps]-predFLM)^2)), 4))
    })
#    cat(predFLM, predFGAM, sep='\n')
  })
  
  output$testres <- renderPrint({
    if(!input$fitaction)
      return(cat('Select tests to use and hit \'Fit Models\' button'))
    isolate({
    X <- datasetInput()
    y <- response()
    N <- nrow(X)
    J <- ncol(X)
    
    if (input$htest==0 || length(input$tests)==0)
      return(cat('Select tests to use using check boxes on left'))
  
    if('RLRT1'%in%input$tests){
      res <- TestKnowEqualVC(y, X, tvals=seq(0, 1, l=J), family=gaussian(),  # 
                   splinepars=list(k=c(input$FGAMxnbf, input$FGAMtnbf), m=list(c(2,2),c(2,2), extendXrange=.0)),
                   REML=TRUE,oracle=0,formratio=FALSE,sig2FLM=NULL,ratio=NULL,
                   psplines=TRUE)[1]
      if(is.na(res))
        res <- 1
      cat('Using RLRT1:\np-value is approximately ', res,'\n')
    }
    if ('RLRT2'%in%input$tests){
      res <- testfun(y, X, tvals=seq(0, 1, l=J), family=gaussian(),  # 
                             splinepars=list(k=c(input$FGAMxnbf, input$FGAMtnbf), m=list(c(2,2),c(2,2), extendXrange=.0)),
                             REML=TRUE,oracle=0,formratio=FALSE,sig2FLM=NULL,ratio=NULL,
                             psplines=TRUE)
      if (is.na(res[1]) && is.na(res[2])){
        res <- 1
      }else if (is.na(res[1])){
        res <- res[2]
      }else if (is.na(res[2])){
        res <- res[3]
      }else{
        res <- signif(max(res[1]*2, 2*res[2]),2)
      }
      cat('Using RLRT2:\np-value is approximately ', res,'\n')
    }
    if ('Bootstrap' %in%input$tests) {
      bootres <- TestTwoVCBOOT(y,X,tvals=seq(0,1,l=ncol(X)), family=gaussian(),
                                         splinepars=list(k=c(input$FGAMxnbf, input$FGAMtnbf), 
                                                         m=list(c(2,2),c(2,2)), extendXrange=.01),
                                         REML=TRUE, oracle=0, formratio=FALSE, sig2FLM=NULL, ratio=NULL,
                                         psplines=TRUE, removeCon = FALSE, BOOT=FALSE, n.sims.boot=5000)[1]
    
        cat('Using bootstrap:\np-value is approximately ',bootres)
    }
    })
  })
  
})