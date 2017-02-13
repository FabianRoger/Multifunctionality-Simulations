
library(multifunc)
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(RColorBrewer)
library(corrplot)



source("helpers.R")


##############################################################################################################################
                                                        # ui # 
##############################################################################################################################

# Define UI 
ui <- shinyUI(fluidPage(
   
#### Application title #####

  fluidRow(
    column(width = 6, offset = 3,
   titlePanel("Multifunctionality simulations")
   )
   ),
  
  br(),
   
############################ Sidebar ############################

  sidebarPanel(
    
#### Slider for number od species ####
    fluidRow(    
         sliderInput("specnum",
                     "Number of species",
                     min = 1,
                     max = 50,
                     value = 10)),

#### Slider for number of functions ####
   fluidRow(  
         sliderInput("funcnum",
                     "Number of functions",
                     min = 1,
                     max = 50,
                     value = 10)),
hr(),
br(),

#### choose distribution and specify parameters ####

  fluidRow(
    
#### choose distribution 

     column(6,
            selectInput("distribution",
                               "probability distribution",
                               list("uniform" = "runif",
                                    "normal" = "rnorm",
                                    "binary" = "rbinom",
                                    "beta" = "rbeta"),
                               selected = "runif",
                        multiple = FALSE),
            
#### specify standardization method
            selectInput("standardization",
                        "standardization method",
                        list("by maximum" = "max",
                             "between [0,1]" = "unit"),
                        selected = "unit",
                        multiple = FALSE),
            
            
#### specify seed
            textInput("seed",
                        "set seed",
                        value = "")),


#### specify parameters conditional on distribution 

     column(6,
            
            # parameters for uniform distribution
            conditionalPanel(
              condition = "input.distribution == 'runif'",
              textInput("min","minimum",
                           value = 0),
              textInput("max","maximum",
                           value = 1)),
            
            # parameters for normal distribution
            conditionalPanel(
              condition = "input.distribution == 'rnorm'",
              textInput("mean","mean",
                           value = 0.5),
              textInput("sd","standard deviation",
                           value = 0.1)),
            
            # parameters for binary distribution
            conditionalPanel(
              condition = "input.distribution == 'rbinom'",
              textInput("size","size",
                        value = 1),
              textInput("prob","probability",
                        value = 0.4)),
            
            # parameters for beta distribution
            conditionalPanel(
              condition = "input.distribution == 'rbeta'",
              textInput("shape1","shape1",
                        value = 3),
              textInput("shape2","shape2",
                        value = 2)))),

br(),
  
#### plot specified distribution ####

  fluidRow( plotOutput("ProbDens")),

hr(),
br(),

#### choose method for diveristy effect and specify conditional parameters ####
  
  fluidRow(
   
#### choose method
   column(6,
          selectInput("method",
                      "diversity effect",
                      list("none" = "av",
                           "complementarity" = "comp"),
#                           "selection" = "sel"),
                      selected = "av",
                      multiple = FALSE)),

#### choose conditional parameters

          column(6,
          
          # parameters method complementarity
          conditionalPanel(
            condition = "input.method == 'comp'",
            textInput("CF","complementarity factor",
                      value = 2),
            textInput("r","complementarity rate",
                      value = 1),
            uiOutput("functionlist_comp")
                        ))),
          
          # parameters method selection
         # conditionalPanel(
        #    condition = "input.method == 'sel'",
        #    uiOutput("functionlist_sel"),
        #    textInput("selfac","selection factor",
        #              value = 1.02))
        #  ),

#### plot complementarity factor if method = comp ####

  fluidRow(
    conditionalPanel(
      condition = "input.method == 'comp'",
    # Show a plot of the chosen distribution
          plotOutput("compfac")
   ))
),
  
  
############################ Main Panel ############################

  mainPanel(
    
#### Action button to draw function values ####
    
    fluidRow(
      column(12, offset = 3,
      actionButton("sample.func", "draw function values", width = 400))),

br(), #line break

#### plots of chosen function values ####
      
   fluidRow(
     
#### plot species - function - matrix
     column(6,
         plotOutput("SpecFuncMat")
         ),
     
#### plot function - correlation - matrix
     column(6,
         plotOutput("FuncCor")
         )
     ),

#### Action button to calculate results ####

fluidRow(
  column(12, offset = 3,
         actionButton("results", "calculate diversity ~ multifunctionality", width = 400))),

br(), #line break

#### plots for Average approach ####
   
   fluidRow(
     
### plots for diveristy single function vlaues 
     column(6,
            plotOutput("SingleFunc")
            ),

### plots for diveristy average function vlaues 
     column(6,
            plotOutput("AvFunc")
            )
     ),

hr(), #horizontal line

#### plots for multi-threshold approach ####

   fluidRow(
     
### plots for distinct thresholds
     column(6,
            plotOutput("singleThresh")
            ),
     
### plots for multiple thresholds
     column(6,
            plotOutput("multipleThresh")
            )
          )

    ) #main panle
  ) #fluid Page
) #Shiny UI
      
   




##############################################################################################################################
                                          # server # 
##############################################################################################################################





# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  
#### function list for method complementarity ####
  
  output$functionlist_comp <- renderUI({
    funclist <- c("all", FunctionList(input$funcnum))
    selectInput("compfunc",
                "complementarity function",
                multiple = TRUE,
                funclist,
                "all")
  })
  
#### function list for method selection ####
  
#  output$functionlist_sel <- renderUI({
#    funclist <- FunctionList(input$funcnum)
#    selectInput("selfunc",
#                "selection function",
#                funclist,
#                "Func_01")
#  })
  
  
#### plot for complementarity factor (if method = comp) ####
  output$compfac <- renderPlot({
    DF <- data.frame(A = 1: input$specnum)
    DF$B <-  as.numeric(input$CF) * ( 1 - ( 1 - 1/as.numeric(input$CF )) * exp(1-DF$A^as.numeric(input$r)))
    
    ggplot(DF, aes(A,B))+
      geom_point()+
      labs(x = "diversity", y = "complementarity factor")+
      scale_x_continuous(breaks = c(1:(input$specnum +1 )))+
      scale_y_continuous(limits = c(0, as.numeric(input$CF)))+
      theme_bw()
    })
  
  
  
#### plot of probability density function ####
  ProbDens <- reactive({
    
    if (input$distribution=="rnorm") {
      
      mean <- as.numeric(input$mean)
      sd <- as.numeric(input$sd)
      
      Dist <- rnorm(10000, mean, sd)
      Lim <- quantile(Dist, c(0.001, 0.999))
      DF <- data.frame(B = signif(Lim,2))
    
    
      ggplot(DF, aes(x = B))+
        stat_function(fun = dnorm, args = list( mean = mean, sd = sd), 
                      colour = "red", size = 2)+
        labs( title = "probablity density function",
              x = "range of function values",
              y = "probablity")
    } else if (input$distribution=="rbeta") {
        
        shape1 <- as.numeric(input$shape1)
        shape2 <- as.numeric(input$shape2)
        
        Dist <- rbeta(10000, shape1, shape2)
        Lim <- quantile(Dist, c(0.001, 0.999))
        DF <- data.frame(B = round(Lim))
        
        
        ggplot(DF, aes(x = B))+
          stat_function(fun = dbeta, args = list( shape1, shape2), 
                        colour = "red", size = 2)+
          labs( title = "probablity density function",
                x = "range of function values",
                y = "probablity")
      } else if (input$distribution=="runif") {
        
        min <- as.numeric(input$min)
        max <- as.numeric(input$max)
        
        DF <- data.frame(B = c(0,1))
        
        
        ggplot(DF, aes(x = B))+
          stat_function(fun = dunif, args = list( min, max), 
                        colour = "red", size = 2)+
          labs( title = "probablity density function",
                x = "range of function values",
                y = "probablity")
      } else if (input$distribution=="rbinom") {
      
      size <- as.numeric(input$size)
      prob <- as.numeric(input$prob)
      
      DF <- data.frame(B = c(0,1))
      
      ggplot(DF, aes(x = B))+
        stat_function(fun = dbinom, args = list( size, prob), 
                      colour = "red", size = 2)+
        labs( title = "probablity density function",
              x = "range of function values",
              y = "probablity")
    }
    
    })
  
  output$ProbDens <- renderPlot({ProbDens()})
  

#################### calculate MF and plot results ###########################
  
  
  #set seed to recompute
  #observeEvent(input$calculate,{
  
  
  
######### calculate reactive values #########
  
  
  
#### draw function values from specified distribution ####
  
  FuncMat <- eventReactive(input$sample.func, {
    
    if(nchar(input$seed) > 0) {set.seed(as.numeric(input$seed))}
    #genrate function matrix
    
    if (input$distribution=="runif") {
      min <- as.numeric(input$min)
      max <- as.numeric(input$max)
      FuncMat <- FunctionValue(input$specnum,input$funcnum, "runif", min = min, max = max)
    } 
    
    if (input$distribution=="rnorm") {
      mean <- as.numeric(input$mean)
      sd <- as.numeric(input$sd)
      FuncMat <- FunctionValue(input$specnum,input$funcnum, "rnorm", mean = mean, sd = sd)
    } 
    
    if (input$distribution=="rbeta") {
      shape1 <- as.numeric(input$shape1)
      shape2 <- as.numeric(input$shape2)
      FuncMat <- FunctionValue(input$specnum,input$funcnum, "rbeta", shape1, shape2)
    } 
    
    if (input$distribution=="rbinom") {
      size <- as.numeric(input$size)
      prob <- as.numeric(input$prob)
      FuncMat <- FunctionValue(input$specnum,input$funcnum, "rbinom", size, prob)
    }
    
    FuncMat <- FuncMat %>% 
        group_by(Functions)
    
  })
    

  
#### calculate species Matrix with specified number of species ###
  
  SpecMat <- reactive({
    SpeciesMatrix(specnum = input$specnum, maxrep = 50)
  })
  
### calculate average multifunctionality with specified method ###
  
  AvFunc <- eventReactive( input$results ,{
    
    AvFunc <- AverageFunction(SpecMat(), 
                    FuncMat(),
                    method = input$method, 
                    CF = as.numeric(input$CF),
                    r = as.numeric(input$r),
                    compfunc = input$compfunc,
                    selfunc = input$selfunc, 
                    selfac = as.numeric(input$selfac))
    
    func.names <- as.character( unique( FuncMat()$Functions))
    
    if (input$standardization == "unit"){
      AvFunc[,func.names] <- apply(AvFunc[,func.names], 2, function(x) {(x - min(x)) / (max(x) - min(x))})
      } 
    
    if (input$standardization == "max"){AvFunc[,func.names] <- apply(AvFunc[,func.names], 2, function(x) {x/max(x)})}
    
    AvFunc$meanFunction <- rowMeans(AvFunc[,func.names])
    
    
        
        return(AvFunc)
    
        })
  
#### calculate slopes for multithreshold approach ###
  
  mixedThresh <- eventReactive( input$results, {
    
    # extract function names
    func.names <- as.character( unique( FuncMat()$Functions))
    
    AvFunc <- AvFunc()
    
    #AvFunc$meanFunction <- rowMeans(AvFunc[,func.names])
    
    getFuncsMaxed(AvFunc, func.names, threshmin=0.05, threshmax=0.99, 
                  prepend=c("Richness"), maxN=1)
  })

  
######### render Plots #########
  
#### plot species - function - matrix ####
   
   output$SpecFuncMat <- renderPlot({
     
     FuncMat <- FuncMat()
     
     if (input$standardization == "unit"){
       FuncMat <- FuncMat %>% 
       group_by(Functions) %>% 
       mutate(Funcval = (Funcval - min(Funcval)) / (max(Funcval) - min(Funcval)))}
     
     if (input$standardization == "max"){
       FuncMat <- FuncMat %>%
         group_by(Functions) %>% 
         mutate(Funcval = Funcval / (max(Funcval)))}
     
     # define colours
     col <- colorRampPalette(c("#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
        
     #plot function matrix
     SF_G <- ggplot(FuncMat, aes(Functions, Species, fill = Funcval))+
       geom_tile(colour = "black", size = 0.7)+
       annotate(geom = "text", x = (1:isolate(input$funcnum))+0.5, y = isolate(input$specnum+isolate(input$specnum)/10), label = unique(FuncMat()$Functions), angle = 30)+
       annotate(geom = "text", x = -1, y = 1:isolate(input$specnum) , label = unique(FuncMat()$Species), colour = "#323232")+
       scale_fill_gradientn(colours = col(length(unique(FuncMat()$Funcval))), limits=c(0, 1))+
       theme_bw()+
       theme(
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         panel.grid.major = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.ticks = element_blank(), 
         legend.position = "right",
         plot.margin=unit(c(1, 1, 0, 0), "cm"))+
      guides(fill = guide_colorbar(barwidth = 2, barheight = 15,
                                    title.position = "top",
                                    direction = "vertical"))
       labs(x = "", y = "")
     
     
     
     SF_G <- ggplot_gtable(ggplot_build(SF_G))
     SF_G$layout$clip[SF_G$layout$name == "panel"] <- "off"
     plot(SF_G)
     })
  
#### plot function correlation matrix ####
   
   output$FuncCor <- renderPlot({
     
     FuncMat_wide <- FuncMat() %>% spread(Functions, Funcval)
     
     C_mat <- cor(FuncMat_wide[,-1])
     corrplot(C_mat, type = "lower", tl.col = "black", cl.ratio = 0.2, 
              cl.length = 11, number.cex = 0.6, addCoef.col = "#323232", diag = F, method="ellipse")})

  
       ###### Average approach ######
  
#### plot diveristy ~ single function values ####
   
   output$SingleFunc <- renderPlot({
     
     AvFunc_long <- gather(AvFunc()[, -c(1:isolate(input$specnum))], Function, FuncVal, -Richness)
     
     ggplot(AvFunc_long, aes(x = Richness, y = FuncVal))+
       geom_point(alpha = 0.2)+
       facet_wrap(~Function) +
       theme_bw(base_size=15)+
       stat_smooth(method="lm", colour="black", size=2) +
       xlab("\nSpecies Richness") +
       ylab("Value of Function\n") +
       theme(panel.grid = element_blank())
     })
  
#### plot diveristy ~ average function values ####
   
   output$AvFunc <- renderPlot({
     # extract function names
     func.names <- as.character( unique( FuncMat()$Functions))
     
     # add on the new (standardized) functions along with the averaged multifunctional index
    # AvFunc <- cbind(AvFunc(), getStdAndMeanFunctions(AvFunc(), func.names))
     
     #plot it
     ggplot(AvFunc(), aes(x=Richness, y=meanFunction))+
       geom_point(size=3, alpha =0.3)+
       theme_bw(base_size=15)+
       stat_smooth(method="lm", colour="black", size=2) +
       xlab("\nSpecies Richness") +
       ylab("Average Value of Standardized Functions\n")+
       scale_y_continuous(limits = c(0,1))
   })
     

#### plot diveristy ~ #functions for single thresholds ####
  
   output$singleThresh <- renderPlot({filter(mixedThresh(), as.character(thresh) %in% as.character(seq(0,1,0.1))) %>% 
       mutate(prct = paste(thresholds * 100, "%")) %>% 
       ggplot(., aes(x = Richness, y = funcMaxed))+
       geom_point(alpha = 0.1)+
       stat_smooth(method = "lm", colour = "red", se = F)+
       facet_wrap(~prct)+
       labs(x = "Species richness", y = "Number of function ≥ Threshold")+
       theme_bw()
     })
   
#### plot diveristy ~ #functions for multiple thresholds ####
  
   output$multipleThresh <- renderPlot({
     
     mixedLinearSlopes <- getCoefTab(funcMaxed ~ Richness, fun = lm,  data=mixedThresh(), 
                                   coefVar="Richness")
     
     colnames(mixedLinearSlopes) <- c("thresholds", "Estimate",  "Std. Error", "t value", "Pr(>|t|)")
     
     SlSum <- SlopeSummary(mixedLinearSlopes)
     #SlSum <- lapply(SlSum, function (y) median(y))
     
     SlSum <- data.frame(thresholds = unlist(SlSum))  %>%
       add_rownames(var = "label")  %>% 
       left_join(.,mixedLinearSlopes) %>% 
       mutate(label = ifelse(grepl("max", label), "max" ,ifelse (
         grepl("sign", label), "sign change", "min")))
   
     SlSum_label <- SlSum %>% group_by(label) %>% 
       summarise(minthresh = min(thresholds),maxthresh = max(thresholds)) %>% 
       left_join(SlSum) %>% 
       select(label, minthresh, maxthresh, Estimate) %>% 
       distinct()
     
     SlSum_label_same <- 
       SlSum_label %>% 
       mutate(same = ifelse(minthresh == maxthresh, "yes", "no")) %>% 
       filter(same == "yes")
     
     SlSum_label_not_same <- 
       SlSum_label %>% 
       mutate(same = ifelse(minthresh == maxthresh, "yes", "no")) %>% 
       filter(same == "no")
     
     
     
   p <-   ggplot(mixedLinearSlopes, aes(x=thresholds)) +
       geom_ribbon(fill="grey50", aes(x=thresholds*100, ymin=Estimate- 1.96*mixedLinearSlopes[["Std. Error"]],
                                      ymax=Estimate+1.96*mixedLinearSlopes[["Std. Error"]])) + 
       geom_point(aes(x=thresholds*100, y=Estimate)) +
       geom_point(data = SlSum, aes(x = thresholds*100, y = Estimate), size = 2, colour = "red")+
       ylab("Change in Number of Functions per Addition of 1 Species\n") + xlab("\nThreshold (%)") +
       geom_abline(intercept=0, slope=0, lwd=1, linetype=2) + 
       theme_bw(base_size=14)
   
   if(nrow(SlSum_label_not_same) > 0) {
    p <-  p +
   geom_label(data = SlSum_label_not_same, aes(x = maxthresh*100, y = Estimate, 
                                               label = paste(label," ", round(minthresh*100), 
                                                             "-",round(maxthresh*100), " %", sep ="" )),
              nudge_y = 0.1, nudge_x = -5)}
   
   if(nrow(SlSum_label_same) > 0) {
     
     p <- p +
     geom_label(data = SlSum_label_same, aes(x = maxthresh*100, y = Estimate, 
                                             label = paste(label," ",round(maxthresh*100), " %", sep ="" )),
                nudge_y = 0.1, nudge_x = -5)}
   
   p
   
  
     })
  
})

shinyApp(ui = ui, server = server)

