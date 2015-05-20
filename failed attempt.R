library('shiny')
# read this: http://alexanderetz.com/2015/04/15/understanding-bayes-a-look-at-the-likelihood/


shinyApp(
  ui = shinyUI(fluidPage(
    sidebarLayout(
      sidebarPanel(
        numericInput("PS", label= "Prior Successes",
                     value=1),
        numericInput("PF", label="Prior Failures",
                     value=1),
        sliderInput("k", label = "k",
                    min = 1, max = 50, value = 7, step = 1),
        sliderInput("N", label = "N",
                    min = 1, max = 50, value = 10, step = 1),

      mainPanel(
        plotOutput("LRplot", height="500px")
        )
      )
   )),
   server = function(input, output) {
      plot.LR <- function(PS, PF) {
        curve(dbeta(x, PS, PF), xlim=c(0,1), ylab= "Density", xlab= "Probability of success", las=1, 
              main="Prior distribution for binomial data",lwd=3, cex.axis=1.5,
              cex.lab=1.5, cex.main=1.5)
        
       
      }
      
      #get.LR <- function(k, N, p1, p2) {
       # MLE <- dbinom(k, N, k/N)
        #L1 <- dbinom(k, N, prob = p1) / MLE
        #L2 <- dbinom(k, N, prob = p2) / MLE
        #L1 / L2
      #}
      
      output$LRplot <- renderPlot({
        PS <- input$PS
        PF <- input$PF
        plot.LR(PS, PF)
      })
      
      #output$LRatio1 <- renderText({
       # LR_12 <- get.LR(input$k, input$N, input$p1, input$p2)
        #paste("L1 / L2: ", round(LR_12, 3))
      #})
      
      #output$LRatio2 <- renderText({
       # LR_21 <- 1 / get.LR(input$k, input$N, input$p1, input$p2)
        #paste("L2 / L1: ", round(LR_21, 3))
      #})
    }
)