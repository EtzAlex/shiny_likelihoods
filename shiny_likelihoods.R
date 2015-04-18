library('shiny')
# read this: http://alexanderetz.com/2015/04/15/understanding-bayes-a-look-at-the-likelihood/


shinyApp(
  ui = shinyUI(fluidPage(
    sidebarLayout(
      sidebarPanel(
        sliderInput("p1", label = "p1",
                    min = 0, max = 1, value = 0.7, step = 0.01),
        sliderInput("p2", label = "p2",
                    min = 0, max = 1, value = 0.5, step = 0.01),
        sliderInput("k", label = "k",
                    min = 1, max = 50, value = 7, step = 1),
        sliderInput("N", label = "N",
                    min = 1, max = 50, value = 10, step = 1),
        br(),
        br(),
        textOutput("LRatio1"),
        br(),
        textOutput("LRatio2"),
        tags$head(tags$style("#LRatio1, #LRatio2 { font-size: 25px; }"))
        ),
      mainPanel(
        plotOutput("LRplot", height="500px")
        )
      )
   )),
   server = function(input, output) {
      plot.LR <- function(k, N, p1, p2) {
        # adapted from blog mentioned above
        MLE <- dbinom(k, N, k/N)
        L1 <- dbinom(k, N, prob = p1) / MLE
        L2 <- dbinom(k, N, prob = p2) / MLE
        
        curve((dbinom(k, N, x) / max(dbinom(k, N, x))), xlim = c(0, 1),
              ylab = "Likelihood", xlab = "Probability of correct answer",
              las = 1, main = "Likelihood function for binomials", lwd = 3,
              cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
        
        points(p1, L1, cex = 2, pch = 21, bg = "cyan")
        points(p2, L2, cex = 2, pch = 21, bg = "cyan")
        lines(c(p1, p2), c(L1, L1), lwd = 3, lty = 2, col = "cyan")
        lines(c(p2, p2), c(L1, L2), lwd = 3, lty = 2, col = "cyan")
        abline(v = k/N, lty = 5, lwd = 1, col = "grey73")
      }
      
      get.LR <- function(k, N, p1, p2) {
        MLE <- dbinom(k, N, k/N)
        L1 <- dbinom(k, N, prob = p1) / MLE
        L2 <- dbinom(k, N, prob = p2) / MLE
        L1 / L2
      }
      
      output$LRplot <- renderPlot({
        k <- input$k
        N <- input$N
        if (k <= N) plot.LR(k, N, input$p1, input$p2)
      })
      
      output$LRatio1 <- renderText({
        LR_12 <- get.LR(input$k, input$N, input$p1, input$p2)
        paste("L1 / L2: ", round(LR_12, 3))
      })
      
      output$LRatio2 <- renderText({
        LR_21 <- 1 / get.LR(input$k, input$N, input$p1, input$p2)
        paste("L2 / L1: ", round(LR_21, 3))
      })
    }
)