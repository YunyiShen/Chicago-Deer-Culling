#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#install_github("https://github.com/YunyiShen/ReCAP",ref = "vital-rate-as-data")
library(shiny)
library(coda)
library(svMisc)
library(RcppArmadillo)
library(Rcpp)

if (!file.exists("R-lib")) {
    dir.create("R-lib")
}
# Unfortunately, there's no way to get deployapp() to ignore this directory, so
# make sure to remove it locally before you call deployapp(). This can be done
# with:
#   unlink("pkgInst/R-lib", recursive = TRUE)

# You may also need to restart R before calling deployapp(), because calling
# runApp() will modify your libpath (below), which can confuse deployapp().

# Add ./R-lib/ to the libPaths
.libPaths( c(normalizePath("R-lib/"), .libPaths()) )

# Install the package if needed.
if (!do.call(require, list("ReCAP"))) {
    install.packages("ReCAP_0.1.1.tar.gz", repos = NULL, type = "source")
}

# Instead of `library(myPackage)`, we'll use do.call, to evade deployapp's
# checks for packages installed locally from source.
do.call(library, list("ReCAP"))



# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Complex 1 Population Dynamics under Different Harvest Weight"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            sliderInput("weight_doe",
                        "Weight on adult female",
                        min = 0.1,
                        max = 3,
                        value = 1),
            sliderInput("weight_other",
                        "Weight on other classes",
                        min = 0.1,
                        max = 3,
                        value = 1),
            radioButtons("quota", "Fixed:",
                         choices = list("Proportion" = FALSE, "Quota" = TRUE),
                         selected = FALSE),
            checkboxGroupInput("sumover",
                               ("Age classes to sum over?"),
                               choices = c("Female Fawn" = 1,
                                              "Female Yearling" = 2,
                                              "Female Adult" = 3,
                                              "Male Fawn" = 4,
                                              "Male Yearling" = 5,
                                              "Male Adult" = 6),
                               selected = 1)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("popdyn")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    load("4Harv_8Fec_11Surv_non_equal.RData")
    #load("4harv_3fec_6surv_equalvar.RData")
    output$popdyn <- renderPlot({
        # generate bins based on input$bins from ui.R
        others = input$weight_other
        weight = input$weight_doe
        Harvest_weight_vector = c(others,others,rep(weight,6),rep(others,3))

        Harvest_matrix_at_this_level = matrix(Harvest_weight_vector,nrow = 11,ncol = 17)

        Scheme_temp = analysisScheme(Chicago_RES$mcmc.objs,
                                            Assumptions,c(8,3),
                                            16,Harvest_matrix_at_this_level,quota = input$quota)

        sum_temp = lapply(Scheme_temp,function(w){
            temp = w
            temp[3,] = colSums(w[3:8,])
            temp = temp[-(4:8),]
        })


        sum_temp = lapply(sum_temp,function(w,sumover){
            colSums( matrix( w[ as.numeric( sumover),],ncol = 17))

        },input$sumover)
        sum_temp = Reduce(rbind,sum_temp)

        sum_temp_mean = colMeans(sum_temp,na.rm = T)
        CI_low = apply(sum_temp,2,quantile,0.025, na.rm = T)
        CI_high = apply(sum_temp,2,quantile,0.975, na.rm = T)

        plot_data = data.frame(year = 1:17+1991,
                               pop_mean = sum_temp_mean,
                               CI_low,
                               CI_high)
        require(ggplot2)
        ggplot(data = plot_data,aes(x = year,y=pop_mean))+
            geom_line() +
            geom_point() +
            geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=.1) +
            ylab("Reconstructed Population") +
            theme(text = element_text(size=18), axis.text.x = element_text(size = 16))

    })
}

# Run the application
shinyApp(ui = ui, server = server)
