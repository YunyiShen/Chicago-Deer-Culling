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
library(ggplot2)
library(jsonlite)

source("helper.R")

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
    install.packages("ReCAP_0.1.2.tar.gz", repos = NULL, type = "source")
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
                        "Weight on adult+yearling female",
                        min = 0.01,
                        max = 1,
                        value = 0.5),
            sliderInput("weight_buck",
                        "Weight on adult+yearling male",
                        min = 0.01,
                        max = 1,
                        value = .5),
            sliderInput("weight_fawn",
                        "Weight on fawns",
                        min = 0.01,
                        max = 1,
                        value = .5),
            sliderInput("K",
                        "Carrying capacity",
                        min = 1000,
                        max = 1800,
                        value = 1200),
            sliderInput("goal",
                        "Population goal",
                        min = 100,
                        max = 500,
                        value = 250),
            sliderInput("alpha",
                        "CI level",
                        min = 0.01,
                        max = 0.99,
                        value = 0.90),
            radioButtons("quota", "Fixed:",
                         choices = list("Proportion" = FALSE),#, "Quota" = TRUE),
                         selected = FALSE),
            checkboxGroupInput("sumover",
                               ("Age classes to sum over?"),
                               choices = c("Female Fawn" = 1,
                                              "Female Yearling" = 2,
                                              "Female Adult" = 3,
                                              "Male Fawn" = 4,
                                              "Male Yearling" = 5,
                                              "Male Adult" = 6),
                               selected = 1:6)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            navbarPage("Showing:",
                       tabPanel("Population",
                                    plotOutput("popdyn"),
                                downloadButton("downloadpop", "Download raw numerical result")
                                ),
                       tabPanel("Chance Reaching the Goal",
                                        plotOutput("p_ctrl"),
                                downloadButton("downloadp_ctr", "Download raw numerical result")
                                ),
                       tabPanel("More detailed numerical results",
                                selectInput("result", "Choose a result:",
                                            choices = c("harvest_median", 
                                                        "harvest_lower_ci", 
                                                        "harvest_upper_ci",
                                                        "post_harvest_median",
                                                        "post_harvest_lower_ci",
                                                        "post_harvest_upper_ci",
                                                        "non_harvest_change_median",
                                                        "non_harvest_change_lower_ci",
                                                        "non_harvest_change_upper_ci")),
                                
                                # Button
                                downloadButton("downloadres", "Download"),
                                tableOutput("table")
                                )
     
                ),
                       
           
           br(),
           br(),
           checkboxGroupInput("noskip",
                              ("Which year to cull?"),
                              choices = c("1992" = 1,"1993" = 2,
                                          "1994" = 3,"1995" = 4,
                                          "1996" = 5,"1997" = 6,
                                          "1998" = 7,"1999" = 8,
                                          "2000" = 9,"2001" = 10,
                                          "2002" = 11,"2003" = 12,
                                          "2004" = 13,"2005" = 14,
                                          "2006" = 15,"2007" = 16,
                                          "2008" = 17),
                              selected = 1:17,inline = TRUE),
           downloadButton("downloadsetting", "Download raw settings in JSON")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    load("4harv_8fec_6surv_equal.RData")
    rm(list = lsf.str())
    output$downloadsetting <- downloadHandler(
        filename = "settings.json",
        content = function(file) {
            write_json(toJSON(reactiveValuesToList(input)), file)
        }
    )


    #load("4harv_3fec_6surv_equalvar.RData")
    output$popdyn <- renderPlot({
        noskip = input$noskip
        skip = 1 - ((1:17) %in% noskip)
        buck = input$weight_buck
        doe = input$weight_doe
        fawn = input$weight_fawn
        K = input$K
        goal = input$goal
        Harvest_weight_vector = c(fawn,doe,rep(doe,6),fawn,rep(buck,2))
        
        Harvest_matrix_at_this_level = matrix(Harvest_weight_vector,nrow = 11,ncol = 17)
        
        Scheme_temp = analysisScheme_SimpleDD(Chicago_RES$mcmc.objs,
                                              Assumptions,c(8,3),
                                              16,Harvest_matrix_at_this_level,quota = input$quota,skip = skip,K=K)
        
        sum_temp = lapply(Scheme_temp,function(w){
            temp = w
            temp[is.na(temp)] = 0
            temp[3,] = colSums(w[3:8,],na.rm = T)
            temp = temp[-(4:8),]
        })
        
        sum_total = sum_temp
        sum_temp = lapply(sum_temp,function(w,sumover){
            colSums( matrix( w[ as.numeric( sumover),],ncol = 17),na.rm = T)
            
        },input$sumover)
        sum_total = Reduce(rbind,sum_total)
        sum_total[is.na(sum_total)] = 0
        sum_temp = Reduce(rbind,sum_temp)
        sum_temp[is.na(sum_temp)] = 0
        
        lev <- (1-input$alpha)/2
        
        sum_temp_mean = apply(sum_temp,2,median,na.rm = T)
        CI_low = apply(sum_temp,2,quantile,lev, na.rm = T)
        CI_high = apply(sum_temp,2,quantile,1-lev, na.rm = T)
        poster_prob_control = colMeans(sum_total<=goal,na.rm = T)
        plot_data = data.frame(year = 1:17+1991,
                               pop_mean = sum_temp_mean,
                               p_ctrl = poster_prob_control,
                               CI_low,
                               CI_high)
        output$downloadpop <- downloadHandler(
            filename = "population_est.csv",
            content = function(file) {
                write.csv(plot_data[,-3], file, row.names = FALSE)
            }
        )
        
        # generate bins based on input$bins from ui.R
        ggplot(data = plot_data,aes(x = year,y=pop_mean))+
            geom_line() +
            geom_point() +
            geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=.1) +
            ylab("Reconstructed Population") +
            theme(text = element_text(size=18), axis.text.x = element_text(size = 16))+
            geom_hline(yintercept=goal, linetype="dashed")

    })
    output$p_ctrl <- renderPlot({
        noskip = input$noskip
        skip = 1 - ((1:17) %in% noskip)
        buck = input$weight_buck
        doe = input$weight_doe
        fawn = input$weight_fawn
        K = input$K
        goal = input$goal
        Harvest_weight_vector = c(fawn,doe,rep(doe,6),fawn,rep(buck,2))
        
        Harvest_matrix_at_this_level = matrix(Harvest_weight_vector,nrow = 11,ncol = 17)
        
        Scheme_temp = analysisScheme_SimpleDD(Chicago_RES$mcmc.objs,
                                              Assumptions,c(8,3),
                                              16,Harvest_matrix_at_this_level,quota = input$quota,skip = skip,K=K)
        
        sum_total = lapply(Scheme_temp,function(w){
            temp = w
            temp[is.na(temp)] = 0
            temp[3,] = colSums(w[3:8,])
            temp = temp[-(4:8),]
        })
        
        sum_total = lapply(sum_total,colSums)
        sum_total = Reduce(rbind,sum_total)
        
        poster_prob_control = colMeans(sum_total<=goal,na.rm = T)
        plot_data = data.frame(year = 1:17+1991,
                               p_ctrl = poster_prob_control)
        output$downloadp_ctr <- downloadHandler(
            filename = "prob_reaching_goal.csv",
            content = function(file) {
                write.csv(plot_data, file, row.names = FALSE)
            }
        )
        ggplot(data = plot_data,aes(x = year,y=p_ctrl))+
            geom_line() +
            geom_point() +
            ylab("Posterior Probability Reaching the Goal")+
            ylim(0,1)+
            theme(text = element_text(size=18), axis.text.x = element_text(size = 16))
    })
    output$table <- renderTable({
        lev <- (1-input$alpha)/2
        what_we_want <- #reactive({
            switch(input$result,
                   "harvest_median" = list("Harvest",0.5), 
                   "harvest_lower_ci" = list("Harvest",lev), 
                   "harvest_upper_ci" = list("Harvest",1-lev),
                   "post_harvest_median" = list("Living",0.5),
                   "post_harvest_lower_ci" = list("Living", lev),
                   "post_harvest_upper_ci" = list("Living", 1-lev),
                   "non_harvest_change_median" = list("non_havest_change", .5),
                   "non_harvest_change_lower_ci" = list("non_havest_change", lev),
                   "non_harvest_change_upper_ci" = list("non_havest_change", 1-lev))
        #})
        
        ## run the simulation
        noskip = input$noskip
        skip = 1 - ((1:17) %in% noskip)
        buck = input$weight_buck
        doe = input$weight_doe
        fawn = input$weight_fawn
        K = input$K
        goal = input$goal
        Harvest_weight_vector = c(fawn,doe,rep(doe,6),fawn,rep(buck,2))
        
        Harvest_matrix_at_this_level = matrix(Harvest_weight_vector,nrow = 11,ncol = 17)
        Scheme_temp = analysisportion_simpleDDScheme_full(Chicago_RES$mcmc.objs,
                                              Assumptions,c(8,3),
                                              16,Harvest_matrix_at_this_level,skip = skip,K=K)
        the_table = extract_things_we_want(Scheme_temp, what_we_want)
        output$downloadres <- downloadHandler(
            filename = function() {
                paste(input$result, ".csv", sep = "")
            },
            content = function(file) {
                write.csv(the_table, file, row.names = TRUE)
            }
        )
        
        the_table
        
    })
}

# Run the application
shinyApp(ui = ui, server = server)
