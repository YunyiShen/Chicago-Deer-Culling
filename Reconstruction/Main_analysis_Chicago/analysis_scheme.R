#ttt = ReCAP::analysisScheme(Chicago_RES$mcmc.objs,Assumptions,c(8,3),16,matrix(1,11,17))

require(dplyr)
library(ggplot2)

others = 1
Female_adult_weight = seq(0.5,2,0.5)

mean_dynamic = matrix(NA,length(Female_adult_weight),17)
sd_dynamic = matrix(NA,length(Female_adult_weight),17)

pp <- list()

for(i in 1:length(Female_adult_weight)){
  weight = Female_adult_weight[i]
  Harvest_weight_vector = c(others,others,rep(weight,6),rep(others,3))
  
  Harvest_matrix_at_this_level = matrix(Harvest_weight_vector,nrow = 11,ncol = 17)
  
  Scheme_temp = ReCAP::analysisScheme(Chicago_RES$mcmc.objs,
                                      Assumptions,c(8,3),
                                      16,Harvest_matrix_at_this_level, skip = c(0,0,0,0,0,0,0,0,0,0, 0,0,0, 0,0,0,0),quota = F)
  
  #Scheme_temp = ReCAP::analysisScheme_SimpleDD(Chicago_RES$mcmc.objs,
  #                                    Assumptions,c(8,3),
  #                                    16,Harvest_matrix_at_this_level, skip = c(0,0,0,0,0,0,1,0,0,1, 0,0,1, 0,0,1,0),
  #                                    quota = F, K=2000)
  
  sum_temp = sapply(Scheme_temp,function(w){
    colSums( w,na.rm = T)
  })
  
  sum_temp <- t(sum_temp)
  
  mean_dynamic_temp = Reduce("+",Scheme_temp)/length(Scheme_temp)
  
  mean_dynamic_temp_sum = colSums(mean_dynamic_temp)
  
  mean_dynamic[i,] = mean_dynamic_temp_sum
  
  plot_temp <- data.frame(year = 1992:2008, 
                          population = mean_dynamic_temp_sum, 
                          lw = apply(sum_temp,2,quantile, 0.05),
                          hi = apply(sum_temp,2,quantile, 0.95))
  
  pp[[i]] <- ggplot(data = plot_temp,aes(x = year,y=population))+
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin=lw, ymax=hi), width=.1) +
    ylab("Reconstructed Population") +
    #theme(text = element_text(size=18), axis.text.x = element_text(size = 16))+
    geom_hline(yintercept=500, linetype="dashed")+
    geom_hline(yintercept=250, linetype="dashed")
  
}

persp(0:16+1992, Female_adult_weight, t(mean_dynamic), phi = 30, theta = 135,
      xlab = "year", ylab = "weight on female adult", zlab = "Population size"
      #main = "Popoulation size after harvest"
)
contour(1:17+1991, Female_adult_weight, t(mean_dynamic))

last_five_years = mean_dynamic[,17-(0:4)] %>%
  rowMeans()

first_three_years = mean_dynamic[,1:3] %>%
  rowMeans()

plot(Female_adult_weight,last_five_years)

overall_lambda = last_five_years/first_three_years
plot(Female_adult_weight,overall_lambda)


lattice::wireframe(mean_dynamic, shade = TRUE,
                   xlab = "Weight", ylab = "Year",zlab = "Population",
                   aspect = c(16/17, 0.8),
                   light.source = c(10,10,10),
                   screen = list(z = 250, x = -60 , y=-0),
                   drape = T,
                   colorkey = T)


ggarrange(plotlist = pp, labels = c("weight 0.5","weight 1","weight 1.5","weight 2"), label.x = 0.2)
ggsave("./noDD_4case.pdf", width = 8, height = 5, scale = 0.9)
