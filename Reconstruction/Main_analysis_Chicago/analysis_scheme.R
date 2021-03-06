#ttt = ReCAP::analysisScheme(Chicago_RES$mcmc.objs,Assumptions,c(8,3),16,matrix(1,11,17))

require(dplyr)


others = 1
Female_adult_weight = seq(0.5,2.5,0.05)

mean_dynamic = matrix(NA,length(Female_adult_weight),17)
sd_dynamic = matrix(NA,length(Female_adult_weight),17)


for(i in 1:length(Female_adult_weight)){
  weight = Female_adult_weight[i]
  Harvest_weight_vector = c(others,others,rep(weight,6),rep(others,3))
  
  Harvest_matrix_at_this_level = matrix(Harvest_weight_vector,nrow = 11,ncol = 17)
  
  Scheme_temp = ReCAP::analysisScheme(Chicago_RES$mcmc.objs,
                                      Assumptions,c(8,3),
                                      16,Harvest_matrix_at_this_level,quota = F)

  mean_dynamic_temp = Reduce("+",Scheme_temp)/length(Scheme_temp)
  
  mean_dynamic_temp_sum = colSums(mean_dynamic_temp)
  
  mean_dynamic[i,] = mean_dynamic_temp_sum
  
  
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
