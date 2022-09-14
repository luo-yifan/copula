####################################################
#### Load libraries
####################################################
library(copula)
library(RCurl)
library(data.table)
library(stats)
library(dplyr)
library(VineCopula)
library(dglm)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(GGally)
library(tidyr)
library(zoo)
library(FAdist)
library(brms)
library(plot.matrix)
#library(BRugs)
#library(R2WinBUGS)

####################################################
#### Load data
####################################################
main.dir 	<- getwd() 
l2s.dir 	<- paste(main.dir, "/l2s_posterior/", sep = "") 	# L2SWBM simulations
data.dir 	<- paste(main.dir, "/historical_data/", sep = "") 	# historical data
plot.dir 	<- paste(main.dir, "/plots/", sep = "") 			# outputed plots

####################################################
#### Constants
####################################################

months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
lakes <- c("Lake Superior", "Lake Michigan-Huron", "Lake Erie")
S_SF 	<- 82103	 	# Lake Superior Surface Area
MH_SF 	<- 117250	 	# Lake Michigan & Huron Surface Area
E_SF 	<- 25740	 	# Lake Erie Surface Area
SPM31 	<- 2678400 		# Seconds per month with 31 days
SPM30 	<- 2592000 		# Seconds per month with 30 days
SPMF 	<- 2419200 		# Seconds in February
SPMFL 	<- 2505600 		# Seconds in February during Leap year 
horizon <- 28 			# copula simulation time horizon in years
gen 	<- 999 			# number of copula predictions per variable; BUGS gives vector length 999

## Monthly Climate Perturbations for each Component
#               |  Jan ~ Dec  
# Sup (factor)  |_____________
# MH  (factor)  |_____________
# Eri (factor)  |_____________
#monthly_climate_mean_perturb <- rep(list(rep(list(matrix(nrow = 3,ncol = 12,dimnames = list(c("Superior","Michigan-Huron","Erie"),months))),3)),4)

## 1 - no climate change  2 - baseline(self-generated)  3 - high-warming(via papers)  4 - Copula Validation Plots
# climate_data <- read.csv(file = paste(main.dir,"/Config.csv", sep = ""),header=F)
# for(s in 1:4){ ## Scenario
#   for(i in 1:3){ ## lake
#     for(j in 1:3){ ## component
#       for(m in 1:12){ ## month
#         if(s==4){
#           monthly_climate_mean_perturb[[s]][[j]][i,m] <- climate_data[36*(i-1)+12*(j-1)+m,1]
#         }else{
#           monthly_climate_mean_perturb[[s]][[j]][i,m] <- climate_data[36*(i-1)+12*(j-1)+m,s]
#         }
#       }
#     }
#   }
# }

################################################################
## Read in median water balance values from L2SWBM output     ##
## L2SWBM output runs from 1950 to 2019                       ##
## Dimension of each record is 840 x 1 (70 years x 12 months) ##
################################################################

####################################################
## Create empty list of vectors ('copulaData') to store L2SWBM output
# 1 - Lake Superior precipitation     4 - Lake Michigan & Huron precipitation     7 - Lake Erie precipitation
# 2 - Lake Superior evaporation       5 - Lake Michigan & Huron evaporation       8 - Lake Erie evaporation
# 3 - Lake Superior runoff            6 - Lake Michigan & Huron runoff            9 - Lake Erie runoff
copulaData <- vector(mode = "list", length = 9)
####################################################

## Set up file names that contain L2SWBM output
glwbDataFiles <- c(
  "superiorPrecip_GLWBData.csv",	"superiorEvap_GLWBData.csv",	"superiorRunoff_GLWBData.csv",
  "miHuronPrecip_GLWBData.csv",	"miHuronEvap_GLWBData.csv",		"miHuronRunoff_GLWBData.csv",
  "eriePrecip_GLWBData.csv",		"erieEvap_GLWBData.csv",		"erieRunoff_GLWBData.csv")

## Read in median values of P, E, R for each lake from L2SWBM 
## ....and transform vectors into matrices (col -- month & row -- year)
for(i in 1:9) copulaData[[i]] <- fread(paste(l2s.dir, glwbDataFiles[[i]], sep = ""), select = "Median")

for(i in 1:9) copulaData[[i]] = as.data.frame(matrix(as.matrix(copulaData[[i]]), ncol=12, byrow=TRUE))

######################################################
######################################################
#### Fitting the L2SWBM data values into a copula ####
######################################################
######################################################

## Combine all L2SWBM median values into single 70 x 108 matrix (for copula fitting)
## Organized where first 3 columns are January P, E, R for superior...
##          ....the next 3 columns are January P, E, R for MHU...
##			....and then 3 columns for January P, E, R for Erie
## [so, the first nine colums are for the month of January] 9 * 12 = 108

## Generate number sequence for this structure ([1,13,25]...[2,14,26]...[12,24...108])
## ....and generate new matrix PERyear
seq.m <- c(); for (i in 1:12) seq.m <- c(seq.m, seq(i, i+96, by=12))
PERyear <- bind_cols(copulaData)[,seq.m]
pseudo_70 <- pobs(PERyear) ## 70 year pseudo-observations: converts data points into rank normal ~ [0,1]
pseudo_50 <- pobs(PERyear[c(1:50),]) ## 50 year
psuedo_1978 <- pobs(PERyear[c(1:28),]) ##1950-1978


## Assume we will use the RVine copula for this study. Other choices, 
## such as Gaussian, Mix, Nested Archimedean were evaluated in preliminary studies.

#struct_vine_50 <- RVineStructureSelect(pseudo_50,type=0,cores=2)
#struct_vine_70 <- RVineStructureSelect(pseudo_70,type=0,cores=2)
#saveRDS(struct_vine_50,"rvine_copula_50.rds")
#saveRDS(struct_vine_70,"rvine_copula_70.rds")

## Generates & Organizes copula generated samples into 108 column matrix by component & lake
generate_predictions <- function(struct_vine,num){
  prediction <- RVineSim(num,struct_vine)
  temp <- vector(mode = "list", length = 9)
  for (m in 1:12) { for(i in 1:9) temp[[i]] <- cbind(temp[[i]], prediction[,((m-1)*9)+i]) }
  return(temp)
}
# copulaPrediction_50 <- generate_predictions(struct_vine_50,gen*horizon)
# copulaPrediction_70 <- generate_predictions(struct_vine_70,gen*horizon)
# saveRDS(copulaPrediction_50,"50_copula_prediction.rds")
# saveRDS(copulaPrediction_70,"70_copula_prediction.rds")

## Extract Predicted/Simulated data points according to components & Lake
# 1 - Lake Superior precipitation     4 - Lake Michigan & Huron precipitation     7 - Lake Erie precipitation
# 2 - Lake Superior evaporation       5 - Lake Michigan & Huron evaporation       8 - Lake Erie evaporation
# 3 - Lake Superior runoff            6 - Lake Michigan & Huron runoff            9 - Lake Erie runoff
copulaPrediction <- list(readRDS("50_copula_prediction.rds"),
                         readRDS("70_copula_prediction.rds"),
                         readRDS("70_copula_prediction.rds"),
                         generate_predictions(readRDS("rvine_copula_70.rds"),1000))

## Convert into real data given climate scenario
for(s in 1:4){
  if(s==4) gen <- 20
  for(year in 1:horizon){
    for(i in 1:108){
      if(s != 1) original <- PERyear[,i]
      else original <- PERyear[c(1:horizon),i]
      lake <- (floor((i-1)/3) %% 3) + 1 # 1 - Sup   2 - MH   3 - Eri
      month <- floor((i-1)/9) + 1  # 1 - Jan ..... 12 - Dec
      component <- ((i-1) %% 3) + 1 # 1 - precip   2 - evap   3 - run
      if(component == 1){ # precipitation (3-parameter Gamma)
        gamma_temp <- dglm(original~1,family = Gamma(link = "log"), mustart = mean(original))
        copulaPrediction[[s]][[3*(lake-1)+component]][c(((year-1)*gen+1):(year*gen)),month] <- qgamma3(copulaPrediction[[s]][[3*(lake-1)+component]][c(((year-1)*gen+1):(year*gen)),month],
                                                                                                       thres = monthly_climate_mean_perturb[[s]][[component]][lake,month] * year,
                                                                                                       shape = exp(-1*gamma_temp$dispersion.fit$coefficients[1]),
                                                                                                       scale = exp(gamma_temp$coefficients[1])/exp(-1*gamma_temp$dispersion.fit$coefficients[1]))
      } else if (component == 2){ # evaporation (normal)
        copulaPrediction[[s]][[3*(lake-1)+component]][c(((year-1)*gen+1):(year*gen)),month] <- qnorm(copulaPrediction[[s]][[3*(lake-1)+component]][c(((year-1)*gen+1):(year*gen)),month], 
                                                                                                     mean = mean(original) + monthly_climate_mean_perturb[[s]][[component]][lake,month] * year,
                                                                                                     sd = sd(original))
      } else if (component == 3){ # runoff (shifted lognormal)
        mu_h <- sum(original)/length(original)
        sig_2_h <- sum((original-mu_h)^2)/length(original)
        copulaPrediction[[s]][[3*(lake-1)+component]][c(((year-1)*gen+1):(year*gen)),month] <-qshifted_lnorm(copulaPrediction[[s]][[3*(lake-1)+component]][c(((year-1)*gen+1):(year*gen)),month],
                                                                                                             shift = monthly_climate_mean_perturb[[s]][[component]][lake,month] * year,
                                                                                                             meanlog = log(mu_h^2) - log(sig_2_h + (mu_h^2))/2, 
                                                                                                             sdlog = sqrt(log(1 + sig_2_h/(mu_h^2))))
      } 
    }
  } 
}
gen <- 999

#NBSsup = (copulaPrediction[[s]][[1]][c(((y-1)*gen+1):(y*gen)),m] - copulaPrediction[[s]][[2]][c(((y-1)*gen+1):(y*gen)),m] + copulaPrediction[[s]][[3]][c(((y-1)*gen+1):(y*gen)),m]
    

###################
### Export Water Level Files
###################
#save(total_wl, file = "wl_sim_data_CIA.Rdata") ##Please Update name for run here, do not overwrite
###################
# 
# ## Climate Perturbation factor Matrix
# for(i in 1:1){
#   mean_perturb_names <- c("Precipitation Mean Perturbation","Evaporation Mean Perturbation","Runoff Mean Perturbation")
#   pdf(paste(plot.dir,"Climate_perturb_factors.pdf", sep=""), width = 14, height = 20, paper = "special",onefile = TRUE) # create pdf
#   par(mfrow = c(1, 3)) # fit 3 plots in 1 page
#   for(i in 1:3){
#     for(s in 1:3){
#       plot(t(monthly_climate_mean_perturb[[s]][[i]]),xlab="Lake",ylab="Month",main=NA,digits=5,col=NA,key=NULL,axis.col=list(side=3),cex=1.5)
#       title(paste(mean_perturb_names[i],"Scenario",s),line = 3)
#     }
#   }
#   dev.off()
# }
# 
# ## Spatial Scatter Correlation Plots
# corr_pdf_names <- list(c("Precip_scatter_corr_Sup_vs_MH.pdf","Evap_scatter_corr_Sup_vs_MH.pdf","Run_scatter_corr_Sup_vs_MH.pdf"),
#                        c("Precip_scatter_corr_MH_vs_Eri.pdf","Evap_scatter_corr_MH_vs_Eri.pdf","Run_scatter_corr_MH_vs_Eri.pdf")
# )
# for(l in 1:2){ # lakes
#   for(c in 1:3){ # component
#     pdf(paste(plot.dir,corr_pdf_names[[l]][c], sep=""), width = 14, height = 20, paper = "special",onefile = FALSE) # create pdf
#     par(mfrow = c(4, 3))
#     for(m in 1:12){ # months
#       plot(t(copulaPrediction[[4]][[3*(l-1)+c]])[m,],t(copulaPrediction[[4]][[3*(l-1)+c+3]])[m,],xlab="",ylab='',type='p',col='blue',main=months[m],pch=16,cex.axis=1.5,cex.main=2)
#       title(ylab = paste('mm/month (',lakes[l+1],')',sep=""), line = 2.5,cex.lab=2)
#       title(xlab = paste('mm/month (',lakes[l],')',sep=""), cex.lab=2)
#       points(t(copulaData[[3*(l-1)+c]])[m,],t(copulaData[[3*(l-1)+c+3]])[m,],col='darkorange2',type='p',pch=16,cex=1.2)
#       abline(lm(t(copulaPrediction[[4]][[3*(l-1)+c+3]])[m,]~t(copulaPrediction[[4]][[3*(l-1)+c]])[m,]), col="blue",lwd=2)
#       abline(lm(t(copulaData[[3*(l-1)+c+3]])[m,]~t(copulaData[[3*(l-1)+c]])[m,]),col="darkorange2",lwd=2)
#     }
#     dev.off()
#   }
# }
# 
# ## ACF Plots
# acf_plot_names <- c("Superior Precipitation","Superior Evaporation","Superior Runoff","Michigan-Huron Precipitation","Michigan-Huron Evaporation","Michigan-Huron Runoff","Erie Precipitation","Erie Evaporation","Erie Runoff")
# pdf(paste(plot.dir, "acf.pdf", sep=""), width = 14, height = 20, paper = "special",onefile = TRUE) # create pdf
# par(mfrow = c(3, 1))
# grid_arrangement <- vector(mode="list",length=9)
# for(i in 1:9){
#   m <- acf(as.vector(t(copulaPrediction[[4]][[i]])),plot=FALSE,lag.max = 12)
#   n <- acf(as.vector(t(copulaData[[i]])),plot=FALSE,lag.max = 12)
#   acf_plot <- data.frame(lag = m$lag,acf1=m$acf,acf2=n$acf)
#   colnames(acf_plot)<-c("Lag","Prediction","Historical")
#   grid_arrangement[[i]] <- acf_plot %>% gather(key = Data, value = Correlation, -Lag) %>% ggplot(aes(x = Lag, y = Correlation, fill = Data)) + geom_col(position = "dodge") + ggtitle(acf_plot_names[i]) + scale_x_continuous(breaks = seq(0, 12, by = 1)) + theme(plot.title = element_text(size=12,hjust=0.5),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10))
# }
# grid.arrange(grid_arrangement[[1]],grid_arrangement[[2]],grid_arrangement[[3]],nrow=3)
# grid.arrange(grid_arrangement[[4]],grid_arrangement[[5]],grid_arrangement[[6]],nrow=3)
# grid.arrange(grid_arrangement[[7]],grid_arrangement[[8]],grid_arrangement[[9]],nrow=3)
# dev.off()
# 
# ## FROM: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
# GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
#                            draw_group = function(self, data, ..., draw_quantiles = NULL) {
#                              data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
#                              grp <- data[1, "group"]
#                              newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
#                              newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
#                              newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
#                              
#                              if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
#                                stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
#                                                                          1))
#                                quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
#                                aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
#                                aesthetics$alpha <- rep(1, nrow(quantiles))
#                                both <- cbind(quantiles, aesthetics)
#                                quantile_grob <- GeomPath$draw_panel(both, ...)
#                                ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
#                              }
#                              else {
#                                ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
#                              }
#                            })
# 
# geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
#                               draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
#                               show.legend = NA, inherit.aes = TRUE) {
#   layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
#         position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
#         params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
# }
# 
# ## Violin Plots
# violin_pdf_names <- c("superior_violin.pdf","michigan_huron_violin.pdf","erie_violin.pdf")
# for(l in 1:3){
#   pdf(paste(plot.dir,violin_pdf_names[l], sep=""), width = 14, height = 20, paper = "special",onefile = TRUE) # create pdf
#   par(mfrow = c(3, 1)) # fit 3 plots in 1 page
#   violin_data <- vector(mode="list",length=3)
#   grid_arrangement <- vector(mode="list",length=3)
#   for(j in 1:3){
#     violin_data[[j]] <- vector(mode = "list", length = 24)
#     for(i in 1:12){
#       violin_data[[j]][[i]] <- cbind(months[i],as.data.frame(copulaPrediction[[4]][[j+3*(l-1)]])[i],"prediction")
#       violin_data[[j]][[i+12]] <- cbind(months[i],as.data.frame(copulaData[[j+3*(l-1)]])[i],"sample")
#       names(violin_data[[j]][[i]])[2] <- "Value"
#       names(violin_data[[j]][[i]])[1] <- "Month"
#       names(violin_data[[j]][[i]])[3] <- "Type"
#       names(violin_data[[j]][[i+12]])[2] <- "Value"
#       names(violin_data[[j]][[i+12]])[1] <- "Month"
#       names(violin_data[[j]][[i+12]])[3] <- "Type"
#     }
#     violin_plot_temp <- rbind(violin_data[[j]][[1]],violin_data[[j]][[2]])
#     for(k in 3:24) violin_plot_temp <- rbind(violin_plot_temp,violin_data[[j]][[k]])
#     grid_arrangement[[j]] <- ggplot(as.data.frame(violin_plot_temp),aes(x=Month,y=Value,fill=Type)) +  geom_split_violin(trim=FALSE) + xlim(months) + stat_summary(fun=median, geom="point", size=2, color="red")
#   }
#   grid.arrange(grid_arrangement[[1]],grid_arrangement[[2]],grid_arrangement[[3]],nrow=3)
#   dev.off()
# }
# 
# 
# 
# 
# # If Statistics.txt already exists, please delete file before running the following line of code.
# sink(file = "Statistics.txt")
# for(i in 1:1){
#   cat(noquote(paste("Observed Annual Average NBS (1950-2019)\n")))
#   cat(noquote(paste("\tLake Superior:",mean(unlist(copulaData[[1]]-copulaData[[2]]+copulaData[[3]])),"\n")))
#   cat(noquote(paste("\tLake Michigan-Huron:",mean(unlist(copulaData[[4]]-copulaData[[5]]+copulaData[[6]])),"\n")))
#   cat(noquote(paste("\tLake Erie:",mean(unlist(copulaData[[7]]-copulaData[[8]]+copulaData[[9]])),"\n")))
#   
#   cat(noquote(paste("Simulated Average Annual NBS Scenario 1 (temp:SC1-RM)\n")))
#   cat(noquote(paste("\tLake Superior:",mean(unlist(copulaPrediction[[1]][[1]]-copulaPrediction[[1]][[2]]+copulaPrediction[[1]][[3]])),"\n")))
#   cat(noquote(paste("\tLake Michigan-Huron:",mean(unlist(copulaPrediction[[1]][[4]]-copulaPrediction[[1]][[5]]+copulaPrediction[[1]][[6]])),"\n")))
#   cat(noquote(paste("\tLake Erie:",mean(unlist(copulaPrediction[[1]][[7]]-copulaPrediction[[1]][[8]]+copulaPrediction[[1]][[9]])),"\n")))
#   
#   cat(noquote(paste("Simulated Average Annual NBS Scenario 2 (temp:SC1-PM)\n")))
#   cat(noquote(paste("\tLake Superior:",mean(unlist(copulaPrediction[[2]][[1]]-copulaPrediction[[2]][[2]]+copulaPrediction[[2]][[3]])),"\n")))
#   cat(noquote(paste("\tLake Michigan-Huron:",mean(unlist(copulaPrediction[[2]][[4]]-copulaPrediction[[2]][[5]]+copulaPrediction[[2]][[6]])),"\n")))
#   cat(noquote(paste("\tLake Erie:",mean(unlist(copulaPrediction[[2]][[7]]-copulaPrediction[[2]][[8]]+copulaPrediction[[2]][[9]])),"\n")))
#   
#   cat(noquote(paste("Simulated Average Annual NBS Scenario 3 (temp:SC1-asdsaM)\n")))
#   cat(noquote(paste("\tLake Superior:",mean(unlist(copulaPrediction[[3]][[1]]-copulaPrediction[[3]][[2]]+copulaPrediction[[3]][[3]])),"\n")))
#   cat(noquote(paste("\tLake Michigan-Huron:",mean(unlist(copulaPrediction[[3]][[4]]-copulaPrediction[[3]][[5]]+copulaPrediction[[3]][[6]])),"\n")))
#   cat(noquote(paste("\tLake Erie:",mean(unlist(copulaPrediction[[3]][[7]]-copulaPrediction[[3]][[8]]+copulaPrediction[[3]][[9]])),"\n")))
# 
# }
# sink(file = NULL)

