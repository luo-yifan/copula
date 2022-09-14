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
library(BRugs)
library(R2WinBUGS)

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
horizon <- 50 			# copula simulation time horizon in years
gen 	<- 999 			# number of copula predictions per variable; BUGS gives vector length 999

## Monthly Climate Perturbations for each Component
#               |  Jan ~ Dec  
# Sup (factor)  |_____________
# MH  (factor)  |_____________
# Eri (factor)  |_____________
monthly_climate_mean_perturb <- rep(list(rep(list(matrix(nrow = 3,ncol = 12,dimnames = list(c("Superior","Michigan-Huron","Erie"),months))),3)),4)

## 1 - no climate change  2 - baseline(self-generated)  3 - high-warming(via papers)  4 - Copula Validation Plots
climate_data <- read.csv(file = paste(main.dir,"/Config.csv", sep = ""),header=F)
for(s in 1:4){ ## Scenario
  for(i in 1:3){ ## lake
    for(j in 1:3){ ## component
      for(m in 1:12){ ## month
        if(s==4){
          monthly_climate_mean_perturb[[s]][[j]][i,m] <- climate_data[36*(i-1)+12*(j-1)+m,1]
        }else{
          monthly_climate_mean_perturb[[s]][[j]][i,m] <- climate_data[36*(i-1)+12*(j-1)+m,s]
        }
      }
    }
  }
}

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
      else original <- PERyear[c(1:50),i]
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

####################################################
####################################################
### Simulate outflow and water levels
####################################################
####################################################

## We want to utilize historical water level data from federal agencies
##    (e.g. coordinating committee)...but we want to use lake-to-lake outflow
##    data from the L2SWBM; these are values that we know close the water balance

## Names of files to be read from
outflowFiles <- c( "SUP_BOM_MM.csv",
                   "superiorOutflow_GLWLRuns_MHDivCons_2M.csv",   
                   "MHG_BOM_MM.csv",
                   "miHuronOutflow_GLWLRuns_MHDivCons_2M.csv",
                   "ERI_BOM_MM.csv",
                   "erieOutflow_GLWLRuns_MHDivCons_2M.csv")

## Emply lists for historical water level data and l2s outflows 
# 1 - Lake Superior water level   3 - Lake Michigan & Huron water level   5 - Lake Erie water level
# 2 - Lake Superior outflow       4 - Lake Michigan & Huron outflow       6 - Lake Erie outflow
outflowData 	<- vector(mode = "list", length = 6)
logOutflowData 	<- vector(mode = "list", length = 6)

for(i in 1:6){
  if(i %% 2 == 0){
    outflowData[[i]] <- fread(paste(data.dir, outflowFiles[[i]], sep = ""), select = "Median")
    outflowData[[i]] <- data.frame(matrix(unlist(outflowData[[i]]),ncol=12,byrow = T)[c(2:51),]) 
  } else{
    outflowData[[i]] <- fread(paste(data.dir, outflowFiles[[i]], sep = ""), select = "Beginning of Month")
    outflowData[[i]] <- data.frame(matrix(unlist(outflowData[[i]]),ncol=12,byrow = T)[c(52:101),])
  }
  logOutflowData[[i]] <- log(outflowData[[i]])
}

## BUGS MODEL
intercept_consts 	<- list(length=3) # seasonal intercepts
stage_consts 		<- list(length=3) # constant stage
vol_consts 			<- list(length=3) # sigma
fall_consts 		<- list(length=1) # constant fall
sill_adjustments 	<- c(181.43,166.55,169.94)
bugs_file_path 		<- paste(getwd(),"/sfd_bug.bug",sep="")

modelCheck(bugs_file_path)

for(i in 1:3){
  # log water level with sill elevation adjustments
  out 	<- as.vector(unlist(t(log(outflowData[[(2*i)]])))) # log outflow
  wl 	<- as.vector(unlist(t(log(outflowData[[(2*i-1)]] - sill_adjustments[i])))) 
  fall 	<- as.vector(unlist(t(log(outflowData[[3]] - outflowData[[5]])))) # log of fall between MH & ERI
  J 	<- length(out)
  m 	<- rep(1:12,length.out=J)
  
  parameters <- c("a","c","sigma","d")
  
  if(i==2){isMH <- 1}else{isMH <- -1}
  data 	<- c("J","m","out","wl","fall","isMH")
  
  sfd.bug <- openbugs(data, inits=NULL, parameters, bugs_file_path, n.chains = 3, n.iter = 20000)
  intercept_consts[[i]] <- sfd.bug$sims.list$a
  stage_consts[[i]] <- sfd.bug$sims.list$c
  vol_consts[[i]] <- sfd.bug$sims.list$sigma
  if(isMH == 1) fall_consts[[1]] <- sfd.bug$sims.list$d
}  

# Plug in diversions 
diversion <- vector(mode = "list",length = 3) 	# IN & OUT in m^3/s

#Sup_IN <- 159;
Sup_IN <- fread(paste(l2s.dir, "superiorDiversion_GLWBData.csv", sep = ""), select = "Median")
Sup_OUT <- 0
diversion[[1]] <- (Sup_IN[14:613,] - Sup_OUT)  	#*(1000^3)/(82097*10^12) # in mm*month/sec

MH_IN <- 0
MH_OUT <- fread(paste(l2s.dir, "miHuronDiversion_GLWBData.csv", sep = ""), select = "Median")
diversion[[2]] <- (MH_IN - MH_OUT[14:613,])		#*(1000^3)/(117318*10^12) # in mm*month/sec

#### Not sure where this came from; I think it's the Welland Canal....water outflow
####   above and beyond what is measured in the Niagara River proper
Eri_IN <- 0; 
Eri_OUT <- fread(paste(l2s.dir, "erieDiversion_GLWBData.csv", sep = ""), select = "Median")
diversion[[3]] <- (Eri_IN - Eri_OUT[14:613,])	#*(1000^3)/(25655*10^12) # in mm*month/sec

## Generate Water level predictions
total_wl <- vector(mode="list",length=3)
total_out <- vector(mode="list",length=3)
for(s in 1:3){
  # year of starting water level = 1900 + addYear (50 ~ 100)
  if(s != 1) addYear <- 50
  else addYear <- 1
  current_level <- list(rep(as.numeric(outflowData[[1]][addYear, 1]),gen),
                        rep(as.numeric(outflowData[[3]][addYear, 1]),genx),
                        rep(as.numeric(outflowData[[5]][addYear, 1]),gen))
  water_level <- vector(mode="list",length=3)
  outflow <- list(length=3)
  outflow_calculation <- list(c(),c(),c()) # holds outflow data from previous predictions
  outflow_plot <- list(c(),c(),c()) # outflow data for plotting
  for(i in 1:3) water_level[[i]] <- rep(current_level[[i]][1],gen) # starting water level
  for(i in 1:(12*horizon)){
    m <- if(i %% 12 == 0) 12 else (i %% 12) # current month
    y <- floor((i-1)/12) + 1 # current year
    
    # Lake Superior
    outflow[[1]] <- exp(intercept_consts[[1]][c(1:gen),m] + stage_consts[[1]][c(1:gen)]*log(current_level[[1]]-sill_adjustments[1]) + rnorm(gen,0,vol_consts[[1]][c(1:gen)]))
    outflow_calculation[[1]] <- cbind(outflow_calculation[[1]],outflow[[1]])
    outflow_plot[[1]] <- cbind(outflow_plot[[1]],outflow[[1]])
    # Convert cms --> mm per month
    if (m %in% c(1,3,5,7,8,10,12)){
      outflow[[1]] <- (outflow[[1]]/S_SF)*(SPM31)/(1000) # months with 31 days
      temp_div_sup <- (as.ts(diversion[[1]])/S_SF)*(SPM31)/(1000)
    } else if (m %in% c(4,6,9,11)){
      outflow[[1]] <- (outflow[[1]]/S_SF)*(SPM30)/(1000) # months with 30 days
      temp_div_sup <- (as.ts(diversion[[1]])/S_SF)*(SPM30)/(1000)
    } else{
      outflow[[1]] <- (outflow[[1]]/S_SF)*(SPMF)/(1000) # February
      temp_div_sup <- (as.ts(diversion[[1]])/S_SF)*(SPMF)/(1000)
    }
    # add water level change to update current simulated water level
    current_level[[1]] <- current_level[[1]] + temp_div_sup[i]/1000 + (copulaPrediction[[s]][[1]][c(((y-1)*gen+1):(y*gen)),m] - copulaPrediction[[s]][[2]][c(((y-1)*gen+1):(y*gen)),m] + copulaPrediction[[s]][[3]][c(((y-1)*gen+1):(y*gen)),m] - outflow[[1]])/1000
    water_level[[1]] <- cbind(water_level[[1]],current_level[[1]])
    
    # Lake Michigan-Huron
    outflow[[2]] <- exp(intercept_consts[[2]][c(1:gen),m] + stage_consts[[2]][c(1:gen)]*log(current_level[[2]]-sill_adjustments[2]) + fall_consts[[1]][c(1:gen)]*log(current_level[[2]]-current_level[[3]]) + rnorm(gen,0,vol_consts[[2]][c(1:gen)]))
    outflow_calculation[[2]] <- cbind(outflow_calculation[[2]],outflow[[2]])
    outflow_plot[[2]] <- cbind(outflow_plot[[2]],outflow[[2]])
    # Convert cubic cm per sec of flow --> water level per month
    if (m %in% c(1,3,5,7,8,10,12)){
      outflow[[2]] <- (outflow[[2]]/MH_SF)*(SPM31)/1000 # months with 31 days
      outflow_calculation[[1]][,ncol(outflow_calculation[[1]])] <- (outflow_calculation[[1]][,ncol(outflow_calculation[[1]])]/MH_SF)*SPM31/1000
      temp_div_mh <- (as.ts(diversion[[2]])/MH_SF)*SPM31/1000
    } else if (m %in% c(4,6,9,11)){
      outflow[[2]] <- (outflow[[2]]/MH_SF)*(SPM30)/1000 # months with 30 days
      outflow_calculation[[1]][,ncol(outflow_calculation[[1]])] <- (outflow_calculation[[1]][,ncol(outflow_calculation[[1]])]/MH_SF)*SPM30/1000
      temp_div_mh <- (as.ts(diversion[[2]])/MH_SF)*SPM30/1000
    } else{
      outflow[[2]] <- (outflow[[2]]/MH_SF)*(SPMF)/1000 # February
      outflow_calculation[[1]][,ncol(outflow_calculation[[1]])] <- (outflow_calculation[[1]][,ncol(outflow_calculation[[1]])]/MH_SF)*SPMF/1000
      temp_div_mh <- (as.ts(diversion[[2]])/MH_SF)*SPMF/1000
    }
    # add water level change to update current simulated water level
    current_level[[2]] <- current_level[[2]] + temp_div_mh[i]/1000 + (outflow_calculation[[1]][,ncol(outflow_calculation[[1]])]*0.8 + copulaPrediction[[s]][[4]][c(((y-1)*gen+1):(y*gen)),m] - copulaPrediction[[s]][[5]][c(((y-1)*gen+1):(y*gen)),m] + copulaPrediction[[s]][[6]][c(((y-1)*gen+1):(y*gen)),m] - outflow[[2]])/1000
    water_level[[2]] <- cbind(water_level[[2]],current_level[[2]])
    
    # Lake Erie
    outflow[[3]] <- exp(intercept_consts[[3]][c(1:gen),m] + stage_consts[[3]][c(1:gen)]*log(current_level[[3]]-sill_adjustments[3]) + rnorm(gen,0,vol_consts[[3]][c(1:gen)]))
    outflow_calculation[[3]] <- cbind(outflow_calculation[[3]],outflow[[3]])
    outflow_plot[[3]] <- cbind(outflow_plot[[3]],outflow[[3]])
    # Convert cubic cm per sec of flow --> water level per month
    if (m %in% c(1,3,5,7,8,10,12)){
      outflow[[3]] <- (outflow[[3]]/E_SF)*(SPM31)/1000 # months with 31 days
      outflow_calculation[[2]][,ncol(outflow_calculation[[2]])] <- (outflow_calculation[[2]][,ncol(outflow_calculation[[2]])]/E_SF)*SPM31/1000
      temp_div_e <- (as.ts(diversion[[3]])/E_SF)*SPM31/1000
    } else if (m %in% c(4,6,9,11)){
      outflow[[3]] <- (outflow[[3]]/E_SF)*(SPM30)/1000 # months with 30 days
      outflow_calculation[[2]][,ncol(outflow_calculation[[2]])] <- (outflow_calculation[[2]][,ncol(outflow_calculation[[2]])]/E_SF)*SPM30/1000
      temp_div_e <- (as.ts(diversion[[3]])/E_SF)*SPM30/1000
    } else{
      outflow[[3]] <- (outflow[[3]]/E_SF)*(SPMF)/1000 # February
      outflow_calculation[[2]][,ncol(outflow_calculation[[2]])] <- (outflow_calculation[[2]][,ncol(outflow_calculation[[2]])]/E_SF)*SPMF/1000
      temp_div_e <- (as.ts(diversion[[3]])/E_SF)*SPMF/1000
    }
    # add water level change to update current simulated water level
    current_level[[3]] <- current_level[[3]] + temp_div_e[i]/1000 + (outflow_calculation[[2]][,ncol(outflow_calculation[[2]])] + copulaPrediction[[s]][[7]][c(((y-1)*gen+1):(y*gen)),m] - copulaPrediction[[s]][[8]][c(((y-1)*gen+1):(y*gen)),m] + copulaPrediction[[s]][[9]][c(((y-1)*gen+1):(y*gen)),m] - outflow[[3]])/1000
    water_level[[3]] <- cbind(water_level[[3]],current_level[[3]])
  }
  total_out[[s]] <- outflow_plot
  total_wl[[s]] <- water_level
}

## Water level Plots
Wl_plot_name <- c("scenario1_wl.pdf","scenario2_wl.pdf","scenario3_wl.pdf")
for(s in 1:3){
  # water level plot
  pdf(paste(plot.dir, Wl_plot_name[[s]], sep=""), width = 14, height = 20, paper = "special",onefile = FALSE) # create pdf
  par(mfrow = c(3, 1))
  
  plot(total_wl[[s]][[1]][1,],xlab="Month",ylab="Water level (m)",ylim=c(182,185),type='l',col='blue',main='Superior Water levels')
  for(i in 2:gen) points(total_wl[[s]][[1]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
  if(s==1) points(as.vector(t(outflowData[[1]][c(1:50),])),type='l',lwd=3)
  
  plot(total_wl[[s]][[2]][1,],xlab="Month",ylab="Water level (m)",ylim=c(175,178),type='l',col='blue',main='Michigan-Huron Water levels')
  for(i in 2:gen) points(total_wl[[s]][[2]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
  if(s==1) points(as.vector(t(outflowData[[3]][c(1:50),])),type='l',lwd=3)
  
  plot(total_wl[[s]][[3]][1,],xlab="Month",ylab="Water level (m)",ylim=c(173,176),type='l',col='blue',main='Erie Water levels')
  for(i in 2:gen) points(total_wl[[s]][[3]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
  if(s==1) points(as.vector(t(outflowData[[5]][c(1:50),])),type='l',lwd=3)
  
  dev.off()
}

## Outflow Plots 
out_plot_name <- c("scenario1_out.pdf","scenario2_out.pdf","scenario3_out.pdf")
for(s in 1:3){
  pdf(paste(plot.dir, out_plot_name[s], sep=""), width = 14, height = 18, paper = "special",onefile = TRUE) # create pdf
  par(mfrow = c(3,1))
  
  par(mar=c(0,5,5,5))
  plot(total_out[[s]][[1]][1,],xlab="",ylab="",ylim=c(0,8000),type='l',col='red',main='',cex.main=2,cex.axis=1.5,cex.lab=2,xaxt='n')
  title(main="Superior Outflow",line = -2,cex.main=2.5)
  for(i in 2:gen) points(total_out[[s]][[1]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
  if(s==1) points(as.vector(t(outflowData[[2]][c(1:50),])),type='l',lwd=3)
  
  par(mar=c(0,5,0,5))
  plot(total_out[[s]][[2]][1,],xlab="",ylab="",ylim=c(2000,10000),type='l',col='red',main='',cex.main=2,cex.axis=1.5,cex.lab=2,xaxt='n')
  title(ylab=expression(Outflow ~ (cm^"3")),line=2,cex.lab=2)
  title(main="Michigan-Huron Outflow",line = -2,cex.main=2.5)
  for(i in 2:gen) points(total_out[[s]][[2]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
  if(s == 1) points(as.vector(t(outflowData[[4]][c(1:50),])),type='l',lwd=3)
  
  par(mar=c(5,5,0,5))
  plot(total_out[[s]][[3]][1,],xlab="Month",ylab="",ylim=c(2000,10000),type='l',col='red',main='',cex.main=2,cex.axis=1.5,cex.lab=2)
  title(main="Erie Outflow",line = -2,cex.main=2.5)
  for(i in 2:gen) points(total_out[[s]][[3]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
  if(s==1 ) points(as.vector(t(outflowData[[6]][c(1:50),])),type='l',lwd=3)
  
  dev.off()
} 

## Long term average
# 1 - Superior  2 - Michigan-Huron   3 - Erie
long_wl_means <- list(
  list(rep(0,gen),rep(0,gen),rep(0,gen)),
  list(rep(0,gen),rep(0,gen),rep(0,gen)),
  list(rep(0,gen),rep(0,gen),rep(0,gen))
)
for(s in 1:3){
  for(i in 1:gen){
    long_wl_means[[s]][[1]][i] <- mean(total_wl[[s]][[1]][i,])
    long_wl_means[[s]][[2]][i] <- mean(total_wl[[s]][[2]][i,])
    long_wl_means[[s]][[3]][i] <- mean(total_wl[[s]][[3]][i,])
  }
}
long_wl_means_mean <- list(
  c(mean(long_wl_means[[1]][[1]]),mean(long_wl_means[[1]][[2]]),mean(long_wl_means[[1]][[3]])),
  c(mean(long_wl_means[[2]][[1]]),mean(long_wl_means[[2]][[2]]),mean(long_wl_means[[2]][[3]])),
  c(mean(long_wl_means[[3]][[1]]),mean(long_wl_means[[3]][[2]]),mean(long_wl_means[[3]][[3]]))
)
long_term_avg_plot_names <- c("Long_term_average_Page_1.jpg","Long_term_average_Page_2.jpg","Long_term_average_Page_3.jpg")
for(i in 1:3){
  jpeg(paste(plot.dir,long_term_avg_plot_names[i],sep=""),width=20,height=7,units = "in", res = 72*4)
  mean_plot <- data.frame()
  for(s in 1:3) mean_plot <- rbind(mean_plot,data.frame(Mean = long_wl_means[[s]][[i]],Scenario = paste("Scenario",s)))
  grid.arrange(ggplot(mean_plot,aes(fill=Scenario)) + geom_density(aes(x=Mean,y=100*..count../sum(..count..)),alpha=0.5) + scale_y_continuous(labels=scales::percent_format(accuracy=0.1)) + ggtitle(lakes[i]) + xlab("Mean") + ylab("Density (%)") + theme(text=element_text(size=23),axis.title.x = element_text(size=23),axis.title.y = element_text(size=23),legend.text=element_text(size=23),axis.text = element_text(size=18)))
  dev.off()
}

## Extrema Frequency Density
annual_extrema_plot_names <- c("Annual Extrema_Page_1.jpg","Annual Extrema_Page_2.jpg","Annual Extrema_Page_3.jpg")
extrema_temp_plot <- vector(mode="list",length=3)
total_annual_extrema <- vector(mode="list",length=3)
for(s in 1:3){
  annual_extrema <- list(list(rep(0,12),rep(0,12)),
                         list(rep(0,12),rep(0,12)),
                         list(rep(0,12),rep(0,12)))
  for(l in 1:3){
    temp_wl <- matrix(total_wl[[s]][[l]][1,][2:601],ncol=12,byrow=TRUE)
    for(i in 2:gen) temp_wl<- rbind(temp_wl,matrix(total_wl[[s]][[l]][i,][2:601],ncol=12,byrow=TRUE))
    for(i in 1:(gen*horizon)){
      for(j in which(temp_wl[i,]==max(temp_wl[i,]))) annual_extrema[[l]][[1]][[j]] <- annual_extrema[[l]][[1]][[j]] + 1
      for(j in which(temp_wl[i,]==min(temp_wl[i,]))) annual_extrema[[l]][[2]][[j]] <- annual_extrema[[l]][[2]][[j]] + 1
    }
    extrema_temp_plot[[l]] <- rbind(extrema_temp_plot[[l]],data.frame(high=as.vector(rep(1:12,annual_extrema[[l]][[1]])),low=as.vector(rep(1:12,annual_extrema[[l]][[2]])),Scenario=paste("Scenario",s)))
  }
  total_annual_extrema[[s]] <- annual_extrema
}
## Extrema for Observed
observed_extrema <- list(list(rep(0,12),rep(0,12)),
                         list(rep(0,12),rep(0,12)),
                         list(rep(0,12),rep(0,12)))
for(l in 1:3){
  temp_wl <- matrix(unlist(outflowData[[2*l-1]]),ncol=12,byrow=F)
  for(i in 1:50){
    for(j in which(temp_wl[i,]==max(temp_wl[i,]))) observed_extrema[[l]][[1]][[j]] <- observed_extrema[[l]][[1]][[j]] + 1
    for(j in which(temp_wl[i,]==min(temp_wl[i,]))) observed_extrema[[l]][[2]][[j]] <- observed_extrema[[l]][[2]][[j]] + 1
  }
}
for(i in 1:3){
  jpeg(paste(plot.dir,annual_extrema_plot_names[i],sep=""),width=20,height=7,units = "in", res = 72*4)
  grid.arrange(ggplot(extrema_temp_plot[[i]],aes(fill=Scenario)) + geom_histogram(aes(x=high,y=(..count..)/sum(..count..)),binwidth = 0.5,position="dodge") + geom_histogram(aes(x=low,y=-(..count..)/sum(..count..)),binwidth = 0.5,position="dodge") + scale_y_continuous(labels=scales::percent) + xlab("Month") + ylab("Density (%)") + scale_x_continuous(breaks=1:12,labels=months) + geom_hline(yintercept = 0, color = "black") + ggtitle(lakes[i])+ theme(text=element_text(size=20),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.text=element_text(size=20)))
  dev.off()
}

## Water Level Density
grid_arrangement <- vector(mode="list",length=6)
plot_temp <- list(length=2)
for(l in 1:3){
  for(s in 2:3){
    Sim_diff <- as.vector(total_wl[[s]][[l]]) - mean(total_wl[[s]][[l]])
    Comp_diff <- as.vector(total_wl[[1]][[l]]) - mean(total_wl[[1]][[l]])
    plot_temp[[1]] <- data.frame(a=Comp_diff,name="Scenario 1")
    plot_temp[[2]] <- data.frame(b=Sim_diff,name=paste("Scenario",s))
    grid_arrangement[[(2*(l-1)+s-1)]] <- ggplot() + geom_density(data=plot_temp[[1]],aes(x=a,y=100*..count../sum(..count..),fill=name),alpha=0.5) + geom_density(data=plot_temp[[2]],aes(x=b,y=100*..count../sum(..count..),fill=name),alpha=0.5) + theme(axis.title.y=element_blank(),axis.title.x = element_blank(),axis.title = element_blank()) + scale_y_continuous(labels=scales::percent_format(accuracy=0.1),limits = c(0,1)) + guides(fill=guide_legend(title="Scenarios"),) + theme(text=element_text(size=20),legend.text = element_text(size=20)) + xlim(-1.0,1.0)
  }
}
jpeg(paste(plot.dir,"WLDensity_Page_1.jpg",sep=""),width=20,height=7,units = "in", res = 72*4)
grid.arrange(arrangeGrob(grid_arrangement[[1]]+theme(legend.position="none",plot.title=element_text(margin=margin(t=40,b=-30)))+ggtitle("Lake Superior"),grid_arrangement[[3]]+theme(legend.position="none",plot.title=element_text(margin=margin(t=40,b=-30)))+ggtitle("Lake Michigan-Huron"),grid_arrangement[[5]]+theme(legend.position="none",plot.title=element_text(margin=margin(t=40,b=-30)))+ggtitle("Lake Erie"),ncol=4,gtable_filter(ggplotGrob(grid_arrangement[[1]]), "guide-box"),widths=c(2.3, 2.3, 2.3, 0.8)),left=textGrob("Density %", gp=gpar(fontsize=20,font=8),rot=90),bottom=textGrob("Difference from mean (m)",gp=gpar(fontsize=20,font=8)))
dev.off()
jpeg(paste(plot.dir,"WLDensity_Page_2.jpg",sep=""),width=20,height=7,units = "in", res = 72*4)
grid.arrange(arrangeGrob(grid_arrangement[[2]]+theme(legend.position="none",plot.title=element_text(margin=margin(t=40,b=-30)))+ggtitle("Lake Superior"),grid_arrangement[[4]]+theme(legend.position="none",plot.title=element_text(margin=margin(t=40,b=-30)))+ggtitle("Lake Michigan-Huron"),grid_arrangement[[6]]+theme(legend.position="none",plot.title=element_text(margin=margin(t=40,b=-30)))+ggtitle("Lake Erie"),ncol=4,gtable_filter(ggplotGrob(grid_arrangement[[4]]), "guide-box"),widths=c(2.3, 2.3, 2.3, 0.8)),left=textGrob("Density %", gp=gpar(fontsize=20,font=8),rot=90),bottom=textGrob("Difference from mean (m)",gp=gpar(fontsize=20,font=8)))
dev.off()

## Climate Perturbation factor Matrix
for(i in 1:1){
  mean_perturb_names <- c("Precipitation Mean Perturbation","Evaporation Mean Perturbation","Runoff Mean Perturbation")
  pdf(paste(plot.dir,"Climate_perturb_factors.pdf", sep=""), width = 14, height = 20, paper = "special",onefile = TRUE) # create pdf
  par(mfrow = c(1, 3)) # fit 3 plots in 1 page
  for(i in 1:3){
    for(s in 1:3){
      plot(t(monthly_climate_mean_perturb[[s]][[i]]),xlab="Lake",ylab="Month",main=NA,digits=5,col=NA,key=NULL,axis.col=list(side=3),cex=1.5)
      title(paste(mean_perturb_names[i],"Scenario",s),line = 3)
    }
  }
  dev.off()
}

## Spatial Scatter Correlation Plots
corr_pdf_names <- list(c("Precip_scatter_corr_Sup_vs_MH.pdf","Evap_scatter_corr_Sup_vs_MH.pdf","Run_scatter_corr_Sup_vs_MH.pdf"),
                       c("Precip_scatter_corr_MH_vs_Eri.pdf","Evap_scatter_corr_MH_vs_Eri.pdf","Run_scatter_corr_MH_vs_Eri.pdf")
                      )
for(l in 1:2){ # lakes
  for(c in 1:3){ # component
    pdf(paste(plot.dir,corr_pdf_names[[l]][c], sep=""), width = 14, height = 20, paper = "special",onefile = FALSE) # create pdf
    par(mfrow = c(4, 3))
    for(m in 1:12){ # months
      plot(t(copulaPrediction[[4]][[3*(l-1)+c]])[m,],t(copulaPrediction[[4]][[3*(l-1)+c+3]])[m,],xlab="",ylab='',type='p',col='blue',main=months[m],pch=16,cex.axis=1.5,cex.main=2)
      title(ylab = paste('mm/month (',lakes[l+1],')',sep=""), line = 2.5,cex.lab=2)
      title(xlab = paste('mm/month (',lakes[l],')',sep=""), cex.lab=2)
      points(t(copulaData[[3*(l-1)+c]])[m,],t(copulaData[[3*(l-1)+c+3]])[m,],col='darkorange2',type='p',pch=16,cex=1.2)
      abline(lm(t(copulaPrediction[[4]][[3*(l-1)+c+3]])[m,]~t(copulaPrediction[[4]][[3*(l-1)+c]])[m,]), col="blue",lwd=2)
      abline(lm(t(copulaData[[3*(l-1)+c+3]])[m,]~t(copulaData[[3*(l-1)+c]])[m,]),col="darkorange2",lwd=2)
    }
    dev.off()
  }
}

## ACF Plots
acf_plot_names <- c("Superior Precipitation","Superior Evaporation","Superior Runoff","Michigan-Huron Precipitation","Michigan-Huron Evaporation","Michigan-Huron Runoff","Erie Precipitation","Erie Evaporation","Erie Runoff")
pdf(paste(plot.dir, "acf.pdf", sep=""), width = 14, height = 20, paper = "special",onefile = TRUE) # create pdf
par(mfrow = c(3, 1))
grid_arrangement <- vector(mode="list",length=9)
for(i in 1:9){
  m <- acf(as.vector(t(copulaPrediction[[4]][[i]])),plot=FALSE,lag.max = 12)
  n <- acf(as.vector(t(copulaData[[i]])),plot=FALSE,lag.max = 12)
  acf_plot <- data.frame(lag = m$lag,acf1=m$acf,acf2=n$acf)
  colnames(acf_plot)<-c("Lag","Prediction","Historical")
  grid_arrangement[[i]] <- acf_plot %>% gather(key = Data, value = Correlation, -Lag) %>% ggplot(aes(x = Lag, y = Correlation, fill = Data)) + geom_col(position = "dodge") + ggtitle(acf_plot_names[i]) + scale_x_continuous(breaks = seq(0, 12, by = 1)) + theme(plot.title = element_text(size=12,hjust=0.5),axis.text.y = element_text(size=10),axis.text.x = element_text(size=10))
}
grid.arrange(grid_arrangement[[1]],grid_arrangement[[2]],grid_arrangement[[3]],nrow=3)
grid.arrange(grid_arrangement[[4]],grid_arrangement[[5]],grid_arrangement[[6]],nrow=3)
grid.arrange(grid_arrangement[[7]],grid_arrangement[[8]],grid_arrangement[[9]],nrow=3)
dev.off()

## FROM: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

## Violin Plots
violin_pdf_names <- c("superior_violin.pdf","michigan_huron_violin.pdf","erie_violin.pdf")
for(l in 1:3){
  pdf(paste(plot.dir,violin_pdf_names[l], sep=""), width = 14, height = 20, paper = "special",onefile = TRUE) # create pdf
  par(mfrow = c(3, 1)) # fit 3 plots in 1 page
  violin_data <- vector(mode="list",length=3)
  grid_arrangement <- vector(mode="list",length=3)
  for(j in 1:3){
    violin_data[[j]] <- vector(mode = "list", length = 24)
    for(i in 1:12){
      violin_data[[j]][[i]] <- cbind(months[i],as.data.frame(copulaPrediction[[4]][[j+3*(l-1)]])[i],"prediction")
      violin_data[[j]][[i+12]] <- cbind(months[i],as.data.frame(copulaData[[j+3*(l-1)]])[i],"sample")
      names(violin_data[[j]][[i]])[2] <- "Value"
      names(violin_data[[j]][[i]])[1] <- "Month"
      names(violin_data[[j]][[i]])[3] <- "Type"
      names(violin_data[[j]][[i+12]])[2] <- "Value"
      names(violin_data[[j]][[i+12]])[1] <- "Month"
      names(violin_data[[j]][[i+12]])[3] <- "Type"
    }
    violin_plot_temp <- rbind(violin_data[[j]][[1]],violin_data[[j]][[2]])
    for(k in 3:24) violin_plot_temp <- rbind(violin_plot_temp,violin_data[[j]][[k]])
    grid_arrangement[[j]] <- ggplot(as.data.frame(violin_plot_temp),aes(x=Month,y=Value,fill=Type)) +  geom_split_violin(trim=FALSE) + xlim(months) + stat_summary(fun=median, geom="point", size=2, color="red")
  }
  grid.arrange(grid_arrangement[[1]],grid_arrangement[[2]],grid_arrangement[[3]],nrow=3)
  dev.off()
}

## Total Statistics
Weighted_extrema <- list(
  list(c(),c(),c()),
  list(c(),c(),c()),
  list(c(),c(),c()),
  list(c(),c(),c())
)
for(l in 1:3){
  Weighted_extrema[[4]][[l]][1] <- sum(seq(1,12)*observed_extrema[[l]][[1]]/sum(observed_extrema[[l]][[1]]))
  Weighted_extrema[[4]][[l]][2] <- sum(seq(1,12)*observed_extrema[[l]][[2]]/sum(observed_extrema[[l]][[2]]))
}
for(s in 1:3){
  for(l in 1:3){
    Weighted_extrema[[s]][[l]][1] <- sum(seq(1,12)*total_annual_extrema[[s]][[l]][[1]]/sum(total_annual_extrema[[s]][[l]][[1]]))
    Weighted_extrema[[s]][[l]][2] <- sum(seq(1,12)*total_annual_extrema[[s]][[l]][[2]]/sum(total_annual_extrema[[s]][[l]][[2]]))
  }
}
observed_0_thresh <- list(length=3)
Comp_0_thresh <- list(length=3)
percent_past_thresh <- list(c(),c(),c(),c())
for(l in 1:3){
  Observed_diff <- as.vector(unlist(outflowData[[2*l-1]])) - mean(as.vector(unlist(outflowData[[2*l-1]])))
  observed_0_thresh[[l]] <- c(mean(Observed_diff)-2*sd(Observed_diff),mean(Observed_diff)+2*sd(Observed_diff))
  Comp_diff <- as.vector(total_wl[[1]][[l]]) - mean(total_wl[[1]][[l]])
  Comp_0_thresh[[l]] <- c(mean(Comp_diff)-2*sd(Comp_diff),mean(Comp_diff)+2*sd(Comp_diff))
  # observed vs sc1
  percent_past_thresh[[1]][l] <- 100*(length(Comp_diff[Comp_diff>observed_0_thresh[[l]][2]])/length(Comp_diff) + length(Comp_diff[Comp_diff<observed_0_thresh[[l]][1]])/length(Comp_diff))
  # sc1 vs sc1
  percent_past_thresh[[2]][l] <- 100*(length(Comp_diff[Comp_diff>Comp_0_thresh[[l]][2]])/length(Comp_diff) + length(Comp_diff[Comp_diff<Comp_0_thresh[[l]][1]])/length(Comp_diff))
  Sim_diff_2 <- as.vector(total_wl[[2]][[l]]) - mean(total_wl[[2]][[l]])
  # sc 1 vs sc 2
  percent_past_thresh[[3]][l] <- 100*(length(Sim_diff_2[Sim_diff_2>Comp_0_thresh[[l]][2]])/length(Sim_diff_2) + length(Sim_diff_2[Sim_diff_2<Comp_0_thresh[[l]][1]])/length(Sim_diff_2))
  Sim_diff_3 <- as.vector(total_wl[[3]][[l]]) - mean(total_wl[[3]][[l]])
  # sc 1 vs sc3
  percent_past_thresh[[4]][l] <- 100*(length(Sim_diff_3[Sim_diff_3>Comp_0_thresh[[l]][2]])/length(Sim_diff_3) + length(Sim_diff_3[Sim_diff_3<Comp_0_thresh[[l]][1]])/length(Sim_diff_3))
}
# If Statistics.txt already exists, please delete file before running the following line of code.
sink(file = "Statistics.txt")
for(i in 1:1){
  cat(noquote(paste("Observed Long Term Average Outflow\n")))
  cat(noquote(paste("\tLake Superior:",mean(unlist(outflowData[[2]])),"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",mean(unlist(outflowData[[4]])),"\n")))
  cat(noquote(paste("\tLake Erie:",mean(unlist(outflowData[[6]])),"\n")))
  
  cat(noquote(paste("Simulated Long Term Average Outflow\n")))
  cat(noquote(paste("\tLake Superior:",mean(total_out[[1]][[1]]),"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",mean(total_out[[1]][[2]]),"\n")))
  cat(noquote(paste("\tLake Erie:",mean(total_out[[1]][[3]]),"\n")))
  
  cat(noquote(paste("Observed Annual Average NBS (1950-2019)\n")))
  cat(noquote(paste("\tLake Superior:",mean(unlist(copulaData[[1]]-copulaData[[2]]+copulaData[[3]])),"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",mean(unlist(copulaData[[4]]-copulaData[[5]]+copulaData[[6]])),"\n")))
  cat(noquote(paste("\tLake Erie:",mean(unlist(copulaData[[7]]-copulaData[[8]]+copulaData[[9]])),"\n")))
  
  cat(noquote(paste("Simulated Average Annual NBS Scenario 1 (temp:SC1-RM)\n")))
  cat(noquote(paste("\tLake Superior:",mean(unlist(copulaPrediction[[1]][[1]]-copulaPrediction[[1]][[2]]+copulaPrediction[[1]][[3]])),"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",mean(unlist(copulaPrediction[[1]][[4]]-copulaPrediction[[1]][[5]]+copulaPrediction[[1]][[6]])),"\n")))
  cat(noquote(paste("\tLake Erie:",mean(unlist(copulaPrediction[[1]][[7]]-copulaPrediction[[1]][[8]]+copulaPrediction[[1]][[9]])),"\n")))
  
  cat(noquote(paste("Simulated Average Annual NBS Scenario 2 (temp:SC1-PM)\n")))
  cat(noquote(paste("\tLake Superior:",mean(unlist(copulaPrediction[[2]][[1]]-copulaPrediction[[2]][[2]]+copulaPrediction[[2]][[3]])),"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",mean(unlist(copulaPrediction[[2]][[4]]-copulaPrediction[[2]][[5]]+copulaPrediction[[2]][[6]])),"\n")))
  cat(noquote(paste("\tLake Erie:",mean(unlist(copulaPrediction[[2]][[7]]-copulaPrediction[[2]][[8]]+copulaPrediction[[2]][[9]])),"\n")))
  
  cat(noquote(paste("Simulated Average Annual NBS Scenario 3 (temp:SC1-asdsaM)\n")))
  cat(noquote(paste("\tLake Superior:",mean(unlist(copulaPrediction[[3]][[1]]-copulaPrediction[[3]][[2]]+copulaPrediction[[3]][[3]])),"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",mean(unlist(copulaPrediction[[3]][[4]]-copulaPrediction[[3]][[5]]+copulaPrediction[[3]][[6]])),"\n")))
  cat(noquote(paste("\tLake Erie:",mean(unlist(copulaPrediction[[3]][[7]]-copulaPrediction[[3]][[8]]+copulaPrediction[[3]][[9]])),"\n")))
  
  cat(noquote(paste("Observed Long Term Average Water Level (1950-2019)\n")))
  cat(noquote(paste("\tLake Superior:",mean(unlist(outflowData[[1]])),"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",mean(unlist(outflowData[[3]])),"\n")))
  cat(noquote(paste("\tLake Erie:",mean(unlist(outflowData[[5]])),"\n")))
  
  cat(noquote(paste("Simulated Long Term Average Water Level Scenario 1 (temp:SC1-RM)\n")))
  cat(noquote(paste("\tLake Superior:",long_wl_means_mean[[1]][[1]],"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",long_wl_means_mean[[1]][[2]],"\n")))
  cat(noquote(paste("\tLake Erie:",long_wl_means_mean[[1]][[3]],"\n")))
  
  cat(noquote(paste("Simulated Long Term Average Water Level Scenario 2 (temp:SC1-PM)\n")))
  cat(noquote(paste("\tLake Superior:",long_wl_means_mean[[2]][[1]],"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",long_wl_means_mean[[2]][[2]],"\n")))
  cat(noquote(paste("\tLake Erie:",long_wl_means_mean[[2]][[3]],"\n")))
  
  cat(noquote(paste("Simulated Long Term Average Water Level Scenario 3 (temp:SC1-asdsaM)\n")))
  cat(noquote(paste("\tLake Superior:",long_wl_means_mean[[3]][[1]],"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",long_wl_means_mean[[3]][[2]],"\n")))
  cat(noquote(paste("\tLake Erie:",long_wl_means_mean[[3]][[3]],"\n")))
  
  cat(noquote(paste("Observed Peak Month\n")))
  cat(noquote(paste("\tLake Superior:",months[which(observed_extrema[[1]][[1]]==max(observed_extrema[[1]][[1]]))],Weighted_extrema[[4]][[1]][1],"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",months[which(observed_extrema[[2]][[1]]==max(observed_extrema[[2]][[1]]))],Weighted_extrema[[4]][[2]][1],"\n")))
  cat(noquote(paste("\tLake Erie:",months[which(observed_extrema[[3]][[1]]==max(observed_extrema[[3]][[1]]))],Weighted_extrema[[4]][[3]][1],"\n")))
  
  cat(noquote(paste("Observed Trough Month\n")))
  cat(noquote(paste("\tLake Superior:",months[which(observed_extrema[[1]][[2]]==max(observed_extrema[[1]][[2]]))],Weighted_extrema[[4]][[1]][2],"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",months[which(observed_extrema[[2]][[2]]==max(observed_extrema[[2]][[2]]))],Weighted_extrema[[4]][[2]][2],"\n")))
  cat(noquote(paste("\tLake Erie:",months[which(observed_extrema[[3]][[2]]==max(observed_extrema[[3]][[2]]))],Weighted_extrema[[4]][[3]][2],"\n")))
  
  cat(noquote(paste("Simulated Peak Month Scenario 1\n")))
  cat(noquote(paste("\tLake Superior:",months[which(total_annual_extrema[[1]][[1]][[1]]==max(total_annual_extrema[[1]][[1]][[1]]))],Weighted_extrema[[1]][[1]][1],"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",months[which(total_annual_extrema[[1]][[2]][[1]]==max(total_annual_extrema[[1]][[2]][[1]]))],Weighted_extrema[[1]][[2]][1],"\n")))
  cat(noquote(paste("\tLake Erie:",months[which(total_annual_extrema[[1]][[3]][[1]]==max(total_annual_extrema[[1]][[3]][[1]]))],Weighted_extrema[[1]][[3]][1],"\n")))
  
  cat(noquote(paste("Simulated Trough Month Scenario 1\n")))
  cat(noquote(paste("\tLake Superior:",months[which(total_annual_extrema[[1]][[1]][[2]]==max(total_annual_extrema[[1]][[1]][[2]]))],Weighted_extrema[[1]][[1]][2],"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",months[which(total_annual_extrema[[1]][[2]][[2]]==max(total_annual_extrema[[1]][[2]][[2]]))],Weighted_extrema[[1]][[2]][2],"\n")))
  cat(noquote(paste("\tLake Erie:",months[which(total_annual_extrema[[1]][[3]][[2]]==max(total_annual_extrema[[1]][[3]][[2]]))],Weighted_extrema[[1]][[3]][2],"\n")))
  
  cat(noquote(paste("Simulated Peak Month Scenario 2\n")))
  cat(noquote(paste("\tLake Superior:",months[which(total_annual_extrema[[2]][[1]][[1]]==max(total_annual_extrema[[2]][[1]][[1]]))],Weighted_extrema[[2]][[1]][1],"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",months[which(total_annual_extrema[[2]][[2]][[1]]==max(total_annual_extrema[[2]][[2]][[1]]))],Weighted_extrema[[2]][[2]][1],"\n")))
  cat(noquote(paste("\tLake Erie:",months[which(total_annual_extrema[[2]][[3]][[1]]==max(total_annual_extrema[[2]][[3]][[1]]))],Weighted_extrema[[2]][[3]][1],"\n")))
  
  cat(noquote(paste("Simulated Trough Month Scenario 2\n")))
  cat(noquote(paste("\tLake Superior:",months[which(total_annual_extrema[[2]][[1]][[2]]==max(total_annual_extrema[[2]][[1]][[2]]))],Weighted_extrema[[2]][[1]][2],"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",months[which(total_annual_extrema[[2]][[2]][[2]]==max(total_annual_extrema[[2]][[2]][[2]]))],Weighted_extrema[[2]][[2]][2],"\n")))
  cat(noquote(paste("\tLake Erie:",months[which(total_annual_extrema[[2]][[3]][[2]]==max(total_annual_extrema[[2]][[3]][[2]]))],Weighted_extrema[[2]][[3]][2],"\n")))
  
  cat(noquote(paste("Simulated Peak Month Scenario 3\n")))
  cat(noquote(paste("\tLake Superior:",months[which(total_annual_extrema[[3]][[1]][[1]]==max(total_annual_extrema[[3]][[1]][[1]]))],Weighted_extrema[[3]][[1]][1],"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",months[which(total_annual_extrema[[3]][[2]][[1]]==max(total_annual_extrema[[3]][[2]][[1]]))],Weighted_extrema[[3]][[2]][1],"\n")))
  cat(noquote(paste("\tLake Erie:",months[which(total_annual_extrema[[3]][[3]][[1]]==max(total_annual_extrema[[3]][[3]][[1]]))],Weighted_extrema[[3]][[3]][1],"\n")))
  
  cat(noquote(paste("Simulated Trough Month Scenario 3\n")))
  cat(noquote(paste("\tLake Superior:",months[which(total_annual_extrema[[3]][[1]][[2]]==max(total_annual_extrema[[3]][[1]][[2]]))],Weighted_extrema[[3]][[1]][2],"\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",months[which(total_annual_extrema[[3]][[2]][[2]]==max(total_annual_extrema[[3]][[2]][[2]]))],Weighted_extrema[[3]][[2]][2],"\n")))
  cat(noquote(paste("\tLake Erie:",months[which(total_annual_extrema[[3]][[3]][[2]]==max(total_annual_extrema[[3]][[3]][[2]]))],Weighted_extrema[[3]][[3]][2],"\n")))
  
  cat(noquote(paste("Scenario 1 Percentage Past Observed Threshold\n")))
  cat(noquote(paste("\tLake Superior:",percent_past_thresh[[1]][1],"%\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",percent_past_thresh[[1]][2],"%\n")))
  cat(noquote(paste("\tLake Erie:",percent_past_thresh[[1]][3],"%\n")))
  
  cat(noquote(paste("Scenario 1 Percentage Past Scenario 1 Threshold\n")))
  cat(noquote(paste("\tLake Superior:",percent_past_thresh[[2]][1],"%\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",percent_past_thresh[[2]][2],"%\n")))
  cat(noquote(paste("\tLake Erie:",percent_past_thresh[[2]][3],"%\n")))
  
  cat(noquote(paste("Scenario 2 Percentage Past Scenario 1 Threshold\n")))
  cat(noquote(paste("\tLake Superior:",percent_past_thresh[[3]][1],"%\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",percent_past_thresh[[3]][2],"%\n")))
  cat(noquote(paste("\tLake Erie:",percent_past_thresh[[3]][3],"%\n")))
  
  cat(noquote(paste("Scenario 3 Percentage Past Scenario 1 Threshold\n")))
  cat(noquote(paste("\tLake Superior:",percent_past_thresh[[4]][1],"%\n")))
  cat(noquote(paste("\tLake Michigan-Huron:",percent_past_thresh[[4]][2],"%\n")))
  cat(noquote(paste("\tLake Erie:",percent_past_thresh[[4]][3],"%\n")))
}
sink(file = NULL)

## Outflow & Water Level Dual Plot
jpeg(paste(plot.dir,"new_figure_1.jpg", sep=""), height = 5, width = 7, units = "in", res = 72*4)
par(mfrow=c(3,2))
par(mar = c(1,1,1.5,1))
par(oma = c(2,5.0,0,4.5))

## Outflow ensemble plotted needs to be the same as the water level ensemble member

## Lake Superior
plot(apply(total_wl[[1]][[1]],2,mean), type = "n", ylim = c(182.5, 184.5), axes = F); box(); axis(2)
for(i in seq(1,gen,10)){lines(total_wl[[1]][[1]][i,], col = "grey")}
lines(as.vector(t(outflowData[[1]])))
lines(total_wl[[1]][[1]][400,], col = "red")
axis(1, labels = F)
mtext("Lake Superior", side = 2, line = 3, cex = 0.75)

plot(apply(outflow_plot[[1]],2,mean), type = "n", ylim = c(0, 6000), axes = F); box(); axis(4)
for(i in seq(1,gen,10)){lines(outflow_plot[[1]][i,], col = "grey")}
lines(outflow_plot[[1]][400,], col = "red")
lines(as.vector(t(outflowData[[2]])))
axis(1, labels = F)

## Lake Michigan-Huron

plot(apply(total_wl[[1]][[2]],2,mean), type = "n", ylim = c(175, 178), axes = F); box()
for(i in seq(1,gen,10)){lines(total_wl[[1]][[2]][i,], col = "grey")}
lines(as.vector(t(outflowData[[3]])))
lines(total_wl[[1]][[2]][400,], col = "red")
axis(1, labels = F); axis(2)
mtext("Lake Michigan-Huron", side = 2, line = 3, cex = 0.7)
mtext("Water surface elevation (m)", side = 2, line = 4.5, cex = 0.85)

plot(apply(outflow_plot[[2]],2,mean), type = "n", ylim = c(3000, 8000), axes = F); box(); axis(4)
for(i in seq(1,gen,10)){lines(outflow_plot[[2]][i,], col = "grey")}
lines(outflow_plot[[2]][400,], col = "red")
lines(as.vector(t(outflowData[[4]])))
axis(1, labels = F)
mtext("Lake outflow (cms)", side = 4, line = 3.5, cex = 0.85)

## Lake Erie

plot(apply(total_wl[[1]][[3]],2,mean), type = "n", ylim = c(173, 175.5))
for(i in seq(1,gen,10)){lines(total_wl[[1]][[3]][i,], col = "grey")}
lines(as.vector(t(outflowData[[5]])))
lines(total_wl[[1]][[3]][400,], col = "red")
mtext("Lake Erie", side = 2, line = 3, cex = 0.7)

plot(apply(outflow_plot[[3]],2,mean), type = "n", ylim = c(4000, 8000), axes = F); box()
for(i in seq(1,gen,10)){lines(outflow_plot[[3]][i,], col = "grey")}
lines(outflow_plot[[3]][400,], col = "red")
lines(as.vector(t(outflowData[[6]])))
axis(4); axis(1)
dev.off()