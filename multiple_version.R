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
mv.dir 	<- paste(main.dir, "/multiple_version/", sep = "") 	# L2SWBM simulations
data.dir 	<- paste(main.dir, "/historical_data/", sep = "") 	# historical data
plot.dir 	<- paste(main.dir, "/plots/", sep = "") 			# outputed plots

ITERNUM = 1000
total_out_list = list()
total_wl_list = list()

####################################################
#### Constants
####################################################
dice <- c(0)
for(dice_index in dice){
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
    paste('superiorPrecip_',dice_index,'.csv',sep = ""),
    paste('superiorEvap_',dice_index,'.csv',sep = ""),
    paste('superiorRunoff_',dice_index,'.csv',sep = ""),
    paste('miHuronPrecip_',dice_index,'.csv',sep = ""),
    paste('miHuronEvap_',dice_index,'.csv',sep = ""),
    paste('miHuronRunoff_',dice_index,'.csv',sep = ""),
    paste('eriePrecip_',dice_index,'.csv',sep = ""),
    paste('erieEvap_',dice_index,'.csv',sep = ""),
    paste('erieRunoff_',dice_index,'.csv',sep = ""))
  
  ## Read in median values of P, E, R for each lake from L2SWBM 
  ## ....and transform vectors into matrices (col -- month & row -- year)
  for(i in 1:9) copulaData[[i]] <- fread(paste(mv.dir, glwbDataFiles[[i]], sep = ""), select = "Median")
  
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
                           generate_predictions(readRDS("rvine_copula_70.rds"),ITERNUM))
  
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
                          rep(as.numeric(outflowData[[3]][addYear, 1]),gen),
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
        outflow[[1]] <- (outflow[[1]]/S_SF)*(SPM31)/(ITERNUM) # months with 31 days
        temp_div_sup <- (as.ts(diversion[[1]])/S_SF)*(SPM31)/(ITERNUM)
      } else if (m %in% c(4,6,9,11)){
        outflow[[1]] <- (outflow[[1]]/S_SF)*(SPM30)/(ITERNUM) # months with 30 days
        temp_div_sup <- (as.ts(diversion[[1]])/S_SF)*(SPM30)/(ITERNUM)
      } else{
        outflow[[1]] <- (outflow[[1]]/S_SF)*(SPMF)/(ITERNUM) # February
        temp_div_sup <- (as.ts(diversion[[1]])/S_SF)*(SPMF)/(ITERNUM)
      }
      # add water level change to update current simulated water level
      current_level[[1]] <- current_level[[1]] + temp_div_sup[i]/ITERNUM + (copulaPrediction[[s]][[1]][c(((y-1)*gen+1):(y*gen)),m] - copulaPrediction[[s]][[2]][c(((y-1)*gen+1):(y*gen)),m] + copulaPrediction[[s]][[3]][c(((y-1)*gen+1):(y*gen)),m] - outflow[[1]])/ITERNUM
      water_level[[1]] <- cbind(water_level[[1]],current_level[[1]])
      
      # Lake Michigan-Huron
      outflow[[2]] <- exp(intercept_consts[[2]][c(1:gen),m] + stage_consts[[2]][c(1:gen)]*log(current_level[[2]]-sill_adjustments[2]) + fall_consts[[1]][c(1:gen)]*log(current_level[[2]]-current_level[[3]]) + rnorm(gen,0,vol_consts[[2]][c(1:gen)]))
      outflow_calculation[[2]] <- cbind(outflow_calculation[[2]],outflow[[2]])
      outflow_plot[[2]] <- cbind(outflow_plot[[2]],outflow[[2]])
      # Convert cubic cm per sec of flow --> water level per month
      if (m %in% c(1,3,5,7,8,10,12)){
        outflow[[2]] <- (outflow[[2]]/MH_SF)*(SPM31)/ITERNUM # months with 31 days
        outflow_calculation[[1]][,ncol(outflow_calculation[[1]])] <- (outflow_calculation[[1]][,ncol(outflow_calculation[[1]])]/MH_SF)*SPM31/ITERNUM
        temp_div_mh <- (as.ts(diversion[[2]])/MH_SF)*SPM31/ITERNUM
      } else if (m %in% c(4,6,9,11)){
        outflow[[2]] <- (outflow[[2]]/MH_SF)*(SPM30)/ITERNUM # months with 30 days
        outflow_calculation[[1]][,ncol(outflow_calculation[[1]])] <- (outflow_calculation[[1]][,ncol(outflow_calculation[[1]])]/MH_SF)*SPM30/ITERNUM
        temp_div_mh <- (as.ts(diversion[[2]])/MH_SF)*SPM30/ITERNUM
      } else{
        outflow[[2]] <- (outflow[[2]]/MH_SF)*(SPMF)/ITERNUM # February
        outflow_calculation[[1]][,ncol(outflow_calculation[[1]])] <- (outflow_calculation[[1]][,ncol(outflow_calculation[[1]])]/MH_SF)*SPMF/ITERNUM
        temp_div_mh <- (as.ts(diversion[[2]])/MH_SF)*SPMF/ITERNUM
      }
      # add water level change to update current simulated water level
      current_level[[2]] <- current_level[[2]] + temp_div_mh[i]/ITERNUM + (outflow_calculation[[1]][,ncol(outflow_calculation[[1]])]*0.8 + copulaPrediction[[s]][[4]][c(((y-1)*gen+1):(y*gen)),m] - copulaPrediction[[s]][[5]][c(((y-1)*gen+1):(y*gen)),m] + copulaPrediction[[s]][[6]][c(((y-1)*gen+1):(y*gen)),m] - outflow[[2]])/ITERNUM
      water_level[[2]] <- cbind(water_level[[2]],current_level[[2]])
      
      # Lake Erie
      outflow[[3]] <- exp(intercept_consts[[3]][c(1:gen),m] + stage_consts[[3]][c(1:gen)]*log(current_level[[3]]-sill_adjustments[3]) + rnorm(gen,0,vol_consts[[3]][c(1:gen)]))
      outflow_calculation[[3]] <- cbind(outflow_calculation[[3]],outflow[[3]])
      outflow_plot[[3]] <- cbind(outflow_plot[[3]],outflow[[3]])
      # Convert cubic cm per sec of flow --> water level per month
      if (m %in% c(1,3,5,7,8,10,12)){
        outflow[[3]] <- (outflow[[3]]/E_SF)*(SPM31)/ITERNUM # months with 31 days
        outflow_calculation[[2]][,ncol(outflow_calculation[[2]])] <- (outflow_calculation[[2]][,ncol(outflow_calculation[[2]])]/E_SF)*SPM31/ITERNUM
        temp_div_e <- (as.ts(diversion[[3]])/E_SF)*SPM31/ITERNUM
      } else if (m %in% c(4,6,9,11)){
        outflow[[3]] <- (outflow[[3]]/E_SF)*(SPM30)/ITERNUM # months with 30 days
        outflow_calculation[[2]][,ncol(outflow_calculation[[2]])] <- (outflow_calculation[[2]][,ncol(outflow_calculation[[2]])]/E_SF)*SPM30/ITERNUM
        temp_div_e <- (as.ts(diversion[[3]])/E_SF)*SPM30/ITERNUM
      } else{
        outflow[[3]] <- (outflow[[3]]/E_SF)*(SPMF)/ITERNUM # February
        outflow_calculation[[2]][,ncol(outflow_calculation[[2]])] <- (outflow_calculation[[2]][,ncol(outflow_calculation[[2]])]/E_SF)*SPMF/ITERNUM
        temp_div_e <- (as.ts(diversion[[3]])/E_SF)*SPMF/ITERNUM
      }
      # add water level change to update current simulated water level
      current_level[[3]] <- current_level[[3]] + temp_div_e[i]/ITERNUM + (outflow_calculation[[2]][,ncol(outflow_calculation[[2]])] + copulaPrediction[[s]][[7]][c(((y-1)*gen+1):(y*gen)),m] - copulaPrediction[[s]][[8]][c(((y-1)*gen+1):(y*gen)),m] + copulaPrediction[[s]][[9]][c(((y-1)*gen+1):(y*gen)),m] - outflow[[3]])/ITERNUM
      water_level[[3]] <- cbind(water_level[[3]],current_level[[3]])
    }
    total_out[[s]] <- outflow_plot
    total_wl[[s]] <- water_level
  }

  total_out_list <- c(total_out_list, total_out)
  total_wl_list <- c(total_wl_list, total_wl)
  
}
save(total_out_list, total_wl_list, file = "data10.RData")

## Water level Plots
Wl_plot_name <- c("scenario1_wl_10.pdf","scenario2_wl_10.pdf","scenario3_wl_10.pdf")
for(s in 1:3){
  # water level plot
  pdf(paste(plot.dir, Wl_plot_name[[s]], sep=""), width = 14, height = 20, paper = "special",onefile = FALSE) # create pdf
  par(mfrow = c(3, 1))
    plot(total_wl_list[[s]][[1]][1,],xlab="Month",ylab="Water level (m)",ylim=c(182,185),type='l',col='blue',main='Superior Water levels')
    for(i in 2:gen) {
      for(index_i in 0:1){
        current = index_i*3 + s
        points(total_wl_list[[current]][[1]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
      }
    }
    if(s==1) points(as.vector(t(outflowData[[1]][c(1:50),])),type='l',lwd=3)
    
    plot(total_wl_list[[s]][[2]][1,],xlab="Month",ylab="Water level (m)",ylim=c(175,178),type='l',col='blue',main='Michigan-Huron Water levels')
    for(i in 2:gen) {
      for(index_i in 0:1){
        current = index_i*3 + s
        points(total_wl_list[[current]][[2]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
      }
    }
    if(s==1) points(as.vector(t(outflowData[[3]][c(1:50),])),type='l',lwd=3)
    
    plot(total_wl_list[[s]][[3]][1,],xlab="Month",ylab="Water level (m)",ylim=c(173,176),type='l',col='blue',main='Erie Water levels')
    for(i in 2:gen) {
      for(index_i in 0:1){
        current = index_i*3 + s
        points(total_wl_list[[current]][[3]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
      }
    }
    if(s==1) points(as.vector(t(outflowData[[5]][c(1:50),])),type='l',lwd=3)
    dev.off()
}

## Outflow Plots 
out_plot_name <- c("scenario1_out_10.pdf","scenario2_out_10.pdf","scenario3_out_10.pdf")
for(s in 1:3){
  pdf(paste(plot.dir, out_plot_name[s], sep=""), width = 14, height = 18, paper = "special",onefile = TRUE) # create pdf
  par(mfrow = c(3,1))
    par(mar=c(0,5,5,5))
    plot(total_out_list[[s]][[1]][1,],xlab="",ylab="",ylim=c(0,8000),type='l',col='red',main='',cex.main=2,cex.axis=1.5,cex.lab=2,xaxt='n')
    title(main="Superior Outflow",line = -2,cex.main=2.5)
    for(i in 2:gen) {
      for(index_i in 0:1){
        current = index_i*3 + s
        points(total_out_list[[current]][[1]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
      }
    }
    if(s==1) points(as.vector(t(outflowData[[2]][c(1:50),])),type='l',lwd=3)
    
    par(mar=c(0,5,0,5))
    plot(total_out_list[[current]][[2]][1,],xlab="",ylab="",ylim=c(2000,10000),type='l',col='red',main='',cex.main=2,cex.axis=1.5,cex.lab=2,xaxt='n')
    title(ylab=expression(Outflow ~ (cm^"3")),line=2,cex.lab=2)
    title(main="Michigan-Huron Outflow",line = -2,cex.main=2.5)
    for(i in 2:gen) {
      for(index_i in 0:1){
        current = index_i*3 + s
        points(total_out_list[[current]][[2]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
      }
    }
     
    if(s == 1) points(as.vector(t(outflowData[[4]][c(1:50),])),type='l',lwd=3)

    par(mar=c(5,5,0,5))
    plot(total_out_list[[current]][[3]][1,],xlab="Month",ylab="",ylim=c(2000,10000),type='l',col='red',main='',cex.main=2,cex.axis=1.5,cex.lab=2)
    title(main="Erie Outflow",line = -2,cex.main=2.5)
    for(i in 2:gen) {
      for(index_i in 0:1){
        current = index_i*3 + s
        points(total_out_list[[current]][[3]][i,],type='l',col=rainbow(gen)[sample.int(gen,1)])
      }
    }
    if(s==1 ) points(as.vector(t(outflowData[[6]][c(1:50),])),type='l',lwd=3)

  dev.off()
} 
gen=999
