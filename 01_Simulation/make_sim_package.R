# make sim package
rm(list=ls())

# library includes function to determine number of cores on machine
library(parallel)
setwd("")

###################################################################
# 1) PARAMETERS TO VARY
# dispersal rate
disp_rate <- c(0.001,0.01,0.1)
# weighting of upstream movement probability (w_up; 0: no upstream movement; 1: upstream and downstream movements are equally likely)
w_up <- c(0,0.5,1)
# use average carrying capacity (K_scale=0) or scale with catchement size (K_scale = 1)
K_scale <- c(0,1)

###################################################################
# 2) FIXED PARAMATERS
# simulation time
sim_time <- 10000
# number of replicates
max_runs <- 10
# no. neural alleles
no_neutral_alleles <- 100
# mutation rate neutral allels
mut_rate_neutral <- 0.0001
# habitat capacity
capacity <- 1000
# fertility
lambda_null <- 2
# environmental stochasticity
sigma <- 0
# Allee effect strength
alpha <- 0
# random patch extinction probability
epsilon <- 0
# dispersal mortality downstream
mu0_down <- 0
# dispersal mortality upstream
mu0_up <- 0
# mutational model used (0 --> random; 1 --> stepwise)
mut_model_smm <- 1

##################################################################
# name of simulation program (is abbreviated as in the "ps" output!)
SIMNAME <- "pop_gen_gammaridae_v0"
SIMNAME_PS <- "pop_gen_gammari"

##################################################################
# save the present working directory
HOMEWD <- getwd()
setwd(HOMEWD)
# get the directory with the individual simulation
SIMWD <- paste(HOMEWD,"/indiv_sim",sep="")

##################################################################
# number of cores on machine
NCORES <- detectCores()-1
# calculate number of simulations
NTOTALSIMS <- length(disp_rate) * length(w_up) * length(K_scale)

##################################################################
# loop over all parameters and start the simulations

# intialize simulation counter
sim_cnt <- 0

# loop 1: dispersal rates
for (disp_rate_cnt in 1:length(disp_rate)){
  # set home as actwd
  ACTWD_0 <- HOMEWD
  setwd(ACTWD_0)
  # make the directory
  system(paste("mkdir disp_",disp_rate[disp_rate_cnt],sep=""))
  
  # loop 2: w_up
  for (w_up_cnt in 1:length(w_up)){
    # set home as actwd
    ACTWD_1 <- paste(ACTWD_0,"/disp_",disp_rate[disp_rate_cnt],sep="")
    setwd(ACTWD_1)
    # make the directory
    system(paste("mkdir w_up_",w_up[w_up_cnt],sep=""))
    
    # loop 3: K scale
    for (K_scale_cnt in 1:length(K_scale)){
      # change into local extinctions directory
      ACTWD_2 <- paste(ACTWD_1,"/w_up_",w_up[w_up_cnt],sep="")
      setwd(ACTWD_2)
      
      # make new directory
      system(paste("mkdir K_scale_",K_scale[K_scale_cnt],sep=""))
      
      # copy the source code and all dependecies into this folder
      system(paste("cp -r ",SIMWD,"/. K_scale_",K_scale[K_scale_cnt],sep=""))
      
      # change into that folder
      ACTWD_3 <- paste(ACTWD_2,"/K_scale_",K_scale[K_scale_cnt],sep="")
      setwd(ACTWD_3)
      
      # write parameter input file
      parameter_values <- matrix(ncol=1,nrow=32)
      parameter_values[1,]  <- "inputfile for parameters"
      parameter_values[2,]  <- ""
      parameter_values[3,]  <- "simulation time (sim_time)"
      parameter_values[4,]  <- sim_time
      parameter_values[5,]  <- "number of repeats (max_runs)"
      parameter_values[6,]  <- max_runs
      parameter_values[7,]  <- "dispersal rate (disp_rate)"
      parameter_values[8,]  <- disp_rate[disp_rate_cnt]
      parameter_values[9,]  <- "no. of alleles (no_neutral_alleles)"
      parameter_values[10,]  <- no_neutral_alleles
      parameter_values[11,]  <- "mutation rate neural loci (mut_rate_neutral)"
      parameter_values[12,]  <- mut_rate_neutral
      parameter_values[13,]  <- "habitat capacity (capacity)"
      parameter_values[14,]  <- capacity
      parameter_values[15,]  <- "fertility (lambda_null)"
      parameter_values[16,]  <- lambda_null
      parameter_values[17,]  <- "environmental stochasticity (sigma)"
      parameter_values[18,]  <- sigma
      parameter_values[19,]  <- "Allee effect strength (alpha)"
      parameter_values[20,]  <- alpha
      parameter_values[21,]  <- "random patch extinction probability (epsilon)"
      parameter_values[22,]  <- epsilon
      parameter_values[23,]  <- "per step (distance in meters) dispersal mortality downstream (mu0_down)"
      parameter_values[24,]  <- mu0_down
      parameter_values[25,]  <- "per step (distance in meters) dispersal mortality upstream (mu0_up)"
      parameter_values[26,]  <- mu0_up
      parameter_values[27,]  <- "weighting of upstream movement probability (w_up; 0: no upstream movement; 1: upstream and downstream movements are equally likely)"
      parameter_values[28,]  <- w_up[w_up_cnt]
      parameter_values[29,]  <- "use average carrying capacity (0) or scale with catchement size (1)"
      parameter_values[30,]  <- K_scale[K_scale_cnt]
      parameter_values[31,]  <- "mutational model used (0 --> random; 1 --> stepwise)"
      parameter_values[32,]  <- mut_model_smm
      write.table(file=paste(ACTWD_3,"/input/parameters.in",sep=""),parameter_values,col.names=F,row.names=F,quote=F)
      
      # now compile the code
      system("make")
      
      # check how many simulations are running and continue only if a core is available
      repeat{
        # get the number of running simulations
        n_act_sims <- length(system(paste("ps -A | grep ",SIMNAME_PS,sep=""),intern=T))
        
        # check number of sims vs number of cores
        if(n_act_sims < NCORES){
          # if less sims are running than we have cores
          break
        }else{
          # else wait for a second
          Sys.sleep(1)
          #print("in loop!")
        }
      }
      
      # start simulation
      system(paste("nohup ./",SIMNAME," &",sep=""))
      # increas simulation counter
      sim_cnt <- sim_cnt + 1
      # print this info
      print(paste(sim_cnt," of ",NTOTALSIMS))
      
    } # end loop: K_scale
  } # end loop: w_ip
} # end loop: dispersal rates
##################################################################