# PopGenNet data preparation
# 2020-12-10
# Author: Roman Alther, Eawag, Duebendorf, Switzerland
# Works in R ver. 3.6.1 and 4.0.2 (tested on GNU/Linux, MacOS, Windows)
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

#### PREPARATION -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
rm(list=ls())

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
WD <- getwd() # save current directory (should correspond to 03_Analysis, source script from there)

##** Settings ####
# Threshold of individuals to include
n_ind <- 15 # minimal number of sequenced individuals for microsatellite data
n <- 10 # number of replicates for each parameter combination
gen_thres <- n_ind # minimal number of sequenced individuals for microsatellite data

tot_time <- Sys.time() # start timer

##** Packages ####
library(hierfstat) # works with version 0.04-22
library(igraph) # works with version 1.2.4.2

##** Functions ####

# Function to standardize a variable to the range [0,1]
stand <- function(variable){
  return((variable-min(variable))/ (max(variable)-min(variable)))
}

#### DATA LOADING-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
##** Microsat data ####
data <- read.table("Gammarus_all_microsat_data_Rhein_2018.txt", header=TRUE)
##** Spatial data ####
load("Gfos_Network.RData")

#### EMPIRICAL DATA PREPARATION -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
##** Microsat data corrections ####
data <- data[!is.na(data$species),]
data <- data[data$species=="fosA" | data$species=="fosB",]
data <- data[data$catchment=="RHEIN",]
# ATTENTION: Anja23 is same location as Anja6
data$locality[data$locality=="Anja_23"] <- "Anja_6"
# ATTENTION: Manebach_1000m has wrong coordinates
data$ykoord[data$locality=="Manebach_1000m"] <- 280926
data <- droplevels(data)

# Prepare dataframe for map plotting
matrix <- data.frame(tapply(data$species, list(data$locality, data$species), length))
matrix$xkoord <- tapply(data$xkoord, list(data$locality), mean) # using mean to extract coordinate per site (all identical). A resulting value with decimals would indicate an issue...
matrix$ykoord <- tapply(data$ykoord, list(data$locality), mean)
matrix$locality <- levels(data$locality)

#### Sites/rows in the matrix to be considered in analysis (containing data and sufficient individuals)
sites_A <- which(!is.na(matrix$fosA) & matrix$fosA>n_ind)
sites_B <- which(!is.na(matrix$fosB) & matrix$fosB>n_ind)

##** Microsat dataset for each species and total with at least threshold number of individuals ####
A <- data[data$species=="fosA" & is.element(data$locality, matrix$locality[sites_A]),c(37,26:35)]
A <- droplevels(A)
B <- data[data$species=="fosB" & is.element(data$locality, matrix$locality[sites_B]),c(37,26:35)]
B <- droplevels(B)
#### Two catchments contain more than one site and are reduced to one site
A_red <- A[A$locality!="Hornbach_Horn_1000m500m"&A$locality!="Manebach_0m"&A$locality!="Eschlibach_0m",] # Hornbach_Horn_1000m500m is the same subcatchment as Anja_16, also Manebach_0m and Eschlibach_0m are the same subcatchment as Speckbach_0m50m
A_red <- droplevels(A_red)
sites_A_red <-  match(levels(A_red$locality),matrix$locality)
rm(list=c("data","A","sites_A","sites_B","sites_A_red"))

#### ANALYSIS NETWORK DATA -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
#### Which network to use

net <- Gfos_Network
rm(list=c("Gfos_Network"))

##** Calculate network metrics ####
betw_undir <- betweenness(net, v=V(net), directed = F) #betweenness, undirected
betw_dir <- betweenness(net, v=V(net), directed = T) #betweenness, directed
clos_undir <- closeness(net, v=V(net), mode="all", normalized = T) #closeness, undirected
oldw <- getOption("warn")
options(warn=-1) # since calculating directed closeness will result in warning 'closeness centrality is not well-defined for disconnected graphs'
clos_dir <- closeness(net, v=V(net), mode="out", normalized = T) #closeness, directed; "out" makes most downstream node most central = highest value
options(warn=oldw)
rm(list=c("oldw"))
degree_undir <- centr_degree(net, mode="all")

#### Undirected network metrics for each species
betw_undir_B <- betw_undir[match(microsite_B,V(net)$name)]
clos_undir_B <- clos_undir[match(microsite_B,V(net)$name)]
degree_undir_B <- degree_undir$res[match(microsite_B,V(net)$name)]
betw_undir_A_red <- betw_undir[match(microsite_A_red,V(net)$name)]
clos_undir_A_red <- clos_undir[match(microsite_A_red,V(net)$name)]
degree_undir_A_red <- degree_undir$res[match(microsite_A_red,V(net)$name)]

#### Directed network metrics for each species
betw_dir_B <- betw_dir[match(microsite_B,V(net)$name)]
clos_dir_B <- clos_dir[match(microsite_B,V(net)$name)]
betw_dir_A_red <- betw_dir[match(microsite_A_red,V(net)$name)]
clos_dir_A_red <- clos_dir[match(microsite_A_red,V(net)$name)]
rm(list=c("betw_dir","clos_dir"))

#### Spatial distance between populations with Fos A and between populations with Fos B
DIST_B <- distances(net, v=V(net)[match(microsite_B,V(net)$name)], to=V(net)[match(microsite_B,V(net)$name)], weights=E(net))
DIST_A_red <- distances(net, v=V(net)[match(microsite_A_red,V(net)$name)], to=V(net)[match(microsite_A_red,V(net)$name)], weights=E(net))

#### Distance matrix to vector
dist_B <- DIST_B[upper.tri(DIST_B)]
dist_A_red <- DIST_A_red[upper.tri(DIST_A_red)]
rm(list=c("DIST_A_red","DIST_B"))

#### Calculate upstream distance from outlet vertex
outlet_vertex <- which.min(centr_degree(net, mode = "out")$res)

#### Total catchment area
catch_B <- V(net)$Total_Catch[match(microsite_B,V(net)$name)]
catch_A_red <- V(net)$Total_Catch[match(microsite_A_red,V(net)$name)]

##** Upstream dist ####
updist <- distances(net, v = V(net)[outlet_vertex], mode = "in", weights = E(net)$weight)

#updist_A <- updist[match(microsite_A,V(net)$name)]
updist_B <- updist[match(microsite_B,V(net)$name)]
updist_A_red <- updist[match(microsite_A_red,V(net)$name)]

####*************************************************************####

#### ANALYSIS EMPIRICAL DATA -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
##** Allelic richness per locus per site ####
Ar_B <- allelic.richness(B, min.n=n_ind)
Ar_A_red <- allelic.richness(A_red, min.n=n_ind)

##** Mean allelic richness per site ####
meanAr_B <- apply(Ar_B$Ar, 2, mean, na.rm=TRUE)
meanAr_A_red <- apply(Ar_A_red$Ar, 2, mean, na.rm=TRUE)

##** Heterozygosity ####
Basic_B <- basic.stats(B)
Basic_A_red <- basic.stats(A_red)

##** Mean observed heterozygosity per site ####
meanHo_B <- apply(Basic_B$Ho, 2, mean, na.rm=TRUE)
meanHo_A_red <- apply(Basic_A_red$Ho, 2, mean, na.rm=TRUE)

#### COMBINED EMPIRICAL DATAFRAME -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
##** Site dataframe preparation ####
# DATA table construction
DATA_B <- cbind(meanAr_B, meanHo_B, betw_undir_B, betw_dir_B, clos_undir_B, clos_dir_B, degree_undir_B, updist_B, catch_B, match_B)
colnames(DATA_B) <- c("meanAr","meanHo","betw_undir","betw_dir","clos_undir","clos_dir","degree","updist","catch","network_match")
DATA_B <- data.frame(DATA_B)

DATA_A_red <- cbind(meanAr_A_red, meanHo_A_red, betw_undir_A_red, betw_dir_A_red, clos_undir_A_red, clos_dir_A_red, degree_undir_A_red, updist_A_red, catch_A_red, match_A_red)
colnames(DATA_A_red) <- c("meanAr","meanHo","betw_undir","betw_dir","clos_undir","clos_dir","degree","updist","catch","network_match")
DATA_A_red <- data.frame(DATA_A_red)

DATA_B$spec <- "B"
DATA_A_red$spec <- "A"

DATA <- rbind(DATA_B,DATA_A_red)
DATA$spec <- as.factor(DATA$spec)

rm(list=c("DATA_A_red","DATA_B","Ar_A_red","Ar_B", "Basic_A_red","Basic_B","catch_A_red","catch_B"))

# Prepare ordered DATA with transformed values
DATA <- DATA[ order(DATA$updist), ]
DATA$log_updist <- log(DATA$updist)
DATA$log_catch <- log(DATA$catch)
DATA$std_clos_undir <- stand(DATA$clos_undir)
DATA$std_clos_dir <- stand(DATA$clos_dir)

##** Network data check ####
# Outlet vertex
if(V(net)$name[outlet_vertex]!="107471"){stop("The outlet vertex name is not '107471', so something went wrong!")}
# Headwater vertex (most upstream vertex)
# dimnames(updist)[[2]][which.max(updist)]

# Network data  preparation ####
#+++ Full network ####
# upstream distance of most downstream site (should hence be zero)
if(updist[which(V(net)$name==107471)]!=0){stop("Something is wrong, since node EZGNR 107471 has usptream distance bigger than 0.")}
# Do the EZGNR node names of full network correspond to the preloaded adjacency matrix node names?
if(!all(V(net)$name==adj.mat.names)){stop("Something is wrong, since EZGNR node names do not correspond to preloaded adjacency matrix.")}

#+++ Gfos A network ####
emp_GfosA <- which(site_coord$gendata==1|site_coord$gendata==3)
modsite_GfosA <- site_coord$EZGNR[emp_GfosA] # what are the EZGNR IDs of those?
if(!all(as.character(microsite_A_red)%in%modsite_GfosA)){stop("Something is wrong. Not all Gfos A microsat sites are represented in modelled sites.")}
SS_GfosA <- formatC(emp_GfosA-1, width=4, flag="0") # create the corresponding "x" number to compare to model output

#+++ Gfos B network ####
emp_GfosB <- which(site_coord$gendata==2|site_coord$gendata==3)
modsite_GfosB <- site_coord$EZGNR[emp_GfosB] # what are the EZGNR IDs of those?
if(!all(as.character(microsite_B)%in%modsite_GfosB)){stop("Something is wrong. Not all Gfos B microsat sites are represented in modelled sites.")}
SS_GfosB <- formatC(emp_GfosB-1, width=4, flag="0") # create the corresponding "x" number to compare to model output

rm(list=c("emp_GfosA","emp_GfosB","SS_GfosA","SS_GfosB"))

#+++ Simulation network ####
empiricaldata <- which(!is.na(site_coord$gendata)) # for which sites from all 2401 sites do we have genetic data?
modsite <- site_coord$EZGNR[empiricaldata] # what are the EZGNR IDs of those?
SS_Mod <- formatC(empiricaldata-1, width=4, flag="0") # create the corresponding "x" number to compare to model output

match_Mod <- match(modsite,site_coord$EZGNR) # Relate microsat data to site_coord

# are all modsites matching a site with microsat data?
if(!all(match(modsite,V(net)$name)%in%which(V(net)$microsat!=0))){stop("Something is wrong, since not all simulation sites match an empirical site.")}

####*************************************************************####

#### ANALYSIS MODEL DATA -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
# create color palette
# pal <- viridis(n)

# define parameter space
D <- c("disp_0.001","disp_0.01","disp_0.1")
D_label <- c("0.001","0.01","0.1")
W <- c("w_up_0","w_up_0.5","w_up_1")
W_label <- c("0","0.5","1")
K <- c("K_scale_0","K_scale_1")
K_label <- c("0","1")

# create folder to save output
dir.create(paste0(WD,"/Output/"), showWarnings=F)

# Looping over parameter space ####
r <- 0

for (d in 1:length(D)){ # looping over dispersal rates 
  
  for (w in 1:length(W)){ # looping over dispersal directionalities
    
    for (k in 1:length(K)){ # looping over carrying capacities
      r <- r+1
      setwd(paste0("../01_Simulation/",D[d],"/",W[w],"/",K[k]))
      
      parameter_time <- Sys.time()
      cat("\n[Parameter space ",r," of ",length(D)*length(W)*length(K),"]: Analyzing popgen data for ",D[d]," ; ",W[w]," ; ",K[k],": 0%   ", sep="")
      for (i in 1:n){ # looping over replicates
        d2 <- read.table(file=paste0("./output/output_individuals_run",i-1,".out"),header=T)
        # define patches
        patches <- formatC(d2$x, width=4, flag="0")
        # reformat diploid genetic data
        genotypes <- matrix(ncol=10,nrow=length(d2$x))
        for (g in 0:9){
          genotypes[,g+1] <- d2[,5+g*2]*100 + d2[,6+g*2]
        }
        # create genetic data dataframe
        genetics_data <- data.frame(patches,genotypes)
        rm(d2) # remove d2 to free memory
        
        gendata <- subset(genetics_data, genetics_data$patches %in% SS_Mod) # create subset of model output only containing nodes with empirical data as well
        gendata <- droplevels(gendata) # drop the uneccessary patch IDs
        
        missing <- which(is.na(match(SS_Mod,levels(genetics_data$patches)))) # which nodes are missing in the model output (happens if carrying capacity is varied)
        
        #*** Calculations ####
        #+++ Basic stats ####
        assign(paste0("calctime_",i-1),Sys.time())
        basic_time <- Sys.time()
        assign(paste0("basicstats_",i-1), basic.stats(gendata,diploid=T))
        basic_time <- Sys.time()-basic_time
        
        #+++ Allelic richness ####
        Ar_time <- Sys.time()
        assign(paste0("Ar_Mod_",i-1),allelic.richness(gendata, min.n=gen_thres))
        assign(paste0("meanAr_Mod_",i-1),apply(eval(parse(text=paste0("Ar_Mod_",i-1,"$Ar"))), 2, mean, na.rm=TRUE))
        if (length(missing)!=0){
          for (a in 1:length(missing)){
            assign(paste0("meanAr_Mod_",i-1),append(eval(parse(text=paste0("meanAr_Mod_",i-1))), NA, after=missing[a]-1))
          }
        }
        Ar_time <- Sys.time()-Ar_time
        
        #+++ Observed heterozygosity ####
        Ho_time <- Sys.time()
        assign(paste0("meanHo_Mod_",i-1),apply(eval(parse(text=paste0("basicstats_",i-1,"$Ho"))), 2, mean, na.rm=TRUE))
        if (length(missing)!=0){
          for (a in 1:length(missing)){
            assign(paste0("meanHo_Mod_",i-1),append(eval(parse(text=paste0("meanHo_Mod_",i-1))), NA, after=missing[a]-1))
          }
        }
        Ho_time <- Sys.time()-Ho_time
        
        assign(paste0("calctime_",i-1),Sys.time()-eval(parse(text=paste0("calctime_",i-1))))
        cat("\r[Parameter space ",r," of ",length(D)*length(W)*length(K),"]: Analyzing popgen data for ",D[d]," ; ",W[w]," ; ",K[k],": ",formatC(i/n*100,digits=3),"%   ", sep = "")
        
      } # end looping over replicates
      parameter_time <- Sys.time()-parameter_time
      cat("\n[Parameter space ",r," of ",length(D)*length(W)*length(K),"]: Analyzed popgen data for ",D[d]," ; ",W[w]," ; ",K[k],": Run time was ",formatC(parameter_time,digits=3)," ",attr(formatC(parameter_time,digits=3),"units"),".",sep="")
        
      meanAr_list <- mget(paste0("meanAr_Mod_", seq(0,n-1)))
      meanAr <- do.call('rbind', meanAr_list)
      meanAr_Mod <- apply(meanAr,2,mean,na.rm = TRUE)
        
      meanHo_list <- mget(paste0("meanHo_Mod_", seq(0,n-1)))
      meanHo <- do.call('rbind', meanHo_list)
      meanHo_Mod <- apply(meanHo,2,mean,na.rm = TRUE)
        
      save(parameter_time, meanAr, meanAr_Mod, meanHo, meanHo_Mod,
             calctime_0, basicstats_0, Ar_Mod_0, meanAr_Mod_0,
             calctime_1, basicstats_1, Ar_Mod_1, meanAr_Mod_1,
             calctime_2, basicstats_2, Ar_Mod_2, meanAr_Mod_2,
             calctime_3, basicstats_3, Ar_Mod_3, meanAr_Mod_3,
             calctime_4, basicstats_4, Ar_Mod_4, meanAr_Mod_4,
             calctime_5, basicstats_5, Ar_Mod_5, meanAr_Mod_5,
             calctime_6, basicstats_6, Ar_Mod_6, meanAr_Mod_6,
             calctime_7, basicstats_7, Ar_Mod_7, meanAr_Mod_7,
             calctime_8, basicstats_8, Ar_Mod_8, meanAr_Mod_8,
             calctime_9, basicstats_9, Ar_Mod_9, meanAr_Mod_9,
             file=paste0(WD,"/Output/IndPopGenData_",D[d],"_",W[w],"_",K[k],".Rdata"))
      
      cat("Done!\n")
      cat("+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+\n")
      
      setwd(WD)
    } # end looping over carrying capacities
    
  } # end looping over dispersal directionalities
  
} # end looping over dispersal rates

tot_time <- Sys.time()-tot_time # stop timer
cat("\n[All ",length(D)*length(W)*length(K)," parameter combinations]: Analyzing popgen data: Overall run time was ",formatC(tot_time,digits=3)," ",attr(formatC(tot_time,digits=3),"units"),".",sep="")

rm(list=c("a","i","g","d","w","k","r","D","W","K","D_label","W_label","K_label","missing","n","n_ind","gen_thres",
          "matrix","patches","gendata","genetics_data","genotypes","meanAr_list","meanHo_list","stand",
          "Ar_Mod_0", "Ar_Mod_1", "Ar_Mod_2", "Ar_Mod_3", "Ar_Mod_4", "Ar_Mod_5", "Ar_Mod_6", "Ar_Mod_7", "Ar_Mod_8", "Ar_Mod_9",
          "meanHo_Mod_0","meanHo_Mod_1","meanHo_Mod_2","meanHo_Mod_3","meanHo_Mod_4","meanHo_Mod_5","meanHo_Mod_6","meanHo_Mod_7","meanHo_Mod_8","meanHo_Mod_9",
          "basicstats_0", "basicstats_1", "basicstats_2", "basicstats_3", "basicstats_4", "basicstats_5", "basicstats_6", "basicstats_7", "basicstats_8", "basicstats_9",
          "calctime_0", "calctime_1", "calctime_2", "calctime_3", "calctime_4", "calctime_5" ,"calctime_6", "calctime_7", "calctime_8", "calctime_9",
          "meanAr", "meanAr_Mod", "meanAr_Mod_0", "meanAr_Mod_1", "meanAr_Mod_2", "meanAr_Mod_3", "meanAr_Mod_4", "meanAr_Mod_5", "meanAr_Mod_6", "meanAr_Mod_7", "meanAr_Mod_8", "meanAr_Mod_9", 
          "meanHo", "meanHo_Mod", "parameter_time","WD","tot_time","basic_time","Ar_time","Ho_time","SS_Mod","outlet_vertex"))

save.image("PopGenNet_data.RData")
