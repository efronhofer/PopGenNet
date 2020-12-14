# PopGenNet analysis
# 2020-12-10
# Author: Roman Alther, Eawag, Duebendorf, Switzerland
# Works in R ver. 3.6.1 and 4.0.3 (tested on GNU/Linux, MacOS, Windows) with RStudio ver. 1.3.1093
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

#### PREPARATION -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
rm(list=ls()) # start from scratch
existing_data=T # Should the analysis run with pre-collated and attached data from the authors? Otherwise the required raw data need to be prepared using the script 'PopGenNet_Prep_20201210.R', creating an output folder in '02_Data_prep'.
internal=F # defaults to FALSE, but was set to TRUE for publication preparation (lab internal use)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
WD <- getwd() # save current directory (should correspond to 03_Analysis, source script from there)

#* Load data ####
setwd("..")
DF <- getwd()
load("02_Data_prep/PopGenNet_data.RData")
setwd(WD)
if (internal){
  load("Gfos_data_20200925.RData")
}
output <- format(Sys.time(), "%Y%m%d")

##* Packages ####
# Check if custom-made packages OpenSwissRiverPlot and MultiPanel are already installed, else install from SWITCHdrive
if (!"OpenSwissRiverPlot" %in% installed.packages()[,"Package"]){
  ulib <- NULL # option to define path for user defined library -> change NULL to desired path
  source("https://drive.switch.ch/index.php/s/kgoAUIbqxYc92YP/download") # install OpenSwissRiverPlot
  # alternatively you may install from the included .tgz file
  # install.packages("OpenSwissRiverPlot.tgz", repos = NULL, type="source", lib=ulib, INSTALL_opts="--no-multiarch")
}
if (!"MultiPanel" %in% installed.packages()[,"Package"]){
  ulib <- NULL # option to define path for user defined library -> change NULL to desired path
  source("https://drive.switch.ch/index.php/s/tdaTpPUM7rH1P4X/download") # install MultiPanel
  # alternatively you may install from the included .tgz file
  # install.packages("MultiPanel.tgz", repos = NULL, type="source", lib=ulib, INSTALL_opts="--no-multiarch")
}

# Load CRAN packages
library(igraph) # to do network "stuff", works with version 1.2.4.2
library(jtools) # for "easy" model output (function summ()), works with version 2.1.1
library(rsq) # for R-squared (function rsq()), works with version 2.0

# Load custom packages and check for updates
if (internal){
  # Load lab internal package for publication figure preparation
  library(SwissRiverPlot) # to plot maps of Switzerland, works with version 0.4-2
  update_SRP()
}else{
  library(OpenSwissRiverPlot) # to plot maps of Switzerland, works with version 0.4-0
  update_OSRP()
}
library(MultiPanel) # to plot multipanel figures, works with version 0.6-3
update_MultiPanel()

#* Figures preparation ####

##*** Figure format ####
pdf=T # set to TRUE if figures should be prepared as PDF

##*** Size ####
fig.width=12 # standard figure width in inches
fig.height=9 # standard figure height in inches

##*** Smoothing span ####
loessspan <- 0.50 # the parameter ?? which controls the degree of smoothing in loess()

##*** Labels ####
label_A <- expression(italic(G.~fossarum)~"type A") # italic G. fossarum A
label_B <- expression(italic(G.~fossarum)~"type B") # italic G. fossarum B
label_mod <- "Model data"

D_label <- c("0.001","0.01","0.1")
W_label <- c("0","0.5","1")
K_label <- c("0","1")
D <- paste0("disp_",D_label)
W <- paste0("w_up_",W_label)
K <- paste0("K_scale_",K_label)
dlab <- c(as.expression(bquote(italic(d)~"="~.(D_label[[1]]))), # d = 000.1
          as.expression(bquote(italic(d)~"="~.(D_label[[2]]))), # d = 00.1
          as.expression(bquote(italic(d)~"="~.(D_label[[3]])))) # d = 0.1
wlab <- c(as.expression(bquote(italic(W)~"="~.(W_label[[1]]))), # W = 0.0
          as.expression(bquote(italic(W)~"="~.(W_label[[2]]))), # W = 0.5
          as.expression(bquote(italic(W)~"="~.(W_label[[3]])))) # W = 1.0
klab_short <- c(as.expression(bquote(italic(K)~"="~.(K_label[[1]]))), # K = 0
                as.expression(bquote(italic(K)~"="~.(K_label[[2]])))) # K = 1
klab <- c("Scaling: No","Scaling: Yes")
labs_short <- expand.grid(K_label, W_label, D_label)
# labs <- expand.grid(klab_short, wlab, dlab)
# labs_comb <- sprintf('%s; %s; %s', labs[,3], labs[,2], labs[,1])lab_Ar <- "Mean allelic richness"
labs_comb <- c(as.expression(bquote(italic(d)~"="~.(D_label[[1]])*";"~italic(W)~"="~.(W_label[[1]])*";"~italic(K)~"="~.(K_label[[1]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[1]])*";"~italic(W)~"="~.(W_label[[1]])*";"~italic(K)~"="~.(K_label[[2]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[1]])*";"~italic(W)~"="~.(W_label[[2]])*";"~italic(K)~"="~.(K_label[[1]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[1]])*";"~italic(W)~"="~.(W_label[[2]])*";"~italic(K)~"="~.(K_label[[2]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[1]])*";"~italic(W)~"="~.(W_label[[3]])*";"~italic(K)~"="~.(K_label[[1]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[1]])*";"~italic(W)~"="~.(W_label[[3]])*";"~italic(K)~"="~.(K_label[[2]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[2]])*";"~italic(W)~"="~.(W_label[[1]])*";"~italic(K)~"="~.(K_label[[1]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[2]])*";"~italic(W)~"="~.(W_label[[1]])*";"~italic(K)~"="~.(K_label[[2]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[2]])*";"~italic(W)~"="~.(W_label[[2]])*";"~italic(K)~"="~.(K_label[[1]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[2]])*";"~italic(W)~"="~.(W_label[[2]])*";"~italic(K)~"="~.(K_label[[2]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[2]])*";"~italic(W)~"="~.(W_label[[3]])*";"~italic(K)~"="~.(K_label[[1]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[2]])*";"~italic(W)~"="~.(W_label[[3]])*";"~italic(K)~"="~.(K_label[[2]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[3]])*";"~italic(W)~"="~.(W_label[[1]])*";"~italic(K)~"="~.(K_label[[1]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[3]])*";"~italic(W)~"="~.(W_label[[1]])*";"~italic(K)~"="~.(K_label[[2]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[3]])*";"~italic(W)~"="~.(W_label[[2]])*";"~italic(K)~"="~.(K_label[[1]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[3]])*";"~italic(W)~"="~.(W_label[[2]])*";"~italic(K)~"="~.(K_label[[2]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[3]])*";"~italic(W)~"="~.(W_label[[3]])*";"~italic(K)~"="~.(K_label[[1]]))),
               as.expression(bquote(italic(d)~"="~.(D_label[[3]])*";"~italic(W)~"="~.(W_label[[3]])*";"~italic(K)~"="~.(K_label[[2]]))))
lab_Ar <- "Mean allelic richness"
lab_Ho <- "Mean observed heterozygosity"
measure1 <- "SPO = Sum of perpendicular offsets"
measure1_short <- "SPO"
measure2 <- "MPO = Median of perpendicular offsets"
measure2_short <- "MPO"
measure3 <- "Perpendicular offset"
measure3a <- "perpendicular offsets"
measure4 <- "MDPO"
lab_sub <- c("(a)","(b)","(c)","(d)")
sum_digits <- 1
sum_cex <- 1.4
median_digits <- 3
median_cex <- 1.4

##*** Colors ####
# Colors for two species
col_Gfos_A <- "#bf812d" # RGB 191, 129, 45; yellowish or orange
col_Gfos_B <- "#35978f" # RGB 53, 151, 143; turquoise
# Color for Gammarus fossarum complex
col_Gfos <- "#9970ab"
# Color for model data
colMod <- rgb(70,70,70,max=255)
# Colors for waterways
col_water <- "darkgrey"
col_rhine <- "lightgrey"
alpha1 <- "CC" # 80% alpha, check https://gist.github.com/lopspower/03fb1cc0ac9f32ef38f4 for transparency code
alpha2 <- "80" # 50% alpha, check https://gist.github.com/lopspower/03fb1cc0ac9f32ef38f4 for transparency code
white_transparent <- paste0("#FFFFFF",alpha2) # Transparent white
# Preparation for heatmap style color palette
col_pal <- c('dark red','white','navy blue')
col_fun <- colorRamp(col_pal)
rgb2hex <- function(r,g,b) rgb(r, g, b, maxColorValue = 255)
paletteLength <- 100 # how finely should the color ramp be divided
my_palette <- colorRampPalette(col_pal)(n = paletteLength)
# Color palette for Local Polynomial Regression figures
loess_pal_short <- c("#54278f","#807dba","#bcbddc")
loess_pal <- rep(loess_pal_short, each=length(W)*length(K))
loess_lty <- rep(c(1,1,2,2,3,3), length(D))
loess_lwd <- rep(c(2,4), length(D)*length(W))

#* Functions ####

#### Function for GLM and prediction appending to data
# The function requires a model formulation, the data, a name for the output, the model family, a significance level, and if selection should be implemented
glm.bind <- function(model, data, name, family, sign, step=F){
  assign(paste0("glm_",name),glm(formula(model), data, family=family))
  b <- get(paste0("glm_",name))
  if (step){
    assign(paste0("sglm_",name), step(get(paste0("glm_",name)))) # backward selection
    c <- get(paste0("sglm_",name))
    assign(paste0("Predict_sglm_",name),predict(get(paste0("sglm_",name)), type="response", se.fit=T))
    assign(paste0("sglm_",name,"_upr"),get(paste0("Predict_sglm_",name))$fit + (sign * get(paste0("Predict_sglm_",name))$se.fit))
    assign(paste0("sglm_",name,"_lwr"),get(paste0("Predict_sglm_",name))$fit - (sign * get(paste0("Predict_sglm_",name))$se.fit))
    assign(paste0("sglm_",name,"_fit"),get(paste0("Predict_sglm_",name))$fit)
    data <- cbind(data, cbind(get(paste0("sglm_",name,"_fit")), get(paste0("sglm_",name,"_lwr")), get(paste0("sglm_",name,"_upr"))))
    colnames(data)[c(ncol(data)-2,ncol(data)-1,ncol(data))] <- c(paste0("sglm_",name,"_fit"),paste0("sglm_",name,"_lwr"),paste0("sglm_",name,"_upr"))
  }else{
    assign(paste0("Predict_glm_",name),predict(get(paste0("glm_",name)), type="response", se.fit=T))
    assign(paste0("glm_",name,"_upr"),get(paste0("Predict_glm_",name))$fit + (sign * get(paste0("Predict_glm_",name))$se.fit))
    assign(paste0("glm_",name,"_lwr"),get(paste0("Predict_glm_",name))$fit - (sign * get(paste0("Predict_glm_",name))$se.fit))
    assign(paste0("glm_",name,"_fit"),get(paste0("Predict_glm_",name))$fit)
    data <- cbind(data, cbind(get(paste0("glm_",name,"_fit")), get(paste0("glm_",name,"_lwr")), get(paste0("glm_",name,"_upr"))))
    colnames(data)[c(ncol(data)-2,ncol(data)-1,ncol(data))] <- c(paste0("glm_",name,"_fit"),paste0("glm_",name,"_lwr"),paste0("glm_",name,"_upr"))
  }
  a <- data
  if (step){
    return(list(a,b,c))
  }else{
    return(list(a,b))
  }
}

#### Function for GLM plots
# wrapper function to prepare the figures of the GLM output, including confidence intervals (as defined in glm.bind(sign=x))
GLMplot <- function(x,y,dat,model,xlabel,ylabel,axislog="", col1=col_Gfos_A, col2=col_Gfos_B, pt.cex=1, CI_border=T, trans=0.3, xrev=F, xax=NULL, pointtrans=F, cex.lab=2, cex.axis=1.5, legend=T, cex.legend=1.5){
  col1trans <- rgb(col2rgb(col1)[1,]/255,col2rgb(col1)[2,]/255,col2rgb(col1)[3,]/255,trans)
  col2trans <- rgb(col2rgb(col2)[1,]/255,col2rgb(col2)[2,]/255,col2rgb(col2)[3,]/255,trans)
  form <- reformulate(y, response = x)
  formfit <- reformulate(y, response = paste0(model,"_fit"))
  formupr <- reformulate(y, response = paste0(model,"_upr"))
  formlwr <- reformulate(y, response = paste0(model,"_lwr"))
  xcol <- which(colnames(dat)==x)
  ycol <- which(colnames(dat)==y)
  lwrcol <- which(colnames(dat)==paste0(model,"_lwr"))
  uprcol <- which(colnames(dat)==paste0(model,"_upr"))
  DATAordered <- dat[order(dat[,ycol]),]
  left <- min(DATAordered[,ycol])
  right <- max(DATAordered[,ycol])
  if (xrev==T){
    xrange <- c(right,left)
  }else{
    xrange <- c(left,right)
  }
  if (pointtrans==T){
    col1point <- col1trans
    col2point <- col2trans
  }else{
    col1point <- col1
    col2point <- col2
  }
  par(mar=c(3.1+cex.lab, 3.1+cex.lab, 0.5, 0.5))
  plot(form, dat, type = "n", las = 1, bty = "l",
       xlab=xlabel,
       ylab=ylabel,
       xlim=xrange,
       log=axislog,
       xaxt=xax, cex.lab=cex.lab, cex.axis=cex.axis)
  polygon(c(rev(DATAordered[,ycol][DATAordered$spec=="A"]), DATAordered[,ycol][DATAordered$spec=="A"]),
          c(rev(DATAordered[,lwrcol][DATAordered$spec=="A"]), DATAordered[,uprcol][DATAordered$spec=="A"]),
          col = col1trans, border = NA)
  polygon(c(rev(DATAordered[,ycol][DATAordered$spec=="B"]), DATAordered[,ycol][DATAordered$spec=="B"]),
          c(rev(DATAordered[,lwrcol][DATAordered$spec=="B"]), DATAordered[,uprcol][DATAordered$spec=="B"]),
          col = col2trans, border = NA)
  points(form, data = subset(DATAordered, spec == "A"), pch = 16, col = col1point, cex=pt.cex)
  points(form, data = subset(DATAordered, spec == "B"), pch = 16, col = col2point, cex=pt.cex)
  lines(formfit, data = subset(DATAordered, spec == "A"), lwd = 2.5, col=col1)
  if(CI_border){lines(formupr, data = subset(DATAordered, spec == "A"), lwd = 2, lty=2, col=col1)}
  if(CI_border){lines(formlwr, data = subset(DATAordered, spec == "A"), lwd = 2, lty=2, col=col1)}
  lines(formfit, data = subset(DATAordered, spec == "B"), lwd = 2.5, col=col2)
  if(CI_border){lines(formupr, data = subset(DATAordered, spec == "B"), lwd = 2, lty=2, col=col2)}
  if(CI_border){lines(formlwr, data = subset(DATAordered, spec == "B"), lwd = 2, lty=2, col=col2)}
  if (legend){
    legend("bottomright",c(label_A,label_B),pch = 16, col = c(col1,col2), bty="n", cex=cex.legend)
  }
}

#### Function to calculate Euclidean distance
# As proposed by user Shambho (https://stackoverflow.com/users/3547167) here: https://stackoverflow.com/a/24747155
euc.dist <- function(x1, x2){sqrt(sum((x1 - x2) ^ 2))}

#### Function to find endpoint for a perpendicular segment from the point (x0,y0) to the line
# As proposed by user MrFlick (https://stackoverflow.com/users/2372064) here: https://stackoverflow.com/a/30399576
perp.segment.coord <- function(x0, y0, a=0,b=1){
  # defined by lm.mod as y=a+b*x
  x1 <- (x0+b*y0-a*b)/(1+b^2)
  y1 <- a + b*x1
  return(list(x0=x0, y0=y0, x1=x1, y1=y1))
}

#### DATA PREPARATION -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
updist_Mod <- updist[match(modsite,V(net)$name)]
catch_Mod <- V(net)$Total_Catch[match(modsite,V(net)$name)]
catch_A_red <- V(net)$Total_Catch[match(microsite_A_red,V(net)$name)] # prepare for Gfos A as well (not preexisting as such)
catch_B <- V(net)$Total_Catch[match(microsite_B,V(net)$name)] # prepare for Gfos B as well (not preexisting as such)
betw_undir_Mod <- betw_undir[match(modsite,V(net)$name)]
clos_undir_Mod <- clos_undir[match(modsite,V(net)$name)]
degree_undir_Mod <- degree_undir$res[match(modsite,V(net)$name)]

if(!existing_data){
  ##*** PopGen table preparation ####
  Ar_modelled <- data.frame(matrix(nrow=length(empiricaldata), ncol=length(D)*length(W)*length(K)+1))
  rownames(Ar_modelled) <- adj.mat.names[empiricaldata]
  Ar_modelled_smooth <- data.frame(matrix(nrow=length(empiricaldata), ncol=length(D)*length(W)*length(K)+1))
  rownames(Ar_modelled_smooth) <- adj.mat.names[empiricaldata]
  Ho_modelled <- data.frame(matrix(nrow=length(empiricaldata), ncol=length(D)*length(W)*length(K)+1))
  rownames(Ho_modelled) <- adj.mat.names[empiricaldata]
  Ho_modelled_smooth <- data.frame(matrix(nrow=length(empiricaldata), ncol=length(D)*length(W)*length(K)+1))
  rownames(Ho_modelled_smooth) <- adj.mat.names[empiricaldata]
  
  r <- 0
  for (d in 1:length(D)){ # looping over dispersal rates 
    for (w in 1:length(W)){ # looping over dispersal directionalities
      for (k in 1:length(K)){ # looping over carrying capacities
        r <- r+1
        
        label_Mod <- paste0("D=",D_label[d],", W_up=",W_label[w],", K=",K_label[k])
        
        load(paste0(DF,"/02_Data_prep/Output/IndPopGenData_",D[d],"_",W[w],"_",K[k],".Rdata"))
        
        orderMod <- order(updist_Mod, decreasing = T)
        
        if (r==1){
          match_Mod <- match_Mod[orderMod]
        }
        
        #*** Mean Ar by upstream distance ####
        updist_Mod_plot <- updist_Mod[orderMod]
        meanAr_Mod_updist <- meanAr_Mod[orderMod]
        loessAr_Mod_updist <- loess(meanAr_Mod_updist~updist_Mod_plot, span=loessspan)
        smoothedAr_Mod_updist <- predict(loessAr_Mod_updist)
        orderA_red <- order(updist_A_red, decreasing = T)
        updist_A_red_plot <- updist_A_red[orderA_red]
        meanAr_A_red_updist <- meanAr_A_red[orderA_red]
        loessAr_A_red_updist <- loess(meanAr_A_red_updist~updist_A_red_plot, span=loessspan)
        smoothedAr_A_red_updist <- predict(loessAr_A_red_updist)
        orderB <- order(updist_B, decreasing = T)
        updist_B_plot <- updist_B[orderB]
        meanAr_B_updist <- meanAr_B[orderB]
        loessAr_B_updist <- loess(meanAr_B_updist~updist_B_plot, span=loessspan)
        smoothedAr_B_updist <- predict(loessAr_B_updist)
        
        #*** Mean Ho by upstream distance ####
        updist_Mod_plot <- updist_Mod[orderMod]
        meanHo_Mod_updist <- meanHo_Mod[orderMod]
        loessHo_Mod_updist <- loess(meanHo_Mod_updist~updist_Mod_plot, span=loessspan)
        smoothedHo_Mod_updist <- predict(loessHo_Mod_updist)
        orderA_red <- order(updist_A_red, decreasing = T)
        updist_A_red_plot <- updist_A_red[orderA_red]
        meanHo_A_red_updist <- meanHo_A_red[orderA_red]
        loessHo_A_red_updist <- loess(meanHo_A_red_updist~updist_A_red_plot, span=loessspan)
        smoothedHo_A_red_updist <- predict(loessHo_A_red_updist)
        orderB <- order(updist_B, decreasing = T)
        updist_B_plot <- updist_B[orderB]
        meanHo_B_updist <- meanHo_B[orderB]
        loessHo_B_updist <- loess(meanHo_B_updist~updist_B_plot, span=loessspan)
        smoothedHo_B_updist <- predict(loessHo_B_updist)
        
        # save modelled popgen values to table
        if (r==1){
          Ar_modelled[,1] <- updist_Mod_plot
          colnames(Ar_modelled)[1] <- "upstream_distance"
        }
        Ar_modelled[,1+r] <- meanAr_Mod_updist
        colnames(Ar_modelled)[r+1] <- paste0(D[d],"_",W[w],"_",K[k])
        
        # save smoothed modelled popgen values to table
        if (r==1){
          Ar_modelled_smooth[,1] <- updist_Mod_plot
          colnames(Ar_modelled_smooth)[1] <- "upstream_distance"
        }
        Ar_modelled_smooth[,1+r] <- smoothedAr_Mod_updist
        colnames(Ar_modelled_smooth)[r+1] <- paste0(D[d],"_",W[w],"_",K[k])
        
        # save modelled popgen values to table
        if (r==1){
          Ho_modelled[,1] <- updist_Mod_plot
          colnames(Ho_modelled)[1] <- "upstream_distance"
        }
        Ho_modelled[,1+r] <- meanHo_Mod_updist
        colnames(Ho_modelled)[r+1] <- paste0(D[d],"_",W[w],"_",K[k])
        
        # save smoothed modelled popgen values to table
        if (r==1){
          Ho_modelled_smooth[,1] <- updist_Mod_plot
          colnames(Ho_modelled_smooth)[1] <- "upstream_distance"
        }
        Ho_modelled_smooth[,1+r] <- smoothedHo_Mod_updist
        colnames(Ho_modelled_smooth)[r+1] <- paste0(D[d],"_",W[w],"_",K[k])
        
        # setwd(WD)
        
      } # end looping over carrying capacities
    } # end looping over dispersal directionalities
  } # end looping over dispersal rates
}else{
  load("PopGenData.RData")
}

#### CORRELATION BETWEEN EXPLANATORY VARS -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
expl.var1 <- DATA[,!(colnames(DATA) %in% c("meanAr","meanHo","network_match","spec","log_updist","log_catch"))]
cor.var1 <- cor(expl.var1, method = "kendall")

expl.var <- DATA[,!(colnames(DATA) %in% c("meanAr","meanHo","degree","betw_dir","clos_dir","network_match","spec","log_updist","log_catch","clos_undir","std_clos_dir"))]
cor.var <- cor(expl.var, method = "kendall")

#### MODELLING EMPIRICAL DATA -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
##*** Model settings ####
#### Define ranges for predict function
range_betw_undir_B <- seq(min(betw_undir_B),max(betw_undir_B),(max(betw_undir_B)-min(betw_undir_B))/1000)
range_betw_undir_A_red <- seq(min(betw_undir_A_red),max(betw_undir_A_red),(max(betw_undir_A_red)-min(betw_undir_A_red))/1000)
range_clos_undir_B <- seq(min(clos_undir_B),max(clos_undir_B),(max(clos_undir_B)-min(clos_undir_B))/1000)
range_clos_undir_A_red <- seq(min(clos_undir_A_red),max(clos_undir_A_red),(max(clos_undir_A_red)-min(clos_undir_A_red))/1000)
range_degree_undir_B <- seq(min(degree_undir_B),max(degree_undir_B),(max(degree_undir_B)-min(degree_undir_B))/1000)
range_degree_undir_A_red <- seq(min(degree_undir_A_red),max(degree_undir_A_red),(max(degree_undir_A_red)-min(degree_undir_A_red))/1000)
range_betw_dir_B <- seq(min(betw_dir_B),max(betw_dir_B),(max(betw_dir_B)-min(betw_dir_B))/1000)
range_betw_dir_A_red <- seq(min(betw_dir_A_red),max(betw_dir_A_red),(max(betw_dir_A_red)-min(betw_dir_A_red))/1000)
range_clos_dir_B <- seq(min(clos_dir_B),max(clos_dir_B),(max(clos_dir_B)-min(clos_dir_B))/1000)
range_clos_dir_A_red <- seq(min(clos_dir_A_red),max(clos_dir_A_red),(max(clos_dir_A_red)-min(clos_dir_A_red))/1000)
range_updist_B <- seq(1,max(updist_B),1000)
range_updist_A_red <- seq(1,max(updist_A_red),1000)
range_dist_B <- seq(1,max(dist_B),1000)
range_dist_A_red <- seq(1,max(dist_A_red),1000)

 ##*** Models meanAr ####
#### Response type and confidence interval
meanAr_family <- "Gamma"
critval <- 1.96 ## approx 95% CI
sigval <- 0.05 ## significance value
test.stat <- "F"

#### Initial model exploration
glm(meanAr ~ updist*betw_dir*std_clos_undir*catch*spec, DATA, family=meanAr_family)
step(glm(meanAr ~ updist*betw_dir*std_clos_undir*catch*spec, DATA, family=meanAr_family))
glm(meanAr ~ updist+betw_dir+std_clos_undir+catch+spec, DATA, family=meanAr_family)
step(glm(meanAr ~ updist+betw_dir+std_clos_undir+catch+spec, DATA, family=meanAr_family))

#### Full GLM of allelic richness without interactions
model <- glm.bind(meanAr ~ updist+betw_dir+std_clos_undir+catch+spec, DATA, "Ar_full", meanAr_family, critval, step=T)
DATA <- model[[1]]
glm_Ar_full <- model[[2]]
summary(glm_Ar_full)
sglm_Ar_full <- model[[3]]
summary(sglm_Ar_full)
anova(sglm_Ar_full, test=test.stat)
# Output for paper
model <- glm(meanAr ~ updist+catch+spec, DATA, family=meanAr_family)
null <- glm(formula = meanAr ~ 1, family = meanAr_family, data = DATA)
attr(summ(model),"rsq")
rsq(model, type="v")
anova(null, model, test="F") # F test statistics for Gamma
anova(model, test=test.stat)
df.residual(model)

#### GLM of allelic richness by upstream distance * species
model <- glm.bind(meanAr ~ updist+spec, DATA, "Ar_updist", meanAr_family, critval, step=T)
DATA <- model[[1]]
glm_Ar_updist <- model[[2]]
summary(glm_Ar_updist)
sglm_Ar_updist <- model[[3]]
summary(sglm_Ar_updist)
# Output for paper
model <- glm(formula = meanAr ~ updist + spec, family = meanAr_family, data = DATA)
attr(summ(model),"rsq")
rsq(model, type="v")
anova(null, model, test="F") # F test statistics for Gamma
anova(model, test=test.stat) # F test statistics for Gamma
df.residual(model)

#### GLM of allelic richness by catchment size * species
model <- glm.bind(meanAr ~ catch+spec, DATA, "Ar_catch", meanAr_family, critval, step=T)
DATA <- model[[1]]
glm_Ar_catch <- model[[2]]
summary(glm_Ar_catch)
sglm_Ar_catch <- model[[3]]
summary(sglm_Ar_catch)
# Output for paper
model <- glm(formula = meanAr ~ catch, family = meanAr_family, data = DATA)
attr(summ(model),"rsq")
rsq(model, type="v")
anova(null, model, test="F") # F test statistics for Gamma
anova(model, test=test.stat) # F test statistics for Gamma
df.residual(model)

#### GLM of allelic richness by betweeness centrality * species
model <- glm.bind(meanAr ~ betw_undir+spec, DATA, "Ar_betw_undir", meanAr_family, critval, step=T)
DATA <- model[[1]]
glm_Ar_betw_undir <- model[[2]]
summary(glm_Ar_betw_undir)
sglm_Ar_betw_undir <- model[[3]]
summary(sglm_Ar_betw_undir)
# Output for paper
model <- glm(formula = meanAr ~ betw_undir, family = meanAr_family, data = DATA)
attr(summ(model),"rsq")
rsq(model, type="v")
anova(model, test=test.stat) # F test statistics for Gamma
df.residual(model)

#### GLM of allelic richness by undirected closeness centrality * species
model <- glm.bind(meanAr ~ std_clos_undir+spec, DATA, "Ar_clos", meanAr_family, critval, step=T)
DATA <- model[[1]]
glm_Ar_clos <- model[[2]]
summary(glm_Ar_clos)
sglm_Ar_clos <- model[[3]]
summary(sglm_Ar_clos)
# Output for paper
model <- glm(formula = meanAr ~ std_clos_undir, family = meanAr_family, data = DATA)
attr(summ(model),"rsq")
rsq(model, type="v")
anova(null, model, test=test.stat) # F test statistics for Gamma
anova(model, test=test.stat) # F test statistics for Gamma
df.residual(model)

##*** Models meanHo ####
#### Response type and confidence interval
meanHo_family <- "quasibinomial"
critval <- 1.96 ## approx 95% CI
test.stat <- "F"

#### Initial model exploration
int.mod <- glm(meanHo ~ updist*betw_dir*std_clos_undir*catch*spec, family=meanHo_family, DATA)
summary(int.mod)
lin.mod <- glm(meanHo ~ updist+betw_dir+std_clos_undir+catch+spec, family=meanHo_family, DATA)
summary(lin.mod)

M1 <- lin.mod
drop1(M1, test=test.stat)
M2 <- update(M1,~.-updist)
drop1(M2, test=test.stat)
M3 <- update(M2,~.-catch)
drop1(M3, test=test.stat)
summary(M3)

#### Combined GLM of observed heterozygosity
model <- glm.bind(meanHo ~ updist+betw_dir+std_clos_undir+catch+spec, DATA, "Ho_full", meanHo_family, critval)
DATA <- model[[1]]
glm_Ho_full <- model[[2]]
summary(glm_Ho_full)
# Output for paper
null <- glm(formula = meanHo ~ 1, data=DATA, family=meanHo_family)
model <- glm(meanHo ~ updist+betw_dir+std_clos_undir+catch+spec, family=meanHo_family, DATA)
attr(summ(model),"rsq")
anova(null, model, test=test.stat) # F test statistics for Quasibinomial

#### Combined GLM of observed heterozygosity
model <- glm.bind(meanHo ~ betw_dir+std_clos_undir+spec, DATA, "Ho_sel", meanHo_family, critval)
DATA <- model[[1]]
glm_Ho_sel <- model[[2]]
summary(glm_Ho_sel)
# Output for paper
model <- glm(meanHo ~ betw_dir+std_clos_undir+spec, DATA, family= meanHo_family)
rsq(model, type="v")
anova(null, model, test=test.stat) # F test statistics for Quasibinomial

#### Combined GLM of observed heterozygosity by closeness centrality * species
model <- glm.bind(meanHo ~ std_clos_undir+spec, DATA, "Ho_clos", meanHo_family, critval)
DATA <- model[[1]]
glm_Ho_clos <- model[[2]]
summary(glm_Ho_clos)
# Output for paper
model <- glm(meanHo ~ betw_dir+std_clos_undir+spec, DATA, family= meanHo_family)
rsq(model, type="v")
anova(null, model, test=test.stat) # F test statistics for Quasibinomial

#### Combined GLM of observed heterozygosity by upstream distance * species
model <- glm.bind(meanHo ~ updist+spec, DATA, "Ho_updist", meanHo_family, critval)
DATA <- model[[1]]
glm_Ho_updist <- model[[2]]
summary(glm_Ho_updist)
# Output for paper
model <- glm(meanHo ~ updist+spec, DATA, family=meanHo_family)
rsq(model, type="v")
anova(null, model, test=test.stat) # F test statistics for Quasibinomial

#### Combined GLM of observed heterozygosity by betweenness centrality * species
model <- glm.bind(meanHo ~ betw_dir+spec, DATA, "Ho_betw", meanHo_family, critval)
DATA <- model[[1]]
glm_Ho_betw <- model[[2]]
summary(glm_Ho_betw)
# Output for paper
model <- glm(meanHo ~ betw_dir+spec, DATA, family=meanHo_family)
rsq(model, type="v")
anova(null, model, test=test.stat) # F test statistics for Quasibinomial

#### Combined GLM of observed heterozygosity by catchment size * species
model <- glm.bind(meanHo ~ catch+spec, DATA, "Ho_catch", meanHo_family, critval)
DATA <- model[[1]]
glm_Ho_catch <- model[[2]]
summary(model[[2]])
# Output for paper
model <- glm(meanHo ~ catch+spec, DATA, family=meanHo_family)
rsq(model, type="v")
anova(null, model, test=test.stat) # F test statistics for Quasibinomial

#### Combined GLM of observed heterozygosity by upstream distance * clos * betw_dir * species
model <- glm.bind(meanHo ~ (updist+std_clos_undir+betw_dir+spec)^2, DATA, "Ho_ext", meanHo_family, critval)
DATA <- model[[1]]
glm_Ho_ext <- model[[2]]
summary(glm_Ho_ext)
# Output for paper
model <- glm(meanHo ~ (updist+std_clos_undir+betw_dir+spec)^2, DATA, family=meanHo_family)
rsq(model, type="v")
anova(null, model, test=test.stat) # F test statistics for Quasibinomial

#### PERPENDICULAR DISTANCES -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
##*** Preparing matrices####
orthdist_Ar_A <- vector()
orthdist_Ar_A_directed <- vector()
sum_orthdist_Ar_A <- vector()
hist_Ar_A <- matrix(nrow=length(modsite_GfosA), ncol=ncol(Ar_modelled)-1)
hist_Ar_A_directed <- matrix(nrow=length(modsite_GfosA), ncol=ncol(Ar_modelled)-1)

orthdist_Ar_B <- vector()
orthdist_Ar_B_directed <- vector()
sum_orthdist_Ar_B <- vector()
hist_Ar_B <- matrix(nrow=length(modsite_GfosB), ncol=ncol(Ar_modelled)-1)
hist_Ar_B_directed <- matrix(nrow=length(modsite_GfosB), ncol=ncol(Ar_modelled)-1)

orthdist_Ho_A <- vector()
orthdist_Ho_A_directed <- vector()
sum_orthdist_Ho_A <- vector()
hist_Ho_A <- matrix(nrow=length(modsite_GfosA), ncol=ncol(Ho_modelled)-1)
hist_Ho_A_directed <- matrix(nrow=length(modsite_GfosA), ncol=ncol(Ho_modelled)-1)

orthdist_Ho_B <- vector()
orthdist_Ho_B_directed <- vector()
sum_orthdist_Ho_B <- vector()
hist_Ho_B <- matrix(nrow=length(modsite_GfosB), ncol=ncol(Ho_modelled)-1)
hist_Ho_B_directed <- matrix(nrow=length(modsite_GfosB), ncol=ncol(Ho_modelled)-1)

##*** Calculating distances####
for (i in 2:ncol(Ar_modelled)){
  Ar_Mod_A <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosA]
  Ar_Mod_B <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosB]
  Ho_Mod_A <- Ho_modelled[,i][rownames(Ho_modelled)%in%modsite_GfosA]
  Ho_Mod_B <- Ho_modelled[,i][rownames(Ho_modelled)%in%modsite_GfosB]
  
  for (j in 1:length(meanAr_A_red_updist)){
    point <- cbind(meanAr_A_red_updist,Ar_Mod_A)[j,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    orthdist_Ar_A[j] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
    orthdist_Ar_A_directed[j] <- sign(seg[2]-seg[4])*orthdist_Ar_A[j]
  }
  for (k in 1:length(meanAr_B_updist)){
    point <- cbind(meanAr_B_updist,Ar_Mod_B)[k,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    orthdist_Ar_B[k] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
    orthdist_Ar_B_directed[k] <- sign(seg[2]-seg[4])*orthdist_Ar_B[k]
  }
  for (j in 1:length(meanHo_A_red_updist)){
    point <- cbind(meanHo_A_red_updist,Ho_Mod_A)[j,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    orthdist_Ho_A[j] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
    orthdist_Ho_A_directed[j] <- sign(seg[2]-seg[4])*orthdist_Ho_A[j]
  }
  for (k in 1:length(meanHo_B_updist)){
    point <- cbind(meanHo_B_updist,Ho_Mod_B)[k,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    orthdist_Ho_B[k] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
    orthdist_Ho_B_directed[k] <- sign(seg[2]-seg[4])*orthdist_Ho_B[k]
  }
  sum_orthdist_Ar_A[[i-1]] <- sum(orthdist_Ar_A)
  sum_orthdist_Ar_B[[i-1]] <- sum(orthdist_Ar_B)
  sum_orthdist_Ho_A[[i-1]] <- sum(orthdist_Ho_A)
  sum_orthdist_Ho_B[[i-1]] <- sum(orthdist_Ho_B)
  
  hist_Ar_A[,i-1] <- orthdist_Ar_A
  hist_Ar_B[,i-1] <- orthdist_Ar_B
  hist_Ho_A[,i-1] <- orthdist_Ho_A
  hist_Ho_B[,i-1] <- orthdist_Ho_B
  
  hist_Ar_A_directed[,i-1] <- orthdist_Ar_A_directed
  hist_Ar_B_directed[,i-1] <- orthdist_Ar_B_directed
  hist_Ho_A_directed[,i-1] <- orthdist_Ho_A_directed
  hist_Ho_B_directed[,i-1] <- orthdist_Ho_B_directed
}
names(sum_orthdist_Ar_A) <- colnames(Ar_modelled)[-1]
names(sum_orthdist_Ar_B) <- colnames(Ar_modelled)[-1]
names(sum_orthdist_Ho_A) <- colnames(Ho_modelled)[-1]
names(sum_orthdist_Ho_B) <- colnames(Ho_modelled)[-1]

#### CORRELATION SMOOTHED VALUES -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
##*** Correlation test Ar####
cortest_Ar <- data.frame(parameter=character(), d=integer(), w=integer(), k=integer(), tau_A=integer(), pval_A=integer(), tau_B=integer(), pval_B=integer(), row.names = 1)
for (i in 2:ncol(Ar_modelled_smooth)){
  smoothAr_Mod_A <- Ar_modelled_smooth[,i][rownames(Ar_modelled_smooth)%in%modsite_GfosA]
  test_A <- cor.test(smoothAr_Mod_A, smoothedAr_A_red_updist, method = "kendall")
  cortest_Ar[i-1,4] <- test_A$estimate
  cortest_Ar[i-1,5] <- test_A$p.value
  smoothAr_Mod_B <- Ar_modelled_smooth[,i][rownames(Ar_modelled_smooth)%in%modsite_GfosB]
  test_B <- cor.test(smoothAr_Mod_B, smoothedAr_B_updist, method = "kendall")
  cortest_Ar[i-1,6] <- test_B$estimate
  cortest_Ar[i-1,7] <- test_B$p.value
}
row.names(cortest_Ar) <- colnames(Ar_modelled_smooth)[-1]
cortest_Ar$d <- rep(c(0.001,0.01,0.1),each=6)
cortest_Ar$w <- rep(c(0,0.5,1),3,each=2)
cortest_Ar$k <- rep(c(0,1),9)

# Matrices for GfosA
cortest_Ar_k1 <- cortest_Ar[cortest_Ar$k==1,]
cortest_Ar_k1_matrix_A <- list()
cortest_Ar_k1_matrix_A[[1]] <- cortest_Ar_k1$tau_A[cortest_Ar_k1$d==0.001]
cortest_Ar_k1_matrix_A[[2]] <- cortest_Ar_k1$tau_A[cortest_Ar_k1$d==0.01]
cortest_Ar_k1_matrix_A[[3]] <- cortest_Ar_k1$tau_A[cortest_Ar_k1$d==0.1]
cortest_Ar_k1_matrix_A = do.call(cbind, cortest_Ar_k1_matrix_A)
colnames(cortest_Ar_k1_matrix_A) <- c("d 0.001","d 0.01","d 0.1")
row.names(cortest_Ar_k1_matrix_A) <- c("w 0","w 0.5","w 1")
cortest_Ar_k0 <- cortest_Ar[cortest_Ar$k==0,]
cortest_Ar_k0_matrix_A <- list()
cortest_Ar_k0_matrix_A[[1]] <- cortest_Ar_k0$tau_A[cortest_Ar_k0$d==0.001]
cortest_Ar_k0_matrix_A[[2]] <- cortest_Ar_k0$tau_A[cortest_Ar_k0$d==0.01]
cortest_Ar_k0_matrix_A[[3]] <- cortest_Ar_k0$tau_A[cortest_Ar_k0$d==0.1]
cortest_Ar_k0_matrix_A = do.call(cbind, cortest_Ar_k0_matrix_A)
colnames(cortest_Ar_k0_matrix_A) <- c("d 0.001","d 0.01","d 0.1")
row.names(cortest_Ar_k0_matrix_A) <- c("w 0","w 0.5","w 1")

# Matrices for GfosB
cortest_Ar_k1_matrix_B <- list()
cortest_Ar_k1_matrix_B[[1]] <- cortest_Ar_k1$tau_B[cortest_Ar_k1$d==0.001]
cortest_Ar_k1_matrix_B[[2]] <- cortest_Ar_k1$tau_B[cortest_Ar_k1$d==0.01]
cortest_Ar_k1_matrix_B[[3]] <- cortest_Ar_k1$tau_B[cortest_Ar_k1$d==0.1]
cortest_Ar_k1_matrix_B = do.call(cbind, cortest_Ar_k1_matrix_B)
colnames(cortest_Ar_k1_matrix_B) <- c("d 0.001","d 0.01","d 0.1")
row.names(cortest_Ar_k1_matrix_B) <- c("w 0","w 0.5","w 1")
cortest_Ar_k0_matrix_B <- list()
cortest_Ar_k0_matrix_B[[1]] <- cortest_Ar_k0$tau_B[cortest_Ar_k0$d==0.001]
cortest_Ar_k0_matrix_B[[2]] <- cortest_Ar_k0$tau_B[cortest_Ar_k0$d==0.01]
cortest_Ar_k0_matrix_B[[3]] <- cortest_Ar_k0$tau_B[cortest_Ar_k0$d==0.1]
cortest_Ar_k0_matrix_B = do.call(cbind, cortest_Ar_k0_matrix_B)
colnames(cortest_Ar_k0_matrix_B) <- c("d 0.001","d 0.01","d 0.1")
row.names(cortest_Ar_k0_matrix_B) <- c("w 0","w 0.5","w 1")

##*** Correlation test Ho####
cortest_Ho <- data.frame(parameter=character(), d=integer(), w=integer(), k=integer(), tau_A=integer(), pval_A=integer(), tau_B=integer(), pval_B=integer(), row.names = 1)
for (i in 2:ncol(Ho_modelled_smooth)){
  smoothAr_Mod_A <- Ho_modelled_smooth[,i][rownames(Ho_modelled_smooth)%in%modsite_GfosA]
  test_A <- cor.test(smoothAr_Mod_A, smoothedAr_A_red_updist, method = "kendall")
  cortest_Ho[i-1,4] <- test_A$estimate
  cortest_Ho[i-1,5] <- test_A$p.value
  smoothAr_Mod_B <- Ho_modelled_smooth[,i][rownames(Ho_modelled_smooth)%in%modsite_GfosB]
  test_B <- cor.test(smoothAr_Mod_B, smoothedAr_B_updist, method = "kendall")
  cortest_Ho[i-1,6] <- test_B$estimate
  cortest_Ho[i-1,7] <- test_B$p.value
}
row.names(cortest_Ho) <- colnames(Ho_modelled_smooth)[-1]
cortest_Ho$d <- rep(c(0.001,0.01,0.1),each=6)
cortest_Ho$w <- rep(c(0,0.5,1),3,each=2)
cortest_Ho$k <- rep(c(0,1),9)

# Matrices for GfosA
cortest_Ho_k1 <- cortest_Ho[cortest_Ho$k==1,]
cortest_Ho_k1_matrix_A <- list()
cortest_Ho_k1_matrix_A[[1]] <- cortest_Ho_k1$tau_A[cortest_Ho_k1$d==0.001]
cortest_Ho_k1_matrix_A[[2]] <- cortest_Ho_k1$tau_A[cortest_Ho_k1$d==0.01]
cortest_Ho_k1_matrix_A[[3]] <- cortest_Ho_k1$tau_A[cortest_Ho_k1$d==0.1]
cortest_Ho_k1_matrix_A = do.call(cbind, cortest_Ho_k1_matrix_A)
colnames(cortest_Ho_k1_matrix_A) <- c("d 0.001","d 0.01","d 0.1")
row.names(cortest_Ho_k1_matrix_A) <- c("w 0","w 0.5","w 1")
cortest_Ho_k0 <- cortest_Ho[cortest_Ho$k==0,]
cortest_Ho_k0_matrix_A <- list()
cortest_Ho_k0_matrix_A[[1]] <- cortest_Ho_k0$tau_A[cortest_Ho_k0$d==0.001]
cortest_Ho_k0_matrix_A[[2]] <- cortest_Ho_k0$tau_A[cortest_Ho_k0$d==0.01]
cortest_Ho_k0_matrix_A[[3]] <- cortest_Ho_k0$tau_A[cortest_Ho_k0$d==0.1]
cortest_Ho_k0_matrix_A = do.call(cbind, cortest_Ho_k0_matrix_A)
colnames(cortest_Ho_k0_matrix_A) <- c("d 0.001","d 0.01","d 0.1")
row.names(cortest_Ho_k0_matrix_A) <- c("w 0","w 0.5","w 1")

# Matrices for GfosB
cortest_Ho_k1_matrix_B <- list()
cortest_Ho_k1_matrix_B[[1]] <- cortest_Ho_k1$tau_B[cortest_Ho_k1$d==0.001]
cortest_Ho_k1_matrix_B[[2]] <- cortest_Ho_k1$tau_B[cortest_Ho_k1$d==0.01]
cortest_Ho_k1_matrix_B[[3]] <- cortest_Ho_k1$tau_B[cortest_Ho_k1$d==0.1]
cortest_Ho_k1_matrix_B = do.call(cbind, cortest_Ho_k1_matrix_B)
colnames(cortest_Ho_k1_matrix_B) <- c("d 0.001","d 0.01","d 0.1")
row.names(cortest_Ho_k1_matrix_B) <- c("w 0","w 0.5","w 1")
cortest_Ho_k0_matrix_B <- list()
cortest_Ho_k0_matrix_B[[1]] <- cortest_Ho_k0$tau_B[cortest_Ho_k0$d==0.001]
cortest_Ho_k0_matrix_B[[2]] <- cortest_Ho_k0$tau_B[cortest_Ho_k0$d==0.01]
cortest_Ho_k0_matrix_B[[3]] <- cortest_Ho_k0$tau_B[cortest_Ho_k0$d==0.1]
cortest_Ho_k0_matrix_B = do.call(cbind, cortest_Ho_k0_matrix_B)
colnames(cortest_Ho_k0_matrix_B) <- c("d 0.001","d 0.01","d 0.1")
row.names(cortest_Ho_k0_matrix_B) <- c("w 0","w 0.5","w 1")

#### VALUES -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
nrow(A_red)+nrow(B) # Number of used individuals
if (internal){
  nrow(Gfos[which(Gfos$type!="K"),])
}

min(meanAr_A_red)
max(meanAr_A_red)
mean(meanAr_A_red)
median(meanAr_A_red)
sd(meanAr_A_red)
min(meanAr_B)
max(meanAr_B)
mean(meanAr_B)
median(meanAr_B)
sd(meanAr_B)

min(meanHo_A_red)
max(meanHo_A_red)
mean(meanHo_A_red)
median(meanHo_A_red)
sd(meanHo_A_red)
min(meanHo_B)
max(meanHo_B)
mean(meanHo_B)
median(meanHo_B)
sd(meanHo_B)

# Get range of scaled carrying capacities (from C++ file, lines 166, 215)
# Extract total catchments sizes from graph object (since included here, otherwise use catch_area.in directly)
total_catch <- V(net)$Total_Catch
# K scaling according to sqrt(catchment size)
sqrt(total_catch)
sum_sqrt_catch_size <- sum(sqrt(total_catch))
K_scaled <- round(sqrt(total_catch) * (1000*length(total_catch)/sum_sqrt_catch_size ))
range(K_scaled)

#### FIGURES -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
dir.create(paste0(WD,"/Analysis_",output), showWarnings=F)

##** Preparing multipanel figures ####
pla <- list(
  a=2, # columns of subplots
  b=1, # rows of subplots
  x=3, # number of rows
  y=3, # number of cols
  sub1=F,
  sub2=T,
  main=F,
  main_t=NULL,
  sub1_t=NULL,
  sub2_t=klab,
  x_t=wlab,
  y_t=dlab,
  main_c=NULL,
  sub1_c=3,
  sub2_c=2.5,
  x_c=2,
  y_c=2)
lab <- c(pla$sub1_t,rep(c(pla$sub2_t),ifelse(pla$sub2,pla$b,0)),rep(pla$x_t,pla$b),rep(pla$y_t,pla$a))
lab_cex <- c(rep(pla$sub1_c,length(pla$sub1_t)),rep(c(pla$sub2_c),ifelse(pla$sub2,pla$a*pla$b,0)),rep(pla$x_c,pla$b*length(pla$x_t)),rep(pla$y_c,pla$a*length(pla$y_t)))

##** FIG 1 ####
## Empirical data maps
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/Fig1.pdf"), width=fig.width, height=(2*fig.width)/3)
}else{
  png(paste0(WD,"/Analysis_",output,"/Fig1.png"), width=fig.width, height=(2*fig.width)/3, units="in", res=300)
}
nf <- layout(matrix(c(9, 5, 6,
                      7, 1, 2,
                      8, 3, 4), nrow=3, byrow=T),
             widths=c(0.5,3,3),
             heights=c(0.5,2,2),
             respect=T)
# layout.show(nf)
figmar <- c(4.5,4.5,0.5,0)
figmgp <- c(1.7,0.3,0)
col_switch <- 1/2
# #### Map of mean allelic richness in G. fossarum A
par(mar=figmar, mgp=figmgp, tcl=0.2, xaxs="i", yaxs="i")
if (internal){
  river_plot(north_arrow = F, overview_map = F, scalebar=T, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, lines_rhone=F, plot_ticino=F, lines_ticino=F, plot_inn=F, lines_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}else{
  river_plot(north_arrow = F, overview_map = F, scalebar=T, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, plot_ticino=F, plot_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}
values <- (meanAr_A_red-min(meanAr_A_red))/ (max(meanAr_A_red)-min(meanAr_A_red)) # transform meanAr values to range [0,1] for heatmap plotting
for(i in 1:length(match_A_red)){
  points(site_coord$x[match_A_red[i]], site_coord$y[match_A_red[i]], col=rgb2hex(col_fun(values))[i], pch=19, cex=2.5)
  l1 <- mean(meanAr_A_red)+(col_switch)*(max(meanAr_A_red)-median(meanAr_A_red))
  l2 <- mean(meanAr_A_red)-(col_switch)*(median(meanAr_A_red)-min(meanAr_A_red))
    text(site_coord$x[match_A_red[i]], site_coord$y[match_A_red[i]], round(meanAr_A_red[i],1), col=ifelse(meanAr_A_red[i]>l1|meanAr_A_red[i]<l2,"white","black"), cex=0.8)
}
overview_map(xl = 495000,xr = 545000,yt = 280000, yb = 230000)
north_arrow(x=825000, y=80000)
gradient.legend(meanAr_A_red, val.cex = 1, palette=col_pal)
mtext(side = 3, text = lab_sub[1], line = 0.5, adj=0, cex = 1.5)

#### Map of mean observed heterozygosity in G. fossarum A
par(mar=figmar, mgp=figmgp, tcl=0.2, xaxs="i", yaxs="i")
if (internal){
  river_plot(north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, lines_rhone=F, plot_ticino=F, lines_ticino=F, plot_inn=F, lines_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}else{
  river_plot(north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, plot_ticino=F, plot_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}
mtext(side = 2, text = expression(paste(bold("CH1903"), " / ", bold("LV03"))), line = 1 + 1, cex = 1)
values <- (meanHo_A_red-min(meanHo_A_red))/ (max(meanHo_A_red)-min(meanHo_A_red)) # transform meanHo values to range [0,1] for heatmap plotting
for(i in 1:length(match_A_red)){
  points(site_coord$x[match_A_red[i]], site_coord$y[match_A_red[i]], col=rgb2hex(col_fun(values))[i], pch=19, cex=2.5)
  l1 <- mean(meanHo_A_red)+(col_switch)*(max(meanHo_A_red)-median(meanHo_A_red))
  l2 <- mean(meanHo_A_red)-(col_switch)*(median(meanHo_A_red)-min(meanHo_A_red))
  text(site_coord$x[match_A_red[i]], site_coord$y[match_A_red[i]], round(meanHo_A_red[i],1), col=ifelse(meanHo_A_red[i]>l1|meanHo_A_red[i]<l2,"white","black"), cex=0.8)
}
gradient.legend(meanHo_A_red, val.cex=1, palette=col_pal)
mtext(side = 3, text = lab_sub[3], line = 0.5, adj=0, cex = 1.5)

#### Map of mean allelic richness in G. fossarum B
par(mar=figmar, mgp=figmgp, tcl=0.2, xaxs="i", yaxs="i")
if (internal){
  river_plot(north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, lines_rhone=F, plot_ticino=F, lines_ticino=F, plot_inn=F, lines_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}else{
  river_plot(north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, plot_ticino=F, plot_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}
mtext(side = 1, text = expression(paste(bold("CH1903"), " / ", bold("LV03"))), line = 1 + 1, cex = 1)
values <- (meanAr_B-min(meanAr_B))/ (max(meanAr_B)-min(meanAr_B)) # transform meanAr values to range [0,1] for heatmap plotting
for(i in 1:length(match_B)){
  points(site_coord$x[match_B[i]], site_coord$y[match_B[i]], col=rgb2hex(col_fun(values))[i], pch=19, cex=2.5)
  l1 <- mean(meanAr_B)+(col_switch)*(max(meanAr_B)-median(meanAr_B))
  l2 <- mean(meanAr_B)-(col_switch)*(median(meanAr_B)-min(meanAr_B))
  text(site_coord$x[match_B[i]], site_coord$y[match_B[i]], round(meanAr_B[i],1), col=ifelse(meanAr_B[i]>l1|meanAr_B[i]<l2,"white","black"), cex=0.8)
}
gradient.legend(meanAr_B, val.cex = 1, palette=col_pal)
mtext(side = 3, text = lab_sub[2], line = 0.5, adj=0, cex = 1.5)

#### Map of mean observed heterozygosity in G. fossarum B
par(mar=figmar, mgp=figmgp, tcl=0.2, xaxs="i", yaxs="i")
if (internal){
  river_plot(north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, lines_rhone=F, plot_ticino=F, lines_ticino=F, plot_inn=F, lines_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}else{
  river_plot(north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, plot_ticino=F, plot_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}
mtext(side = 1, text = expression(paste(bold("CH1903"), " / ", bold("LV03"))), line = 1 + 1, cex = 1)
mtext(side = 2, text = expression(paste(bold("CH1903"), " / ", bold("LV03"))), line = 1 + 1, cex = 1)
values <- (meanHo_B-min(meanHo_B))/ (max(meanHo_B)-min(meanHo_B)) # transform meanHo values to range [0,1] for heatmap plotting
for(i in 1:length(match_B)){
  points(site_coord$x[match_B[i]], site_coord$y[match_B[i]], col=rgb2hex(col_fun(values))[i], pch=19, cex=2.5)
  l1 <- mean(meanHo_B)+(col_switch)*(max(meanHo_B)-median(meanHo_B))
  l2 <- mean(meanHo_B)-(col_switch)*(median(meanHo_B)-min(meanHo_B))
  text(site_coord$x[match_B[i]], site_coord$y[match_B[i]], round(meanHo_B[i],1), col=ifelse(meanHo_B[i]>l1|meanHo_B[i]<l2,"white","black"), cex=0.8)
}
gradient.legend(meanHo_B, val.cex = 1, palette=col_pal)
mtext(side = 3, text = lab_sub[4], line = 0.5, adj=0, cex = 1.5)

par(mar=c(0,0,0,0),mgp=c(3,1,0))
plot.new()
text(0.5,0.5, lab_Ar ,adj=c(0.5,0.5), cex=2.5)
plot.new()
text(0.5,0.5, lab_Ho ,adj=c(0.5,0.5), cex=2.5)
plot.new()
text(1,0.5,label_A,adj=c(0.5,0), cex=2.5, srt=90)
plot.new()
text(1,0.5,label_B,adj=c(0.5,0), cex=2.5, srt=90)
dev.off()

##** FIG 2 ####
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/Fig2.pdf"), width=15, height=10)
}else{
  png(paste0(WD,"/Analysis_",output,"/Fig2.png"), width=15, height=10, units="in", res=300)
}
par(mfrow=c(2,2))
#### meanAr~updist
GLMplot("meanAr","updist",DATA,model="sglm_Ar_updist",pt.cex=2, CI_border=F,xlabel="Upstream distance [km]",ylabel="Mean allelic richness", xrev=T, xax="n", cex.lab=2, cex.axis=1.5, cex.legend=1.5)
axis(1, c(0,50,100,150,200,250,300), at=c(0,50000,100000,150000,200000,250000,300000), cex.axis=1.5)
text(300000, par("usr")[4], lab_sub[1], cex=2, adj=c(0.5,1))
#### meanHo~updist
GLMplot("meanHo","updist",DATA,model="glm_Ho_updist",pt.cex=2, CI_border=F,xlabel="Upstream distance [km]",ylabel="Mean observed heterozygosity", xrev=T, xax="n", cex.lab=2, cex.axis=1.5, legend=F)
axis(1, c(0,50,100,150,200,250,300), at=c(0,50000,100000,150000,200000,250000,300000), cex.axis=1.5)
text(300000, par("usr")[4], lab_sub[3], cex=2, adj=c(0.5,1))
#### meanAr~std_clos_undir
GLMplot("meanAr","std_clos_undir",DATA,model="sglm_Ar_clos",pt.cex=2, CI_border=F,xlabel="Standardized closeness centrality",ylabel="Mean allelic richness", legend=F)
text(0, par("usr")[4], lab_sub[2], cex=2, adj=c(0.5,1))
#### meanHo~std_clos_undir
GLMplot("meanHo","std_clos_undir",DATA,model="glm_Ho_clos",pt.cex=2, CI_border=F,xlabel="Standardized closeness centrality",ylabel="Mean observed heterozygosity", legend=F)
text(0, par("usr")[4], lab_sub[4], cex=2, adj=c(0.5,1))
dev.off()

##** FIG 3 ####
##** Modelled data maps
#*** Mean Ar maps
main_3=F
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,sub.row=pla$x,sub.col=pla$y,sub1=pla$sub1,sub2=pla$sub2, main=main_3, h.main=0.5, w.legend=0, h.sub2=0.3, w.axis=0.7, h.axis=0, spacer.sub.col=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/Fig3.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/Fig3.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
par(mar=c(0,0,0,0))
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,sub.row=pla$x,sub.col=pla$y,sub1=pla$sub1,sub2=pla$sub2, main=main_3, h.main=0.5, w.legend=0, h.sub2=0.3, w.axis=0.7, h.axis=0, spacer.sub.col=0.5)

for (j in 2:ncol(Ar_modelled)){
  if (internal){
    river_plot(width_country=0.5, lwd_rivers=0.5, lwd_lakes=0.5, xlimit=c(495000,825000), col_water=col_water, north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, plot_rhone=F, lines_rhone=F, plot_ticino=F, lines_ticino=F, plot_inn=F, lines_inn=F, lakes=TRUE, rivers=TRUE, axes="none", river_nr=T)
  }else{
    river_plot(width_country=0.5, lwd_rivers=0.5, lwd_lakes=0.5, xlimit=c(495000,825000), col_water=col_water, north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, plot_rhone=F, plot_ticino=F, plot_inn=F, lakes=TRUE, rivers=TRUE, axes="none", river_nr=F)
  }
  values <- (Ar_modelled[,j]-min(Ar_modelled[,j]))/ (max(Ar_modelled[,j])-min(Ar_modelled[,j])) # transform meanAr values to range [0,1] for heatmap plotting
  for(i in 1:length(match_Mod)){
    points(site_coord$x[match_Mod[i]], site_coord$y[match_Mod[i]], bg=rgb2hex(col_fun(values))[i], pch=21, lwd=0.5, cex=1)
  }
  if (j %in% c(2,3)){
    mtext(lab[6], side=3, line=-1, cex=1.3)
  }
  if (j %in% c(8,9)){
    mtext(lab[7], side=3, line=-1, cex=1.3)
  }
  if (j %in% c(14,15)){
    mtext(lab[8], side=3, line=-1, cex=1.3)
  }
  gradient.legend(Ar_modelled[,j],alpha=1, val.midpoint=F, round=1, val.cex=1, val.gap = 0.5, title.gap=0.1, xl = 505000, xr = 635000, yb = 20000, yt = 40000, horizontal=T)
}
if(main_3){
  plot.new()
  text(0.5,0.5,"Simulated mean allelic richness",adj=c(0.5,0.5), cex=3)
}
for (i in 1:5){
  plot.new()
  if (i %in% c(1:2)){
    text(0.5,1,lab[i],adj=c(0.5,1), cex=lab_cex[i])
  }
  if (i %in% c(3:5)){
    text(0.5,0.5,lab[i],adj=c(0.5,0.5), cex=lab_cex[i])
  }
}
dev.off()

##** FIG 4 ####
##** Histogram of orthogonal distance to 1:1 line
main_4=F
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_4, h.main=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/Fig4.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/Fig4.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_4, h.main=0.5)

for (i in 1:ncol(hist_Ar_A)){
  par(mar=c(0,0,0,0))
  yclip <- 30
  ylim <- 35
  hist(hist_Ar_A[,i],
       breaks=seq(0,ceiling(max(hist_Ar_A,hist_Ar_B))-0.5,0.5),
       xlim=c(0,ceiling(max(hist_Ar_A,hist_Ar_B))),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="")
  clip(0,ceiling(max(hist_Ar_A,hist_Ar_B)), 0, yclip)
  abline(v=median(hist_Ar_A[,i]),col=col_Gfos_A)
  abline(v=median(hist_Ar_B[,i]),col=col_Gfos_B)
  clip(0,ceiling(max(hist_Ar_A,hist_Ar_B)), 0, ylim)
  textbox(ceiling(max(hist_Ar_A,hist_Ar_B)),30,paste0(measure2_short," = ",formatC(round(median(hist_Ar_A[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_A, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  textbox(ceiling(max(hist_Ar_A,hist_Ar_B)),25,paste0(measure2_short," = ",formatC(round(median(hist_Ar_B[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_B, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  hist(hist_Ar_A[,i],
       breaks=seq(0,ceiling(max(hist_Ar_A,hist_Ar_B))-0.5,0.5),
       xlim=c(0,ceiling(max(hist_Ar_A,hist_Ar_B))),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="", add=T)
  hist(hist_Ar_B[,i],
       breaks=seq(0,ceiling(max(hist_Ar_A,hist_Ar_B))-0.5,0.5),
       xlim=c(0,ceiling(max(hist_Ar_A,hist_Ar_B))),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_B,
       main="",
       add=T)
 
  if (i %in% c(1,2)){
    mtext(lab[6], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(7,8)){
    mtext(lab[7], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(13,14)){
    mtext(lab[8], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(3)){
    mtext("Frequencies [counts]", side=2, line=1, cex=1)
  }
  if (i %in% c(11,12)){
    mtext(measure3, side=1, line=3, cex=1)
  }
  if (i %in% c(5:6,11:12,17:18)){
    axis(1)
  }
  if (i %in% c(1,3,5)){
    axis(2, at=seq(0,yclip,5), labels=F)
  }
  if (i %in% c(2,4,6)){
    axis(2, at=seq(0,yclip,5))
  }
  if (i %in% c(13:18)){
    axis(4, at=seq(0,yclip,5), labels=F)
  }
}

if(main_4){
  plot.new()
  text(0.5,0.7,paste0("Mean allelic richness: Distribution of ",measure3a),adj=c(0.5,0.5),cex=3)
}

for (i in 1:length(lab)){
  plot.new()
  if (i %in% c(1:2)){
    text(0.5,1,lab[i],adj=c(0.5,1), cex=lab_cex[i])
  }
  if (i %in% c(3:5)){
    text(0,0.5,lab[i],adj=c(0,0.5), cex=lab_cex[i])
  }
}
# Add example plot
par(mar=c(0,6,16,0.5), pty="s")
i <- ncol(Ar_modelled)
Ar_Mod_A <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosA]
Ar_Mod_B <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosB]
plot(Ar_Mod_A~meanAr_A_red_updist, type="n",
     xlim=c(min(Ar_modelled[,-1]),max(Ar_modelled[,-1])),
     ylim=c(min(Ar_modelled[,-1]),max(Ar_modelled[,-1])),
     xaxt="s",
     yaxt="n",
     xlab="",
     ylab="",
     asp=1)

par(xpd=T)
polygon(c(-4.5,-4.5,par("usr")[3],par("usr")[3]),c(par("usr")[1],9,par("usr")[4],par("usr")[1]), col="lightgrey", border=NA)
segments(-4.5,par("usr")[1],par("usr")[3],par("usr")[1], lty=1, col="grey", lwd=1.5)
segments(-4.5,9,par("usr")[3],par("usr")[4], lty=1, col="grey", lwd=1.5)

text(-2.5,18,"Example plot for\nperpendicular offsets", cex=2, adj=0, col="darkgrey")
text(-2.5,15,"See Fig. S6 for all plots",cex=1, adj=0, col="darkgrey")

par(xpd=F)

mtext(expression(bold("Model data:")*" Ar"), side=2, line=2, cex=0.8)
mtext(expression(bold("Empirical data:")*" Ar"), side=1, line=2.5, cex=0.8)
box()
axis(3,labels=F)
axis(2)

abline(0,1, lwd=1, lty=2) # add 1:1 line

points(Ar_Mod_A~meanAr_A_red_updist,col=col_Gfos_A, pch=16)
points(Ar_Mod_B~meanAr_B_updist,col=col_Gfos_B, pch=16)

for (j in 1:length(meanAr_A_red_updist)){
  point <- cbind(meanAr_A_red_updist,Ar_Mod_A)[j,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  segments(seg[1],seg[2],seg[3],seg[4], col=col_Gfos_A, lty=2, lwd=0.5)
}
for (k in 1:length(meanAr_B_updist)){
  point <- cbind(meanAr_B_updist,Ar_Mod_B)[k,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  segments(seg[1],seg[2],seg[3],seg[4], col=col_Gfos_B, lty=2, lwd=0.5)
}
text(max(Ar_modelled[,-1]),min(Ar_modelled[,-1])+1.5,paste0(measure1_short, " = ",formatC(round(sum(orthdist_Ar_A),sum_digits),digits=sum_digits, format="f")), adj=1, col=col_Gfos_A)
text(max(Ar_modelled[,-1]),min(Ar_modelled[,-1])+0.5,paste0(measure1_short," = ",formatC(round(sum(orthdist_Ar_B),sum_digits),digits=sum_digits, format="f")), adj=1, col=col_Gfos_B)
par(xpd=T)
legend(-5,30,c(label_A,label_B),pch = 16, col = c(col_Gfos_A,col_Gfos_B), bty="n", cex=2)

# text(0,0,"SOSO = Sum of squared orthogonals", adj=0)
par(xpd=F)
dev.off()

##** FIG 5 ####
##** Histogram of orthogonal distance to 1:1 line
main_5=F
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_5, h.main=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/Fig5.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/Fig5.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_5, h.main=0.5)

for (i in 1:ncol(hist_Ho_A)){
  par(mar=c(0,0,0,0))
  yclip <- 30
  ylim <- 35
  hist(hist_Ho_A[,i],
       breaks=seq(0,round(max(hist_Ho_A,hist_Ho_B),1),0.05),
       xlim=c(0,round(max(hist_Ho_A,hist_Ho_B),1)),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="")
  clip(0,round(max(hist_Ho_A,hist_Ho_B),1), 0, yclip)
  abline(v=median(hist_Ho_A[,i]),col=col_Gfos_A)
  abline(v=median(hist_Ho_B[,i]),col=col_Gfos_B)
  clip(0,round(max(hist_Ho_A,hist_Ho_B),1), 0, ylim)
  textbox(round(max(hist_Ho_A,hist_Ho_B),1),30,paste0(measure2_short," = ",formatC(round(median(hist_Ho_A[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_A, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  textbox(round(max(hist_Ho_A,hist_Ho_B),1),25,paste0(measure2_short," = ",formatC(round(median(hist_Ho_B[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_B, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  hist(hist_Ho_A[,i],
       breaks=seq(0,round(max(hist_Ho_A,hist_Ho_B),1),0.05),
       xlim=c(0,round(max(hist_Ho_A,hist_Ho_B),1)),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="", add=T)
  hist(hist_Ho_B[,i],
       breaks=seq(0,round(max(hist_Ho_A,hist_Ho_B),1),0.05),
       xlim=c(0,round(max(hist_Ho_A,hist_Ho_B),1)),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_B,
       main="",
       add=T)
  
  if (i %in% c(1,2)){
    mtext(lab[6], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(7,8)){
    mtext(lab[7], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(13,14)){
    mtext(lab[8], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(3)){
    mtext("Frequencies [counts]", side=2, line=1, cex=1)
  }
  if (i %in% c(11,12)){
    mtext(measure3, side=1, line=3, cex=1)
  }
  if (i %in% c(5:6,11:12,17:18)){
    axis(1, at=seq(0,round(max(hist_Ho_A,hist_Ho_B),1),0.1), labels=c("0.0","","0.2","","0.4","",""))
  }
  if (i %in% c(1,3,5)){
    axis(2, at=seq(0,yclip,5), labels=F)
  }
  if (i %in% c(2,4,6)){
    axis(2, at=seq(0,yclip,5))
  }
  if (i %in% c(13:18)){
    axis(4, at=seq(0,yclip,5), labels=F)
  }
}

if(main_5){
  plot.new()
  text(0.5,0.7,paste0("Mean observed heterozygosity: Distribution of ",measure3a),adj=c(0.5,0.5),cex=3)
}

for (i in 1:length(lab)){
  plot.new()
  if (i %in% c(1:2)){
    text(0.5,1,lab[i],adj=c(0.5,1), cex=lab_cex[i])
  }
  if (i %in% c(3:5)){
    text(0,0.5,lab[i],adj=c(0,0.5), cex=lab_cex[i])
  }
}
par(mar=c(0,6,16,0.5), pty="s")
i <- ncol(Ar_modelled)
Ar_Mod_A <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosA]
Ar_Mod_B <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosB]
plot(Ar_Mod_A~meanAr_A_red_updist, type="n",
     xlim=c(min(Ar_modelled[,-1]),max(Ar_modelled[,-1])),
     ylim=c(min(Ar_modelled[,-1]),max(Ar_modelled[,-1])),
     xaxt="n",
     yaxt="n",
     xlab="",
     ylab="",
     bty="n",
     asp=1)

par(xpd=T)
legend(-5,yclip,c(label_A,label_B),pch = 16, col = c(col_Gfos_A,col_Gfos_B), bty="n", cex=2)
par(xpd=F)
dev.off()

#### SUPP INFO -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
dir.create(paste0(WD,"/Analysis_",output,"/SuppFigs"), showWarnings=F)

##** FIG S1 ####
# create an overview map (once at the end of looping over parameter space)
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS1.pdf"), width=8, height=6)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS1.png"), width=8, height=6, units="in", res=300)
}
par(mar=c(0,0,0,0))
V(net)$modsite <- ifelse(V(net)$name%in%modsite,2,1)
if (internal){
  river_plot(overview_map = F, col_rhine=col_rhine, col_rhone=NA, col_ticino = NA, col_inn=NA, col_water = col_water, axes="none")
}else{
  river_plot(overview_map = F, col_rhine=col_rhine, col_rhone=NA, col_ticino = NA, col_inn=NA, col_water = col_water, axes="none")
}
plot(net, layout=net_layout, edge.arrow.size=0, edge.width=2.5, edge.color="blue",
     vertex.color=c(adjustcolor("white",0),"dark red")[V(net)$modsite], vertex.size=ifelse(V(net)$modsite==2,1,1), vertex.frame.color=col_water, vertex.shape="none",
     vertex.label=NA,
     rescale=F, xlim=c(min(net_layout[,1]), max(net_layout[,1])), ylim = c(min(net_layout[,2]), max(net_layout[,2])), asp = 0, add=T) # look at the existing/preloaded data
cex_small=0.5
cex_big=0.75
if (internal){
  points(Rhine$chx[Rhine$sample==0], Rhine$chy[Rhine$sample==0], pch=21, col="black", bg="white", cex=cex_small, lwd=0.75) # plot empty samples
  points(Rhine$chx[Rhine$sample==1], Rhine$chy[Rhine$sample==1], pch=21, col="black", bg="black", cex=cex_small, lwd=0.75) # plot amphipod data
  points(Gfos$chx[Gfos$type=="K"], Gfos$chy[Gfos$type=="K"], pch=21, col="black", bg=col_Gfos, cex=cex_big) # plot Gammarus fossarum complex
  points(Gfos$chx[Gfos$type=="A"], Gfos$chy[Gfos$type=="A"], pch=21, col="black", bg=col_Gfos_A, cex=cex_big) # plot Gammarus fossarum type A
  points(Gfos$chx[Gfos$type=="B"], Gfos$chy[Gfos$type=="B"], pch=21, col="black", bg=col_Gfos_B, cex=cex_big) # plot Gammarus fossarum type B
  points(Gfos$chx[Gfos$type=="C"], Gfos$chy[Gfos$type=="C"], pch=21, col="black", bg=col_Gfos, cex=cex_big) # plot Gammarus fossarum type C
}
points(site_coord$x, site_coord$y, col=c(1,1,0)[site_coord$gendata], bg=c(col_Gfos_A,col_Gfos_B,0)[site_coord$gendata], cex=1.5, pch=24)
points(site_coord$x[site_coord$gendata==3], site_coord$y[site_coord$gendata==3], col=1, bg=col_Gfos_B, cex=1.5, pch=24)
points(site_coord$x[site_coord$gendata==3], site_coord$y[site_coord$gendata==3], col=col_Gfos_A, cex=0.75, pch=17)
text(465000,315000,"Microsat data (> 15 ind.)", adj=0)
legend(465000,315000, c(label_A,label_B,"Both"), col=c(1,1,1), pt.bg=c(col_Gfos_A,col_Gfos_B,col_Gfos_B), pt.cex=1.5, pch=24, bty="n")
legend(465000,315000, c("","",""), col=c(col_Gfos_A,col_Gfos_B,col_Gfos_A), pt.cex=0.75, pch=17, bty="n")
text(465000,265000,"Presence/Absence data", adj=0, cex=0.8)
legend(465000,265000, cex=0.8, c(expression(italic(G. ~ fossarum) ~ "complex"), label_A,label_B, "Amphipods present", "Amphipods absent"), col=c(1,1,1,1,1), pt.bg=c(col_Gfos,col_Gfos_A,col_Gfos_B,"black","white"), pt.cex=c(cex_big,cex_big,cex_big,cex_small,cex_small), pch=21, bty="n")
dev.off()

##** FIG S2 ####
#*** Mean Ho maps
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,sub.row=pla$x,sub.col=pla$y,sub1=pla$sub1,sub2=pla$sub2, main=T, h.main=0.5, w.legend=0, h.sub2=0.3, w.axis=0.7, h.axis=0, spacer.sub.col=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS2.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS2.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
par(mar=c(0,0,0,0))
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,sub.row=pla$x,sub.col=pla$y,sub1=pla$sub1,sub2=pla$sub2, main=T, h.main=0.5, w.legend=0, h.sub2=0.3, w.axis=0.7, h.axis=0, spacer.sub.col=0.5)

for (j in 2:ncol(Ho_modelled)){
  if (internal){
    river_plot(width_country=0.5, lwd_rivers=0.5, xlimit=c(495000,825000), col_water=col_water, north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, plot_rhone=F, lines_rhone=F, plot_ticino=F, lines_ticino=F, plot_inn=F, lines_inn=F, lakes=TRUE, rivers=TRUE, axes="none", river_nr=T)
  }else{
    river_plot(width_country=0.5, lwd_rivers=0.5, xlimit=c(495000,825000), col_water=col_water, north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, plot_rhone=F, plot_ticino=F, plot_inn=F, lakes=TRUE, rivers=TRUE, axes="none", river_nr=F)
  }
  values <- (Ho_modelled[,j]-min(Ho_modelled[,j]))/ (max(Ho_modelled[,j])-min(Ho_modelled[,j])) # transform meanAr values to range [0,1] for heatmap plotting
  for(i in 1:length(match_Mod)){
    points(site_coord$x[match_Mod[i]], site_coord$y[match_Mod[i]], bg=rgb2hex(col_fun(values))[i], pch=21, lwd=0.5, cex=1)
  }
  if (j %in% c(2,3)){
    mtext(lab[6], side=3, line=-1, cex=1.3)
  }
  if (j %in% c(8,9)){
    mtext(lab[7], side=3, line=-1, cex=1.3)
  }
  if (j %in% c(14,15)){
    mtext(lab[8], side=3, line=-1, cex=1.3)
  }
  gradient.legend(Ho_modelled[,j],alpha=1, val.midpoint=F, round=2, val.cex=1, val.gap = 0.5, title.gap=0.1, xl = 505000, xr = 635000, yb = 20000, yt = 40000, horizontal=T)
}
plot.new()
text(0.5,0.5,"Simulated heterozygosity",adj=c(0.5,0.5), cex=3)
for (i in 1:5){
  plot.new()
  if (i %in% c(1:2)){
    text(0.5,1,lab[i],adj=c(0.5,1), cex=lab_cex[i])
  }
  if (i %in% c(3:5)){
    text(0.5,0.5,lab[i],adj=c(0.5,0.5), cex=lab_cex[i])
  }
}
dev.off()

##** FIG S3 ####
##** Perpendicular offset histogram comparison for d
maxhist_Ar <- ceiling(max(hist_Ar_A,hist_Ar_B))
maxhist_Ho <- ceiling(10*max(hist_Ho_A,hist_Ho_B))/10
d0001_Ar_A <- hist_Ar_A[,(1:6)]
d001_Ar_A <- hist_Ar_A[,c(7:12)]
d01_Ar_A <- hist_Ar_A[,c(13:18)]
d0001_Ar_B <- hist_Ar_B[,(1:6)]
d001_Ar_B <- hist_Ar_B[,c(7:12)]
d01_Ar_B <- hist_Ar_B[,c(13:18)]
d0001_Ho_A <- hist_Ho_A[,(1:6)]
d001_Ho_A <- hist_Ho_A[,c(7:12)]
d01_Ho_A <- hist_Ho_A[,c(13:18)]
d0001_Ho_B <- hist_Ho_B[,(1:6)]
d001_Ho_B <- hist_Ho_B[,c(7:12)]
d01_Ho_B <- hist_Ho_B[,c(13:18)]

diffSPO_d0001both_Ar_A <- c(apply(d0001_Ar_A,2,sum)-apply(d001_Ar_A,2,sum),apply(d0001_Ar_A,2,sum)-apply(d01_Ar_A,2,sum))
diffSPO_d0001both_Ar_B <- c(apply(d0001_Ar_B,2,sum)-apply(d001_Ar_B,2,sum),apply(d0001_Ar_B,2,sum)-apply(d01_Ar_B,2,sum))
diffSPO_d0001both_Ho_A <- c(apply(d0001_Ho_A,2,sum)-apply(d001_Ho_A,2,sum),apply(d0001_Ho_A,2,sum)-apply(d01_Ho_A,2,sum))
diffSPO_d0001both_Ho_B <- c(apply(d0001_Ho_B,2,sum)-apply(d001_Ho_B,2,sum),apply(d0001_Ho_B,2,sum)-apply(d01_Ho_B,2,sum))
diffMPO_d0001both_Ar_A <- c(apply(d0001_Ar_A,2,median)-apply(d001_Ar_A,2,median),apply(d0001_Ar_A,2,median)-apply(d01_Ar_A,2,median))
diffMPO_d0001both_Ar_B <- c(apply(d0001_Ar_B,2,median)-apply(d001_Ar_B,2,median),apply(d0001_Ar_B,2,median)-apply(d01_Ar_B,2,median))
diffMPO_d0001both_Ho_A <- c(apply(d0001_Ho_A,2,median)-apply(d001_Ho_A,2,median),apply(d0001_Ho_A,2,median)-apply(d01_Ho_A,2,median))
diffMPO_d0001both_Ho_B <- c(apply(d0001_Ho_B,2,median)-apply(d001_Ho_B,2,median),apply(d0001_Ho_B,2,median)-apply(d01_Ho_B,2,median))

total_comparison_d0001both <- length(c(diffMPO_d0001both_Ar_A,diffMPO_d0001both_Ar_B,diffMPO_d0001both_Ho_A,diffMPO_d0001both_Ho_B,diffSPO_d0001both_Ar_A,diffSPO_d0001both_Ar_B,diffSPO_d0001both_Ho_A,diffSPO_d0001both_Ho_B))
improved_fit_d0001both <- sum(c(diffMPO_d0001both_Ar_A,diffMPO_d0001both_Ar_B,diffMPO_d0001both_Ho_A,diffMPO_d0001both_Ho_B,diffSPO_d0001both_Ar_A,diffSPO_d0001both_Ar_B,diffSPO_d0001both_Ho_A,diffSPO_d0001both_Ho_B)<0)
improved_fit_d0001both/total_comparison_d0001both


if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS3.pdf"), width=8, height=6)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS3.png"), width=8, height=6, units="in", res=300)
}
op <- par(mfrow = c(3,2),
          oma = c(5,5,2,0) + 0.1,
          mar = c(0,0,2,1) + 0.1)
hist(d0001_Ar_A, col=col_Gfos_A, cex.axis=2, xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ar), breaks=seq(0,maxhist_Ar,0.25))
hist(d0001_Ar_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ar,0.25))
axis(1, labels=F, cex.axis=2)
mtext(side = 3, text = "Mean allelic richness", line = 1, adj=0.5, cex = 1.5)
abline(v=median(d0001_Ar_A), col=white_transparent, lwd=2)
abline(v=median(d0001_Ar_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(d0001_Ar_B), col=white_transparent, lwd=2)
abline(v=median(d0001_Ar_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(4,65, txt = lab[6], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)
legend(2,50, c(expression("Median" ~ italic(G. ~ fossarum) ~ "type A"),expression("Median" ~ italic(G. ~ fossarum) ~ "type B")),
       lty=2, col=c(col_Gfos_A,col_Gfos_B), bty="n", cex=1.3, lwd=2)

hist(d0001_Ho_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ho), breaks=seq(0,maxhist_Ho,0.025))
hist(d0001_Ho_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ho,0.025))
axis(1, labels=F, cex.axis=2)
axis(2, labels=F, cex.axis=2)
mtext(side = 3, text = "Mean observed heterozygosity", line = 1, adj=0.5, cex = 1.5)
abline(v=median(d0001_Ho_A), col=white_transparent, lwd=2)
abline(v=median(d0001_Ho_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(d0001_Ho_B), col=white_transparent, lwd=2)
abline(v=median(d0001_Ho_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,65, txt = lab[6], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

hist(d001_Ar_A, col=col_Gfos_A, cex.axis=2, xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ar), breaks=seq(0,maxhist_Ar,0.25))
hist(d001_Ar_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ar,0.25))
axis(1, labels=F, cex.axis=2)
abline(v=median(d001_Ar_A), col=white_transparent, lwd=2)
abline(v=median(d001_Ar_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(d001_Ar_B), col=white_transparent, lwd=2)
abline(v=median(d001_Ar_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(4,65, txt = lab[7], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

hist(d001_Ho_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ho), breaks=seq(0,maxhist_Ho,0.025))
hist(d001_Ho_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ho,0.025))
axis(1, labels=F, cex.axis=2)
axis(2, labels=F, cex.axis=2)
abline(v=median(d001_Ho_A), col=white_transparent, lwd=2)
abline(v=median(d001_Ho_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(d001_Ho_B), col=white_transparent, lwd=2)
abline(v=median(d001_Ho_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,65, txt = lab[7], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

hist(d01_Ar_A, col=col_Gfos_A, cex.axis=2, xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ar), breaks=seq(0,maxhist_Ar,0.25))
hist(d01_Ar_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ar,0.25))
axis(1, labels=T, cex.axis=2, padj=0.5)
abline(v=median(d01_Ar_A), col=white_transparent, lwd=2)
abline(v=median(d01_Ar_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(d01_Ar_B), col=white_transparent, lwd=2)
abline(v=median(d01_Ar_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(4,65, txt = lab[8], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

hist(d01_Ho_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ho), breaks=seq(0,maxhist_Ho,0.025))
hist(d01_Ho_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ho,0.025))
axis(1, labels=T, cex.axis=2, padj=0.5)
axis(2, labels=F, cex.axis=2)
abline(v=median(d01_Ho_A), col=white_transparent, lwd=2)
abline(v=median(d01_Ho_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(d01_Ho_B), col=white_transparent, lwd=2)
abline(v=median(d01_Ho_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,65, txt = lab[8], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

title(xlab = "Perpendicular offset",
      ylab = "Frequency",
      outer = TRUE, line = 3.5, cex.lab=2)
dev.off()

##** FIG S4 ####
##** Perpendicular offset histogram comparison for W
w00_Ar_A <- hist_Ar_A[,c(1,2,7,8,13,14)]
w05_Ar_A <- hist_Ar_A[,c(3,4,9,10,15,16)]
w10_Ar_A <- hist_Ar_A[,c(5,6,11,12,17,18)]
w00_Ar_B <- hist_Ar_B[,c(1,2,7,8,13,14)]
w05_Ar_B <- hist_Ar_B[,c(3,4,9,10,15,16)]
w10_Ar_B <- hist_Ar_B[,c(5,6,11,12,17,18)]
w00_Ho_A <- hist_Ho_A[,c(1,2,7,8,13,14)]
w05_Ho_A <- hist_Ho_A[,c(3,4,9,10,15,16)]
w10_Ho_A <- hist_Ho_A[,c(5,6,11,12,17,18)]
w00_Ho_B <- hist_Ho_B[,c(1,2,7,8,13,14)]
w05_Ho_B <- hist_Ho_B[,c(3,4,9,10,15,16)]
w10_Ho_B <- hist_Ho_B[,c(5,6,11,12,17,18)]

diffSPO_w00both_Ar_A <- c(apply(w00_Ar_A,2,sum)-apply(w05_Ar_A,2,sum),apply(w00_Ar_A,2,sum)-apply(w10_Ar_A,2,sum))
diffSPO_w00both_Ar_B <- c(apply(w00_Ar_B,2,sum)-apply(w05_Ar_B,2,sum),apply(w00_Ar_B,2,sum)-apply(w10_Ar_B,2,sum))
diffSPO_w00both_Ho_A <- c(apply(w00_Ho_A,2,sum)-apply(w05_Ho_A,2,sum),apply(w00_Ho_A,2,sum)-apply(w10_Ho_A,2,sum))
diffSPO_w00both_Ho_B <- c(apply(w00_Ho_B,2,sum)-apply(w05_Ho_B,2,sum),apply(w00_Ho_B,2,sum)-apply(w10_Ho_B,2,sum))
diffMPO_w00both_Ar_A <- c(apply(w00_Ar_A,2,median)-apply(w05_Ar_A,2,median),apply(w00_Ar_A,2,median)-apply(w10_Ar_A,2,median))
diffMPO_w00both_Ar_B <- c(apply(w00_Ar_B,2,median)-apply(w05_Ar_B,2,median),apply(w00_Ar_B,2,median)-apply(w10_Ar_B,2,median))
diffMPO_w00both_Ho_A <- c(apply(w00_Ho_A,2,median)-apply(w05_Ho_A,2,median),apply(w00_Ho_A,2,median)-apply(w10_Ho_A,2,median))
diffMPO_w00both_Ho_B <- c(apply(w00_Ho_B,2,median)-apply(w05_Ho_B,2,median),apply(w00_Ho_B,2,median)-apply(w10_Ho_B,2,median))

total_comparison_w00both <- length(c(diffMPO_w00both_Ar_A,diffMPO_w00both_Ar_B,diffMPO_w00both_Ho_A,diffMPO_w00both_Ho_B,diffSPO_w00both_Ar_A,diffSPO_w00both_Ar_B,diffSPO_w00both_Ho_A,diffSPO_w00both_Ho_B))
improved_fit_w00both <- sum(c(diffMPO_w00both_Ar_A,diffMPO_w00both_Ar_B,diffMPO_w00both_Ho_A,diffMPO_w00both_Ho_B,diffSPO_w00both_Ar_A,diffSPO_w00both_Ar_B,diffSPO_w00both_Ho_A,diffSPO_w00both_Ho_B)<0)
improved_fit_w00both/total_comparison_w00both

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS4.pdf"), width=8, height=6)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS4.png"), width=8, height=6, units="in", res=300)
}
op <- par(mfrow = c(3,2),
          oma = c(5,5,2,0) + 0.1,
          mar = c(0,0,2,1) + 0.1)
hist(w00_Ar_A, col=col_Gfos_A, cex.axis=2, xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ar), breaks=seq(0,maxhist_Ar,0.25))
hist(w00_Ar_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ar,0.25))
axis(1, labels=F, cex.axis=2)
mtext(side = 3, text = "Mean allelic richness", line = 1, adj=0.5, cex = 1.5)
abline(v=median(w00_Ar_A), col=white_transparent, lwd=2)
abline(v=median(w00_Ar_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(w00_Ar_B), col=white_transparent, lwd=2)
abline(v=median(w00_Ar_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(4,65, txt = lab[3], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)
legend(2,50, c(expression("Median" ~ italic(G. ~ fossarum) ~ "type A"),expression("Median" ~ italic(G. ~ fossarum) ~ "type B")),
       lty=2, col=c(col_Gfos_A,col_Gfos_B), bty="n", cex=1.3, lwd=2)

hist(w00_Ho_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ho), breaks=seq(0,maxhist_Ho,0.025))
hist(w00_Ho_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ho,0.025))
axis(1, labels=F, cex.axis=2)
axis(2, labels=F, cex.axis=2)
mtext(side = 3, text = "Mean observed heterozygosity", line = 1, adj=0.5, cex = 1.5)
abline(v=median(w00_Ho_A), col=white_transparent, lwd=2)
abline(v=median(w00_Ho_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(w00_Ho_B), col=white_transparent, lwd=2)
abline(v=median(w00_Ho_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,65, txt = lab[3], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

hist(w05_Ar_A, col=col_Gfos_A, cex.axis=2, xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ar), breaks=seq(0,maxhist_Ar,0.25))
hist(w05_Ar_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ar,0.25))
axis(1, labels=F, cex.axis=2)
abline(v=median(w05_Ar_A), col=white_transparent, lwd=2)
abline(v=median(w05_Ar_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(w05_Ar_B), col=white_transparent, lwd=2)
abline(v=median(w05_Ar_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(4,65, txt = lab[4], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

hist(w05_Ho_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ho), breaks=seq(0,maxhist_Ho,0.025))
hist(w05_Ho_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ho,0.025))
axis(1, labels=F, cex.axis=2)
axis(2, labels=F, cex.axis=2)
abline(v=median(w05_Ho_A), col=white_transparent, lwd=2)
abline(v=median(w05_Ho_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(w05_Ho_B), col=white_transparent, lwd=2)
abline(v=median(w05_Ho_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,65, txt = lab[4], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

hist(w10_Ar_A, col=col_Gfos_A, cex.axis=2, xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ar), breaks=seq(0,maxhist_Ar,0.25))
hist(w10_Ar_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ar,0.25))
axis(1, labels=T, cex.axis=2, padj=0.5)
abline(v=median(w10_Ar_A), col=white_transparent, lwd=2)
abline(v=median(w10_Ar_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(w10_Ar_B), col=white_transparent, lwd=2)
abline(v=median(w10_Ar_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(4,65, txt = lab[5], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

hist(w10_Ho_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ho), breaks=seq(0,maxhist_Ho,0.025))
hist(w10_Ho_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ho,0.025))
axis(1, labels=T, cex.axis=2, padj=0.5)
axis(2, labels=F, cex.axis=2)
abline(v=median(w10_Ho_A), col=white_transparent, lwd=2)
abline(v=median(w10_Ho_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(w10_Ho_B), col=white_transparent, lwd=2)
abline(v=median(w10_Ho_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,65, txt = lab[5], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)
strwidth(lab[5]) * 2
title(xlab = "Perpendicular offset",
      ylab = "Frequency",
      outer = TRUE, line = 3.5, cex.lab=2)
dev.off()

##** FIG S5 ####
##** Perpendicular offset histogram comparison for K
k0_Ar_A <- hist_Ar_A[,c(1,3,5,7,9,11,13,15,17)]
k1_Ar_A <- hist_Ar_A[,c(2,4,6,8,10,12,14,16,18)]
k0_Ar_B <- hist_Ar_B[,c(1,3,5,7,9,11,13,15,17)]
k1_Ar_B <- hist_Ar_B[,c(2,4,6,8,10,12,14,16,18)]
k0_Ho_A <- hist_Ho_A[,c(1,3,5,7,9,11,13,15,17)]
k1_Ho_A <- hist_Ho_A[,c(2,4,6,8,10,12,14,16,18)]
k0_Ho_B <- hist_Ho_B[,c(1,3,5,7,9,11,13,15,17)]
k1_Ho_B <- hist_Ho_B[,c(2,4,6,8,10,12,14,16,18)]

diffSPO_k0k1_Ar_A <- apply(k0_Ar_A,2,sum)-apply(k1_Ar_A,2,sum)
diffSPO_k0k1_Ar_B <- apply(k0_Ar_B,2,sum)-apply(k1_Ar_B,2,sum)
diffSPO_k0k1_Ho_A <- apply(k0_Ho_A,2,sum)-apply(k1_Ho_A,2,sum)
diffSPO_k0k1_Ho_B <- apply(k0_Ho_B,2,sum)-apply(k1_Ho_B,2,sum)
diffMPO_k0k1_Ar_A <- apply(k0_Ar_A,2,median)-apply(k1_Ar_A,2,median)
diffMPO_k0k1_Ar_B <- apply(k0_Ar_B,2,median)-apply(k1_Ar_B,2,median)
diffMPO_k0k1_Ho_A <- apply(k0_Ho_A,2,median)-apply(k1_Ho_A,2,median)
diffMPO_k0k1_Ho_B <- apply(k0_Ho_B,2,median)-apply(k1_Ho_B,2,median)

total_comparison_k0k1 <- length(c(diffMPO_k0k1_Ar_A,diffMPO_k0k1_Ar_B,diffMPO_k0k1_Ho_A,diffMPO_k0k1_Ho_B,diffSPO_k0k1_Ar_A,diffSPO_k0k1_Ar_B,diffSPO_k0k1_Ho_A,diffSPO_k0k1_Ho_B))
improved_fit_k0k1 <- sum(c(diffMPO_k0k1_Ar_A,diffMPO_k0k1_Ar_B,diffMPO_k0k1_Ho_A,diffMPO_k0k1_Ho_B,diffSPO_k0k1_Ar_A,diffSPO_k0k1_Ar_B,diffSPO_k0k1_Ho_A,diffSPO_k0k1_Ho_B)<0)
improved_fit_k0k1/total_comparison_k0k1

scal.fact <- 0.75
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS5.pdf"), width=8, height=6*scal.fact)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS5.png"), width=8, height=6*scal.fact, units="in", res=300)
}
op <- par(mfrow = c(3*scal.fact,2),
          oma = c(5*scal.fact,5*scal.fact,2*scal.fact,0) + 0.1,
          mar = c(0,0,2*scal.fact,1*scal.fact) + 0.1)
hist(k0_Ar_A, col=col_Gfos_A, tcl=-0.5*scal.fact, cex.axis=2*scal.fact, xaxt="n", yaxt="n", main="", ylim=c(0,90), xlim = c(0,maxhist_Ar), breaks=seq(0,maxhist_Ar,0.25))
hist(k0_Ar_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ar,0.25))
axis(1, labels=F, cex.axis=2*scal.fact, tcl=-0.5*scal.fact)
axis(2, labels=T, cex.axis=2*scal.fact, tcl=-0.5*scal.fact, padj=0.5)
mtext(side = 3, text = "Mean allelic richness", line = 1, adj=0.5, cex = 1.5)
abline(v=median(k0_Ar_A), col=white_transparent, lwd=2)
abline(v=median(k0_Ar_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(k0_Ar_B), col=white_transparent, lwd=2)
abline(v=median(k0_Ar_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(4,85, txt = lab[1], txt.adj=0.5, txt.cex = 2*scal.fact, frm.brd = NA, frm.col = white_transparent)
legend(2,70, c(expression("Median" ~ italic(G. ~ fossarum) ~ "type A"),expression("Median" ~ italic(G. ~ fossarum) ~ "type B")),
       lty=2, col=c(col_Gfos_A,col_Gfos_B), bty="n", cex=1, lwd=2)

hist(k0_Ho_A, col=col_Gfos_A, tcl=-0.5*scal.fact, cex.axis=2*scal.fact, xaxt="n", yaxt="n", main="", ylim=c(0,90), xlim = c(0,maxhist_Ho), breaks=seq(0,maxhist_Ho,0.025))
hist(k0_Ho_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ho,0.025))
axis(1, labels=F, cex.axis=2*scal.fact, tcl=-0.5*scal.fact)
axis(2, labels=F, cex.axis=2*scal.fact, tcl=-0.5*scal.fact)
mtext(side = 3, text = "Mean observed heterozygosity", line = 1, adj=0.5, cex = 1.5)
abline(v=median(k0_Ho_A), col=white_transparent, lwd=2)
abline(v=median(k0_Ho_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(k0_Ho_B), col=white_transparent, lwd=2)
abline(v=median(k0_Ho_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,85, txt = lab[1], txt.adj=0.5, txt.cex = 2*scal.fact, frm.brd = NA, frm.col = white_transparent)

hist(k1_Ar_A, col=col_Gfos_A, tcl=-0.5*scal.fact, cex.axis=2*scal.fact, xaxt="n", yaxt="n", main="", ylim=c(0,90), xlim = c(0,maxhist_Ar), breaks=seq(0,maxhist_Ar,0.25))
hist(k1_Ar_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ar,0.25))
axis(1, labels=T, cex.axis=2*scal.fact, tcl=-0.5*scal.fact, padj=0)
axis(2, labels=T, cex.axis=2*scal.fact, tcl=-0.5*scal.fact, padj=0.5)
abline(v=median(k1_Ar_A), col=white_transparent, lwd=2)
abline(v=median(k1_Ar_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(k1_Ar_B), col=white_transparent, lwd=2)
abline(v=median(k1_Ar_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(4,85, txt = lab[2], txt.adj=0.5, txt.cex = 2*scal.fact, frm.brd = NA, frm.col = white_transparent)

hist(k1_Ho_A, col=col_Gfos_A, tcl=-0.5*scal.fact, cex.axis=2*scal.fact, xaxt="n", yaxt="n", main="", ylim=c(0,90), xlim = c(0,maxhist_Ho), breaks=seq(0,maxhist_Ho,0.025))
hist(k1_Ho_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ho,0.025))
axis(1, labels=T, cex.axis=2*scal.fact, tcl=-0.5*scal.fact, padj=0)
axis(2, labels=F, cex.axis=2*scal.fact, tcl=-0.5*scal.fact)
abline(v=median(k1_Ho_A), col=white_transparent, lwd=2)
abline(v=median(k1_Ho_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(k1_Ho_B), col=white_transparent, lwd=2)
abline(v=median(k1_Ho_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,85, txt = lab[2], txt.adj=0.5, txt.cex = 2*scal.fact, frm.brd = NA, frm.col = white_transparent)

title(xlab = "Perpendicular offset",
      ylab = "Frequency",
      outer = TRUE, line = 3.5*scal.fact, cex.lab=2*scal.fact)
dev.off()

##** FIG S6 ####
##** Orthogonal distance to 1:1 line
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=pla$main)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS6.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS6.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=pla$main)

for (i in 2:ncol(Ar_modelled)){
  Ar_Mod_A <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosA]
  Ar_Mod_B <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosB]
  par(mar=c(0,0,0,0))
  plot(Ar_Mod_A~meanAr_A_red_updist, type="n",
       xlim=c(min(Ar_modelled[,-1]),max(Ar_modelled[,-1])),
       ylim=c(min(Ar_modelled[,-1]),max(Ar_modelled[,-1])),
       xaxt=ifelse(i%in%c(6,7,12,13,18,19), "s", "n"),
       yaxt=ifelse(i%in%c(3,5,7), "s", "n"),
       asp=1) # create empty plot
  if (i %in% c(14,16,18)){
    axis(4, labels=F)
  }
  if (i %in% c(2,4,6)){
    axis(2, labels=F)
  }
  if (i %in% c(2,3)){
    mtext(lab[6], side=3, line=0.8, cex=1.4)
  }
  if (i %in% c(8,9)){
    mtext(lab[7], side=3, line=0.8, cex=1.4)
  }
  if (i %in% c(14,15)){
    mtext(lab[8], side=3, line=0.8, cex=1.4)
  }
  if (i %in% c(4)){
    mtext(expression(bold("Model data:")*" Mean allelic richness"), side=2, line=1, cex=1)
  }
  if (i %in% c(12,13)){
    mtext(expression(bold("Empirical data:")*" Mean allelic richness"), side=1, line=3, cex=1)
  }
  abline(0,1, lwd=1, lty=2) # add 1:1 line
  
  points(Ar_Mod_A~meanAr_A_red_updist,col=col_Gfos_A, pch=16)
  points(Ar_Mod_B~meanAr_B_updist,col=col_Gfos_B, pch=16)
  
  for (j in 1:length(meanAr_A_red_updist)){
    point <- cbind(meanAr_A_red_updist,Ar_Mod_A)[j,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    segments(seg[1],seg[2],seg[3],seg[4], col=col_Gfos_A, lty=2, lwd=0.5)
  }
  for (k in 1:length(meanAr_B_updist)){
    point <- cbind(meanAr_B_updist,Ar_Mod_B)[k,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    segments(seg[1],seg[2],seg[3],seg[4], col=col_Gfos_B, lty=2, lwd=0.5)
  }
  textbox(max(Ar_modelled[,-1]),min(Ar_modelled[,-1])+2,paste0(measure1_short, " = ",formatC(round(sum_orthdist_Ar_A[[i-1]],sum_digits),digits=sum_digits, format="f")), txt.cex=sum_cex, txt.adj=1, txt.col=col_Gfos_A, frm.col=white_transparent, frm.brd = NA, frm.siz = 0.2)
  textbox(max(Ar_modelled[,-1]),min(Ar_modelled[,-1])+0.5,paste0(measure1_short, " = ",formatC(round(sum_orthdist_Ar_B[[i-1]],sum_digits),digits=sum_digits, format="f")), txt.cex=sum_cex, txt.adj=1, txt.col=col_Gfos_B, frm.col=white_transparent, frm.brd = NA, frm.siz = 0.2)
}

for (i in 1:length(lab)){
  plot.new()
  if (i %in% c(1:2)){
    text(0.5,1,lab[i],adj=c(0.5,1), cex=lab_cex[i])
  }
  if (i %in% c(3:5)){
    text(0,0.5,lab[i],adj=c(0,0.5), cex=lab_cex[i])
  }
}
plot.new()
legend("topleft",c(label_A,label_B),pch = 16, col = c(col_Gfos_A,col_Gfos_B), bty="n", cex=2)
text(0,0,measure1, adj=0)
dev.off()


##** FIG S7 ####
##** Orthogonal distance to 1:1 line
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=pla$main)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS7.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS7.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=pla$main)

for (i in 2:ncol(Ho_modelled)){
  Ho_Mod_A <- Ho_modelled[,i][rownames(Ho_modelled)%in%modsite_GfosA]
  Ho_Mod_B <- Ho_modelled[,i][rownames(Ho_modelled)%in%modsite_GfosB]
  par(mar=c(0,0,0,0))
  plot(Ho_Mod_A~meanHo_A_red_updist, type="n",
       xlim=c(min(Ho_modelled[,-1]),max(Ho_modelled[,-1])),
       ylim=c(min(Ho_modelled[,-1]),max(Ho_modelled[,-1])),
       xaxt=ifelse(i%in%c(6,7,12,13,18,19), "s", "n"),
       yaxt=ifelse(i%in%c(3,5,7), "s", "n"),
       asp=1) # create empty plot
  if (i %in% c(14,16,18)){
    axis(4, labels=F)
  }
  if (i %in% c(2,4,6)){
    axis(2, labels=F)
  }
  if (i %in% c(2,3)){
    mtext(lab[6], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(8,9)){
    mtext(lab[7], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(14,15)){
    mtext(lab[8], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(4)){
    mtext(expression(bold("Model data:")*" Mean observed heterozygosity"), side=2, line=1, cex=1)
  }
  if (i %in% c(12,13)){
    mtext(expression(bold("Empirical data:")*" Mean observed heterozygosity"), side=1, line=3, cex=1)
  }
  abline(0,1, lwd=1, lty=2) # add 1:1 line
  
  points(Ho_Mod_A~meanHo_A_red_updist,col=col_Gfos_A, pch=16)
  points(Ho_Mod_B~meanHo_B_updist,col=col_Gfos_B, pch=16)
  
  for (j in 1:length(meanHo_A_red_updist)){
    point <- cbind(meanHo_A_red_updist,Ho_Mod_A)[j,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    segments(seg[1],seg[2],seg[3],seg[4], col=col_Gfos_A, lty=2, lwd=0.5)
  }
  for (k in 1:length(meanHo_B_updist)){
    point <- cbind(meanHo_B_updist,Ho_Mod_B)[k,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    segments(seg[1],seg[2],seg[3],seg[4], col=col_Gfos_B, lty=2, lwd=0.5)
  }
  textbox(max(Ho_modelled[,-1]),min(Ho_modelled[,-1])+0.18,paste0(measure1_short," = ",formatC(round(sum_orthdist_Ho_A[[i-1]],sum_digits),digits=sum_digits, format="f")), txt.cex=sum_cex, txt.adj=1, txt.col=col_Gfos_A, frm.col=white_transparent, frm.brd = NA, frm.siz = 0.2)
  textbox(max(Ho_modelled[,-1]),min(Ho_modelled[,-1])+0.05,paste0(measure1_short," = ",formatC(round(sum_orthdist_Ho_B[[i-1]],sum_digits),digits=sum_digits, format="f")), txt.cex=sum_cex, txt.adj=1, txt.col=col_Gfos_B, frm.col=white_transparent, frm.brd = NA, frm.siz = 0.2)
}
for (i in 1:length(lab)){
  plot.new()
  if (i %in% c(1:2)){
    text(0.5,1,lab[i],adj=c(0.5,1), cex=lab_cex[i])
  }
  if (i %in% c(3:5)){
    text(0,0.5,lab[i],adj=c(0,0.5), cex=lab_cex[i])
  }
}
plot.new()
legend("topleft",c(label_A,label_B),pch = 16, col = c(col_Gfos_A,col_Gfos_B), bty="n", cex=2)
text(0,0,measure1, adj=0)
dev.off()

##** TABLE S8 ####
SPOrank_Ar_A <- labs_comb[order(apply(hist_Ar_A,2,sum))]
SPOrank_Ar_B <- labs_comb[order(apply(hist_Ar_B,2,sum))]

SPOrank_spo_Ar_A <- apply(hist_Ar_A,2,sum)[order(apply(hist_Ar_A,2,sum))]
SPOrank_d_Ar_A <- labs_short[,3][order(apply(hist_Ar_A,2,sum))]
SPOrank_W_Ar_A <- labs_short[,2][order(apply(hist_Ar_A,2,sum))]
SPOrank_K_Ar_A <- labs_short[,1][order(apply(hist_Ar_A,2,sum))]
SPOrank_spo_Ar_B <- apply(hist_Ar_B,2,sum)[order(apply(hist_Ar_B,2,sum))]
SPOrank_d_Ar_B <- labs_short[,3][order(apply(hist_Ar_B,2,sum))]
SPOrank_W_Ar_B <- labs_short[,2][order(apply(hist_Ar_B,2,sum))]
SPOrank_K_Ar_B <- labs_short[,1][order(apply(hist_Ar_B,2,sum))]

SPOrank_Ar <- data.frame("Ar_A_SPO"=SPOrank_spo_Ar_A,"Ar_A_d"=SPOrank_d_Ar_A,"Ar_A_W"=SPOrank_W_Ar_A,"Ar_A_K"=SPOrank_K_Ar_A,
                      "Ar_B_SPO"=SPOrank_spo_Ar_B,"Ar_B_d"=SPOrank_d_Ar_B,"Ar_B_W"=SPOrank_W_Ar_B,"Ar_B_K"=SPOrank_K_Ar_B)
write.csv2(SPOrank_Ar, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S8.csv"))

##** TABLE S9 ####
SPOrank_Ho_A <- labs_comb[order(apply(hist_Ho_A,2,sum))]
SPOrank_Ho_B <- labs_comb[order(apply(hist_Ho_B,2,sum))]

SPOrank_spo_Ho_A <- apply(hist_Ho_A,2,sum)[order(apply(hist_Ho_A,2,sum))]
SPOrank_d_Ho_A <- labs_short[,3][order(apply(hist_Ho_A,2,sum))]
SPOrank_W_Ho_A <- labs_short[,2][order(apply(hist_Ho_A,2,sum))]
SPOrank_K_Ho_A <- labs_short[,1][order(apply(hist_Ho_A,2,sum))]
SPOrank_spo_Ho_B <- apply(hist_Ho_B,2,sum)[order(apply(hist_Ho_B,2,sum))]
SPOrank_d_Ho_B <- labs_short[,3][order(apply(hist_Ho_B,2,sum))]
SPOrank_W_Ho_B <- labs_short[,2][order(apply(hist_Ho_B,2,sum))]
SPOrank_K_Ho_B <- labs_short[,1][order(apply(hist_Ho_B,2,sum))]

SPOrank_Ho <- data.frame("Ho_A_SPO"=SPOrank_spo_Ho_A,"Ho_A_d"=SPOrank_d_Ho_A,"Ho_A_W"=SPOrank_W_Ho_A,"Ho_A_K"=SPOrank_K_Ho_A,
                         "Ho_B_SPO"=SPOrank_spo_Ho_B,"Ho_B_d"=SPOrank_d_Ho_B,"Ho_B_W"=SPOrank_W_Ho_B,"Ho_B_K"=SPOrank_K_Ho_B)
write.csv2(SPOrank_Ho, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S9.csv"))

##** TABLE S10 ####
MPOrank_Ar_A <- labs_comb[order(apply(hist_Ar_A,2,median))]
MPOrank_Ar_B <- labs_comb[order(apply(hist_Ar_B,2,median))]

MPOrank_mpo_Ar_A <- apply(hist_Ar_A,2,median)[order(apply(hist_Ar_A,2,median))]
MPOrank_d_Ar_A <- labs_short[,3][order(apply(hist_Ar_A,2,median))]
MPOrank_W_Ar_A <- labs_short[,2][order(apply(hist_Ar_A,2,median))]
MPOrank_K_Ar_A <- labs_short[,1][order(apply(hist_Ar_A,2,median))]
MPOrank_mpo_Ar_B <- apply(hist_Ar_B,2,median)[order(apply(hist_Ar_B,2,median))]
MPOrank_d_Ar_B <- labs_short[,3][order(apply(hist_Ar_B,2,median))]
MPOrank_W_Ar_B <- labs_short[,2][order(apply(hist_Ar_B,2,median))]
MPOrank_K_Ar_B <- labs_short[,1][order(apply(hist_Ar_B,2,median))]

MPOrank_Ar <- data.frame("Ar_A_MPO"=MPOrank_mpo_Ar_A,"Ar_A_d"=MPOrank_d_Ar_A,"Ar_A_W"=MPOrank_W_Ar_A,"Ar_A_K"=MPOrank_K_Ar_A,
                         "Ar_B_MPO"=MPOrank_mpo_Ar_B,"Ar_B_d"=MPOrank_d_Ar_B,"Ar_B_W"=MPOrank_W_Ar_B,"Ar_B_K"=MPOrank_K_Ar_B)
write.csv2(MPOrank_Ar, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S10.csv"))

##** TABLE S11 ####
MPOrank_Ho_A <- labs_comb[order(apply(hist_Ho_A,2,median))]
MPOrank_Ho_B <- labs_comb[order(apply(hist_Ho_B,2,median))]

MPOrank_mpo_Ho_A <- apply(hist_Ho_A,2,median)[order(apply(hist_Ho_A,2,median))]
MPOrank_d_Ho_A <- labs_short[,3][order(apply(hist_Ho_A,2,median))]
MPOrank_W_Ho_A <- labs_short[,2][order(apply(hist_Ho_A,2,median))]
MPOrank_K_Ho_A <- labs_short[,1][order(apply(hist_Ho_A,2,median))]
MPOrank_mpo_Ho_B <- apply(hist_Ho_B,2,median)[order(apply(hist_Ho_B,2,median))]
MPOrank_d_Ho_B <- labs_short[,3][order(apply(hist_Ho_B,2,median))]
MPOrank_W_Ho_B <- labs_short[,2][order(apply(hist_Ho_B,2,median))]
MPOrank_K_Ho_B <- labs_short[,1][order(apply(hist_Ho_B,2,median))]

MPOrank_Ho <- data.frame("Ho_A_MPO"=MPOrank_mpo_Ho_A,"Ho_A_d"=MPOrank_d_Ho_A,"Ho_A_W"=MPOrank_W_Ho_A,"Ho_A_K"=MPOrank_K_Ho_A,
                         "Ho_B_MPO"=MPOrank_mpo_Ho_B,"Ho_B_d"=MPOrank_d_Ho_B,"Ho_B_W"=MPOrank_W_Ho_B,"Ho_B_K"=MPOrank_K_Ho_B)
write.csv2(MPOrank_Ho, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S11.csv"))

##** FIG S12 ####
##** Histogram of directed perpendicular offset to 1:1 line
main_s12=F
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s12, h.main=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS12.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS12.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s12, h.main=0.5)

min(abs(apply(hist_Ar_A_directed, 2, median)))

for (i in 1:ncol(hist_Ar_A)){
  par(mar=c(0,0,0,0))
  yclip <- 30
  ylim <- 35
  hist(hist_Ar_A_directed[,i],
       breaks=seq(floor(min(hist_Ar_A_directed,hist_Ar_B_directed)),ceiling(max(hist_Ar_A_directed,hist_Ar_B_directed))-0.5,0.5),
       xlim=c(floor(min(hist_Ar_A_directed,hist_Ar_B_directed)),ceiling(max(hist_Ar_A_directed,hist_Ar_B_directed))),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="")
  clip(floor(min(hist_Ar_A_directed,hist_Ar_B_directed)),ceiling(max(hist_Ar_A_directed,hist_Ar_B_directed)), 0, yclip)
  abline(v=0, lwd=1.5, lty=2)
  abline(v=median(hist_Ar_A_directed[,i]),col=col_Gfos_A)
  abline(v=median(hist_Ar_B_directed[,i]),col=col_Gfos_B)
  clip(floor(min(hist_Ar_A_directed,hist_Ar_B_directed)),ceiling(max(hist_Ar_A_directed,hist_Ar_B_directed)), 0, ylim)
  textbox(ceiling(max(hist_Ar_A_directed,hist_Ar_B_directed)),30,paste0(measure4," = ",formatC(round(median(hist_Ar_A_directed[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_A, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  textbox(ceiling(max(hist_Ar_A_directed,hist_Ar_B_directed)),25,paste0(measure4," = ",formatC(round(median(hist_Ar_B_directed[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_B, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  hist(hist_Ar_A_directed[,i],
       breaks=seq(floor(min(hist_Ar_A_directed,hist_Ar_B_directed)),ceiling(max(hist_Ar_A_directed,hist_Ar_B_directed))-0.5,0.5),
       xlim=c(floor(min(hist_Ar_A_directed,hist_Ar_B_directed)),ceiling(max(hist_Ar_A_directed,hist_Ar_B_directed))),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="", add=T)
  hist(hist_Ar_B_directed[,i],
       breaks=seq(floor(min(hist_Ar_A_directed,hist_Ar_B_directed)),ceiling(max(hist_Ar_A_directed,hist_Ar_B_directed))-0.5,0.5),
       xlim=c(floor(min(hist_Ar_A_directed,hist_Ar_B_directed)),ceiling(max(hist_Ar_A_directed,hist_Ar_B_directed))),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_B,
       main="",
       add=T)
  
  if (i %in% c(1,2)){
    mtext(lab[6], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(7,8)){
    mtext(lab[7], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(13,14)){
    mtext(lab[8], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(3)){
    mtext("Frequencies [counts]", side=2, line=1, cex=1)
  }
  if (i %in% c(11,12)){
    mtext(measure3, side=1, line=3, cex=1)
  }
  if (i %in% c(5:6,11:12,17:18)){
    axis(1)
  }
  if (i %in% c(1,3,5)){
    axis(2, at=seq(0,yclip,5), labels=F)
  }
  if (i %in% c(2,4,6)){
    axis(2, at=seq(0,yclip,5))
  }
  if (i %in% c(13:18)){
    axis(4, at=seq(0,yclip,5), labels=F)
  }
}

if(main_s12){
  plot.new()
  text(0.5,0.7,paste0("Mean allelic richness: Distribution of directed ",measure3a),adj=c(0.5,0.5),cex=3)
}

for (i in 1:length(lab)){
  plot.new()
  if (i %in% c(1:2)){
    text(0.5,1,lab[i],adj=c(0.5,1), cex=lab_cex[i])
  }
  if (i %in% c(3:5)){
    text(0,0.5,lab[i],adj=c(0,0.5), cex=lab_cex[i])
  }
}
par(mar=c(0,6,16,0.5), pty="s")
i <- ncol(Ar_modelled)
Ar_Mod_A <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosA]
Ar_Mod_B <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosB]
plot(Ar_Mod_A~meanAr_A_red_updist, type="n",
     xlim=c(min(Ar_modelled[,-1]),max(Ar_modelled[,-1])),
     ylim=c(min(Ar_modelled[,-1]),max(Ar_modelled[,-1])),
     xaxt="n",
     yaxt="n",
     xlab="",
     ylab="",
     bty="n",
     asp=1)

par(xpd=T)
legend(-5,yclip,c(label_A,label_B),pch = 16, col = c(col_Gfos_A,col_Gfos_B), bty="n", cex=2)
par(xpd=F)
dev.off()

##** FIG S13 ####
##** Histogram of directed perpendicular offset to 1:1 line
main_s13=F
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s13, h.main=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS13.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS13.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s13, h.main=0.5)

for (i in 1:ncol(hist_Ho_A_directed)){
  par(mar=c(0,0,0,0))
  yclip <- 30
  ylim <- 35
  xlim <- round(max(abs(min(hist_Ho_A_directed,hist_Ho_B_directed)),max(hist_Ho_A_directed,hist_Ho_B_directed)),1)
  hist(hist_Ho_A_directed[,i],
       breaks=seq(-xlim,xlim,0.05),
       xlim=c(-xlim,xlim),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="")
  clip(-xlim,xlim, 0, yclip)
  abline(v=0, lwd=1.5, lty=2)
  abline(v=median(hist_Ho_A_directed[,i]),col=col_Gfos_A)
  abline(v=median(hist_Ho_B_directed[,i]),col=col_Gfos_B)
  clip(-xlim,xlim, 0, ylim)
  textbox(xlim,30,paste0(measure4," = ",formatC(round(median(hist_Ho_A_directed[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_A, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  textbox(xlim,25,paste0(measure4," = ",formatC(round(median(hist_Ho_B_directed[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_B, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  hist(hist_Ho_A_directed[,i],
       breaks=seq(-xlim,xlim,0.05),
       xlim=c(-xlim,xlim),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="", add=T)
  hist(hist_Ho_B_directed[,i],
       breaks=seq(-xlim,xlim,0.05),
       xlim=c(-xlim,xlim),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_B,
       main="",
       add=T)
  
  if (i %in% c(1,2)){
    mtext(lab[6], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(7,8)){
    mtext(lab[7], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(13,14)){
    mtext(lab[8], side=3, line=0.8, cex=1.3)
  }
  if (i %in% c(3)){
    mtext("Frequencies [counts]", side=2, line=1, cex=1)
  }
  if (i %in% c(11,12)){
    mtext(measure3, side=1, line=3, cex=1)
  }
  if (i %in% c(5:6,11:12,17:18)){
    axis(1, at=as.numeric(formatC(seq(-xlim+0.1,xlim-0.1,0.1),digits=1, format="f")), labels=T)
  }
  if (i %in% c(1,3,5)){
    axis(2, at=seq(0,yclip,5), labels=F)
  }
  if (i %in% c(2,4,6)){
    axis(2, at=seq(0,yclip,5))
  }
  if (i %in% c(13:18)){
    axis(4, at=seq(0,yclip,5), labels=F)
  }
}

if(main_s13){
  plot.new()
  text(0.5,0.7,paste0("Mean observed heterozygosity: Distribution of directed ",measure3a),adj=c(0.5,0.5),cex=3)
}

for (i in 1:length(lab)){
  plot.new()
  if (i %in% c(1:2)){
    text(0.5,1,lab[i],adj=c(0.5,1), cex=lab_cex[i])
  }
  if (i %in% c(3:5)){
    text(0,0.5,lab[i],adj=c(0,0.5), cex=lab_cex[i])
  }
}
par(mar=c(0,6,16,0.5), pty="s")
i <- ncol(Ar_modelled)
Ar_Mod_A <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosA]
Ar_Mod_B <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosB]
plot(Ar_Mod_A~meanAr_A_red_updist, type="n",
     xlim=c(min(Ar_modelled[,-1]),max(Ar_modelled[,-1])),
     ylim=c(min(Ar_modelled[,-1]),max(Ar_modelled[,-1])),
     xaxt="n",
     yaxt="n",
     xlab="",
     ylab="",
     bty="n",
     asp=1)

par(xpd=T)
legend(-5,yclip,c(label_A,label_B),pch = 16, col = c(col_Gfos_A,col_Gfos_B), bty="n", cex=2)
par(xpd=F)
dev.off()

##** FIG S14 ####
##** Smoothed scatterplot Ar by updist
main_s14=F
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS14.pdf"), width=14, height=8)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS14.png"), width=14, height=8, units="in", res=300)
}
if(main_s14){
  par(mfrow=c(1,1), mar=c(4.6, 4.6, 1.6, 14), xpd=TRUE)
  plot(Ar_modelled_smooth[,2]~updist_Mod_plot, col=colMod, type="n",
       ylab="Mean allelic richness", xlab="Upstream distance [km]",
       xlim=rev(range(c(updist_B_plot,updist_A_red_plot,updist_Mod_plot))), ylim=range(Ar_modelled_smooth[,-1]), pch=16, bty = "l", xaxt="n",
       main=paste0("Mean allelic richness by upstream distance"), cex.main=2, cex.axis=1.5, cex.lab=1.5)
}else{
  par(mfrow=c(1,1), mar=c(4.6, 4.6, 0.6, 14), xpd=TRUE)
  plot(Ar_modelled_smooth[,2]~updist_Mod_plot, col=colMod, type="n",
       ylab="Mean allelic richness", xlab="Upstream distance [km]",
       xlim=rev(range(c(updist_B_plot,updist_A_red_plot,updist_Mod_plot))), ylim=range(Ar_modelled_smooth[,-1]), pch=16, bty = "l", xaxt="n",
       main="", cex.axis=1.5, cex.lab=1.5)
}

for (i in 2:ncol(Ar_modelled_smooth)){
  lines(Ar_modelled_smooth[,i], x=updist_Mod_plot, col=loess_pal[i-1], lwd=loess_lwd[i-1], lty=loess_lty[i-1])
}
lines(smoothedAr_A_red_updist, x=updist_A_red_plot, col=paste0(col_Gfos_A,alpha2), lwd=9)
lines(smoothedAr_A_red_updist, x=updist_A_red_plot, col=col_Gfos_A, lwd=3)
lines(smoothedAr_B_updist, x=updist_B_plot, col=paste0(col_Gfos_B,alpha2), lwd=9)
lines(smoothedAr_B_updist, x=updist_B_plot, col=col_Gfos_B, lwd=3)
axis(1, c(0,50,100,150,200,250,300), at=c(0,50000,100000,150000,200000,250000,300000), cex.axis=1.5, cex.lab=1.5)
legend("topright", cex=1.2, inset=c(-0.27,0), c(label_A,label_B,"",labs_comb[1:6],"",labs_comb[7:12],"",labs_comb[13:18]),lty = c(1,1,0,loess_lty[1:6],0,loess_lty[1:6],0,loess_lty[1:6]), col = c(col_Gfos_A,col_Gfos_B,0,loess_pal[1:6],0,loess_pal[7:12],0,loess_pal[13:18]), lwd = c(6,6,0,loess_lwd[1:6],0,loess_lwd[7:12],0,loess_lwd[13:18]), bty="n")
dev.off()

##** FIG S15 ####
##** Smoothed scatterplot Ho by updist
main_s15=F
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS15.pdf"), width=14, height=8)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS15.png"), width=14, height=8, units="in", res=300)
}
if(main_s15){
  par(mfrow=c(1,1), mar=c(4.6, 4.6, 1.6, 14), xpd=TRUE)
  plot(Ho_modelled_smooth[,2]~updist_Mod_plot, col=colMod, type="n",
       ylab="Mean observed heterozygosity", xlab="Upstream distance [km]",
       xlim=rev(range(c(updist_B_plot,updist_A_red_plot,updist_Mod_plot))), ylim=range(Ho_modelled_smooth[,-1]), pch=16, bty = "l", xaxt="n",
       main=paste0("Mean observed heterozygosity by upstream distance"), cex.main=2, cex.axis=1.5, cex.lab=1.5)
}else{
  par(mfrow=c(1,1), mar=c(4.6, 4.6, 0.6, 14), xpd=TRUE)
  plot(Ho_modelled_smooth[,2]~updist_Mod_plot, col=colMod, type="n",
       ylab="Mean observed heterozygosity", xlab="Upstream distance [km]",
       xlim=rev(range(c(updist_B_plot,updist_A_red_plot,updist_Mod_plot))), ylim=range(Ho_modelled_smooth[,-1]), pch=16, bty = "l", xaxt="n",
       main="", cex.axis=1.5, cex.lab=1.5)
}

for (i in 2:ncol(Ho_modelled_smooth)){
  lines(Ho_modelled_smooth[,i], x=updist_Mod_plot, col=loess_pal[i-1], lwd=loess_lwd[i-1], lty=loess_lty[i-1])
}
lines(smoothedHo_A_red_updist, x=updist_A_red_plot, col=paste0(col_Gfos_A,alpha2), lwd=9)
lines(smoothedHo_A_red_updist, x=updist_A_red_plot, col=col_Gfos_A, lwd=3)
lines(smoothedHo_B_updist, x=updist_B_plot, col=paste0(col_Gfos_B,alpha2), lwd=9)
lines(smoothedHo_B_updist, x=updist_B_plot, col=col_Gfos_B, lwd=3)
axis(1, c(0,50,100,150,200,250,300), at=c(0,50000,100000,150000,200000,250000,300000), cex.axis=1.5, cex.lab=1.5)
legend("topright", cex=1.2, inset=c(-0.27,0), c(label_A,label_B,"",labs_comb[1:6],"",labs_comb[7:12],"",labs_comb[13:18]),lty = c(1,1,0,loess_lty[1:6],0,loess_lty[1:6],0,loess_lty[1:6]), col = c(col_Gfos_A,col_Gfos_B,0,loess_pal[1:6],0,loess_pal[7:12],0,loess_pal[13:18]), lwd = c(6,6,0,loess_lwd[1:6],0,loess_lwd[7:12],0,loess_lwd[13:18]), bty="n")
dev.off()