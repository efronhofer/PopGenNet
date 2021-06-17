# PopGenNet analysis
# 2021-06-11
# Author: Roman Alther, Eawag, Duebendorf, Switzerland
# Works in R ver. 3.6.1 and 4.0.3 (tested on GNU/Linux, MacOS, Windows) with RStudio ver. 1.3.1093
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#

#### PREPARATION -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
rm(list=ls()) # start from scratch
existing_data=T # Should the analysis run with pre-collated and attached data from the authors? Otherwise the required raw data need to be prepared using the script 'PopGenNet_Prep_20201210.R', creating an output folder in '02_Data_prep'.
internal=F # defaults to FALSE, but was set to TRUE for publication preparation (lab internal use)
log_trans=T # log-transform some explanatory variables ("betw_undir","betw_dir","degree","catch")
critval <- 1.96 ## approx 95% confidence interval (plus/minus x times SD)

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
WD <- getwd() # save current directory (should correspond to 03_Analysis, source script from there)

#* Load data ####
setwd("..")
DF <- getwd()
if (internal){
  prep_folder <- "Output20210528_wo_gf21"
  rdata <- "PopGenNet_data20210528_wo_gf21"
}else{
  prep_folder <- "Output"
  rdata <- "PopGenNet_data"
}
load(paste0("02_Data_prep/",rdata,".RData"))
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

# Load custom packages and check for updates
if (internal){
  # Load lab internal package for publication figure preparation
  library(SwissRiverPlot) # to plot maps of Switzerland, works with version 0.6-7
  update_SRP()
}else{
  library(OpenSwissRiverPlot) # to plot maps of Switzerland, works with version 0.4-0
  update_OSRP()
}
library(MultiPanel) # to plot multipanel figures, works with version 0.7-0
update_MultiPanel()

#* Figures preparation ####

##*** Figure format ####
pdf=T # set to TRUE if figures should be prepared as PDF

##*** Size ####
fig.width=12 # standard figure width in inches
fig.height=9 # standard figure height in inches

##*** Smoothing span ####
# loessspan <- 0.50 # the parameter ?? which controls the degree of smoothing in loess()

##*** Labels ####
label_A <- expression(italic(G.~fossarum)~"type A") # italic G. fossarum A
label_B <- expression(italic(G.~fossarum)~"type B") # italic G. fossarum B
label_mod <- "Simulation data"

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
lab_Hs <- "Expecteded heterozygosity"
lab_Ar_short <- "Mean Ar"
lab_Ho_short <- "Mean Ho"
lab_Hs_short <- "He"
measure1 <- "SPO = Sum of perpendicular offsets"
measure1_short <- "SPO"
measure2 <- "MPO = Median of perpendicular offsets"
measure2_short <- "MPO"
measure3 <- "Perpendicular offset"
measure3a <- "perpendicular offsets"
measure4 <- "DMPO"
lab_sub <- c("(a)","(b)","(c)","(d)","(e)","(f)")
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
# loess_pal_short <- c("#54278f","#807dba","#bcbddc")
# loess_pal <- rep(loess_pal_short, each=length(W)*length(K))
# loess_lty <- rep(c(1,1,2,2,3,3), length(D))
# loess_lwd <- rep(c(2,4), length(D)*length(W))

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

# The function requires a model formulation, the data, a name for the output, the model family, a significance level, and if selection should be implemented
lm.bind <- function(model, data, name, sign, step=F){
  assign(paste0("lm_",name),lm(formula(model), data))
  b <- get(paste0("lm_",name))
  if (step){
    assign(paste0("slm_",name), step(get(paste0("lm_",name)))) # backward selection
    c <- get(paste0("slm_",name))
    assign(paste0("Predict_slm_",name),predict(get(paste0("slm_",name)), type="response", se.fit=T))
    assign(paste0("slm_",name,"_upr"),get(paste0("Predict_slm_",name))$fit + (sign * get(paste0("Predict_slm_",name))$se.fit))
    assign(paste0("slm_",name,"_lwr"),get(paste0("Predict_slm_",name))$fit - (sign * get(paste0("Predict_slm_",name))$se.fit))
    assign(paste0("slm_",name,"_fit"),get(paste0("Predict_slm_",name))$fit)
    data <- cbind(data, cbind(get(paste0("slm_",name,"_fit")), get(paste0("slm_",name,"_lwr")), get(paste0("slm_",name,"_upr"))))
    colnames(data)[c(ncol(data)-2,ncol(data)-1,ncol(data))] <- c(paste0("slm_",name,"_fit"),paste0("slm_",name,"_lwr"),paste0("slm_",name,"_upr"))
  }else{
    assign(paste0("Predict_lm_",name),predict(get(paste0("lm_",name)), type="response", se.fit=T))
    assign(paste0("lm_",name,"_upr"),get(paste0("Predict_lm_",name))$fit + (sign * get(paste0("Predict_lm_",name))$se.fit))
    assign(paste0("lm_",name,"_lwr"),get(paste0("Predict_lm_",name))$fit - (sign * get(paste0("Predict_lm_",name))$se.fit))
    assign(paste0("lm_",name,"_fit"),get(paste0("Predict_lm_",name))$fit)
    data <- cbind(data, cbind(get(paste0("lm_",name,"_fit")), get(paste0("lm_",name,"_lwr")), get(paste0("lm_",name,"_upr"))))
    colnames(data)[c(ncol(data)-2,ncol(data)-1,ncol(data))] <- c(paste0("lm_",name,"_fit"),paste0("lm_",name,"_lwr"),paste0("lm_",name,"_upr"))
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
GLMplot <- function(x,y,dat,model,xlabel,ylabel,ylim=NULL,axislog="", col1=col_Gfos_A, col2=col_Gfos_B, pt.cex=1, CI_border=T, trans=0.3, xrev=F, xax=NULL, yax=NULL, pointtrans=F, cex.lab=2, cex.axis=1.5, legend=T, cex.legend=1.5){
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
       ylim=ylim,
       log=axislog,
       xaxt=xax, yaxt=yax, cex.lab=cex.lab, cex.axis=cex.axis)
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

#### Function to calculate decimal variant of ceiling()
# As proposed by user Ferroao (https://stackoverflow.com/users/6388753/ferroao) here: https://stackoverflow.com/a/59861612/6380381
ceiling_dec <- function(x, decimals=1) {
  x2<-x*10^decimals
  ceiling(x2)/10^decimals
}

#### Function to specify decimals
# Copied from https://stackoverflow.com/a/12135122/6380381
specify_decimal <- function(x, k, formatC=TRUE){
  if (formatC){
    trimws(formatC(round(x, digits=k),digits=k, format="f"))
  }else{
    trimws(format(round(x, k), nsmall=k))
  }
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
  Ho_modelled <- data.frame(matrix(nrow=length(empiricaldata), ncol=length(D)*length(W)*length(K)+1))
  rownames(Ho_modelled) <- adj.mat.names[empiricaldata]
  Hs_modelled <- data.frame(matrix(nrow=length(empiricaldata), ncol=length(D)*length(W)*length(K)+1))
  rownames(Hs_modelled) <- adj.mat.names[empiricaldata]
  
  r <- 0
  for (d in 1:length(D)){ # looping over dispersal rates 
    for (w in 1:length(W)){ # looping over dispersal directionalities
      for (k in 1:length(K)){ # looping over carrying capacities
        r <- r+1
        
        label_Mod <- paste0("D=",D_label[d],", W_up=",W_label[w],", K=",K_label[k])
        
        load(paste0(DF,"/02_Data_prep/",prep_folder,"/IndPopGenData_",D[d],"_",W[w],"_",K[k],".Rdata"))
        
        orderMod <- order(updist_Mod, decreasing = T)
        
        if (r==1){
          match_Mod <- match_Mod[orderMod]
        }
        
        #*** Mean Ar by upstream distance ####
        updist_Mod_plot <- updist_Mod[orderMod]
        meanAr_Mod_updist <- meanAr_Mod[orderMod]
        orderA_red <- order(updist_A_red, decreasing = T)
        updist_A_red_plot <- updist_A_red[orderA_red]
        meanAr_A_red_updist <- meanAr_A_red[orderA_red]
        orderB <- order(updist_B, decreasing = T)
        updist_B_plot <- updist_B[orderB]
        meanAr_B_updist <- meanAr_B[orderB]
        
        #*** Mean Ho by upstream distance ####
        updist_Mod_plot <- updist_Mod[orderMod]
        meanHo_Mod_updist <- meanHo_Mod[orderMod]
        orderA_red <- order(updist_A_red, decreasing = T)
        updist_A_red_plot <- updist_A_red[orderA_red]
        meanHo_A_red_updist <- meanHo_A_red[orderA_red]
        orderB <- order(updist_B, decreasing = T)
        updist_B_plot <- updist_B[orderB]
        meanHo_B_updist <- meanHo_B[orderB]
        
        #*** He by upstream distance ####
        updist_Mod_plot <- updist_Mod[orderMod]
        meanHs_Mod_updist <- meanHs_Mod[orderMod]
        orderA_red <- order(updist_A_red, decreasing = T)
        updist_A_red_plot <- updist_A_red[orderA_red]
        meanHs_A_red_updist <- meanHs_A_red[orderA_red]
        orderB <- order(updist_B, decreasing = T)
        updist_B_plot <- updist_B[orderB]
        meanHs_B_updist <- meanHs_B[orderB]
        
        # save modelled popgen values to table
        if (r==1){
          Ar_modelled[,1] <- updist_Mod_plot
          colnames(Ar_modelled)[1] <- "upstream_distance"
        }
        Ar_modelled[,1+r] <- meanAr_Mod_updist
        colnames(Ar_modelled)[r+1] <- paste0(D[d],"_",W[w],"_",K[k])
        
        # save modelled popgen values to table
        if (r==1){
          Ho_modelled[,1] <- updist_Mod_plot
          colnames(Ho_modelled)[1] <- "upstream_distance"
        }
        Ho_modelled[,1+r] <- meanHo_Mod_updist
        colnames(Ho_modelled)[r+1] <- paste0(D[d],"_",W[w],"_",K[k])
        
        # save modelled popgen values to table
        if (r==1){
          Hs_modelled[,1] <- updist_Mod_plot
          colnames(Hs_modelled)[1] <- "upstream_distance"
        }
        Hs_modelled[,1+r] <- meanHs_Mod_updist
        colnames(Hs_modelled)[r+1] <- paste0(D[d],"_",W[w],"_",K[k])
        
        # setwd(WD)
        
      } # end looping over carrying capacities
    } # end looping over dispersal directionalities
  } # end looping over dispersal rates
}else{
  load("PopGenData.RData")
}

#### EXPLANATORY VARIABLES -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ ####
# check for correlations
expl.var1 <- DATA[,!(colnames(DATA) %in% c("meanAr","meanHo","meanHs","network_match","spec","log_updist","log_catch","clos_undir","clos_dir"))]
cor.var1 <- cor(expl.var1, method = "kendall")

# transform explanatory variables
if(log_trans){
  DATA[,c("betw_undir","betw_dir","degree","catch")] <- log(DATA[c("betw_undir","betw_dir","degree","catch")])
  DATA[DATA==-Inf] <- 0
}

expl.var <- DATA[,!(colnames(DATA) %in% c("meanAr","meanHo","meanHs","degree","betw_undir","clos_dir","network_match","spec","log_updist","log_catch","clos_undir","std_clos_dir"))]
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
#### Initial model exploration
shapiro.test(DATA$meanAr)

# Interaction model
int.mod.Ar <- lm(meanAr ~ updist*betw_dir*std_clos_undir*catch*spec, DATA, na.action = "na.fail")
summary(int.mod.Ar)
opar <- par(mfrow=c(2,2))
plot(int.mod.Ar)
par(opar)
step(int.mod.Ar)
# best model keeps some interactions
if(log_trans){
  int.sel.mod.Ar <- lm(meanAr ~ updist + betw_dir + std_clos_undir + catch + 
                         spec + updist:betw_dir + updist:std_clos_undir + betw_dir:std_clos_undir + 
                         updist:catch + betw_dir:catch + std_clos_undir:catch + updist:spec + 
                         betw_dir:spec + std_clos_undir:spec + catch:spec + updist:betw_dir:std_clos_undir + 
                         updist:betw_dir:catch + updist:std_clos_undir:catch + betw_dir:std_clos_undir:catch + 
                         updist:betw_dir:spec + updist:catch:spec + betw_dir:catch:spec + 
                         updist:betw_dir:std_clos_undir:catch + updist:betw_dir:catch:spec, 
                       data = DATA, na.action = "na.fail")
}else{
  int.sel.mod.Ar <- lm(meanAr ~ updist + betw_dir + std_clos_undir + catch + 
                     spec + updist:betw_dir + updist:std_clos_undir + betw_dir:std_clos_undir + 
                     updist:catch + betw_dir:catch + std_clos_undir:catch + updist:spec + 
                     betw_dir:spec + std_clos_undir:spec + catch:spec + updist:betw_dir:std_clos_undir + 
                     updist:betw_dir:catch + betw_dir:std_clos_undir:catch + updist:betw_dir:spec + 
                     updist:std_clos_undir:spec + betw_dir:std_clos_undir:spec + 
                     updist:catch:spec + betw_dir:catch:spec + std_clos_undir:catch:spec + 
                     updist:betw_dir:catch:spec + betw_dir:std_clos_undir:catch:spec, 
                   data = DATA, na.action = "na.fail")
}

# Linear model
lin.mod.Ar <- lm(meanAr ~ updist+betw_dir+std_clos_undir+catch+spec, DATA, na.action = "na.fail")
summary(lin.mod.Ar)
opar <- par(mfrow=c(2,2))
plot(lin.mod.Ar)
par(opar)
MuMIn::dredge(lin.mod.Ar)
step(lin.mod.Ar)
# best model: ctc+upd (untransformed); btw_dir+spc+upd (log-transformed)
if(log_trans){
  sel.mod.Ar <- lm(meanAr ~ updist+betw_dir+std_clos_undir, DATA, na.action = "na.fail")
}else{
  sel.mod.Ar <- lm(meanAr ~ updist+catch, DATA, na.action = "na.fail")
}  

# Comparison
AIC(int.sel.mod.Ar)
AIC(sel.mod.Ar)
AIC(glm(meanAr ~ updist+catch+spec, DATA, family="Gamma"))
# linear model outperforms interaction model
summary(sel.mod.Ar)
opar <- par(mfrow=c(2,2))
plot(sel.mod.Ar)
par(opar)
car::vif(sel.mod.Ar)

#### Full LM of allelic richness without interactions
model <- lm.bind(meanAr ~ updist+betw_dir+std_clos_undir+catch+spec, DATA, "Ar_full", critval, step=T)
DATA <- model[[1]]
lm_Ar_full <- model[[2]]
summary(lm_Ar_full)
car::vif(lm_Ar_full)
slm_Ar_full <- model[[3]]
summary(slm_Ar_full)
car::vif(slm_Ar_full) # should be the same as above (car::vif(sel.mod.Ar))

#### LM of allelic richness by upstream distance * species
AIC(lm(meanAr ~ updist+spec, DATA))
AIC(lm(meanAr ~ updist, DATA))
model <- lm.bind(meanAr ~ updist, DATA, "Ar_updist", critval, step=T)
DATA <- model[[1]]
lm_Ar_updist <- model[[2]]
slm_Ar_updist <- model[[3]]

#### LM of allelic richness by undirected closeness centrality
AIC(lm(meanAr ~ std_clos_undir+spec, DATA))
AIC(lm(meanAr ~ std_clos_undir, DATA))
model <- lm.bind(meanAr ~ std_clos_undir, DATA, "Ar_clos", critval, step=T)
DATA <- model[[1]]
lm_Ar_clos <- model[[2]]
slm_Ar_clos <- model[[3]]

#### LM of allelic richness by directed betweenness centrality
AIC(lm(meanAr ~ betw_dir+spec, DATA))
AIC(lm(meanAr ~ betw_dir, DATA))
model <- lm.bind(meanAr ~ betw_dir, DATA, "Ar_betw_dir", critval, step=T)
DATA <- model[[1]]
lm_Ar_betw_dir <- model[[2]]
slm_Ar_betw_dir <- model[[3]]

##*** Models meanHo ####
shapiro.test(DATA$meanHo)

# Interaction model
int.mod.Ho <- lm(meanHo ~ updist*betw_dir*std_clos_undir*catch*spec, DATA, na.action = "na.fail")
summary(int.mod.Ho)
opar <- par(mfrow=c(2,2))
plot(int.mod.Ho)
par(opar)
step(int.mod.Ho)
# best model keeps some interactions
if(log_trans){
  int.sel.mod.Ho <- lm(meanHo ~ updist * betw_dir * std_clos_undir * catch * 
                         spec, data = DATA, na.action = "na.fail")
}else{
  int.sel.mod.Ho <- lm(meanHo ~ updist + betw_dir + std_clos_undir + catch + 
                     spec + updist:betw_dir + updist:std_clos_undir + betw_dir:std_clos_undir + 
                     updist:catch + betw_dir:catch + std_clos_undir:catch + updist:spec + 
                     betw_dir:spec + std_clos_undir:spec + catch:spec + updist:betw_dir:std_clos_undir + 
                     betw_dir:std_clos_undir:catch + updist:betw_dir:spec + std_clos_undir:catch:spec, 
                   data = DATA, na.action = "na.fail")
}

# Linear model
lin.mod.Ho <- lm(meanHo ~ updist+betw_dir+std_clos_undir+catch+spec, DATA, na.action = "na.fail")
summary(lin.mod.Ho)
opar <- par(mfrow=c(2,2))
plot(lin.mod.Ho)
par(opar)
MuMIn::dredge(lin.mod.Ho)
step(lin.mod.Ho)
# best model: btw_dir+spc+std_cls_und (untransformed); btw_dir+spc+std_cls_und (log-transformed)
if(log_trans){
  sel.mod.Ho <- lm(meanHo ~ betw_dir+std_clos_undir, DATA, na.action = "na.fail")
}else{
  sel.mod.Ho <- lm(meanHo ~ betw_dir+std_clos_undir+spec, DATA, na.action = "na.fail")
}  

# Comparison
AIC(int.sel.mod.Ho)
AIC(sel.mod.Ho)
# AIC(glm(meanHo ~ betw_dir+std_clos_undir+spec, DATA, family= "quasibinomial"))
# linear model outperforms interaction model
summary(int.sel.mod.Ho)
opar <- par(mfrow=c(2,2))
plot(int.sel.mod.Ho)
par(opar)
car::vif(int.sel.mod.Ho)
summary(sel.mod.Ho)
opar <- par(mfrow=c(2,2))
plot(sel.mod.Ho)
par(opar)
car::vif(sel.mod.Ho)

#### Full LM of mean observed heterozygosity without interactions
model <- lm.bind(meanHo ~ updist+betw_dir+std_clos_undir+catch+spec, DATA, "Ho_full", critval, step=T)
DATA <- model[[1]]
lm_Ho_full <- model[[2]]
summary(lm_Ho_full)
car::vif(lm_Ho_full)
slm_Ho_full <- model[[3]]
summary(slm_Ho_full)
car::vif(slm_Ho_full) # should be the same as above (car::vif(sel.mod.Ho))

#### LM of observed heterozygosity by upstream distance * species
AIC(lm(meanHo ~ updist+spec, DATA))
AIC(lm(meanHo ~ updist, DATA))
model <- lm.bind(meanHo ~ updist, DATA, "Ho_updist", critval, step=F)
DATA <- model[[1]]
lm_Ho_updist <- model[[2]]
# slm_Ho_updist <- model[[3]]

#### LM of observed heterozygosity by closeness centrality * species
AIC(lm(meanHo ~ std_clos_undir+spec, DATA))
AIC(lm(meanHo ~ std_clos_undir, DATA))
model <- lm.bind(meanHo ~ std_clos_undir, DATA, "Ho_clos", critval, step=T)
DATA <- model[[1]]
lm_Ho_clos <- model[[2]]
slm_Ho_clos <- model[[3]]

#### LM of observed heterozygosity by directed betweenness centrality * species
AIC(lm(meanHo ~ betw_dir+spec, DATA))
AIC(lm(meanHo ~ betw_dir, DATA))
model <- lm.bind(meanHo ~ betw_dir, DATA, "Ho_betw_dir", critval, step=T)
DATA <- model[[1]]
lm_Ho_betw_dir <- model[[2]]
slm_Ho_betw_dir <- model[[3]]

##*** Models He ####
shapiro.test(DATA$meanHs)

# Interaction model
int.mod.Hs <- lm(meanHs ~ updist*betw_dir*std_clos_undir*catch*spec, DATA, na.action = "na.fail")
summary(int.mod.Hs)
opar <- par(mfrow=c(2,2))
plot(int.mod.Hs)
par(opar)
step(int.mod.Hs)
# best model keeps some interactions
if(log_trans){
  int.sel.mod.Hs <- lm(meanHs ~ betw_dir + spec + betw_dir:spec, data = DATA, na.action = "na.fail")
}else{
  int.sel.mod.Hs <- lm(meanHs ~ catch, data = DATA, na.action = "na.fail")
}

# Linear model
lin.mod.Hs <- lm(meanHs ~ updist+betw_dir+std_clos_undir+catch+spec, DATA, na.action = "na.fail")
summary(lin.mod.Hs)
opar <- par(mfrow=c(2,2))
plot(lin.mod.Hs)
par(opar)
MuMIn::dredge(lin.mod.Hs)
step(lin.mod.Hs)
# best model: ctc (untransformed); (Int) (log-transformed); ctc+sps (log-transformed, step)
if(log_trans){
  sel.mod.Hs <- lm(meanHs ~ 1, DATA, na.action = "na.fail")
}else{
  sel.mod.Hs <- lm(meanHs ~ catch, DATA, na.action = "na.fail")
}

# Comparison
AIC(int.sel.mod.Hs)
AIC(sel.mod.Hs)
# AIC(glm(meanHs ~ catch, DATA, family= "quasibinomial"))
# linear model outperforms interaction model
summary(sel.mod.Hs)
opar <- par(mfrow=c(2,2))
plot(sel.mod.Hs)
par(opar)
car::vif(sel.mod.Hs)

#### Full LM of expected heterozygosity without interactions
model <- lm.bind(meanHs ~ updist+betw_dir+std_clos_undir+catch+spec, DATA, "Hs_full", critval, step=T)
DATA <- model[[1]]
lm_Hs_full <- model[[2]]
summary(lm_Hs_full)
car::vif(lm_Hs_full)
slm_Hs_full <- model[[3]]
summary(slm_Hs_full)
car::vif(slm_Hs_full) # should be the same as above (car::vif(sel.mod.Hs))

#### LM of expected heterozygosity by upstream distance * species
AIC(lm(meanHs ~ updist+spec, DATA))
AIC(lm(meanHs ~ updist, DATA))
model <- lm.bind(meanHo ~ updist, DATA, "Hs_updist", critval, step=F)
DATA <- model[[1]]
lm_Hs_updist <- model[[2]]
# slm_Hs_updist <- model[[3]]

#### LM of expected heterozygosity by closeness centrality * species
AIC(lm(meanHs ~ std_clos_undir+spec, DATA))
AIC(lm(meanHs ~ std_clos_undir, DATA))
model <- lm.bind(meanHs ~ std_clos_undir, DATA, "Hs_clos", critval, step=F)
DATA <- model[[1]]
lm_Hs_clos <- model[[2]]
# slm_Hs_clos <- model[[3]]

#### LM of expected heterozygosity by directed betweenness centrality * species
AIC(lm(meanHs ~ betw_dir+spec, DATA))
AIC(lm(meanHs ~ betw_dir, DATA))
model <- lm.bind(meanHs ~ betw_dir, DATA, "Hs_betw_dir", critval, step=T)
DATA <- model[[1]]
lm_Hs_betw_dir <- model[[2]]
slm_Hs_betw_dir <- model[[3]]

##*** Models Fst ####
#### Spatial distance between populations with Fos A and between populations with Fos B
DIST_B <- distances(net, v=V(net)[match(microsite_B,V(net)$name)], to=V(net)[match(microsite_B,V(net)$name)], weights=E(net))
DIST_A_red <- distances(net, v=V(net)[match(microsite_A_red,V(net)$name)], to=V(net)[match(microsite_A_red,V(net)$name)], weights=E(net))

#### Mantel test of genetic differentiation by instream distance
mantel_A_red <- vegan::mantel(FST_A_red, DIST_A_red, method="pearson", permutations=1000)
mantel_B <- vegan::mantel(FST_B, DIST_B, method="pearson", permutations=1000)

#### LM of genetic differentiation by spatial distance
model_dist_fst_B <- lm(fst_B~dist_B)
summary(model_dist_fst_B)
opar <- par(mfrow=c(2,2))
plot(model_dist_fst_B)
par(opar)
resp_dist_fst_B <- predict(model_dist_fst_B, list(dist_B = range_dist_B), type="response")
model_dist_fst_A_red <- lm(fst_A_red~dist_A_red)
summary(model_dist_fst_A_red)
opar <- par(mfrow=c(2,2))
plot(model_dist_fst_A_red)
par(opar)
resp_dist_fst_A_red <- predict(model_dist_fst_A_red, list(dist_A_red = range_dist_A_red), type="response")

shapiro.test(DISTDATA$nonneg_fst)

int.mod.Fst <- lm(nonneg_fst ~ dist*spec, DISTDATA, na.action = "na.fail")
summary(int.mod.Fst)
opar <- par(mfrow=c(2,2))
plot(int.mod.Fst)
par(opar)

MuMIn::dredge(int.mod.Fst)
step(int.mod.Fst)
# interaction model outperforms linear model
car::vif(int.mod.Fst)

lin.mod.Fst <- lm(nonneg_fst ~ dist+spec, DISTDATA, na.action = "na.fail")
summary(lin.mod.Fst)
opar <- par(mfrow=c(2,2))
plot(lin.mod.Fst)
par(opar)
MuMIn::dredge(lin.mod.Fst)
car::vif(lin.mod.Fst)

log.mod.Fst <- lm(nonneg_fst ~ log(dist)*spec, DISTDATA, na.action = "na.fail")

power <- seq(0,1,0.01)
AICpower <- c()
for (i in 1:length(power)){
  pow.mod.Fst <- lm(nonneg_fst ~ I(dist^power[i])*spec, DISTDATA, na.action = "na.fail")
  AICpower[i] <- AIC(pow.mod.Fst)
}
plot(AICpower~power, xlab="power term", ylab="AIC")
pow.mod.Fst <- lm(nonneg_fst ~ I(dist^power[which.min(AICpower)])*spec, DISTDATA, na.action = "na.fail")

AIC(int.mod.Fst)
AIC(lin.mod.Fst)
AIC(log.mod.Fst)
AIC(pow.mod.Fst)
summary(pow.mod.Fst)

#### Full LM of expected heterozygosity without interactions
model <- lm.bind(nonneg_fst ~ I(dist^power[which.min(AICpower)])*spec, DISTDATA, "fst_power", critval, step=T)
DISTDATA <- model[[1]]
lm_fst_power <- model[[2]]
summary(lm_fst_power)
car::vif(lm_fst_power)
slm_fst_power <- model[[3]]
summary(slm_fst_power)
car::vif(slm_fst_power)

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

orthdist_Hs_A <- vector()
orthdist_Hs_A_directed <- vector()
sum_orthdist_Hs_A <- vector()
hist_Hs_A <- matrix(nrow=length(modsite_GfosA), ncol=ncol(Hs_modelled)-1)
hist_Hs_A_directed <- matrix(nrow=length(modsite_GfosA), ncol=ncol(Hs_modelled)-1)

orthdist_Hs_B <- vector()
orthdist_Hs_B_directed <- vector()
sum_orthdist_Hs_B <- vector()
hist_Hs_B <- matrix(nrow=length(modsite_GfosB), ncol=ncol(Hs_modelled)-1)
hist_Hs_B_directed <- matrix(nrow=length(modsite_GfosB), ncol=ncol(Hs_modelled)-1)

##*** Calculating distances####
for (i in 2:ncol(Ar_modelled)){
  Ar_Mod_A <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosA]
  Ar_Mod_B <- Ar_modelled[,i][rownames(Ar_modelled)%in%modsite_GfosB]
  Ho_Mod_A <- Ho_modelled[,i][rownames(Ho_modelled)%in%modsite_GfosA]
  Ho_Mod_B <- Ho_modelled[,i][rownames(Ho_modelled)%in%modsite_GfosB]
  Hs_Mod_A <- Hs_modelled[,i][rownames(Hs_modelled)%in%modsite_GfosA]
  Hs_Mod_B <- Hs_modelled[,i][rownames(Hs_modelled)%in%modsite_GfosB]
  
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
  for (j in 1:length(meanHs_A_red_updist)){
    point <- cbind(meanHs_A_red_updist,Hs_Mod_A)[j,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    orthdist_Hs_A[j] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
    orthdist_Hs_A_directed[j] <- sign(seg[2]-seg[4])*orthdist_Ho_A[j]
  }
  for (k in 1:length(meanHs_B_updist)){
    point <- cbind(meanHs_B_updist,Hs_Mod_B)[k,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    orthdist_Hs_B[k] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
    orthdist_Hs_B_directed[k] <- sign(seg[2]-seg[4])*orthdist_Ho_B[k]
  }
  sum_orthdist_Ar_A[[i-1]] <- sum(orthdist_Ar_A)
  sum_orthdist_Ar_B[[i-1]] <- sum(orthdist_Ar_B)
  sum_orthdist_Ho_A[[i-1]] <- sum(orthdist_Ho_A)
  sum_orthdist_Ho_B[[i-1]] <- sum(orthdist_Ho_B)
  sum_orthdist_Hs_A[[i-1]] <- sum(orthdist_Hs_A)
  sum_orthdist_Hs_B[[i-1]] <- sum(orthdist_Hs_B)
  
  hist_Ar_A[,i-1] <- orthdist_Ar_A
  hist_Ar_B[,i-1] <- orthdist_Ar_B
  hist_Ho_A[,i-1] <- orthdist_Ho_A
  hist_Ho_B[,i-1] <- orthdist_Ho_B
  hist_Hs_A[,i-1] <- orthdist_Hs_A
  hist_Hs_B[,i-1] <- orthdist_Hs_B
  
  hist_Ar_A_directed[,i-1] <- orthdist_Ar_A_directed
  hist_Ar_B_directed[,i-1] <- orthdist_Ar_B_directed
  hist_Ho_A_directed[,i-1] <- orthdist_Ho_A_directed
  hist_Ho_B_directed[,i-1] <- orthdist_Ho_B_directed
  hist_Hs_A_directed[,i-1] <- orthdist_Hs_A_directed
  hist_Hs_B_directed[,i-1] <- orthdist_Hs_B_directed
}
names(sum_orthdist_Ar_A) <- colnames(Ar_modelled)[-1]
names(sum_orthdist_Ar_B) <- colnames(Ar_modelled)[-1]
names(sum_orthdist_Ho_A) <- colnames(Ho_modelled)[-1]
names(sum_orthdist_Ho_B) <- colnames(Ho_modelled)[-1]
names(sum_orthdist_Hs_A) <- colnames(Hs_modelled)[-1]
names(sum_orthdist_Hs_B) <- colnames(Hs_modelled)[-1]

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

min(meanHs_A_red)
max(meanHs_A_red)
mean(meanHs_A_red)
median(meanHs_A_red)
sd(meanHs_A_red)
min(meanHs_B)
max(meanHs_B)
mean(meanHs_B)
median(meanHs_B)
sd(meanHs_B)

min(fst_A_red)
max(fst_A_red)
mean(fst_A_red)
median(fst_A_red)
sd(fst_A_red)
min(fst_B)
max(fst_B)
mean(fst_B)
median(fst_B)
sd(fst_B)

mantel_A_red
mantel_B

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
  pdf(paste0(WD,"/Analysis_",output,"/Fig1.pdf"), width=1.5*fig.width, height=(3*fig.width)/4.5)
}else{
  png(paste0(WD,"/Analysis_",output,"/Fig1.png"), width=1.5*fig.width, height=(3*fig.width)/4.5, units="in", res=300)
}
nf <- layout(matrix(c(12, 7, 8, 9,
                      10, 1, 2, 3,
                      11, 4, 5, 6), nrow=3, byrow=T),
             widths=c(0.5,3,3,3),
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

#### Map of expected heterozygosity in G. fossarum A
par(mar=figmar, mgp=figmgp, tcl=0.2, xaxs="i", yaxs="i")
if (internal){
  river_plot(north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, lines_rhone=F, plot_ticino=F, lines_ticino=F, plot_inn=F, lines_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}else{
  river_plot(north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, plot_ticino=F, plot_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}
mtext(side = 2, text = expression(paste(bold("CH1903"), " / ", bold("LV03"))), line = 1 + 1, cex = 1)
values <- (meanHs_A_red-min(meanHs_A_red))/ (max(meanHs_A_red)-min(meanHs_A_red)) # transform meanHs values to range [0,1] for heatmap plotting
for(i in 1:length(match_A_red)){
  points(site_coord$x[match_A_red[i]], site_coord$y[match_A_red[i]], col=rgb2hex(col_fun(values))[i], pch=19, cex=2.5)
  l1 <- mean(meanHs_A_red)+(col_switch)*(max(meanHs_A_red)-median(meanHs_A_red))
  l2 <- mean(meanHs_A_red)-(col_switch)*(median(meanHs_A_red)-min(meanHs_A_red))
  text(site_coord$x[match_A_red[i]], site_coord$y[match_A_red[i]], round(meanHs_A_red[i],1), col=ifelse(meanHs_A_red[i]>l1|meanHs_A_red[i]<l2,"white","black"), cex=0.8)
}
gradient.legend(meanHs_A_red, val.cex=1, palette=col_pal)
mtext(side = 3, text = lab_sub[5], line = 0.5, adj=0, cex = 1.5)

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

#### Map of expected heterozygosity in G. fossarum B
par(mar=figmar, mgp=figmgp, tcl=0.2, xaxs="i", yaxs="i")
if (internal){
  river_plot(north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, lines_rhone=F, plot_ticino=F, lines_ticino=F, plot_inn=F, lines_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}else{
  river_plot(north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, col_water=col_water, plot_rhone=F, plot_ticino=F, plot_inn=F, lakes=TRUE, rivers=TRUE, axes="no_label", river_nr=FALSE)
}
mtext(side = 1, text = expression(paste(bold("CH1903"), " / ", bold("LV03"))), line = 1 + 1, cex = 1)
mtext(side = 2, text = expression(paste(bold("CH1903"), " / ", bold("LV03"))), line = 1 + 1, cex = 1)
values <- (meanHs_B-min(meanHs_B))/ (max(meanHs_B)-min(meanHs_B)) # transform meanHs values to range [0,1] for heatmap plotting
for(i in 1:length(match_B)){
  points(site_coord$x[match_B[i]], site_coord$y[match_B[i]], col=rgb2hex(col_fun(values))[i], pch=19, cex=2.5)
  l1 <- mean(meanHs_B)+(col_switch)*(max(meanHs_B)-median(meanHs_B))
  l2 <- mean(meanHs_B)-(col_switch)*(median(meanHs_B)-min(meanHs_B))
  text(site_coord$x[match_B[i]], site_coord$y[match_B[i]], round(meanHs_B[i],1), col=ifelse(meanHs_B[i]>l1|meanHs_B[i]<l2,"white","black"), cex=0.8)
}
gradient.legend(meanHs_B, val.cex = 1, palette=col_pal)
mtext(side = 3, text = lab_sub[6], line = 0.5, adj=0, cex = 1.5)

par(mar=c(0,0,0,0),mgp=c(3,1,0))
plot.new()
text(0.5,0.5, lab_Ar ,adj=c(0.5,0.5), cex=2.5)
plot.new()
text(0.5,0.5, lab_Ho ,adj=c(0.5,0.5), cex=2.5)
plot.new()
text(0.5,0.5, lab_Hs ,adj=c(0.5,0.5), cex=2.5)
plot.new()
text(1,0.5,label_A,adj=c(0.5,0), cex=2.5, srt=90)
plot.new()
text(1,0.5,label_B,adj=c(0.5,0), cex=2.5, srt=90)
dev.off()

##** FIG 2 ####
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/Fig2.pdf"), width=22.5, height=10)
}else{
  png(paste0(WD,"/Analysis_",output,"/Fig2.png"), width=22.5, height=10, units="in", res=300)
}
par(mfrow=c(2,3))
ylim <- c(min(meanHo_A_red,meanHo_B,meanHs_A_red,meanHs_B),max(meanHo_A_red,meanHo_B,meanHs_A_red,meanHs_B))
#### meanAr~updist
GLMplot("meanAr","updist",DATA,model="slm_Ar_updist",pt.cex=2, CI_border=F,xlabel="Upstream distance [km]",ylabel=lab_Ar, xrev=T, xax="n", cex.lab=2, cex.axis=1.5, cex.legend=1.5)
axis(1, c(0,50,100,150,200,250,300), at=c(0,50000,100000,150000,200000,250000,300000), cex.axis=1.5)
text(300000, par("usr")[4], lab_sub[1], cex=2, adj=c(0.5,1))
#### meanHo~updist
GLMplot("meanHo","updist",DATA,model="lm_Ho_updist",ylim=ylim,pt.cex=2, CI_border=F,xlabel="Upstream distance [km]",ylabel=lab_Ho, xrev=T, xax="n", cex.lab=2, cex.axis=1.5, legend=F)
axis(1, c(0,50,100,150,200,250,300), at=c(0,50000,100000,150000,200000,250000,300000), cex.axis=1.5)
text(300000, par("usr")[4], lab_sub[3], cex=2, adj=c(0.5,1))
#### meanHs~updist
GLMplot("meanHs","updist",DATA,model="lm_Hs_updist",ylim=ylim,pt.cex=2, CI_border=F,xlabel="Upstream distance [km]",ylabel=lab_Hs, xrev=T, xax="n", cex.lab=2, cex.axis=1.5, legend=F)
axis(1, c(0,50,100,150,200,250,300), at=c(0,50000,100000,150000,200000,250000,300000), cex.axis=1.5)
text(300000, par("usr")[4], lab_sub[5], cex=2, adj=c(0.5,1))
#### meanAr~std_clos_undir
GLMplot("meanAr","std_clos_undir",DATA,model="slm_Ar_clos",pt.cex=2, CI_border=F,xlabel="Standardized closeness centrality",ylabel=lab_Ar, legend=F)
text(0, par("usr")[4], lab_sub[2], cex=2, adj=c(0.5,1))
#### meanHo~std_clos_undir
GLMplot("meanHo","std_clos_undir",DATA,model="slm_Ho_clos",ylim=ylim,pt.cex=2, CI_border=F,xlabel="Standardized closeness centrality",ylabel=lab_Ho, legend=F)
text(0, par("usr")[4], lab_sub[4], cex=2, adj=c(0.5,1))
# dev.off()
#### meanHs~std_clos_undir
GLMplot("meanHs","std_clos_undir",DATA,model="lm_Hs_clos",ylim=ylim,pt.cex=2, CI_border=F,xlabel="Standardized closeness centrality",ylabel=lab_Hs, legend=F)
text(0, par("usr")[4], lab_sub[6], cex=2, adj=c(0.5,1))
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
  ylim <- 40
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
  text(0.5,0.7,paste0(lab_Ar,": Distribution of ",measure3a),adj=c(0.5,0.5),cex=3)
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
text(-2.5,15,"See Fig. S10 for all plots",cex=1, adj=0, col="darkgrey")

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
##** Fst by instream distance (power function) by species, showing IBD
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/Fig5.pdf"), width=8, height=6)
}else{
  png(paste0(WD,"/Analysis_",output,"/Fig5.png"), width=8, height=6, units="in", res=300)
}
GLMplot("fst","dist",dat=DISTDATA,model="slm_fst_power", CI_border = F, xlabel="Instream distance [km]",ylabel=expression('Genetic diff. [Pairwise Nei F'[ST]*']'),xax="n",pointtrans = T)
axis(1, c(0,50,100,150,200,250), at=c(0,50000,100000,150000,200000,250000), cex.axis=1.5)
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
##** Correlation plot all explanatory variables
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS2.pdf"), width=8, height=6)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS2.png"), width=8, height=6, units="in", res=300)
}
PerformanceAnalytics::chart.Correlation(expl.var1, method = "kendall")
dev.off()

##** FIG S3 ####
##** Correlation plot selected and transformed explanatory variables
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS3.pdf"), width=8, height=6)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS3.png"), width=8, height=6, units="in", res=300)
}
PerformanceAnalytics::chart.Correlation(expl.var, method = "kendall")
dev.off()

##** FIG S4 ####
#Perpendicular offset example
x <- seq(1,10,1)
y1 <- seq(1,10,1)
y2 <- seq(3,12,1)
y3 <- seq(-1,8,1)
y4 <- seq(2,20,2)
y5 <- seq(20,2,-2)
y6 <- sample(seq(5.499,5.501,0.00001),10)
y6 <- y6[order(y6)[c(1,10,3,8,5,6,7,4,9,2)]]

cor1 <- cor(x,y1)
cor2 <- cor(x,y2)
cor3 <- cor(x,y3)
cor4 <- cor(x,y4)
cor5 <- cor(x,y5)
cor6 <- abs(cor(x,y6))

orthdist_y1 <- c()
orthdist_y1_directed <- c()
orthdist_y2 <- c()
orthdist_y2_directed <- c()
orthdist_y3 <- c()
orthdist_y3_directed <- c()
orthdist_y4 <- c()
orthdist_y4_directed <- c()
orthdist_y5 <- c()
orthdist_y5_directed <- c()
orthdist_y6 <- c()
orthdist_y6_directed <- c()

# all combined in one plot
col1 <- "darkgreen"
col2 <- "steelblue"
col3 <- "red"
col4 <- "orange"
col5 <- "purple"
col6 <- "gold"
plot(x,y1, xlim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)), ylim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)),
     asp=1, type="n", xlab="Empirical data", ylab="Simulated data")
abline(0,1, lwd=1, lty=2)

points(x,y1,col=col1, pch=16)
for (i in 1:length(x)){
  point <- cbind(x,y1)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y1[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y1_directed[i] <- sign(seg[2]-seg[4])*orthdist_y1[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col1, lty=2, lwd=0.5)
}

points(x,y2,col=col2, pch=16)
for (i in 1:length(x)){
  point <- cbind(x,y2)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y2[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y2_directed[i] <- sign(seg[2]-seg[4])*orthdist_y2[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col2, lty=2, lwd=0.5)
}

points(x,y3,col=col3, pch=16)
for (i in 1:length(x)){
  point <- cbind(x,y3)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y3[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y3_directed[i] <- sign(seg[2]-seg[4])*orthdist_y3[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col3, lty=2, lwd=0.5)
}

points(x,y4,col=col4, pch=16)
for (i in 1:length(x)){
  point <- cbind(x,y4)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y4[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y4_directed[i] <- sign(seg[2]-seg[4])*orthdist_y4[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col4, lty=2, lwd=0.5)
}

points(x,y5,col=col5, pch=16)
for (i in 1:length(x)){
  point <- cbind(x,y5)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y5[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y5_directed[i] <- sign(seg[2]-seg[4])*orthdist_y5[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col5, lty=2, lwd=0.5)
}

points(x,y6,col=col6, pch=16)
for (i in 1:length(x)){
  point <- cbind(x,y6)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y6[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y6_directed[i] <- sign(seg[2]-seg[4])*orthdist_y6[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col6, lty=2, lwd=0.5)
}

text(8,0,paste0("Cor: ", specify_decimal(cor1,1)), col=col1)
text(15,0,paste0("SPO: ", specify_decimal(sum(orthdist_y1),1)), col=col1)
text(20,0,paste0("MPO: ", specify_decimal(median(orthdist_y1),1)),col=col1)
text(25,0,paste0("DMPO: ", specify_decimal(median(orthdist_y1_directed),1)), col=col1)
text(8,1,paste0("Cor: ", specify_decimal(cor2,1)), col=col2)
text(15,1,paste0("SPO: ", specify_decimal(sum(orthdist_y2),1)), col=col2)
text(20,1,paste0("MPO: ", specify_decimal(median(orthdist_y2),1)),col=col2)
text(25,1,paste0("DMPO: ", specify_decimal(median(orthdist_y2_directed),1)), col=col2)
text(8,2,paste0("Cor: ", specify_decimal(cor3,1)), col=col3)
text(15,2,paste0("SPO: ", specify_decimal(sum(orthdist_y3),1)), col=col3)
text(20,2,paste0("MPO: ", specify_decimal(median(orthdist_y3),1)),col=col3)
text(25,2,paste0("DMPO: ", specify_decimal(median(orthdist_y3_directed),1)), col=col3)
text(8,3,paste0("Cor: ", specify_decimal(cor4,1)), col=col4)
text(15,3,paste0("SPO: ", specify_decimal(sum(orthdist_y4),1)), col=col4)
text(20,3,paste0("MPO: ", specify_decimal(median(orthdist_y4),1)),col=col4)
text(25,3,paste0("DMPO: ", specify_decimal(median(orthdist_y4_directed),1)), col=col4)
text(8,4,paste0("Cor: ", specify_decimal(cor5,1)), col=col5)
text(15,4,paste0("SPO: ", specify_decimal(sum(orthdist_y5),1)), col=col5)
text(20,4,paste0("MPO: ", specify_decimal(median(orthdist_y5),1)),col=col5)
text(25,4,paste0("DMPO: ", specify_decimal(median(orthdist_y5_directed),1)), col=col5)
text(8,5,paste0("Cor: ", specify_decimal(cor6,1)), col=col6)
text(15,5,paste0("SPO: ", specify_decimal(sum(orthdist_y6),1)), col=col6)
text(20,5,paste0("MPO: ", specify_decimal(median(orthdist_y6),1)),col=col6)
text(25,5,paste0("DMPO: ", specify_decimal(median(orthdist_y6_directed),1)), col=col6)

# six plots
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS4.pdf"), width=9, height=6)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS4.png"), width=9, height=6, units="in", res=300)
}
nf <- multipanel.layout(main.col=1,main.row=1,sub.col=3,sub.row=2,sub1=F,sub2=F,main=pla$main, w.legend=0.05, w.axis=0.3, h.axis=0.3)

x <- seq(1,10,1)
y1 <- seq(1,10,1)
y2 <- seq(3,12,1)
y3 <- seq(-1,8,1)
y4 <- seq(2,20,2)
y5 <- seq(20,2,-2)

while(cor6>=0.05){
  y6 <- sample(seq(5.499,5.501,0.00001),10)
  y6 <- y6[order(y6)[c(1,10,3,8,5,6,7,4,9,2)]]
  cor6 <- abs(cor(x,y6))
}

par(mar=c(0,0,0,0))
x_text <- 20
y_text_1 <- 6
y_text_2 <- 4
y_text_3 <- 2
y_text_4 <- 0
cex_perp <- 2
col_perp <- "steelblue"
col1 <- col_perp
col2 <- col_perp
col3 <- col_perp
col4 <- col_perp
col5 <- col_perp
col6 <- col_perp
lwd_perp <- 1
pt.cex <- 2

plot(x,y1, xlim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)), ylim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)),
     asp=1, type="n", xlab="Empirical data", ylab="Simulated data", xaxt="n", yaxt="n")
abline(0,1, lwd=1, lty=2)
points(x,y1,col=col1, pch=16, cex=pt.cex)
for (i in 1:length(x)){
  point <- cbind(x,y1)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y1[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y1_directed[i] <- sign(seg[2]-seg[4])*orthdist_y1[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col1, lty=2, lwd=lwd_perp)
}
text(par("usr")[1]+((par("usr")[2]-par("usr")[1])/40),par("usr")[4]-((par("usr")[4]-par("usr")[3])/40), "(a)", adj=c(0,1), cex=cex_perp)
text(x_text,y_text_1,paste0("Cor: ", specify_decimal(cor1,1)), col=col1, adj=1, cex=cex_perp)
text(x_text,y_text_2,paste0("SPO: ", specify_decimal(sum(orthdist_y1),1)), col=col1, adj=1, cex=cex_perp)
text(x_text,y_text_3,paste0("MPO: ", specify_decimal(median(orthdist_y1),1)),col=col1, adj=1, cex=cex_perp)
text(x_text,y_text_4,paste0("DMPO: ", specify_decimal(median(orthdist_y1_directed),1)), col=col1, adj=1, cex=cex_perp)
axis(2)

plot(x,y4, xlim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)), ylim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)),
     asp=1, type="n", xlab="Empirical data", ylab="Simulated data", xaxt="n", yaxt="n")
abline(0,1, lwd=1, lty=2)
points(x,y4,col=col4, pch=16, cex=pt.cex)
for (i in 1:length(x)){
  point <- cbind(x,y4)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y4[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y4_directed[i] <- sign(seg[2]-seg[4])*orthdist_y4[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col4, lty=2, lwd=lwd_perp)
}
text(par("usr")[1]+((par("usr")[2]-par("usr")[1])/40),par("usr")[4]-((par("usr")[4]-par("usr")[3])/40), "(d)", adj=c(0,1), cex=cex_perp)
text(x_text,y_text_1,paste0("Cor: ", specify_decimal(cor4,1)), col=col4, adj=1, cex=cex_perp)
text(x_text,y_text_2,paste0("SPO: ", specify_decimal(sum(orthdist_y4),1)), col=col4, adj=1, cex=cex_perp)
text(x_text,y_text_3,paste0("MPO: ", specify_decimal(median(orthdist_y4),1)),col=col4, adj=1, cex=cex_perp)
text(x_text,y_text_4,paste0("DMPO: ", specify_decimal(median(orthdist_y4_directed),1)), col=col4, adj=1, cex=cex_perp)
axis(1)
axis(2)

plot(x,y2, xlim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)), ylim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)),
     asp=1, type="n", xlab="Empirical data", ylab="Simulated data", xaxt="n", yaxt="n")
abline(0,1, lwd=1, lty=2)
points(x,y2,col=col2, pch=16, cex=pt.cex)
for (i in 1:length(x)){
  point <- cbind(x,y2)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y2[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y2_directed[i] <- sign(seg[2]-seg[4])*orthdist_y2[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col2, lty=2, lwd=lwd_perp)
}
text(par("usr")[1]+((par("usr")[2]-par("usr")[1])/40),par("usr")[4]-((par("usr")[4]-par("usr")[3])/40), "(b)", adj=c(0,1), cex=cex_perp)
text(x_text,y_text_1,paste0("Cor: ", specify_decimal(cor2,1)), col=col2, adj=1, cex=cex_perp)
text(x_text,y_text_2,paste0("SPO: ", specify_decimal(sum(orthdist_y2),1)), col=col2, adj=1, cex=cex_perp)
text(x_text,y_text_3,paste0("MPO: ", specify_decimal(median(orthdist_y2),1)),col=col2, adj=1, cex=cex_perp)
text(x_text,y_text_4,paste0("DMPO: ", specify_decimal(median(orthdist_y2_directed),1)), col=col2, adj=1, cex=cex_perp)

plot(x,y5, xlim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)), ylim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)),
     asp=1, type="n", xlab="Empirical data", ylab="Simulated data", xaxt="n", yaxt="n")
abline(0,1, lwd=1, lty=2)
points(x,y5,col=col5, pch=16, cex=pt.cex)
for (i in 1:length(x)){
  point <- cbind(x,y5)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y5[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y5_directed[i] <- sign(seg[2]-seg[4])*orthdist_y5[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col5, lty=2, lwd=lwd_perp)
}
text(par("usr")[1]+((par("usr")[2]-par("usr")[1])/40),par("usr")[4]-((par("usr")[4]-par("usr")[3])/40), "(e)", adj=c(0,1), cex=cex_perp)
text(x_text,y_text_1,paste0("Cor: ", specify_decimal(cor5,1)), col=col5, adj=1, cex=cex_perp)
text(x_text,y_text_2,paste0("SPO: ", specify_decimal(sum(orthdist_y5),1)), col=col5, adj=1, cex=cex_perp)
text(x_text,y_text_3,paste0("MPO: ", specify_decimal(median(orthdist_y5),1)),col=col5, adj=1, cex=cex_perp)
text(x_text,y_text_4,paste0("DMPO: ", specify_decimal(median(orthdist_y5_directed),1)), col=col5, adj=1, cex=cex_perp)
axis(1)

plot(x,y3, xlim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)), ylim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)),
     asp=1, type="n", xlab="Empirical data", ylab="Simulated data", xaxt="n", yaxt="n")
abline(0,1, lwd=1, lty=2)
points(x,y3,col=col3, pch=16, cex=pt.cex)
for (i in 1:length(x)){
  point <- cbind(x,y3)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y3[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y3_directed[i] <- sign(seg[2]-seg[4])*orthdist_y3[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col3, lty=2, lwd=lwd_perp)
}
text(par("usr")[1]+((par("usr")[2]-par("usr")[1])/40),par("usr")[4]-((par("usr")[4]-par("usr")[3])/40), "(c)", adj=c(0,1), cex=cex_perp)
text(x_text,y_text_1,paste0("Cor: ", specify_decimal(cor3,1)), col=col3, adj=1, cex=cex_perp)
text(x_text,y_text_2,paste0("SPO: ", specify_decimal(sum(orthdist_y3),1)), col=col3, adj=1, cex=cex_perp)
text(x_text,y_text_3,paste0("MPO: ", specify_decimal(median(orthdist_y3),1)),col=col3, adj=1, cex=cex_perp)
text(x_text,y_text_4,paste0("DMPO: ", specify_decimal(median(orthdist_y3_directed),1)), col=col3, adj=1, cex=cex_perp)

plot(x,y6, xlim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)), ylim=c(min(x,y1,y2,y3,y4,y5,y6),max(x,y1,y2,y3,y4,y5,y6)),
     asp=1, type="n", xlab="Empirical data", ylab="Simulated data", xaxt="n", yaxt="n")
abline(0,1, lwd=1, lty=2)
points(x,y6,col=col6, pch=16, cex=pt.cex)
for (i in 1:length(x)){
  point <- cbind(x,y6)[i,]
  seg <- unlist(perp.segment.coord(point[1],point[2]))
  orthdist_y6[i] <- euc.dist(c(seg[1],seg[2]),c(seg[3],seg[4]))
  orthdist_y6_directed[i] <- sign(seg[2]-seg[4])*orthdist_y6[i]
  segments(seg[1],seg[2],seg[3],seg[4], col=col6, lty=2, lwd=lwd_perp)
}
text(par("usr")[1]+((par("usr")[2]-par("usr")[1])/40),par("usr")[4]-((par("usr")[4]-par("usr")[3])/40), "(f)", adj=c(0,1), cex=cex_perp)
text(x_text,y_text_1,paste0("Cor: ", specify_decimal(cor6,1)), col=col6, adj=1, cex=cex_perp)
text(x_text,y_text_2,paste0("SPO: ", specify_decimal(sum(orthdist_y6),1)), col=col6, adj=1, cex=cex_perp)
text(x_text,y_text_3,paste0("MPO: ", specify_decimal(median(orthdist_y6),1)),col=col6, adj=1, cex=cex_perp)
text(x_text,y_text_4,paste0("DMPO: ", specify_decimal(median(orthdist_y6_directed),1)), col=col6, adj=1, cex=cex_perp)
axis(1)

plot.new()
text(0.5,0.5,"Simulated data", cex=2, srt = 90)

plot.new()
text(0.5,0.5,"Simulated data", cex=2, srt = 90)

plot.new()
text(0.5,0.5,"Empirical data", cex=2)

plot.new()
text(0.5,0.5,"Empirical data", cex=2)

plot.new()
text(0.5,0.5,"Empirical data", cex=2)

dev.off()

##** FIG S5 ####
#*** Mean Ho maps
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,sub.row=pla$x,sub.col=pla$y,sub1=pla$sub1,sub2=pla$sub2, main=T, h.main=0.5, w.legend=0, h.sub2=0.3, w.axis=0.7, h.axis=0, spacer.sub.col=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS5.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS5.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
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
text(0.5,0.5,"Simulated observed heterozygosity",adj=c(0.5,0.5), cex=3)
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

##** FIG S6 ####
#*** Mean Hs maps
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,sub.row=pla$x,sub.col=pla$y,sub1=pla$sub1,sub2=pla$sub2, main=T, h.main=0.5, w.legend=0, h.sub2=0.3, w.axis=0.7, h.axis=0, spacer.sub.col=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS6.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS6.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
par(mar=c(0,0,0,0))
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,sub.row=pla$x,sub.col=pla$y,sub1=pla$sub1,sub2=pla$sub2, main=T, h.main=0.5, w.legend=0, h.sub2=0.3, w.axis=0.7, h.axis=0, spacer.sub.col=0.5)

for (j in 2:ncol(Hs_modelled)){
  if (internal){
    river_plot(width_country=0.5, lwd_rivers=0.5, xlimit=c(495000,825000), col_water=col_water, north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, plot_rhone=F, lines_rhone=F, plot_ticino=F, lines_ticino=F, plot_inn=F, lines_inn=F, lakes=TRUE, rivers=TRUE, axes="none", river_nr=T)
  }else{
    river_plot(width_country=0.5, lwd_rivers=0.5, xlimit=c(495000,825000), col_water=col_water, north_arrow = F, overview_map = F, scalebar=F, arrows = F, border_outline=F, width_border=2, col_rhine = col_rhine, plot_rhone=F, plot_ticino=F, plot_inn=F, lakes=TRUE, rivers=TRUE, axes="none", river_nr=F)
  }
  values <- (Hs_modelled[,j]-min(Hs_modelled[,j]))/ (max(Hs_modelled[,j])-min(Hs_modelled[,j])) # transform meanAr values to range [0,1] for heatmap plotting
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
  gradient.legend(Hs_modelled[,j],alpha=1, val.midpoint=F, round=2, val.cex=1, val.gap = 0.5, title.gap=0.1, xl = 505000, xr = 635000, yb = 20000, yt = 40000, horizontal=T)
}
plot.new()
text(0.5,0.5,"Simulated expected heterozygosity",adj=c(0.5,0.5), cex=3)
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

##** d: Model performance ####
maxhist_Ar <- ceiling(max(hist_Ar_A,hist_Ar_B))
maxhist_Ho <- ceiling(10*max(hist_Ho_A,hist_Ho_B))/10
maxhist_Hs <- ceiling(10*max(hist_Hs_A,hist_Hs_B))/10
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
d0001_Hs_A <- hist_Hs_A[,(1:6)]
d001_Hs_A <- hist_Hs_A[,c(7:12)]
d01_Hs_A <- hist_Hs_A[,c(13:18)]
d0001_Hs_B <- hist_Hs_B[,(1:6)]
d001_Hs_B <- hist_Hs_B[,c(7:12)]
d01_Hs_B <- hist_Hs_B[,c(13:18)]

diffSPO_d0001both_Ar_A <- c(apply(d0001_Ar_A,2,sum)-apply(d001_Ar_A,2,sum),apply(d0001_Ar_A,2,sum)-apply(d01_Ar_A,2,sum))
diffSPO_d0001both_Ar_B <- c(apply(d0001_Ar_B,2,sum)-apply(d001_Ar_B,2,sum),apply(d0001_Ar_B,2,sum)-apply(d01_Ar_B,2,sum))
diffSPO_d0001both_Ho_A <- c(apply(d0001_Ho_A,2,sum)-apply(d001_Ho_A,2,sum),apply(d0001_Ho_A,2,sum)-apply(d01_Ho_A,2,sum))
diffSPO_d0001both_Ho_B <- c(apply(d0001_Ho_B,2,sum)-apply(d001_Ho_B,2,sum),apply(d0001_Ho_B,2,sum)-apply(d01_Ho_B,2,sum))
diffSPO_d0001both_Hs_A <- c(apply(d0001_Hs_A,2,sum)-apply(d001_Hs_A,2,sum),apply(d0001_Hs_A,2,sum)-apply(d01_Hs_A,2,sum))
diffSPO_d0001both_Hs_B <- c(apply(d0001_Hs_B,2,sum)-apply(d001_Hs_B,2,sum),apply(d0001_Hs_B,2,sum)-apply(d01_Hs_B,2,sum))
diffMPO_d0001both_Ar_A <- c(apply(d0001_Ar_A,2,median)-apply(d001_Ar_A,2,median),apply(d0001_Ar_A,2,median)-apply(d01_Ar_A,2,median))
diffMPO_d0001both_Ar_B <- c(apply(d0001_Ar_B,2,median)-apply(d001_Ar_B,2,median),apply(d0001_Ar_B,2,median)-apply(d01_Ar_B,2,median))
diffMPO_d0001both_Ho_A <- c(apply(d0001_Ho_A,2,median)-apply(d001_Ho_A,2,median),apply(d0001_Ho_A,2,median)-apply(d01_Ho_A,2,median))
diffMPO_d0001both_Ho_B <- c(apply(d0001_Ho_B,2,median)-apply(d001_Ho_B,2,median),apply(d0001_Ho_B,2,median)-apply(d01_Ho_B,2,median))
diffMPO_d0001both_Hs_A <- c(apply(d0001_Hs_A,2,median)-apply(d001_Hs_A,2,median),apply(d0001_Hs_A,2,median)-apply(d01_Hs_A,2,median))
diffMPO_d0001both_Hs_B <- c(apply(d0001_Hs_B,2,median)-apply(d001_Hs_B,2,median),apply(d0001_Hs_B,2,median)-apply(d01_Hs_B,2,median))

list_d0001both <- c(diffMPO_d0001both_Ar_A,diffMPO_d0001both_Ar_B,
                    diffMPO_d0001both_Ho_A,diffMPO_d0001both_Ho_B,
                    diffMPO_d0001both_Hs_A,diffMPO_d0001both_Hs_B,
                    diffSPO_d0001both_Ar_A,diffSPO_d0001both_Ar_B,
                    diffSPO_d0001both_Ho_A,diffSPO_d0001both_Ho_B,
                    diffSPO_d0001both_Hs_A,diffSPO_d0001both_Hs_B)
total_comparison_d0001both <- length(list_d0001both)
improved_fit_d0001both <- sum(list_d0001both<0)
# Model performance d
improved_fit_d0001both/total_comparison_d0001both

##** FIG S7 ####
##** Perpendicular offset histogram comparison for d
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS7.pdf"), width=8, height=6)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS7.png"), width=8, height=6, units="in", res=300)
}
op <- par(mfrow = c(3,3),
          oma = c(5,5,2,0) + 0.1,
          mar = c(0,0,2,1) + 0.1)
hist(d0001_Ar_A, col=col_Gfos_A, cex.axis=2, xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ar), breaks=seq(0,maxhist_Ar,0.25))
hist(d0001_Ar_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ar,0.25))
axis(1, labels=F, cex.axis=2)
mtext(side = 3, text = lab_Ar_short, line = 1, adj=0.5, cex = 1.5)
abline(v=median(d0001_Ar_A), col=white_transparent, lwd=2)
abline(v=median(d0001_Ar_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(d0001_Ar_B), col=white_transparent, lwd=2)
abline(v=median(d0001_Ar_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(4,65, txt = lab[6], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)
legend(2,50, c(expression("Median" ~ italic(G. ~ fossarum) ~ "type A"),expression("Median" ~ italic(G. ~ fossarum) ~ "type B")),
       lty=2, col=c(col_Gfos_A,col_Gfos_B), bty="n", cex=0.85, lwd=2)

hist(d0001_Ho_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ho), breaks=seq(0,maxhist_Ho,0.025))
hist(d0001_Ho_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ho,0.025))
axis(1, labels=F, cex.axis=2)
axis(2, labels=F, cex.axis=2)
mtext(side = 3, text = lab_Ho_short, line = 1, adj=0.5, cex = 1.5)
abline(v=median(d0001_Ho_A), col=white_transparent, lwd=2)
abline(v=median(d0001_Ho_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(d0001_Ho_B), col=white_transparent, lwd=2)
abline(v=median(d0001_Ho_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,65, txt = lab[6], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

hist(d0001_Hs_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Hs), breaks=seq(0,maxhist_Hs,0.025))
hist(d0001_Hs_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Hs,0.025))
axis(1, labels=F, cex.axis=2)
axis(2, labels=F, cex.axis=2)
mtext(side = 3, text = lab_Hs_short, line = 1, adj=0.5, cex = 1.5)
abline(v=median(d0001_Hs_A), col=white_transparent, lwd=2)
abline(v=median(d0001_Hs_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(d0001_Hs_B), col=white_transparent, lwd=2)
abline(v=median(d0001_Hs_B), col=col_Gfos_B, lty=2, lwd=2)
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

hist(d001_Hs_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Hs), breaks=seq(0,maxhist_Hs,0.025))
hist(d001_Hs_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Hs,0.025))
axis(1, labels=F, cex.axis=2)
axis(2, labels=F, cex.axis=2)
abline(v=median(d001_Hs_A), col=white_transparent, lwd=2)
abline(v=median(d001_Hs_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(d001_Hs_B), col=white_transparent, lwd=2)
abline(v=median(d001_Hs_B), col=col_Gfos_B, lty=2, lwd=2)
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

hist(d01_Hs_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Hs), breaks=seq(0,maxhist_Hs,0.025))
hist(d01_Hs_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Hs,0.025))
axis(1, labels=T, cex.axis=2, padj=0.5)
axis(2, labels=F, cex.axis=2)
abline(v=median(d01_Hs_A), col=white_transparent, lwd=2)
abline(v=median(d01_Hs_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(d01_Hs_B), col=white_transparent, lwd=2)
abline(v=median(d01_Hs_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,65, txt = lab[8], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

title(xlab = "Perpendicular offset",
      ylab = "Frequency",
      outer = TRUE, line = 3.5, cex.lab=2)
dev.off()

##** W: Model performance ####
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
w00_Hs_A <- hist_Hs_A[,c(1,2,7,8,13,14)]
w05_Hs_A <- hist_Hs_A[,c(3,4,9,10,15,16)]
w10_Hs_A <- hist_Hs_A[,c(5,6,11,12,17,18)]
w00_Hs_B <- hist_Hs_B[,c(1,2,7,8,13,14)]
w05_Hs_B <- hist_Hs_B[,c(3,4,9,10,15,16)]
w10_Hs_B <- hist_Hs_B[,c(5,6,11,12,17,18)]

diffSPO_w00both_Ar_A <- c(apply(w00_Ar_A,2,sum)-apply(w05_Ar_A,2,sum),apply(w00_Ar_A,2,sum)-apply(w10_Ar_A,2,sum))
diffSPO_w00both_Ar_B <- c(apply(w00_Ar_B,2,sum)-apply(w05_Ar_B,2,sum),apply(w00_Ar_B,2,sum)-apply(w10_Ar_B,2,sum))
diffSPO_w00both_Ho_A <- c(apply(w00_Ho_A,2,sum)-apply(w05_Ho_A,2,sum),apply(w00_Ho_A,2,sum)-apply(w10_Ho_A,2,sum))
diffSPO_w00both_Ho_B <- c(apply(w00_Ho_B,2,sum)-apply(w05_Ho_B,2,sum),apply(w00_Ho_B,2,sum)-apply(w10_Ho_B,2,sum))
diffSPO_w00both_Hs_A <- c(apply(w00_Hs_A,2,sum)-apply(w05_Hs_A,2,sum),apply(w00_Hs_A,2,sum)-apply(w10_Hs_A,2,sum))
diffSPO_w00both_Hs_B <- c(apply(w00_Hs_B,2,sum)-apply(w05_Hs_B,2,sum),apply(w00_Hs_B,2,sum)-apply(w10_Hs_B,2,sum))
diffMPO_w00both_Ar_A <- c(apply(w00_Ar_A,2,median)-apply(w05_Ar_A,2,median),apply(w00_Ar_A,2,median)-apply(w10_Ar_A,2,median))
diffMPO_w00both_Ar_B <- c(apply(w00_Ar_B,2,median)-apply(w05_Ar_B,2,median),apply(w00_Ar_B,2,median)-apply(w10_Ar_B,2,median))
diffMPO_w00both_Ho_A <- c(apply(w00_Ho_A,2,median)-apply(w05_Ho_A,2,median),apply(w00_Ho_A,2,median)-apply(w10_Ho_A,2,median))
diffMPO_w00both_Ho_B <- c(apply(w00_Ho_B,2,median)-apply(w05_Ho_B,2,median),apply(w00_Ho_B,2,median)-apply(w10_Ho_B,2,median))
diffMPO_w00both_Hs_A <- c(apply(w00_Hs_A,2,median)-apply(w05_Hs_A,2,median),apply(w00_Hs_A,2,median)-apply(w10_Hs_A,2,median))
diffMPO_w00both_Hs_B <- c(apply(w00_Hs_B,2,median)-apply(w05_Hs_B,2,median),apply(w00_Hs_B,2,median)-apply(w10_Hs_B,2,median))

list_w00both <- c(diffMPO_w00both_Ar_A,diffMPO_w00both_Ar_B,
                  diffMPO_w00both_Ho_A,diffMPO_w00both_Ho_B,
                  diffMPO_w00both_Hs_A,diffMPO_w00both_Hs_B,
                  diffSPO_w00both_Ar_A,diffSPO_w00both_Ar_B,
                  diffSPO_w00both_Ho_A,diffSPO_w00both_Ho_B,
                  diffSPO_w00both_Hs_A,diffSPO_w00both_Hs_B)
total_comparison_w00both <- length(list_w00both)
improved_fit_w00both <- sum(list_w00both<0)
# Model performance W
improved_fit_w00both/total_comparison_w00both

##** FIG S8 ####
##** Perpendicular offset histogram comparison for W
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS8.pdf"), width=8, height=6)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS8.png"), width=8, height=6, units="in", res=300)
}
op <- par(mfrow = c(3,3),
          oma = c(5,5,2,0) + 0.1,
          mar = c(0,0,2,1) + 0.1)
hist(w00_Ar_A, col=col_Gfos_A, cex.axis=2, xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ar), breaks=seq(0,maxhist_Ar,0.25))
hist(w00_Ar_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ar,0.25))
axis(1, labels=F, cex.axis=2)
mtext(side = 3, text = lab_Ar_short, line = 1, adj=0.5, cex = 1.5)
abline(v=median(w00_Ar_A), col=white_transparent, lwd=2)
abline(v=median(w00_Ar_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(w00_Ar_B), col=white_transparent, lwd=2)
abline(v=median(w00_Ar_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(4,65, txt = lab[3], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)
legend(2,50, c(expression("Median" ~ italic(G. ~ fossarum) ~ "type A"),expression("Median" ~ italic(G. ~ fossarum) ~ "type B")),
       lty=2, col=c(col_Gfos_A,col_Gfos_B), bty="n", cex=0.85, lwd=2)

hist(w00_Ho_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Ho), breaks=seq(0,maxhist_Ho,0.025))
hist(w00_Ho_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ho,0.025))
axis(1, labels=F, cex.axis=2)
axis(2, labels=F, cex.axis=2)
mtext(side = 3, text = lab_Ho_short, line = 1, adj=0.5, cex = 1.5)
abline(v=median(w00_Ho_A), col=white_transparent, lwd=2)
abline(v=median(w00_Ho_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(w00_Ho_B), col=white_transparent, lwd=2)
abline(v=median(w00_Ho_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,65, txt = lab[3], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

hist(w00_Hs_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Hs), breaks=seq(0,maxhist_Hs,0.025))
hist(w00_Hs_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Hs,0.025))
axis(1, labels=F, cex.axis=2)
axis(2, labels=F, cex.axis=2)
mtext(side = 3, text = lab_Hs_short, line = 1, adj=0.5, cex = 1.5)
abline(v=median(w00_Hs_A), col=white_transparent, lwd=2)
abline(v=median(w00_Hs_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(w00_Hs_B), col=white_transparent, lwd=2)
abline(v=median(w00_Hs_B), col=col_Gfos_B, lty=2, lwd=2)
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

hist(w05_Hs_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Hs), breaks=seq(0,maxhist_Hs,0.025))
hist(w05_Hs_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Hs,0.025))
axis(1, labels=F, cex.axis=2)
axis(2, labels=F, cex.axis=2)
abline(v=median(w05_Hs_A), col=white_transparent, lwd=2)
abline(v=median(w05_Hs_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(w05_Hs_B), col=white_transparent, lwd=2)
abline(v=median(w05_Hs_B), col=col_Gfos_B, lty=2, lwd=2)
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

hist(w10_Hs_A, col=col_Gfos_A, yaxt="n", xaxt="n", main="", ylim=c(0,70), xlim = c(0,maxhist_Hs), breaks=seq(0,maxhist_Hs,0.025))
hist(w10_Hs_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Hs,0.025))
axis(1, labels=T, cex.axis=2, padj=0.5)
axis(2, labels=F, cex.axis=2)
abline(v=median(w10_Hs_A), col=white_transparent, lwd=2)
abline(v=median(w10_Hs_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(w10_Hs_B), col=white_transparent, lwd=2)
abline(v=median(w10_Hs_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,65, txt = lab[5], txt.adj=0.5, txt.cex = 2, frm.brd = NA, frm.col = white_transparent)

strwidth(lab[5]) * 2
title(xlab = "Perpendicular offset",
      ylab = "Frequency",
      outer = TRUE, line = 3.5, cex.lab=2)
dev.off()

##** K: Model performance ####
k0_Ar_A <- hist_Ar_A[,c(1,3,5,7,9,11,13,15,17)]
k1_Ar_A <- hist_Ar_A[,c(2,4,6,8,10,12,14,16,18)]
k0_Ar_B <- hist_Ar_B[,c(1,3,5,7,9,11,13,15,17)]
k1_Ar_B <- hist_Ar_B[,c(2,4,6,8,10,12,14,16,18)]
k0_Ho_A <- hist_Ho_A[,c(1,3,5,7,9,11,13,15,17)]
k1_Ho_A <- hist_Ho_A[,c(2,4,6,8,10,12,14,16,18)]
k0_Ho_B <- hist_Ho_B[,c(1,3,5,7,9,11,13,15,17)]
k1_Ho_B <- hist_Ho_B[,c(2,4,6,8,10,12,14,16,18)]
k0_Hs_A <- hist_Hs_A[,c(1,3,5,7,9,11,13,15,17)]
k1_Hs_A <- hist_Hs_A[,c(2,4,6,8,10,12,14,16,18)]
k0_Hs_B <- hist_Hs_B[,c(1,3,5,7,9,11,13,15,17)]
k1_Hs_B <- hist_Hs_B[,c(2,4,6,8,10,12,14,16,18)]

diffSPO_k0k1_Ar_A <- apply(k0_Ar_A,2,sum)-apply(k1_Ar_A,2,sum)
diffSPO_k0k1_Ar_B <- apply(k0_Ar_B,2,sum)-apply(k1_Ar_B,2,sum)
diffSPO_k0k1_Ho_A <- apply(k0_Ho_A,2,sum)-apply(k1_Ho_A,2,sum)
diffSPO_k0k1_Ho_B <- apply(k0_Ho_B,2,sum)-apply(k1_Ho_B,2,sum)
diffSPO_k0k1_Hs_A <- apply(k0_Hs_A,2,sum)-apply(k1_Hs_A,2,sum)
diffSPO_k0k1_Hs_B <- apply(k0_Hs_B,2,sum)-apply(k1_Hs_B,2,sum)
diffMPO_k0k1_Ar_A <- apply(k0_Ar_A,2,median)-apply(k1_Ar_A,2,median)
diffMPO_k0k1_Ar_B <- apply(k0_Ar_B,2,median)-apply(k1_Ar_B,2,median)
diffMPO_k0k1_Ho_A <- apply(k0_Ho_A,2,median)-apply(k1_Ho_A,2,median)
diffMPO_k0k1_Ho_B <- apply(k0_Ho_B,2,median)-apply(k1_Ho_B,2,median)
diffMPO_k0k1_Hs_A <- apply(k0_Hs_A,2,median)-apply(k1_Hs_A,2,median)
diffMPO_k0k1_Hs_B <- apply(k0_Hs_B,2,median)-apply(k1_Hs_B,2,median)

list_k0k1 <- c(diffMPO_k0k1_Ar_A,diffMPO_k0k1_Ar_B,
               diffMPO_k0k1_Ho_A,diffMPO_k0k1_Ho_B,
               diffMPO_k0k1_Hs_A,diffMPO_k0k1_Hs_B,
               diffSPO_k0k1_Ar_A,diffSPO_k0k1_Ar_B,
               diffSPO_k0k1_Ho_A,diffSPO_k0k1_Ho_B,
               diffSPO_k0k1_Hs_A,diffSPO_k0k1_Hs_B)
total_comparison_k0k1 <- length(list_k0k1)
improved_fit_k0k1 <- sum(list_k0k1<0)
# Model performance K
improved_fit_k0k1/total_comparison_k0k1

##** FIG S9 ####
##** Perpendicular offset histogram comparison for K
scal.fact <- 0.75
if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS9.pdf"), width=8, height=6*scal.fact)
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS9.png"), width=8, height=6*scal.fact, units="in", res=300)
}
op <- par(mfrow = c(3*scal.fact,3),
          oma = c(5*scal.fact,5*scal.fact,2*scal.fact,0) + 0.1,
          mar = c(0,0,2*scal.fact,1*scal.fact) + 0.1)
hist(k0_Ar_A, col=col_Gfos_A, tcl=-0.5*scal.fact, cex.axis=2*scal.fact, xaxt="n", yaxt="n", main="", ylim=c(0,90), xlim = c(0,maxhist_Ar), breaks=seq(0,maxhist_Ar,0.25))
hist(k0_Ar_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ar,0.25))
axis(1, labels=F, cex.axis=2*scal.fact, tcl=-0.5*scal.fact)
axis(2, labels=T, cex.axis=2*scal.fact, tcl=-0.5*scal.fact, padj=0.5)
mtext(side = 3, text = lab_Ar_short, line = 1, adj=0.5, cex = 1.5)
abline(v=median(k0_Ar_A), col=white_transparent, lwd=2)
abline(v=median(k0_Ar_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(k0_Ar_B), col=white_transparent, lwd=2)
abline(v=median(k0_Ar_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(4,85, txt = lab[1], txt.adj=0.5, txt.cex = 2*scal.fact, frm.brd = NA, frm.col = white_transparent)
legend(2,70, c(expression("Median" ~ italic(G. ~ fossarum) ~ "type A"),expression("Median" ~ italic(G. ~ fossarum) ~ "type B")),
       lty=2, col=c(col_Gfos_A,col_Gfos_B), bty="n", cex=0.85, lwd=2)

hist(k0_Ho_A, col=col_Gfos_A, tcl=-0.5*scal.fact, cex.axis=2*scal.fact, xaxt="n", yaxt="n", main="", ylim=c(0,90), xlim = c(0,maxhist_Ho), breaks=seq(0,maxhist_Ho,0.025))
hist(k0_Ho_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Ho,0.025))
axis(1, labels=F, cex.axis=2*scal.fact, tcl=-0.5*scal.fact)
axis(2, labels=F, cex.axis=2*scal.fact, tcl=-0.5*scal.fact)
mtext(side = 3, text = lab_Ho_short, line = 1, adj=0.5, cex = 1.5)
abline(v=median(k0_Ho_A), col=white_transparent, lwd=2)
abline(v=median(k0_Ho_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(k0_Ho_B), col=white_transparent, lwd=2)
abline(v=median(k0_Ho_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,85, txt = lab[1], txt.adj=0.5, txt.cex = 2*scal.fact, frm.brd = NA, frm.col = white_transparent)

hist(k0_Hs_A, col=col_Gfos_A, tcl=-0.5*scal.fact, cex.axis=2*scal.fact, xaxt="n", yaxt="n", main="", ylim=c(0,90), xlim = c(0,maxhist_Hs), breaks=seq(0,maxhist_Hs,0.025))
hist(k0_Hs_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Hs,0.025))
axis(1, labels=F, cex.axis=2*scal.fact, tcl=-0.5*scal.fact)
axis(2, labels=F, cex.axis=2*scal.fact, tcl=-0.5*scal.fact)
mtext(side = 3, text = lab_Hs_short, line = 1, adj=0.5, cex = 1.5)
abline(v=median(k0_Hs_A), col=white_transparent, lwd=2)
abline(v=median(k0_Hs_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(k0_Hs_B), col=white_transparent, lwd=2)
abline(v=median(k0_Hs_B), col=col_Gfos_B, lty=2, lwd=2)
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

hist(k1_Hs_A, col=col_Gfos_A, tcl=-0.5*scal.fact, cex.axis=2*scal.fact, xaxt="n", yaxt="n", main="", ylim=c(0,90), xlim = c(0,maxhist_Hs), breaks=seq(0,maxhist_Hs,0.025))
hist(k1_Hs_B, col=col_Gfos_B, add=T, breaks=seq(0,maxhist_Hs,0.025))
axis(1, labels=T, cex.axis=2*scal.fact, tcl=-0.5*scal.fact, padj=0)
axis(2, labels=F, cex.axis=2*scal.fact, tcl=-0.5*scal.fact)
abline(v=median(k1_Hs_A), col=white_transparent, lwd=2)
abline(v=median(k1_Hs_A), col=col_Gfos_A, lty=2, lwd=2)
abline(v=median(k1_Hs_B), col=white_transparent, lwd=2)
abline(v=median(k1_Hs_B), col=col_Gfos_B, lty=2, lwd=2)
textbox(0.3,85, txt = lab[2], txt.adj=0.5, txt.cex = 2*scal.fact, frm.brd = NA, frm.col = white_transparent)

title(xlab = "Perpendicular offset",
      ylab = "Frequency",
      outer = TRUE, line = 3.5*scal.fact, cex.lab=2*scal.fact)
dev.off()

##** FIG S10 ####
##** Orthogonal distance to 1:1 line
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=pla$main)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS10.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS10.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
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
    mtext(expression(bold("Model data: ")*" Mean allelic richness"), side=2, line=1, cex=1)
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

##** TABLE S11 ####
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
write.csv2(SPOrank_Ar, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S11.csv"))

##** FIG S12 ####
###** Histogram of orthogonal distance to 1:1 line
main_5=F
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_5, h.main=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS12.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS12.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_5, h.main=0.5)

for (i in 1:ncol(hist_Ho_A)){
  par(mar=c(0,0,0,0))
  yclip <- 30
  ylim <- 35
  hist(hist_Ho_A[,i],
       breaks=seq(0,ceiling(max(hist_Ho_A,hist_Ho_B)/0.1)*0.1,0.05),
       xlim=c(0,ceiling(max(hist_Ho_A,hist_Ho_B)/0.1)*0.1),
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
       breaks=seq(0,ceiling(max(hist_Ho_A,hist_Ho_B)/0.1)*0.1,0.05),
       xlim=c(0,ceiling(max(hist_Ho_A,hist_Ho_B)/0.1)*0.1),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="", add=T)
  hist(hist_Ho_B[,i],
       breaks=seq(0,ceiling(max(hist_Ho_A,hist_Ho_B)/0.1)*0.1,0.05),
       xlim=c(0,ceiling(max(hist_Ho_A,hist_Ho_B)/0.1)*0.1),
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
  text(0.5,0.7,paste0(lab_Ho,": Distribution of ",measure3a),adj=c(0.5,0.5),cex=3)
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
##** Orthogonal distance to 1:1 line
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=pla$main)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS13.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS13.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
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

##** TABLE S14 ####
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
write.csv2(SPOrank_Ho, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S14.csv"))

##** FIG S15 ####
##** Histogram of orthogonal distance to 1:1 line
main_s15=F
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s15, h.main=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS15.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS15.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s15, h.main=0.5)

for (i in 1:ncol(hist_Hs_A)){
  par(mar=c(0,0,0,0))
  yclip <- 30
  ylim <- 35
  hist(hist_Hs_A[,i],
       breaks=seq(0,ceiling_dec(max(hist_Hs_A,hist_Hs_B),1),0.05),
       xlim=c(0,ceiling_dec(max(hist_Hs_A,hist_Hs_B),1)),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="")
  clip(0,ceiling_dec(max(hist_Hs_A,hist_Hs_B),1), 0, yclip)
  abline(v=median(hist_Hs_A[,i]),col=col_Gfos_A)
  abline(v=median(hist_Hs_B[,i]),col=col_Gfos_B)
  clip(0,ceiling_dec(max(hist_Hs_A,hist_Hs_B),1), 0, ylim)
  textbox(ceiling_dec(max(hist_Hs_A,hist_Hs_B),1),30,paste0(measure2_short," = ",formatC(round(median(hist_Hs_A[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_A, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  textbox(ceiling_dec(max(hist_Hs_A,hist_Hs_B),1),25,paste0(measure2_short," = ",formatC(round(median(hist_Hs_B[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_B, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  hist(hist_Hs_A[,i],
       breaks=seq(0,ceiling_dec(max(hist_Hs_A,hist_Hs_B),1),0.05),
       xlim=c(0,ceiling_dec(max(hist_Hs_A,hist_Hs_B),1)),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="", add=T)
  hist(hist_Hs_B[,i],
       breaks=seq(0,ceiling_dec(max(hist_Hs_A,hist_Hs_B),1),0.05),
       xlim=c(0,ceiling_dec(max(hist_Hs_A,hist_Hs_B),1)),
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
    axis(1, at=seq(0,ceiling_dec(max(hist_Hs_A,hist_Hs_B),1),0.1), labels=c("0.0","","0.2","","0.4","",""))
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

if(main_s15){
  plot.new()
  text(0.5,0.7,paste0(lab_Hs,": Distribution of ",measure3a),adj=c(0.5,0.5),cex=3)
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

##** FIG S16 ####
##** Orthogonal distance to 1:1 line
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=pla$main)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS16.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS16.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=pla$main)

for (i in 2:ncol(Hs_modelled)){
  Hs_Mod_A <- Hs_modelled[,i][rownames(Hs_modelled)%in%modsite_GfosA]
  Hs_Mod_B <- Hs_modelled[,i][rownames(Hs_modelled)%in%modsite_GfosB]
  par(mar=c(0,0,0,0))
  plot(Hs_Mod_A~meanHs_A_red_updist, type="n",
       xlim=c(min(Hs_modelled[,-1],meanHs_A_red_updist,meanHs_B_updist),max(Hs_modelled[,-1],meanHs_A_red_updist,meanHs_B_updist)),
       ylim=c(min(Hs_modelled[,-1],meanHs_A_red_updist,meanHs_B_updist),max(Hs_modelled[,-1],meanHs_A_red_updist,meanHs_B_updist)),
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
    mtext(expression(bold("Model data:")*" Expected heterozygosity"), side=2, line=1, cex=1)
  }
  if (i %in% c(12,13)){
    mtext(expression(bold("Empirical data:")*" Expected heterozygosity"), side=1, line=3, cex=1)
  }
  abline(0,1, lwd=1, lty=2) # add 1:1 line
  
  points(Hs_Mod_A~meanHs_A_red_updist,col=col_Gfos_A, pch=16)
  points(Hs_Mod_B~meanHs_B_updist,col=col_Gfos_B, pch=16)
  
  for (j in 1:length(meanHs_A_red_updist)){
    point <- cbind(meanHs_A_red_updist,Hs_Mod_A)[j,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    segments(seg[1],seg[2],seg[3],seg[4], col=col_Gfos_A, lty=2, lwd=0.5)
  }
  for (k in 1:length(meanHs_B_updist)){
    point <- cbind(meanHs_B_updist,Hs_Mod_B)[k,]
    seg <- unlist(perp.segment.coord(point[1],point[2]))
    segments(seg[1],seg[2],seg[3],seg[4], col=col_Gfos_B, lty=2, lwd=0.5)
  }
  textbox(max(Hs_modelled[,-1],meanHs_A_red_updist,meanHs_B_updist),min(Hs_modelled[,-1],meanHs_A_red_updist,meanHs_B_updist)+0.18,paste0(measure1_short," = ",formatC(round(sum_orthdist_Hs_A[[i-1]],sum_digits),digits=sum_digits, format="f")), txt.cex=sum_cex, txt.adj=1, txt.col=col_Gfos_A, frm.col=white_transparent, frm.brd = NA, frm.siz = 0.2)
  textbox(max(Hs_modelled[,-1],meanHs_A_red_updist,meanHs_B_updist),min(Hs_modelled[,-1],meanHs_A_red_updist,meanHs_B_updist)+0.05,paste0(measure1_short," = ",formatC(round(sum_orthdist_Hs_B[[i-1]],sum_digits),digits=sum_digits, format="f")), txt.cex=sum_cex, txt.adj=1, txt.col=col_Gfos_B, frm.col=white_transparent, frm.brd = NA, frm.siz = 0.2)
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

##** TABLE S17 ####
SPOrank_Hs_A <- labs_comb[order(apply(hist_Hs_A,2,sum))]
SPOrank_Hs_B <- labs_comb[order(apply(hist_Hs_B,2,sum))]

SPOrank_spo_Hs_A <- apply(hist_Hs_A,2,sum)[order(apply(hist_Hs_A,2,sum))]
SPOrank_d_Hs_A <- labs_short[,3][order(apply(hist_Hs_A,2,sum))]
SPOrank_W_Hs_A <- labs_short[,2][order(apply(hist_Hs_A,2,sum))]
SPOrank_K_Hs_A <- labs_short[,1][order(apply(hist_Hs_A,2,sum))]
SPOrank_spo_Hs_B <- apply(hist_Hs_B,2,sum)[order(apply(hist_Hs_B,2,sum))]
SPOrank_d_Hs_B <- labs_short[,3][order(apply(hist_Hs_B,2,sum))]
SPOrank_W_Hs_B <- labs_short[,2][order(apply(hist_Hs_B,2,sum))]
SPOrank_K_Hs_B <- labs_short[,1][order(apply(hist_Hs_B,2,sum))]

SPOrank_Hs <- data.frame("He_A_SPO"=SPOrank_spo_Hs_A,"He_A_d"=SPOrank_d_Hs_A,"He_A_W"=SPOrank_W_Hs_A,"He_A_K"=SPOrank_K_Hs_A,
                         "He_B_SPO"=SPOrank_spo_Hs_B,"He_B_d"=SPOrank_d_Hs_B,"He_B_W"=SPOrank_W_Hs_B,"He_B_K"=SPOrank_K_Hs_B)
write.csv2(SPOrank_Hs, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S17.csv"))

##** TABLE S18 ####
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
write.csv2(MPOrank_Ar, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S18.csv"))

##** TABLE S19 ####
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
write.csv2(MPOrank_Ho, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S19.csv"))

##** TABLE S20 ####
MPOrank_Hs_A <- labs_comb[order(apply(hist_Hs_A,2,median))]
MPOrank_Hs_B <- labs_comb[order(apply(hist_Hs_B,2,median))]

MPOrank_mpo_Hs_A <- apply(hist_Hs_A,2,median)[order(apply(hist_Hs_A,2,median))]
MPOrank_d_Hs_A <- labs_short[,3][order(apply(hist_Hs_A,2,median))]
MPOrank_W_Hs_A <- labs_short[,2][order(apply(hist_Hs_A,2,median))]
MPOrank_K_Hs_A <- labs_short[,1][order(apply(hist_Hs_A,2,median))]
MPOrank_mpo_Hs_B <- apply(hist_Hs_B,2,median)[order(apply(hist_Hs_B,2,median))]
MPOrank_d_Hs_B <- labs_short[,3][order(apply(hist_Hs_B,2,median))]
MPOrank_W_Hs_B <- labs_short[,2][order(apply(hist_Hs_B,2,median))]
MPOrank_K_Hs_B <- labs_short[,1][order(apply(hist_Hs_B,2,median))]

MPOrank_Hs <- data.frame("He_A_MPO"=MPOrank_mpo_Hs_A,"He_A_d"=MPOrank_d_Hs_A,"He_A_W"=MPOrank_W_Hs_A,"He_A_K"=MPOrank_K_Hs_A,
                         "He_B_MPO"=MPOrank_mpo_Hs_B,"He_B_d"=MPOrank_d_Hs_B,"He_B_W"=MPOrank_W_Hs_B,"He_B_K"=MPOrank_K_Hs_B)
write.csv2(MPOrank_Hs, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S20.csv"))

##** FIG S21 ####
##** Histogram of directed perpendicular offset to 1:1 line
main_s21=F
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s21, h.main=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS21.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS21.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s21, h.main=0.5)

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

if(main_s21){
  plot.new()
  text(0.5,0.7,paste0(lab_Ar,": Distribution of directed ",measure3a),adj=c(0.5,0.5),cex=3)
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

##** FIG S22 ####
##** Histogram of directed perpendicular offset to 1:1 line
main_s22=F
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s22, h.main=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS22.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS22.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s22, h.main=0.5)

for (i in 1:ncol(hist_Ho_A_directed)){
  par(mar=c(0,0,0,0))
  yclip <- 30
  ylim <- 35
  # xlim <- round(max(abs(min(hist_Ho_A_directed,hist_Ho_B_directed)),max(hist_Ho_A_directed,hist_Ho_B_directed)),1)
  xlim <- ceiling(max(hist_Ho_A_directed,hist_Ho_B_directed)/0.1)*0.1
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

if(main_s22){
  plot.new()
  text(0.5,0.7,paste0(lab_Ho,": Distribution of directed ",measure3a),adj=c(0.5,0.5),cex=3)
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

##** FIG S23 ####
##** Histogram of directed perpendicular offset to 1:1 line
main_s23=F
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s23, h.main=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS23.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS23.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s23, h.main=0.5)

for (i in 1:ncol(hist_Hs_A_directed)){
  par(mar=c(0,0,0,0))
  yclip <- 30
  ylim <- 35
  # xlim <- round(max(abs(min(hist_Hs_A_directed,hist_Hs_B_directed)),max(hist_Hs_A_directed,hist_Hs_B_directed)),1)
  xlim <- ceiling(max(hist_Hs_A_directed,hist_Hs_B_directed)/0.1)*0.1
  hist(hist_Hs_A_directed[,i],
       breaks=seq(-xlim,xlim,0.05),
       xlim=c(-xlim,xlim),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="")
  clip(-xlim,xlim, 0, yclip)
  abline(v=0, lwd=1.5, lty=2)
  abline(v=median(hist_Hs_A_directed[,i]),col=col_Gfos_A)
  abline(v=median(hist_Hs_B_directed[,i]),col=col_Gfos_B)
  clip(-xlim,xlim, 0, ylim)
  textbox(xlim,30,paste0(measure4," = ",formatC(round(median(hist_Hs_A_directed[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_A, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  textbox(xlim,25,paste0(measure4," = ",formatC(round(median(hist_Hs_B_directed[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_B, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  hist(hist_Hs_A_directed[,i],
       breaks=seq(-xlim,xlim,0.05),
       xlim=c(-xlim,xlim),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="", add=T)
  hist(hist_Hs_B_directed[,i],
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

if(main_s23){
  plot.new()
  text(0.5,0.7,paste0(lab_Hs,": Distribution of directed ",measure3a),adj=c(0.5,0.5),cex=3)
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

##** FIG S24 ####
# Fst perpendicular offset
##*** Preparing matrices####
orthdist_Fst_A <- vector()
orthdist_Fst_A_directed <- vector()
sum_orthdist_Fst_A <- vector()
hist_Fst_A <- matrix(nrow=length(fst_A_red), ncol=ncol(Ar_modelled)-1)
hist_Fst_A_directed <- matrix(nrow=length(fst_A_red), ncol=ncol(Ar_modelled)-1)

orthdist_Fst_B <- vector()
orthdist_Fst_B_directed <- vector()
sum_orthdist_Fst_B <- vector()
hist_Fst_B <- matrix(nrow=length(fst_B), ncol=ncol(Ar_modelled)-1)
hist_Fst_B_directed <- matrix(nrow=length(fst_B), ncol=ncol(Ar_modelled)-1)

# dir.create(paste0(WD,"/Analysis_",output,"/SuppFigs/Fst"), showWarnings=F)

r <- 0
for (d in 1:length(D)){ # looping over dispersal rates 
  for (w in 1:length(W)){ # looping over dispersal directionalities
    for (k in 1:length(K)){ # looping over carrying capacities
      r <- r+1
      
      if(!existing_data){
        load(paste0(DF,"/02_Data_prep/",prep_folder,"/IndPopGenData_",D[d],"_",W[w],"_",K[k],".Rdata"))
      }else{
        load(paste0(DF,"/02_Data_prep/Fst_data/FstData_",D[d],"_",W[w],"_",K[k],".Rdata"))
      }
      
      # Distance matrix to vector
      MEANFST_Mod <- meanFst_Mod
      rownames(MEANFST_Mod) <- modsite
      colnames(MEANFST_Mod) <- modsite
      
      fst_match_B <- match(modsite[modsite%in%microsite_B],microsite_B)
      MEANFST_Mod_B <- MEANFST_Mod[fst_match_B,fst_match_B]
      meanFst_Mod_B <- MEANFST_Mod_B[upper.tri(MEANFST_Mod_B)]
      
      fst_match_A_red <- match(modsite[modsite%in%microsite_A],microsite_A_red)
      MEANFST_Mod_A_red <- MEANFST_Mod[fst_match_A_red,fst_match_A_red]
      meanFst_Mod_A_red <- MEANFST_Mod_A_red[upper.tri(MEANFST_Mod_A_red)]
      
      seg <- matrix(nrow=length(fst_A_red),ncol=4)
      point <- cbind(fst_A_red,meanFst_Mod_A_red)
      for (j in 1:nrow(point)){
        seg[j,] <- unlist(perp.segment.coord(point[j,1],point[j,2]))
        orthdist_Fst_A[j] <- euc.dist(c(seg[j,1],seg[j,2]),c(seg[j,3],seg[j,4]))
        orthdist_Fst_A_directed[j] <- sign(seg[j,2]-seg[j,4])*orthdist_Fst_A[j]
      }
      
      seg <- matrix(nrow=length(fst_B),ncol=4)
      point <- cbind(fst_B,meanFst_Mod_B)
      for (j in 1:nrow(point)){
        seg[j,] <- unlist(perp.segment.coord(point[j,1],point[j,2]))
        orthdist_Fst_B[j] <- euc.dist(c(seg[j,1],seg[j,2]),c(seg[j,3],seg[j,4]))
        orthdist_Fst_B_directed[j] <- sign(seg[j,2]-seg[j,4])*orthdist_Fst_B[j]
      }
      
      sum_orthdist_Fst_A[[r]] <- sum(orthdist_Fst_A)
      sum_orthdist_Fst_B[[r]] <- sum(orthdist_Fst_B)
      
      hist_Fst_A[,r] <- orthdist_Fst_A
      hist_Fst_B[,r] <- orthdist_Fst_B
      
      hist_Fst_A_directed[,r] <- orthdist_Fst_A_directed
      hist_Fst_B_directed[,r] <- orthdist_Fst_B_directed
      
      label_Mod_short <- paste0("D",D_label[d],"_W",W_label[w],"_K",K_label[k])
      # if(pdf){
      #   pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/Fst/Fig_FstMod_",label_Mod_short,".pdf"), width=6, height=6)
      # }else{
      #   png(paste0(WD,"/Analysis_",output,"/SuppFigs/Fst/Fig_FstMod_",label_Mod_short,".png"), width=6, height=6, units="in", res=300)
      # }
      plot(fst_A_red,meanFst_Mod_A_red,
           xlim=c(min(fst_A_red,meanFst_Mod_A_red,fst_B,meanFst_Mod_B, na.rm=T),
                  max(fst_A_red,meanFst_Mod_A_red,fst_B,meanFst_Mod_B, na.rm=T)),
           ylim=c(min(fst_A_red,meanFst_Mod_A_red,fst_B,meanFst_Mod_B, na.rm=T),
                  max(fst_A_red,meanFst_Mod_A_red,fst_B,meanFst_Mod_B, na.rm=T)),
           col=col_Gfos_A, asp=1)
      points(fst_B,meanFst_Mod_B, col=col_Gfos_B)
      abline(0,1,col="red")
      
      mtext(label_Mod_short)
      # dev.off()
    } # end looping over carrying capacities
  } # end looping over dispersal directionalities
} # end looping over dispersal rates

names(sum_orthdist_Fst_A) <- colnames(Ar_modelled)[-1]
names(sum_orthdist_Fst_B) <- colnames(Ar_modelled)[-1]

#### Spatial distance between populations in simulations
mod_vertices <- match(modsite,V(net)$name)
DIST_Mod <- distances(net, v=V(net)[mod_vertices], to=V(net)[mod_vertices], weights=E(net))

#### Distance matrix to vector
dist_Mod <- DIST_Mod[upper.tri(DIST_Mod)]

mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=pla$main)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS24.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS24.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=pla$main)

r <- 0
for (d in 1:length(D)){ # looping over dispersal rates 
  for (w in 1:length(W)){ # looping over dispersal directionalities
    for (k in 1:length(K)){ # looping over carrying capacities
      r <- r+1
      
      if(!existing_data){
        load(paste0(DF,"/02_Data_prep/",prep_folder,"/IndPopGenData_",D[d],"_",W[w],"_",K[k],".Rdata"))
      }else{
        load(paste0(DF,"/02_Data_prep/Fst_data/FstData_",D[d],"_",W[w],"_",K[k],".Rdata"))
      }
      
      # Distance matrix to vector
      MEANFST_Mod <- meanFst_Mod
      meanFst_Mod <- meanFst_Mod[upper.tri(meanFst_Mod)]
      
      # DISTDATA table construction
      DISTDATA_Mod <- cbind(meanFst_Mod,dist_Mod)
      colnames(DISTDATA_Mod) <- c("fst","dist")
      DISTDATA_Mod <- data.frame(DISTDATA_Mod)
      
      # Order DISTDATA according to dist
      DISTDATA_Mod <- DISTDATA_Mod[ order(DISTDATA_Mod$dist), ]
      
      # Remove NAs
      DISTDATA_Mod <- DISTDATA_Mod[which(!is.na(DISTDATA_Mod$fst)),]
      
      # Prepare DISTDATA with non-zero Fst values
      DISTDATA_Mod$nonneg_fst <- DISTDATA_Mod$fst
      DISTDATA_Mod$nonneg_fst[which(DISTDATA_Mod$fst<0)] <- 0
      
      # Prepare DISTDATA with log-transformed dist values
      DISTDATA_Mod$log_dist <- log(DISTDATA_Mod$dist)
      
      #### Combined GLM of genetic differentiation by instream distance * species (power term)
      power <- seq(0,1,0.01)
      AICpower <- c()
      for (i in 1:length(power)){
        pow.mod.Fst <- lm(nonneg_fst ~ I(dist^power[i]), DISTDATA_Mod, na.action = "na.fail")
        AICpower[i] <- AIC(pow.mod.Fst)
      }
      model <- lm.bind(nonneg_fst ~ I(dist^power[which.min(AICpower)]), DISTDATA_Mod, "fst_power", critval, step=T)
      DISTDATA_Mod <- model[[1]]
      lm_fst_power <- model[[2]]
      slm_fst_power <- model[[3]]
      
      DISTDATA_Mod$spec <- "Mod"
      
      DISTDATA_Mod <- DISTDATA_Mod[,-which(colnames(DISTDATA_Mod)=="log_dist")]
      DISTDATA_temp <- DISTDATA[,match(colnames(DISTDATA_Mod),colnames(DISTDATA))]
      DISTDATA_temp <- rbind(DISTDATA_Mod,DISTDATA_temp)
      
      label_Mod <- paste0("D=",D_label[d],", W_up=",W_label[w],", K=",K_label[k])
      label_Mod_short <- paste0("D",D_label[d],"_W",W_label[w],"_K",K_label[k])
      
      par(mar=c(0,0,0,0))
      x <- "fst"
      y <- "dist"
      dat=DISTDATA_temp
      model="slm_fst_power"
      CI_border = F
      xlabel="Instream distance [km]"
      ylabel=expression('Genetic diff. [Pairwise Nei F'[ST]*']')
      xax="n"
      yax="n"
      bty="o"
      yrange=c(0,1)
      pointtrans = T
      trans=0.2
      trans_mod=0.1
      xrev=F
      axislog=""
      pt.cex=0.5
      lwd=1
      cex.lab=2
      cex.axis=1.5
      legend=F
      cex.legend=1.5
      main=F
      col1=col_Gfos_A
      col2=col_Gfos_B
      col3=colMod
      col1trans <- rgb(col2rgb(col1)[1,]/255,col2rgb(col1)[2,]/255,col2rgb(col1)[3,]/255,trans)
      col2trans <- rgb(col2rgb(col2)[1,]/255,col2rgb(col2)[2,]/255,col2rgb(col2)[3,]/255,trans)
      col3trans <- rgb(col2rgb(col3)[1,]/255,col2rgb(col3)[2,]/255,col2rgb(col3)[3,]/255,trans_mod)
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
        col3point <- col3trans
      }else{
        col1point <- col1
        col2point <- col2
        col3point <- col3
      }
      plot(form, dat, type = "n", las = 1, bty = bty,
           xlab=xlabel,
           ylab=ylabel,
           xlim=xrange,
           ylim=yrange,
           log=axislog,
           xaxt=xax, yaxt=yax, cex.lab=cex.lab, cex.axis=cex.axis)
      polygon(c(rev(DATAordered[,ycol][DATAordered$spec=="A"]), DATAordered[,ycol][DATAordered$spec=="A"]),
              c(rev(DATAordered[,lwrcol][DATAordered$spec=="A"]), DATAordered[,uprcol][DATAordered$spec=="A"]),
              col = col1trans, border = NA)
      polygon(c(rev(DATAordered[,ycol][DATAordered$spec=="B"]), DATAordered[,ycol][DATAordered$spec=="B"]),
              c(rev(DATAordered[,lwrcol][DATAordered$spec=="B"]), DATAordered[,uprcol][DATAordered$spec=="B"]),
              col = col2trans, border = NA)
      polygon(c(rev(DATAordered[,ycol][DATAordered$spec=="Mod"]), DATAordered[,ycol][DATAordered$spec=="Mod"]),
              c(rev(DATAordered[,lwrcol][DATAordered$spec=="Mod"]), DATAordered[,uprcol][DATAordered$spec=="Mod"]),
              col = col3trans, border = NA)
      points(form, data = subset(DATAordered, spec == "A"), pch = 16, col = col1point, cex=pt.cex)
      points(form, data = subset(DATAordered, spec == "B"), pch = 16, col = col2point, cex=pt.cex)
      points(form, data = subset(DATAordered, spec == "Mod"), pch = 16, col = col3point, cex=pt.cex)
      lines(formfit, data = subset(DATAordered, spec == "A"), lwd = lwd, col=col1)
      if(CI_border){lines(formupr, data = subset(DATAordered, spec == "A"), lwd = 2, lty=2, col=col1)}
      if(CI_border){lines(formlwr, data = subset(DATAordered, spec == "A"), lwd = 2, lty=2, col=col1)}
      lines(formfit, data = subset(DATAordered, spec == "B"), lwd = lwd, col=col2)
      if(CI_border){lines(formupr, data = subset(DATAordered, spec == "B"), lwd = 2, lty=2, col=col2)}
      if(CI_border){lines(formlwr, data = subset(DATAordered, spec == "B"), lwd = 2, lty=2, col=col2)}
      lines(formfit, data = subset(DATAordered, spec == "Mod"), lwd = lwd, col=col3)
      if(CI_border){lines(formupr, data = subset(DATAordered, spec == "Mod"), lwd = 2, lty=2, col=col3)}
      if(CI_border){lines(formlwr, data = subset(DATAordered, spec == "Mod"), lwd = 2, lty=2, col=col3)}
      
      if (r%in%c(2,4,6)){
        axis(2, cex.axis=1.5)
      }
      if (r%in%c(5,6,11,12,17,18)){
        axis(1, c(0,50,100,150,200,250), at=c(0,50000,100000,150000,200000,250000), cex.axis=1.5)
      }
      if (r %in% c(13,15,17)){
        axis(4, labels=F)
      }
      if (r %in% c(1,3,5)){
        axis(2, labels=F)
      }
      if (r %in% c(1,2)){
        mtext(lab[6], side=3, line=0.8, cex=1.4)
      }
      if (r %in% c(7,8)){
        mtext(lab[7], side=3, line=0.8, cex=1.4)
      }
      if (r %in% c(13,14)){
        mtext(lab[8], side=3, line=0.8, cex=1.4)
      }
      if (r %in% c(3)){
        mtext(ylabel, side=2, line=1, cex=1)
      }
      if (r %in% c(11,12)){
        mtext(xlabel, side=1, line=3, cex=1)
      }
    } # end looping over carrying capacities
  } # end looping over dispersal directionalities
} # end looping over dispersal rates
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
legend("topleft",c(label_A,label_B, label_mod),pch = 16, col = c(col_Gfos_A,col_Gfos_B, colMod), bty="n", cex=2)
dev.off()

##** FIG S25 ####
##** Histogram of directed perpendicular offset to 1:1 line
main_s25=F
median_cex <- 0.7
mp_dim <- multipanel.dimensions(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s25, h.main=0.5)

if(pdf){
  pdf(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS25.pdf"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]))
}else{
  png(paste0(WD,"/Analysis_",output,"/SuppFigs/FigS25.png"), width=fig.width, height=fig.width*(mp_dim[2]/mp_dim[1]), units="in", res=600)
}
nf <- multipanel.layout(main.col=pla$a,main.row=pla$b,pla$x,pla$y,sub1=pla$sub1,sub2=pla$sub2,main=main_s25, h.main=0.5)

# min(abs(apply(hist_Fst_A_directed, 2, median)))

for (i in 1:ncol(hist_Fst_A)){
  par(mar=c(0,0,0,0))
  yclip <- 450
  ylim <- 500
  hist(hist_Fst_A_directed[,i],
       breaks=seq(floor(min(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),0.05),
       xlim=c(floor(min(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T))),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="")
  clip(floor(min(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)), 0, yclip)
  abline(v=0, lwd=1.5, lty=2)
  abline(v=median(hist_Fst_A_directed[,i]),col=col_Gfos_A)
  abline(v=median(hist_Fst_B_directed[,i]),col=col_Gfos_B)
  clip(floor(min(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)), 0, ylim)
  # textbox(ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),150,paste0(measure4," = ",formatC(round(median(hist_Fst_A_directed[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_A, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  # textbox(ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),50,paste0(measure4," = ",formatC(round(median(hist_Fst_B_directed[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_B, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  hist(hist_Fst_A_directed[,i],
       breaks=seq(floor(min(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),0.05),
       xlim=c(floor(min(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T))),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_A,
       main="", add=T)
  hist(hist_Fst_B_directed[,i],
       breaks=seq(floor(min(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),0.05),
       xlim=c(floor(min(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T))),
       ylim=c(0,ylim),
       xaxt="n",
       yaxt="n",
       col=col_Gfos_B,
       main="",
       add=T)
  textbox(ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),yclip,paste0(measure4," = ",formatC(round(median(hist_Fst_A_directed[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_A, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  textbox(ceiling(max(hist_Fst_A_directed,hist_Fst_B_directed, na.rm=T)),yclip-50,paste0(measure4," = ",formatC(round(median(hist_Fst_B_directed[,i]),median_digits),digits=median_digits, format="f")), txt.cex=median_cex, txt.col=col_Gfos_B, txt.adj=1, frm.col="white", frm.brd = "white", frm.siz = 0.2)
  
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
  '%notin%' <- Negate('%in%')
  if (i %notin% c(5:6,11:12,17:18)){
    axis(1, labels = F)
  }
  if (i %in% c(1,3,5)){
    axis(2, at=seq(0,yclip,200), labels=F)
  }
  if (i %in% c(2,4,6)){
    axis(2, at=seq(0,yclip,200))
  }
  if (i %in% c(13:18)){
    axis(4, at=seq(0,yclip,200), labels=F)
  }
}

if(main_s25){
  plot.new()
  text(0.5,0.7,paste0(lab_Ar,": Distribution of directed ",measure3a),adj=c(0.5,0.5),cex=3)
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
dev.off()

##** TABLE S26 ####
SPOrank_Fst_A <- labs_comb[order(apply(hist_Fst_A,2,sum))]
SPOrank_Fst_B <- labs_comb[order(apply(hist_Fst_B,2,sum))]

SPOrank_spo_Fst_A <- apply(hist_Fst_A,2,sum)[order(apply(hist_Fst_A,2,sum))]
SPOrank_d_Fst_A <- labs_short[,3][order(apply(hist_Fst_A,2,sum))]
SPOrank_W_Fst_A <- labs_short[,2][order(apply(hist_Fst_A,2,sum))]
SPOrank_K_Fst_A <- labs_short[,1][order(apply(hist_Fst_A,2,sum))]
SPOrank_spo_Fst_B <- apply(hist_Fst_B,2,sum)[order(apply(hist_Fst_B,2,sum))]
SPOrank_d_Fst_B <- labs_short[,3][order(apply(hist_Fst_B,2,sum))]
SPOrank_W_Fst_B <- labs_short[,2][order(apply(hist_Fst_B,2,sum))]
SPOrank_K_Fst_B <- labs_short[,1][order(apply(hist_Fst_B,2,sum))]

SPOrank_Fst <- data.frame("Fst_A_SPO"=SPOrank_spo_Fst_A,"Fst_A_d"=SPOrank_d_Fst_A,"Fst_A_W"=SPOrank_W_Fst_A,"Fst_A_K"=SPOrank_K_Fst_A,
                         "Fst_B_SPO"=SPOrank_spo_Fst_B,"Fst_B_d"=SPOrank_d_Fst_B,"Fst_B_W"=SPOrank_W_Fst_B,"Fst_B_K"=SPOrank_K_Fst_B)
write.csv2(SPOrank_Fst, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S26.csv"))

##** TABLE S27 ####
MPOrank_Fst_A <- labs_comb[order(apply(hist_Fst_A,2,median))]
MPOrank_Fst_B <- labs_comb[order(apply(hist_Fst_B,2,median))]

MPOrank_mpo_Fst_A <- apply(hist_Fst_A,2,median)[order(apply(hist_Fst_A,2,median))]
MPOrank_d_Fst_A <- labs_short[,3][order(apply(hist_Fst_A,2,median))]
MPOrank_W_Fst_A <- labs_short[,2][order(apply(hist_Fst_A,2,median))]
MPOrank_K_Fst_A <- labs_short[,1][order(apply(hist_Fst_A,2,median))]
MPOrank_mpo_Fst_B <- apply(hist_Fst_B,2,median)[order(apply(hist_Fst_B,2,median))]
MPOrank_d_Fst_B <- labs_short[,3][order(apply(hist_Fst_B,2,median))]
MPOrank_W_Fst_B <- labs_short[,2][order(apply(hist_Fst_B,2,median))]
MPOrank_K_Fst_B <- labs_short[,1][order(apply(hist_Fst_B,2,median))]

MPOrank_Fst <- data.frame("Fst_A_MPO"=MPOrank_mpo_Fst_A,"Fst_A_d"=MPOrank_d_Fst_A,"Fst_A_W"=MPOrank_W_Fst_A,"Fst_A_K"=MPOrank_K_Fst_A,
                         "Fst_B_MPO"=MPOrank_mpo_Fst_B,"Fst_B_d"=MPOrank_d_Fst_B,"Fst_B_W"=MPOrank_W_Fst_B,"Fst_B_K"=MPOrank_K_Fst_B)
write.csv2(MPOrank_Fst, paste0(WD,"/Analysis_",output,"/SuppFigs/Table_S27.csv"))

##** d: Model performance ####
d0001_Fst_A <- hist_Fst_A[,(1:6)]
d001_Fst_A <- hist_Fst_A[,c(7:12)]
d01_Fst_A <- hist_Fst_A[,c(13:18)]
d0001_Fst_B <- hist_Fst_B[,(1:6)]
d001_Fst_B <- hist_Fst_B[,c(7:12)]
d01_Fst_B <- hist_Fst_B[,c(13:18)]
# Here we calculate the difference to d=0.01 (instead of d=0.001)
diffSPO_d001both_Fst_A <- c(apply(d001_Fst_A,2,sum)-apply(d0001_Fst_A,2,sum),apply(d001_Fst_A,2,sum)-apply(d01_Fst_A,2,sum))
diffSPO_d001both_Fst_B <- c(apply(d001_Fst_B,2,sum)-apply(d0001_Fst_B,2,sum),apply(d001_Fst_B,2,sum)-apply(d01_Fst_B,2,sum))
diffMPO_d001both_Fst_A <- c(apply(d001_Fst_A,2,median)-apply(d0001_Fst_A,2,median),apply(d001_Fst_A,2,median)-apply(d01_Fst_A,2,median))
diffMPO_d001both_Fst_B <- c(apply(d001_Fst_B,2,median)-apply(d0001_Fst_B,2,median),apply(d001_Fst_B,2,median)-apply(d01_Fst_B,2,median))
list_d001both <- c(diffMPO_d001both_Fst_A,diffMPO_d001both_Fst_B,
                    diffSPO_d001both_Fst_A,diffSPO_d001both_Fst_B)
total_comparison_d001both <- length(na.omit(list_d001both))
improved_fit_d001both <- sum(na.omit(list_d001both<0))
# Model performance d
improved_fit_d001both/total_comparison_d001both

##** W: Model performance ####
w00_Fst_A <- hist_Fst_A[,c(1,2,7,8,13,14)]
w05_Fst_A <- hist_Fst_A[,c(3,4,9,10,15,16)]
w10_Fst_A <- hist_Fst_A[,c(5,6,11,12,17,18)]
w00_Fst_B <- hist_Fst_B[,c(1,2,7,8,13,14)]
w05_Fst_B <- hist_Fst_B[,c(3,4,9,10,15,16)]
w10_Fst_B <- hist_Fst_B[,c(5,6,11,12,17,18)]
diffSPO_w00both_Fst_A <- c(apply(w00_Fst_A,2,sum)-apply(w05_Fst_A,2,sum),apply(w00_Fst_A,2,sum)-apply(w10_Fst_A,2,sum))
diffSPO_w00both_Fst_B <- c(apply(w00_Fst_B,2,sum)-apply(w05_Fst_B,2,sum),apply(w00_Fst_B,2,sum)-apply(w10_Fst_B,2,sum))
diffMPO_w00both_Fst_A <- c(apply(w00_Fst_A,2,median)-apply(w05_Fst_A,2,median),apply(w00_Fst_A,2,median)-apply(w10_Fst_A,2,median))
diffMPO_w00both_Fst_B <- c(apply(w00_Fst_B,2,median)-apply(w05_Fst_B,2,median),apply(w00_Fst_B,2,median)-apply(w10_Fst_B,2,median))

list_w00both <- c(diffMPO_w00both_Fst_A,diffMPO_w00both_Fst_B,
                  diffSPO_w00both_Fst_A,diffSPO_w00both_Fst_B)
total_comparison_w00both <- length(na.omit(list_w00both))
improved_fit_w00both <- sum(na.omit(list_w00both<0))
# Model performance W
improved_fit_w00both/total_comparison_w00both

##** K: Model performance ####
k0_Fst_A <- hist_Fst_A[,c(1,3,5,7,9,11,13,15,17)]
k1_Fst_A <- hist_Fst_A[,c(2,4,6,8,10,12,14,16,18)]
k0_Fst_B <- hist_Fst_B[,c(1,3,5,7,9,11,13,15,17)]
k1_Fst_B <- hist_Fst_B[,c(2,4,6,8,10,12,14,16,18)]
diffSPO_k0k1_Fst_A <- apply(k0_Fst_A,2,sum)-apply(k1_Fst_A,2,sum)
diffSPO_k0k1_Fst_B <- apply(k0_Fst_B,2,sum)-apply(k1_Fst_B,2,sum)
diffMPO_k0k1_Fst_A <- apply(k0_Fst_A,2,median)-apply(k1_Fst_A,2,median)
diffMPO_k0k1_Fst_B <- apply(k0_Fst_B,2,median)-apply(k1_Fst_B,2,median)
list_k0k1 <- c(diffMPO_k0k1_Fst_A,diffMPO_k0k1_Fst_B,
               diffSPO_k0k1_Fst_A,diffSPO_k0k1_Fst_B)
total_comparison_k0k1 <- length(na.omit(list_k0k1))
improved_fit_k0k1 <- sum(na.omit(list_k0k1)<0)
# Model performance K
improved_fit_k0k1/total_comparison_k0k1
