start_time = format(Sys.time(), '%X') # get starttime of analysis

### USER SETTINGS ###

VERSION_NR <- 001                    # change for multiple analysis (only for documentation)

work_dir <- ("/home/daniel/Raman_script/")     # Working directory (name of project folder)

data_folder <- "spectra/Data/"              # Name of your Data folder in your working directory
groups <- c("C", "E")                 # Filenames (unique part)
ending <- ".TXT"                     # Filetype
legend <- c("Control", "Etoposide")     # Naming in publication

plot_original <- T                   # plot original spectra: yes = T | no = F
Legend_pos <- "topright"             # position of legend in plots "topleft" or "topright"
num_x <- 1                           # number of plots in one single window: 1 or 2

max_range <- "1700"                  # reduce maximum wavenumber to...
min_range <- "440"                   # reduce minimum wavenumber to...

remove_area <- F                     # remove Area: yes = T | no = F
area <- c("1250", "1450")            # exclude values from-to wavenumber

spectra_per_cell <- 1                # How many spectra were measured per cell. 
                                     # The average spectrum for each cell will be calculated 
                                     # and used for further analysis.

vector_normalize <- T                # Normalize each spectrum to a vector of length 1
select_outlier <- T                  # select outlier during PCA: yes = T |no = F
PC <- 8                              # Number of components in variance plots
PCs <- 4                             # Number of PCs in overview plots
PC_PCA <- 10                         # number of PCs calculated in PCA

plot_loadings <- T                   # show loading plots
loadings <- 3                        # number of loadings plotted simultaneously
x_PC <- 1                            # PC on x-axis of score plots
y_PC <- 2                            # PC on y-axis of score plots

perform_LDA <- T                     # perform linear discriminant analysis
repetitions <- 100                   # Number of repetitions of the cross-validation
segments.out <- 3                    # Number of segments in the outer loop
segments.in <- 5                     # Number of segments in the inner loop
max.PCs <- 25                        # Max number of PCs to use for LDA
pos.class <- "E"

### Choose Pretreatments

# Available Methods:
# ================== 
#  1. Standard Normal Variate
#  2. Detrend
#  3. Asymmetric Least Squares
#  4. FillPeaks
#  5. Iterative Restricted Least Squares
#  6. Median Window
#  7. Modified Polynomial Fitting
#  8. Simultaneous Peak Detection and Baseline Correction
#  9. Robust Baseline Estimation
# 10. Rolling Ball

# Enter the numbers of the treatments that should be performed here:
pretreatments <- c(4)

### Load packages
Packages <- c("prospectr","baseline","knitr","ggplot2","ggpubr","ggpmisc",
              "matrixStats","stringr","MASS","caret","ROCR", "binom", "cvAUC")

for (p in Packages) {
  library(p, character.only = TRUE)
}

### Function for import
# work_dir: Path to files
# data_folder: folder with data
# groups: String vector of unique groupnames
# ending: Filetype
# Maximum 8 Groups possible 

Import.data <- function(work_dir, data_folder, groups, ending) {
  
  location <- paste(work_dir, data_folder, sep="")
  
  Files_full <- list.files(path = location, pattern = ending)
  Files_full <- str_replace(Files_full, ending, "")
  x <- c()
  for (group in groups) {
    x <- c(x,grep(pattern = group, x = Files_full))
  }
  Files <- paste(Files_full[x], ending, sep="")
  
  import <- function(data) {
    setwd(location)
    df <- try(read.csv(data, header = FALSE, 
                       sep = ",",              
                       dec = "." )) 
    setwd(work_dir)
    return(df)  
  }
  
  # read files
  Raw.list <- lapply(Files,import)
  #combine list of dataframes
  Raw.data <- do.call("cbind", Raw.list)
  # remove wavenumbers
  Spectra <- as.data.frame(t(Raw.data[,-c(which(colnames(Raw.data) == "V1"))]))
  # extract wavenumbers and set to column names
  Wavenumber <- as.numeric(Raw.data[,1])
  colnames(Spectra) <- Wavenumber
  # check for groupvector
  if(groups[1] != Files[1]){
    # create groupvector
    groups_v <- c(1:length(Files))
    for (i in 1:length(groups)){
      Pos <- c(grep(pattern = as.character(groups[i]), x = str_replace(Files, ending, ""))) # check for parts of filenames
      groups_v[Pos] <- groups[i] # fill in vector
    }
  }else{
    groups_v <- sub(pattern=ending, x=groups, replacement="")
  }
  
  #create a list of groups, wavenumbers and spectra
  OriginalData <- list(Wavenumber,Spectra,as.factor(groups_v),Files)
  names(OriginalData) <- c("Wavenumber","Spectra","Groups","Files")
  rownames(OriginalData$Spectra) <- 1:length(OriginalData$Groups)
  
  return(OriginalData)
}

### Import data
Data <- Import.data(work_dir = work_dir, data_folder = data_folder, groups = groups, ending = ending)

if (spectra_per_cell > 1) {
  # Consistency check groups 
  if (any(tapply(Data$Groups, 
             c(gl(length(Data$Groups), spectra_per_cell, length(Data$Groups))), 
             function(x) {(length(unique(x)) != 1)}))) {
    warning("At least one average spectrum contains multiple groups. Only the first group is retained. Please check the sample groups again")
  }
  
  
  Data$Spectra <- aggregate(Data$Spectra, 
                            by=list(gl(nrow(Data$Spectra), 
                                       spectra_per_cell, 
                                       nrow(Data$Spectra))),
                            FUN=sum)[-1]
  Data$Groups <- Data$Groups[seq(1, length(Data$Groups), spectra_per_cell)]
  
  Data$Files <- Data$Files[seq(1, length(Data$Files), spectra_per_cell)]
}

### Function for Spectra plots
# Liste = list with Wavenumber, Spectra, Groups
# Spektren = which Spectra should be plotted
# Wellenzahl = Vektor with wavenumbers
# area = which area should be plotted
# Code = vector for colours 
# Bereich = select wavenumbers

plot.spectra <- function(Liste,
                         Spektren,
                         Wellenzahl = Liste$Wavenumber,
                         area = c(max(Wellenzahl),min(Wellenzahl)), 
                         Code = Liste$Groups, 
                         Bereich = Sample,
                         main = main){
  
  if(Spektren == "Spectra_d"){          # check for derivatives
    Wellenzahl = Liste$Wavenumber_d     
    Bereich = Sample_d
  }
  
  matplot(x = Wellenzahl[Bereich],
          y = t(as.matrix(Liste[[Spektren]][,Bereich])),
          lty = 1, 
          type = "l", 
          col = Code,
          main = main,
          xlab = "Wavenumber [1/cm]",
          ylab = "Intensity [-]",
          font = 2, 
          font.lab = 2,  ###font = 2 ist Fett gedruckt
          lab = c(20,15,10), 
          xlim = area,
          ylim = c(min(Liste[[Spektren]][,Bereich]),
                   max(Liste[[Spektren]][,Bereich])), 
          bty = "l", 
          family = "sans", 
          xaxs = "i")
  grid(lwd = 0.8)
  
  
}

### Plot of the original spectra
if (plot_original == T) {
  
  plot.spectra(Liste = Data,
               Spektren = "Spectra",
               Bereich = c(1:length(Data$Wavenumber)),
               area = c(max(Data$Wavenumber),min(Data$Wavenumber)),
               main = paste("Original Spectra: ","max Wavenumber = ",
                            max(Data$Wavenumber),", min Wavenumber = ",
                            min(Data$Wavenumber), 
                            sep=""))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
}

# Plot with reduced range
plot.spectra(Liste = Data,
             Spektren = "Spectra", 
             Bereich = c(1:length(Data$Wavenumber)),
             main = "Reduce Wavenumber Range")

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")

abline(v = as.numeric(max_range), lty = "dashed")
abline(v = as.numeric(max(Data$Wavenumber)), lty = "dashed")
abline(v = as.numeric(min(Data$Wavenumber)), lty = "dashed")
abline(v = as.numeric(min_range), lty = "dashed")

arrows(x0 = as.numeric(max_range), 
       y0 = 0, 
       x1 = as.numeric(max(Data$Wavenumber)), 
       y1 = 0, 
       col = "blue", 
       lwd = 2, 
       length = 0.1, 
       code = 3)
arrows(x0 = as.numeric(min(Data$Wavenumber)), 
       y0 = 0, 
       x1 = as.numeric(min_range), 
       y1 = 0, 
       col = "blue", 
       lwd = 2, 
       length = 0.1, 
       code = 3)

Data$Spectra_min_max <- subset(Data$Spectra, 
                               select = -c(which(colnames(Data$Spectra)==max_range):
                                           which(colnames(Data$Spectra)==max(Data$Wavenumber)),
                                           which(colnames(Data$Spectra)==min(Data$Wavenumber)):
                                           which(colnames(Data$Spectra)==min_range)))

Data$Wavenumber_min_max <- as.numeric(colnames(Data$Spectra_min_max))

plot.spectra(Liste = Data,
             Spektren = "Spectra_min_max",
             Bereich = c(1:length(Data$Wavenumber_min_max)),
             Wellenzahl = Data$Wavenumber_min_max, 
             area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)), 
             main = paste("Reduced Wavenumber Range from", max_range, "to", min_range))

legend(Legend_pos, 
       legend = legend, 
       pch = 16, 
       col = unique(Data$Groups), 
       inset = 0.05, 
       bty = "n")


### Plot with excluded area
if (remove_area == T) {
  plot.spectra(Liste = Data, 
               Spektren = "Spectra_min_max",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               Wellenzahl = Data$Wavenumber_min_max,
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = "Exclude Wavenumber Area")
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  abline(v = as.numeric(area[1]), 
         lty = "dashed")
  abline(v = as.numeric(area[2]), 
         lty = "dashed")
  
  arrows(x0 = as.numeric(area[1]), 
         y0 = 0, 
         x1 = as.numeric(area[2]), 
         y1 = 0, 
         col = "blue", 
         lwd = 2, 
         length = 0.1, 
         code = 3)
  
  Data$Spectra_area <- subset(Data$Spectra_min_max, 
                                     select = -c(which(colnames(Data$Spectra_min_max)==area[1]):
                                                   which(colnames(Data$Spectra_min_max)==area[2])))
  
  Data$Wavenumber_area <- as.numeric(colnames(Data$Spectra_area))
  
  plot.spectra(Liste = Data,
               Spektren = "Spectra_area",
               Bereich = c(1:length(Data$Wavenumber_area)),
               Wellenzahl = Data$Wavenumber_area, 
               area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)), 
               main = paste("Spectra with excluded Area from", area[2], "to", area[1]))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
}

### Data pre-treatment
# Data types
Type_rw <- c("Spectra_SNV_rw", "Spectra_det_rw", "Spectra_als_rw", "Spectra_fil_rw","Spectra_irl_rw", 
          "Spectra_med_rw", "Spectra_mod_rw", "Spectra_pea_rw", "Spectra_rbe_rw", "Spectra_rol_rw")

Type_a <- c("Spectra_SNV_area", "Spectra_det_area", "Spectra_als_area", "Spectra_fil_area","Spectra_irl_area", 
          "Spectra_med_area", "Spectra_mod_area", "Spectra_pea_area", "Spectra_rbe_area", "Spectra_rol_area")
# Data names for plots
Names <- c("SNV", "Detrend", "Asymmetric Least Squares", "FillPeaks", "Iterative Restricted Least Squares", 
           "Median window", "Modified polynomial fitting", "Simultaneous Peak Detection and Baseline Correction", 
           "Robust Baseline Estimation", "Rolling ball")

### Standard Normal Variate - REDUCED WAVENUMBER

if (1 %in% pretreatments) {
  
  par(mfrow = c(num_x,1))
  Data$Spectra_SNV_rw <- standardNormalVariate(Data$Spectra_min_max)
  
  # calculate mean of spectra
  Spectra_SNV_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_SNV_mean_rw <- rbind(Spectra_SNV_mean_rw, t(colMeans(as.matrix(Data$Spectra_SNV_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_SNV_mean_rw) <- groups
  Data$Spectra_SNV_mean_rw <- Spectra_SNV_mean_rw
  
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_SNV_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (mean) corrected with Standard Normal Variate, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_SNV_mean_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  # calculate median of spectra
  Spectra_SNV_median_rw <- data.frame()
  for (name in groups) {
    Spectra_SNV_median_rw <- rbind(Spectra_SNV_median_rw, t(colMedians(as.matrix(Data$Spectra_SNV_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_SNV_median_rw) <- groups
  Data$Spectra_SNV_median_rw <- Spectra_SNV_median_rw
  
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_SNV_median_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (median) corrected with Standard Normal Variate, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_SNV_median_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### Standard Normal Variate - REDUCED WAVENUMBER/EXCLUDED AREA
  if (remove_area == T) {
    
    Data$Spectra_SNV_area <- standardNormalVariate(Data$Spectra_area)
    
    # calculate mean of spectra
    Spectra_SNV_mean_area <- data.frame()
    for (name in groups) {
      Spectra_SNV_mean_area <- rbind(Spectra_SNV_mean_area, t(colMeans(as.matrix(Data$Spectra_SNV_area[Data$Groups==name,]))))
    }
    rownames(Spectra_SNV_mean_area) <- groups
    Data$Spectra_SNV_mean_area <- Spectra_SNV_mean_area
    
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_SNV_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (mean) corrected with Standard Normal Variate, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_SNV_mean_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    # calculate median of spectra
    Spectra_SNV_median_area <- data.frame()
    for (name in groups) {
      Spectra_SNV_median_area <- rbind(Spectra_SNV_median_area, t(colMedians(as.matrix(Data$Spectra_SNV_area[Data$Groups==name,]))))
    }
    rownames(Spectra_SNV_median_area) <- groups
    Data$Spectra_SNV_median_area <- Spectra_SNV_median_area
    
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_SNV_median_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (median) corrected with Standard Normal Variate, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_SNV_median_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    par(mfrow = c(1,1))
  }
  
}


### Detrend - REDUCED WAVENUMBER

if (2 %in% pretreatments) {
  
  par(mfrow = c(num_x,1))
  
  Data$Spectra_det_rw <- detrend(Data$Spectra_min_max, wav = Data$Wavenumber_min_max)
  
  # calculate mean of spectra
  Spectra_det_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_det_mean_rw <- rbind(Spectra_det_mean_rw, t(colMeans(as.matrix(Data$Spectra_det_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_det_mean_rw) <- groups
  Data$Spectra_det_mean_rw <- Spectra_det_mean_rw
  
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_det_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (mean) corrected with Detrend, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_det_mean_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  # calculate median of spectra
  Spectra_det_median_rw <- data.frame()
  for (name in groups) {
    Spectra_det_median_rw <- rbind(Spectra_det_median_rw, t(colMeans(as.matrix(Data$Spectra_det_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_det_median_rw) <- groups
  Data$Spectra_det_median_rw <- Spectra_det_median_rw
  
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_det_median_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (median) corrected with Detrend, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_det_median_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  rm(Spectra_det_mean_rw, Spectra_det_median_rw)
  
  ### Detrend - REDUCED WAVENUMBER/EXCLUDED AREA
  if (remove_area == T) {
    
    
    Data$Spectra_det_area <- detrend(Data$Spectra_area, wav = Data$Wavenumber_area)
    
    # calculate mean of spectra
    Spectra_det_mean_area <- data.frame()
    for (name in groups) {
      Spectra_det_mean_area <- rbind(Spectra_det_mean_area, t(colMeans(as.matrix(Data$Spectra_det_area[Data$Groups==name,]))))
    }
    rownames(Spectra_det_mean_area) <- groups
    Data$Spectra_det_mean_area <- Spectra_det_mean_area
    
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_det_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (mean) corrected with Detrend, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_det_mean_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    # calculate median of spectra
    Spectra_det_median_area <- data.frame()
    for (name in groups) {
      Spectra_det_median_area <- rbind(Spectra_det_median_area, t(colMeans(as.matrix(Data$Spectra_det_area[Data$Groups==name,]))))
    }
    rownames(Spectra_det_median_area) <- groups
    Data$Spectra_det_median_area <- Spectra_det_median_area
    
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_det_median_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (median) corrected with Detrend, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_det_median_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    par(mfrow = c(1,1))
  }
  
}


### Asymmetric Least Squares - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN

if (3 %in% pretreatments) {
  
  # calculate baseline with "baseline" package
  bc_als_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                        method = "als",  
                        lambda = 5, 
                        p = 0.05, 
                        maxit = 20)
  
  # extract data from baseline object
  Data$Spectra_als_orig_rw <- as.data.frame(getSpectra(bc_als_rw))
  colnames(bc_als_rw@baseline) <- colnames(bc_als_rw@spectra)
  Data$Spectra_als_bl_rw <- as.data.frame(getBaseline(bc_als_rw))
  Data$Spectra_als_rw <- as.data.frame(getCorrected(bc_als_rw))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced spectra
  
  Spectra_orig_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_als_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_rw) <- groups
  Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw
  
  ### calculate median of reduced spectra
  
  Spectra_orig_median_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_als_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_rw) <- groups
  Data$Spectra_orig_median_rw <- Spectra_orig_median_rw
  
  ### calculate mean of baseline
  
  Spectra_als_mean_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_als_mean_bl_rw <- rbind(Spectra_als_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_als_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_als_mean_bl_rw) <- groups
  Data$Spectra_als_mean_bl_rw <- Spectra_als_mean_bl_rw
  
  
  ### calculate median of baseline
  
  Spectra_als_median_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_als_median_bl_rw <- rbind(Spectra_als_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_als_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_als_median_bl_rw) <- groups
  Data$Spectra_als_median_bl_rw <- Spectra_als_median_bl_rw
  
  
  
  ### combine mean of reduced spectra and baseline
  Data$Spectra_orig_and_bl_mean_rw <- data.frame()
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_als_mean_bl_rw)
  
  # plot mean of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Mean of original spectra and calculated baseline with Asymmetric Least Squares, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced spectra and baseline
  Data$Spectra_orig_and_bl_med_rw <- data.frame()
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_als_median_bl_rw)
  
  # plot median of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_med_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Median of original spectra and calculated baseline with Asymmetric Least Squares, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))#rownames(Data$Spectra_orig_and_bl_rw)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced spectra
  
  Spectra_als_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_als_mean_rw <- rbind(Spectra_als_mean_rw, t(colMeans(as.matrix(Data$Spectra_als_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_als_mean_rw) <- groups
  Data$Spectra_als_mean_rw <- Spectra_als_mean_rw
  
  # plot mean of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_als_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (mean) corrected with Asymmetric Least Squares, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_als_mean_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced spectra
  
  Spectra_als_median_rw <- data.frame()
  for (name in groups) {
    Spectra_als_median_rw <- rbind(Spectra_als_median_rw, t(colMedians(as.matrix(Data$Spectra_als_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_als_median_rw) <- groups
  Data$Spectra_als_median_rw <- Spectra_als_median_rw
  
  # plot median of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_als_median_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (median) corrected with Asymmetric Least Squares, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_als_median_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
  
  ### Asymmetric Least Squares - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
  if (remove_area == T) {
    
    # calculate baseline with "baseline" package
    bc_als_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                            method = "als",  
                            lambda = 7, 
                            p = 0.05, 
                            maxit = 20)
    
    # extract data from baseline object
    Data$Spectra_als_orig_area <- as.data.frame(getSpectra(bc_als_area))
    colnames(bc_als_area@baseline) <- colnames(bc_als_area@spectra)
    Data$Spectra_als_bl_area <- as.data.frame(getBaseline(bc_als_area))
    Data$Spectra_als_area <- as.data.frame(getCorrected(bc_als_area))
    
    par(mfrow = c(num_x,1))
    
    ### calculate mean of reduced/excluded area spectra
    
    Spectra_orig_mean_area <- data.frame()
    for (name in groups) {
      Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_als_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_mean_area) <- groups
    Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
    
    
    ### calculate median of reduced/excluded area spectra
    
    Spectra_orig_median_area <- data.frame()
    for (name in groups) {
      Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_als_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_median_area) <- groups
    Data$Spectra_orig_median_area <- Spectra_orig_median_area
    
    
    ### calculate mean of baseline
    
    Spectra_als_mean_bl_area <- data.frame()
    for (name in groups) {
      Spectra_als_mean_bl_area <- rbind(Spectra_als_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_als_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_als_mean_bl_area) <- groups
    Data$Spectra_als_mean_bl_area <- Spectra_als_mean_bl_area
    
    
    ### calculate median of baseline
    
    Spectra_als_median_bl_area <- data.frame()
    for (name in groups) {
      Spectra_als_median_bl_area <- rbind(Spectra_als_median_bl_area, t(colMedians(as.matrix(Data$Spectra_als_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_als_median_bl_area) <- groups
    Data$Spectra_als_median_bl_area <- Spectra_als_median_bl_area
    
    
    
    ### combine mean of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_mean_area <- data.frame()
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_als_mean_bl_area)
    
    # plot mean of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Mean of original spectra and calculated baseline with Asymmetric Least Squares, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))#rownames(Data$Spectra_orig_and_bl_area)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### combine median of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_med_area <- data.frame()
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_als_median_bl_area)
    
    # plot median of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_med_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Median of original spectra and calculated baseline with Asymmetric Least Squares, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))#rownames(Data$Spectra_orig_and_bl_area)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    
    ### calculate mean of corrected reduced/excluded area spectra
    
    Spectra_als_mean_area <- data.frame()
    for (name in groups) {
      Spectra_als_mean_area <- rbind(Spectra_als_mean_area, t(colMeans(as.matrix(Data$Spectra_als_area[Data$Groups==name,]))))
    }
    rownames(Spectra_als_mean_area) <- groups
    Data$Spectra_als_mean_area <- Spectra_als_mean_area
    
    # plot mean of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_als_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (mean) corrected with Asymmetric Least Squares, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_als_mean_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### calculate median of corrected reduced/excluded area spectra
    
    Spectra_als_median_area <- data.frame()
    for (name in groups) {
      Spectra_als_median_area <- rbind(Spectra_als_median_area, t(colMedians(as.matrix(Data$Spectra_als_area[Data$Groups==name,]))))
    }
    rownames(Spectra_als_median_area) <- groups
    Data$Spectra_als_median_area <- Spectra_als_median_area
    
    # plot median of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_als_median_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (median) corrected with Asymmetric Least Squares, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_als_median_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    par(mfrow = c(1,1))
  }
  
}

### FillPeaks - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN

if (4 %in% pretreatments) {
  
  # calculate baseline with "baseline" package
  bc_fil_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                        method = "fillPeaks",  
                        lambda=1, 
                        hwi=8, 
                        it=5, 
                        int=500)
  
  # extract data from baseline object
  Data$Spectra_fil_orig_rw <- as.data.frame(getSpectra(bc_fil_rw))
  colnames(bc_fil_rw@baseline) <- colnames(bc_fil_rw@spectra)
  Data$Spectra_fil_bl_rw <- as.data.frame(getBaseline(bc_fil_rw))
  Data$Spectra_fil_rw <- as.data.frame(getCorrected(bc_fil_rw))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced spectra
  
  Spectra_orig_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_fil_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_rw) <- groups
  Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw
  
  ### calculate median of reduced spectra
  
  Spectra_orig_median_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_fil_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_rw) <- groups
  Data$Spectra_orig_median_rw <- Spectra_orig_median_rw
  
  ### calculate mean of baseline
  
  Spectra_fil_mean_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_fil_mean_bl_rw <- rbind(Spectra_fil_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_fil_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_fil_mean_bl_rw) <- groups
  Data$Spectra_fil_mean_bl_rw <- Spectra_fil_mean_bl_rw
  
  
  ### calculate median of baseline
  
  Spectra_fil_median_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_fil_median_bl_rw <- rbind(Spectra_fil_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_fil_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_fil_median_bl_rw) <- groups
  Data$Spectra_fil_median_bl_rw <- Spectra_fil_median_bl_rw
  
  
  
  ### combine mean of reduced spectra and baseline
  Data$Spectra_orig_and_bl_mean_rw <- data.frame()
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_fil_mean_bl_rw)
  
  # plot mean of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Mean of original spectra and calculated baseline with FillPeaks, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced spectra and baseline
  Data$Spectra_orig_and_bl_med_rw <- data.frame()
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_fil_median_bl_rw)
  
  # plot median of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_med_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Median of original spectra and calculated baseline with FillPeaks, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced spectra
  
  Spectra_fil_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_fil_mean_rw <- rbind(Spectra_fil_mean_rw, t(colMeans(as.matrix(Data$Spectra_fil_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_fil_mean_rw) <- groups
  Data$Spectra_fil_mean_rw <- Spectra_fil_mean_rw
  
  # plot mean of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_fil_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (mean) corrected with FillPeaks, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_fil_mean_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  # calculate median of corrected reduced spectra
  
  Spectra_fil_median_rw <- data.frame()
  for (name in groups) {
    Spectra_fil_median_rw <- rbind(Spectra_fil_median_rw, t(colMedians(as.matrix(Data$Spectra_fil_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_fil_median_rw) <- groups
  Data$Spectra_fil_median_rw <- Spectra_fil_median_rw
  
  # plot median of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_fil_median_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (median) corrected with FillPeaks, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_fil_median_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
  
  rm(Spectra_fil_mean_bl_rw, Spectra_fil_mean_rw, 
     Spectra_fil_median_bl_rw, Spectra_fil_median_rw, 
     Spectra_orig_mean_rw, Spectra_orig_median_rw,
     bc_fil_rw)
  
  ### FillPeaks - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
  if (remove_area == T) {
    
    # calculate baseline with "baseline" package
    bc_fil_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                            method = "als",  
                            lambda = 7, 
                            p = 0.05, 
                            maxit = 20)
    
    # extract data from baseline object
    Data$Spectra_fil_orig_area <- as.data.frame(getSpectra(bc_fil_area))
    colnames(bc_fil_area@baseline) <- colnames(bc_fil_area@spectra)
    Data$Spectra_fil_bl_area <- as.data.frame(getBaseline(bc_fil_area))
    Data$Spectra_fil_area <- as.data.frame(getCorrected(bc_fil_area))
    
    par(mfrow = c(num_x,1))
    
    ### calculate mean of reduced/excluded area spectra
    
    Spectra_orig_mean_area <- data.frame()
    for (name in groups) {
      Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_fil_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_mean_area) <- groups
    Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
    
    
    ### calculate median of reduced/excluded area spectra
    
    Spectra_orig_median_area <- data.frame()
    for (name in groups) {
      Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_fil_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_median_area) <- groups
    Data$Spectra_orig_median_area <- Spectra_orig_median_area
    
    
    ### calculate mean of baseline
    
    Spectra_fil_mean_bl_area <- data.frame()
    for (name in groups) {
      Spectra_fil_mean_bl_area <- rbind(Spectra_fil_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_fil_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_fil_mean_bl_area) <- groups
    Data$Spectra_fil_mean_bl_area <- Spectra_fil_mean_bl_area
    
    
    ### calculate median of baseline
    
    Spectra_fil_median_bl_area <- data.frame()
    for (name in groups) {
      Spectra_fil_median_bl_area <- rbind(Spectra_fil_median_bl_area, t(colMedians(as.matrix(Data$Spectra_fil_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_fil_median_bl_area) <- groups
    Data$Spectra_fil_median_bl_area <- Spectra_fil_median_bl_area
    
    
    
    ### combine mean of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_mean_area <- data.frame()
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_fil_mean_bl_area)
    
    # plot mean of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Mean of original spectra and calculated baseline with FillPeaks, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### combine median of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_med_area <- data.frame()
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_fil_median_bl_area)
    
    # plot median of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_med_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Median of original spectra and calculated baseline with FillPeaks, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    
    ### calculate mean of corrected reduced/excluded area spectra
    
    Spectra_fil_mean_area <- data.frame()
    for (name in groups) {
      Spectra_fil_mean_area <- rbind(Spectra_fil_mean_area, t(colMeans(as.matrix(Data$Spectra_fil_area[Data$Groups==name,]))))
    }
    rownames(Spectra_fil_mean_area) <- groups
    Data$Spectra_fil_mean_area <- Spectra_fil_mean_area
    
    # plot mean of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_fil_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (mean) corrected with FillPeaks, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_fil_mean_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### calculate median of corrected reduced/excluded area spectra
    
    Spectra_fil_median_area <- data.frame()
    for (name in groups) {
      Spectra_fil_median_area <- rbind(Spectra_fil_median_area, t(colMedians(as.matrix(Data$Spectra_fil_area[Data$Groups==name,]))))
    }
    rownames(Spectra_fil_median_area) <- groups
    Data$Spectra_fil_median_area <- Spectra_fil_median_area
    
    # plot median of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_fil_median_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (median) corrected with FillPeaks, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_fil_median_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    par(mfrow = c(1,1))
  }
  
}


### Iterative Restricted Least Squares - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN                                #

if (5 %in% pretreatments) {
  
  # calculate baseline with "baseline" package
  bc_irl_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                        method = "irls", 
                        lambda1 = 1, 
                        lambda2 = 5, 
                        maxit = 200, 
                        wi = 0.01)
  
  # extract data from baseline object
  Data$Spectra_irl_orig_rw <- as.data.frame(getSpectra(bc_irl_rw))
  colnames(bc_irl_rw@baseline) <- colnames(bc_irl_rw@spectra)
  Data$Spectra_irl_bl_rw <- as.data.frame(getBaseline(bc_irl_rw))
  Data$Spectra_irl_rw <- as.data.frame(getCorrected(bc_irl_rw))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced spectra
  
  Spectra_orig_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_irl_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_rw) <- groups
  Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw
  
  ### calculate median of reduced spectra
  
  Spectra_orig_median_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_irl_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_rw) <- groups
  Data$Spectra_orig_median_rw <- Spectra_orig_median_rw
  
  ### calculate mean of baseline
  
  Spectra_irl_mean_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_irl_mean_bl_rw <- rbind(Spectra_irl_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_irl_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_irl_mean_bl_rw) <- groups
  Data$Spectra_irl_mean_bl_rw <- Spectra_irl_mean_bl_rw
  
  
  ### calculate median of baseline
  
  Spectra_irl_median_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_irl_median_bl_rw <- rbind(Spectra_irl_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_irl_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_irl_median_bl_rw) <- groups
  Data$Spectra_irl_median_bl_rw <- Spectra_irl_median_bl_rw
  
  
  
  ### combine mean of reduced spectra and baseline
  Data$Spectra_orig_and_bl_mean_rw <- data.frame()
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_irl_mean_bl_rw)
  
  # plot mean of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Mean of original spectra and calculated baseline with Iterative Restricted Least Squares, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced spectra and baseline
  Data$Spectra_orig_and_bl_med_rw <- data.frame()
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_irl_median_bl_rw)
  
  # plot median of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_med_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Median of original spectra and calculated baseline with Iterative Restricted Least Squares, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced spectra
  
  Spectra_irl_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_irl_mean_rw <- rbind(Spectra_irl_mean_rw, t(colMeans(as.matrix(Data$Spectra_irl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_irl_mean_rw) <- groups
  Data$Spectra_irl_mean_rw <- Spectra_irl_mean_rw
  
  # plot mean of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_irl_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (mean) corrected with Iterative Restricted Least Squares, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_irl_mean_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced spectra
  
  Spectra_irl_median_rw <- data.frame()
  for (name in groups) {
    Spectra_irl_median_rw <- rbind(Spectra_irl_median_rw, t(colMedians(as.matrix(Data$Spectra_irl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_irl_median_rw) <- groups
  Data$Spectra_irl_median_rw <- Spectra_irl_median_rw
  
  # plot median of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_irl_median_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (median) corrected with Iterative Restricted Least Squares, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_irl_median_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
  
  ### Iterative Restricted Least Squares - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
  if (remove_area == T) {
    
    # calculate baseline with "baseline" package
    bc_irl_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                            method = "irls", 
                            lambda1 = 5, 
                            lambda2 = 9, 
                            maxit = 200, 
                            wi = 0.02)
    
    # extract data from baseline object
    Data$Spectra_irl_orig_area <- as.data.frame(getSpectra(bc_irl_area))
    colnames(bc_irl_area@baseline) <- colnames(bc_irl_area@spectra)
    Data$Spectra_irl_bl_area <- as.data.frame(getBaseline(bc_irl_area))
    Data$Spectra_irl_area <- as.data.frame(getCorrected(bc_irl_area))
    
    par(mfrow = c(num_x,1))
    
    ### calculate mean of reduced/excluded area spectra
    
    Spectra_orig_mean_area <- data.frame()
    for (name in groups) {
      Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_irl_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_mean_area) <- groups
    Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
    
    
    ### calculate median of reduced/excluded area spectra
    
    Spectra_orig_median_area <- data.frame()
    for (name in groups) {
      Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_irl_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_median_area) <- groups
    Data$Spectra_orig_median_area <- Spectra_orig_median_area
    
    
    ### calculate mean of baseline
    
    Spectra_irl_mean_bl_area <- data.frame()
    for (name in groups) {
      Spectra_irl_mean_bl_area <- rbind(Spectra_irl_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_irl_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_irl_mean_bl_area) <- groups
    Data$Spectra_irl_mean_bl_area <- Spectra_irl_mean_bl_area
    
    
    ### calculate median of baseline
    
    Spectra_irl_median_bl_area <- data.frame()
    for (name in groups) {
      Spectra_irl_median_bl_area <- rbind(Spectra_irl_median_bl_area, t(colMedians(as.matrix(Data$Spectra_irl_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_irl_median_bl_area) <- groups
    Data$Spectra_irl_median_bl_area <- Spectra_irl_median_bl_area
    
    
    
    ### combine mean of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_mean_area <- data.frame()
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_irl_mean_bl_area)
    
    # plot mean of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Mean of original spectra and calculated baseline with Iterative Restricted Least Squares, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### combine median of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_med_area <- data.frame()
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_irl_median_bl_area)
    
    # plot median of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_med_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Median of original spectra and calculated baseline with Iterative Restricted Least Squares, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    
    ### calculate mean of corrected reduced/excluded area spectra
    
    Spectra_irl_mean_area <- data.frame()
    for (name in groups) {
      Spectra_irl_mean_area <- rbind(Spectra_irl_mean_area, t(colMeans(as.matrix(Data$Spectra_irl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_irl_mean_area) <- groups
    Data$Spectra_irl_mean_area <- Spectra_irl_mean_area
    
    # plot mean of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_irl_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (mean) corrected with Iterative Restricted Least Squares, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_irl_mean_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### calculate median of corrected reduced/excluded area spectra
    
    Spectra_irl_median_area <- data.frame()
    for (name in groups) {
      Spectra_irl_median_area <- rbind(Spectra_irl_median_area, t(colMedians(as.matrix(Data$Spectra_irl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_irl_median_area) <- groups
    Data$Spectra_irl_median_area <- Spectra_irl_median_area
    
    # plot median of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_irl_median_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (median) corrected with Iterative Restricted Least Squares, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_irl_median_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    par(mfrow = c(1,1))
  }
  
}


### Median window - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN

if (6 %in% pretreatments) {
 
  # calculate baseline with "baseline" package
  bc_med_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                        method = "medianWindow", 
                        hwm = 100)
  
  # extract data from baseline object
  Data$Spectra_med_orig_rw <- as.data.frame(getSpectra(bc_med_rw))
  colnames(bc_med_rw@baseline) <- colnames(bc_med_rw@spectra)
  Data$Spectra_med_bl_rw <- as.data.frame(getBaseline(bc_med_rw))
  Data$Spectra_med_rw <- as.data.frame(getCorrected(bc_med_rw))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced spectra
  
  Spectra_orig_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_med_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_rw) <- groups
  Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw
  
  ### calculate median of reduced spectra
  
  Spectra_orig_median_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_med_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_rw) <- groups
  Data$Spectra_orig_median_rw <- Spectra_orig_median_rw
  
  ### calculate mean of baseline
  
  Spectra_med_mean_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_med_mean_bl_rw <- rbind(Spectra_med_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_med_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_med_mean_bl_rw) <- groups
  Data$Spectra_med_mean_bl_rw <- Spectra_med_mean_bl_rw
  
  
  ### calculate median of baseline
  
  Spectra_med_median_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_med_median_bl_rw <- rbind(Spectra_med_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_med_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_med_median_bl_rw) <- groups
  Data$Spectra_med_median_bl_rw <- Spectra_med_median_bl_rw
  
  
  
  ### combine mean of reduced spectra and baseline
  Data$Spectra_orig_and_bl_mean_rw <- data.frame()
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_med_mean_bl_rw)
  
  # plot mean of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Mean of original spectra and calculated baseline with Median Window, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced spectra and baseline
  Data$Spectra_orig_and_bl_med_rw <- data.frame()
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_med_median_bl_rw)
  
  # plot median of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_med_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Median of original spectra and calculated baseline with Median Window, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced spectra
  
  Spectra_med_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_med_mean_rw <- rbind(Spectra_med_mean_rw, t(colMeans(as.matrix(Data$Spectra_med_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_med_mean_rw) <- groups
  Data$Spectra_med_mean_rw <- Spectra_med_mean_rw
  
  # plot mean of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_med_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (mean) corrected with Median Window, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_med_mean_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced spectra
  
  Spectra_med_median_rw <- data.frame()
  for (name in groups) {
    Spectra_med_median_rw <- rbind(Spectra_med_median_rw, t(colMedians(as.matrix(Data$Spectra_med_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_med_median_rw) <- groups
  Data$Spectra_med_median_rw <- Spectra_med_median_rw
  
  # plot median of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_med_median_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (median) corrected with Median Window, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_med_median_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
  
  ### Median Window - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
  if (remove_area == T) {
    
    # calculate baseline with "baseline" package
    bc_med_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                            method = "medianWindow", 
                            hwm = 300)
    
    # extract data from baseline object
    Data$Spectra_med_orig_area <- as.data.frame(getSpectra(bc_med_area))
    colnames(bc_med_area@baseline) <- colnames(bc_med_area@spectra)
    Data$Spectra_med_bl_area <- as.data.frame(getBaseline(bc_med_area))
    Data$Spectra_med_area <- as.data.frame(getCorrected(bc_med_area))
    
    par(mfrow = c(num_x,1))
    
    ### calculate mean of reduced/excluded area spectra
    
    Spectra_orig_mean_area <- data.frame()
    for (name in groups) {
      Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_med_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_mean_area) <- groups
    Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
    
    
    ### calculate median of reduced/excluded area spectra
    
    Spectra_orig_median_area <- data.frame()
    for (name in groups) {
      Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_med_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_median_area) <- groups
    Data$Spectra_orig_median_area <- Spectra_orig_median_area
    
    
    ### calculate mean of baseline
    
    Spectra_med_mean_bl_area <- data.frame()
    for (name in groups) {
      Spectra_med_mean_bl_area <- rbind(Spectra_med_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_med_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_med_mean_bl_area) <- groups
    Data$Spectra_med_mean_bl_area <- Spectra_med_mean_bl_area
    
    
    ### calculate median of baseline
    
    Spectra_med_median_bl_area <- data.frame()
    for (name in groups) {
      Spectra_med_median_bl_area <- rbind(Spectra_med_median_bl_area, t(colMedians(as.matrix(Data$Spectra_med_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_med_median_bl_area) <- groups
    Data$Spectra_med_median_bl_area <- Spectra_med_median_bl_area
    
    
    
    ### combine mean of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_mean_area <- data.frame()
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_med_mean_bl_area)
    
    # plot mean of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Mean of original spectra and calculated baseline with Median Window, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### combine median of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_med_area <- data.frame()
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_med_median_bl_area)
    
    # plot median of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_med_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Median of original spectra and calculated baseline with Median Window, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    
    ### calculate mean of corrected reduced/excluded area spectra
    
    Spectra_med_mean_area <- data.frame()
    for (name in groups) {
      Spectra_med_mean_area <- rbind(Spectra_med_mean_area, t(colMeans(as.matrix(Data$Spectra_med_area[Data$Groups==name,]))))
    }
    rownames(Spectra_med_mean_area) <- groups
    Data$Spectra_med_mean_area <- Spectra_med_mean_area
    
    # plot mean of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_med_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (mean) corrected with Median Window, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_med_mean_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### calculate median of corrected reduced/excluded area spectra
    
    Spectra_med_median_area <- data.frame()
    for (name in groups) {
      Spectra_med_median_area <- rbind(Spectra_med_median_area, t(colMedians(as.matrix(Data$Spectra_med_area[Data$Groups==name,]))))
    }
    rownames(Spectra_med_median_area) <- groups
    Data$Spectra_med_median_area <- Spectra_med_median_area
    
    # plot median of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_med_median_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (median) corrected with Median Window, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_med_median_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    par(mfrow = c(1,1))
  }
   
}

### Modified Polynomial Fitting - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN

if (7 %in% pretreatments) {
 
  # calculate baseline with "baseline" package
  bc_mod_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                        method = "modpolyfit", 
                        deg = 15,
                        tol = 0.001,
                        rep = 100)
  
  # extract data from baseline object
  Data$Spectra_mod_orig_rw <- as.data.frame(getSpectra(bc_mod_rw))
  colnames(bc_mod_rw@baseline) <- colnames(bc_mod_rw@spectra)
  Data$Spectra_mod_bl_rw <- as.data.frame(getBaseline(bc_mod_rw))
  Data$Spectra_mod_rw <- as.data.frame(getCorrected(bc_mod_rw))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced spectra
  
  Spectra_orig_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_mod_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_rw) <- groups
  Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw
  
  ### calculate median of reduced spectra
  
  Spectra_orig_median_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_mod_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_rw) <- groups
  Data$Spectra_orig_median_rw <- Spectra_orig_median_rw
  
  ### calculate mean of baseline
  
  Spectra_mod_mean_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_mod_mean_bl_rw <- rbind(Spectra_mod_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_mod_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_mod_mean_bl_rw) <- groups
  Data$Spectra_mod_mean_bl_rw <- Spectra_mod_mean_bl_rw
  
  
  ### calculate median of baseline
  
  Spectra_mod_median_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_mod_median_bl_rw <- rbind(Spectra_mod_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_mod_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_mod_median_bl_rw) <- groups
  Data$Spectra_mod_median_bl_rw <- Spectra_mod_median_bl_rw
  
  
  
  ### combine mean of reduced spectra and baseline
  Data$Spectra_orig_and_bl_mean_rw <- data.frame()
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_mod_mean_bl_rw)
  
  # plot mean of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Mean of original spectra and calculated baseline with Modified Polynomial Fitting, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced spectra and baseline
  Data$Spectra_orig_and_bl_med_rw <- data.frame()
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_mod_median_bl_rw)
  
  # plot median of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_med_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Median of original spectra and calculated baseline with Modified Polynomial Fitting, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced spectra
  
  Spectra_mod_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_mod_mean_rw <- rbind(Spectra_mod_mean_rw, t(colMeans(as.matrix(Data$Spectra_mod_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_mod_mean_rw) <- groups
  Data$Spectra_mod_mean_rw <- Spectra_mod_mean_rw
  
  # plot mean of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_mod_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (mean) corrected with Modified Polynomial Fitting, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_mod_mean_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced spectra
  
  Spectra_mod_median_rw <- data.frame()
  for (name in groups) {
    Spectra_mod_median_rw <- rbind(Spectra_mod_median_rw, t(colMedians(as.matrix(Data$Spectra_mod_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_mod_median_rw) <- groups
  Data$Spectra_mod_median_rw <- Spectra_mod_median_rw
  
  # plot median of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_mod_median_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (median) corrected with Modified Polynomial Fitting, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_mod_median_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
  
  ### Modified Polynomial Fitting - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
  if (remove_area == T) {
    
    # calculate baseline with "baseline" package
    bc_mod_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                            method = "modpolyfit", 
                            deg = 8)
    
    # extract data from baseline object
    Data$Spectra_mod_orig_area <- as.data.frame(getSpectra(bc_mod_area))
    colnames(bc_mod_area@baseline) <- colnames(bc_mod_area@spectra)
    Data$Spectra_mod_bl_area <- as.data.frame(getBaseline(bc_mod_area))
    Data$Spectra_mod_area <- as.data.frame(getCorrected(bc_mod_area))
    
    par(mfrow = c(num_x,1))
    
    ### calculate mean of reduced/excluded area spectra
    
    Spectra_orig_mean_area <- data.frame()
    for (name in groups) {
      Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_mod_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_mean_area) <- groups
    Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
    
    
    ### calculate median of reduced/excluded area spectra
    
    Spectra_orig_median_area <- data.frame()
    for (name in groups) {
      Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_mod_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_median_area) <- groups
    Data$Spectra_orig_median_area <- Spectra_orig_median_area
    
    
    ### calculate mean of baseline
    
    Spectra_mod_mean_bl_area <- data.frame()
    for (name in groups) {
      Spectra_mod_mean_bl_area <- rbind(Spectra_mod_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_mod_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_mod_mean_bl_area) <- groups
    Data$Spectra_mod_mean_bl_area <- Spectra_mod_mean_bl_area
    
    
    ### calculate median of baseline
    
    Spectra_mod_median_bl_area <- data.frame()
    for (name in groups) {
      Spectra_mod_median_bl_area <- rbind(Spectra_mod_median_bl_area, t(colMedians(as.matrix(Data$Spectra_mod_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_mod_median_bl_area) <- groups
    Data$Spectra_mod_median_bl_area <- Spectra_mod_median_bl_area
    
    
    
    ### combine mean of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_mean_area <- data.frame()
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_mod_mean_bl_area)
    
    # plot mean of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Mean of original spectra and calculated baseline with Modified Polynomial Fitting, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### combine median of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_med_area <- data.frame()
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_mod_median_bl_area)
    
    # plot median of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_med_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Median of original spectra and calculated baseline with Modified Polynomial Fitting, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    
    ### calculate mean of corrected reduced/excluded area spectra
    
    Spectra_mod_mean_area <- data.frame()
    for (name in groups) {
      Spectra_mod_mean_area <- rbind(Spectra_mod_mean_area, t(colMeans(as.matrix(Data$Spectra_mod_area[Data$Groups==name,]))))
    }
    rownames(Spectra_mod_mean_area) <- groups
    Data$Spectra_mod_mean_area <- Spectra_mod_mean_area
    
    # plot mean of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_mod_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (mean) corrected with Modified Polynomial Fitting, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_mod_mean_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### calculate median of corrected reduced/excluded area spectra
    
    Spectra_mod_median_area <- data.frame()
    for (name in groups) {
      Spectra_mod_median_area <- rbind(Spectra_mod_median_area, t(colMedians(as.matrix(Data$Spectra_mod_area[Data$Groups==name,]))))
    }
    rownames(Spectra_mod_median_area) <- groups
    Data$Spectra_mod_median_area <- Spectra_mod_median_area
    
    # plot median of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_mod_median_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (median) corrected with Modified Polynomial Fitting, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_mod_median_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    par(mfrow = c(1,1))
  }
   
}


### Simultaneous Peak Detection and Baseline Correction - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN

if (8 %in% pretreatments) {
  
  options(warn=-1) # suppress errors
  # calculate baseline with "baseline" package
  bc_pea_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                        method = "peakDetection", 
                        left=30,
                        right=30,
                        lwin=10, 
                        rwin=10)
  
  # extract data from baseline object
  Data$Spectra_pea_orig_rw <- as.data.frame(getSpectra(bc_pea_rw))
  colnames(bc_pea_rw@baseline) <- colnames(bc_pea_rw@spectra)
  Data$Spectra_pea_bl_rw <- as.data.frame(getBaseline(bc_pea_rw))
  Data$Spectra_pea_rw <- as.data.frame(getCorrected(bc_pea_rw))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced spectra
  
  Spectra_orig_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_pea_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_rw) <- groups
  Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw
  
  ### calculate median of reduced spectra
  
  Spectra_orig_median_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_pea_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_rw) <- groups
  Data$Spectra_orig_median_rw <- Spectra_orig_median_rw
  
  ### calculate mean of baseline
  
  Spectra_pea_mean_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_pea_mean_bl_rw <- rbind(Spectra_pea_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_pea_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_pea_mean_bl_rw) <- groups
  Data$Spectra_pea_mean_bl_rw <- Spectra_pea_mean_bl_rw
  
  
  ### calculate median of baseline
  
  Spectra_pea_median_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_pea_median_bl_rw <- rbind(Spectra_pea_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_pea_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_pea_median_bl_rw) <- groups
  Data$Spectra_pea_median_bl_rw <- Spectra_pea_median_bl_rw
  
  
  
  ### combine mean of reduced spectra and baseline
  Data$Spectra_orig_and_bl_mean_rw <- data.frame()
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_pea_mean_bl_rw)
  
  # plot mean of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Mean of original spectra and calculated baseline with Simultaneous Peak Detection and Baseline Correction, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced spectra and baseline
  Data$Spectra_orig_and_bl_med_rw <- data.frame()
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_pea_median_bl_rw)
  
  # plot median of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_med_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Median of original spectra and calculated baseline with Simultaneous Peak Detection and Baseline Correction, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced spectra
  
  Spectra_pea_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_pea_mean_rw <- rbind(Spectra_pea_mean_rw, t(colMeans(as.matrix(Data$Spectra_pea_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_pea_mean_rw) <- groups
  Data$Spectra_pea_mean_rw <- Spectra_pea_mean_rw
  
  # plot mean of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_pea_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (mean) corrected with Simultaneous Peak Detection and Baseline Correction, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_pea_mean_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced spectra
  
  Spectra_pea_median_rw <- data.frame()
  for (name in groups) {
    Spectra_pea_median_rw <- rbind(Spectra_pea_median_rw, t(colMedians(as.matrix(Data$Spectra_pea_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_pea_median_rw) <- groups
  Data$Spectra_pea_median_rw <- Spectra_pea_median_rw
  
  # plot median of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_pea_median_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (median) corrected with Simultaneous Peak Detection and Baseline Correction, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_pea_median_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
  
  ### Simultaneous Peak Detection and Baseline Correction - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
  if (remove_area == T) {
    
    # calculate baseline with "baseline" package
    bc_pea_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                            method = "peakDetection")#, 
    #left=300, 
    #right=300, 
    #lwin=50, 
    #rwin=50)
    
    # extract data from baseline object
    Data$Spectra_pea_orig_area <- as.data.frame(getSpectra(bc_pea_area))
    colnames(bc_pea_area@baseline) <- colnames(bc_pea_area@spectra)
    Data$Spectra_pea_bl_area <- as.data.frame(getBaseline(bc_pea_area))
    Data$Spectra_pea_area <- as.data.frame(getCorrected(bc_pea_area))
    
    par(mfrow = c(num_x,1))
    
    ### calculate mean of reduced/excluded area spectra
    
    Spectra_orig_mean_area <- data.frame()
    for (name in groups) {
      Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_pea_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_mean_area) <- groups
    Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
    
    
    ### calculate median of reduced/excluded area spectra
    
    Spectra_orig_median_area <- data.frame()
    for (name in groups) {
      Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_pea_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_median_area) <- groups
    Data$Spectra_orig_median_area <- Spectra_orig_median_area
    
    
    ### calculate mean of baseline
    
    Spectra_pea_mean_bl_area <- data.frame()
    for (name in groups) {
      Spectra_pea_mean_bl_area <- rbind(Spectra_pea_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_pea_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_pea_mean_bl_area) <- groups
    Data$Spectra_pea_mean_bl_area <- Spectra_pea_mean_bl_area
    
    
    ### calculate median of baseline
    
    Spectra_pea_median_bl_area <- data.frame()
    for (name in groups) {
      Spectra_pea_median_bl_area <- rbind(Spectra_pea_median_bl_area, t(colMedians(as.matrix(Data$Spectra_pea_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_pea_median_bl_area) <- groups
    Data$Spectra_pea_median_bl_area <- Spectra_pea_median_bl_area
    
    
    
    ### combine mean of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_mean_area <- data.frame()
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_pea_mean_bl_area)
    
    # plot mean of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Mean of original spectra and calculated baseline with Simultaneous Peak Detection and Baseline Correction, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### combine median of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_med_area <- data.frame()
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_pea_median_bl_area)
    
    # plot median of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_med_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Median of original spectra and calculated baseline with Simultaneous Peak Detection and Baseline Correction, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    
    ### calculate mean of corrected reduced/excluded area spectra
    
    Spectra_pea_mean_area <- data.frame()
    for (name in groups) {
      Spectra_pea_mean_area <- rbind(Spectra_pea_mean_area, t(colMeans(as.matrix(Data$Spectra_pea_area[Data$Groups==name,]))))
    }
    rownames(Spectra_pea_mean_area) <- groups
    Data$Spectra_pea_mean_area <- Spectra_pea_mean_area
    
    # plot mean of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_pea_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (mean) corrected with Simultaneous Peak Detection and Baseline Correction, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_pea_mean_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### calculate median of corrected reduced/excluded area spectra
    
    Spectra_pea_median_area <- data.frame()
    for (name in groups) {
      Spectra_pea_median_area <- rbind(Spectra_pea_median_area, t(colMedians(as.matrix(Data$Spectra_pea_area[Data$Groups==name,]))))
    }
    rownames(Spectra_pea_median_area) <- groups
    Data$Spectra_pea_median_area <- Spectra_pea_median_area
    
    # plot median of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_pea_median_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (median) corrected with Simultaneous Peak Detection and Baseline Correction, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_pea_median_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    par(mfrow = c(1,1))
  }
  options(warn=0)
  
}

### Robust Baseline Estimation - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN

if (9 %in% pretreatments) {

  # calculate baseline with "baseline" package
  bc_rbe_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                        method = "rfbaseline", 
                        span=0.1)
  
  # extract data from baseline object
  Data$Spectra_rbe_orig_rw <- as.data.frame(getSpectra(bc_rbe_rw))
  colnames(bc_rbe_rw@baseline) <- colnames(bc_rbe_rw@spectra)
  Data$Spectra_rbe_bl_rw <- as.data.frame(getBaseline(bc_rbe_rw))
  Data$Spectra_rbe_rw <- as.data.frame(getCorrected(bc_rbe_rw))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced spectra
  
  Spectra_orig_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_rbe_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_rw) <- groups
  Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw
  
  ### calculate median of reduced spectra
  
  Spectra_orig_median_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_rbe_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_rw) <- groups
  Data$Spectra_orig_median_rw <- Spectra_orig_median_rw
  
  ### calculate mean of baseline
  
  Spectra_rbe_mean_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_rbe_mean_bl_rw <- rbind(Spectra_rbe_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_rbe_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_rbe_mean_bl_rw) <- groups
  Data$Spectra_rbe_mean_bl_rw <- Spectra_rbe_mean_bl_rw
  
  
  ### calculate median of baseline
  
  Spectra_rbe_median_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_rbe_median_bl_rw <- rbind(Spectra_rbe_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_rbe_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_rbe_median_bl_rw) <- groups
  Data$Spectra_rbe_median_bl_rw <- Spectra_rbe_median_bl_rw
  
  ### combine mean of reduced spectra and baseline
  Data$Spectra_orig_and_bl_mean_rw <- data.frame()
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_rbe_mean_bl_rw)
  
  # plot mean of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Mean of original spectra and calculated baseline with Robust Baseline Estimation, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced spectra and baseline
  Data$Spectra_orig_and_bl_med_rw <- data.frame()
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_rbe_median_bl_rw)
  
  # plot median of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_med_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Median of original spectra and calculated baseline with Robust Baseline Estimation, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  
  ### calculate mean of corrected reduced spectra
  
  Spectra_rbe_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_rbe_mean_rw <- rbind(Spectra_rbe_mean_rw, t(colMeans(as.matrix(Data$Spectra_rbe_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_rbe_mean_rw) <- groups
  Data$Spectra_rbe_mean_rw <- Spectra_rbe_mean_rw
  
  # plot mean of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_rbe_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (mean) corrected with Robust Baseline Estimation, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_rbe_mean_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced spectra
  
  Spectra_rbe_median_rw <- data.frame()
  for (name in groups) {
    Spectra_rbe_median_rw <- rbind(Spectra_rbe_median_rw, t(colMedians(as.matrix(Data$Spectra_rbe_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_rbe_median_rw) <- groups
  Data$Spectra_rbe_median_rw <- Spectra_rbe_median_rw
  
  # plot median of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_rbe_median_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (median) corrected with Robust Baseline Estimation, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_rbe_median_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
  
  ### Robust Baseline Estimation - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
  if (remove_area == T) {
    
    # calculate baseline with "baseline" package
    bc_rbe_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                            method = "rfbaseline", 
                            span=NULL, 
                            NoXP=1000)
    
    # extract data from baseline object
    Data$Spectra_rbe_orig_area <- as.data.frame(getSpectra(bc_rbe_area))
    colnames(bc_rbe_area@baseline) <- colnames(bc_rbe_area@spectra)
    Data$Spectra_rbe_bl_area <- as.data.frame(getBaseline(bc_rbe_area))
    Data$Spectra_rbe_area <- as.data.frame(getCorrected(bc_rbe_area))
    
    par(mfrow = c(num_x,1))
    
    ### calculate mean of reduced/excluded area spectra
    
    Spectra_orig_mean_area <- data.frame()
    for (name in groups) {
      Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_rbe_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_mean_area) <- groups
    Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
    
    
    ### calculate median of reduced/excluded area spectra
    
    Spectra_orig_median_area <- data.frame()
    for (name in groups) {
      Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_rbe_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_median_area) <- groups
    Data$Spectra_orig_median_area <- Spectra_orig_median_area
    
    
    ### calculate mean of baseline
    
    Spectra_rbe_mean_bl_area <- data.frame()
    for (name in groups) {
      Spectra_rbe_mean_bl_area <- rbind(Spectra_rbe_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_rbe_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_rbe_mean_bl_area) <- groups
    Data$Spectra_rbe_mean_bl_area <- Spectra_rbe_mean_bl_area
    
    ### calculate median of baseline
    
    Spectra_rbe_median_bl_area <- data.frame()
    for (name in groups) {
      Spectra_rbe_median_bl_area <- rbind(Spectra_rbe_median_bl_area, t(colMedians(as.matrix(Data$Spectra_rbe_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_rbe_median_bl_area) <- groups
    Data$Spectra_rbe_median_bl_area <- Spectra_rbe_median_bl_area
    
    ### combine mean of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_mean_area <- data.frame()
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_rbe_mean_bl_area)
    
    # plot mean of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Mean of original spectra and calculated baseline with Robust Baseline Estimation, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### combine median of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_med_area <- data.frame()
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_rbe_median_bl_area)
    
    # plot median of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_med_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Median of original spectra and calculated baseline with Robust Baseline Estimation, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### calculate mean of corrected reduced/excluded area spectra
    
    Spectra_rbe_mean_area <- data.frame()
    for (name in groups) {
      Spectra_rbe_mean_area <- rbind(Spectra_rbe_mean_area, t(colMeans(as.matrix(Data$Spectra_rbe_area[Data$Groups==name,]))))
    }
    rownames(Spectra_rbe_mean_area) <- groups
    Data$Spectra_rbe_mean_area <- Spectra_rbe_mean_area
    
    # plot mean of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_rbe_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (mean) corrected with Robust Baseline Estimation, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_rbe_mean_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### calculate median of corrected reduced/excluded area spectra
    
    Spectra_rbe_median_area <- data.frame()
    for (name in groups) {
      Spectra_rbe_median_area <- rbind(Spectra_rbe_median_area, t(colMedians(as.matrix(Data$Spectra_rbe_area[Data$Groups==name,]))))
    }
    rownames(Spectra_rbe_median_area) <- groups
    Data$Spectra_rbe_median_area <- Spectra_rbe_median_area
    
    # plot median of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_rbe_median_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (median) corrected with Robust Baseline Estimation, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_rbe_median_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    par(mfrow = c(1,1))
  }
    
}


### Rolling ball - REDUCED WAVENUMBER SPECTRA MEAN AND MEDIAN

if (10 %in% pretreatments) {
  
  # calculate baseline with "baseline" package
  bc_rol_rw <- baseline(spectra = as.matrix(Data$Spectra_min_max), 
                        method = "rollingBall", 
                        wm=40, 
                        ws=40)
  
  # extract data from baseline object
  Data$Spectra_rol_orig_rw <- as.data.frame(getSpectra(bc_rol_rw))
  colnames(bc_rol_rw@baseline) <- colnames(bc_rol_rw@spectra)
  Data$Spectra_rol_bl_rw <- as.data.frame(getBaseline(bc_rol_rw))
  Data$Spectra_rol_rw <- as.data.frame(getCorrected(bc_rol_rw))
  
  par(mfrow = c(num_x,1))
  
  ### calculate mean of reduced spectra
  
  Spectra_orig_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_mean_rw <- rbind(Spectra_orig_mean_rw, t(colMeans(as.matrix(Data$Spectra_rol_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_mean_rw) <- groups
  Data$Spectra_orig_mean_rw <- Spectra_orig_mean_rw
  
  ### calculate median of reduced spectra
  
  Spectra_orig_median_rw <- data.frame()
  for (name in groups) {
    Spectra_orig_median_rw <- rbind(Spectra_orig_median_rw, t(colMedians(as.matrix(Data$Spectra_rol_orig_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_orig_median_rw) <- groups
  Data$Spectra_orig_median_rw <- Spectra_orig_median_rw
  
  ### calculate mean of baseline
  
  Spectra_rol_mean_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_rol_mean_bl_rw <- rbind(Spectra_rol_mean_bl_rw, t(colMeans(as.matrix(Data$Spectra_rol_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_rol_mean_bl_rw) <- groups
  Data$Spectra_rol_mean_bl_rw <- Spectra_rol_mean_bl_rw
  
  
  ### calculate median of baseline
  
  Spectra_rol_median_bl_rw <- data.frame()
  for (name in groups) {
    Spectra_rol_median_bl_rw <- rbind(Spectra_rol_median_bl_rw, t(colMedians(as.matrix(Data$Spectra_rol_bl_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_rol_median_bl_rw) <- groups
  Data$Spectra_rol_median_bl_rw <- Spectra_rol_median_bl_rw
  
  ### combine mean of reduced spectra and baseline
  Data$Spectra_orig_and_bl_mean_rw <- data.frame()
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_orig_mean_rw)
  Data$Spectra_orig_and_bl_mean_rw <- rbind(Data$Spectra_orig_and_bl_mean_rw,Data$Spectra_rol_mean_bl_rw)
  
  # plot mean of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Mean of original spectra and calculated baseline with Rolling ball, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### combine median of reduced spectra and baseline
  Data$Spectra_orig_and_bl_med_rw <- data.frame()
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_orig_median_rw)
  Data$Spectra_orig_and_bl_med_rw <- rbind(Data$Spectra_orig_and_bl_med_rw,Data$Spectra_rol_median_bl_rw)
  
  # plot median of reduced spectra and baseline
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_orig_and_bl_med_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Median of original spectra and calculated baseline with Rolling ball, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(c(groups,groups)))
  
  legend(Legend_pos,
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate mean of corrected reduced spectra
  
  Spectra_rol_mean_rw <- data.frame()
  for (name in groups) {
    Spectra_rol_mean_rw <- rbind(Spectra_rol_mean_rw, t(colMeans(as.matrix(Data$Spectra_rol_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_rol_mean_rw) <- groups
  Data$Spectra_rol_mean_rw <- Spectra_rol_mean_rw
  
  # plot mean of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_rol_mean_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (mean) corrected with Rolling ball, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_rol_mean_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  ### calculate median of corrected reduced spectra
  
  Spectra_rol_median_rw <- data.frame()
  for (name in groups) {
    Spectra_rol_median_rw <- rbind(Spectra_rol_median_rw, t(colMedians(as.matrix(Data$Spectra_rol_rw[Data$Groups==name,]))))
  }
  rownames(Spectra_rol_median_rw) <- groups
  Data$Spectra_rol_median_rw <- Spectra_rol_median_rw
  
  # plot median of corrected reduced spectra
  plot.spectra(Liste = Data,
               Wellenzahl = Data$Wavenumber_min_max,
               Spektren = "Spectra_rol_median_rw",
               Bereich = c(1:length(Data$Wavenumber_min_max)),
               area = c(max(Data$Wavenumber_min_max),min(Data$Wavenumber_min_max)),
               main = paste("Spectra (median) corrected with Rolling ball, 
                          reduced range from",max_range,"to",min_range),
               Code = as.factor(rownames(Data$Spectra_rol_median_rw)))
  
  legend(Legend_pos, 
         legend = legend, 
         pch = 16, 
         col = unique(Data$Groups), 
         inset = 0.05, 
         bty = "n")
  
  par(mfrow = c(1,1))
  
  # Rolling ball - REDUCED WAVENUMBER/EXCLUDED AREA SPECTRA MEAN AND MEDIAN
  if (remove_area == T) {
    
    # calculate baseline with "baseline" package
    bc_rol_area <- baseline(spectra = as.matrix(Data$Spectra_area), 
                            method = "rollingBall", 
                            wm=100, 
                            ws=100)
    
    # extract data from baseline object
    Data$Spectra_rol_orig_area <- as.data.frame(getSpectra(bc_rol_area))
    colnames(bc_rol_area@baseline) <- colnames(bc_rol_area@spectra)
    Data$Spectra_rol_bl_area <- as.data.frame(getBaseline(bc_rol_area))
    Data$Spectra_rol_area <- as.data.frame(getCorrected(bc_rol_area))
    
    par(mfrow = c(num_x,1))
    
    ### calculate mean of reduced/excluded area spectra
    
    Spectra_orig_mean_area <- data.frame()
    for (name in groups) {
      Spectra_orig_mean_area <- rbind(Spectra_orig_mean_area, t(colMeans(as.matrix(Data$Spectra_rol_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_mean_area) <- groups
    Data$Spectra_orig_mean_area <- Spectra_orig_mean_area
    
    ### calculate median of reduced/excluded area spectra
    
    Spectra_orig_median_area <- data.frame()
    for (name in groups) {
      Spectra_orig_median_area <- rbind(Spectra_orig_median_area, t(colMedians(as.matrix(Data$Spectra_rol_orig_area[Data$Groups==name,]))))
    }
    rownames(Spectra_orig_median_area) <- groups
    Data$Spectra_orig_median_area <- Spectra_orig_median_area
    
    ### calculate mean of baseline
    
    Spectra_rol_mean_bl_area <- data.frame()
    for (name in groups) {
      Spectra_rol_mean_bl_area <- rbind(Spectra_rol_mean_bl_area, t(colMeans(as.matrix(Data$Spectra_rol_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_rol_mean_bl_area) <- groups
    Data$Spectra_rol_mean_bl_area <- Spectra_rol_mean_bl_area
    
    ### calculate median of baseline
    
    Spectra_rol_median_bl_area <- data.frame()
    for (name in groups) {
      Spectra_rol_median_bl_area <- rbind(Spectra_rol_median_bl_area, t(colMedians(as.matrix(Data$Spectra_rol_bl_area[Data$Groups==name,]))))
    }
    rownames(Spectra_rol_median_bl_area) <- groups
    Data$Spectra_rol_median_bl_area <- Spectra_rol_median_bl_area
    
    ### combine mean of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_mean_area <- data.frame()
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_orig_mean_area)
    Data$Spectra_orig_and_bl_mean_area <- rbind(Data$Spectra_orig_and_bl_mean_area,Data$Spectra_rol_mean_bl_area)
    
    # plot mean of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Mean of original spectra and calculated baseline with Rolling ball, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### combine median of reduced/excluded area spectra and baseline
    Data$Spectra_orig_and_bl_med_area <- data.frame()
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_orig_median_area)
    Data$Spectra_orig_and_bl_med_area <- rbind(Data$Spectra_orig_and_bl_med_area,Data$Spectra_rol_median_bl_area)
    
    # plot median of reduced/excluded area spectra and baseline
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_orig_and_bl_med_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Median of original spectra and calculated baseline with Rolling ball, 
                            reduced range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(c(groups,groups)))
    
    legend(Legend_pos,
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### calculate mean of corrected reduced/excluded area spectra
    
    Spectra_rol_mean_area <- data.frame()
    for (name in groups) {
      Spectra_rol_mean_area <- rbind(Spectra_rol_mean_area, t(colMeans(as.matrix(Data$Spectra_rol_area[Data$Groups==name,]))))
    }
    rownames(Spectra_rol_mean_area) <- groups
    Data$Spectra_rol_mean_area <- Spectra_rol_mean_area
    
    # plot mean of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_rol_mean_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (mean) corrected with Rolling ball, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_rol_mean_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    ### calculate median of corrected reduced/excluded area spectra
    
    Spectra_rol_median_area <- data.frame()
    for (name in groups) {
      Spectra_rol_median_area <- rbind(Spectra_rol_median_area, t(colMedians(as.matrix(Data$Spectra_rol_area[Data$Groups==name,]))))
    }
    rownames(Spectra_rol_median_area) <- groups
    Data$Spectra_rol_median_area <- Spectra_rol_median_area
    
    # plot median of corrected reduced/excluded area spectra
    plot.spectra(Liste = Data,
                 Wellenzahl = Data$Wavenumber_area,
                 Spektren = "Spectra_rol_median_area",
                 Bereich = c(1:length(Data$Wavenumber_area)),
                 area = c(max(Data$Wavenumber_area),min(Data$Wavenumber_area)),
                 main = paste("Spectra (median) corrected with Rolling ball, 
                            reduced wavenumber range from",max_range,"to",min_range,"and excluded range from",area[2],"to",area[1]),
                 Code = as.factor(rownames(Data$Spectra_rol_median_area)))
    
    legend(Legend_pos, 
           legend = legend, 
           pch = 16, 
           col = unique(Data$Groups), 
           inset = 0.05, 
           bty = "n")
    
    par(mfrow = c(1,1))
  }
  
}



# PCA - reduced wavenumber 
par(mfrow = c(1,1))

Type = Type_rw
Sample <- c(1:length(Data$Wavenumber_min_max))

for (Stats in pretreatments) {
  print("---------")
  print(Names[Stats])
  
  pca.data <- Data[[Type[Stats]]][!is.na(Data$Groups),Sample]
  
  if (vector_normalize == T) {
    
    # Normalization
    norm.factors <- apply(pca.data, 1, function(x){sqrt(sum(x^2))})
    pca.data <- sweep(pca.data, 1, norm.factors, "/")
    
  }

  
  PCA <- prcomp(pca.data,
                center = TRUE,
                rank. = PC_PCA, 
                scale. = FALSE )
  
  # remove "NAs" for correct colours in the plot
  PCA$Groups <- Data$Groups[!is.na(Data$Groups)]
  
  print(summary(PCA))
  
  Data[[paste("Var_", Type[Stats])]] <- PCA$sdev^2/sum(PCA$sdev^2)
  
  plot(x = 1:PC, 
       y = cumsum(Data[[paste("Var_", Type[Stats])]][1:PC]*100),
       type = "b", 
       main = paste("Method: ",Names[Stats], sep=""),
       xlab = "Number of Principal Components",
       ylab = "Explained Variance [%]",
       col = "dimgrey", 
       pch = 21, 
       bg = "darkgrey",
       font = 2, 
       font.lab = 2,
       lab = c(10,20,20),
       xaxt = "n")
  axis(1,xaxp=c(1,PC,PC-1))
  
  grid(lwd = 0.8)
  
  pairs(PCA$x[,1:PCs],
        col = Data$Groups,
        main = paste("Method: ",Names[Stats], sep=""),
        pch = 19)
  
  par(xpd = T) 
  
  legend("bottomright", 
         legend = legend, 
         pch = 19, 
         col = unique(Data$Groups), 
         inset = 0.05,
         cex = 0.7,
         horiz = T, 
         bty = "n")
  
  par(xpd = F)
  
  plot(x = PCA$x[,x_PC], 
       y = PCA$x[,y_PC],
       xlab = paste("Scores of","PC" ,x_PC, "[",round(Data[[paste("Var_", Type[Stats])]][x_PC]*100,1),"%]", sep = " "),
       ylab = paste("Scores of","PC" ,y_PC, "[",round(Data[[paste("Var_", Type[Stats])]][y_PC]*100,1),"%]", sep = " "),
       main = paste("Method: ",Names[Stats], sep=""),
       col = Data$Groups[!is.na(Data$Groups)], 
       pch = 19, 
       font.lab = 2)
  
  grid(lwd = 0.8)
  
  abline(h = 0, v = 0)
  
  legend("bottomright", 
         legend = legend, 
         pch = 20, 
         col = unique(Data$Groups[!is.na(Data$Groups)]), 
         inset = 0.01,
         bty = "n")
  
  # remove outlier
  
  Data[paste(Type[Stats],"_outl",sep="")] <- Data[paste(Type[Stats])]
  Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]] <- Data$Groups
  
  if(select_outlier == T){
    print(paste(Stats, "Select outlier: ESC to exit"))
    outlier = identify( x = PCA$x[,x_PC], y = PCA$x[,y_PC])
    print(outlier)
    print(Data$Files[outlier])
    Data[[paste(Type[Stats],"_outl",sep="")]][outlier,] = NA
    Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]][outlier] = NA
    #print(Data$Groups)
    #print(Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]])
  }
  
  if (plot_loadings == T) {
    
    DataLoadings <- data.frame(WaveNum=Data$Wavenumber_min_max)
    
    for (loadingNR in (1:loadings)) {
      #print(paste("PC",loadingNR,sep=""))
      DataLoadings <- cbind(DataLoadings,PCA$rotation[,loadingNR])
      colnames(DataLoadings)[(loadingNR+1)] <- paste("PC",loadingNR,sep="")
    }
    
    DataLoadingsRev <- DataLoadings[dim(DataLoadings)[1]:1,]
    
    matplot(x = DataLoadingsRev[1],
            y = DataLoadingsRev[-1],
            type = "l", 
            main = paste("Loadings: ",Names[Stats], sep=""), 
            ylab = "",
            xlab = "Wavenumber [1/cm]",
            xlim = c(max(DataLoadings$WaveNum),min(DataLoadings$WaveNum)),
            xaxs = "i")
    abline(0,0, col="black", lwd=1.5)
    legend("topleft", 
           legend = colnames(DataLoadings)[2:(loadings+1)], 
           pch = 20, 
           col = c(1:loadings), 
           inset = 0.01,
           bty = "n")
    grid(lwd = 0.8)
    
  }
}

### PCA - reduced wavenumber - without outlier
if (select_outlier == T) {
  for (Stats in pretreatments) {
    print(Names[Stats])
    
    pca.data <- na.omit(Data[[paste(Type[Stats],"_outl", sep="")]][!is.na(Data$Groups),Sample])
    
    if (vector_normalize == T) {
      
      # Normalization
      norm.factors <- apply(pca.data, 1, function(x){sqrt(sum(x^2))})
      pca.data <- sweep(pca.data, 1, norm.factors, "/")
      
    }
    
    PCA <- prcomp(pca.data,
                  center = TRUE,
                  rank. = PC_PCA, 
                  scale. = FALSE )
    
    # remove "NAs" for correct colours in the plot
    PCA$Groups <- Data$Groups[!is.na(Data$Groups)]
    #print(Data$Groups)
    #print(Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]])])
    
    print(summary(PCA))
    
    Data[[paste("Var_", Type[Stats],"_outl", sep="")]] <- PCA$sdev^2/sum(PCA$sdev^2)
    
    plot(x = 1:PC, 
         y = cumsum(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][1:PC]*100),
         type = "b", 
         main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
         xlab = "Number of Principal Components",
         ylab = "Explained Variance [%]",
         col = "dimgrey", 
         pch = 21, 
         bg = "darkgrey",
         font = 2, 
         font.lab = 2,
         lab = c(10,20,20),
         xaxt = "n")
    axis(1,xaxp=c(1,PC,PC-1))
    
    grid(lwd = 0.8)
    
    pairs(PCA$x[,1:PCs],
          #col = Data$Groups,
          col = Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]])],
          main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
          pch = 19)
    
    par(xpd = T) 
    
    legend("bottomright", 
           legend = legend, 
           pch = 19, 
           col = unique(Data$Groups), 
           inset = 0.05,
           cex = 0.7,
           horiz = T, 
           bty = "n")
    
    par(xpd = F)
    
    plot(x = PCA$x[,x_PC], 
         y = PCA$x[,y_PC],
         xlab = paste("Scores of","PC" ,x_PC, "[",round(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][x_PC]*100,1),"%]", sep = " "),
         ylab = paste("Scores of","PC" ,y_PC, "[",round(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][y_PC]*100,1),"%]", sep = " "),
         main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
         #col = Data$Groups[!is.na(Data$Groups)],
         col = Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_red_Groups",sep="")]])],
         pch = 19, 
         font.lab = 2)
    
    grid(lwd = 0.8)
    
    abline(h = 0, v = 0)
    
    legend("bottomleft", 
           legend = legend, 
           pch = 20, 
           col = unique(Data$Groups[!is.na(Data$Groups)]), 
           inset = 0.01,
           bty = "n")
    
    if (plot_loadings == T) {
      #plot.load(main=paste("Method: ",Names[Stats], " (Outl)", sep=""))
      
      DataLoadings <- data.frame(WaveNum=Data$Wavenumber_min_max)
      for (loadingNR in (1:loadings)) {
        #print(paste("PC",loadingNR,sep=""))
        DataLoadings <- cbind(DataLoadings,PCA$rotation[,loadingNR])
        colnames(DataLoadings)[(loadingNR+1)] <- paste("PC",loadingNR,sep="")
      }
      
      DataLoadingsRev <- DataLoadings[dim(DataLoadings)[1]:1,]
      
      matplot(x = DataLoadingsRev[1],
              y = DataLoadingsRev[-1],
              type = "l", 
              main=paste("Loadings: ",Names[Stats], " (Outl)", sep=""),
              xlab = "Wavenumber [1/cm]",
              ylab = "", 
              #xaxt ="n",
              xaxs = "i",
              xlim = c(max(DataLoadings$WaveNum),min(DataLoadings$WaveNum)))
      
      abline(0,0, col="black", lwd=1.5)
      legend("topleft", 
             legend = colnames(DataLoadings)[2:(loadings+1)], 
             pch = 20, 
             col = c(1:loadings), 
             inset = 0.01,
             bty = "n")
      grid(lwd = 0.8)
      
    }
  }
}


### PCA - reduced wavenumber/area 
if (remove_area == T) {

  Type = Type_a
  Sample <- c(1:length(Data$Wavenumber_area))
  
  for (Stats in pretreatments) {
    print(Names[Stats])
    
    pca.data <- Data[[Type[Stats]]][!is.na(Data$Groups),Sample]
    
    if (vector_normalize == T) {
      
      # Normalization
      norm.factors <- apply(pca.data, 1, function(x){sqrt(sum(x^2))})
      pca.data <- sweep(pca.data, 1, norm.factors, "/")
      
    }
    
    PCA <- prcomp(pca.data,
                  center = TRUE,
                  rank. = PC_PCA, 
                  scale. = FALSE )
    
    # remove "NAs" for correct colours in the plot
    PCA$Groups <- Data$Groups[!is.na(Data$Groups)]
    
    print(summary(PCA))
    
    Data[[paste("Var_", Type[Stats])]] <- PCA$sdev^2/sum(PCA$sdev^2)
    
    plot(x = 1:PC, 
         y = cumsum(Data[[paste("Var_", Type[Stats])]][1:PC]*100),
         type = "b", 
         main = paste("Method: ",Names[Stats], sep=""),
         xlab = "Number of Principal Components",
         ylab = "Explained Variance [%]",
         col = "dimgrey", 
         pch = 21, 
         bg = "darkgrey",
         font = 2, 
         font.lab = 2,
         lab = c(10,20,20),
         xaxt = "n")
    axis(1,xaxp=c(1,PC,PC-1))
    
    grid(lwd = 0.8)
    
    pairs(PCA$x[,1:PCs],
          col = Data$Groups,
          main = paste("Method: ",Names[Stats], sep=""),
          pch = 19)
    
    par(xpd = T) 
    
    legend("bottomright", 
           legend = legend, 
           pch = 19, 
           col = unique(Data$Groups), 
           inset = 0.05,
           cex = 0.7,
           horiz = T, 
           bty = "n")
    
    par(xpd = F)
    
    plot(x = PCA$x[,x_PC], 
         y = PCA$x[,y_PC],
         xlab = paste("Scores of","PC" ,x_PC, "[",round(Data[[paste("Var_", Type[Stats])]][x_PC]*100,1),"%]", sep = " "),
         ylab = paste("Scores of","PC" ,y_PC, "[",round(Data[[paste("Var_", Type[Stats])]][y_PC]*100,1),"%]", sep = " "),
         main = paste("Method: ",Names[Stats], sep=""),
         col = Data$Groups[!is.na(Data$Groups)], 
         pch = 19, 
         font.lab = 2)
    
    grid(lwd = 0.8)
    
    abline(h = 0, v = 0)
    
    legend("bottomleft", 
           legend = legend, 
           pch = 20, 
           col = unique(Data$Groups[!is.na(Data$Groups)]), 
           inset = 0.01,
           bty = "n")
    
    # remove outlier
    Data[paste(Type[Stats],"_outl",sep="")] <- Data[paste(Type[Stats])]
    Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]] <- Data$Groups
    
    if(select_outlier == T){
      print(paste(Stats, "Select outlier: ESC to exit"))
      outlier = identify( x = PCA$x[,x_PC], y = PCA$x[,y_PC])
      print(outlier)
      print(Data$Files[outlier])
      Data[[paste(Type[Stats],"_outl",sep="")]][outlier,] = NA
      Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]][outlier] = NA
      #print(Data$Groups)
      #print(Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]])
      
    }
    
    if (plot_loadings == T) {
      
      DataLoadings <- data.frame(WaveNum=Data$Wavenumber_area)
      
      for (loadingNR in (1:loadings)) {
        #print(paste("PC",loadingNR,sep=""))
        DataLoadings <- cbind(DataLoadings,PCA$rotation[,loadingNR])
        colnames(DataLoadings)[(loadingNR+1)] <- paste("PC",loadingNR,sep="")
      }
      
      DataLoadingsRev <- DataLoadings[dim(DataLoadings)[1]:1,]
      
      matplot(x = DataLoadingsRev[1],
              y = DataLoadingsRev[-1],
              type = "l", 
              main = paste("Loadings: ",Names[Stats], sep=""), 
              ylab = "",
              xlab = "Wavenumber [1/cm]",
              xlim = c(max(DataLoadings$WaveNum),min(DataLoadings$WaveNum)),
              xaxs = "i")
      abline(0,0, col="black", lwd=1.5)
      legend("topleft", 
             legend = colnames(DataLoadings)[2:(loadings+1)], 
             pch = 20, 
             col = c(1:loadings), 
             inset = 0.01,
             bty = "n")
      grid(lwd = 0.8)
      
    }
  }
}

### PCA - reduced wavenumber/area - without outlier
if (remove_area == T) {
  
  if (select_outlier == T) {
    for (Stats in pretreatments) {
      print(Names[Stats])
      
      pca.data <- na.omit(Data[[paste(Type[Stats],"_outl", sep="")]][!is.na(Data$Groups),Sample])
      
      if (vector_normalize == T) {
       
        # Normalization
        norm.factors <- apply(pca.data, 1, function(x){sqrt(sum(x^2))})
        pca.data <- sweep(pca.data, 1, norm.factors, "/")
         
      }
      
      PCA <- prcomp(pca.data,
                    center = TRUE,
                    rank. = PC_PCA, 
                    scale. = FALSE )
      
      # remove "NAs" for correct colours in the plot
      PCA$Groups <- Data$Groups[!is.na(Data$Groups)]
      #print(Data$Groups)
      #print(Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]])])
      
      print(summary(PCA))
      
      Data[[paste("Var_", Type[Stats],"_outl", sep="")]] <- PCA$sdev^2/sum(PCA$sdev^2)
      
      plot(x = 1:PC, 
           y = cumsum(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][1:PC]*100),
           type = "b", 
           main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
           xlab = "Number of Principal Components",
           ylab = "Explained Variance [%]",
           col = "dimgrey", 
           pch = 21, 
           bg = "darkgrey",
           font = 2, 
           font.lab = 2,
           lab = c(10,20,20),
           xaxt = "n")
      axis(1,xaxp=c(1,PC,PC-1))
      
      grid(lwd = 0.8)
      
      pairs(PCA$x[,1:PCs],
            #col = Data$Groups,
            col = Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]])],
            main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
            pch = 19)
      
      par(xpd = T) 
      
      legend("bottomright", 
             legend = legend, 
             pch = 19, 
             col = unique(Data$Groups), 
             inset = 0.05,
             cex = 0.7,
             horiz = T, 
             bty = "n")
      
      par(xpd = F)
      
      plot(x = PCA$x[,x_PC], 
           y = PCA$x[,y_PC],
           xlab = paste("Scores of","PC" ,x_PC, "[",round(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][x_PC]*100,1),"%]", sep = " "),
           ylab = paste("Scores of","PC" ,y_PC, "[",round(Data[[paste("Var_", Type[Stats],"_outl", sep="")]][y_PC]*100,1),"%]", sep = " "),
           main = paste("Method: ",Names[Stats], " (Outl)", sep=""),
           #col = Data$Groups[!is.na(Data$Groups)],
           col = Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]][!is.na(Data[[paste(Type[Stats],"_outl_ar_Groups",sep="")]])],
           pch = 19, 
           font.lab = 2)
      
      grid(lwd = 0.8)
      
      abline(h = 0, v = 0)
      
      legend("bottomleft", 
             legend = legend, 
             pch = 20, 
             col = unique(Data$Groups[!is.na(Data$Groups)]), 
             inset = 0.01,
             bty = "n")
      
      if (plot_loadings == T) {
        
        DataLoadings <- data.frame(WaveNum=Data$Wavenumber_area)
        for (loadingNR in (1:loadings)) {
          #print(paste("PC",loadingNR,sep=""))
          DataLoadings <- cbind(DataLoadings,PCA$rotation[,loadingNR])
          colnames(DataLoadings)[(loadingNR+1)] <- paste("PC",loadingNR,sep="")
        }
        
        DataLoadingsRev <- DataLoadings[dim(DataLoadings)[1]:1,]
        
        matplot(x = DataLoadingsRev[1],
                y = DataLoadingsRev[-1],
                type = "l", 
                main=paste("Loadings: ",Names[Stats], " (Outl)", sep=""),
                xlab = "Wavenumber [1/cm]",
                ylab = "", 
                #xaxt ="n",
                xaxs = "i",
                xlim = c(max(DataLoadings$WaveNum),min(DataLoadings$WaveNum)))
        
        abline(0,0, col="black", lwd=1.5)
        legend("topleft", 
               legend = colnames(DataLoadings)[2:(loadings+1)], 
               pch = 20, 
               col = c(1:loadings), 
               inset = 0.01,
               bty = "n")
        grid(lwd = 0.8)
        
      }
    }
  }
}

rm(PCA, pca.data)


if (perform_LDA == T) {
  
  if (remove_area == TRUE) {
    Type = Type_a
  } else {
    Type = Type_rw
  }
  
  
  for (Stats in pretreatments) {
    
    print("---------")
    print(Names[Stats])
    
    #extract data
    if (select_outlier == TRUE) {
      data.lda <- na.omit(Data[[paste(Type[Stats],"_outl", sep="")]][!is.na(Data$Groups),Sample])
      groups.lda <- na.omit(Data[[paste(Type[Stats],"_outl_red_Groups", sep="")]])
    } else {
      data.lda <- Data[[Type[Stats]]][!is.na(Data$Groups),Sample]
      groups.lda <- Data$Groups
    }
    colnames(data.lda) <- Data$Wavenumber_min_max
    
    if (vector_normalize == T) {
      # normalization
      norm.factors <- apply(data.lda, 1, function(x){sqrt(sum(x^2))})
      data.lda <- sweep(data.lda, 1, norm.factors, "/")
      rm(norm.factors)
      
    }
    
    i <- 0
    j <- 0
    
    max.LDs <- length(levels(groups.lda))-1
    
    # Prepare matrices and arrays for results
    CV.results <- array(NA, c(max.PCs, 3, segments.out*repetitions))
    colnames(CV.results) <- c("PCs", "mean_error", "sd_error")
    opt.PCs <- rep(NA, repetitions*segments.out)
    PCA.scores <- array(NA, c(nrow(data.lda), max.PCs, segments.out*repetitions))
    PCA.loadings <- array(NA, c(ncol(data.lda), max.PCs, segments.out*repetitions))
    
    LDA.scores <- matrix(NA, nrow=nrow(data.lda), ncol=repetitions)
    LDA.loadings <- matrix(NA, nrow=ncol(data.lda), ncol=repetitions) 
    predict.probs <- matrix(NA, nrow=nrow(data.lda), ncol=repetitions)
    
    cm.results <- array(NA, c(length(levels(groups.lda)),length(levels(groups.lda)), 
                              segments.out*repetitions), dimnames = list(levels(groups.lda),levels(groups.lda))) 
    predict.results <- array(NA, c(ceiling(length(groups.lda)/(segments.out-1)), 2, segments.out*repetitions))
    colnames(predict.results) <- c("reference", "prediction")
    
    roc.results <- array(NA, c(ceiling(length(groups.lda)/(segments.out-1)), 2, segments.out*repetitions))
    colnames(roc.results) <- c("FPR", "TPR")
    
    auc.results <- rep(NA, repetitions*segments.out)
    
    for (r in c(1:repetitions)) { # START Cross Validation
      tfolds <- createFolds(groups.lda, k=segments.out, list=T) # split data into segments
      
      for (t in c(1:segments.out)) {
      
        i <- i + 1
        # 1 Segment for testing, rest for calibration
        test.set <- sort(unlist(tfolds[t], use.names=F))
        calib.set <- sort(unlist(tfolds[-t], use.names=F))
        
        kfolds <- createFolds(groups.lda[calib.set], k=segments.in, list=T) # split calibration data into new segments
        CV.err <- matrix(NA, nrow=segments.in, ncol=max.PCs) # prepare matrix for misclassification errors
        
        for (k in c(1:segments.in)) {
          
          # 1 Segment far validation, rest for training
          valid.set <- unlist(kfolds[k], use.names=F)
          train.set <- unlist(kfolds[-k], use.names=F)
          
          # PCA of training set
          PCA <- prcomp(data.lda[calib.set,][train.set,], center=T, scale=F, rank=max.PCs)
          scores.train <- PCA$x
          scores.valid <- predict(PCA, data.lda[calib.set,][valid.set,])
          
          for (p in c(1:max.PCs)) { # Different amount of PCs
            
            train.df <- data.frame(scores.train[, 1:p])
            colnames(train.df) <- c(1:p)
            valid.df <- data.frame(scores.valid[, 1:p])
            colnames(valid.df) <- c(1:p)
            
            #LDA of training set
            model.LDA <- lda(train.df, grouping=groups.lda[calib.set][train.set], CV=F)
            valid.predict <- predict(model.LDA, valid.df)
            
            valid.groups <- as.matrix(table(valid.predict$class, groups.lda[calib.set][valid.set]))
            
            # Calculate misclassification rate
            delta <- row(valid.groups) - col(valid.groups)
            CV.err[k,p] <- sum(valid.groups[delta < 0 | delta > 0]) / sum(valid.groups)
            
          } # PCs loop
          
        } # inner CV loop
          
        # Add error rates to results array
        CV.results[,1,i] <- c(1:max.PCs)
        CV.results[,2,i] <- apply(CV.err, 2, mean)
        CV.results[,3,i] <- apply(CV.err, 2, sd)
        
        CV.res.min <- CV.results[CV.results[,2,i] == min(CV.results[,2,i]),,i]
        if (is.vector(CV.res.min)) {CV.res.min <- t(as.matrix(CV.res.min))}
        
        parsimony <- 0 # How many stdevs of error rate the "optimal" number of 
                       # components can be above the minimum one
        threshold <- CV.res.min[1,2] + parsimony * CV.res.min[1,3]/sqrt(segments.in)
        
        CV.res.opt <- CV.results[(CV.results[,2,i] <= threshold),,i]
        if (is.vector(CV.res.opt)) {CV.res.opt <- t(as.matrix(CV.res.opt))}
        opt.PCs[i] <- CV.res.opt[1,1]
        
        # PCA of the whole calibration set 
        PCA <- prcomp(data.lda[calib.set,], center=T, scale=F, rank=max.PCs)
        scores.calib <- PCA$x
        PCA.scores[,,i][1:nrow(PCA$x),] <- PCA$x
        PCA.loadings[,,i] <- PCA$rotation
        
        # Save scores for later
        scores.calib.df <- data.frame(groups.lda[calib.set], PCA$x)
        scores.test <- predict(PCA, data.lda[test.set,])
        PCA.scores[,,i][test.set,] <- scores.test
        
        # LDA with optimal number of PCs
        model.LDA <- lda(data.frame(scores.calib[,1:opt.PCs[i]]), grouping=groups.lda[calib.set], CV=F)
        options(warn=-1)
        test.predict <- predict(model.LDA, data.frame(scores.test[,1:opt.PCs[i]]))
        options(warn=0)
        LDA.scores[test.set, r] <- test.predict$x
        LDA.loadings[, r] <- PCA$rotation[,1:opt.PCs[i]] %*% model.LDA$scaling
        
        predict.probs[test.set, r] <- test.predict$posterior[,colnames(test.predict$posterior)==pos.class]
        
        cm <- confusionMatrix(data=as.factor(test.predict$class), reference=groups.lda[test.set], positive=pos.class)
        cm.results[,,i] <- cm$table
        
        predict.results[,,i][1:length(groups.lda[test.set]),1] <- groups.lda[test.set] # store reference labels of TEST
        predict.results[,,i][1:length(test.predict$class),2] <- as.factor(test.predict$class) # store predictions for TEST
        
        # Receiver Operator Characteristic Curve
        roc.prob <- prediction(predict.probs[test.set, r], groups.lda[test.set])
        roc.perf <- performance(roc.prob, "tpr","fpr")
        auc.results[i] <- performance(roc.prob, "auc")@y.values[1]
        
        roc.results[,,i][1:length(roc.perf@x.values[[1]]),1] <-roc.perf@x.values[[1]]
        roc.results[,,i][1:length(roc.perf@y.values[[1]]),2] <-roc.perf@y.values[[1]]
        
        print (paste("repetition =",r,"/", repetitions, "  fold =",t, "/", segments.out, "  model",i))
        
      } # outer CV loop
      
    } # repetitions loop
    
    rm(calib.set, cm, CV.err, CV.res.min, delta, 
       kfolds, model.LDA, PCA, roc.perf, roc.prob, scores.calib, scores.calib.df, scores.test, 
       scores.train, scores.valid, test.predict, tfolds, train.df, valid.df, 
       valid.predict)
    
    # Postprocessing
    
    # Add together confusion matrices of each outer loop
    cm.results.merged <- array(NA, c(2,2,repetitions))
    
    j <- 1
    for (i in seq(1, (segments.out*repetitions), by=segments.out)) {
      cm.results.merged[,,j] <- apply(cm.results[,,c(i:(i+(segments.out-1)))], c(1,2), sum)
      j <- j + 1
    }
    cm.results <- cm.results.merged
    
    #########################
    # Average ROC curves of each outer loop
    roc.results.avg <- array(NA, c(nrow(roc.results), 2, repetitions))
    
    j <- 1
    for (i in seq(1, (segments.out*repetitions), by=segments.out)) {
      roc.results.avg[,,j] <- apply(roc.results[,,c(i:(i+(segments.out-1)))], c(1,2), mean)
      j <- j + 1
    }
    
    roc.results <- roc.results.avg
    colnames(roc.results) <- c("FPR", "TPR")
    
    rm(cm.results.merged, roc.results.avg)
    
    ############
    # Plotting #
    ############
    
    # Mean Error rate by number of PCs
    CV.results.df <- as.data.frame(CV.results)
    CV.summary.df <- data.frame(PCs=c(1:max.PCs))

    CV.summary.df$mean_error <- apply(CV.results[,2,], 1, mean)
    CV.summary.df$sdlow <- apply(CV.results[,2,], 1, function(x) {ifelse(mean(x)-sd(x)<0, 0, mean(x)-sd(x))})                        
    CV.summary.df$sdhigh <- apply(CV.results[,2,], 1, function(x) {mean(x)+sd(x)})
    
    p <- ggplot(data=CV.results.df)
    
    for (i in c(0:(segments.out*repetitions-1))) {
      p <- p + geom_line(aes_string(x=names(CV.results.df)[i*3+1], y=names(CV.results.df)[i*3+2], col=shQuote("lg")))
    }
    
    p <- p + geom_line(data=CV.summary.df, mapping=aes(x=PCs, y=mean_error, col="black"), size=1) +
             geom_line(data=CV.summary.df, mapping=aes(x=PCs, y=sdlow, col="black"), size=1, linetype="dashed") +
             geom_line(data=CV.summary.df, mapping=aes(x=PCs, y=sdhigh, col="black"), size=1, linetype="dashed") + 
             scale_color_manual("Legend", values=c(lg="#d3d3d3", black="black")) + 
             theme_bw() + theme(legend.position = "none") + 
             labs(x="Principal Components", y="Mean Error Rate") +
             scale_x_continuous(breaks=c(1:max.PCs)) +
             coord_cartesian(expand=F)
    
    print(p)
    
    ############################################
    
    # Frequency of each number of PCs being determined as optimal
    p <- ggplot(mapping=aes(opt.PCs)) + geom_bar(aes(y= ..count..*100/sum(..count..))) + theme_bw() + 
         labs(x="Optimal number of PCs", y="Frequency [%]") +
         scale_x_continuous(expand = expansion(mult=0.01), breaks=c(1:max.PCs)) +
         scale_y_continuous(expand = expansion(mult=c(0, 0.05))) +
         coord_cartesian(xlim=c(0.5,max.PCs+0.5))
    
    print(p)
    
    ################################################
    
    # Confusion matrix with histograms
    
    pos.class.name <- legend[groups == pos.class]
    neg.class.name <- legend[groups != pos.class]
    
    cm.class <- factor(rep(c(neg.class.name, neg.class.name, pos.class.name, pos.class.name),
                           times=repetitions),
                       levels=c(neg.class.name, pos.class.name))
    
    cm.pred <- factor(rep(c(neg.class.name, pos.class.name, neg.class.name, pos.class.name),
                          times=repetitions),
                      levels=c(neg.class.name, pos.class.name))
    
    cm.quadrant <- factor(rep(c("TN", "FP", "FN", "TP"), times=repetitions), levels=c("TN", "FP", "FN", "TP"))
    
    cm.df <- data.frame(c(cm.results), cm.class, cm.pred, cm.quadrant)
    colnames(cm.df) <- c("Count", "Group", "Prediction", "Quadrant")
    
    
    x_max <- length(groups.lda) / length(levels(groups.lda))
    y_max <- max(as.data.frame(table(cm.df$Count))$Freq) * 1.2
    
    cm.annotations <- as.data.frame(aggregate(cm.df$Count, by=list(cm.df$Quadrant), median))
    cm.annotations <- as.data.frame(aggregate(cm.df$Count, by=list(cm.df$Group, cm.df$Prediction), median))
    colnames(cm.annotations) <- c("Group", "Prediction", "Median")
    
    cm.annotations$x <- x_max/2
    cm.annotations$y <- y_max/2
    
    cm.labs <- c(neg.class.name, pos.class.name)
    
    cm.xlabs <- paste("Prediction:", cm.labs, sep=" ")
    names(cm.xlabs) <- c(neg.class.name, pos.class.name)
    
    cm.ylabs <- paste("Group:", cm.labs, sep=" ")
    names(cm.ylabs) <- c(neg.class.name, pos.class.name)
    
    p <- ggplot(data=cm.df, 
                aes(x=Count)) +
      
         coord_cartesian(xlim=c(0,x_max), ylim=c(0,y_max)) +
      
         geom_bar(width=1) + 
      
         geom_vline(data=cm.annotations, size=1,
                    mapping=aes(xintercept=Median, col="red")) +
      
         geom_text(data=cm.annotations, 
                   aes(label=Median, 
                       x=x, y=y), 
                   color="red",
                   size=5) +
      
         facet_grid(rows=vars(Group), cols=vars(Prediction), switch="y",
                    labeller=labeller(.rows = cm.ylabs, .cols = cm.xlabs)) +
      
         scale_y_continuous(position="right") +
    
         labs(x="No. of elements", 
              y="Frequency (%)") +
      
         theme_bw() + 
         theme(legend.position = "none",
               strip.text.x = element_text(size=12, face="bold"),
               strip.text.y = element_text(size=12, face="bold"))
    
    print(p)
    
    
    ###############################
    
    # Quality measures of classification
    cm.values <- data.frame(matrix(NA, nrow=repetitions, ncol=6))
    colnames(cm.values) <- c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "AUROC")
    
    TN <- cm.results[1,1,]
    FN <- cm.results[1,2,]
    FP <- cm.results[2,1,]
    TP <- cm.results[2,2,]
    
    cm.values$Accuracy <- (TN+TP) / (TN+FP+FN+TP)
    cm.values$Sensitivity <- TP / (TP+FN)
    cm.values$Specificity <- TN / (TN+FP)
    cm.values$PPV <- TP / (TP+FP)
    cm.values$NPV <- TN / (TN+FN)
    
    cm.values.ci <- data.frame(matrix(NA, nrow=6, ncol=3))
    colnames(cm.values.ci) <- c("Mean", "Low", "High")
    rownames(cm.values.ci) <- c("Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "AUROC")
    
    n_tot <- length(groups.lda)
    n_pos <- sum(groups.lda == pos.class)
    n_neg <- n_tot - n_pos
    n_pos_pred <- TP+FP
    n_neg_pred <- TN+FN
    
    
    cm.values.ci["Accuracy",] <- binom.confint(x=n_tot*mean(cm.values$Accuracy),
                                               n=n_tot,
                                               methods = "exact")[,4:6]
    
    cm.values.ci["Sensitivity",] <- binom.confint(x=n_pos*mean(cm.values$Sensitivity),
                                               n=n_pos,
                                               methods = "exact")[,4:6]
    
    cm.values.ci["Specificity",] <- binom.confint(x=n_neg*mean(cm.values$Specificity),
                                               n=n_neg,
                                               methods = "exact")[,4:6]
    
    cm.values.ci["PPV",] <- colMeans(binom.confint(x=n_pos_pred*mean(cm.values$PPV),
                                                   n=n_pos_pred,
                                                   methods = "exact")[,4:6])
    
    cm.values.ci["NPV",] <- colMeans(binom.confint(x=n_neg_pred*mean(cm.values$NPV),
                                                   n=n_neg_pred,
                                                   methods = "exact")[,4:6])
    
    
    auc.labels <- matrix(rep(groups.lda, repetitions), nrow(predict.probs), repetitions)
    cm.values.ci["AUROC",] <- unlist(ci.cvAUC(predict.probs, auc.labels)[c("cvAUC", "ci")])
    
    cm.values.ci <- round(cm.values.ci*100, 2)
    
    print(cm.values.ci)
     
     ########################
     
     # ROC curves
     roc.results.df <- as.data.frame(roc.results)
     roc.results.df$median_FPR <- apply(roc.results[,1,], 1, median, na.rm=T)
     roc.results.df$median_TPR <- apply(roc.results[,2,], 1, median, na.rm=T)
     
     roc.results.df <- rbind(roc.results.df, rep(1, ncol(roc.results.df)))
     
     p <- ggplot(data=roc.results.df)
     
     for (i in c(0:(repetitions-1))) {
       p <- p + geom_line(aes_string(x=names(roc.results.df)[i*2+1], y=names(roc.results.df)[i*2+2], col=shQuote("lg")))
     }
     
     
     p <- p + geom_line(aes(x=median_FPR, y=median_TPR, col="black"), size=1) + 
       geom_segment(x=0, y=0, xend=1, yend=1, linetype="dashed") +
       scale_color_manual("Legend", values=c(lg="#d3d3d3", black="black")) + 
       theme_bw() + theme(legend.position = "none") + 
       labs(x="False Positive Rate", y="True Positive Rate") + 
       scale_x_continuous(expand=expansion(mult=0.01)) + 
       scale_y_continuous(expand=expansion(mult=0.01))
     
     options(warn=-1)
     print(p)
     options(warn=0)
     
     ############################
     
     # LDA scores by group
     LDA.scores.df <- data.frame(groups.lda, apply(LDA.scores, 1 , function(x) median(x)))
     colnames(LDA.scores.df) <- c("Group", "Score")
     
     p <- ggplot(data=LDA.scores.df, aes(x=Group, y=Score)) + 
       stat_boxplot(geom="errorbar", width=0.3) +
       geom_boxplot(outlier.shape = NA, width=0.7) + 
       geom_jitter(aes(col=Group), width=0.1, show.legend=F) +
       geom_hline(yintercept=0, linetype="dashed") +
       scale_x_discrete(breaks=groups, labels=legend) + 
       ggtitle("LD scores by group") + 
       theme_bw() + 
       theme(plot.title = element_text(hjust = 0.5))
     
     print(p) 
     
     ###############################################
     
     # LDA loadings
     LDA.loadings.df <- data.frame(
       Data$Wavenumber_min_max,
       apply(LDA.loadings, 1, median),
       apply(LDA.loadings, 1, function(x) {quantile(x, 0.25)}),
       apply(LDA.loadings, 1, function(x) {quantile(x, 0.75)})
     )
     
     colnames(LDA.loadings.df) <- c("Wavenumber", "Median", "FirstQuart", "ThirdQuart")
     
     p <- ggplot(data=LDA.loadings.df, aes(x=Wavenumber, y=Median)) + 
          geom_ribbon(aes(ymin=FirstQuart, ymax=ThirdQuart), col="lightgray", fill="lightgray") +
          geom_line(size=0.5) +
          geom_hline(yintercept=0, linetype="dashed") + theme_bw() +
          stat_peaks(geom="text", span=41, color="black", x.label.fmt="%.0f", ignore_threshold=0.7, angle=90, vjust=0.5, hjust=-0.75) + 
          stat_valleys(geom="text", span=41, color="black", x.label.fmt="%.0f", ignore_threshold=0.35, angle=90, vjust=0.5, hjust=1.75) +
          scale_x_continuous(name=bquote(bold("Raman shift" ~(cm^-1))), expand=expansion(add=0)) + 
          scale_y_continuous(name=bquote(bold("Loadings of PCA-LDA")), expand=expansion(mult=0.25))
     
     print(p)
     
     
     ##############################
     
     # Median spectra by group with interquart range
     spectra_median.df <- data.frame(matrix(NA, nrow=length(Data$Wavenumber_min_max), ncol=1+3*length(groups)))
     
     colnames(spectra_median.df) <- c("Wavenumber", paste(rep(c("Median","Q1","Q3"), times=length(groups)), rep(groups, each=3), sep="_"))
     
     spectra_median.df$Wavenumber <- Data$Wavenumber_min_max
     
     for (group in groups) {
       spectra_median.df[[paste0("Median_", group)]] <- apply(data.lda[groups.lda == group,], 2, median)
       spectra_median.df[[paste0("Q1_", group)]] <- apply(data.lda[groups.lda == group,], 2, function(x) {quantile(x, 0.25)})
       spectra_median.df[[paste0("Q3_", group)]] <- apply(data.lda[groups.lda == group,], 2, function(x) {quantile(x, 0.75)})
     }
     
     spectra.combined <- data.frame(matrix(NA, nrow=length(Data$Wavenumber_min_max), ncol=2))
     colnames(spectra.combined) <- c("WN", "INT")
     
     spectra.combined$WN <- Data$Wavenumber_min_max
     spectra.combined$INT <- apply(spectra_median.df[,paste0("Median_", groups)], 1, max)
     
     spectra_median.df <- reshape(spectra_median.df, varying=2:ncol(spectra_median.df), sep="_", direction = "long", timevar="Group")
     
     p <- ggplot(spectra_median.df, aes(x=Wavenumber, y=Median)) + 
          geom_line(aes(col=Group), size=0.5) +
          geom_ribbon(aes(ymin=Q1, ymax=Q3, fill=Group), alpha=0.2, show.legend=FALSE) +
          # stat_peaks(data=spectra.combined, aes(x=WN, y=INT), geom="text", span=33, 
          #            color="black", x.label.fmt="%.0f", ignore_threshold=0.02, angle=90, 
          #            vjust=0.5, hjust=-0.75) +
          scale_color_discrete(label=legend) +
          scale_fill_discrete(label=legend) +
          scale_x_continuous(expand=expansion(0)) + 
          scale_y_continuous(expand=expansion(c(0,0.25))) +
          theme_bw() + theme(legend.position = "bottom") + 
          labs(x=bquote(bold("Raman shift" ~(cm^-1))), y=bquote(bold("Intensity")))
     
     print(p)
     
     #############################
     
     # Difference spectrum
     if (length(groups) > 2) {
       print("Difference Spectrum is only possible with exactly two groups")
     } else {
       
       diff.spectrum.df <- data.frame(Wavenumber=Data$Wavenumber_min_max, Median=NA, Q1=NA, Q3=NA)
       
       diff.spectrum.df$Median <- 
         unlist(subset(spectra_median.df, Group==pos.class, select=Median) - 
         subset(spectra_median.df, Group!=pos.class, select=Median))
       
       diff.spectrum.df$Q1 <- 
         unlist(subset(spectra_median.df, Group==pos.class, select=Q1) - 
         subset(spectra_median.df, Group!=pos.class, select=Q3))
       
       diff.spectrum.df$Q3 <- 
         unlist(subset(spectra_median.df, Group==pos.class, select=Q3) - 
         subset(spectra_median.df, Group!=pos.class, select=Q1))
       
       p <- ggplot(diff.spectrum.df, aes(x=Wavenumber, y=Median)) + 
            geom_line(size=0.5) + 
            geom_ribbon(aes(ymin=Q1, ymax=Q3), alpha=0.25) +
            geom_hline(yintercept=0, linetype="dashed") +
            stat_peaks(geom="text", span=31, color="black", 
                       x.label.fmt="%.0f", ignore_threshold=0.65, 
                       angle=90, vjust=0.5, hjust=-0.75) + 
            stat_valleys(geom="text", span=31, color="black", 
                         x.label.fmt="%.0f", ignore_threshold=0.4, 
                         angle=90, vjust=0.5, hjust=1.75) +
            scale_x_continuous(expand=expansion(0)) +
            theme_bw() + 
            labs(x=bquote(bold("Wavenumber" ~(cm^-1))), y=bquote(bold(~Delta~ "Intensity")))
       
       print(p)
       
       
     }
     
  }
  
}


end_time = format(Sys.time(), '%X') # get end time of analysis
# print time difference between start time and end time
print(paste("Time used for analysis:", round(as.difftime(end_time, units = "mins")-as.difftime(start_time, units = "mins"),digits=2),"minutes"))

# y1lim <- c(0,0.2)
# y2lim <- c(-1.5,1)
# b <- diff(y1lim)/diff(y2lim)
# a <- y1lim[1] - b*y2lim[1]
# 
# 
# p <- ggplot(spectra_median.df, aes(x=Wavenumber, y=Median)) + 
#   geom_line(aes(col=Group), size=0.7) +
#   geom_line(data=LDA.loadings.df, mapping=aes(x=Wavenumber, y=a+Median*b), color="gray") + 
#   geom_hline(yintercept=a, linetype="dashed") +
#   geom_ribbon(aes(ymin=Q1, ymax=Q3, fill=Group), alpha=0.2, show.legend=FALSE) +
#   stat_peaks(data=spectra.combined, aes(x=WN, y=INT), geom="text", span=35, 
#              color="black", x.label.fmt="%.0f", ignore_threshold=0.05, angle=90, 
#              vjust=0.5, hjust=-0.75) +
#   scale_color_discrete(label=legend) +
#   scale_fill_discrete(label=legend) +
#   scale_x_continuous(expand=expansion(0)) + 
#   scale_y_continuous(expand=expansion(c(0,0.25)), name="Intensity",
#                      sec.axis=sec_axis(~ (. - a)/b, name="Loadings")) +
#   theme_bw() + theme(legend.position = "bottom") + 
#   labs(x=bquote("Raman shift" ~(cm^-1)))
# 
# print(p)
# 
