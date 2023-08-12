### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
### Differential Methylation across Ceccareli's glioma subtype ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

#@@@ TCGA glioma samples
#@@@ Aim: tune DNAmet probes (features) to classify those 6 glioma subtypes based on the label establish by Ceccarelli et al 2016
#@@@ probes from 450k and 27k illumina array 
#**IMPORTANT**#
#@@@ If we don't get a good fine tuned model from that we should try to add EPIC array. Maybe only intersections within EPIC and 450k ? 

# Load packages 
library(dplyr)


# Load DNAmtx 
load("/media/hd/maycon/Glioma_classifier/DNAmet_mtx_no_maskprobes_neither_chrprobes.rda")
DNAmet_mtx_sub[1:4, 1:4]
dim(DNAmet_mtx_sub) #23493   932

# Load metadata
library(readr)
metadata <- read_csv("/media/hd/maycon/Glioma_classifier/mmc2.csv")
metadata = data.frame(metadata)
metadata %>% head()
dim(metadata) #1122 51

metadata <- metadata %>% 
  filter(Case %in% colnames(DNAmet_mtx_sub)) #keep only DNAmet samples
dim(metadata) #932  51



# Anova (Tathi's code) -------------------------

# Aggregating 'PA-like' subtype into 'LGm6-GBM' (because now there are 6 subtypes instead of 7)
metadata$Ceccarelli_six_subtypes <- as.character(metadata$Supervised.DNA.Methylation.Cluster)
metadata[metadata$Ceccarelli_six_subtypes %in% 'PA-like', ]$Ceccarelli_six_subtypes <- 'LGm6-GBM'
metadata$Ceccarelli_six_subtypes %>% table()

# Handling NAs 
table(is.na(DNAmet_mtx_sub))
# FALSE     TRUE 
# 20091130  1804346 
DNAmet_mtx_sub <- na.omit(DNAmet_mtx_sub) #removing probes (rows) within at least one beta-value NA
dim(DNAmet_mtx_sub) #18789 932
table(is.na(DNAmet_mtx_sub))
# FALSE 
# 17511348

# Ordering datasets
meta <- metadata
rownames(meta) <- meta$Case
data <- DNAmet_mtx_sub
data <- data[, rownames(meta)]
head(data)

identical(colnames(data), rownames(meta)) #T (samples in the same order. data = your data. meta = meta data)

require(parallel)
values <- as.data.frame(t(data)) #[samples, features]
w.p.values <- unlist(mclapply(values,
                              function(probe) {
                                probe <- data.frame(probe)
                                probe$Ceccarelli_six_subtypes <- meta$Ceccarelli_six_subtypes #column with your groups
                                colnames(probe)[1] <- "value"
                                if(nrow(na.omit(probe)) > 1){ #se você só tiver NA no seu objeto ele não realiza o teste
                                  test <- summary(aov(value ~ Ceccarelli_six_subtypes, data=probe)) #faz o teste anova
                                  if(ncol(test[[1]]) > 3)  #dimensão onde o p-value está armazenado
                                    return(test[[1]][[5]][[1]])
                                  else
                                    return(NA)
                                }
                                else
                                  return(NA)
                                
                              }, mc.cores=8))
w.p.values.adj <- p.adjust(w.p.values, method = "BH")
save(w.p.values, w.p.values.adj, file="/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/anova.pvalue.Rda")
# what about the tukey to correct ANOVA test? ??? Answer: BH method has already been choose

identical(rownames(data),
          names(w.p.values.adj)) #TRUE
data <- as.data.frame(data)
data$p_val_adj <- NA
data$p_val_adj <- as.vector(w.p.values.adj)
data[data$p_val_adj < 0.01, ]

meta$Ceccarelli_six_subtypes %>% table()
# Classic-like            Codel      G-CIMP-high       G-CIMP-low 
# 148                     174         249               25 
# LGm6-GBM                Mesenchymal-like 
# 67                      215

Classic_like <-  meta[meta$Ceccarelli_six_subtypes %in% 'Classic-like', ]$Case
Codel <-  meta[meta$Ceccarelli_six_subtypes %in% 'Codel', ]$Case
G_CIMP_high <- meta[meta$Ceccarelli_six_subtypes %in% 'G-CIMP-high', ]$Case
G_CIMP_low <- meta[meta$Ceccarelli_six_subtypes %in% 'G-CIMP-low', ]$Case
LGm6_GBM <- meta[meta$Ceccarelli_six_subtypes %in% 'LGm6-GBM', ]$Case
Mesenchymal_like <- meta[meta$Ceccarelli_six_subtypes %in% 'Mesenchymal-like', ]$Case

# Classic_like Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_Codel_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_G_CIMP_high_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_G_CIMP_low_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_LGm6_GBM_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_Mesenchymal_like_ClassicOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation > 0.3, ] %>% dim() # 4 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_Classic_Codel_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation > 0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_high_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation > 0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_low_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation > 0.3, ] %>% rownames(),
                 
                  DiffMean_Classic_LGm6_GBM_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation > 0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_Mesenchymal_like_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation > 0.3, ] %>% rownames()
                  )

upset(fromList(listInput_hyper), order.by = "freq", nsets = 6) 
# 4/5 overlapped comparisons probeset = 85 probes; 6/6 overlapped groups probeset = 4. Take them all (85+4 = 89 probes) to be the Classic-like subtype probe signature/features 


data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation < -0.3, ] %>% dim() # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Classic_Codel_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_high_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_low_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_LGm6_GBM_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_Mesenchymal_like_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 6) # to many groups overlapping on same probes within DeafMean < -0.3


# Extract the 85+4 hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #413 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #413 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # only probes overlapped on 5/5 comparisons 
x1[ rowSums(x$New_data) == 4] # only probes overlapped on 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

Classic_probeset_89 <- c(x1[ rowSums(x$New_data) == 5], x1[ rowSums(x$New_data) == 4])



# Codel Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_Classic_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_G_CIMP_high_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_G_CIMP_low_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_LGm6_GBM_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_Mesenchymal_like_CodelOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


# Hyper methylated probes
data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation > 0.3 &
           data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation > 0.3 &
           data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation > 0.3 &
           data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation > 0.3 &
           data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation > 0.3, ] %>% dim() # 14 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_Codel_Classic_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Codel_G_CIMP_high_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Codel_G_CIMP_low_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Codel_LGm6_GBM_CodelOrientation = data_sig[data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Codel_Mesenchymal_like_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation > 0.3, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# 4/5 overlapped comparisons probeset = 365 probes; 5/5 overlapped comparisons probeset = 14. Take them all (85+4 = 89 probes) to be the Classic-like subtype probe signature/features 


# Hypo methylated probes
data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation < -0.3, ] %>% dim() #  # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Codel_Classic_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_G_CIMP_high_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_G_CIMP_low_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_LGm6_GBM_CodelOrientation = data_sig[data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_Mesenchymal_like_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # to many groups overlapping on same probes within DeafMean < -0.3


# Extract the 365+14 hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #2460 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #2460 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # only probes overlapped on 5/5 comparisons 
x1[ rowSums(x$New_data) == 4] # only probes overlapped on 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

Codel_probeset_379p <- c(x1[ rowSums(x$New_data) == 5], x1[ rowSums(x$New_data) == 4])




# G_CIMP_high Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


# Hyper methylated probes
data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation > 0.3, ] %>% dim() # 0 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_G_CIMP_high_Classic_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_Codel_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation > 0.3, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# 2 4/5 overlapped comparisons probeset = i) 177 probes and ii) 4 probes. Take them all (177+4 = 181 probes) to be the Classic-like subtype probe signature/features 


# Hypo methylated probes
data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation < -0.3, ] %>% dim() # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_G_CIMP_high_Classic_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation < -0.3 , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_Codel_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation < -0.3 , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation < -0.3 , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation < -0.3 , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation < -0.3 , ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # to many groups overlapping on same probes within DeafMean < -0.3


# Extract the 365+14 hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #2012 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #2012 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none 
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 groups

# length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret. - don't need to do this

GcimpHigh_probeset_181p <- c(#x1[ rowSums(x$New_data) == 5], - don't need to do this
                         x1[ rowSums(x$New_data) == 4])




# G_CIMP_low Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


# Hyper methylated probes
data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation> 0.2 &
           data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation> 0.2, ] %>% dim() # 5 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_G_CIMP_low_Classic_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_Codel_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation > 0.2, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# We had to set the threshold to 0.2 otherwise we couldn't find enough probes (~ 90 probes as in the other comparisons we've done)


# Hypo methylated probes
data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation < -0.2 &
           data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation < -0.2 , ] %>% dim() # 69 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_G_CIMP_low_Classic_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_Codel_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation < -0.2, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # In this one we got some  hypo methylated probes ! So for the GcimpLow probeset we're getting both hyper and hypo methylated probes 


# Extract hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #1798 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #1798 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # all probes present it only 5/5 comparisons
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

hyper_probes <- c(x1[ rowSums(x$New_data) == 5],
                  x1[ rowSums(x$New_data) == 4])

# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #2426 probes somehow overlapped
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #2426 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # all probes present it only 5/5 comparisons
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

hypo_probes <- c(x1[ rowSums(x$New_data) == 5],
                  x1[ rowSums(x$New_data) == 4])

length(intersect(hyper_probes, hypo_probes)) # zero. Corret.
length(hyper_probes) #64
length(hypo_probes) #115
#  115+64 = 179 probes at total
GcimpLow_probeset_179p <- c(hyper_probes, hypo_probes)





# LGm6_GBM Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


# Hyper methylated probes
data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation> 0.2 &
           data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation > 0.2, ] %>% dim() # 0 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation > 0.2, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# We had to set the threshold to 0.2 otherwise we couldn't find enough probes (~ 90 probes as in the other comparisons we've done)


# Hypo methylated probes
data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation < -0.2, ] %>% dim() # 28 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation < -0.2, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # In this one (LGm6) we got ONLY  hypo methylated probes !   


# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #4216 probes somehow overlapped
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #4216 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # all probes present it only 5/5 comparisons
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.
length(x1[ rowSums(x$New_data) == 5]) #28
length(x1[ rowSums(x$New_data) == 4]) #206
# 28 + 206 = 234

LGm6_probeset_234p <- c(x1[ rowSums(x$New_data) == 5],
                        x1[ rowSums(x$New_data) == 4])





# Mesenchymal_like Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_Classic_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_Codel_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]


#**Alert**#
dim(data) #18789   probes
dim(data_sig) #17701   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 17701 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


# Hyper methylated probes
data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation> 0.3 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation, ] %>% dim() # 0 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_Mesenchymal_Classic_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_Codel_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation > 0.2, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# too few probes and one comparison fell out for not having anyprobe in common to the other comparisons

# Hypo methylated probes
data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation, ] %>% dim() # 19 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Mesenchymal_Classic_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_Codel_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation < -0.2, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # In this one (Mesenchymal) we got ONLY  hypo methylated probes but still too few probes (20 probes) ! I found more intresting keep few specific probes than picking a lower threshold (eg 0.1) because it might increase non-specific probes into the model


# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #3828 probes somehow overlapped
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #3828 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

# length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret. - not necessary 
length(x1[ rowSums(x$New_data) == 4]) #20
# 20 probes at total

Mesenchymal_probeset_20p <- c(x1[ rowSums(x$New_data) == 4])


DMP_anova_diffmean_subtypes_probeset_list <- list(Classic_probeset_89 = Classic_probeset_89,
                                               Codel_probeset_379p = Codel_probeset_379p,
                                               GcimpHigh_probeset_181p= GcimpHigh_probeset_181p,
                                               GcimpLow_probeset_179p = GcimpLow_probeset_179p,
                                               LGm6_probeset_234p = LGm6_probeset_234p,
                                               Mesenchymal_probeset_20p = Mesenchymal_probeset_20p)

saveRDS(DMP_anova_diffmean_subtypes_probeset_list, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DMP_anova_diffmean_subtypes_probeset_list.rds')


# Heatmap visualization -----------
head(metadata) # from the beginning of this code 

my_sample_col <- data.frame(row.names = metadata$Case, Subtypes_six = metadata$Ceccarelli_six_subtypes, Subtypes_seven = metadata$Supervised.DNA.Methylation.Cluster, DNAmet_cluster = metadata$IDH.specific.DNA.Methylation.Cluster) # patient ID in row names; any columns are for sample information

data_sig # DNAmtx 
DNAmtx <- data_sig[, colnames(data_sig) %in% metadata$Case] #keep only beta-values

length(DMP_anova_diffmean_subtypes_probeset_list) 
tune_probset <- c(DMP_anova_diffmean_subtypes_probeset_list[[1]],
  DMP_anova_diffmean_subtypes_probeset_list[[2]],
  DMP_anova_diffmean_subtypes_probeset_list[[3]],
  DMP_anova_diffmean_subtypes_probeset_list[[4]],
  DMP_anova_diffmean_subtypes_probeset_list[[5]],
  DMP_anova_diffmean_subtypes_probeset_list[[6]])
  


library(pheatmap)
p = pheatmap(DNAmtx[rownames(DNAmtx) %in% tune_probset, ], 
             #annotation_row = my_probe_col, 
             annotation_col = my_sample_col,
             show_rownames = FALSE,
             main = 'Glioma subtypes tuned probes'); p

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p, "/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/finetunedprobes_heatmap.pdf")

# Ceccareli pan_glioma_probes
library(openxlsx)
pan_glioma_probes <- read.xlsx("/media/hd/maycon/Glioma_classifier/PanGlioma_MethylationSignatures.xlsx", sheet = 1) 
pan_glioma_probes <- pan_glioma_probes$`1,300.pan-glioma.tumor.specific.probes.(Figure.2A)`


library(pheatmap)
p_cecca = pheatmap(DNAmtx[rownames(DNAmtx) %in% pan_glioma_probes, ], 
             #annotation_row = my_probe_col, 
             annotation_col = my_sample_col,
             show_rownames = FALSE,
             main = 'Glioma subtypes tuned probes'); p_cecca

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p_cecca, "/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/ceccarelli_pangliomaprobes_heatmap.pdf")

# We still need to clean up intersected probes across probesets 
# ATTENTION: maybe the intesect probes are methylated in different orientations. Plots a headmap to evaluate the 6 probesets before start removing any probes
# 1. plot upsetplot of all 6 probesets
# 2. extract only probes in x1[ rowSums(x$New_data) == 1]. Rembemer now it's 6 groups

# Unique probes from tuned probeset
library(UpSetR)
listInput_allprobeset <- list(Classic_probeset_89 = DMP_anova_diffmean_subtypes_probeset_list[[1]],
                       
                       Codel_probeset_379p = DMP_anova_diffmean_subtypes_probeset_list[[2]],
                       
                       GcimpHigh_probeset_181p = DMP_anova_diffmean_subtypes_probeset_list[[3]],
                       
                       GcimpLow_probeset_179p =  DMP_anova_diffmean_subtypes_probeset_list[[4]],
                       
                       LGm6_probeset_234p = DMP_anova_diffmean_subtypes_probeset_list[[5]],
                       
                       Mesenchymal_probeset_20p = DMP_anova_diffmean_subtypes_probeset_list[[6]]
                       
)

upset(fromList(listInput_allprobeset), order.by = "freq", nsets = 6)

# Extract non-intersected probes
x <- upset(fromList(listInput_allprobeset), nsets = 6)
x$New_data[1:5, 1:5]
dim(x$New_data) #901 
x1 <- unlist(listInput_allprobeset, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #901 

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 1] # unique probes for each group (glioma subtype)


tune_probset_unique <- x1[ rowSums(x$New_data) == 1] 


library(pheatmap)
p_uniq = pheatmap(DNAmtx[rownames(DNAmtx) %in% tune_probset_unique, ], 
             #annotation_row = my_probe_col, 
             annotation_col = my_sample_col,
             show_rownames = FALSE,
             main = 'Glioma subtypes tuned probes'); p_uniq

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(p_uniq, "/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/finetunedprobes_uniqueprobes_heatmap.pdf")



library(ComplexHeatmap)
library(matlab)
#label 'my_sample_col'
#                    Subtypes_six   Subtypes_seven DNAmet_cluster
# TCGA-CS-4938      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4941 Mesenchymal-like Mesenchymal-like       IDHwt-K2
# TCGA-CS-4942      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4943      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4944      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-5390            Codel            Codel      IDHmut-K3
# Current label

### Pallete 
library(viridis)
n_colors <- length(table(my_sample_col$Subtypes_six))
Subtypes_six_pal <- viridis(n = n_colors, option = "plasma", direction = -1)
n_colors <- length(table(my_sample_col$Subtypes_seven))
Subtypes_seven_pal <- viridis(n = n_colors, option = "inferno", direction = -1)
n_colors <- length(table(my_sample_col$DNAmet_cluster))
DNAmet_cluster_pal <- viridis(n = n_colors, option = "cividis", direction = -1)

Subtypes_six_lvs <- names(table(my_sample_col$Subtypes_six))
Subtypes_seven_lvs <- names(table(my_sample_col$Subtypes_seven))
DNAmet_cluster_lvs <- names(table(my_sample_col$DNAmet_cluster))

# Name yout columns and it's levels
column_name_1 <- "Subtypes_six"
levels_1 <- Subtypes_six_lvs

column_name_2 <- "Subtypes_seven"
levels_2 <- Subtypes_seven_lvs

column_name_3 <- "DNAmet_cluster"
levels_3 <- DNAmet_cluster_lvs

# Generate color vector
colors_1 <- Subtypes_six_pal
colors_2 <- Subtypes_seven_pal
colors_3 <- DNAmet_cluster_pal
# Assign names to color vector
names(colors_1) <- levels_1
names(colors_2) <- levels_2
names(colors_3) <- levels_3
# Create named list
color_mapping <- list(column_name_1 = colors_1, 
                      column_name_2 = colors_2,
                      column_name_3 = colors_3)


top.anno = HeatmapAnnotation(df = my_sample_col, 
                             col= color_mapping,
                             show_annotation_name = T, annotation_name_gp = gpar(fontsize=7),
                             na_col= "white")


#hm_tun_pset <- Heatmap(as.matrix(DNAmtx[rownames(DNAmtx) %in% tune_probset, ]),
#hm_uniq_tun_pset <- Heatmap(as.matrix(DNAmtx[rownames(DNAmtx) %in% tune_probset_unique, ]),
hm_pan_glioma <- Heatmap(as.matrix(DNAmtx[rownames(DNAmtx) %in% pan_glioma_probes, ]),
                               cluster_columns=T,
                               cluster_rows = T,
                               clustering_method_rows="complete",  
                               show_row_names = F,
                               row_names_gp = gpar(fontsize = 7),
                               show_column_names = F,
                               name = "CpG probes methylation",
                               col= jet.colors(75),
                               row_title = "CpG probes",
                               #row_names_gp = gpar(fontsize = 12),
                               #column_title = "Glioma subtypes tuned probes",
                               #column_title = "Glioma subtypes unique tuned probes",
                               column_title = "Glioma subtypes pan_glioma probes (Ceccarelli)",
                               #split=label.EMT.genes$E_ou_M,
                               #row_names_side="left",
                               
                               top_annotation = top.anno,
                               #row_names_gp = gpar(fontsize = 8),
                               #column_names_gp = gpar(fontsize = 12),
                               heatmap_legend_param = list(
                                 color_bar = 'continuous',
                                 legend_direction = 'vertical',
                                 legend_width = unit(8, 'cm'),
                                 legend_height = unit(5.0, 'cm'),
                                 title_position = 'leftcenter-rot',
                                 title_gp=gpar(fontsize = 12, fontface = 'bold'),
                                 labels_gp=gpar(fontsize = 12, fontface = 'bold')))

#print(hm_tun_pset)
#print(hm_uniq_tun_pset)
print(hm_pan_glioma)

### For left side annotation
# RowAnn <- HeatmapAnnotation(df = label, col=list("Epi.Mes.CellC" = c("Epi" = "blue", "Mes" = "black", "c.Cycle" = "cadetblue1")),
#                             show_annotation_name = T, annotation_name_gp = gpar(fontsize=7),
#                             na_col= "white",which = "row",show_legend =T)


### For left side annotation
# draw(hm_EMT.cCycle.MALTA + RowAnn) 

pdf("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/panglioma_probeset.pdf",width = 5, height = 10)
### For left side annotation
# draw(hm_EMT.cCycle.MALTA, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()




#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#**Tathi's advice**#
# - use knn to replace NA 
# - be MORE stringest on DiffMean threshold. We don't want to many probes
# - run Random Forest with the probe list we got until know including Ceccarelli's (pan glioma probset)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Prepare the data again ------------
# Load DNAmtx 
load("/media/hd/maycon/Glioma_classifier/DNAmet_mtx_no_maskprobes_neither_chrprobes.rda")
DNAmet_mtx_sub[1:4, 1:4]
dim(DNAmet_mtx_sub) #23493   932

# Load metadata
library(readr)
metadata <- read_csv("/media/hd/maycon/Glioma_classifier/mmc2.csv")
metadata = data.frame(metadata)
metadata %>% head()
dim(metadata) #1122 51

metadata <- metadata %>% 
  filter(Case %in% colnames(DNAmet_mtx_sub)) #keep only DNAmet samples
dim(metadata) #932  51


# Aggregating 'PA-like' subtype into 'LGm6-GBM' (because now there are 6 subtypes instead of 7)
metadata$Ceccarelli_six_subtypes <- as.character(metadata$Supervised.DNA.Methylation.Cluster)
metadata[metadata$Ceccarelli_six_subtypes %in% 'PA-like', ]$Ceccarelli_six_subtypes <- 'LGm6-GBM'
metadata$Ceccarelli_six_subtypes %>% table()

# Handling NAs 
table(is.na(DNAmet_mtx_sub))

library(impute)
DNAmet_mtx_sub <- impute.knn(as.matrix(DNAmet_mtx_sub), k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)[[1]]
# Warning message:
#   In knnimp(x, k, maxmiss = rowmax, maxp = maxp) :
#   1919 rows with more than 80 % entries missing;
# mean imputation used for these rows

table(is.na(DNAmet_mtx_sub)) # all NAs have been replaced


# Ordering datasets
meta <- metadata
rownames(meta) <- meta$Case
data <- DNAmet_mtx_sub
data <- data[, rownames(meta)]
head(data)

identical(colnames(data), rownames(meta)) #T (samples in the same order. data = your data. meta = meta data)

require(parallel)
values <- as.data.frame(t(data)) #[samples, features]
w.p.values <- unlist(mclapply(values,
                              function(probe) {
                                probe <- data.frame(probe)
                                probe$Ceccarelli_six_subtypes <- meta$Ceccarelli_six_subtypes #column with your groups
                                colnames(probe)[1] <- "value"
                                if(nrow(na.omit(probe)) > 1){ #se você só tiver NA no seu objeto ele não realiza o teste
                                  test <- summary(aov(value ~ Ceccarelli_six_subtypes, data=probe)) #faz o teste anova
                                  if(ncol(test[[1]]) > 3)  #dimensão onde o p-value está armazenado
                                    return(test[[1]][[5]][[1]])
                                  else
                                    return(NA)
                                }
                                else
                                  return(NA)
                                
                              }, mc.cores=8))
w.p.values.adj <- p.adjust(w.p.values, method = "BH")
save(w.p.values, w.p.values.adj, file="/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/KNN_NAreplace_anova.pvalue.Rda")
# what about the tukey to correct ANOVA test? ??? Answer: BH method has already been choose

identical(rownames(data),
          names(w.p.values.adj)) #TRUE
data <- as.data.frame(data)
data$p_val_adj <- NA
data$p_val_adj <- as.vector(w.p.values.adj)
data[data$p_val_adj < 0.01, ] %>% dim()

meta$Ceccarelli_six_subtypes %>% table()
# Classic-like            Codel      G-CIMP-high       G-CIMP-low 
# 148                     174         249               25 
# LGm6-GBM                Mesenchymal-like 
# 67                      215

Classic_like <-  meta[meta$Ceccarelli_six_subtypes %in% 'Classic-like', ]$Case
Codel <-  meta[meta$Ceccarelli_six_subtypes %in% 'Codel', ]$Case
G_CIMP_high <- meta[meta$Ceccarelli_six_subtypes %in% 'G-CIMP-high', ]$Case
G_CIMP_low <- meta[meta$Ceccarelli_six_subtypes %in% 'G-CIMP-low', ]$Case
LGm6_GBM <- meta[meta$Ceccarelli_six_subtypes %in% 'LGm6-GBM', ]$Case
Mesenchymal_like <- meta[meta$Ceccarelli_six_subtypes %in% 'Mesenchymal-like', ]$Case


# Classic_like Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_Codel_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_G_CIMP_high_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_G_CIMP_low_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_LGm6_GBM_ClassicOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Classic_Mesenchymal_like_ClassicOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #23493   probes
dim(data_sig) #22180   probes - it seems there is too much significant probes... 
# Obs: it's ok to get the same 22180 significant probes each time we're subseting probes for each glioma subtype because we're just filering the same statistics we had made before within all groups together 


data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation > 0.3 &
       data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation > 0.3, ] %>% dim() # 4 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_Classic_Codel_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation > 0.4, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_high_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation > 0.4, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_low_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation > 0.4, ] %>% rownames(),
                 
                  DiffMean_Classic_LGm6_GBM_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation > 0.4, ] %>% rownames(),
                  
                  DiffMean_Classic_Mesenchymal_like_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation > 0.4, ] %>% rownames()
                  )

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 


data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation < -0.3 &
       data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation < -0.3, ] %>% dim() # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Classic_Codel_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Codel_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_high_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_high_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_G_CIMP_low_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_G_CIMP_low_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_LGm6_GBM_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_LGm6_GBM_ClassicOrientation < -0.3, ] %>% rownames(),
                  
                  DiffMean_Classic_Mesenchymal_like_ClassicOrientation = data_sig[data_sig$DiffMean_Classic_Mesenchymal_like_ClassicOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5)  


# Extract the hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #167 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #167 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none
x1[ rowSums(x$New_data) == 4] # only probes overlapped on 4/5 comparisons

Classic_probeset_36p <- c(x1[ rowSums(x$New_data) == 4])
# only hyper probes; diffmean > 0.4; pvalue < 0.05; 4/5 comparisons had this probeset


# Codel Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_Classic_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_G_CIMP_high_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_G_CIMP_low_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_LGm6_GBM_CodelOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Codel) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Codel_Mesenchymal_like_CodelOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #23493
dim(data_sig) #22180


# Hyper methylated probes
data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation > 0.5 &
           data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation > 0.5 &
           data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation > 0.5 &
           data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation > 0.5 &
           data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation > 0.5, ] %>% dim()

library(UpSetR)
listInput_hyper <- list(DiffMean_Codel_Classic_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation > 0.5, ] %>% rownames(),
                        
                        DiffMean_Codel_G_CIMP_high_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation > 0.5, ] %>% rownames(),
                        
                        DiffMean_Codel_G_CIMP_low_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation > 0.5, ] %>% rownames(),
                        
                        DiffMean_Codel_LGm6_GBM_CodelOrientation = data_sig[data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation > 0.5, ] %>% rownames(),
                        
                        DiffMean_Codel_Mesenchymal_like_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation > 0.5, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 



# Hypo methylated probes
data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation < -0.3 &
           data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation < -0.3, ] %>% dim() #  # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Codel_Classic_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Classic_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_G_CIMP_high_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_high_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_G_CIMP_low_CodelOrientation = data_sig[data_sig$DiffMean_Codel_G_CIMP_low_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_LGm6_GBM_CodelOrientation = data_sig[data_sig$DiffMean_Codel_LGm6_GBM_CodelOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Codel_Mesenchymal_like_CodelOrientation = data_sig[data_sig$DiffMean_Codel_Mesenchymal_like_CodelOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # to few comparisons overlapped


# Extract the hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #605 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #605 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # only probes overlapped on 5/5 comparisons 
x1[ rowSums(x$New_data) == 4] # only probes overlapped on 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

Codel_probeset_10p <- c(x1[ rowSums(x$New_data) == 5], x1[ rowSums(x$New_data) == 4])
# only hyper probes; diffmean > 0.5; pvalue < 0.05; 5/5 and 4/5 comparisons had this probeset




# G_CIMP_high Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_high) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #23493   
dim(data_sig) #22180  


# Hyper methylated probes
data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation > 0.3 &
           data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation > 0.3, ] %>% dim() 
library(UpSetR)
listInput_hyper <- list(DiffMean_G_CIMP_high_Classic_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_Codel_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation > 0.4, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 

# Hypo methylated probes
data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation < -0.3 &
           data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation < -0.3, ] %>% dim() # 0 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_G_CIMP_high_Classic_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Classic_GcimpHighOrientation < -0.4 , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_Codel_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Codel_GcimpHighOrientation < -0.4  , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_G_CIMP_low_GcimpHighOrientation < -0.4  , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_LGm6_GBM_GcimpHighOrientation < -0.4  , ] %>% rownames(),
                       
                       DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation = data_sig[data_sig$DiffMean_G_CIMP_high_Mesenchymal_like_GcimpHighOrientation < -0.4  , ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) # 


# Extract the hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #1079 
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #1079 

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none 
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 groups

# length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret. - don't need to do this

GcimpHigh_probeset_18p <- c(x1[ rowSums(x$New_data) == 4])
# only hyper probes; diffmean > 0.4; pvalue < 0.05; 4/5 comparisons had this probeset


# G_CIMP_low Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(G_CIMP_low) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #23493   
dim(data_sig) #22180   


# Hyper methylated probes
data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation > 0.2 &
           data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation> 0.2 &
           data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation> 0.2, ] %>% dim() 
library(UpSetR)
listInput_hyper <- list(DiffMean_G_CIMP_low_Classic_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_Codel_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation > 0.4, ] %>% rownames(),
                        
                        DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation > 0.4, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 
# We had to set the threshold to 0.2 otherwise we couldn't find enough probes (~ 90 probes as in the other comparisons we've done)


# Hypo methylated probes
data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation < -0.2 &
           data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation < -0.2  &
           data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation < -0.2 , ] %>% dim() # 69 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_G_CIMP_low_Classic_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Classic_GcimpLowOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_Codel_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Codel_GcimpLowOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_G_CIMP_high_GcimpLowOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_LGm6_GBM_GcimpLowOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation = data_sig[data_sig$DiffMean_G_CIMP_low_Mesenchymal_like_GcimpLowOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5)

# Extract hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #684 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #684 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

hyper_probes <- c(x1[ rowSums(x$New_data) == 5],
                  x1[ rowSums(x$New_data) == 4])

# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #1321 
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #1321

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # all probes present it only 5/5 comparisons
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

hypo_probes <- c(x1[ rowSums(x$New_data) == 5],
                  x1[ rowSums(x$New_data) == 4])

length(intersect(hyper_probes, hypo_probes)) # zero. Corret.
GcimpLow_probeset_22p <- c(hyper_probes, hypo_probes)
# hyper probes; diffmean > 0.4; pvalue < 0.05; 4/5 comparisons had this probeset
# hypo probes; diffmean < -0.3; pvalue < 0.05; 5/5 and 4/5 comparisons had this probeset





# LGm6_GBM Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(LGm6_GBM) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
data$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]



#**Alert**#
dim(data) #23493   
dim(data_sig) #22180   


# Hyper methylated probes
data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation > 0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation> 0.2 &
           data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation > 0.2, ] %>% dim() # 0 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation > 0.3, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 


# Hypo methylated probes
data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation < -0.2 &
           data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation < -0.2, ] %>% dim() 

library(UpSetR)
listInput_hypo <- list(DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Classic_LGm6GBMOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Codel_LGm6GBMOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_high_LGm6GBMOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_G_CIMP_low_LGm6GBMOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation = data_sig[data_sig$DiffMean_LGm6_GBM_Mesenchymal_like_LGm6GBMOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) 


# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #2661 probes somehow overlapped
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #2661 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # all probes present it only 5/5 comparisons
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret.

LGm6_probeset_45p <- c(x1[ rowSums(x$New_data) == 5],
                        x1[ rowSums(x$New_data) == 4])
# hypo probes; diffmean < -0.4; pvalue < 0.05; 5/5 and 4/5 comparisons had this probeset



# Mesenchymal_like Fold Change from all mean comparisons  --------------
data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Classic_like)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_Classic_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(Codel)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_Codel_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation <- data$meanM1 - data$meanM2

data$meanM1 <- apply(data[,as.character(Mesenchymal_like) ],1,mean,na.rm=T)  #group 1 
data$meanM2 <- apply(data[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
data$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation <- data$meanM1 - data$meanM2

data_sig <- data[data$p_val_adj < 0.05, ]


#**Alert**#
dim(data) #23493   
dim(data_sig) #22180   

# Hyper methylated probes
data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation> 0.3 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation > 0.3 &
           data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation, ] %>% dim() # 0 probes
library(UpSetR)
listInput_hyper <- list(DiffMean_Mesenchymal_Classic_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_Codel_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation > 0.2, ] %>% rownames(),
                        
                        DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation > 0.2, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 

# Hypo methylated probes
data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation < -0.2 &
           data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation, ] %>% dim() # 19 probes

library(UpSetR)
listInput_hypo <- list(DiffMean_Mesenchymal_Classic_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Classic_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_Codel_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_Codel_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_high_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_G_CIMP_low_MesenchymaOrientation < -0.2, ] %>% rownames(),
                       
                       DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation = data_sig[data_sig$DiffMean_Mesenchymal_LGm6_GBM_MesenchymaOrientation < -0.2, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5) 



# Extract hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #327 probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #327 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

hyper_probes <- c(x1[ rowSums(x$New_data) == 4])

# Extract hypo methylated probes 
x <- upset(fromList(listInput_hypo), nsets = 5)
x$New_data[1:5, 1:5]
dim(x$New_data) #4370 probes somehow overlapped
x1 <- unlist(listInput_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ]
length(x1) #4370 probes somehow overlapped

# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # none
x1[ rowSums(x$New_data) == 4] # all probes present it only 4/5 comparisons

# length(intersect(x1[ rowSums(x$New_data) == 5] , x1[ rowSums(x$New_data) == 4])) # zero. Corret. - not necessary 

hypo_probes <- c(x1[ rowSums(x$New_data) == 4])

Mesenchymal_probeset_47p <- c(hyper_probes, hypo_probes)
# hyoer probes; diffmean > 0.2; pvalue < 0.05; 4/5 comparisons had this probeset
# hypo probes; diffmean < -0.2; pvalue < 0.05; 4/5 comparisons had this probeset




DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list <- list(Classic_probeset_36p = Classic_probeset_36p,
                                                                         Codel_probeset_10p = Codel_probeset_10p,
                                                                         GcimpHigh_probeset_18p= GcimpHigh_probeset_18p,
                                                                         GcimpLow_probeset_22p = GcimpLow_probeset_22p,
                                                                         LGm6_probeset_45p = LGm6_probeset_45p,
                                                                         Mesenchymal_probeset_47p = Mesenchymal_probeset_47p)

saveRDS(DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DMP_KNN_StringestThresholdanova_diffmean_subtypes_probeset_list.rds')



tune_probset <- c(DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[1]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[2]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[3]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[4]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[5]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[6]])

tune_probset_unique <- unique(tune_probset) # we dont need to use that 'upset strategy' to get only unique probes. We just need to remove duplicates like this.


# Heatmap visualization ---------------

my_sample_col <- data.frame(row.names = metadata$Case, Subtypes_six = metadata$Ceccarelli_six_subtypes, Subtypes_seven = metadata$Supervised.DNA.Methylation.Cluster, DNAmet_cluster = metadata$IDH.specific.DNA.Methylation.Cluster) # patient ID in row names; any columns are for sample information

DNAmtx <- data_sig[, colnames(data_sig) %in% metadata$Case] #keep only beta-values

library(ComplexHeatmap)
library(matlab)
#label 'my_sample_col'
#                    Subtypes_six   Subtypes_seven DNAmet_cluster
# TCGA-CS-4938      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4941 Mesenchymal-like Mesenchymal-like       IDHwt-K2
# TCGA-CS-4942      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4943      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-4944      G-CIMP-high      G-CIMP-high      IDHmut-K2
# TCGA-CS-5390            Codel            Codel      IDHmut-K3
# Current label

### Pallete 
library(viridis)
n_colors <- length(table(my_sample_col$Subtypes_six))
Subtypes_six_pal <- viridis(n = n_colors, option = "plasma", direction = -1)
n_colors <- length(table(my_sample_col$Subtypes_seven))
Subtypes_seven_pal <- viridis(n = n_colors, option = "inferno", direction = -1)
n_colors <- length(table(my_sample_col$DNAmet_cluster))
DNAmet_cluster_pal <- viridis(n = n_colors, option = "cividis", direction = -1)

Subtypes_six_lvs <- names(table(my_sample_col$Subtypes_six))
Subtypes_seven_lvs <- names(table(my_sample_col$Subtypes_seven))
DNAmet_cluster_lvs <- names(table(my_sample_col$DNAmet_cluster))

# Name yout columns and it's levels
column_name_1 <- "Subtypes_six"
levels_1 <- Subtypes_six_lvs

column_name_2 <- "Subtypes_seven"
levels_2 <- Subtypes_seven_lvs

column_name_3 <- "DNAmet_cluster"
levels_3 <- DNAmet_cluster_lvs

# Generate color vector
colors_1 <- Subtypes_six_pal
colors_2 <- Subtypes_seven_pal
colors_3 <- DNAmet_cluster_pal
# Assign names to color vector
names(colors_1) <- levels_1
names(colors_2) <- levels_2
names(colors_3) <- levels_3
# Create named list
color_mapping <- list(column_name_1 = colors_1, 
                      column_name_2 = colors_2,
                      column_name_3 = colors_3)


top.anno = HeatmapAnnotation(df = my_sample_col, 
                             col= color_mapping,
                             show_annotation_name = T, annotation_name_gp = gpar(fontsize=7),
                             na_col= "white")



hm_unique_probset <- Heatmap(as.matrix(DNAmtx[rownames(DNAmtx) %in% tune_probset_unique, ]),
                         cluster_columns=T,
                         cluster_rows = T,
                         clustering_method_rows="complete",  
                         show_row_names = F,
                         row_names_gp = gpar(fontsize = 7),
                         show_column_names = F,
                         name = "CpG probes methylation",
                         col= jet.colors(75),
                         row_title = "CpG probes",
                         #row_names_gp = gpar(fontsize = 12),
                         #column_title = "Glioma subtypes tuned probes",
                         #column_title = "Glioma subtypes unique tuned probes",
                         #column_title = "Glioma subtypes unique probset (knn + stringest DiffMean)",
                         column_title = "Glioma subtypes unique probset",
                         #split=label.EMT.genes$E_ou_M,
                         #row_names_side="left",
                         
                         top_annotation = top.anno,
                         #row_names_gp = gpar(fontsize = 8),
                         #column_names_gp = gpar(fontsize = 12),
                         heatmap_legend_param = list(
                           color_bar = 'continuous',
                           legend_direction = 'vertical',
                           legend_width = unit(12, 'cm'),
                           legend_height = unit(10, 'cm'),
                           title_position = 'leftcenter-rot',
                           title_gp=gpar(fontsize = 16, fontface = 'bold'),
                           labels_gp=gpar(fontsize = 16, fontface = 'bold')))

print(hm_unique_probset)

### For left side annotation
# RowAnn <- HeatmapAnnotation(df = label, col=list("Epi.Mes.CellC" = c("Epi" = "blue", "Mes" = "black", "c.Cycle" = "cadetblue1")),
#                             show_annotation_name = T, annotation_name_gp = gpar(fontsize=7),
#                             na_col= "white",which = "row",show_legend =T)


### For left side annotation
# draw(hm_EMT.cCycle.MALTA + RowAnn) 

pdf("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/knn_stringestDiffmean_unique_probeset.pdf",width = 5, height = 10)
### For left side annotation
# draw(hm_EMT.cCycle.MALTA, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()




# Silhouette coef. for a clustering comparison metric -----------
library(cluster)
library(factoextra)
# Didn't worked out ...


# Cifarp Isabela abstract objective information ---------------
# Probeset intersection 
pan_glioma_probes <- read.xlsx("/media/hd/maycon/Glioma_classifier/PanGlioma_MethylationSignatures.xlsx", sheet = 1) 
pan_glioma_probes <- pan_glioma_probes$`1,300.pan-glioma.tumor.specific.probes.(Figure.2A)`

DMP_KNN_StringestThresholdanova_diffmean_subtypes_probeset_list <- readRDS("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DMP_KNN_StringestThresholdanova_diffmean_subtypes_probeset_list.rds")
tune_probset <- c(DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[1]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[2]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[3]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[4]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[5]],
                  DMP_knn_stringestThreshold_anova_diffmean_subtypes_probeset_list[[6]])

tune_probset_unique <- unique(tune_probset)


length(pan_glioma_probes) #1300
length(tune_probset_unique) #143
length(intersect(pan_glioma_probes,
                 tune_probset_unique)) #63 


# Catching how many gliomas IDH mut have been miss clustered next to the IDHwt in our probset *143probes

my_sample_col$index <- 1:dim(my_sample_col)[1]
my_sample_col$Case <- rownames(my_sample_col)
rownames(my_sample_col) <- my_sample_col$index
my_sample_col[column_order(hm_unique_probset), ][700:932, ] #part of the dendogram on heatmap which I've seen gliomas IDHmut next to IDHwt
# 12 G-CIMP-high / 249 (4.8%)
# 12 codel / 174 (6.8 %)





### Install ChAMP ### - damn it, lots of dependencies laying in the way
#** it installed ChAMP on my local machine !! **
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("minfi")
# BiocManager::install("minfi", version = "1.34.0")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocGenerics")

library(BiocGenerics)
packageVersion('BiocGenerics')
packageVersion('BiocManager')


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.17")
library(BiocManager)



# Check package version - is it the latest one?
packageVersion("BiocGenerics") #‘0.46.0’
BiocManager::install("Biobase", version =)


###########################################################################
# dmpFinder() function - I'm working on it in my local machine
# https://www.bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html#identifying-dmrs-and-dmps
###########################################################################

# Geting probes contained in both EPIC and 450k arrays (it an object from my master's codes) ----------
load('/media/hd/maycon/MESTRADO/projeto.principal/Matrix_metilacao/data.complete.sub.GSC.Rda')
dim(data.complete.sub.GSC) #453093    153
df_probes_epic_450k <- data.frame(rownames(data.complete.sub.GSC))
names(df_probes_epic_450k) <- 'epic_450k_probeID'

saveRDS(df_probes_epic_450k,"/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/df_probes_epic_450k.rds")



# Saving processed TCGA gliomas DNAmet samples within most variable feature ----------
# I processed these files more in my masters because the processed files from TCGA were with a lot of NAs on their probes 
load("/media/hd/maycon/MESTRADO/SAVE.point/matrix.methy_pdata_global/dataset.predicao.25.08.2021/pData.prediction_4801_atualizado_30.11.2021.Rda")
dim(pData.prediction_4801) #4801   15

# load matrix meth
load("/media/hd/maycon/MESTRADO/SAVE.point/matrix.methy_pdata_global/dataset.predicao.25.08.2021/matrix.meth.global_tcgaNormalized.TCGA_ID_recovered.24.01.2022.GBM.tcga.CORRETO.Rda")
dim(matrix.meth.global_tcgaNormalized_4801) #452832   4801

glioma_TCGA_proc <- pData.prediction_4801[pData.prediction_4801$Class.added %in% c('LGG_tcga.proc',
                                                               'GBM_tcga.proc'), ]
meta_glioma_TCGA_proc <- glioma_TCGA_proc[, c('barcode', 'Class.added')]

DNAmtx_glioma_TCGA_proc <- matrix.meth.global_tcgaNormalized_4801[, colnames(matrix.meth.global_tcgaNormalized_4801) %in% glioma_TCGA_proc$barcode]
dim(DNAmtx_glioma_TCGA_proc) #452832    689

### Most variable feature 
beta_values <- DNAmtx_glioma_TCGA_proc
# Removing mask/chr probes
load("/media/hd/maycon/Glioma_classifier/hm450.anno.Rda")
probes_retain <- subset(hm450.anno, !chrm_A %in% c("chrX","chrY", "chrM") & MASK_general == FALSE)$probeID
beta_values = beta_values[rownames(beta_values) %in% probes_retain, ]
# Replacing NAs
library(impute)
beta_values <- impute.knn(as.matrix(beta_values), k = 10, rowmax = 0.8, colmax = 0.8, maxp = 1500, rng.seed=362436069)[[1]] # it takes some time (like 5 min)
saveRDS(beta_values, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DNAmtx_glioma_TCGA_proc_NAreplacedKNN.rds')

# Calculate the variance of beta values for each CpG site
dim(beta_values) #383940    689
variance_values <- apply(beta_values, 1, var)

# Create a data frame with CpG sites and their variance values
variance_data <- data.frame(CpG_Site = rownames(beta_values), Variance = variance_values)

# Visualize all probes on the elbow plot
library(ggplot2); theme_set(theme_classic())
elbow_plot <- variance_data %>%
  arrange(desc(Variance)) %>%
  mutate(Rank = row_number()) %>%
  ggplot(aes(x = Rank, y = Variance)) +
  geom_line() +
  geom_point() +
  labs(title = "Elbow Plot for Variance Cutoff",
       x = "Number of CpG Sites",
       y = "Variance") +
  theme_minimal(); elbow_plot # 0.125 seems a good cutoff 

variance_data[variance_data$Variance >=  0.125, ] %>% dim() #310 probes 
variance_data[variance_data$Variance >=  0.100, ] %>% dim() #1821    probes
variance_data[variance_data$Variance >=  0.050, ] %>% dim() #28934     probes. It's a decent amount of probes to go through a Differential DNAmet

# Select the top N variable features (e.g., top 100)
library(dplyr)
top_n <- 28934
selected_features <- variance_data %>%
  arrange(desc(Variance)) %>%
  slice(1:top_n)

dim(selected_features)
head(selected_features)
names(selected_features) <- c('probeID', 'Variance')
beta_values <- beta_values[rownames(beta_values) %in% rownames(selected_features), ]
beta_values <- as.data.frame(beta_values)
beta_values$probeID <- rownames(beta_values)
DNAmtx <- merge(beta_values, selected_features, by.y = 'probeID')
DNAmtx[1:4, 1:4]

saveRDS(DNAmtx, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DNAmtx_glioma_TCGA_proc_NAreplaced_MostVariableFeature.rds')



# PCA - Most Variable Feature --------------
# It has been done on my local machine
# Plot only PCA with whole array (~ 380k probes, mask/chr/NAs probes solved)
DNAmtx_whoArray <- readRDS('/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DNAmtx_glioma_TCGA_proc_NAreplacedKNN.rds')
dim(DNAmtx_whoArray) #383940    689
library(stringr)
colnames(DNAmtx_whoArray) <- str_split_fixed(as.character(colnames(DNAmtx_whoArray)), "[-][0-9][0-9][A-Z]", 2)[,1]

load("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/from_local_mc/metadata.rda")
dim(metadata) # 878  52



# Checking on duplicates 
table(duplicated(colnames(DNAmtx_whoArray)))
# FALSE  TRUE 
#   658    31
DNAmtx_whoArray <- DNAmtx_whoArray[, !duplicated(colnames(DNAmtx_whoArray))]
table(duplicated(colnames(DNAmtx_whoArray))) # no duplcates anymore
# Adjusting samples on both metada and DNAmtx_whoArray data
metadata <- metadata[metadata$Case %in% colnames(DNAmtx_whoArray), ]
length(metadata$Case) #641
DNAmtx_whoArray <- DNAmtx_whoArray[, colnames(DNAmtx_whoArray) %in% metadata$Case]
length(colnames(DNAmtx_whoArray)) #641



pca_data <- DNAmtx_whoArray
dim(pca_data) # 383940   641
table(is.na(pca_data))
# FALSE      
# 264534660

pca2 <- prcomp(t(pca_data)) # this is the part where you calculate the PCA
aux <- as.data.frame(pca2$x[, 1:3]) #get only the info you need
pca_metadata <- merge(metadata, aux, by.y=0, by.x="Case", all.x=T)
saveRDS(pca_metadata, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/glioma_TCGA_wholearray_EPIC450k_PCA_out.rds')

library(viridis)
n_colors <- length(table(pca_metadata$Ceccarelli_six_subtypes))
pal <- viridis(n = n_colors, option = "K", direction = -1)


library(ggplot2); theme_set(theme_classic())
ggplot(pca_metadata, aes(x=PC1, y=PC2, colour=Ceccarelli_six_subtypes)) +
  geom_point() +
  scale_color_manual(values=c("orange","green", "purple", "red","darkgreen","darkblue"), name="Legend") +
  #scale_fill_manual(values=c(pal), name="Legend") + 
  #scale_color_manual(values=c(pal), name="Legend") +
  #scale_color_viridis(discrete=TRUE) +
  xlab(paste0("PC1 (",prettyNum(summary(pca2)$importance[2,1]*100, digits = 2),"%)")) +
  ylab(paste0("PC2 (",prettyNum(summary(pca2)$importance[2,2]*100, digits = 2 ),"%)")) +
  ggtitle("Gliomas TCGA - EPIC/450k genome wide") +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
  ) 

ggsave(file='/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/glioma_TCGA_wholearray_EPIC450k_PCA_out.pdf',width = 12, height = 6, dpi = 300, 
) 


# t-SNE 
library(Rtsne)
tsne_realData <- Rtsne(t(pca_data), perplexity=30, check_duplicates = FALSE) # #function to run t-sne
### Plotando o tnse 
# add as 2 colunas (V1 e V2) com os dados do tsne
pdata.teste.tsne <- metadata #pData inicial
pdata.teste.tsne$V1 <- tsne_realData$Y[,1]
pdata.teste.tsne$V2 <- tsne_realData$Y[,2]

library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne[, ], aes(x=V1, V2, colour=Ceccarelli_six_subtypes)) +
  geom_point() +
  scale_color_manual(values=c("orange","green", "purple", "red","darkgreen","darkblue"), name="Legend") +
  # scale_color_viridis(discrete=TRUE) +
  
  ggtitle(label = "t-SNE Genome wide (380k probes)",
          subtitle = "perplexity=30")

ggsave("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/glioma_TCGA_wholearray_EPIC450k_tSNE_out.pdf",width = 6, height = 4)


library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne[, ], aes(x=V1, V2, colour=IDH.status)) +
  geom_point() +
  scale_color_manual(values=c("orange","green", "purple", "red","darkgreen","darkblue"), name="Legend") +
  # scale_color_viridis(discrete=TRUE) +
  
  ggtitle(label = "t-SNE Genome wide (380k probes)",
          subtitle = "perplexity=30")






###############################################################
###** Fine tune the MostVariableFeature (MVF) by DiffMean **###

library(dplyr)
library(stringr)
library(UpSetR)

# 1. load DNAmtx with MVF
DNAmtx <- readRDS("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DNAmtx_glioma_TCGA_proc_NAreplaced_MostVariableFeature.rds")
dim(DNAmtx) # 28934   691
DNAmtx %>% dim() #28934   691
rownames(DNAmtx) <- DNAmtx$probeID
DNAmtx$probeID <- NULL
DNAmtx$Variance <- NULL
DNAmtx <- as.matrix(DNAmtx)
colnames(DNAmtx) <- str_split_fixed(as.character(colnames(DNAmtx)), "[-][0-9][0-9][A-Z]", 2)[,1] #to match to metadata sample ID


# 2. laod metadata
load("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/from_local_mc/metadata.rda")
dim(metadata) #878  52
metadata[1:4, 1:4]


# 3. correct samples intersected in both metadata and DNAmtx
table(duplicated(colnames(DNAmtx)))
# FALSE  TRUE 
#   658    31
DNAmtx <- DNAmtx[, !duplicated(colnames(DNAmtx))]
table(duplicated(colnames(DNAmtx))) # no duplcates anymore
# Adjusting samples on both metada and DNAmtx data
metadata <- metadata[metadata$Case %in% colnames(DNAmtx), ]
length(metadata$Case) #641
DNAmtx <- DNAmtx[, colnames(DNAmtx) %in% metadata$Case]
length(colnames(DNAmtx)) #641

#** NOTE **#
# All of these MVF seems to be significant different accordingly to dmpFinder() - from {minfi} - & this set of probes made a great clustering on t-SNE over the 6 subtypes. That's why we're moving foward without a "Differential DNA methylation" 


# 4. Find DiffMean in a group by group way (as we've done before)

metadata$Ceccarelli_six_subtypes %>% table() # there're fewer samples because now we're working only with 450k samples

# Classic-like    Codel            G-CIMP-high       G-CIMP-low 
# 76              173              240               12 
# LGm6-GBM        Mesenchymal-like 
# 39              101 

Classic_like <-  metadata[metadata$Ceccarelli_six_subtypes %in% 'Classic-like', ]$Case
Codel <-  metadata[metadata$Ceccarelli_six_subtypes %in% 'Codel', ]$Case
G_CIMP_high <- metadata[metadata$Ceccarelli_six_subtypes %in% 'G-CIMP-high', ]$Case
G_CIMP_low <- metadata[metadata$Ceccarelli_six_subtypes %in% 'G-CIMP-low', ]$Case
LGm6_GBM <- metadata[metadata$Ceccarelli_six_subtypes %in% 'LGm6-GBM', ]$Case
Mesenchymal_like <- metadata[metadata$Ceccarelli_six_subtypes %in% 'Mesenchymal-like', ]$Case

# Check NAs in DNAmtx
table(is.na(DNAmtx))
# FALSE 
# 18546694

# Turn DNAmtx to a "data frame" 
DNAmtx <- DNAmtx %>% as.data.frame()


# Classic_like Fold Change from all mean comparisons  --------------

# Calculate DiffMean 
DNAmtx$meanM1 <- apply(DNAmtx[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
DNAmtx$meanM2 <- apply(DNAmtx[,as.character(Codel)],1,mean,na.rm=T) #group n 
DNAmtx$DiffMean_Classic_Codel_ClassicOrientation <- DNAmtx$meanM1 - DNAmtx$meanM2

DNAmtx$meanM1 <- apply(DNAmtx[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
DNAmtx$meanM2 <- apply(DNAmtx[,as.character(G_CIMP_high)],1,mean,na.rm=T) #group n 
DNAmtx$DiffMean_Classic_G_CIMP_high_ClassicOrientation <- DNAmtx$meanM1 - DNAmtx$meanM2

DNAmtx$meanM1 <- apply(DNAmtx[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
DNAmtx$meanM2 <- apply(DNAmtx[,as.character(G_CIMP_low)],1,mean,na.rm=T) #group n 
DNAmtx$DiffMean_Classic_G_CIMP_low_ClassicOrientation <- DNAmtx$meanM1 - DNAmtx$meanM2

DNAmtx$meanM1 <- apply(DNAmtx[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
DNAmtx$meanM2 <- apply(DNAmtx[,as.character(LGm6_GBM)],1,mean,na.rm=T) #group n 
DNAmtx$DiffMean_Classic_LGm6_GBM_ClassicOrientation <- DNAmtx$meanM1 - DNAmtx$meanM2

DNAmtx$meanM1 <- apply(DNAmtx[,as.character(Classic_like) ],1,mean,na.rm=T)  #group 1 
DNAmtx$meanM2 <- apply(DNAmtx[,as.character(Mesenchymal_like)],1,mean,na.rm=T) #group n 
DNAmtx$DiffMean_Classic_Mesenchymal_like_ClassicOrientation <- DNAmtx$meanM1 - DNAmtx$meanM2


# Check on probes intersection by DNAmet orientation
# Hyper methylated probes
library(UpSetR)
listInput_hyper <- list(DiffMean_Classic_Codel_ClassicOrientation = DNAmtx[DNAmtx$DiffMean_Classic_Codel_ClassicOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Classic_G_CIMP_high_ClassicOrientation = DNAmtx[DNAmtx$DiffMean_Classic_G_CIMP_high_ClassicOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Classic_G_CIMP_low_ClassicOrientation = DNAmtx[DNAmtx$DiffMean_Classic_G_CIMP_low_ClassicOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Classic_LGm6_GBM_ClassicOrientation = DNAmtx[DNAmtx$DiffMean_Classic_LGm6_GBM_ClassicOrientation > 0.3, ] %>% rownames(),
                        
                        DiffMean_Classic_Mesenchymal_like_ClassicOrientation = DNAmtx[DNAmtx$DiffMean_Classic_Mesenchymal_like_ClassicOrientation > 0.3, ] %>% rownames()
)

upset(fromList(listInput_hyper), order.by = "freq", nsets = 5) 


# Hypo methylated probes
library(UpSetR)
listInput_hypo <- list(DiffMean_Classic_Codel_ClassicOrientation = DNAmtx[DNAmtx$DiffMean_Classic_Codel_ClassicOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Classic_G_CIMP_high_ClassicOrientation = DNAmtx[DNAmtx$DiffMean_Classic_G_CIMP_high_ClassicOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Classic_G_CIMP_low_ClassicOrientation = DNAmtx[DNAmtx$DiffMean_Classic_G_CIMP_low_ClassicOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Classic_LGm6_GBM_ClassicOrientation = DNAmtx[DNAmtx$DiffMean_Classic_LGm6_GBM_ClassicOrientation < -0.3, ] %>% rownames(),
                       
                       DiffMean_Classic_Mesenchymal_like_ClassicOrientation = DNAmtx[DNAmtx$DiffMean_Classic_Mesenchymal_like_ClassicOrientation < -0.3, ] %>% rownames()
)

upset(fromList(listInput_hypo), order.by = "freq", nsets = 5)  

# Extract only probes that went specific for the group been analyzed compared to ALL THE OTHER 5 GROUPS
# In this case, only hyper methylated probes (|0.3|) were suited  

# Extract the hyper methylated probes 
x <- upset(fromList(listInput_hyper), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. Sometimes there're more than one probeset intersected across "4" groups. You gonna get all those probes in every probeset wich contains 4 groups sharing that probeset
length(x1[ rowSums(x$New_data) == 5])

Classic_probeset_69p <- x1[ rowSums(x$New_data) == 5]





Classic_like <-  metadata[metadata$Ceccarelli_six_subtypes %in% 'Classic-like', ]$Case
Codel <-  metadata[metadata$Ceccarelli_six_subtypes %in% 'Codel', ]$Case
G_CIMP_high <- metadata[metadata$Ceccarelli_six_subtypes %in% 'G-CIMP-high', ]$Case
G_CIMP_low <- metadata[metadata$Ceccarelli_six_subtypes %in% 'G-CIMP-low', ]$Case
LGm6_GBM <- metadata[metadata$Ceccarelli_six_subtypes %in% 'LGm6-GBM', ]$Case
Mesenchymal_like <- metadata[metadata$Ceccarelli_six_subtypes %in% 'Mesenchymal-like', ]$Case

# Automating the DiffMean calculation -----------------------
#DNAmtx <- DNAmtx[, -(642:648)] just deleting the last diffmean try
# Testing the loop
i <- 1
condition_one <- 'Classic-like'
data <- DNAmtx
metadata <- metadata
threshold <- 0.3

group_names <- names(table(metadata$Ceccarelli_six_subtypes))
calculate_diff_mean <- function(data, 
                                metadata, 
                                condition_one, 
                                threshold) {
  list_diffmean_dfs <- list()
  data <- as.data.frame(data)
  condition_all_but_one <- group_names[!group_names %in% condition_one]
  ###condition_all_but_one <- 'Mesenchymal-like' #TESTING only. 
  for(i in 1:length(condition_all_but_one)) {
  print(paste0("Processing group: ", condition_all_but_one[i]))
    
  try({ #it ignores any error that might prevent the code of going forward
        #I did it because of the absence of hyper or hypo probes were making the loop stop
  loop_data <- data
  condition_1_sampleID <- metadata[metadata$Ceccarelli_six_subtypes %in% condition_one, ]$Case
  condition_2_sampleID <- metadata[metadata$Ceccarelli_six_subtypes %in% condition_all_but_one[i], ]$Case
 
  print(paste0("Length of condition_1_sampleID: ", length(condition_1_sampleID)))
  print(paste0("Length of condition_2_sampleID: ", length(condition_2_sampleID)))
  
  loop_data$meanM1 <- apply(loop_data[, condition_1_sampleID], 1, mean, na.rm = TRUE)
  loop_data$meanM2 <- apply(loop_data[, condition_2_sampleID], 1, mean, na.rm = TRUE)
  loop_data$DiffMean <- loop_data$meanM1 - loop_data$meanM2
  loop_data$Comparison <- paste0(condition_one, '_', condition_all_but_one[i], '_', condition_one, '_', 'Orientation')
  

  
  loop_data$DNAmet_orientation <- NA
  #loop_data$DNAmet_orientation <- as.character(loop_data$DNAmet_orientation) #this time it's not necessary
  if(dim(loop_data[loop_data$DiffMean > threshold, ])[1] > 0) {
    loop_data[loop_data$DiffMean > threshold, ]$DNAmet_orientation <- 'hyper'
  } else {
    # do nothing
    }
  
  if(dim(loop_data[loop_data$DiffMean < -threshold, ])[1] > 0) {
    loop_data[loop_data$DiffMean < -threshold, ]$DNAmet_orientation <- 'hypo'
  } else {
    # do nothing
  }

  loop_data[loop_data$DNAmet_orientation %in% NA, ]$DNAmet_orientation <- 'not_diff'
  
  
  
  loop_data$probeID <- rownames(loop_data)
  
  list_diffmean_dfs[[i]] <- loop_data[, c('DiffMean', 'Comparison', 'DNAmet_orientation', 'probeID')]
  print(paste0(list_diffmean_dfs[[i]]$Comparison[1], ' has been stored into the list.'))
  },  silent = FALSE)
  }
  return(list_diffmean_dfs) 
  }


# Classic-like DNAmet markers ---------------
list_diffmean_dfs <-calculate_diff_mean(#data = DNAmtx,
                                        data = DNAmtx_merged,
                    metadata = metadata,
                    condition_one = group_names[1], #"Classic-like" 
                    threshold = 0.3)

Classiclike_markers <- do.call('rbind', list_diffmean_dfs)
head(Classiclike_markers)
length(names(table(Classiclike_markers$Comparison))) # 5 (it should be 5) 
# Attention 2: it's okay to have more probes than the ~28k MVF because in this dataframe should be 5 comparisons. So it's the sum of 5*~28k probes at total.
to_upsetplot <- Classiclike_markers[Classiclike_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  Classic_Codel = to_upsetplot[to_upsetplot$Comparison %in% 'Classic-like_Codel_Classic-like_Orientation', ]$probeID,
  Classic_GCIMPhigh = to_upsetplot[to_upsetplot$Comparison %in% 'Classic-like_G-CIMP-high_Classic-like_Orientation', ]$probeID,
  Classic_GCIMPlow = to_upsetplot[to_upsetplot$Comparison %in% 'Classic-like_G-CIMP-low_Classic-like_Orientation', ]$probeID,
  Classic_LGm6GBM = to_upsetplot[to_upsetplot$Comparison %in% 'Classic-like_LGm6-GBM_Classic-like_Orientation', ]$probeID,
  Classic_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'Classic-like_Mesenchymal-like_Classic-like_Orientation', ]$probeID
)

upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

#Classic_probeset_70p <- x1[ rowSums(x$New_data) == 5]
Classic_probeset_115p <- x1[ rowSums(x$New_data) == 5]


# Codel DNAmet markers ---------------
list_diffmean_dfs <- calculate_diff_mean(# data = DNAmtx,
                                         data = DNAmtx_merged,
                    metadata = metadata,
                    condition_one = group_names[2], #"Codel" 
                    threshold = 0.3)

Codel_markers <- do.call('rbind', list_diffmean_dfs)
head(Codel_markers)
length(names(table(Codel_markers$Comparison))) # 5 (it should be 5)

to_upsetplot <- Codel_markers[Codel_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  Codel_Classic = to_upsetplot[to_upsetplot$Comparison %in% 'Codel_Classic-like_Codel_Orientation', ]$probeID,
  Codel_GCIMPhigh = to_upsetplot[to_upsetplot$Comparison %in% 'Codel_G-CIMP-high_Codel_Orientation', ]$probeID,
  Codel_GCIMPlow = to_upsetplot[to_upsetplot$Comparison %in% 'Codel_G-CIMP-low_Codel_Orientation', ]$probeID,
  Codel_LGm6GBM = to_upsetplot[to_upsetplot$Comparison %in% 'Codel_LGm6-GBM_Codel_Orientation', ]$probeID,
  Codel_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'Codel_Mesenchymal-like_Codel_Orientation', ]$probeID
)

upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

#Codel_probeset_147p <- x1[ rowSums(x$New_data) == 5]
Codel_probeset_147p <- x1[ rowSums(x$New_data) == 5]




# G-CIMP-high DNAmet markers ---------------
list_diffmean_dfs <- calculate_diff_mean(#data = DNAmtx,
                                         data = DNAmtx_merged,
                    metadata = metadata,
                    condition_one = group_names[3], #"G-CIMP-high" 
                    threshold = 0.3)

GCIMPhigh_markers <- do.call('rbind', list_diffmean_dfs)
head(GCIMPhigh_markers)
length(names(table(GCIMPhigh_markers$Comparison))) # 5 (it should be 5)

to_upsetplot <- GCIMPhigh_markers[GCIMPhigh_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  GCIMPhigh_Classic = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-high_Classic-like_G-CIMP-high_Orientation', ]$probeID,
  GCIMPhigh_Codel = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-high_Codel_G-CIMP-high_Orientation', ]$probeID,
  GCIMPhigh_GCIMPlow = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-high_G-CIMP-low_G-CIMP-high_Orientation', ]$probeID,
  GCIMPhigh_LGm6GBM = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-high_LGm6-GBM_G-CIMP-high_Orientation', ]$probeID,
  GCIMPhigh_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-high_Mesenchymal-like_G-CIMP-high_Orientation', ]$probeID
)

upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

#CGIMPhigh_probeset_13p <- x1[ rowSums(x$New_data) == 5]
CGIMPhigh_probeset_13p <- x1[ rowSums(x$New_data) == 5]


# G-CIMP-low DNAmet markers ---------------
list_diffmean_dfs <- calculate_diff_mean(#data = DNAmtx,
                                         data = DNAmtx_merged,
                    metadata = metadata,
                    condition_one = group_names[4], #"G-CIMP-low" 
                    threshold = 0.3)

GCIMPlow_markers <- do.call('rbind', list_diffmean_dfs)
head(GCIMPlow_markers)
length(names(table(GCIMPlow_markers$Comparison))) # 5 (it should be 5)

to_upsetplot <- GCIMPlow_markers[GCIMPlow_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  GCIMPlow_Classic = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-low_Classic-like_G-CIMP-low_Orientation', ]$probeID,
  GCIMPlow_Codel = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-low_Codel_G-CIMP-low_Orientation', ]$probeID,
  GCIMPlow_GCIMPhigh = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-low_G-CIMP-high_G-CIMP-low_Orientation', ]$probeID,
  GCIMPlow_LGm6GBM = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-low_LGm6-GBM_G-CIMP-low_Orientation', ]$probeID,
  GCIMPlow_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'G-CIMP-low_Mesenchymal-like_G-CIMP-low_Orientation', ]$probeID
)

upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

#CGIMPlow_probeset_129p <- x1[ rowSums(x$New_data) == 5]
CGIMPlow_probeset_129p <- x1[ rowSums(x$New_data) == 5]




# LGm6-GBM DNAmet markers ---------------
list_diffmean_dfs <- calculate_diff_mean(#data = DNAmtx,
                                         data = DNAmtx_merged,
                    metadata = metadata,
                    condition_one = group_names[5], #"LGm6-GBM" 
                    threshold = 0.3)

LGm6GBM_markers <- do.call('rbind', list_diffmean_dfs)
head(LGm6GBM_markers)
length(names(table(LGm6GBM_markers$Comparison))) # 5 (it should be 5)

to_upsetplot <- LGm6GBM_markers[LGm6GBM_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  LGm6GBM_Classic = to_upsetplot[to_upsetplot$Comparison %in% 'LGm6-GBM_Classic-like_LGm6-GBM_Orientation', ]$probeID,
  LGm6GBM_Codel = to_upsetplot[to_upsetplot$Comparison %in% 'LGm6-GBM_Codel_LGm6-GBM_Orientation', ]$probeID,
  LGm6GBM_GCIMPhigh = to_upsetplot[to_upsetplot$Comparison %in% 'LGm6-GBM_G-CIMP-high_LGm6-GBM_Orientation', ]$probeID,
  LGm6GBM_LGm6GBM = to_upsetplot[to_upsetplot$Comparison %in% 'LGm6-GBM_G-CIMP-low_LGm6-GBM_Orientation', ]$probeID,
  LGm6GBM_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'LGm6-GBM_Mesenchymal-like_LGm6-GBM_Orientation', ]$probeID
)

upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

#LGm6GBM_probeset_12p <- x1[ rowSums(x$New_data) == 5]
LGm6GBM_probeset_12p <- x1[ rowSums(x$New_data) == 5]





# Mesenchymal-like DNAmet markers ---------------
list_diffmean_dfs <- calculate_diff_mean(#data = DNAmtx,
                                         data = DNAmtx_merged,
                    metadata = metadata,
                    condition_one = group_names[6], #"Mesenchymal-like" 
                    threshold = 0.3)

Mesenchymallike_markers <- do.call('rbind', list_diffmean_dfs)
head(Mesenchymallike_markers)
length(names(table(Mesenchymallike_markers$Comparison))) # 5 (it should be 5)

to_upsetplot <- Mesenchymallike_markers[Mesenchymallike_markers$DNAmet_orientation %in% c('hyper', 'hypo'), ]
table(to_upsetplot$Comparison)

listInput_hyper_hypo <- list(
  Mesenchymal_Classic = to_upsetplot[to_upsetplot$Comparison %in% 'Mesenchymal-like_Classic-like_Mesenchymal-like_Orientation', ]$probeID,
  Mesenchymal_Codel = to_upsetplot[to_upsetplot$Comparison %in% 'Mesenchymal-like_Codel_Mesenchymal-like_Orientation', ]$probeID,
  Mesenchymal_GCIMPhigh = to_upsetplot[to_upsetplot$Comparison %in% 'Mesenchymal-like_G-CIMP-high_Mesenchymal-like_Orientation', ]$probeID,
  Mesenchymal_GCIMPlow = to_upsetplot[to_upsetplot$Comparison %in% 'Mesenchymal-like_G-CIMP-low_Mesenchymal-like_Orientation', ]$probeID,
  Mesenchymal_Mesenchymal = to_upsetplot[to_upsetplot$Comparison %in% 'Mesenchymal-like_LGm6-GBM_Mesenchymal-like_Orientation', ]$probeID
)

upset(fromList(listInput_hyper_hypo), order.by = "freq", nsets = 5) 

# Extract probes 
x <- upset(fromList(listInput_hyper_hypo), nsets = 5)
x$New_data[1:5, 1:5] # dummy df to probes in each group comparison
dim(x$New_data)[1] # n of probes somehow overlapped
x1 <- unlist(listInput_hyper_hypo, use.names = FALSE)
x1 <- x1[ !duplicated(x1) ] # actual probe list 
length(x1) # n of probes somehow overlapped
# x and x1 have the their elements aligned 
x1[ rowSums(x$New_data) == 5] # == <number> recover all probes present in n intersections. 
length(x1[ rowSums(x$New_data) == 5])

#Mesenchymal_probeset_7p <- x1[ rowSums(x$New_data) == 5]
Mesenchymal_probeset_7p <- x1[ rowSums(x$New_data) == 5]

# Extreme_upset_cutoff_Glioma_probes_vector <- c(
#   Classic_probeset_70p, 
#   Codel_probeset_147p,
#   CGIMPhigh_probeset_13p, 
#   CGIMPlow_probeset_129p, 
#   LGm6GBM_probeset_12p, 
#   Mesenchymal_probeset_7p)

Extreme_upset_cutoff_Glioma_probes_vector <- c(
  Classic_probeset_115p,
  Codel_probeset_147p,
  CGIMPhigh_probeset_13p,
  CGIMPlow_probeset_129p,
  LGm6GBM_probeset_12p,
  Mesenchymal_probeset_7p)

# saveRDS(Extreme_upset_cutoff_Glioma_probes_vector, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/Extreme_upset_cutoff_Glioma_probes_vector.rds')


DNAmtx_ext_pcutoff <- DNAmtx[rownames(DNAmtx) %in% Extreme_upset_cutoff_Glioma_probes_vector, ]
dim(DNAmtx_ext_pcutoff) #371 641
dim(metadata) #641  52


# PCA 
DNAmtx_ext_pcutoff$probeID <- NULL
pca_data <- DNAmtx_ext_pcutoff
dim(pca_data) # 371 641
table(is.na(pca_data))
# FALSE      
# 237811

pca2 <- prcomp(t(pca_data)) # this is the part where you calculate the PCA
aux <- as.data.frame(pca2$x[, 1:3]) #get only the info you need
pca_metadata <- merge(metadata, aux, by.y=0, by.x="Case", all.x=T)
#saveRDS(pca_metadata, file = '/.rds')


library(ggplot2); theme_set(theme_classic())
ggplot(pca_metadata, aes(x=PC1, y=PC2, colour=Ceccarelli_six_subtypes)) +
  geom_point() +
  scale_color_manual(values=c("orange","green", "purple", "red","darkgreen","darkblue"), name="Legend") +
  #scale_fill_manual(values=c(pal), name="Legend") + 
  #scale_color_manual(values=c(pal), name="Legend") +
  #scale_color_viridis(discrete=TRUE) +
  xlab(paste0("PC1 (",prettyNum(summary(pca2)$importance[2,1]*100, digits = 2),"%)")) +
  ylab(paste0("PC2 (",prettyNum(summary(pca2)$importance[2,2]*100, digits = 2 ),"%)")) +
  ggtitle("Gliomas TCGA - EPIC/450k 371 fine tuned probes") +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)
  ) 

ggsave(file='/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/glioma_TCGA_wholearray_EPIC450k_371finetunedprobes_PCA_out.pdf',width = 12, height = 6, dpi = 300,
) 


# t-SNE 
library(Rtsne)
tsne_realData <- Rtsne(t(pca_data), perplexity=30, check_duplicates = FALSE) # #function to run t-sne
### Plotando o tnse 
# add as 2 colunas (V1 e V2) com os dados do tsne
pdata.teste.tsne <- metadata #pData inicial
pdata.teste.tsne$V1 <- tsne_realData$Y[,1]
pdata.teste.tsne$V2 <- tsne_realData$Y[,2]

library(ggplot2); theme_set(theme_classic())
ggplot(pdata.teste.tsne[, ], aes(x=V1, V2, colour=Ceccarelli_six_subtypes)) +
  geom_point() +
  scale_color_manual(values=c("orange","green", "purple", "red","darkgreen","darkblue"), name="Legend") +
  # scale_color_viridis(discrete=TRUE) +
  
  ggtitle(label = "t-SNE Gliomas TCGA - EPIC/450k 371 fine tuned probes",
          subtitle = "perplexity=30")

ggsave("/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/glioma_TCGA_wholearray_EPIC450k_371finetunedprobes_tSNE_out.pdf",width = 6, height = 4)








### Remark ### ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ### 
# Classic and Mesenchymal subtypes are too much similar 
# We must go back and get a Most Variable Feature only from these two subtypes
# Then we can add thes probes to our MVF probeset to see if we get a better DeaffMean than 0.1. 
### Remark ### ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ###  ### 

DNAmtx_MVF_Classic_Mesenchymal_Based <- readRDS('/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DNAmtx_glioma_Classic_Mesenchymal_TCGA_proc_NAreplaced_MostVariableFeature.rds') # code from /media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/MVF_classic_vs_mesemchymal.R

dim(DNAmtx) #28934   641 MVF within all the 6 subtypes
dim(DNAmtx_MVF_Classic_Mesenchymal_Based) #963 643 MVF within only Classic and Mesenchymal subtypes

# samples intersected
length(intersect(colnames(DNAmtx),
                 colnames(DNAmtx_MVF_Classic_Mesenchymal_Based))) #641

# probes intersected
length(intersect(rownames(DNAmtx),
                 rownames(DNAmtx_MVF_Classic_Mesenchymal_Based))) #855

# Removing intersected probes
remove_probes <- intersect(rownames(DNAmtx),
                           rownames(DNAmtx_MVF_Classic_Mesenchymal_Based))
DNAmtx_MVF_Classic_Mesenchymal_Based <- DNAmtx_MVF_Classic_Mesenchymal_Based[!rownames(DNAmtx_MVF_Classic_Mesenchymal_Based) %in% 
                                       remove_probes, ]

# probes intersected
length(intersect(rownames(DNAmtx),
                 rownames(DNAmtx_MVF_Classic_Mesenchymal_Based))) #0 okay

DNAmtx <- as.data.frame(DNAmtx)
DNAmtx$probeID <- rownames(DNAmtx)
DNAmtx_merged <- plyr::rbind.fill(as.data.frame(DNAmtx), 
                                  as.data.frame(DNAmtx_MVF_Classic_Mesenchymal_Based))
head(DNAmtx_merged)
dim(DNAmtx_merged)
rownames(DNAmtx_merged) <- DNAmtx_merged$probeID
DNAmtx_merged$probeID <- NULL
DNAmtx_merged$Variance <- NULL

saveRDS(DNAmtx_merged, file = '/media/hd/maycon/Glioma_classifier/just_tune_Cecca_model/DNAmtx_glioma_TCGA_proc_NAreplaced_AllMostVariableFeature_and_ClassicMesenchymalMVF.rds')


# Now, try to apply 'calculate_diff_mean' function into the DNAmtx with the new added features
# Try it on Classic and on Mesenchymal 

# Find Classic probe markers
list_diffmean_dfs <- calculate_diff_mean(data = DNAmtx_merged,
                                        metadata = metadata,
                                        condition_one = group_names[1], #"Classic-like" 
                                        threshold = 0.3)

Classiclike_markers <- do.call('rbind', list_diffmean_dfs)
head(Classiclike_markers)
length(names(table(Classiclike_markers$Comparison))) # 4 (it should be 5) 

# Find Mesenhcymal probe markers
list_diffmean_dfs <- calculate_diff_mean(data = DNAmtx_merged,
                                         metadata = metadata,
                                         condition_one = group_names[6], #"Mesenchymal-like" 
                                         threshold = 0.3)

Mesenchymallike_markers <- do.call('rbind', list_diffmean_dfs)
head(Mesenchymallike_markers)
length(names(table(Mesenchymallike_markers$Comparison))) # 4 (it should be 5)

#** Note **#
#Extracting MVF only from Classic and Mesenhcymal subtypes hasn't brought any gain to find their probes markers within threshold = 0.3 for Differential Mean Methylation











