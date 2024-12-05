# Clear environment
rm(list=ls())

pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(e1071)
library(caTools)
library(dplyr)


## Functions
norm_data <- function(x){
  (x-mean(x))/(sd(x))
}

t_test_col <- function(x,y=NA){
  
  t_test <- matrix(ncol=2, nrow=ncol(x))
  col_names_t <- c("t","p")
  
  for (t_index in 1:ncol(x)){
    col_c <- x[,t_index]
    t_test_c <- t.test(col_c, mu=0.5)
    t_test[t_index,1] <- t_test_c$statistic
    t_test[t_index,2] <- t_test_c$p.value
  }
  t_test <- data.frame(t_test)
  colnames(t_test) <- col_names_t
  if (!is.na(y)){
    rownames(t_test) <- y
  }
  return(t_test)
}


## Prepare
sub_dir <- dir('.')[file.info(dir('.',full.names=T))$isdir]

accuracy_lh <- matrix(ncol=3, nrow=length(sub_dir))
accuracy_rh <- accuracy_lh

set.seed(42)


sub_index <- 1
for (sub in sub_dir){
  
  message(sprintf("Subject %s\n", sub_index))
  
  #sub <- sub_dir[1]
  path_c <- paste(pwd, sub, sep='/')
  roi_list = list.files(path = path_c, pattern = "*02.txt")
  
  roi_index <- 1
  for (roi_sub in roi_list){
    
    #roi_sub <- roi_list[1]
    roi_c <- paste(path_c, roi_sub, sep='/')
    roi_data <- read.table(roi_c,sep=",",header=T)
    
    # 1 = Can/NonCan; 0/1 = can/non-can; 0/1 = side view, tilted/depth view, rotated; 0/1 = metal/non-metal
    # 2 = global shape; 2/8 = scrambling; 0 = filler; 0 = filler
    
    ## Data by experiment
    # Objects can-noncan
    exp_can <- roi_data[roi_data[,1]==1,]
    exp_can <- select(exp_can, -betas1, -betas3, -betas4)
    
    data_can <- exp_can[exp_can[,1]==0,]
    data_noncan <- exp_can[exp_can[,1]==1,]
    
    labels_can <- select(data_can, betas2)
    #labels_can[labels_can=="0"] <- "FALSE"
    labels_can[labels_can=="0"] <- "TRUE"
    labels_can <- as.character(labels_can$betas2)
    
    data_can <- select(data_can, -betas2)

    labels_noncan <- select(data_noncan, betas2)
    labels_noncan[labels_noncan=="1"] <- "TRUE"
    labels_noncan <- as.character(labels_noncan$betas2)
    
    data_noncan <- select(data_noncan, -betas2)
    
    # Global shapes
    exp_global <- roi_data[roi_data[,1]==2,]
    exp_global <- select(exp_global, -betas1, -betas3, -betas4)
    
    exp_global <- exp_global[exp_global[,1]==2,]
    exp_global$betas2 <- 1
    
    labels_global <- select(exp_global, betas2)
    labels_global[labels_global=="1"] <- "TRUE"
    labels_global <- as.character(labels_global$betas2)
    
    data_global <- select(exp_global, -betas2)
    
    # Split global into training and test data
    global_split <- sample.split(data_global, SplitRatio = 0.80)
    
    labels_global_train <- labels_global[global_split]
    data_global_train <- data_global[global_split,]
    
    labels_global_test <- labels_global[-global_split]
    data_global_test <- data_global[-global_split,]
    
    rm(exp_can, exp_global, roi_data)
    
    # Normalize features
    #data_can <- t(apply(data_can, 1, norm_data))
    #data_noncan <- t(apply(data_noncan, 1, norm_data))
    #data_global <- t(apply(data_global, 1, norm_data))
    
    
    ## SVM
    #svm.global <- svm(data_global,y=NULL, type='one-classification', nu = 0.10, scale=TRUE, kernel="radial", cross = 10)
    svm.global <- svm(data_global_train, y=NULL, type='one-classification', nu = 0.05, scale = TRUE, kernel="radial", cross = 10)
    
    svm.predglobal <- predict(svm.global, data_global_test)
    svm.prednoncan <- predict(svm.global, data_noncan)
    svm.predcan <- predict(svm.global, data_can)
    
    conf_global <- table(Predicted = svm.predglobal, Reference = labels_global_test)
    conf_noncan <- table(Predicted = svm.prednoncan, Reference = labels_noncan)
    conf_can <- table(Predicted = svm.predcan, Reference = labels_can)
    
    accuracy_global <- conf_global[2]/(conf_global[1] + conf_global[2])
    accuracy_noncan <- conf_noncan[2]/(conf_noncan[1] + conf_noncan[2])
    accuracy_can <- conf_can[2]/(conf_can[1] + conf_can[2])
    
    
    # Save data
    if (roi_index == 1){
      accuracy_lh[sub_index,] <- c(accuracy_global, accuracy_noncan, accuracy_can)
    } else {
      accuracy_rh[sub_index,] <- c(accuracy_global, accuracy_noncan, accuracy_can)
    }
    
    roi_index <- roi_index + 1
  }
  
  sub_index <- sub_index + 1
  
}


## Prepare results
col_names_tpj <- c("ACC_Global", "ACC_NonCan", "ACC_Can")

accuracy_lh  = data.frame(accuracy_lh )
colnames(accuracy_lh) <- col_names_tpj

accuracy_rh  = data.frame(accuracy_rh )
colnames(accuracy_rh) <- col_names_tpj

write.table(accuracy_lh, file = "SVM_OneClass_TPJ_Global_Obj_accuracy_lh.txt", sep = ",", row.names = FALSE, col.names = TRUE)
write.table(accuracy_rh, file = "SVM_OneClass_TPJ_Global_Obj_accuracy_rh.txt", sep = ",", row.names = FALSE, col.names = TRUE)
