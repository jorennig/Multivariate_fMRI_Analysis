# Clear environment
rm(list=ls())

pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(e1071)
library(caTools)
library(dplyr)
library(stringr)

## Prepare
sub_dir <- dir('.')[file.info(dir('.',full.names=T))$isdir]

set.seed(42)

accuracy = data.frame()
for (sub in sub_dir){
  #sub <- sub_dir[1]
  
  message(sprintf("Subject %s", sub))
  
  path_c <- paste(pwd, sub, sep='/')
  roi_list = list.files(path = path_c, pattern = "*.txt")
  
  for (roi_sub in roi_list){
    #roi_sub <- roi_list[1]
    
    roi_name <- str_split(roi_sub, ".txt", simplify = TRUE)[,1]
    split_name <- str_split(roi_name, "_", simplify = TRUE)
    roi <- split_name[,3]
    hem <- split_name[,2]
    hem <- str_replace(hem, 'h', '')
    
    message(sprintf("Subject: %s, ROI: %s, Hemisphere: %s", sub, roi, hem))
    
    roi_c <- paste(path_c, roi_sub, sep='/')
    roi_data <- read.table(roi_c,sep=",",header=T)
    
    if (ncol(roi_data) > 4){
      
      ## Data by experiment
      # Global (2 = global, 8 = scrambled)
      exp_global <- roi_data[roi_data[,1]==2, ]
      exp_global <- select(exp_global, -betas1, -betas3, -betas4)
      exp_tot <- exp_global

      exp_tot$betas2[exp_tot$betas2 == 2]  <- 1
      exp_tot$betas2[exp_tot$betas2 == 8]  <- 2
      
      y_exp_tot <- factor(exp_tot[,1], levels = c(1, 2))
      x_exp_tot <- select(exp_tot, -betas2)
      
      exp_tot_split <- sample.split(y_exp_tot, SplitRatio = 0.80)
      
      x_train <- subset(x_exp_tot, exp_tot_split == TRUE)
      x_test <- subset(x_exp_tot, exp_tot_split == FALSE)
      
      y_train <- subset(y_exp_tot, exp_tot_split == TRUE)
      y_test <- subset(y_exp_tot, exp_tot_split == FALSE)
      
      
      # Normalize features
      x_exp_tot <- scale(x_exp_tot) 
      x_train <- scale(x_train) 
      x_test <- scale(x_test) 
      
      
      ## SVM
      #svm_tune <- tune(svm, train.x = x_can_noncan, train.y = y_can_noncan, kernel = "radial", gamma = c(0.1,.5,1,2), ranges = list(cost=10^(-1:2)))
      #svm_can_noncan <- svm(x_train, y_train, kernel = "radial", gamma = svm_tune$best.parameters$gamma, cost = svm_tune$best.parameters$cost)
      
      svm_tune <- tune(svm, train.x = x_exp_tot, train.y = y_exp_tot, kernel = "linear", ranges = list(cost=10^(-1:2)))
      svm_exp_tot <- svm(x_train, y_train, kernel = "linear", cost = svm_tune$best.parameters$cost)
      
      svm_pred <- predict(svm_exp_tot, x_test)
      conf_exp_tot <- table(Predicted = svm_pred, Reference = y_test)
      
      # Accuracy total and per condition
      acc_tot <- (conf_exp_tot[1,1] + conf_exp_tot[2,2])/(sum(conf_exp_tot))
      acc_global <- conf_exp_tot[1,1]/(sum(conf_exp_tot[,1]))
      acc_scrambeled <- conf_exp_tot[2,2]/(sum(conf_exp_tot[,2]))
      
      # Save data
      accuracy_c <- c(sub, roi, hem, acc_tot, acc_global, acc_scrambeled)
      parameters_c <- c(sub, roi, hem, svm_tune$best.parameters$gamma, svm_tune$best.parameters$cost)
      
      accuracy <- rbind(accuracy, accuracy_c)
    }
    
  }
  
}


## Prepare results
col_names <- c("sub", "roi", "hem", "acc_tot", "acc_global", "acc_scrambeled")
colnames(accuracy) <- col_names

write.table(accuracy, file = "SVM_Ventral_Global_acc.txt", sep = ",", row.names = FALSE, col.names = TRUE)
