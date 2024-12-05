# Clear environment
rm(list=ls())

pwd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(pwd)

library(e1071)
library(caTools)
library(dplyr)
library(stringr)

## Prepare
conditions <- read.table("conditions_MVPA.csv", sep=',',header=T)

conditions$congruency <- as.numeric(str_replace_all(conditions$congruency, c("con" = "1", "in" = "2")))
conditions$form <- as.numeric(str_replace_all(conditions$form, c("c" = "1", "s" = "2")))
conditions <- conditions %>% select(-scrambling)

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
    
    r <- str_split(roi_sub, "_")[[1]][1]
    h <- str_split(roi_sub, "_")[[1]][2]
    h <- str_split(h, ".txt")[[1]][1]
    
    message(sprintf("%s %s", r, h))
    
    roi_c <- paste(path_c, roi_sub, sep='/')
    roi_data <- read.table(roi_c,sep=",",header=T)
    
    roi_data <- roi_data %>% 
      rename("run" = "betas1",
             "stim" = "betas2",
             "scrambling" = "betas3")
    roi_data <- roi_data[roi_data$run==1 | roi_data$run==2, ]
    
    roi_data <- merge(roi_data, conditions, by=c("stim", "run"), all=TRUE, sort=FALSE) %>% 
      drop_na() 
    
    if (ncol(roi_data) > 9){
      
      ## Data by experiment
      # Scrambling: 1 = global, 2 = scrambled; Congruency: 1 = congruent, 2 = incongruent; Form: 1 = Circle, 2 = Square
      y_exp_tot <- factor(roi_data$form, levels = c(1, 2))
      x_exp_tot <- select(roi_data, contains("betas"))
      
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
      acc_circle <- conf_exp_tot[1,1]/(sum(conf_exp_tot[,1]))
      acc_square <- conf_exp_tot[2,2]/(sum(conf_exp_tot[,2]))

      # Save data
      accuracy_c <- c(sub, r, h, acc_tot, acc_circle, acc_square)
      parameters_c <- c(sub, r, h, svm_tune$best.parameters$gamma, svm_tune$best.parameters$cost)
      
      accuracy <- rbind(accuracy, accuracy_c)
    }
    
  }
  
}


# Prepare results
col_names <- c("sub", "roi", "hem", "acc_tot", "acc_circle", "acc_square")
colnames(accuracy) <- col_names

write.table(accuracy, file = "SVM_Ventral_Global_Form_acc.txt", sep = ",", row.names = FALSE, col.names = TRUE)
