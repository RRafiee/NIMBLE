## Shiny server for MB MassArray Classification
## Machine learning model and classifier code by Dr Reza Rafiee 2014-2017
## Adaptation and Shiny code: Dr Matthew Bashton


## load libraries
library(shiny)
library(e1071) #for SVM classifier
library(parallel) # For mclapply speeds up probability estimation
library(gtools) # Needed for numerically rather than lexicographically sorted strings
library(NanoStringNorm)

##### Threshold setting ####
# NanoString
threshold <- 0.6847321

## MB ##
# Need to load in samples here
# Will also need to call Ed's clean up funtion, get UI working first, before testing lastest NewGene data.

### Get input file name from UI
shinyServer(function(input, output) {
  
  
  #############################################################################
  ######################## Reactive classifier function #######################
  #############################################################################
  
  classifier <- reactive({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    cat("input file is here:\n")
    cat(inFile$datapath, "\n")
    
    withProgress(message = 'Processing data', value = 0, {
      
      # Start the clock
      cat("Timing start\n")
      ptm <- proc.time()
      
      # For raw data use cleanSeq4
      #incProgress(0.10, detail = paste("Processing and cleaning raw data"))
      # File system less data passing
      
      # ## MB changes here for legacy code compatibility and new BS conversion effiency data.
      # #Sample.test <- cleanSeq4(filename=inFile$datapath)
      # returned_data <- cleanSeq4(filename=inFile$datapath)
      # # Re-create old Sample.tests
      # Sample.test <- returned_data[[1]]
      # # New BS_Eff
      # BS_Eff <- returned_data[[2]]
      
      
      #raw_betas <- Sample.test
      #browser()
      
      # For cleaned data
      Sample.test <- as.data.frame(read.csv(inFile$datapath, header=TRUE))#,row.names=0))
      
      # Reza
      ######################################################################################################################
      # Reading an input csv file including raw gene expression data of 19 genes as well as conrtol genes
      
      NanostringSamplesGenes <- Sample.test
      # NanostringSamplesGenes <- as.data.frame(read.csv("~/NanoString/81RawNanoString08112016_TestSample.csv",header=T))
      head(NanostringSamplesGenes)
      # Subgroup medulloblastoma genes 
      # WNT_genes <- c("DKK2","EMX2","GAD1","TNC","WIF1")
      # SHH_genes <- c("ATOH1","EYA1","HHIP","SFRP1")
      # Grp3_genes <- c("GABRA5","IMPG2","MAB21L2","NPR3","NRL")
      # Grp4_genes <- c("EOMES","KCNA1","KHDRBS2","RBM24","UNC5D")
      
      All_genes_ordered <- c("DKK2","EMX2","GAD1","TNC","WIF1","ATOH1","EYA1","HHIP","SFRP1","GABRA5","IMPG2","MAB21L2","NPR3","NRL","EOMES","KCNA1","KHDRBS2","RBM24","UNC5D",
                             "NEG_A", "NEG_B", "NEG_C", "NEG_D", "NEG_E", "NEG_F", "NEG_G", "NEG_H", "POS_A", "POS_B", "POS_C", "POS_D", "POS_E", "POS_F", "GAPDH", "ACTB", "LDHA")
      
      # Order the subgroup-based genes of the input file
      rownames(NanostringSamplesGenes) <- NanostringSamplesGenes$Name
      NanostringSamplesGenes <- NanostringSamplesGenes[All_genes_ordered,]
      
      # Separating data from annotation
      NanostringSamples.anno <- NanostringSamplesGenes[,c(1:3)]
      NanostringSamples.data <- NanostringSamplesGenes[,-c(1:3)]
      
      # Normalising raw data
      NanoStriang_mRNA_norm_tech_bio <- NanoStringNorm(
        x = NanostringSamples.data, 
        anno = NanostringSamples.anno,
        CodeCount = 'geo.mean',
        Background = 'mean.2sd',
        SampleContent = 'housekeeping.geo.mean',
        round.values = TRUE,
        take.log = TRUE,  #log2 transformed
        return.matrix.of.endogenous.probes = TRUE
      );
      
      # Chnage return.matrix.of.endogenous.probes value into FALSE
      # Plot.NanoStringNorm(
      #   x = NanoStriang_mRNA_norm_tech_bio,
      #   label.best.guess = TRUE,
      #   #plot.type = c('volcano'),
      #   plot.type = c('norm.factors'),  #'cv', 'mean.sd', 'positive.controls', 'RNA.estimates'
      #   title = FALSE
      # )
      
      
      # NanoString.mRNA.norm.tech.bio # results
      # Normalising gene expression data into range of [0 1]
      Nano1_log2_Normalised <- (NanoStriang_mRNA_norm_tech_bio-min(NanoStriang_mRNA_norm_tech_bio))/(max(NanoStriang_mRNA_norm_tech_bio)-min(NanoStriang_mRNA_norm_tech_bio))
      Nano1_log2_Normalised <- as.data.frame(Nano1_log2_Normalised)
      
      Sample_test <- Nano1_log2_Normalised
      Sample.test <- Sample_test
      
      #############################################################
      # Quality Control1: checking Normalization performance
      # and rejecting those which have not meet the criterion
      #############################################################
      # Reza
      ######################################################################################################################
      # ## MB ##
      # # Not having 17 probes will cause issues later down stream in the analysis, these are the white listed ones
      # Probes_17 <- c("cg00583535", "cg18788664", "cg08123444", "cg17185060", "cg04541368", "cg25923609", "cg06795768", "cg19336198", "cg05851505", "cg20912770", "cg09190051", "cg01986767", "cg01561259", "cg12373208", "cg24280645", "cg00388871", "cg09923107")
      # # Replacing orginal object in case it's used below
      # Sample.test <- Sample.test[Probes_17,]
      seq.test.BEM.97 <- Sample.test
      Total.No.of.Samples <- ncol(seq.test.BEM.97)
      Original.No.of.Samples <-ncol(seq.test.BEM.97)
      
      #############################################################
      # Quality Control1: checking Normalisation results  
      # and rejecting those which have not meet the criterion
      #############################################################
      ## MB
      #incProgress(0.10, detail = paste("Checking for missing probes"))
      ############################################################
      # Should only run if we have more than 1 sample which has less than 7 missing probes 
      if (Total.No.of.Samples > 0 ) {
        
      ############################################################
      ## MB
        # incProgress(0.10, detail = paste("Loading training sets"))
        ## 13 October 2015, 220 Training set (225 - 5 Grp3 samples: NMB273, NMB376, NMB405, NMB666, NMB717)
         # Trainingset450k17WithSubgroup <- as.matrix(read.csv("220TrainingSet450KforSequenomClassifierwithSubgroupOriginal13Oct2015_ver4.csv",header=T,row.names=1))
         # Trainingset450k17 <- Trainingset450k17WithSubgroup[1:17,]
         # labels220 <- as.character(Trainingset450k17WithSubgroup[18,])
         # subgroup.labels <- factor(labels220)
         # y1 <- subgroup.labels
        
      #############################################################
      # Training set, n=101 (readcount log2 transformed), 
        incProgress(0.10, detail = paste("Loading training sets"))
        TrainingsetRNA_seq19 <- as.data.frame(as.matrix(read.delim("101RNA_seq19gene_log2_normalised01_5thApril2016.csv",sep=",",header=T)))
        sam_names <- as.character(TrainingsetRNA_seq19[,1])
        TrainingsetRNA_seq19 <- TrainingsetRNA_seq19[,-1]
        rownames(TrainingsetRNA_seq19) <- sam_names
        
        subgroup_labels <- as.integer(c(rep.int(1,11),rep.int(2,24),rep.int(3,33),rep.int(4,33)))# WNT=11, SHH=24, Grp3=33, Grp4=33
        y_training <- factor(subgroup_labels)
        y1 <- factor(subgroup_labels)
        
        ############################################################
        
        ## If no NAs in seq.test.BEM.97 then skip this bit too. ####
        # if (anyNA(seq.test.BEM.97 == TRUE)) {
        #   
        #   cat("\nWe have some missing NA probes\nImputing missing vlaues\n")
        #   
        #   ## MB
        #   incProgress(0.10, detail = paste("Imputation Modelling missing probes"))
        #   
        #   ## MB setting bounds for amelia
        #   amelia_bounds <- matrix(c(1:17, rep(0,17), rep(1,17)), nrow = 17, ncol = 3)
        #   
        #   # Multiple Imputation Modelling using Bootstrapping Expectation Maximization Algorithm
        #   # Handling missing probes
        #   # Amelia package
        #   # Using EM algorithm and Bootstrapping for handling missing probes
        #   # using passed samples for multiple imputation modelling :seq.test.BEM.97
        #   # Imputation cohort has now updated based on 101 NMBs instead of 103 NMBs
        #   Cohort101.test <- read.csv("101GoldCohortSeqDataAfterBEM17Probes08Oct2015.csv",header=T,row.names=1)
        #   Combined.cohort2 <- cbind(Cohort101.test,seq.test.BEM.97)
        #   set.seed(1234)
        #   Combined.datasets.models <- amelia(x = t(Combined.cohort2), m = 20, p2s = 1, frontend = FALSE, tolerance = 0.0001, bounds = amelia_bounds, max.resample = 20) # 18/12/2015
        #   summary(Combined.datasets.models)
        #   cl2_11 <- t(Combined.datasets.models$imputations[[1]])  # Saving the third imputed data
        #   cl2_11 <- cl2_11[,-c(1:ncol(Cohort101.test))]
        #   m <- 20 # number of imputation dataset
        #   # computing the average
        #   Bothcohort_combined <- cl2_11
        #   
        #   for (i in 1:nrow(Sample.test))  #should be 17 probes
        #   {
        #     for (j in 1:Total.No.of.Samples) #ncol(seq.test.BEM.97)) #the number of samples
        #     {
        #       if (is.na(seq.test.BEM.97[i,j]))
        #       {
        #         sum.imputed <- 0
        #         for (k in 1:m)
        #         {
        #           cl2_temp <- t(Combined.datasets.models$imputations[[k]])
        #           cl2_temp <- cl2_temp[,-c(1:ncol(Cohort101.test))]
        #           sum.imputed <- sum.imputed + cl2_temp[i,j]
        #         }
        #         avg.imputed <- sum.imputed/m
        #         Bothcohort_combined[i,j] <- avg.imputed
        #       }
        #     }
        #   }
        #   seq.test.BEM.97 <- Bothcohort_combined
        #   Total.No.of.Samples <- ncol(seq.test.BEM.97)
        # } else if (anyNA(seq.test.BEM.97 == FALSE)) {
        #   cat ("\nNo missing NA probes skiped Imputation Modelling\n")
        # }
        
        #############################################################
        ## Further analysis to assess confidence of calls ###########
        
        ## MB
        incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 1"))
        
        cost_1 <- 4.4 
        gamma_1 <- 0.01 
        kernel_1 <- "radial"
        
        x <- 1000 ## Number of iterations
        # Samples in rows 
        train.beta <- TrainingsetRNA_seq19 #mat_data_training  #Trainingset450k17
        amount <- round(0.9*nrow(train.beta))
        
        sel2 <- lapply(1:x, function(i) {
          set.seed(i)
          sample(1:nrow(train.beta), amount, replace=F)
        })
        
        ## MB this bit causes a delay
        incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 2"))
        linear.svms <- mclapply(1:x,
                                mc.cores=4,
                                function(i)  svm(x = TrainingsetRNA_seq19[sel2[[i]],],
                                                 y = y1[sel2[[i]]], scale = F,
                                                 tolerance = 0.00001, type = "C-classification",
                                                 kernel = kernel_1, cost = cost_1,
                                                 gamma=gamma_1, probability = T,
                                                 seed=i)
        )
        
        ## Test on Sequenom cases
        incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 3"))
        linear.tests <- mclapply(1:x,
                                 mc.cores=4,
                                 function(i) predict(linear.svms[[i]],
                                                     newdata=t(seq.test.BEM.97),
                                                     decision.values = T,
                                                     probability = T)
        )
        
        incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 4"))
        linear.calls <-lapply(1:x, function(i) linear.tests[[i]][1:19])
        
        incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 5"))
        prob.test <- (lapply(1:x,
                             function(i) attr(linear.tests[[i]], "probabilities"))
        )
        
        # old code for old conf interval 13/02/2015
        #probs2 <- predictConf(prob.test)
        
        ####################################### New Creating Pobes2 #################################
        # MB code from Reza to replace that of predictConf
        k <- FALSE
        for (j in 1:x) # the number of iterations
        {
          
          predProbTemp <-prob.test[[j]] # j iteration
          predProbTemp <- predProbTemp[,c("1", "2", "3","4")] # order the matrix based on the subgroup orders
          
          
          if (k == FALSE) # Making defult tables
          {
            
            predProbabilities <- matrix(ncol = 4, nrow =nrow(predProbTemp)*x, 0.0)
            predProbabilities <- predProbTemp
            colnames(predProbabilities) <- c("WNT","SHH", "Grp3", "Grp4")
            
            k <- TRUE
          }
          else
          {
            #Adding other iteration probabilities to the created table in the ordered columns
            predProbabilities <- rbind(predProbabilities,predProbTemp)
          }
        }
        
        
        probs2 <- matrix(ncol=nrow(predProbTemp),nrow=x,0.0)
        colnames(probs2) <- rownames(predProbTemp)
        
        for (ttt in 1:nrow(predProbTemp)) # number of samples
        {
          mmm <- matrix(ncol = 4, nrow =x, 0.0)
          colnames(mmm) <- c("WNT","SHH", "Grp3", "Grp4")
          gg <- 0
          for (fftt in 1:x)
          {
            gg <- gg + 1
            mmm[gg,] <- predProbabilities[ttt+nrow(predProbTemp)*(fftt-1),]
            #predProbabilities[1+3*0,1+3*1,1+3*2,1+nrow(predProbTemp)*(x-1)] # for the first sample, n=3
          }
          
          ProbSubgroup <- apply(mmm[,1:4],1,max)
          probs2[,ttt] <- ProbSubgroup
          
        }
        ####################################### End New Creating Pobes2 #############################
        
        #browser()
        presumed.class <- c(rep("Unknown",Total.No.of.Samples))
        
        i=1234
        incProgress(0.10, detail = paste("Assessing confidence of subgroup calls: stage 6"))
        model <- svm(TrainingsetRNA_seq19,y1,scale = F, tolerance = 0.00001, type = "C-classification", kernel = kernel_1, cost = cost_1, gamma=gamma_1, probability = T, seed=i)
        
        test.pred <- predict(object=model, newdata=t(seq.test.BEM.97), probability=TRUE)
        prob.test <- signif(attr(test.pred, "probabilities"), digits=2)
        maxProbs <- apply(prob.test,1,max)
        # New code here 13/02/2015
        # maxProbsWhich <- predict(model,newdata=t(seq.test.BEM.97))
        maxProbsWhich <- factor(test.pred[1:nrow(prob.test)],levels=c("1", "2", "3", "4"))
        maxProbsCol <- ifelse(maxProbsWhich==1,"blue",ifelse(maxProbsWhich==2,"red",
                                                             ifelse(maxProbsWhich==3,"yellow2","darkgreen")))
        
        maxProbsCol2 <- ifelse(maxProbsCol=="yellow2","#EEEE0066", ifelse(maxProbsCol=="blue","#0000FF66",
                                                                          ifelse(maxProbsCol=="darkgreen","#00640066","#FF000066")))
        
        # MB Output for classification table
        levels(maxProbsWhich) <- c("WNT", "SHH", "Grp3", "Grp4")
        #results.df <- data.frame(names(maxProbsWhich), maxProbsWhich, maxProbs, row.names = NULL)
        results.df <- data.frame(names(maxProbsWhich), as.character(maxProbsWhich), maxProbs, row.names = NULL, stringsAsFactors = FALSE)
        colnames(results.df) <- c("Sample", "Subgroup", "Confidence")
        
        # Stop the clock
        time<- (proc.time() - ptm)
        #browser()
        # Return all the things:
        classified_data <- list(results.df,
                                Sample.test,
                                #missing_summary,
                                Total.No.of.Samples,
                                Original.No.of.Samples,
                                #failed.samples,
                                #failed_sample_names,
                                #probe_threshold,
                                probs2,
                                maxProbsWhich,
                                maxProbs,
                                maxProbsCol,
                                maxProbsCol2,
                                time)#,
                                #BS_Eff)
        
        names(classified_data) <- c("results.df",
                                    "Sample.test",
                                    #"missing_summary",
                                    "Total.No.of.Samples",
                                    "Original.No.of.Samples",
                                    #"failed.samples",
                                    #"failed_sample_names",
                                    #"probe_threshold",
                                    "probs2",
                                    "maxProbsWhich",
                                    "maxProbs",
                                    "maxProbsCol",
                                    "maxProbsCol2",
                                    "time")#,
                                    #"BS_Eff")
        
        return(classified_data)
        
        # End if if (Total.No.of.Samples > 0 ) 
      } else if (Total.No.of.Samples == 0) {
        # At this point need to set up empty data if we've got nothing to classify
        
        # Stop the clock
        time<- (proc.time() - ptm)
        
        # Return all the things:
        classified_data <- list(Sample.test,
                                #missing_summary,
                                Total.No.of.Samples,
                                Original.No.of.Samples,
                                #failed.samples,
                                #failed_sample_names,
                                #probe_threshold,
                                time)#,
                                #BS_Eff)
        
        names(classified_data) <- c("Sample.test",
                                    #"missing_summary",
                                    "Total.No.of.Samples",
                                    "Original.No.of.Samples",
                                    #"failed.samples",
                                    #"failed_sample_names",
                                    #"probe_threshold",
                                    "time")#,
                                    #"BS_Eff")
        
        return(classified_data)
        
      }
      
      
      # Let user know we're done
      setProgress(value = 1, message = "Done!")
      cat("\nDone classification\n")
      
      
    }) # End with progress
    
    
    #########################################################################
    ################## End of classifier reactive function ##################
    #########################################################################
    
  }) # End reactive classifier function
  
  
  # Output classification_table #####
  
  output$classification_table <- renderDataTable({
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    
    # Now check for no of samples, if more than 0 run existing code
    Total.No.of.Samples <- classified_data$Total.No.of.Samples
    
    if( Total.No.of.Samples > 0 ) {
      # Now run existing code
      
      results.df <- classified_data$results.df
      # Change name to Subgroup Call of 2nd col
      colnames(results.df)[2] <- "Subgroup Call"
      # Add a Normalization QC column
      results.df[,4] <- "Pass"
      colnames(results.df)[4] <- "Normalization QC"
      failed.samples <- classified_data$failed.samples
      
      # Apply threshold and label samples as unclassifiable
      thresholded_results.df <- results.df
      i <- 1
      for (i in 1:nrow(results.df)) {
        if (!results.df[i,"Confidence"] > threshold) {
          thresholded_results.df[i,"Subgroup Call"] <- "Unclassifiable"
          thresholded_results.df[i,"Confidence"] <- NA
        }
      }
      
      # Inject failed samples that did not pass missing probes threshold test if we have them
      if (length(failed.samples > 0)) {
        failed_sample_names <- classified_data$failed_sample_names
        i <- 1
        new.results.df <- thresholded_results.df
        for (i in i:length(failed.samples)) {
          new.results.df <- rbind(new.results.df, c(failed_sample_names[i], "-", NA, "Fail"))
        }
        # Convert to percentage (because medics)
        new.results.df[,3] <- as.character(as.numeric(new.results.df[,3])*100)
        colnames(new.results.df)[3] <- "Probability %"
        # Now return the df but with NAs replaced by -
        new.results.df[is.na(new.results.df)] <- "-"
        # Sort via sample ID (correctly)
        new.results.df <- new.results.df[mixedorder(new.results.df[,1]),]
        return(new.results.df)
        # Where we don't have any failed samples simply return the df but with NAs
        # replaced by -
      } else {
        # Convert to percentage (because medics)
        thresholded_results.df[,3] <- as.character(as.numeric(thresholded_results.df[,3])*100)
        colnames(thresholded_results.df)[3] <- "Probability %"
        thresholded_results.df[is.na(thresholded_results.df)] <- "-"
        # Sort via sample ID (correctly)
        thresholded_results.df <- thresholded_results.df[mixedorder(thresholded_results.df[,1]),]
        return(thresholded_results.df)
      }
      # End failed sample injector
      
    } else if (Total.No.of.Samples == 0 ) {
      # Recreate results.df
      # Get failed samples and their names
      failed_sample_names <- classified_data$failed_sample_names
      i <- 1
      results.df <- data.frame(matrix(nrow = 0, ncol =4), stringsAsFactors = FALSE)
      for (i in 1:length(failed_sample_names)) {
        results.df <- rbind(results.df, c(failed_sample_names[i], "-", "-", "Fail"), stringsAsFactors = FALSE)
      }
      colnames(results.df) <- c("Sample", "Subgroup Call", "Probability %", "Normalization QC")
      new.results.df <- results.df[mixedorder(results.df[,1]),]
      return(new.results.df)
    }
    
    # End Total.No.of.Samples > 0
    
  })
  
  output$downloadClassification <- downloadHandler(
    filename = "MB_classification.csv",
    content =  function(file) {
      
      classified_data <- classifier()
      if (is.null(classified_data)) return(NULL)
      
      # Now check for no of samples, if more than 0 run existing code
      Total.No.of.Samples <- classified_data$Total.No.of.Samples
      
      # Also need to check Total.No.Of.Samples here
      if( Total.No.of.Samples > 0 ) {
        # Now run existing code
        
        results.df <- classified_data$results.df
        # Change name to Subgroup Call of 2nd col
        colnames(results.df)[2] <- "Subgroup Call"
        # Add a Normalization QC column
        results.df[,4] <- "Pass"
        colnames(results.df)[4] <- "Normalization QC"
        failed.samples <- classified_data$failed.samples
        
        # Apply threshold and label samples as unclassifiable
        thresholded_results.df <- results.df
        i <- 1
        for (i in 1:nrow(results.df)) {
          if (!results.df[i,"Confidence"] > threshold) {
            thresholded_results.df[i,"Subgroup Call"] <- "Unclassifiable"
            thresholded_results.df[i,"Confidence"] <- NA
          }
        }
        
        # Inject failed samples that did not pass missing probes threshold test if we have them
        if (length(failed.samples > 0)) {
          failed_sample_names <- classified_data$failed_sample_names
          i <- 1
          new.results.df <- thresholded_results.df
          for (i in i:length(failed.samples)) {
            new.results.df <- rbind(new.results.df, c(failed_sample_names[i], "-", NA, "Fail"))
          }
          # Convert to percentage (because medics)
          new.results.df[,3] <- as.character(as.numeric(new.results.df[,3])*100)
          colnames(new.results.df)[3] <- "Probability %"
          # Now return the df but with NAs replaced by -
          new.results.df[is.na(new.results.df)] <- "-"
          # Sort via sample ID (correctly)
          new.results.df <- new.results.df[mixedorder(new.results.df[,1]),]
          
          write.csv(new.results.df, file, row.names = FALSE)
          # Where we don't have any failed samples simply return the df but with NAs
          # replaced by -
        } else {
          # Convert to percentage (because medics)
          thresholded_results.df[,3] <- as.character(as.numeric(thresholded_results.df[,3])*100)
          thresholded_results.df[is.na(thresholded_results.df)] <- "-"
          # Sort via sample ID (correctly)
          thresholded_results.df <- thresholded_results.df[mixedorder(thresholded_results.df[,1]),]
          colnames(thresholded_results.df)[3] <- "Probability %"
          write.csv(thresholded_results.df, file, row.names = FALSE)
        }
        # End failed sample injector
        
      } else if (Total.No.of.Samples == 0) {
        # Recreate results.df
        # Get failed samples and their names
        failed_sample_names <- classified_data$failed_sample_names
        i <- 1
        results.df <- data.frame(matrix(nrow = 0, ncol =4), stringsAsFactors = FALSE)
        for (i in 1:length(failed_sample_names)) {
          results.df <- rbind(results.df, c(failed_sample_names[i], "-", "-", "Fail"), stringsAsFactors = FALSE)
        }
        colnames(results.df) <- c("Sample", "Subgroup Call", "Probability %", "Normalization QC")
        results.df <- results.df[mixedorder(results.df[,1]),]
        write.csv(results.df, file, row.names = FALSE)
      }
      
    }
  )
  
  ###################################
  
  
  # Output Missing probe summary ####
  # Note now changed to ouput "informative probes"
  
  output$mp <- renderDataTable({
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    missing_summary <- classified_data$missing_summary
    missing_table <- data.frame(rownames(missing_summary), 17-missing_summary[, 1], row.names=NULL, stringsAsFactors = FALSE)
    missing_table <- missing_table[mixedorder(missing_table[,1]),]
    colnames(missing_table) <- c("Sample", "Number of Informative Probes")
    missing_table
  })
  
  output$downloadMissing <- downloadHandler(
    filename = "MB_informative_probes.csv",
    content =  function(file) {
      classified_data <- classifier()
      if (is.null(classified_data)) return(NULL)
      missing_summary <- classified_data$missing_summary
      missing_table <- data.frame(rownames(missing_summary), 17-missing_summary[, 1], row.names=NULL, stringsAsFactors = FALSE)
      missing_table <- missing_table[mixedorder(missing_table[,1]),]
      colnames(missing_table) <- c("Sample", "Number of Informative Probes")
      write.csv(missing_table, file, row.names = FALSE)
    }
  )
  
  ###################################
  
  
  # Output number of failed samples #
  output$fs <- renderText({
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    
    Total.No.of.Samples <- classified_data$Total.No.of.Samples
    failed.samples <- classified_data$failed.samples
    failed_sample_names <- classified_data$failed_sample_names
    probe_threshold <- classified_data$probe_threshold
    
    if (length(failed.samples) > 0) {
      c(length(failed.samples), "sample(s) failed Normalization QC having", probe_threshold, "or more missing probes:", paste(failed_sample_names, collapse = ", "))
    } else if (length(failed.samples) == 0) {
      "All samples passed Normalization QC"
    }
    
  })
  
  ###################################
  
  
  # Output number of unclassifable samples
  
  # List samples (if any) that could not be classified above the threshold
  output$fc <- renderText({
    
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    
    # Now check for no of samples, if more than 0 run existing code
    Total.No.of.Samples <- classified_data$Total.No.of.Samples
    
    if(Total.No.of.Samples > 0) {
      # Now run existing code
      
      results.df <- classified_data$results.df
      unclassifiable <- results.df[results.df[3] < threshold, 1]
      
      if (length(unclassifiable) > 0) {
        c(length(unclassifiable), "samples(s) passing Normalization QC could not be confidently assigned a subgroup call:", paste(unclassifiable, collapse = ", "))
      } else if (length(unclassifiable) == 0) {
        "All samples passing Normalization QC were successfully assigned a subgroup"
      }
      
      # End Total.No.of.Samples > 0
    } else if (Total.No.of.Samples == 0) {
      "All samples failed Normalization QC, no samples can be classified"
    }
    
  })
  
  ###################################
  
  
  # Output BS con eff ###############
  
  output$BS_Eff <- renderDataTable({
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    
    BS_Eff_table <- classified_data$BS_Eff
    
    # Convert sample name factor to character vector
    BS_Eff_table[,1] <- as.character(BS_Eff_table[,1])
    # Change col name
    colnames(BS_Eff_table)[2] <- "Bisulphite conversion efficiency %"
    # New Order
    BS_Eff_table <- BS_Eff_table[mixedorder(BS_Eff_table[,1]),]
    # Round down percentages
    is.num <- sapply(BS_Eff_table, is.numeric)
    BS_Eff_table[is.num] <- lapply(BS_Eff_table[is.num], round, 1)
    BS_Eff_table
    
  })
  
  
  output$downloadBS_Eff <- downloadHandler(
    filename = "MB_BS_Eff.csv",
    content =  function(file) {
      classified_data <- classifier()
      if (is.null(classified_data)) return(NULL)
      
      BS_Eff_table <- classified_data$BS_Eff
      
      # Convert sample name factor to character vector
      BS_Eff_table[,1] <- as.character(BS_Eff_table[,1])
      # Change col name
      colnames(BS_Eff_table)[2] <- "Bisulphite conversion efficiency %"
      # New Order
      BS_Eff_table <- BS_Eff_table[mixedorder(BS_Eff_table[,1]),]
      # Round down percentages
      is.num <- sapply(BS_Eff_table, is.numeric)
      BS_Eff_table[is.num] <- lapply(BS_Eff_table[is.num], round, 1)
      write.csv(BS_Eff_table, file, row.names = FALSE)
    }
  )
  
  ###################################
  
  
  # Output graph ####################
  ## MB totally reworked to get sane graph of WNT, SHH, Grp3, Grp4
  output$classifierPlot <- renderPlot({
    
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    
    # Now check for no of samples, if more than 0 run existing code
    Total.No.of.Samples <- classified_data$Total.No.of.Samples
    
    if (Total.No.of.Samples > 0) {
      # Now run existing code
      
      probs2 <- classified_data$probs2
      maxProbsWhich <- classified_data$maxProbsWhich
      maxProbs <- classified_data$maxProbs
      maxProbsCol <- classified_data$maxProbsCol
      maxProbsCol2 <- classified_data$maxProbsCol2
      
      # Code to remove samples below threshold from plot
      cat(paste("Removing data points below threshold", threshold, "from graph:\n"))
      index <- maxProbs > threshold
      cat(names(maxProbs[!index]), "\n")
      new.probs2 <- probs2[,index]
      new.maxProbs <- maxProbs[index]
      new.maxProbsWhich <- maxProbsWhich[index]
      new.maxProbsCol <- maxProbsCol[index]
      new.maxProbsCol2 <- maxProbsCol2[index]
      new.Total.No.of.Samples <- length(maxProbs[index])
      
      par(mfrow=c(1,1))
      #par(mar=c(6,4,2,1) + 0.1)
      par(mar=c(6,4,4,1) + 0.1)
      par(cex=1.3)
      par(cex.axis=1)
      
      heading <- paste("Medulloblastoma subgroup call confidence intervals for", new.Total.No.of.Samples, "samples")
      
      boxplot(yaxt="n",xlab="",main=heading,ylab="Probability",new.probs2[,order(new.maxProbsWhich, new.maxProbs)],outpch=NA,ylim=c(0,1),las=2,
              col=new.maxProbsCol2[order(new.maxProbsWhich,new.maxProbs)] )
      
      abline(col="grey",lty = 1, h = threshold)
      
      # How many subgroups of each colour are we plotting
      tmp <- table(new.maxProbsCol)
      desired_col_order <-c("blue", "red", "yellow2", "darkgreen")
      to_sort <- names(tmp)
      # Re order by correct sub group col order using match on the desired_col_order vector
      tmp <- tmp[to_sort[order(match(to_sort,desired_col_order))]]
      # Index of where to draw the sub group deviders via cumsum
      grp.sum <- cumsum(tmp)
      # Add 0.5 to grp.sum for abline
      grp.sum <- grp.sum + 0.5
      # Index out final element of grp.sum to get rid of unwanted final abline
      grp.sum <- grp.sum[1:length(grp.sum)-1]
      # Check
      grp.sum
      abline(v=grp.sum)
      #lines(col="black",lwd=2,new.maxProbs[order(new.maxProbsWhich,new.maxProbs)])
      points(col=new.maxProbsCol[order(new.maxProbsWhich,new.maxProbs)],pch=19, new.maxProbs[order(new.maxProbsWhich,new.maxProbs)])
      legend("bottomleft", legend = c("WNT", "SHH", "Grp3", "Grp4"), col=c("blue", "red", "yellow2", "darkgreen"), pch=19)
      axis(2, las=2)
      
      # End if Total.No.of.Samples > 1
    } else if (Total.No.of.Samples == 0) {
      plot(0,0, xaxt = "n", yaxt = "n", xlab = '', ylab = '', frame.plot = FALSE, pch = 4, cex = 10, col = "red", main = "No samples to classify")
    }
    
  })
  # End output graph ################
  
  
  # Output graph download ####################
  ## MB totally reworked to get sane graph of WNT, SHH, Grp3, Grp4
  output$PlotDownload <- downloadHandler(
    
    filename = "MB_classification.png",
    content = function(file) {
      
      classified_data <- classifier()
      if (is.null(classified_data)) return(NULL)
      
      # Now check for no of samples, if more than 0 run existing code
      Total.No.of.Samples <- classified_data$Total.No.of.Samples
      
      if (Total.No.of.Samples > 0) {
        # Now run existing code
        
        probs2 <- classified_data$probs2
        maxProbsWhich <- classified_data$maxProbsWhich
        maxProbs <- classified_data$maxProbs
        maxProbsCol <- classified_data$maxProbsCol
        maxProbsCol2 <- classified_data$maxProbsCol2
        
        # Code to remove samples below threshold from plot
        cat(paste("Removing data points below threshold", threshold, "from graph:\n"))
        index <- maxProbs > threshold
        cat(names(maxProbs[!index]), "\n")
        new.probs2 <- probs2[,index]
        new.maxProbs <- maxProbs[index]
        new.maxProbsWhich <- maxProbsWhich[index]
        new.maxProbsCol <- maxProbsCol[index]
        new.maxProbsCol2 <- maxProbsCol2[index]
        new.Total.No.of.Samples <- length(maxProbs[index])
        
        heading <- paste("Medulloblastoma subgroup call confidence intervals for", new.Total.No.of.Samples, "samples")
        
        png(file, height = 1280, width = 1440)
        par(mfrow=c(1,1))
        par(mar=c(7,4,4,1) + 0.1)
        par(cex=2)
        par(cex.axis=1)
        boxplot(yaxt="n",xlab="",main=heading,ylab="Probability",new.probs2[,order(new.maxProbsWhich, new.maxProbs)],outpch=NA,ylim=c(0,1),las=2,
                col=new.maxProbsCol2[order(new.maxProbsWhich,new.maxProbs)] )
        
        abline(col="grey",lty = 1, h = threshold)
        
        # How many subgroups of each colour are we plotting
        tmp <- table(new.maxProbsCol)
        desired_col_order <-c("blue", "red", "yellow2", "darkgreen")
        to_sort <- names(tmp)
        # Re order by correct sub group col order using match on the desired_col_order vector
        tmp <- tmp[to_sort[order(match(to_sort,desired_col_order))]]
        # Index of where to draw the sub group deviders via cumsum
        grp.sum <- cumsum(tmp)
        # Add 0.5 to grp.sum for abline
        grp.sum <- grp.sum + 0.5
        # Index out final element of grp.sum to get rid of unwanted final abline
        grp.sum <- grp.sum[1:length(grp.sum)-1]
        abline(v=grp.sum)
        #lines(col="black",lwd=2,new.maxProbs[order(new.maxProbsWhich,new.maxProbs)])
        points(col=new.maxProbsCol[order(new.maxProbsWhich,new.maxProbs)],pch=19, new.maxProbs[order(new.maxProbsWhich,new.maxProbs)])
        legend("bottomleft", legend = c("WNT", "SHH", "Grp3", "Grp4"), col=c("blue", "red", "yellow2", "darkgreen"), pch=19)
        axis(2, las=2)
        dev.off()
        
        # End if Total.No.of.Samples > 1
      } else if (Total.No.of.Samples == 0) {
        png(file, height = 1280, width = 1440)
        plot(0,0, xaxt = "n", yaxt = "n", xlab = '', ylab = '', frame.plot = FALSE, pch = 4, cex = 10, col = "red", main = "No samples to classify")
        dev.off()
      }
      
    })
  # End output graph download ################
  
  # Output time taken ###############
  
  output$time <- renderText({
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    
    time <- classified_data$time
    
    c("Classification took", format(time[3]), "seconds")
  })
  
  ###################################
  
 # Output  gene expression values
  
 Gene_class <- c("DKK2"="WNT",
            "EMX2"="WNT",
            "GAD1"="WNT",
            "TNC"="WNT",
            "WIF1"="WNT",
            "ATOH1"="SHH",
            "EYA1"="SHH",
            "HHIP"="SHH",
            "SFRP1"="SHH",
            "GABRA5"="Grp3",
            "IMPG2"="Grp3",
            "MAB21L2"="Grp3",
            "NPR3"="Grp3",
            "NRL"="Grp3",
            "EOMES"="Grp4",
            "KCNA1"="Grp4",
            "KHDRBS2"="Grp4",
            "RBM24"="Grp4",
            "UNC5D"="Grp4")
            
  
  
    
  output$Beta <- renderDataTable(options = list(
    lengthMenu = list(c(10, 19, -1), c('10', '19', 'All')),
    pageLength = 19
  ), {
    classified_data <- classifier()
    if (is.null(classified_data)) return(NULL)
    betas <- round(classified_data$Sample.test, 2)
    betas[is.na(betas)] <- "-"
    Gene_class <- Gene_class[rownames(betas)]
    betas <- cbind(Gene_class, betas)
    betas <- cbind(rownames(classified_data$Sample.test), betas)
    colnames(betas) <- c("Gene ", colnames(betas)[-1])
    betas <- betas[order(betas[,"Gene_class"], betas[,"Gene "]),]
    # Blue Monday
    neworder <- mixedorder(colnames(betas[,c(-1,-2)]))+2
    betas <- betas[,c(1,2,neworder)]
    betas
  })
  
  output$downloadBeta <- downloadHandler(
    filename = "MB_beta_values.csv",
    content = function(file) {
      classified_data <- classifier()
      if (is.null(classified_data)) return(NULL)
      betas <- round(classified_data$Sample.test, 2)
      Gene_class <- Gene_class[rownames(betas)]
      betas <- cbind(Gene_class, betas)
      betas <- cbind(rownames(classified_data$Sample.test), betas)
      colnames(betas) <- c("Gene ", colnames(betas)[-1])
      betas <- betas[order(betas[,"Gene_class"], betas[,"Gene "]),]
      # Blue Monday
      neworder <- mixedorder(colnames(betas[,c(-1,-2)]))+2
      betas <- betas[,c(1,2,neworder)]
      write.csv(betas, file, row.names = FALSE)
    }
  )
  
}) # End shinyServer
