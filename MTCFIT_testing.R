source(paste0(here::here(), "~/path/MTCFIT.R"))


#---------------- DATA WITH "OR" NOT SINIFICANT
simulated_data_IZ <- readRDS( "~/path/simulated_data_IZ.Rds")

MTCFIT(mydata = simulated_data_IZ, 
                    id = "id", 
                    trt = "treatment", 
                    Y = "outcome", 
                    typeY = "binomial", 
                    vars_keep = c('c1', 'c5', 'c9','c12'),
                    verbose = TRUE, # default is TRUE
                    alpha = 0.05, # default
                    myseed = 200, # default
                    ratio_k = c(1, 2, 3), # default
                    caliper = c(0.1, 0.2, 0.3), # default
                    cut_length = 30 # default, equal to the number of characters 
                    # allowed on the first line. the rest of the string appears 
                    # indented on the next line.
)






#--------------------------------------------------------------------------------

#---------------- DATA WITH "OR" SINIFICANT
more_sample <- readRDS( "~/path/simulated_data_IZ.Rds")
more_sample$treat <- ifelse(more_sample$c10 < 23.1,1,0)
drop <- c("c10")
more_sample =  more_sample[,!(names(more_sample) %in% drop)]


MTCFIT(mydata = more_sample,
                    id = "id",
                    trt = "treat",
                    Y = "outcome",
                    typeY = "binomial",
                    vars_keep = c('c1', 'c6', 'c8','c9'),
                    verbose = TRUE, # default is TRUE
                    alpha = 0.05, # default
                    myseed = 200, # default
                    ratio_k = c(1, 2, 3), # default
                    caliper = c(0.1, 0.2, 0.3), # default
                    cut_length = 30 # default, equal to the number of characters
                    # allowed on the first line. the rest of the string appears
                    # indented on the next line.
)

