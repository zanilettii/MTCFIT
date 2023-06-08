Scenario 1: 

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


Scenario 1 subset n=120

set.seed(100)
rand_df <- simulated_data_IZ[sample(nrow(simulated_data_IZ), size=130), ]

MTCFIT(mydata = rand_df, 
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

Scenario 2:

more_sample <- simulated_data_IZ
more_sample$treat <- ifelse(more_sample$c10 < 24.4,1,0)
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


Scenario 3:

data("lalonde")
lalonde2<- lalonde %>% mutate(id = 1:nrow(lalonde))
lalonde2$black <- ifelse(lalonde2$race == 'black', 1, 0) 
lalonde3<- lalonde2[-grep('race', colnames(lalonde2))] 

MTCFIT(mydata = lalonde3, 
       id = "id", 
       trt = "treat", 
       Y = "re78", 
       typeY = "gaussian",
       vars_keep = c("re74", "age", "educ", "married"),
       verbose = TRUE, # default is TRUE
       alpha = 0.05, # default
       myseed = 200, # default
       ratio_k = c(1, 2, 3), # default
       caliper = c(0.1, 0.2, 0.3), # default
       cut_length = 30 # default, equal to the number of characters 
       # allowed on the first line. the rest of the string appears 
       # indented on the next line.
)


