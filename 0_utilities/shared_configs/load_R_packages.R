

 

{
  
  
  ## Load in dependencies needed for BayesMVP (and other)
  require(Rcpp)
  require(RcppEigen)
  require(RcppParallel)
  require(RcppZiggurat)
  require(RcppClock)
  require(devtools)
  require(githubinstall)
  require(sitmo)  
  require(dqrng)  
  require(beepr)
  require(plyr)
  require(dplyr)
  require(ggplot2)
  require(patchwork)
  require(tictoc) 
  require(remotes)  
  require(jsonlite)
  require(expm)
  require(stringr)
  require(pracma)

  require(R6)
  require(roxygen2)

  require(Rfast)    
  require(LaplacesDemon)  
  require(TruncatedNormal) 
  require(rockchalk)    
  
  # 
  # remotes::install_github("https://github.com/roualdes/bridgestan",
  #                  #       subdir="~/Documents/Work/PhD_work/R_packages/BayesMVP",
  #                      force = TRUE)
  require(bridgestan)
  # # # # # # Install the cmdstanr package from CRAN
  # # we recommend running this in a fresh R session or restarting the current session
  # install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
  #  remotes::install_github("stan-dev/cmdstanr",force = TRUE)
  # # # Load cmdstanr package
  require(cmdstanr)
  # # 
  # # # # Install the latest version of CmdStan
  # install_cmdstan(cores = 16, overwrite = TRUE, cpp_options = list("STAN_MODEL_LDFLAGS" = "-shared",
  #                                              #  "LDFLAGS" = "-shared" ,
  #                                                 "CXXFLAGS" = "-fPIC"))
  #                                                # "CXX" = "/opt/AMD/aocc-compiler-4.2.0/bin/clang++") )
  # #               
  # # 

  # install.packages("StanHeaders")
  
  ## Mplus [Demo] must be installed seperately (i..e, outside of R) first!
  require(MplusAutomation)
  
  
  # ## Now load in BayesMVP 
  # require(BayesMVP)
  
  
  
  
}
 

 



















