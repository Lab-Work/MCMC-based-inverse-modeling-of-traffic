MCMC-based-inverse-modeling-of-traffic
======================================

Supporting source code for the article, "Markov Chain Monte Carlo based inverse modeling of traffic flows using GPS data," by Tossavainen and Work 

This project contains the following:
1)All results presented in the article
2)M-files to reproduce all results presented in the article

How to use the code:
  1) Create the true state (v-field and trajectories) by running the script:
    simulateTrueStateUsingRiemannSolver.m 
    This will save the true velocity field and all vehicle trajectory information on the disk.
  2) Recover parameters:
    a)Using same 3-link configuration as in the true state generation by
      running recoverThreeDiagramsWithRatios.m
    b)Using the 4-link configuration by running 
      recoverFourDiagramsWithRatios.m
    c)Using 3-link configuration with poor mixing of the chain by
      running recoverThreeDiagramsNORATIOS.m
      
