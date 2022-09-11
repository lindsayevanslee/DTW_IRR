# DTW_IRR

This is the code that was used to produce results for my dissertation for the MSc in Applied Statistics program at the University of Oxford. My dissertation mentor was Dr. Órlaith Burke (Nuffield Department of Population Health, University of Oxford). 

The paper is available upon request. This repository is now archived.

## Abstract

Wearable cameras are increasingly being used in research as a means to gather more robust information about people’s daily activity levels. However, the quantity and complexity of the data produced by these cameras requires the development of an image annotation scheme that is sufficiently descriptive yet able to be implemented consistently by multiple researchers. In order to analyze the ability of the image annotation protocol to produce consistent annotations, inter-rater reliability (IRR) needs to be assessed. Traditional methods of assessing IRR like Cohen’s kappa do not utilize the unique features of the image annotation data, so a new method is needed.

In this report we implement a dynamic time warping (DTW) algorithm with the longest common prefix (LCP) distance metric to assess IRR for annotations of images produced by wearable cameras. Two raters each produced image annotations for twelve study participants. We use DTW with four different step patterns to align the two time series and calculate the normalized distance for each warping path. This normalized distance is used as the new metric of inter-rater reliability. We also implement DTW for randomly simulated annotations in order to provide a baseline normalized distance to compare against the normalized distances for our observed data. We conclude that DTW with the LCP distance metric and Rabiner-Juang type 5 step pattern is the most appropriate method to analyze IRR for this type of data. It represents an improvement upon traditional measures of IRR because it utilizes the time series characteristics of the data and the hierarchical nature of the categories in the image annotation protocol. 
