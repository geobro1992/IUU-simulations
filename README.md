## Integrating information on fisher behavior into ecological fisheries models.
Fisheries are coupled social and ecological systems that exemplify the challenges arising from the intricate dynamics and feedback loops between human behavior and the natural environment. Assessments of the status of a fishery rely on data on fish population variability over time and estimates of harvest rate and reported catch. These models are used to understand the regulatory mechanisms underlying the relationship between fishing and population productivity. However, when fishing activity is illegal, unreported, and/or unregulated (IUU), such models can underestimate the vulnerability of fish stocks to collapse. There is a need, therefore, to develop general methods to incorporate information about human decisions regarding when and how much to fish, and compliance with reporting and fishing regulations into population assessments. Here we propose and assess a novel approach to model population demography of a fished population, incorporating ecological estimates of natural mortality, population censuses, and catch data. We describe the ways that knowledge of fisher motivations and behavior in different types of fisheries can inform estimates of reported and unreported fishing in our modeling framework. We consider four alternative modelling parameterizations that reflect real-world scenarios for which the degree of unreported fishing and the relationship between reporting rate and fish abundance varies. We show that in some cases ignoring IUU fishing can severely bias estimates of vital rates and population dynamics. Using social knowledge of fisher behavioral responses to changes in fish abundance (stock status) can inform estimates of IUU in model formulations and improve their predictive accuracy.

### Contents
#### Four scripts to create the model formulations in JAGS format for each type of fishery:
- subsistence_fishing_jags.R
- commercial_fishing_jags.R
- recreational_fishing_jags.R
- exotic_fishing_jags.R

#### Four scripts to simulate data from each fishery and run the JAGS models:
- subsistence_fishing_run.R
- commercial_fishing_run.R
- recreational_fishing_run.R
- exotic_fishing_run.R
