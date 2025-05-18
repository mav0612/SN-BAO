Main files for parameter estimation

Setup_SN_BAO.m : Loads data, covariance matrices

ChiSq_SNBAO.m : computes Log likelihood (-1/2* chisquared) values (using both data sets)  V0V1 model

mcmc3.m : Markov Chain Monte Carlo code

MakeMCMC.m : Generates likelihood distribution. Use plotmatrix(m) to generate triangle plots

columnStats.m : displays means and standard deviations



D2combow0wa.m: : computes Log likelihood (-1/2* chisquared) values (using both data sets)  w0wa model

generateFakeData.m: makes fake dataset based on a V0V1 model using covariance from data

D2comboF.m: computes Log likelihood from this fake data

fakerun.m: generates parameter estimates with ensemble of fake data sets



sweepV1.m: Template for finding best fit parameter values over range of V1 values. Need to convert log likelihood calculator to chi squared calculator
