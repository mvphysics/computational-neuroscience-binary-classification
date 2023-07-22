# computational-neuroscience-binary-classification
* The project needs the MATLAB EEG Lab add-on in order to work.
* The hjorth.m function was found on github https://github.com/donnchadh/biosig/blob/master/biosig/t300_FeatureExtraction/hjorth.m
* fitAR.m, GCI.m, CGCI.m, hjorth.m are the white noise process, Granger causality, conditional Granger causality and Hjorth parameters implemetations respectively used in the project.m.
* project.m is the main project and project_experiments.m is an expansion to the main project and can run only after the project.m has run as it uses variables from the project file. Please do not run other projects in between so that the expansion obtains the correct values for the variables.
* The feature_selection.pdf contains all the features that create each dataset. In all the datasets, except dataset and dataset3, the features are ranked according to their significance. The classification_results.pdf contains the accuracy values of all the models that have been created.
