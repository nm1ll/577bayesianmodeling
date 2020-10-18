# Statistics 577 - Introduction to Bayesian Modeling
## Final Project

The primary aim of the analysis is to quantify the effect of beta-carotene treatment on serum beta-carotene levels over time.

In addition, there are two secondary aims for the analysis. (1) Quantify whether the effect of treatment on serum beta-carotene differs by age, gender, BMI, or cholesterol. (2) Quantify the effect of treatment on serum vitamin E levels over time, and determine if serum vitamin E levels are correlated with serum beta-carotene levels over time.

It was a group project, and my part of the project was to work on model building and model checking. The task: "Discuss any statistical work you did to find and check
your final model. If you considered a range of possible models in your analysis, you should explain how you picked the major models you reported on in the previous section. You should also include discussion of Bayesian convergence metrics and a sensitivity analysis examining the importance of your choice of prior for your final models."

This repository contains the code I've used to complete the task. Languages used are R and OpenBUGS. The latter is primarily accessed through the R package "R2OpenBUGS" to facilitate for easier scripting and model building in R. As such, R scripts access OpenBUGS by running Bayesian models recorded in txt files.

Files starting with 1 answer the first aim of the study, to determine how the beta-carotene levels change over time when given treatment. Files starting with 2 address the task to determine how the vitamin E levels change over time when given treatment for beta-carotene.
