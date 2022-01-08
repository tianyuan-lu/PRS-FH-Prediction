# PRS-FH-Prediction

This is an R toolkit for combining parental history and existing polygenic risk scores for improved prediction of complex traits or disease risks.

The computational scripts in the [R](R/) folder can be directly downloaded. The only prerequisite is installation of the [pracma](https://cran.r-project.org/web/packages/pracma/index.html) R pacakge.

Refer to our [preprint](https://www.medrxiv.org/content/10.1101/2022.01.06.22268853v1) for detailed descriptions.

## Predicting continuous traits

Predicting a continuous trait for the children requires 

(1) an estimate of the proportion of variance explained by a polygenic risk score (*h2_PRS*);

(2) an estimate of the proportion of variance explained by a mid-parental predictor (*h2_midp*);

(3) data of the target population.

Both (1) and (2) can be obtained from large cohort studies that have similar genetic ancestries and demographic characteristics as the target population.

In command line:

```ruby
Rscript /path/to/predict_continuous.R h2_PRS h2_midp /path/to/data
```

The data file should be a tab-delimited file and should contain one column with the header **PRS** corresponding to the children's polygenic risk scores, and at least one column corresponding to the maternal or paternal trait measure, with the header **Mother** or **Father**. These columns should be standardized to have zero mean and unit variance.

```
Mother	Father	PRS
0.28236405580198  0.0478243107177947	-0.89184382599075
0.759665799839955	-0.00721719819976707	-0.237169912244541
1.01630718554635	0.785608541053388	0.863639182659233
-1.89783916306742	0.713878129174335	0.68661572124903
-1.33398782749324	-0.997270913147984	0.271900690073316
-1.9072137451484	-0.501673356594952	0.855450956488146
-0.199858668999688	0.926516478157096	1.23195377799939
0.571884596761782	-0.527252725302699	-1.60609961780322
0.543640566012509	0.711280253861894	0.00859903139352139
```

When both **Mother** and **Father** columns are supplied, prediction will be based on the children's polygenic risk scores and both parents' trait measures; If only one parental trait measure is supplied, prediction will be based on the children's polygenic risk scores and the provided parental trait measure. Additional columns can be included in the data file, and the order of columns can be flexible.

The joint predictor will be appended to the data file, which will then be stored in an RData with the affix **.joint.pred.continuous.RData**.

## Predicting binary disease outcomes

Predicting a binary trait (e.g. a disease outcome) for the children requires 

(1) an estimate of the baseline log-odds (*mu_0*): this should be estimated using a logistic regression model adjusted for covariate effects based on the reference population, but could be approximated by log-[prevalence / (1 - prevalence)] when a regression-based estimate is not available;

(2) an estimate of the log-odds ratio associated with one standard deviation increase in the polygenic risk score based on the reference population (*beta_PRS*);

(3) an estimate of the log-odds ratio associated with having a parental disease history based on the reference population (*beta_parent*): when the maternal history-based and the paternal history-based estimates are both available, this should be the average of the two estimates;

(4) an estimate of the log-men-to-women odds ratio based on the reference population (*beta_men*): this could be set to 0 when only one parent's disease history is provided or when sex has no influence on the disease risk;

(5) data of the target population.

In command line:

```ruby
Rscript /path/to/predict_binary.R mu_0 beta_PRS beta_parent beta_men /path/to/data
```

Again, the data file should a tab-delimited file and should contain one column with the header **PRS** corresponding to the children's polygenic risk scores, and at least one column corresponding to the maternal or paternal disease history (or binary trait measure), with the header **Mother_Disease** or **Father_Disease**. The **PRS** column should be standardized to have zero mean and unit variance.

```
Mother_Disease	Father_Disease	PRS
0	1	-0.89184382599075
0	0	-0.237169912244541
1	0	0.863639182659233
0	0	0.68661572124903
1	0	0.271900690073316
0	0	0.855450956488146
0	0	1.23195377799939
0	0	-1.60609961780322
0	1	0.00859903139352139
```

When both columns are supplied, prediction will be based on the children's polygenic risk scores and both parents' disease history; If only one parental disease history is supplied, prediction will be based on the children's polygenic risk scores and the provided parental disease history. Additional columns can be included in the data file, and the order of columns can be flexible.

The joint predictor will be appended to the data file, which will then be stored in an RData with the affix **.joint.pred.binary.RData**.

## Data example

A simulated toy dataset is provided [here](RData/simulate_phenotype.txt). Scripts used to generate these data are stored [here](R/simulate_phenotype.R). 

Briefly, for a continuous trait, a polygenic risk score was designed to explain 9% of the trait heritability while an orthogonal genetic component was designed to explain 64% of the trait heritability. Based on the data, we can find that the score explained 9.3% of the measured trait variance, and that a mid-parental predictor explained 27.7% of the measured trait variance.

We therefore generate a joint predictor with the following command using command line:

```ruby
Rscript predict_continuous.R 0.093 0.277 simulate_phenotype.txt
```

This should lead to an RData file similar to the [output](RData/simulate_phenotype.txt.joint.pred.continuous.RData).

And the joint predictor should explain 31.2% of the measured trait variance in the given dataset.

A binary outcome was simulated based on the above continuous trait with a baseline log-odds designed to be -1. Based on the data, we can find that per one standard deviation increase in the polygenic risk score is associated with a log-odds ratio of 0.247; having a parental disease history is associated with a log-odds ratio of 0.266; and the estimated baseline log-odds is -1.035. We do not consider sex effect in this example.

We therefore generate a joint predictor with the following command using command line:

```ruby
Rscript predict_binary.R -1.035 0.247 0.266 0 simulate_phenotype.txt
```

This should lead to an RData file similar to the [output](RData/simulate_phenotype.txt.joint.pred.binary.RData).

And the joint predictor should achieve an area under the receiver operating characteristic curve of 0.5744, higher than 0.5702 by the polygenic risk score alone.

## Estimating disease heritability using family disease history

Disease heritability can be estimated using the **estimate_h2** function in this [file](R/estimate_h2.R) with

(1) an estimate of the baseline log-odds (*mu_0*);

(2) an estimate of the log-odds ratio associated with having a parental disease history based on the reference population (*beta_parent*).

In an R console or RStudio:

```ruby
estimate_h2(mu_0, beta_parent, degree = 1)
```

Note that this function allows for *beta_parent* estimates based on 2nd-degree or more distant relatives, which can be specified by the option *degree*.




