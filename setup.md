---
title: Setup
---

This workshop is designed to be conducted on the [NeSI](https://www.nesi.org.nz) compute infrastructure. All software and data is already set up for you to use during the workshop.

Running the workshop locally on your own computer will involve installing the required programs and downloading the example data. Some software packages may not work on all operating systems. Below are the requirements for running locally.


## Software

The software needed to run this lesson are listed below. Versions reported are the ones we use on NeSI but the lesson is likely to work with a range of versions of these softwares.

A Unix shell is required. See [https://carpentries.github.io/workshop-template/#setup](https://carpentries.github.io/workshop-template/#setup) for instructions on installing.

| Software      | Version | Manual      | Description 	|
| ----------- | ----------- | ----------- | ----------- |
|R | 4.1.0 | [link](https://www.r-project.org/) | Statistical computing language.|
|cutadaptQC|||
|QIIME2/2021.4|
|vsearch||https://github.com/torognes/vsearch||


R packages

```r
install.packages(c('phyloseq', 'tibble', 'ggplot2', 'dplyr',
                    'tidyr', 'ape', 'vegan','stringr'))
```


## Data

The data is already set up and available for you on NeSI.




{% include links.md %}
