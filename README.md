# Analysis of DSBs in each sample 

> Please refer to html report for the data analysis code and nextflow pipeline for the processing of the sequence data

Main observations:

* Sample 1, 2, 4, 5, 6, 7, 8 have 0 breaks identified from the intersection with breaks in the predicted breaks bed file provided
* Sample 3 has 1 break which intersects with the predicted breaks
* Sample 9 has 5, Sample 10 has 10, Sample 11 also has 10, Sample 12 has 5, Sample 13 has 2, Sample 14 has 7, Sample 15 has 15, Sample 16 has 18 breaks which intersect with the predicted sites


Control versus treated samples:

According to the data analysis and the README description, control samples would be the samples where there are no intersections with AsiSI enzyme breaks defined in the bed file and treated samples would be the ones where the AsiSI restrictions sites are induced with DSBs. 

* From the BED file analysis, we can conclude that Sample 1, 2, 4, 5, 6, 7 and 8 are definitely control samples as they do not have any AsiSI induced DSBs
* Sample 3 and Sample 13 have 1 and 2 DSBs which can be due to sequencing artefacts or random sampling bias from the reads which are extracted from the actual FASTQ sequences, so we are not entirely sure if they are control samples or treated samples
* Sample 9, 10, 11, 12, 14, 15, 16 are treated samples as they have more number of DSBs overlapping with AsiSI restriction sites in the bed file provided.
* The maximum percentage of possible AsiSI DSBs in a single sample observed in the data is 6.7844 %.
