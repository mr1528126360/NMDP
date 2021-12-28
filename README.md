# NMDP: A novel Multi-omics Drug response predic-tion model based on weighted edge sparse PCA and deep learning

Abstract
Motivation: The goal of precision oncology is to use the genomics information to customize the best treatment plan for the patient. Among them, predicting the drug response of patients based on multi-omics genetic data is the key issue. The latest research shows that building models based on deep learning can effectively predict drug response using multi-omics genetic data. However, the existing deep learning drug response prediction models have three main problems. First, ge-nomics data are characterized by high dimensionality and small sample size, the existing deep learning models have a high risk of overfitting. Second, the existing deep learning prediction mod-els lack biological interpretability. Third, there is still much room for improvement in the results of the existing models.  
Results: We propose NMDP, a novel drug response prediction model based on weighted edge sparse principal component analysis (SPCA), similarity network and deep learning. The NMDP model contains three modules including feature selection module, data fusion module and drug response prediction module. In the feature selection module, we propose a weighted sparse PCA method based on gene networks to extract key multi-group biometric features. In the data fusion module, we use the sample similarity network to solve the problem of sample dimensionality and data fusion. In the drug response prediction module, we use convolutional neural network (CNN) to build a classifier to predict drug response. We conduct experiments on 68 drugs including both targeted therapy drugs and non-specific drugs. Experimental results show that the average sensi-tivity and specificity of the NMDP model can reach 0.98 and 0.99 respectively, which is 11%-54% higher than the comparison models. Bio-enrichment analysis shows that the NMDP model can find target gene probes and pathway information that are highly related to the target drug. In conclusion, the NMDP model can achieve high-performance prediction of drug response while possessing biological interpretability.

The copy_number folder is a processing program for copy_number genomics

The methylation folder is a processing program for methylation genomics

The seq-file folder is a processing program for seq-file genomics

The total folder is a processing program for build network and classification 
