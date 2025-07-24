# The engagement of the cerebellum and basal ganglia enhances expertise in a sensorimotor adaptation task
This repository contains the code used in my project *"The engagement of the cerebellum and basal ganglia enhances expertise in a sensorimotor adaptation task"* https://doi.org/10.1162/imag_a_00271. 

**Code** was written using a combination of MATLAB and Python scripts.

**Data** was downloaded from *OpenNeuro* [ds004021](https://openneuro.org/datasets/ds004021/versions/1.0.0) and originally collected from this [paper](https://academic.oup.com/cercor/article-abstract/33/8/4761/6761518?login=false)
## code
The code folder contains the code required to run the analyses and produce the figures. The following scripts will be described in an order that fits with the manuscript and analyses. Note that these scripts are not functions and should be viewed similar to a notebook.
- [dataprep.m](https://github.com/ShineLabUSYD/VM_Adaptation/blob/main/Code/dataprep.m) includes code for processing the data after pre-processing. Reads in the behavioural data (.tsv) and timeseries data (.mat), and removes outliers, missing values, as well as normalises number of trials per subject (refer to manuscript for more details). **The outputs from this script are used for all analyses described below.** This script contains code to generate Figures 1 and 2A.
- [behavioural_analyses.m](https://github.com/ShineLabUSYD/VM_Adaptation/blob/main/Code/behavioural_analysis.m) compares behavioural performance both within and across task conditions. This script also has code to produce Figures 2B/C, and Figure 4A.

- [glm.m](https://github.com/ShineLabUSYD/VM_Adaptation/blob/main/Code/glm.m) creates all the *general linear models (GLM)* described in the manuscript. The script also contains code for comparing brain maps against [Neurosynth data](https://neurosynth.org/), and produces Figure 3.
- [clustering_energy.m](https://github.com/ShineLabUSYD/VM_Adaptation/blob/main/Code/clustering_energy.m) runs *K-means Clustering*, *inter-trial correlations* and the *energy landscape analysis* ([Munn et al., 2021](https://www-nature-com.ezproxy.library.sydney.edu.au/articles/s41467-021-26268-x)). Original energy landscape code can be found [here](https://github.com/ShineLabUSYD/Brainstem_DTI_Attractor_Paper). This script also produces Figures 4B/C, and Figure 5
- [supp_material.m](https://github.com/ShineLabUSYD/VM_Adaptation/blob/main/Code/supp_material.m) runs the main analyses from [glm.m](https://github.com/ShineLabUSYD/VM_Adaptation/blob/main/Code/glm.m) but uses 1000 Schaefer cortical nodes.

We also provide code used to download and extract meta-analyses from the *Neurosynth database* ([neurosynth_analyses.py](https://github.com/ShineLabUSYD/VM_Adaptation/blob/main/Code/neurosynth_analysis.py)).

## Visualisation
The visualisation folder contains the code required to create the brain visualisations. Details regarding dependencies are described in the README file.
