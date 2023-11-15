# -*- coding: utf-8 -*-
"""
Creating probabilitic activation maps from Neurosynth data

Created on Wed Aug  2 13:41:55 2023

@author: JoshB

- Downloads neurosynth data (version 7 - can be modified) and converts into NiMARE dataset using individual terms (can be modified to topics)
- Conducts coordinate-based meta-analyses using Multiple Kernel density analysis for all Neurosynth terms with Montecarlo multiple comparisons test
- Extracts time-series data from probabilistic maps using voltron_400 parcellation (can be modified to the preferred parcellation)
- Can z-score time-series however data already has zeroes and z-scoring has no effect

"""

# To install NiMARE
# Run "pip install nimare" in terminal
# Refer to "https://nimare.readthedocs.io/en/stable/installation.html" for details on installation

# biopython is unnecessary here, but is required by download_abstracts.
# We import it here only to document the dependency and cause an early failure if it's missing.
import Bio  # pip install biopython

from nimare.extract import download_abstracts, fetch_neurosynth
from nimare.io import convert_neurosynth_to_dataset
from nimare.correct import FWECorrector
from nimare.dataset import Dataset
from nimare.meta.cbma.mkda import MKDADensity

import nibabel as nib
from nilearn.maskers import NiftiLabelsMasker
import os
from pprint import pprint
import pandas as pd
from scipy.io import savemat


#%% Using NiMARE to setup neurosynth

# Download Neurosynth
out_dir = os.path.abspath("C:/PythonPrograms/NiMARE/example_data/") # Where to download data to
os.makedirs(out_dir, exist_ok=True)

# Takes Version 7 of Neurosynth, can be changed to preferred version
files = fetch_neurosynth(
    data_dir=out_dir,
    version="7",
    overwrite=False,
    source="abstract",
    vocab="terms",
)
# Note that the files are saved to a new folder within "out_dir" named "neurosynth".
pprint(files)
neurosynth_db = files[0]

# Convert Neurosynth database to NiMARE dataset file
# In neurosynth_db, you can modify the files (features + vocabulary) used e.g. instead of using terms you can use LDA topics for annotations (features)
# data-neurosynth_version-7_vocab-LDA400_source-abstract_type-weight_features.npz for LDA400 Topics
# data-neurosynth_version-7_vocab-terms_source-abstract_type-tfidf_features.npz for individual terms
ns_dset = convert_neurosynth_to_dataset(
    coordinates_file=neurosynth_db["coordinates"],
    metadata_file=neurosynth_db["metadata"],
    annotations_files=neurosynth_db["features"],
)
# Save dataset without abstracts
ns_dset.save(os.path.join(out_dir, "neurosynth_dataset.pkl.gz"))
print(ns_dset)

# Add article abstracts to dataset
ns_dset = download_abstracts(ns_dset, "jtan5724@uni.sydney.edu.au") # Change to your uni email
# Save dataset with abstracts
ns_dset.save(os.path.join(out_dir, "neurosynth_dataset_with_abstracts.pkl.gz"))

# How many studies in dataset
print(f"There are {len(ns_dset.ids)} studies in the Neurosynth database.")


#%% Neurosynth meta-analysis and parcellation

# Load neurosynth dataset
ns_dset = Dataset.load("C:/PythonPrograms/NiMARE/example_data/neurosynth_dataset_with_abstracts.pkl.gz")

# Load in custom atlas (voltron)
atlas = nib.load('C:/Users/JoshB/OneDrive/Documents/MATLAB_Analysis/MATLAB/Functions/schaefer_parcellation/voltron_1000.nii.gz') # path to custom atlas .nii.gz

# Create masker object - using voltron labels, z-score data
masker = NiftiLabelsMasker(labels_img=atlas)

# Extract list of neurosynth features/terms
annotations_list = list(ns_dset.annotations)
start = 3 # First neurosynth term
start = 46
ns_terms = annotations_list[start:]

# Save neurosynth terms as csv
ns_terms2 = pd.DataFrame(ns_terms)
ns_terms2.to_csv('C:/PythonPrograms/NiMARE/ns_LDA50_terms.csv')

# Run coordinate-based meta-analysis for each neurosynth term
# Multilevel Kernel Density Analysis algorithm
meta_ts = []
i = 1
while i < len(ns_terms):
    # Subsample neurosynth dataset using terms
    print(f"Iteration {i} out of {len(ns_terms)-1}...")
    msg = f'Running meta-analysis for term: {ns_terms[i]}'
    print(msg, end='\r', flush=True)
    temp_ids = ns_dset.get_studies_by_label(ns_terms[i])
    temp_dset = ns_dset.slice(temp_ids)
    print(f"\nThere are {len(temp_ids)} studies labeled with '{ns_terms[i]}'.")
    
    # Run Multilevel Kernel Density Analysis
    meta = MKDADensity()
    results = meta.fit(temp_dset)
    
    # Montecarlo multiple comparisons test
    # Available methods are 'montecarlo', 'bonferroni'
    corr = FWECorrector(method="montecarlo", n_iters=10, n_cores=1)
    cres = corr.transform(results)
    
    meta_map = cres.get_map("z_level-voxel_corr-FWE_method-montecarlo")
    temp_ts = masker.fit_transform(meta_map)
    meta_ts.append(temp_ts)
       
    i = i + 1
    
# Save time-series
savemat('C:/PythonPrograms/NiMARE/ns_LDA50_activity_44-49.mat', {'ts': meta_ts})