# Instructions

This repo contains a jupyter notebook that produces tracking detector figures for the MAIA 10 TeV Muon Collider Detector paper. 

The notebook processes muon gun data with and without BIB, performs the analysis, and produces figures 

Noteboosk assume input files live in `json_data`, you can modify this as needed

For the latest samples please go to the OSG cluster into the directory `/scratch/lrozanov/mucolstudies` and download the following into your working directory:
 
```
# No BIB data
v2_noBIB_merged.json
v0_noBIB_all.json
v0_noBIB_all_5TeV.json

# BIB data
v0_BIB_all.json
v0_BIB_0_50.json
v0_BIB_50_250.json
v0_BIB_250_1000.json
v0_BIB_5TeV.json
```


