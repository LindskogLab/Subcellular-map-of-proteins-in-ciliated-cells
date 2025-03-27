# A high-resolution subcellular map of proteins in ciliated cells

Here, we provide the scripts used for the paper: A high-resolution subcellular map of proteins in ciliated cells
In this project, we investigated protein expression in motile ciliated cells in five different tissue types, using a multiplex immunohistochemistry (mIHC) panel. 
The idea with the panel was to target the motile ciliated cell using five markers targeting the cilium, transition zone, rootlet, cytoplasm, and nucleus, with a sixth position in the panel open for any candidate protein we wanted to study. 
This allowed us to see with which panel marker the candidate protein overlaped with, indicating its subcellular location. 

Scripts for this project:
- Cilia_image_analysis : This is an ImageJ Macros script for the automated image analysis
- Data_analysis_and_Heatmap : The script to process the output thata from the image analysis, and to generate a heatmap with clusters based on those results
- GO_enrichment : Gene ontology enrichment for the proteins in the study using the clusters identified from the heatmap
- CILIA_script_sup_fig_2-3 : Script for supplementary figure 2 and 3 (Fig. 2 shows scRNA levels in ciliated cells compared to different cell tupes, Fig.3 show Human Proteome Project (HPP) data for the clusters, showing a distribution of how well known the proteins are).




