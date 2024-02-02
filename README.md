# imc_tools

Some useful tools to process IMC data from Hyperion.

The pipeline now is composed of:
- Image Extraction (with Steinbock) for filtering
- Segmentation of images from MCD files (with Steinbock, https://bodenmillergroup.github.io/steinbock, and custom hot pixel filtering)
- Sample merge (if necessary) to process combined analyses
- IMC Data Analysis with https://bodenmillergroup.github.io/IMCDataAnalysis/

```
bash imc_extraction.sh -r ROOTFOLDER
```
where ROOTFOLDER if the folder cointaing the dir samples

```
bash imc_segmentation.sh -r ROOTFOLDER -p panel.tsv
```
where ROOTFOLDER if the folder cointaing the dir samples and panels.tsv is the file containg markers

```
python3 merge_samples.py IN_PATH OUT_PATH --has_neighbors
```
where IN_PATH if the folder cointaing the dir samples together and OUT_PATH is the new folder containg the merged samples