# imc_tools

Some useful tools to process IMC data from Hyperion.

The pipeline now is composed of:
- Segmentation of images from MCD files (with Steinbock, https://bodenmillergroup.github.io/steinbock, and custom hot pixel filtering)
- Sample merge (if necessary) to process combined analyses
- IMC Data Analysis with https://bodenmillergroup.github.io/IMCDataAnalysis/
