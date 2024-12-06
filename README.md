# Xenium Analysis
This repository is used to analyze the Xenium data in a basic way. 

## Basic Analsysis of Xenium Transcripts
Version 1.1 (Dec 2024) - adapt to the lastest version of xeniumranger v3.0

`basic_analysis_xenium_transcripts.ipynb` can be directly run for analysis. 2 Scenarios are included:
- whole Xenium sections
- selected area(s) on the Xenium sections

## Xenium Re-segmentation Using Cellpose3
#### Version log
- Version 0.0 - use Cellpose3 to re-seg: DAPI and interior RNA
- Version 1.0 - use original DAPI signals, binary (invalid); add functions to convert mask .tiff to 32-bit
- Version 1.1 - use original DAPI seg from xeniumranger output
- Version 2.0 - wrap
- Version 2.1 - add functions to filter weak RNA signals before re-seg using Cellpose3
#### Required software
- Cellpose3
- Xeniumranger v3.0 (latest version)
#### Usage
```
$ bash reseg_2step_pipe.sh -h

Usage: reseg_2step_pipe.sh -o outdir -x xenium_bundle -s sample_name [-d diameter] [-f filter_value]
  -o outdir         Specify the output directory
  -x xenium_bundle  Specify the PATH of the xenium bundle
  -s sampleid       Specify the name of the sample/slide
  -d diameter       (Optional) Specify the diameter for Cellpose3 cyto3 model (default: 70)
  -f filter_value   (Optional) Run channel intensity filter step with the given thresthold
  -h                Usage
```
#### Example
```
$ bash reseg_2step_pipe.sh \
  -o /dfs3b/ruic20_lab/tingty7/projects/DRG_spatial/DRG_Xenium/DRG_Xenium_241203/segmentation \
  -x /dfs3b/ruic20_lab/tingty7/data/DRG_Xenium/20241126__215026__drg/output-XETG00320__0055259__Region_2__20241126__215119 \
  -s s259_R2 \
  -d 70 \
  -f 2000
```

## Xenium Data Combination and Co-embedding with scRNA Reference
Please see `xenium_coembed.R` as an example. The process includes:
1) extract and integrate xenium data from different slides/regions
2) handle the scRNA reference (integration; if needed)
3) co-embedding xenium data with the scRNA reference using Seurat RCPA

