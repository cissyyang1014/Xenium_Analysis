#!/bin/bash
# 2-step segmentation, version 2.1
# add functions to filter weak RNA signals before re-seg using Cellpose3
# usage: $ bash reseg_2step_pipe.sh -o outdir -x xenium_bundle -s sample_name [-d diameter] [-f filter_value]

diameter=70
filter_value="" 

usage() {
  echo "Usage: $0 -o outdir -x xenium_bundle -s sample_name [-d diameter] [-f filter_value]"
  echo "  -o outdir         Specify the output directory"
  echo "  -x xenium_bundle  Specify the PATH of the xenium bundle"
  echo "  -s sampleid       Specify the name of the sample/slide"
  echo "  -d diameter       (Optional) Specify the diameter for Cellpose3 cyto3 model (default: 70)"
  echo "  -f filter_value   (Optional) Run channel intensity filter step with the given thresthold"
  echo "  -h                Usage"
  exit 0
}

while getopts "o:x:s:d:f:h" opt; do
  case $opt in
    o) outdir=$OPTARG ;;  
    x) original_bundle=$OPTARG ;;  
    s) sampleid=$OPTARG ;;    
    d) diameter=$OPTARG ;;   
    f) filter_value=$OPTARG ;; 
    h) usage ;;          
    *)
      usage             
      ;;
  esac
done

if [ -z "$outdir" ] || [ -z "$original_bundle" ] || [ -z "$sampleid" ]; then
  echo "Error: Missing required arguments"
  usage
fi

# main
echo "Xenium Bundle PATH: $original_bundle"
echo "SampleID: $sampleid"
output_dir=$outdir/$sampleid
mkdir -p $output_dir
echo "Output Directory: $output_dir"
mkdir -p $output_dir/morphology_focus

cd $output_dir
cp $original_bundle/cells.zarr.zip $output_dir/
cp $original_bundle/morphology_focus/morphology_focus_0002.ome.tif $output_dir/morphology_focus/
# mask DAPI signals on the RNA signals
python mask_orgDAPI_from_RNAs.py

cd $output_dir/morphology_focus/

# filter weak signals
if [ -n "$filter_value" ]; then
  python filter_weak_signals.py --file_path modified_morphology_focus_0002.ome.tif --threshold "$filter_value"
  image_PATH=$output_dir/morphology_focus/filtered_channel_3.tiff
  masks_PATH=filtered_channel_3_cp_masks.tif
else
  image_PATH=$output_dir/morphology_focus/modified_morphology_focus_0002.ome.tif
  masks_PATH=modified_morphology_focus_0002.ome_cp_masks.tif
fi


# segment using masked RNA signals using Cellpose3
python -m cellpose \
--image_path $image_PATH \
--pretrained_model cyto3 \
--chan 0 \
--chan2 0 \
--diameter $diameter \
--save_tif \
--verbose

# convert tif to 32-bit (optional, to get required input tiff for xeniumranger)
python convert_tiff.py $masks_PATH modified_morphology_focus_0002.ome_cp_masks_32bit.tif

cd $output_dir
# re-segment using the results of Cellpose3 and the original DAPI segmentation
xeniumranger import-segmentation --id=$sampleid \
--xenium-bundle=$original_bundle \
--nuclei=cells.zarr.zip \
--expansion-distance=0 \
--cells=$output_dir/morphology_focus/modified_morphology_focus_0002.ome_cp_masks_32bit.tif
