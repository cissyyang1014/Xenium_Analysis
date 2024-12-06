#conver 16-bit .tiff to 32-bit .tiff
#usage: python convert_tiff.py path_to_input_16bit.tif path_to_output_32bit.tif

import tifffile
import numpy as np
import argparse

# Function to convert 16-bit TIFF to 32-bit TIFF
def convert_to_32bit(input_path, output_path):
    # Load the 16-bit TIFF image
    img_16bit = tifffile.imread(input_path)
    
    # Convert the image to 32-bit
    img_32bit = img_16bit.astype(np.uint32)
    
    # Save the 32-bit TIFF image
    tifffile.imwrite(output_path, img_32bit)
    print(f"Converted {input_path} to 32-bit and saved as {output_path}.")

# Command line argument parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert 16-bit TIFF to 32-bit TIFF.")
    parser.add_argument('input_path', type=str, help="Path to the input 16-bit TIFF file.")
    parser.add_argument('output_path', type=str, help="Path to save the 32-bit TIFF file.")

    args = parser.parse_args()

    # Convert the provided TIFF file to 32-bit
    convert_to_32bit(args.input_path, args.output_path)


