# filter weak signal based on the channel intensity
# usage: python filter_weak_signals.py --file_path modified_morphology_focus_0002.ome.tif --threshold 2000

import tifffile as tiff
import numpy as np
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Process a 3rd channel of a TIFF file with a threshold.")
parser.add_argument("--file_path", required=True, help="Path to the input TIFF file")
parser.add_argument("--threshold", type=int, required=True, help="Intensity threshold for filtering")
args = parser.parse_args()

file_path = args.file_path
threshold = args.threshold

image = tiff.imread(file_path)
print(f"Image shape: {image.shape}")

channel_3 = image[2]

# Plot the histogram of the intensity 
hist_path = "channel_intensity_histogram.png"
plt.figure(figsize=(8, 6))
plt.hist(channel_3.flatten(), bins=256, range=(0, channel_3.max()), color='blue', alpha=0.7)
plt.xlim(0, 5000)
plt.ylim(0, 200000)
plt.title("Intensity Histogram for 3rd Channel")
plt.xlabel("Pixel Intensity")
plt.ylabel("Frequency")
plt.grid(True)
plt.savefig(hist_path)
print(f"Histogram saved to {hist_path}")
plt.close()

# filter
filtered_channel_3 = np.where(channel_3 > threshold, channel_3, 0)

# Plot the original image
original_channel_path = "original_channel_3.png"
plt.figure(figsize=(10, 10))
plt.imshow(channel_3, cmap='gray')
plt.title("Original 3rd Channel")
plt.axis('off')
plt.savefig(original_channel_path)
print(f"Original channel image saved to {original_channel_path}")
plt.close()

# Plot the filtered image
filtered_channel_path = "filtered_channel_3.png"
plt.figure(figsize=(10, 10))
plt.imshow(filtered_channel_3, cmap='gray')
plt.title(f"Filtered 3rd Channel (Intensity > {threshold})")
plt.axis('off')
plt.savefig(filtered_channel_path)
print(f"Filtered channel image saved to {filtered_channel_path}")
plt.close()


filtered_tiff_path = "filtered_channel_3.tiff"
tiff.imwrite(filtered_tiff_path, filtered_channel_3.astype(np.uint16))
print(f"Filtered image saved to {filtered_tiff_path}")
