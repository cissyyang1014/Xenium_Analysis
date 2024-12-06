import zarr
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tiff

# Load the Zarr file
z = zarr.open("cells.zarr.zip", mode='r')

# Load the masks
mask_0 = np.array(z['masks/0'])
mask_1 = np.array(z['masks/1'])

# Visualize both masks
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.imshow(mask_0, cmap='gray')
plt.title("Mask 0")

plt.subplot(1, 2, 2)
plt.imshow(mask_1, cmap='gray')
plt.title("Mask 1")

# Save the plot with a higher resolution (e.g., 300 DPI)
plt.savefig("nuclei_and_cell_masks_high_res.png", dpi=300)  # 300 DPI for high resolution

# Display the plot
plt.show()

 
# Load the nuclei segmentation (mask_0) from the Zarr file
nuclei_segmentation_mask = mask_0  # Assuming mask_0 is already loaded

# Load the RNA stain image (interior signal from OME-TIFF format)
ome_tiff_path = 'morphology_focus/morphology_focus_0002.ome.tif'  # RNA interior signal image
ome_image = tiff.imread(ome_tiff_path)

# Make a copy of the original OME-TIFF image
original_ome_image = ome_image.copy()

# Check shapes
print(f"Nuclei Segmentation Mask Shape: {nuclei_segmentation_mask.shape}")
print(f"OME Image Shape: {ome_image.shape}")

# Ensure the segmentation mask matches the OME image dimensions
# Resize if necessary, or check dimensions manually
if nuclei_segmentation_mask.shape != ome_image.shape[1:]:
    print("Warning: Segmentation mask and OME image dimensions do not match!")

# Define background label
background_label = 0

# Create a mask to exclude regions corresponding to nuclei (binary mask: 0 outside, 1 inside nuclei)
nuclei_mask = nuclei_segmentation_mask != background_label

# Apply the mask to each channel in the OME-TIFF image (assuming it has multiple channels)
for channel in range(ome_image.shape[0]):  
    ome_image[channel][nuclei_mask] = background_label  # Mask nuclei regions by setting to 0

# Save the modified OME-TIFF image
modified_ome_tiff_path = 'morphology_focus/modified_morphology_focus_0002.ome.tif'
tiff.imwrite(modified_ome_tiff_path, ome_image, photometric='minisblack')

# Optional: Verify the result by printing or visualizing
print(f"Modified OME Image saved at: {modified_ome_tiff_path}")

