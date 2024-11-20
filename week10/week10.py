#!/usr/bin/env python

import numpy as np
import scipy
import matplotlib.pyplot as plt
import imageio
import plotly.express as px
import plotly
import pandas as pd

# helps to reshape things easily 
na = np.newaxis


### Exercise 1 ### 

# test sample image to figure out size
test_img = imageio.v3.imread("~/qbb2024-answers/week10/APEX1_field0_DAPI.tif").astype(np.uint16)

# make dictionary called images
images = {}

# specify variables for loop
genename = ["APEX1", "PIM2", "POLR2B", "SRSF1"]
fields = ["field0", "field1"]
channels = ["DAPI", "PCNA", "nascentRNA"] 

# deleted DAPI_c images because not supposed to be there

# loop through images to load in
for name in genename:
     for field in fields: 
        for channel in channels: 
            file_path = f"~/qbb2024-answers/week10/{name}_{field}_{channel}.tif"
            try:
                img = imageio.v3.imread(file_path).astype(np.uint16)
                images[f"{name}_{field}_{channel}"] = img
                # print(f"{name}_{field}_{channel} processed")
            except FileNotFoundError:
                continue 
                # print(f"{name}_{field}_{channel} not found")

# make list of image arrays
all_img_arrays = []

# create numpy array to load images into
for name in genename:
    for field in fields: 
        img_array =  np.zeros((test_img.shape[0], test_img.shape[1], 3), np.uint16)
        for i, channel in enumerate(channels):
            try: 
                img_array[:, :, i] = images[f"{name}_{field}_{channel}"]
                img_array[:, :, i] -= np.amin(img_array[:, :, i])
                img_array[:, :, i] /= np.amax(img_array[:, :, i])
            except:continue
        # print plots to show array worked
        # plt.imshow(img_array)
        # plt.show()
        all_img_arrays.append(img_array)

all_img_arrays = np.array(all_img_arrays)
            
### Exercise 2 ###

## 2.1 

# make mask where each image will scale to the corresponding DAPI mean
binary_mask = []
for img in all_img_arrays:
    DAPI_channel = img[:, :, 0]
    DAPI_mean = np.mean(DAPI_channel)
    DAPI_mask = DAPI_channel >= DAPI_mean
    binary_mask.append(DAPI_mask)

# test to make sure it worked for first image
# plt.imshow(mask[0])
# plt.show()

## 2.2

# from live coding, find labels based on dapi mask
# Let's use a function to find nuclei bodies
def find_labels(mask):
    # Set initial label
    l = 0
    # Create array to hold labels
    labels = np.zeros(mask.shape, np.int32)
    # Create list to keep track of label associations
    equivalence = [0]
    # Check upper-left corner
    if mask[0, 0]:
        l += 1
        equivalence.append(l)
        labels[0, 0] = l
    # For each non-zero column in row 0, check back pixel label
    for y in range(1, mask.shape[1]):
        if mask[0, y]:
            if mask[0, y - 1]:
                # If back pixel has a label, use same label
                labels[0, y] = equivalence[labels[0, y - 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[0, y] = l
    # For each non-zero row
    for x in range(1, mask.shape[0]):
        # Check left-most column, up  and up-right pixels
        if mask[x, 0]:
            if mask[x - 1, 0]:
                # If up pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 0]]
            elif mask[x - 1, 1]:
                # If up-right pixel has label, use that label
                labels[x, 0] = equivalence[labels[x - 1, 1]]
            else:
                # Add new label
                l += 1
                equivalence.append(l)
                labels[x, 0] = l
        # For each non-zero column except last in nonzero rows, check up, up-right, up-right, up-left, left pixels
        for y in range(1, mask.shape[1] - 1):
            if mask[x, y]:
                if mask[x - 1, y]:
                    # If up pixel has label, use that label
                    labels[x, y] = equivalence[labels[x - 1, y]]
                elif mask[x - 1, y + 1]:
                    # If not up but up-right pixel has label, need to update equivalence table
                    if mask[x - 1, y - 1]:
                        # If up-left pixel has label, relabel up-right equivalence, up-left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x - 1, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x - 1, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    elif mask[x, y - 1]:
                        # If left pixel has label, relabel up-right equivalence, left equivalence, and self with smallest label
                        labels[x, y] = min(equivalence[labels[x, y - 1]], equivalence[labels[x - 1, y + 1]])
                        equivalence[labels[x, y - 1]] = labels[x, y]
                        equivalence[labels[x - 1, y + 1]] = labels[x, y]
                    else:
                        # If neither up-left or left pixels are labeled, use up-right equivalence label
                        labels[x, y] = equivalence[labels[x - 1, y + 1]]
                elif mask[x - 1, y - 1]:
                    # If not up, or up-right pixels have labels but up-left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x - 1, y - 1]]
                elif mask[x, y - 1]:
                    # If not up, up-right, or up-left pixels have labels but left does, use that equivalence label
                    labels[x, y] = equivalence[labels[x, y - 1]]
                else:
                    # Otherwise, add new label
                    l += 1
                    equivalence.append(l)
                    labels[x, y] = l
        # Check last pixel in row
        if mask[x, -1]:
            if mask[x - 1, -1]:
                # if up pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -1]]
            elif mask[x - 1, -2]:
                # if not up but up-left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x - 1, -2]]
            elif mask[x, -2]:
                # if not up or up-left but left pixel is labeled use that equivalence label 
                labels[x, -1] = equivalence[labels[x, -2]]
            else:
                # Otherwise, add new label
                l += 1
                equivalence.append(l)
                labels[x, -1] = l
    equivalence = np.array(equivalence)
    # Go backwards through all labels
    for i in range(1, len(equivalence))[::-1]:
        # Convert labels to the lowest value in the set associated with a single object
        labels[np.where(labels == i)] = equivalence[i]
    # Get set of unique labels
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[np.where(labels == j)] = i
    return labels

# to check it worked
# labels = find_labels(binary_mask[0])
# plt.imshow(labels)
# plt.show()

# this is super long
# for mask in binary_mask:
#     labels = find_labels(mask)
#     plt.imshow(labels)
#     plt.show()

# create label array
label_array = []
for mask in binary_mask:
    labels = find_labels(mask)
    label_array.append(labels)
# print(type(label_array))
# plt.imshow(label_array[0])
# plt.show

# make sure label_array is 8, when DAPI_c was there it was 10 and messing up all the code
# print(len(label_array))
# print(len(all_img_arrays))

## 2.3

# filter by size, from live code 
def filter_by_size(labels, minsize = 100, maxsize = 42000000000):
    # Find label sizes
    sizes = np.bincount(labels.ravel())
    # Iterate through labels, skipping background
    for i in range(1, sizes.shape[0]):
        # If the number of pixels falls outsize the cutoff range, relabel as background
        if sizes[i] < minsize or sizes[i] > maxsize:
            # Find all pixels for label
            where = np.where(labels == i)
            labels[where] = 0
    # Get set of unique labels
    ulabels = np.unique(labels)
    for i, j in enumerate(ulabels):
        # Relabel so labels span 1 to # of labels
        labels[np.where(labels == j)] = i
    return labels


# background label is 0, adjust for more contrast 
label_copy = np.copy(labels)
label_copy[np.where(label_copy == 0)] -= 50
# plt.imshow(label_copy)
# plt.show()

# save to array
for i in range(len(label_array)):
    label_array[i] = filter_by_size(label_array[i])

def filter_by_size_meanstdev(labels):
    sizes = np.bincount(labels.ravel())
    # skip background sizes which are zero
    non_zero_sizes = sizes[1:] 
    mean_size = np.mean(non_zero_sizes)
    std_size = np.std(non_zero_sizes)
    min_size = mean_size - std_size
    max_size = mean_size + std_size
    for i in range(1, sizes.shape[0]):
        if sizes[i] < min_size or sizes[i] > max_size:
            where = np.where(labels == i)
            labels[where] = 0
        
### Exercise 3 ###

# 3.1

# find mean signal for each nucleus from the PCNA and nascent RNA channels

# make lists to store information
genes_ls = []
fields_ls = []
nuclei_ls = []
PCNA_signal = []
nRNA_signal = []
log2_ratio = []

# for loop to go through each label for each array
for genefield in range(len(label_array[:])):
    for nucleus in range(np.amax(label_array[genefield]+1)):
        if nucleus == 0: 
            continue
        else: 
            # label with gene names
            where = np.where(label_array[genefield] == nucleus)
            # APEX1 is first two in list
            if genefield in [0,1]:
                genes_ls.append("APEX1")
            # then are the two PIM2 images, elif = else if
            elif genefield in [2,3]:
                genes_ls.append("PIM2")
            # then two POLR2B images
            elif genefield in [4,5]:
                genes_ls.append("POLR2B")
            # then two SRSF1 images
            elif genefield in [6,7]:
                genes_ls.append("SRSF1")
            if genefield in [0,2,4,6]: fields_ls.append(0)
            else: fields_ls.append(1)
            nuclei_ls.append(nucleus)
            PCNA_signal.append(np.average(all_img_arrays[genefield][where][1]))
            nRNA_signal.append(np.average(all_img_arrays[genefield][where][2]))
            log2_ratio.append(np.log2(np.average(all_img_arrays[genefield][where][2])/np.average(all_img_arrays[genefield][where][1])))

# save as csv
pd.DataFrame({'gene':genes_ls,'viewField':fields_ls,'nucleusIDX':nuclei_ls,'nascentRNA':nRNA_signal,'PCNA':PCNA_signal,'Ratio':log2_ratio}).to_csv("mean_nucleus_signal.csv",sep=",",mode="w",header=True,index=False)


