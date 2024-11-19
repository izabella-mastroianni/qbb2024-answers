#!/usr/bin/env python

import numpy as np
import scipy
import matplotlib.pyplot as plt
import imageio
import plotly.express as px
import plotly

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
channels = ["DAPI", "PCNA", "nascentRNA", "DAPI_c"] 
# DAPI_c only appears for a few genes, will be ignored by the code later

# loop through images to load in
for name in genename:
     for field in fields: 
        for channel in channels: 
            file_path = f"~/qbb2024-answers/week10/{name}_{field}_{channel}.tif"
            try:
                img = imageio.v3.imread(file_path).astype(np.uint16)
                images[f"{name}_{field}_{channel}"] = img
                print(f"{name}_{field}_{channel} processed")
            except FileNotFoundError:
                print(f"{name}_{field}_{channel} not found")

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

