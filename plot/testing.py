import skimage.data as data
import bigfish.plot as plot
import bigfish.stack as stack
import bigfish.segmentation as seg
import matplotlib.pyplot as plt
from pbwrap.plot import plot_labels


image = data.checkerboard()
image = seg.thresholding(image, threshold= 200)
label = seg.label_instances(image)
print(label.max())


plot.plot_images([label])
plot_labels(label)
