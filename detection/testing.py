import bigfish.stack as stack
import pbwrap.detection.detection_wrappers as detection
import os
import time

path = '/media/floricslimani/SSD 4To/SSD_floricslimani/4_centrosome/test_input'
if not path.endswith('/') : path += '/'
files = os.listdir(path)

voxel_size = (300, 103, 103)
spot_size = (150,100,100)


images = [stack.read_image(path + file) for file in files if file.endswith('.tiff')]


clock = time.process_time()
threshold = detection.compute_auto_threshold(images, voxel_size, spot_radius=spot_size, im_number= 5)
print("time for threshold computing with shuffle : {0}".format(time.process_time() - clock))
clock = time.process_time()
# _,bigfish_threshold = detection.detect_spots(images, ndim= 3, threshold_penalty=1, voxel_size=voxel_size, spot_radius=spot_size, return_threshold=True)
bigfish_threshold = 728.0
print("best time for bigFish : {0}".format(time.process_time() - clock))

print("Newly computed threshold : ", threshold)
print("Big computed threshold : ", bigfish_threshold)