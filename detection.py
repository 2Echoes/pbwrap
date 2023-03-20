import bigfish.stack as stack
import bigfish.multistack as multistack
import bigfish.plot as plot
import bigfish.detection as detection
import numpy as np

def spot_decomposition_nobckgrndrmv(image, spots, spot_radius, voxel_size_nm, alpha= 0.5, beta= 1):
    """ Basically same function as bigfish.detection.decompose_dense but without the remove background gaussian.
    
    Detect dense and bright regions with potential clustered spots and
    simulate a more realistic number of spots in these regions.

    #. We build a reference spot by aggregating predetected spots.
    #. We fit gaussian parameters on the reference spots.
    #. We detect dense regions to decompose.
    #. We simulate as many gaussians as possible in the candidate regions.

    Parameters
    ----------
    image : np.ndarray
        Image with shape (z, y, x) or (y, x).
    spots : np.ndarray
        Coordinate of the spots with shape (nb_spots, 3) or (nb_spots, 2)
        for 3-d or 2-d images respectively.
    voxel_size : int, float, Tuple(int, float) or List(int, float)
        Size of a voxel, in nanometer. One value per spatial dimension (zyx or
        yx dimensions). If it's a scalar, the same value is applied to every
        dimensions.
    spot_radius : int, float, Tuple(int, float) or List(int, float)
        Radius of the spot, in nanometer. One value per spatial dimension (zyx
        or yx dimensions). If it's a scalar, the same radius is applied to
        every dimensions.
    kernel_size : int, float, Tuple(float, int), List(float, int) or None
        Standard deviation used for the gaussian kernel (one for each
        dimension), in pixel. If it's a scalar, the same standard deviation is
        applied to every dimensions. If None, we estimate the kernel size from
        'spot_radius', 'voxel_size' and 'gamma'
    alpha : int or float
        Intensity percentile used to compute the reference spot, between 0
        and 1. The higher, the brighter are the spots simulated in the dense
        regions. Consequently, a high intensity score reduces the number of
        spots added. Default is 0.5, meaning the reference spot considered is
        the median spot.
    beta : int or float
        Multiplicative factor for the intensity threshold of a dense region.
        Default is 1. Threshold is computed with the formula:

        .. math::
            \\mbox{threshold} = \\beta * \\mbox{max(median spot)}

        With :math:`\\mbox{median spot}` the median value of all detected spot
        signals.
    gamma : int or float
        Multiplicative factor use to compute the gaussian kernel size:

        .. math::
            \\mbox{kernel size} = \\frac{\\gamma * \\mbox{spot radius}}{\\mbox{
            voxel size}}

        We perform a large gaussian filter with such scale to estimate image
        background and remove it from original image. A large gamma increases
        the scale of the gaussian filter and smooth the estimated background.
        To decompose very large bright areas, a larger gamma should be set.

    Notes
    -----
    If ``gamma = 0`` and ``kernel_size = None``, image is not denoised.

    Returns
    -------
    spots_postdecom : np.ndarray
        Coordinate of the spots detected, with shape (nb_spots, 3) or
        (nb_spots, 2). One coordinate per dimension (zyx or yx coordinates).

    """
    ndim = image.ndim


    #Spot decomposition    

    # reference spot
    reference_spot = detection.build_reference_spot(
        image=image,
        spots=spots,
        voxel_size=voxel_size_nm, 
        spot_radius= spot_radius,
        alpha=alpha)

    # fit a gaussian function on the reference spot
    sigma_z, sigma_yx, amplitude, background = detection.modelize_spot(
        reference_spot=reference_spot, 
        voxel_size= voxel_size_nm, 
        spot_radius= spot_radius)

    # detect dense regions
    regions_to_decompose, spots_out_regions, region_size = detection.get_dense_region(
        image= image, 
        spots= spots,
        voxel_size=voxel_size_nm,
        spot_radius= spot_radius,
        beta= beta)

    # precompute gaussian function values
    max_grid = max(200, region_size + 1)
    precomputed_gaussian = detection.precompute_erf(
        ndim= ndim,
        voxel_size=voxel_size_nm,
        sigma=(sigma_z, sigma_yx, sigma_yx),
        max_grid=max_grid)

    # simulate gaussian mixtures
    spots_in_regions, _ = detection.simulate_gaussian_mixture(
        image= image,
        candidate_regions=regions_to_decompose,
        voxel_size= voxel_size_nm,
        sigma=(sigma_z, sigma_yx, sigma_yx),
        amplitude=amplitude,
        background=background,
        precomputed_gaussian=precomputed_gaussian)

    spots_postdecomp = np.concatenate((spots_out_regions, spots_in_regions[:, :3]), axis=0)

    return spots_postdecomp

"""

## Testing
image_input = "/Users/floric.slimani/Documents/Projets/1_PBodies/RNA-grouped_3D-z-stacked/input/"
output = "/Users/floric.slimani/Documents/Projets/1_PBodies/RNA-grouped_3D-z-stacked/output/quantification/"
segmentation_input = "/Users/floric.slimani/Documents/Projets/1_PBodies/RNA-grouped_3D-z-stacked/seg_res/"

#Input
dapi = stack.read_image(image_input + "r02c06-SPDL1f01-DAPI-sk1fk1fl1.tiff")
cy3 = stack.read_image(image_input + "r02c06-SPDL1f01-Cy3-sk1fk1fl1.tiff")
alexa647 = stack.read_image(image_input + "r02c06-SPDL1f01-Alexa 647-sk1fk1fl1.tiff")
cell_label = stack.read_array(segmentation_input + "cytoplasm_2Dlabel.npy")
nuc_label = stack.read_array(segmentation_input + "nucleus_3Dlabel.npy")
rna_spots = stack.read_array(segmentation_input + "rnaspots.npy")
malat1_spots = stack.read_array(segmentation_input + "malat1spots.npy")

#Parameters
ndim = 3
voxel_sz = (350,103,103)
#spot decomp
alpha=0.8  # alpha impacts the number of spots per candidate region
beta=5  # beta impacts the number of candidate regions to decompose
gamma=2 # gamma the filtering step to denoise the image
rna_spot_radius=(350, 150, 150)
malat1_spot_radius = (350,100,100)

#focus proj
cy3_proj = stack.focus_projection(cy3)
alexa647_proj = stack.focus_projection(alexa647)


#spots decomp
cy3_spots_postdecomp = spot_decomposition_nobckgrndrmv(cy3, spots= rna_spots, spot_radius= rna_spot_radius, voxel_size_nm= voxel_sz, alpha= alpha, beta= beta)
alexa647_spots_postdecomp = spot_decomposition_nobckgrndrmv(alexa647, spots= malat1_spots, spot_radius= malat1_spot_radius, voxel_size_nm= voxel_sz, alpha= alpha, beta= beta)

cy3_spots_postclustering, cy3_clusters = detection.detect_clusters(spots= cy3_spots_postdecomp, voxel_size= voxel_sz)
alexa647_spots_postclustering, alexa647_clusters = detection.detect_clusters(spots= alexa647_spots_postdecomp, radius = 500, voxel_size= voxel_sz)

print("detected spots after clustering")
print("\r shape: {0}".format(alexa647_spots_postclustering.shape))
print("\r dtype: {0}".format(alexa647_spots_postclustering.dtype), "\n")
print("detected clusters")
print("\r shape: {0}".format(alexa647_clusters.shape))
print("\r dtype: {0}".format(alexa647_clusters.dtype))



plot.plot_detection(alexa647_proj, 
                    spots=[alexa647_spots_postdecomp, alexa647_clusters[:, :3]], 
                    shape=["circle", "polygon"], 
                    radius=[3, 6], 
                    color=["red", "blue"],
                    linewidth=[1, 2], 
                    fill=[False, True], 
                    contrast=True)




spots_no_ts, foci, ts = multistack.remove_transcription_site(rna_spots, clusters, nuc_label, ndim=3)




fov_results = multistack.extract_cell(cell_label=cell_label[5], ndim=ndim, nuc_label= nuc_label[5], rna_coord= rna_spots, others_coord= other_coords)
for i, cell_results in enumerate(fov_results):
    print("cell {0}".format(i))
    
    # get cell results
    cell_mask = cell_results["cell_mask"]
    cell_coord = cell_results["cell_coord"]
    nuc_mask = cell_results["nuc_mask"]
    nuc_coord = cell_results["nuc_coord"]
    rna_coord = cell_results["rna_coord"]
    #foci_coord = cell_results["foci"]
    #ts_coord = cell_results["transcription_site"]
    #image_contrasted = cell_results["image"]
    print("\r number of rna {0}".format(len(rna_coord)))
    #print("\r number of foci {0}".format(len(foci_coord)))
    #print("\r number of transcription sites {0}".format(len(ts_coord)))
    
    # plot cell
    plot.plot_cell(
        ndim=3, cell_coord=cell_coord, nuc_coord=nuc_coord, 
        rna_coord=rna_coord, cell_mask=cell_mask, nuc_mask=nuc_mask, 
        title="Cell {0}".format(i))

 """
