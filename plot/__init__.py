# -*- coding: utf-8 -*-
# Author: Floric Slimani <floric.slimani@live.fr>


"""
This pbwrap subpackage groups custom plots tuned to our analysis pipelines.
"""

from .control_plots import plot_labels, plot_detection_steps, plot_cell
from .control_plots import save_plot
from .visuals import output_spot_tiffvisual, nucleus_signal_control, dapi_artifact_overlay, colocalisation_plot
from .utils import get_colors_list, hist_maximum, plot_horizontal_bar, make_color_frame
from .histogram import histogram