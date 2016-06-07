import numpy as np
from astropy.io import fits
from tqdm import tqdm
from ctapipe.utils.datasets import get_path
from ctapipe.io.hessio import hessio_event_source
import os.path
import math
import sys

filepath = get_path("/Users/Jason/Software/outputs/sim_telarray/meudon_gamma/simtel_runmeudon_gamma_30tel_30deg_1.gz")
# filepath = get_path(sys.argv[1])

events_per_byte = 16673/2150294565

est_events = os.path.getsize(filepath) * events_per_byte
source = hessio_event_source(filepath, max_events=1)
event = next(source)
tel = list(event.dl0.tels_with_data)[0]
pos_arr = np.array(event.meta.pixel_pos[tel])
n_pix = pos_arr.shape[1]
for pix in range(n_pix):
    pos_x = pos_arr[0,pix]
    pos_y = pos_arr[1,pix]
    print("({0:.5g}, {1:.5g}),  # {2}".format(pos_x,-pos_y,pix))
