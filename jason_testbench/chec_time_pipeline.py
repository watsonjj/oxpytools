import argparse
import numpy as np
import target_calib
from astropy import log
from ctapipe.calib.camera.integrators import local_peak_integration
from targetpipe.io.files import CHECInputFile as InputFile
from targetpipe.io.pixels import Pixels
from targetpipe.utils.plotting import intensity_to_hex

from IPython import embed
from tqdm import tqdm


params = dict(integration_window=[7, 3])

input_path = '/Users/Jason/Software/outputs/lab/ped_wfs.fits'
input_file = InputFile(input_path, max_events=1000)
n_events = input_file.num_events
telid = 0

first_event = input_file.get_event(0)
n_pix = first_event.dl0.tel[0].adc_samples[telid].shape[0]
n_samples = first_event.dl0.tel[0].adc_samples[telid].shape[1]

ped_obj = target_calib.Pedestal(32)
pixels = Pixels()
pix2tmtmpix = pixels.pix_to_tm_tmpix.astype(np.uint16)

def generate_pedestal():
    source = input_file.read()
    desc = "Filling pedestal"
    with tqdm(total=n_events, desc=desc) as pbar:
        for event in source:
            pbar.update(1)
            first_cell_ids = event.dl0.tel[telid].first_cell_ids
            waveforms = event.dl0.tel[telid].adc_samples[0]
            ped_obj.AddEvent(pix2tmtmpix, waveforms, first_cell_ids)
    return ped_obj
#Initial: 30Hz
#Update to new hessio layout: 35Hz
#Reduced to 96 samples: 36Hz
#(IO)Loop pixels in c++: 73Hz
#(Add Ped)Loop pixels in c++: 280Hz

def calibrate(ped_obj):
    source = input_file.read()

    pedsub = np.empty((n_pix, n_samples), dtype=np.float32)

    with tqdm(total=n_events, desc='Calibrating') as pbar:
        for event in source:
            pbar.update(1)
            waveforms = event.dl0.tel[telid].adc_samples[0]
            fci = event.dl0.tel[telid].first_cell_ids

            pedsub[:] = waveforms[:]
            if ped_obj:
                ped_obj.ApplyEvent(pix2tmtmpix, waveforms, fci, pedsub)

            # TODO: Transfer function

            # TODO: Proper extraction (using components)
            charge = local_peak_integration(pedsub[None, ...], params)

            # TODO: Convert to pe?
#Initial: 28Hz
#Update to new hessio layout: 32Hz
#Reduced to 96 samples: 32Hz
#(IO)Loop pixels in c++: 60Hz
#(Sub Ped)Loop pixels in c++: 170Hz

def colours():
    source = input_file.read()

    pedsub = np.empty((n_pix, n_samples), dtype=np.float32)

    with tqdm(total=n_events, desc='Calibrating') as pbar:
        for event in source:
            pbar.update(1)
            waveforms = event.dl0.tel[telid].adc_samples[0]
            # h = intensity_to_hex(waveforms[:, 0])


embed()
