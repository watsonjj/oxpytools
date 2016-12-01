from ctapipe.core import Tool
from traitlets import Dict, List

from ctapipe.io.files import InputFile
from ctapipe.io.hessio import hessio_event_source
from ctapipe.utils.datasets import get_path
from ctapipe.io import CameraGeometry
import numpy as np

from ctapipe.calib.camera.charge_extraction import *


def get_test_parameters():
    parameters = {"integrator": "nb_peak_integration",
                  "integration_window": [7, 3],
                  "integration_sigamp": [2, 4],
                  "integration_lwt": 0}
    return parameters


def get_test_event():
    filename = get_path('gamma_test.simtel.gz')
    for event in hessio_event_source(filename):
        if event.dl0.event_id == 409:
            return event

def test_full_integration():
    telid = 11
    event = get_test_event()
    nsamples = event.dl0.tel[telid].num_samples
    data = np.array(
        list(event.dl0.tel[telid].adc_samples.values()))
    ped = event.dl0.tel[telid].pedestal
    data_ped = data - np.atleast_3d(ped / nsamples)
    data_ped = np.array(
        [data_ped[0], data_ped[0]])  # Test LG functionality

    integrator = FullIntegrator(data_ped, None)
    integration, window, peakpos = integrator.extract_charge()

    assert integration[0][0] == 149
    assert integration[1][0] == 149
    assert peakpos[0] is None
    assert peakpos[1] is None

test_full_integration()

# telid = 11
# waveforms = get_test_event()
# nsamples = waveforms.dl0.tel[telid].n_samples
# data = np.array(list(waveforms.dl0.tel[telid].adc_samples.values()))
# ped = waveforms.dl0.tel[telid].pedestal
# data_ped = data - np.atleast_3d(ped / nsamples)
# data_ped = np.array([data_ped[0], data_ped[0]])
#
# pixel_pos = waveforms.meta.pixel_pos[telid]
# optical_foclen = waveforms.meta.optical_foclen[telid]
#
# geom = CameraGeometry.guess(*pixel_pos, optical_foclen)
# nei = geom.neighbors
#
# params = get_test_parameters()
#
#
# class MyTool(Tool):
#     name = "mytool"
#     description = "do some things and stuff"
#     aliases = Dict(dict('ChargeExtractorFactory.extractor'))
#
#     # Which classes are registered for configuration
#     classes = List([ChargeExtractorFactory])
#
#     # local configuration parameters
#     # iterations = Integer(5, help="Number of times to run",
#     #                      allow_none=False).tag(config=True)
#
#     def setup_comp(self):
#         self.comp = MyComponent(self, config=self.config)
#         self.comp2 = SecondaryMyComponent(self, config=self.config)
#
#     def setup_advanced(self):
#         self.advanced = AdvancedComponent(self, config=self.config)
#
#     def setup(self):
#         factory = ChargeExtractorFactory()
#         self.setup_comp()
#         self.setup_advanced()
#
#     def start(self):
#         self.log.info("Performing {} iterations...".format(self.iterations))
#         for ii in range(self.iterations):
#             self.log.info("ITERATION {}".format(ii))
#             self.comp.do_thing()
#             self.comp2.do_thing()
#             sleep(0.5)
#
#     def finish(self):
#         self.log.warning("Shutting down.")