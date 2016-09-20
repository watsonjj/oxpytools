from ctapipe.utils.datasets import get_path
from ctapipe.io.hessio import hessio_event_source
from ctapipe.calib.camera.calibrators import calibrate_event, calibrate_source
from ctapipe.calib.camera.integrators import integrator_dict
import numpy as np
from astropy import log


def get_test_parameters():
    parameters = {"integrator": "local_peak_integration",
                  "window": 7,
                  "shift": 3,
                  "sigamp": [2, 4],
                  "lwt": 0,
                  "calib_scale": 1.05}
    return parameters
#
# log.setLevel(10)
# params = get_test_parameters()
# filename = get_path('gamma_test.simtel.gz')
# source = hessio_event_source(filename)
# for event in source:
#     if event.dl0.event_id == 409:
#         telid = 11
#         for i in range(5):
#
#             integrator = i+1
#             integrator_name = integrator_dict()[0][integrator]
#             params['integrator'] = integrator_name
#             calibrated_event = calibrate_event(event, params)
#             peakpos = calibrated_event.dl1.tel[telid].peakpos
#             print(integrator_name, peakpos)
#         break


from ctapipe.analysis.camera.tests.test_chargeresolution import *

test_add_charges()
test_get_charge_resolution()
test_get_binned_charge_resolution()
test_limit_curves()
test_requirement()
test_goal()