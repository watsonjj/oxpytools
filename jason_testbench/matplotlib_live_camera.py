#!python
"""
Demonstration of how a live waveforms viewer might be implemented. This script
repeatedly loops through a CHEC file and displays the camera image and two
waveforms.
This script is only valid for targetio files.
"""

from targetpipe.io.files import CHECInputFile as InputFile
from ctapipe.io import CameraGeometry
from ctapipe.visualization import CameraDisplay

from matplotlib import pyplot as plt
import numpy as np
import argparse
from astropy import log
from astropy import units as u
from copy import deepcopy
from matplotlib.animation import FuncAnimation


def rate(timeit_object):
    timeperevent = timeit_object.best/num_events
    rate = 1/timeperevent * u.Hz
    return rate


def get_data(inputf):
    while True:
        source = inputf.read()
        for event in source:
            for i in range(128):
                data = event.dl0.tel[telid].adc_samples[0][:,i]
                yield data


def camera_image(input_file, ax):
    first_event = input_file.get_event(0)
    geom = CameraGeometry.guess(*first_event.meta.pixel_pos[telid],
                                first_event.meta.optical_foclen[telid])
    camera = CameraDisplay(geom, ax=ax)
    camera.cmap = plt.cm.viridis
    import time

    def update(iframe, source):
        data = next(source)
        print(iframe)
        # data = np.zeros((2048,128))
        # data[iframe,0] = iframe
        camera.image = data
        return camera.pixels,

    source = get_data(input_file)

    anim = FuncAnimation(fig, update, fargs=[source], blit=True,
                         frames=1000, interval=1, repeat=False)

    plt.show()
    # anim.save('/tmp/animation.gif', writer="imagemagick", fps=15)

if __name__ == '__main__':
    input_path = '/Users/Jason/Software/outputs/libCHEC/sky/cameradata_run1596.fits'
    input_file = InputFile(input_path, 'targetio')
    telid = 0
    num_events = len(input_file.get_list_of_event_ids())

    plt.style.use("ggplot")
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    camera_image(input_file, ax)