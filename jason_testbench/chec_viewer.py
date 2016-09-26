from ctapipe.io import CameraGeometry
from ctapipe.io.files import InputFile
from ctapipe.visualization import CameraDisplay

from matplotlib import pyplot as plt
import numpy as np

fp = '/Users/Jason/Software/outputs/libCHEC/sky/cameradata_run1594.fits'
input_file = InputFile(fp, 'targetio')
telid = 0

fig = plt.figure(figsize=(15, 7))
ax0 = fig.add_subplot(1, 2, 1)
ax1 = fig.add_subplot(2, 2, 2)
ax2 = fig.add_subplot(2, 2, 4)
fig.subplots_adjust(hspace=.5)

# Build CameraPlotter from first CHEC event
first_event = input_file.get_event(0)
data = first_event.dl0.tel[telid].adc_samples[0]
ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("Amplitude (ADC)")
title1 = ax1.set_title("Waveform (pixel: {})".format(0))
w1, = ax1.plot(data[0, :])
ax2.set_xlabel("Time (ns)")
ax2.set_ylabel("Amplitude (ADC)")
title2 = ax2.set_title("Waveform (pixel: {})".format(0))
w2, = ax2.plot(data[0, :])
geom = CameraGeometry.guess(*first_event.meta.pixel_pos[telid],
                            first_event.meta.optical_foclen[telid])
camera = CameraDisplay(geom, ax=ax0)
camera.cmap = plt.cm.viridis
camera.add_colorbar(ax=ax0, label="Amplitude (ADC)")
title0 = ax0.set_title("")

while True:
    source = input_file.read()
    for event in source:
        data = event.dl0.tel[telid].adc_samples[0]
        data_min = np.min(data[np.nonzero(data)])
        data_max = np.max(data)
        camera.set_limits_minmax(data_min * 1.05,
                                 data_max * 0.95)
        camera.image = data[:, 0]
        title0.set_text("Event: {}".format(event.count))

        max_charges = np.max(data, axis=1)
        max_pixel = int(np.argmax(max_charges))
        min_pixel = int(np.argmin(max_charges[np.nonzero(max_charges)]))
        w1.set_ydata(data[max_pixel, :])
        ax1.relim()  # make sure all the data fits
        ax1.autoscale()
        title1.set_text("Waveform (pixel: {})".format(max_pixel))
        w2.set_ydata(data[min_pixel, :])
        ax2.relim()  # make sure all the data fits
        ax2.autoscale()
        title2.set_text("Waveform (pixel: {})".format(min_pixel))
        print(event.count)

        plt.pause(0.00001)
