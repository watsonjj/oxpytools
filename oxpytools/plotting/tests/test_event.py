from ...io.files import InputFile
from ...utils.datasets import get_oxpytools_extra_path
from ..event import EventPlotter
import numpy as np


def test_eventplotter():
    dataset = get_oxpytools_extra_path("target_sky.fits")
    file = InputFile(dataset)
    source = file.read()
    event = next(source)
    data = event.dl0.tel[0].adc_samples[0]
    plotter = EventPlotter(event)

    camera = plotter.draw_camera(0, data[:, 0])
    assert camera is not None
    np.testing.assert_array_equal(camera.image, data[:, 0])

    plotter.draw_camera_pixel_ids(0)

    waveform = plotter.draw_waveform(data[0, :])
    assert waveform is not None
    np.testing.assert_array_equal(waveform.get_ydata(), data[0, :])

    line = plotter.draw_waveform_positionline(0)
    assert line is not None
    np.testing.assert_array_equal(line.get_xdata(), [0, 0])


