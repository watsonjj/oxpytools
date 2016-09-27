from targetpipe.io.targetio import targetio_event_source
from ...utils.datasets import get_oxpytools_extra_path


def test_targetio_event_source():
    dataset = get_oxpytools_extra_path("target_sky.fits")
    source = targetio_event_source(dataset)
    event = next(source)
    assert event.dl0.event_id == 0
    assert event.dl0.tel[0].adc_samples[0][0][0] == 287
    event = next(source)
    assert event.dl0.event_id == 1
    assert event.dl0.tel[0].adc_samples[0][0][0] == 288
