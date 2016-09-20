"""Module to handle to storage of events extracted with `targetio` into
containers defined in `ctapipe.io.containers`.
"""

import target_io
from astropy import log
import numpy as np

from .pixels import get_pixel_id, get_all_pixel_pos

from ctapipe.io.containers import RawData
from ctapipe.io.containers import RawCameraData
from ctapipe.core import Container

from astropy import units as u


def get_targetio_event(reader, event_index):
    """
    Get all the adc_samples for an event

    Parameters
    ----------
    reader : target_io.EventFileReader()
        `target_io` object created when reading target files
    event_index : int

    Returns
    -------
    event : array-like
        An array of size [npix,nsamples] containing the adc_samples for each
        pixel and timeslice

    """

    # Get the sample size of the waveforms
    event_packet = reader.GetEventPacket(0, 0)
    packet = target_io.DataPacket()
    packet.Assign(event_packet, reader.GetPacketSize())
    wav = packet.GetWaveform(0)
    num_samples = wav.GetSamples()

    # unsorted_adc_samples = []
    # pixel_ids = []
    event = np.zeros((2048, num_samples))

    for ipack in range(reader.GetNPacketsPerEvent()):
        event_packet = reader.GetEventPacket(event_index, ipack)
        packet = target_io.DataPacket()
        packet.Assign(event_packet, reader.GetPacketSize())
        module = packet.GetSlotID()
        for iwav in range(packet.GetNumberOfWaveforms()):
            wav = packet.GetWaveform(iwav)
            asic = wav.GetASIC()
            channel = wav.GetChannel()
            pixel_id = get_pixel_id(module, asic, channel)
            for isam in range(wav.GetSamples()):
                event[pixel_id, isam] = wav.GetADC(isam)

    return event


def targetio_event_source(url, max_events=None):
    """A generator that streams data from a GCT target file

    Parameters
    ----------
    url : str
        filepath to the target file

    Returns
    -------
    container
        object containing entire single event information

    """

    reader = target_io.EventFileReader(url)
    n_events = reader.GetNEvents()
    log.info("[run][n_events] {}".format(n_events))

    counter = 0
    container = Container("targetio_container")
    container.meta.add_item('targetio_input', url)
    container.meta.add_item('pixel_pos', dict())
    container.meta.add_item('optical_foclen', dict())
    container.add_item("dl0", RawData())
    container.add_item("count")

    tel_id = 0

    for event_index in range(n_events):
        if counter > max_events:
            break

        container.dl0.run_id = 0
        container.dl0.event_id = event_index
        container.dl0.tels_with_data = {tel_id}

        container.count = counter

        container.dl0.tel = dict()  # clear the previous telescopes

        container.meta.pixel_pos[tel_id] = get_all_pixel_pos() * u.m
        container.meta.optical_foclen[tel_id] = 2.283 * u.m

        nchans = 1
        container.dl0.tel[tel_id] = RawCameraData(tel_id)
        container.dl0.tel[tel_id].num_channels = 1

        container.dl0.tel[tel_id].adc_samples[0] = \
            get_targetio_event(reader, event_index)

        yield container
        counter += 1
