"""Module to handle to storage of events extracted with `target_io` into
containers defined in `ctapipe.io.containers`.
"""

import target_io
import numpy as np

from .pixels import get_pixel_id, get_all_pixel_pos

from ctapipe.io.containers import RawData
from ctapipe.io.containers import RawCameraData
from ctapipe.core import Container

from astropy import units as u


class TargetioExtractor:
    """
    Extract waveforms from `target_io` and build them into a camera image

    Attributes
    ----------
    reader : target_io.EventFileReader()
    n_events : int
        number of events in the fits file
    num_samples : int
        number of samples in the waveform
    event : ndarray
        two dimensional array to store the waveform for each pixel

    """
    def __init__(self, url):
        """
        Parameters
        ----------
        url : string
            path to the TARGET fits file
        """
        self.__url = None
        self.__event_index = None

        self.reader = None
        self.n_events = None
        self.num_samples = None
        self.event = None

        self.url = url

    @property
    def url(self):
        return self.__url

    @url.setter
    def url(self, string):
        self.__url = string
        if self.reader:
            self.close_reader()
        self.reader = target_io.EventFileReader(string)
        self.n_events = self.reader.GetNEvents()

        # Get the sample size of the waveforms
        event_packet = self.reader.GetEventPacket(0, 0)
        packet = target_io.DataPacket()
        packet.Assign(event_packet, self.reader.GetPacketSize())
        wav = packet.GetWaveform(0)
        self.num_samples = wav.GetSamples()

    @property
    def event_index(self):
        return self.__event_index

    @event_index.setter
    def event_index(self, val):
        self.__event_index = val
        self.event = np.zeros((2048, self.num_samples))
        for ipack in range(self.reader.GetNPacketsPerEvent()):
            event_packet = self.reader.GetEventPacket(val, ipack)
            packet = target_io.DataPacket()
            packet.Assign(event_packet, self.reader.GetPacketSize())
            module = packet.GetSlotID()
            for iwav in range(packet.GetNumberOfWaveforms()):
                wav = packet.GetWaveform(iwav)
                asic = wav.GetASIC()
                channel = wav.GetChannel()
                pixel_id = get_pixel_id(module, asic, channel)
                self.event[pixel_id] = wav.GetADCArray(wav.GetSamples())

    def close_reader(self):
        self.reader.Close()
        del self.reader
        self.reader = None


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

    targetio_extractor = TargetioExtractor(url)
    n_events = targetio_extractor.n_events

    counter = 0
    container = Container("targetio_container")
    container.meta.add_item('targetio_input', url)
    container.meta.add_item('pixel_pos', dict())
    container.meta.add_item('optical_foclen', dict())
    container.meta.add_item('n_events', n_events)
    container.add_item("dl0", RawData())
    container.add_item("count")

    tel_id = 0

    for targetio_extractor.event_index in range(n_events):
        if max_events is not None:
            if counter > max_events:
                break

        container.dl0.run_id = 0
        container.dl0.event_id = targetio_extractor.event_index
        container.dl0.tels_with_data = {tel_id}

        container.count = counter

        container.dl0.tel = dict()  # clear the previous telescopes

        container.meta.pixel_pos[tel_id] = get_all_pixel_pos() * u.m
        container.meta.optical_foclen[tel_id] = 2.283 * u.m

        container.dl0.tel[tel_id] = RawCameraData(tel_id)
        container.dl0.tel[tel_id].num_channels = 1
        container.dl0.tel[tel_id].adc_samples[0] = targetio_extractor.event

        yield container
        counter += 1

    del targetio_extractor

