from astropy import log
from ctapipe.io.files import InputFile
from ctapipe.io.hessio import hessio_event_source
from ctapipe.utils.datasets import get_path
from oxpytools.io.targetio import targetio_event_source


class CHECInputFile(InputFile):
    def __init__(self, input_path, file_origin=None):
        if not file_origin:
            file_origin = 'targetio'
            if input_path.endswith('.gz'):
                file_origin = 'hessio'

        InputFile.__init__(self, input_path, file_origin)

    def read(self, max_events=None):
        """
        Read the file using the appropriate method depending on the file origin

        Parameters
        ----------
        max_events : int
            Maximum number of events to read

        Returns
        -------
        source : generator
            A generator that can be iterated over to obtain events
        """

        # Obtain relevent source
        log.debug("[file] Reading file...")
        if max_events:
            log.info("[file] Max events being read = {}".format(max_events))
        switch = {
            'hessio':
                lambda: hessio_event_source(get_path(self.input_path),
                                            max_events=max_events),
            'targetio':
                lambda: targetio_event_source(self.input_path,
                                              max_events=max_events)
        }
        try:
            source = switch[self.origin]()
        except KeyError:
            log.exception("unknown file origin '{}'".format(self.origin))
            raise
        log.debug("[file] Reading complete")

        return source
