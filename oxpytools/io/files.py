"""Definitions of file-related classes.

The classes defined in this module are to hold files the will be read or wrote.

"""

from pathlib import Path
from os.path import basename, splitext, dirname
from astropy import log

from .targetio import targetio_event_source
from ctapipe.utils.datasets import get_path
from ctapipe.io.hessio import hessio_event_source


class InputFile:
    """
    Class to handle input files

    The input files intended for this class are either sim_telarray .gz files,
    or target .fits files.

    Attributes
    ----------
    input_path : str
    directory : str
        Automatically set from `input_path`.
    filename : str
        Name of the file without the extension.
        Automatically set from `input_path`.
    extension : str
        Automatically set from `input_path`.
    type : {'simtel', 'target'}
        The type of file, related to its source.
        Automatically set from `input_path`.

    """

    def __init__(self, input_path):
        """
        Parameters
        ----------
        input_path : str
            Full path to the file

        """
        self.__input_path = None
        self.directory = None
        self.filename = None
        self.extension = None
        self.type = None

        self.input_path = input_path

        log.info("[file] {}".format(self.input_path))
        log.info("[file][type] {}".format(self.type))

    @property
    def input_path(self):
        return self.__input_path

    @input_path.setter
    def input_path(self, string):
        path = Path(string)
        try:
            if not path.exists():
                raise FileNotFoundError
        except FileNotFoundError as e:
            log.exception("file path does not exist: '{}'".format(string))

        self.__input_path = path.as_posix()
        self.directory = dirname(self.__input_path)
        self.filename = splitext(basename(self.__input_path))[0]
        self.extension = splitext(self.__input_path)[1]

        try:
            if self.extension == ".gz":
                self.type = "hessio"
            elif self.extension == ".fits":
                self.type = "targetio"
            else:
                raise RuntimeError()
        except RuntimeError as e:
            log.exception("unknown file extension '{}'".format(self.__input_path))

    def read(self):
        """
        Read the file using the appropriate method depending on the file type

        Returns
        -------
        source : generator
            A generator that can be iterated over to obtain events

        """

        try:
            if self.type == "hessio":
                return hessio_event_source(get_path(self.input_path))
            elif self.type == "targetio":
                return targetio_event_source(self.input_path)
            else:
                raise RuntimeError()
        except RuntimeError as e:
            log.exception("unknown file type '{}'".format(self.type))
