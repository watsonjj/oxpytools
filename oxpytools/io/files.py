from pathlib import Path
from os.path import basename, splitext, dirname
from astropy import log


class InputFile:
    def __init__(self, input_path):
        self.__input_path = None
        self.directory = None
        self.filename = None
        self.extension = None
        self.type = None

        self.input_path = input_path

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
        self.filename = basename(self.__input_path)
        self.extension = splitext(self.__input_path)[1]

        try:
            if self.extension == ".gz":
                self.type = "simtel"
            elif self.extension == ".fits":
                self.type = "target"
            else:
                raise RuntimeError("unknown file extension '{}'".format(self.__input_path))
        except RuntimeError as e:
            log.exception("unknown file extension '{}'".format(self.__input_path))