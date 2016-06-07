"""Methods involving the obtaining of datasets from oxpytools-extra

"""

from pathlib import Path
from astropy import log
import os
import oxpytools


def get_oxpytools_extra_path(filename):
    """
    Get path to a dataset file inside `oxpytools-extra`.

    Parameters
    ----------
    filename : str
        Name of the file inside oxpytools-extra/datasets including extension
        e.g. target_sky.fits

    Returns
    -------
    str
        Path to the dataset file

    """

    environ_variable_name = 'OXPYTOOLS_EXTRA_DIR'

    try:
        path = Path(os.environ[environ_variable_name])
    except KeyError:
        path = Path(oxpytools.__file__).parent.parent.joinpath('oxpytools-extra')

    dataset = Path(path.joinpath('datasets').joinpath(filename))
    try:
        if not dataset.exists():
            raise FileNotFoundError
    except FileNotFoundError as e:
        log.exception("file path does not exist: '{}'".format(dataset))

    return dataset.as_posix()
