from pathlib import Path
from astropy import log

def get_oxpytools_extra_path(filename):
    """Get path to a dataset file inside `oxpytools-extra`.
    """

    import oxpytools
    path = Path(oxpytools.__file__).parent.parent.joinpath('oxpytools-extra')
    dataset = Path(path.joinpath('datasets').joinpath(filename))
    try:
        if not dataset.exists():
            raise FileNotFoundError
    except FileNotFoundError as e:
        log.exception("file path does not exist: '{}'".format(dataset))

    return dataset
