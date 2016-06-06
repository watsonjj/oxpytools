from ..files import InputFile
from ...utils.datasets import get_oxpytools_extra_path
from pathlib import Path


def test_inputfile_simtelarray():
    dataset = get_oxpytools_extra_path("sim_telarray_trueq.gz")
    file = InputFile(dataset)
    datasets_path = Path(get_oxpytools_extra_path(""))
    assert file.input_path == datasets_path.joinpath("sim_telarray_trueq.gz").as_posix()
    assert file.directory == datasets_path.as_posix()
    assert file.extension == ".gz"
    assert file.filename == "sim_telarray_trueq"
    assert file.type == "simtel"


def test_inputfile_target():
    dataset = get_oxpytools_extra_path("target_sky.fits")
    file = InputFile(dataset)
    datasets_path = Path(get_oxpytools_extra_path(""))
    assert file.input_path == datasets_path.joinpath("target_sky.fits").as_posix()
    assert file.directory == datasets_path.as_posix()
    assert file.extension == ".fits"
    assert file.filename == "target_sky"
    assert file.type == "target"
