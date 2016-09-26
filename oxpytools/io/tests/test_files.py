from ..files import CHECInputFile
from ...utils.datasets import get_oxpytools_extra_path
from pathlib import Path


def test_inputfile_simtelarray():
    dataset = get_oxpytools_extra_path("sim_telarray_dummy.gz")
    file = CHECInputFile(dataset)
    datasets_path = Path(get_oxpytools_extra_path(""))
    assert file.input_path == datasets_path.joinpath("sim_telarray_dummy.gz").as_posix()
    assert file.directory == datasets_path.as_posix()
    assert file.extension == ".gz"
    assert file.filename == "sim_telarray_dummy"
    assert file.origin == "hessio"


def test_inputfile_target():
    dataset = get_oxpytools_extra_path("target_sky.fits")
    file = CHECInputFile(dataset)
    datasets_path = Path(get_oxpytools_extra_path(""))
    assert file.input_path == datasets_path.joinpath("target_sky.fits").as_posix()
    assert file.directory == datasets_path.as_posix()
    assert file.extension == ".fits"
    assert file.filename == "target_sky"
    assert file.origin == "targetio"
    source = file.read()
    event = next(source)
    assert event.dl0.tels_with_data == {0}