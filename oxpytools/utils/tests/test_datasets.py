from ..datasets import get_oxpytools_extra_path
# from pathlib import Path


def test_get_oxpytools_extra_path():
    path = get_oxpytools_extra_path("")
    assert path is not None
