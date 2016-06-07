from ..pixels import get_pixel_id, get_pixel_pos, get_all_pixel_pos


def test_get_pixel_id():
    assert get_pixel_id(0, 0, 0) == 231


def test_get_pixel_pos():
    assert get_pixel_pos(0)[0] == -0.10286
    assert get_pixel_pos(0)[1] == 0.15714


def test_get_all_pixel_pos():
    x, y = get_all_pixel_pos()
    assert x.size == 2048
    assert y.size == 2048
