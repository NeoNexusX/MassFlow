import math

from massflow.agent_tools.spatial_tools import rasterize_coordinate_values, summarize_raster


def test_rasterize_coordinate_values_builds_expected_image() -> None:
    image, meta = rasterize_coordinate_values(
        [
            (3, 5, 10.0),
            (4, 5, 20.0),
            (3, 6, 30.0),
        ]
    )

    assert image.shape == (2, 2)
    assert image[0, 0] == 10.0
    assert image[0, 1] == 20.0
    assert image[1, 0] == 30.0
    assert math.isnan(image[1, 1])
    assert meta["width"] == 2
    assert meta["height"] == 2
    assert meta["pixel_count"] == 3


def test_summarize_raster_ignores_nan_values() -> None:
    image, _ = rasterize_coordinate_values(
        [
            (1, 1, 5.0),
            (2, 1, 15.0),
        ]
    )
    summary = summarize_raster(image)

    assert summary["count"] == 2
    assert summary["min"] == 5.0
    assert summary["max"] == 15.0
    assert summary["mean"] == 10.0
