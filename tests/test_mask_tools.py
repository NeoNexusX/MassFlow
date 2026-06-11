import numpy as np
import unittest

from massflow.agent_tools.mask_tools import (
    binarize_mask_image,
    select_coordinates_from_mask,
    summarize_selected_coordinates,
)


class TestMaskTools(unittest.TestCase):
    def test_binarize_mask_image_uses_white_foreground_threshold(self) -> None:
        mask_image = np.array([[0.2, 0.5], [0.8, 1.0]], dtype=np.float32)

        binary_mask = binarize_mask_image(
            mask_image,
            threshold=0.5,
            keep_mode="white_foreground",
        )

        self.assertEqual(binary_mask.tolist(), [[False, True], [True, True]])

    def test_binarize_mask_image_uses_black_foreground_threshold(self) -> None:
        mask_image = np.array([[0.2, 0.5], [0.8, 1.0]], dtype=np.float32)

        binary_mask = binarize_mask_image(
            mask_image,
            threshold=0.5,
            keep_mode="black_foreground",
        )

        self.assertEqual(binary_mask.tolist(), [[True, True], [False, False]])

    def test_binarize_mask_image_rejects_unknown_keep_mode(self) -> None:
        with self.assertRaisesRegex(ValueError, "keep_mode"):
            binarize_mask_image(
                np.ones((2, 2), dtype=np.float32),
                threshold=0.5,
                keep_mode="invalid",
            )

    def test_select_coordinates_from_mask_maps_pixels_back_to_original_coordinates(self) -> None:
        selected_coordinates = select_coordinates_from_mask(
            [
                (10, 20, 1.0),
                (11, 20, 1.0),
                (10, 21, 1.0),
                (11, 21, 1.0),
            ],
            np.array([[True, False], [False, True]], dtype=bool),
            {
                "min_x": 10,
                "min_y": 20,
                "width": 2,
                "height": 2,
            },
        )

        self.assertEqual(selected_coordinates, [[10, 20, 1], [11, 21, 1]])

    def test_summarize_selected_coordinates_returns_expected_bbox(self) -> None:
        summary = summarize_selected_coordinates(
            [
                [10, 20, 1],
                [12, 21, 1],
                [11, 24, 1],
            ]
        )

        self.assertEqual(
            summary,
            {
                "selected_count": 3,
                "min_x": 10,
                "max_x": 12,
                "min_y": 20,
                "max_y": 24,
            },
        )

    def test_summarize_selected_coordinates_handles_empty_selection(self) -> None:
        self.assertEqual(
            summarize_selected_coordinates([]),
            {
                "selected_count": 0,
                "min_x": None,
                "max_x": None,
                "min_y": None,
                "max_y": None,
            },
        )


if __name__ == "__main__":
    unittest.main()
