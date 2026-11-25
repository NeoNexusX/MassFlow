"""MassFlow pipeline example: imzML → H5 ion image."""

from pathlib import Path

from massflow.data_manager import MSDataManagerImzML
from massflow.preprocess import Preprocessor


def main():
    input_path = Path("data/example.imzML")
    output_path = Path("data/example.h5")

    with MSDataManagerImzML(filepath=str(input_path)) as dm:
        dm.load_head_data()

        processed_dm = (
            Preprocessor(
                data_manager=dm,
                output_format="ion_image_h5",
                output_path=str(output_path),
            )
            .start()
        )

        print(f"Output: {output_path}")
        print(f"  n_ions: {processed_dm.n_ions}")
        print(f"  n_pixels: {processed_dm.n_pixels}")
        print(f"  image: {processed_dm.height} x {processed_dm.width}")
        print(f"  size: {output_path.stat().st_size / (1024**3):.2f} GB")
        processed_dm.close()


if __name__ == "__main__":
    main()
