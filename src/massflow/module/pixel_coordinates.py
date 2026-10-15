from typing import Union, Sequence
import numpy as np
from massflow.tools.logger import get_logger

logger = get_logger("massflow.module.pixel_coordinates")


class PixelCoordinates:
    """
    Coordinate holder for MSI pixels.

    This class stores raw `x`, `y`, `z`coordinates without automatic base adjustment.

    Attributes:
        x (int): X coordinate.
        y (int): Y coordinate.
        z (int): Z coordinate.
    """

    def __init__(
        self, x: Union[int, Sequence[int], "PixelCoordinates"], y: int = 0, z: int = 0
    ):
        """
        Initialize pixel coordinates.

        Args:
            x (int): X coordinate value.
            y (int): Y coordinate value.
            z (int): Z coordinate value.

        Raises:
            TypeError: If any of x, y, or z cannot be converted to int.
        """

        # Store raw coordinates internally.
        self._x = None
        self._y = None
        self._z = None

        if isinstance(x, (list, tuple, np.ndarray)) and len(x) == 3:
            self._x, self._y, self._z = x

        elif isinstance(x, PixelCoordinates):
            self._x = x._x
            self._y = x._y
            self._z = x._z

        elif isinstance(x, (int, np.integer)):
            self._x = x
            self._y = y
            self._z = z

        else:
            logger.error(
                "List x must have exactly three elements with x, y, z coordinates."
            )
            raise ValueError(
                "List x must have exactly three elements with x, y, z coordinates."
            )

    @staticmethod
    def _require(value, axis: str):
        if value is None:
            logger.error(f"{axis} coordinate is not set.")
            raise ValueError(f"{axis} coordinate is not set.")
        return value

    @property
    def x(self) -> int:
        """
        Get the X coordinate.

        Returns:
            int: X coordinate.
        """
        return self._require(self._x, "X")

    @x.setter
    def x(self, value: int):
        self._x = value

    @property
    def y(self) -> int:
        """
        Get the Y coordinate.

        Returns:
            int: Y coordinate.
        """
        return self._require(self._y, "Y")

    @y.setter
    def y(self, value: int):
        self._y = value

    @property
    def z(self) -> int:
        """
        Get the Z coordinate.
        """
        return self._require(self._z, "Z")

    @z.setter
    def z(self, value: int):
        self._z = value

    def get_tuple(self) -> tuple:
        """
        Get coordinates as a tuple.

        Returns:
            tuple: (x, y, z) coordinates.

        Raises:
            ValueError: If any coordinate is not set.
        """
        return (self.x, self.y, self.z)

    def __eq__(self, other: object) -> bool:
        """
        Compare equality based on adjusted coordinates.

        Args:
            other (PixelCoordinates): Another PixelCoordinates instance to compare.

        Returns:
            bool: True if all adjusted coordinates match; False otherwise.
        """
        if not isinstance(other, PixelCoordinates):
            return NotImplemented

        return self.x == other.x and self.y == other.y and self.z == other.z

    def __hash__(self):
        """
        Compute hash based on adjusted coordinates.

        Returns:
            int: Hash value derived from adjusted (x, y, z).
        """
        return hash((self.x, self.y, self.z))

    def __repr__(self) -> str:
        """
        Return a readable string representation of the coordinates.

        Returns:
            str: A string formatted as '(x, y, z)'.

        Raises:
            None
        """
        return f"({self.x}, {self.y}, {self.z})"

    def __len__(self) -> int:
        """
        Return the number of dimensions.

        Returns:
            int: Returns 3 when all coordinates are initialized; otherwise 0.

        Raises:
            None
        """
        if self._x is None or self._y is None or self._z is None:
            return 0
        return 3
