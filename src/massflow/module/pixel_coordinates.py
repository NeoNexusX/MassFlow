from typing import Optional
from massflow.logger import get_logger

logger = get_logger("ms_module")

class PixelCoordinates:
    """
    Coordinate holder for MSI pixels with optional zero-based view.

    This class stores raw `x`, `y`, `z` coordinates and provides a `zero_based`
    switch that, when enabled, makes the public getters return values offset by `-1`
    (useful when converting 1-based coordinates to 0-based indices).

    Attributes:
        x (int): X coordinate (getter returns adjusted value if `zero_based` is True).
        y (int): Y coordinate (getter returns adjusted value if `zero_based` is True).
        z (int): Z coordinate (getter returns adjusted value if `zero_based` is True).
        zero_based (bool): Whether to expose coordinates in zero-based form.
    """

    def __init__(self, x: int, y: int, z: int, zero_based: bool = False):
        """
        Initialize pixel coordinates with optional zero-based adjustment.

        Args:
            x (int): X coordinate value.
            y (int): Y coordinate value.
            z (int): Z coordinate value.
            zero_based (bool): If True, getters will return (value - 1)
                to convert 1-based coordinates to 0-based.

        Returns:
            None

        Raises:
            TypeError: If any of x, y, or z cannot be converted to int.
        """

        # internal storage
        self._zero_based = None
        # Store raw coordinates internally; getters apply adjustment when required.
        self._x : Optional[int] = None
        self._y : Optional[int] = None
        self._z : Optional[int] = None

        self.x = x
        self.y = y
        self.z = z
        self.zero_based = zero_based

    @property
    def x(self) -> Optional[int]:
        """
        Get the X coordinate.

        Returns:
            int: X coordinate, adjusted by `-1` when `zero_based` is True.
        """
        if self._x is not None:
            return self._x - 1 if self.zero_based else self._x
        else:
            logger.error("X coordinate is not set.")
            return None

    @x.setter
    def x(self, value: int):
        self._x = value

    @property
    def y(self) -> Optional[int]:
        """
        Get the Y coordinate.

        Returns:
            int: Y coordinate, adjusted by `-1` when `zero_based` is True.
        """
        if self._y is  not None:
            return self._y - 1 if self.zero_based else self._y
        else:
            logger.error("Y coordinate is not set.")
            return None

    @y.setter
    def y(self, value: int):
        self._y = value

    @property
    def z(self) -> Optional[int]:
        """
        Get the Z coordinate.

        Returns:
            int: Z coordinate, adjusted by `-1` when `zero_based` is True.
        """
        if self._z is not None:
            return self._z - 1 if self.zero_based else self._z
        else:
            logger.error("Z coordinate is not set.")
            return None

    @z.setter
    def z(self, value: int):
        self._z = value

    @property
    def zero_based(self) -> Optional[bool]:
        """
        Flag indicating whether getters return zero-based coordinates.

        Returns:
            bool: True if coordinates are zero-based; False otherwise.
        """
        if self._zero_based is not None:
            return self._zero_based
        else:
            logger.error("zero_based flag is not set.")
            return None

    @zero_based.setter
    def zero_based(self, value: bool):
        """
        Set the flag controlling zero-based view of coordinates.

        Args:
            value (bool): True for zero-based; False for one-based.
        """
        self._zero_based = bool(value)

    def get_tuple(self) -> Optional[tuple]:
        """
        Get coordinates as a tuple.

        Returns:
            tuple: (x, y, z) coordinates, adjusted by `-1` if `zero_based` is True.
        """
        if self.x is not None and self.y is not None and self.z is not None:
            return (self.x, self.y, self.z)
        else:
            logger.error("One or more coordinates are not set.")
            return None

    def __eq__(self, other: object) -> bool:
        """
        Compare equality based on adjusted coordinates.

        Args:
            other (PixelCoordinates): Another PixelCoordinates instance to compare.

        Returns:
            bool: True if all adjusted coordinates match; False otherwise.

        Raises:
            TypeError: If `other` is not an instance of PixelCoordinates.
        """
        if not isinstance(other, PixelCoordinates):
            raise TypeError("other must be a PixelCoordinates instance")
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __hash__(self):
        """
        Compute hash based on adjusted coordinates.

        Returns:
            int: Hash value derived from adjusted (x, y, z).

        Raises:
            None
        """
        return hash((self.x, self.y, self.z))

    def lefter(self, other: "PixelCoordinates") -> Optional[bool]:
        """
        Determine if this pixel is to the left of another pixel.

        Args:
            other (PixelCoordinates): The other pixel to compare against.

        Returns:
            bool: True if this pixel's X is less than the other's X; otherwise False.
        """
        if self.x is None or other.x is None:
            logger.error("X coordinate is None, cannot compare.")
            return None
        else:
            return self.x < other.x

    def righter(self, other: "PixelCoordinates") -> Optional[bool]:
        """
        Determine if this pixel is to the right of another pixel.

        Args:
            other (PixelCoordinates): The other pixel to compare against.

        Returns:
            bool: True if this pixel's X is greater than the other's X; otherwise False.
        """
        if self.x is None or other.x is None:
            logger.error("X coordinate is None, cannot compare.")
            return None
        else:
            return self.x > other.x

    def upper(self, other: "PixelCoordinates") -> Optional[bool]:
        """
        Determine if this pixel is above another pixel.

        Args:
            other (PixelCoordinates): The other pixel to compare against.

        Returns:
            bool: True if this pixel's Y is greater than the other's Y; otherwise False.
        """
        if self.y is None or other.y is None:
            logger.error("Y coordinate is None, cannot compare.")
            return None
        else:
            return self.y > other.y

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