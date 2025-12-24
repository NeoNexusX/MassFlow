from typing import Optional,Union,Sequence
from massflow.logger import get_logger

logger = get_logger("ms_module")

class PixelCoordinates:
    """
    Coordinate holder for MSI pixels.

    This class stores raw `x`, `y`, `z` coordinates without automatic base adjustment.

    Attributes:
        x (int): X coordinate.
        y (int): Y coordinate.
        z (int): Z coordinate.
    """
    
    def __init__(self, x: Union[int,Sequence[int],"PixelCoordinates"], y: int=0, z: int=0):
        """
        Initialize pixel coordinates.

        Args:
            x (int): X coordinate value.
            y (int): Y coordinate value.
            z (int): Z coordinate value.

        Returns:
            None

        Raises:
            TypeError: If any of x, y, or z cannot be converted to int.
        """

        # internal storage
        # Store raw coordinates internally.
        self._x : Optional[int] = None
        self._y : Optional[int] = None
        self._z : Optional[int] = None

        if isinstance(x, (list, tuple)):
            if len(x) == 3:
                self.x, self.y, self.z = x
            else:
                logger.error("List x must have exactly three elements with x, y, z coordinates.")
                raise ValueError("List x must have exactly three elements with x, y, z coordinates..")
        
        elif isinstance(x, PixelCoordinates):
            self.x = x.x # type: ignore
            self.y = x.y # type: ignore
            self.z = x.z # type: ignore

        else:
            self.x = x # type: ignore
            self.y = y
            self.z = z

    @property
    def x(self) -> Optional[int]:
        """
        Get the X coordinate.

        Returns:
            int: X coordinate.
        """
        if self._x is not None:
            return self._x
        else:
            logger.error("X coordinate is not set.")
            return None

    @x.setter
    def x(self, value: int):
        self._x = int(value)

    @property
    def y(self) -> Optional[int]:
        """
        Get the Y coordinate.

        Returns:
            int: Y coordinate.
        """
        if self._y is  not None:
            return self._y
        else:
            logger.error("Y coordinate is not set.")
            return None

    @y.setter
    def y(self, value: int):
        self._y = int(value)

    @property
    def z(self) -> Optional[int]:
        """
        Get the Z coordinate.

        Returns:
            int: Z coordinate.
        """
        if self._z is not None:
            return self._z
        else:
            logger.error("Z coordinate is not set.")
            return None

    @z.setter
    def z(self, value: int):
        self._z = int(value)

    def get_tuple(self) -> tuple:
        """
        Get coordinates as a tuple.

        Returns:
            tuple: (x, y, z) coordinates.
        """
        if self.x is not None and self.y is not None and self.z is not None:
            return (self.x, self.y, self.z)
        else:
            logger.error("One or more coordinates are not set.")
            raise ValueError("One or more coordinates are not set.")

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