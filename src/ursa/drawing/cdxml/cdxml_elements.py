# based on https://chemapps.stolaf.edu/iupac/cdx/sdk/IntroCDXML.htm specification
# (severely outdated, but still working)
from __future__ import annotations

from dataclasses import dataclass

from lxml import etree  # nosec: B410

# only used for trusted data


@dataclass
class CDXMLBoundingBox:
    """
    A bounding box is a rectangle that contains an element.
    It is defined by its minimum and maximum x and y coordinates.

    :param min_x: The minimum x coordinate of the bounding box;
    :param min_y: The minimum y coordinate of the bounding box;
    :param max_x: The maximum x coordinate of the bounding box;
    :param max_y: The maximum y coordinate of the bounding box;
    """

    min_x: float
    min_y: float
    max_x: float
    max_y: float

    @property
    def width(self) -> float:
        """
        The width of the bounding box.
        :return: The width of the bounding box;
        """

        return self.max_x - self.min_x

    @property
    def height(self) -> float:
        """
        The height of the bounding box.
        :return: The height of the bounding box;
        """

        return self.max_y - self.min_y

    @property
    def center_y(self) -> float:
        """
        The y coordinate of the vertical center of the bounding box.
        :return: The y coordinate of the center of the bounding box;
        """

        return self.min_y + self.height / 2


@dataclass
class CDXMLLine:
    """
    Describes a line.

    :param beg_pos_x: The x coordinate of the tail of the line;
    :param beg_pos_y: The y coordinate of the tail of the line;
    :param end_pos_x: The x coordinate of the head of the line;
    :param end_pos_y: The y coordinate of the head of the line;
    """

    beg_pos_x: float
    beg_pos_y: float
    end_pos_x: float
    end_pos_y: float

    def _prop_dict(self) -> dict[str, str]:
        """
        Returns a dictionary of properties for the line element.
        :return: A dictionary of properties for the line element;
        """

        return {
            "BoundingBox": (
                f"{self.beg_pos_x:.2f} {self.beg_pos_y:.2f} "
                f"{self.end_pos_x:.2f} {self.end_pos_y:.2f}"
            ),
            "GraphicType": "Line",
        }

    def as_cdxml_element(self) -> etree.Element:
        """
        Returns the CDXML element for the line.
        :return: The CDXML element for the line;
        """

        return etree.Element("graphic", **self._prop_dict())


@dataclass
class CDXMLRectangle:
    """
    Describes a rectangle.

    :param beg_pos_x: The x coordinate of the top-left corner of the rectangle;
    :param beg_pos_y: The y coordinate of the top-left corner of the rectangle;
    :param end_pos_x: The x coordinate of the bottom-right corner of the rectangle;
    :param end_pos_y: The y coordinate of the bottom-right corner of the rectangle;
    """

    beg_pos_x: float
    beg_pos_y: float
    end_pos_x: float
    end_pos_y: float

    @classmethod
    def from_bbox(cls, bbox: CDXMLBoundingBox) -> CDXMLRectangle:
        """
        Creates a CDXMLRectangle from a bounding box.
        :param bbox: The bounding box;
        :return: The CDXMLRectangle object;
        """

        return cls(
            beg_pos_x=bbox.min_x,
            beg_pos_y=bbox.min_y,
            end_pos_x=bbox.max_x,
            end_pos_y=bbox.max_y,
        )

    def _prop_dict(self) -> dict[str, str]:
        """
        Returns a dictionary of properties for the rectangle element.
        :return: A dictionary of properties for the rectangle element;
        """

        return {
            "BoundingBox": (
                f"{self.beg_pos_x:.2f} {self.beg_pos_y:.2f} "
                f"{self.end_pos_x:.2f} {self.end_pos_y:.2f}"
            ),
            "GraphicType": "Rectangle",
        }

    def as_cdxml_element(self) -> etree.Element:
        """
        Returns the CDXML element for the rectangle.
        :return: The CDXML element for the rectangle;
        """

        return etree.Element("graphic", **self._prop_dict())


@dataclass
class CDXMLText:
    """
    Describes a text element.

    :param text: The text of the element;
    :param pos_x: The x coordinate of the text;
    :param pos_y: The y coordinate of the text;
    """

    text: str
    pos_x: float
    pos_y: float

    def as_cdxml_element(self) -> etree.Element:
        """
        Returns the CDXML element for the text.
        :return: The CDXML element for the text;
        """

        t_element = etree.Element("t", p=f"{self.pos_x:.2f} {self.pos_y:.2f}")
        s_element = etree.SubElement(
            t_element,
            "s",
            font="1",
            size="10.0",
            face="96",
            color="0",
        )
        s_element.text = self.text

        return t_element
