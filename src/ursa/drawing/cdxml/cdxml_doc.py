import importlib.resources
import itertools

from lxml import etree  # nosec: B410

from .cdxml_tree import CDXMLTree

# only used for trusted data

VERTICAL_SPACING = 100.0


def build_cdxml_document(routes: list[dict]) -> str:
    """
    Builds a CDXML document from a list of route dict representations.
    :param routes: The list of route dicts;
    :return: The CDXML document;
    """

    template_path = importlib.resources.files(__package__).joinpath(
        "acs_1996_template.cdxml"
    )
    with template_path.open() as f:
        # bandit is complaining about the use of etree.parse, but it is safe
        # because we are using a local file integrated into the package
        doc = etree.parse(f)  # nosec: B320
    root = doc.getroot()
    page = etree.SubElement(root, "page")

    obj_id_it = itertools.count(1)
    y_offset = 0.0
    doc_max_x = 0.0
    doc_max_y = 0.0

    for route in routes:
        tree = CDXMLTree.from_route(route, obj_id_it)
        tree.translate(x=0, y=y_offset)
        page.append(tree.as_cdxml_element())

        tree_bbox = tree.get_bounding_box()
        y_offset += tree_bbox.height + VERTICAL_SPACING
        doc_max_x = max(doc_max_x, tree_bbox.max_x)
        doc_max_y = tree_bbox.max_y

    root.set("BoundingBox", f"0 0 {doc_max_x:.2f} {doc_max_y:.2f}")

    return etree.tostring(doc, pretty_print=True).decode("utf-8")
