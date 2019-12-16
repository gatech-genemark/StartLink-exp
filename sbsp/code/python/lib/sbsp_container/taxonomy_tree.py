import copy
import logging
import pandas as pd
from typing import *

from sbsp_general.general import get_value
from sbsp_io.objects import save_obj, load_obj

logger = logging.getLogger(__name__)


class Node:
    """
    Represents a node in the taxonomy tree.
    """

    def __init__(self, tax_id=None, parent=None, attributes=None):
        # type: (int, Union[Node, None], Dict[str, Any]) -> None

        self.tax_id = tax_id
        if self.tax_id is None:
            raise ValueError("Node ID cannot be None")

        self._parent = parent

        self.attributes = copy.copy(attributes)
        if self.attributes is None:
            self.attributes = dict()

        self._children = dict()
        self.attributes["taxid"] = tax_id

    def add_child(self, child):
        # type: (Node) -> None
        self._children[child.tax_id] = child

    def add_children(self, children):
        # type: (Iterable[Node]) -> None
        for child in children:
            self.add_child(child)

    def children(self):
        return self._children.values()

    def set_parent(self, parent):
        # type: (Node) -> None
        self._parent = parent

    def is_leaf(self):
        # type: () -> bool
        return len(self._children) == 0


class TaxonomyTree:

    def __init__(self, root_node):
        # type: (Node) -> None

        self.root = root_node
        self._check_for_cycles()

    def _check_for_cycles(self):
        # type: () -> None
        pass

    def add_attributes(self, dict_id_to_attributes):
        # type: (Dict[int, Dict[str, Any]]) -> None

        def add_attributes_helper(curr_node):
            # type: (Node) -> None

            # add current info
            if curr_node.tax_id in dict_id_to_attributes:
                curr_node.attributes.update(dict_id_to_attributes[curr_node.tax_id])

            # add for each child
            for child in curr_node.children():
                add_attributes_helper(child)

        add_attributes_helper(self.root)


    @staticmethod
    def _get_names_per_taxid(df):
        # type: (pd.DataFrame) -> Dict[int, str]

        taxid_to_name = dict()
        for index, row in df.iterrows():
            node_id = row.iloc[0]
            name_class = row.iloc[3].strip()

            if name_class == "scientific name":

                taxid_to_name[node_id] = row.iloc[1].strip()

        return taxid_to_name

    @classmethod
    def init_from_file(cls, pf_input, pf_names, file_format=None):
        # type: (str, str, str) -> TaxonomyTree

        """
        Create a taxonomy tree from nodes dump file
        :param pf_input:
        :param file_format:
        :return:
        """

        df_nodes = pd.read_csv(pf_input, header=None, delimiter="|")
        df_names = pd.read_csv(pf_names, header=None, delimiter="|")

        # create node for each value
        dict_tax_id_node = TaxonomyTree._create_node_per_dataframe_row(df_nodes)
        dict_taxid_names = TaxonomyTree._get_names_per_taxid(df_names)

        root_nodes = list()

        # add each node to its parent's children
        for tax_id, node in dict_tax_id_node.items():

            parent_id = node.attributes["parent_id"]

            if node.tax_id in dict_taxid_names:
                node.attributes["name_txt"] = dict_taxid_names[node.tax_id]

            if parent_id is not None:

                parent_node = dict_tax_id_node[parent_id]

                parent_node.add_child(node)
                node._parent = parent_node
            else:
                root_nodes.append(node)

        if len(root_nodes) > 1:
            raise ValueError("More than one root node available")
        if len(root_nodes) == 0:
            raise ValueError("No root node detected")

        return TaxonomyTree(root_nodes[0])

    @staticmethod
    def _create_node_per_dataframe_row(df_nodes):
        # type: (pd.DataFrame) -> Dict[int, Node]

        dict_tax_id_node = dict()

        for index, row in df_nodes.iterrows():
            node_id = row.iloc[0]
            parent_id = row.iloc[1]
            rank = row.iloc[2]
            genetic_code = row.iloc[6]

            if parent_id == node_id:
                parent_id = None

            attributes = {
                "parent_id": parent_id,
                "rank": rank,
                "genetic_code": genetic_code
            }

            dict_tax_id_node[node_id] = Node(
                node_id, attributes=attributes
            )

        return dict_tax_id_node

    def save(self, pf_save):
        # type: (str) -> None
        save_obj(self, pf_save)

    @staticmethod
    def load(pf_load):
        # type: (str) -> TaxonomyTree
        return load_obj(pf_load)

    def to_string(self, **kwargs):
        # type: (Dict[str, Any]) -> str
        return TaxonomyTree.to_string_helper(self.root, depth=0, **kwargs)

    @staticmethod
    def to_string_current_level(node, depth, **kwargs):
        # type: (Node, int, Dict[str, Any]) -> str

        tag_name = get_value(kwargs, "tag_name", None)

        attribute_name = get_value(kwargs, "attribute_name", None)

        output = ""

        single_level = "    |"
        depth_level = single_level * depth

        if depth > 0:
            output = depth_level + "__ "

        # get tag
        tag_value = node.tax_id
        if tag_name is not None:
            tag_value = get_value(node.attributes, tag_name, node.tax_id, default_if_none=True)

        output += str(tag_value)

        if attribute_name is not None:
            output += "\t({})".format(node.attributes[attribute_name])

        return output

    @staticmethod
    def to_string_helper(node, depth, **kwargs):
        # type: (Node, int, Dict[str, Any]) -> str

        max_depth = get_value(kwargs, "max_depth", None)

        # print current node level
        output = TaxonomyTree.to_string_current_level(node, depth, **kwargs) + "\n"

        # print for children if not reached max depth
        if max_depth is None or depth < max_depth:
            for child in node.children():
                output += TaxonomyTree.to_string_helper(
                    child, depth + 1, **kwargs
                )

        return output

    def to_string_tree_with_stats(self, attribute_name, accumulator, func_kwargs, **kwargs):
        # type: (str, Callable[List[Dict[str, Any]], Dict[str, Any]], Dict[str, Any], Dict[str, Any]) -> str
        TaxonomyTree.update_tree_with_stats(self.root, attribute_name, accumulator, func_kwargs, **kwargs)

        return TaxonomyTree.to_string_helper(self.root, 0, attribute_name=attribute_name, **kwargs)

    @staticmethod
    def update_tree_with_stats(node, attribute_name, accumulator, func_kwargs, **kwargs):
        # type: (Node, str, Callable[List[Dict[str, Any]], Dict[str, Any]],Dict[str, Any], Dict[str, Any]) -> None

        # clear stats for current node
        # node.attributes["stats"] = dict()

        # update stats for children
        for child in node.children():
            TaxonomyTree.update_tree_with_stats(child, attribute_name, accumulator, func_kwargs, **kwargs)

        # update stats for current node
        val = accumulator([c.attributes for c in node.children()], node.attributes, attribute_name, **func_kwargs)
        node.attributes[attribute_name] = val

    @staticmethod
    def is_node_with_tag(node, tag, tag_type=None):
        # type: (Node, Any, Union[None, str]) -> None

        if tag_type is None:
            return tag == node.tax_id

        return node.attributes[tag_type] == tag

    def get_node_with_tag(self, ancestor_tag, tag_type):

        current_node = self.root

        lifo = list()

        lifo.append(current_node)

        while len(lifo) > 0:
            p = lifo.pop()

            # if is node we're searching for
            if TaxonomyTree.is_node_with_tag(p, ancestor_tag, tag_type):
                return p

            # otherwise add all children
            for child in p.children():
                lifo.append(child)

        return None

    @staticmethod
    def get_leaves_under_node(node):
        # type: (Node) -> Generator[Node]

        lifo = list()

        lifo.append(node)

        while len(lifo) > 0:
            p = lifo.pop()

            # if is leaf
            if p.is_leaf():
                yield p

            # otherwise add all children
            for child in p.children():
                lifo.append(child)

    def get_genomes_under_ancestor(self, ancestor_tag, tag_type):
        # type: (Union[str, int], str) -> Generator[Dict[str, Any]]

        ancestor_node = self.get_node_with_tag(ancestor_tag, tag_type)

        if ancestor_tag is None:  # empty generator
            return (_ for _ in ())

        for curr_node in TaxonomyTree.get_leaves_under_node(ancestor_node):
            yield curr_node.attributes

    @staticmethod
    def get_nodes_under_ancestor(node):
        # type: (Node) -> Generator[Node]

        lifo = list()

        lifo.append(node)

        while len(lifo) > 0:
            p = lifo.pop()

            yield p

            # otherwise add all children
            for child in p.children():
                lifo.append(child)

    def get_possible_genomes_under_ancestor(self, ancestor_tag, tag_type):
        # type: (Union[str, int], str) -> Generator[Dict[str, Any]]

        ancestor_node = self.get_node_with_tag(ancestor_tag, tag_type)

        if ancestor_tag is None:  # empty generator
            return (_ for _ in ())

        for curr_node in TaxonomyTree.get_nodes_under_ancestor(ancestor_node):
            yield curr_node.attributes
