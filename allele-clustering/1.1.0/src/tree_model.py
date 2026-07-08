from dataclasses import dataclass
from typing import Iterator

@dataclass
class Node:
    id: int
    name: str
    branch_length: float
    children: list[int] # list of node IDs
    parent: None|int # node ID

    @classmethod
    def new(
        cls,
        id: int,
        name: str="",
        branch_length: float=0.,
        children: list[int]=None,
        parent: int=None
    ):
        if children is None:
            children = []
        return cls(id, name, branch_length, children, parent)

    @property
    def is_leaf(self) -> bool:
        return len(self.children) == 0

    @property
    def is_internal(self) -> bool:
        return not self.is_leaf

class Tree:
    def __init__(self):
        self.nodes: list[Node] = [Node.new(0)]
        self.root: int = 0

    @classmethod
    def empty(cls):
        tree = cls()
        tree.nodes = []
        tree.root = None
        return tree

    @property
    def leaf_names(self) -> list[str]:
        def traverse(tree: "Tree", node_id: int) -> Iterator[str]:
            this_node = tree.nodes[node_id]
            if this_node.is_leaf:
                yield this_node.name
            else:
                for child in this_node.children:
                    yield from traverse(tree, child)

        return [name for name in traverse(self, self.root)]

    def add_node(
        self,
        parent: int=None,
        name: str="",
        branch_length: float=0.,
    ) -> int:
        """Add a new node to the tree and return its index"""
        if (
            parent is not None
            and (
                parent >= len(self.nodes)
                or parent < 0
                or not isinstance(parent, int)
            )
        ):
            raise ValueError(f"Invalid parent ID provided when adding new node: {parent}")
        new_id = len(self.nodes)
        self.nodes.append(
            Node.new(
                id=new_id,
                name=name,
                branch_length=branch_length,
                children=None,
                parent=parent
            )
        )
        if parent is not None:
            self.nodes[parent].children.append(new_id)
        return new_id

    def add_node_with_children(
        self,
        children: list[int],
        parent: int=None,
        name: str="",
        branch_length: float=0.
    ) -> int:
        new_id = self.add_node(parent, name, branch_length)
        new_node = self.nodes[new_id]
        new_node.children = children
        # add new_node as parent of each child
        for child in children:
            self.nodes[child].parent = new_id
        return new_id

    def add_leaf_polytomy(
        self,
        leaf_names: list[str],
        leaf_brlens: list[float]=None,
        parent: int=None,
        internal_name: str="",
        internal_brlen: float=0.
    ) -> int:
        """Create a polytomy of provided leaves beneath a new internal node

        Assumes provided leaf names and branch lengths (optional) are in the same order.
        If no branch lengths are provided, then branches are set a t 0.0
        """
        if leaf_brlens is None:
            leaf_brlens = [0. for _ in leaf_names]
        internal_id = self.add_node(
            parent=parent,
            name=internal_name,
            branch_length=internal_brlen
        )
        for name, brlen in zip(leaf_names, leaf_brlens):
            self.add_node(
                parent=internal_id,
                name=name,
                branch_length=brlen
            )
        return internal_id

    @classmethod
    def from_newick(cls, nwk_string: str):
        def parse_newick(tree: "Tree", nwk_string: str, idx: int, parent_node: int) -> int:
            """parse substring of newick tree

            From a specified index, parse newick string until a corresponding closing paren is found

            Args:
                tree: Tree instance for accessing and adding nodes
                nwk_string: the newick string being parsed
                idx: where in the newick string should parsing begin for this subtree
                parent_node: the index of the node to which this subtree should be a child

            returns:
                the index of the closing paren terminating this subtree
            """
            node_chars = []
            brlen_digits = []
            in_node_name = True # track what current chars correspond to
            this_node_id = tree.add_node(parent_node)
            while True:
                idx += 1
                char = nwk_string[idx]
                match char:
                    case "(":
                        # subtree found. call recursively
                        idx = parse_newick(tree, nwk_string, idx, this_node_id)
                    case ")":
                        # end of this subtree. Add current node to parent and end
                        try:
                            brlen = float("".join(brlen_digits))
                        except ValueError:
                            raise ValueError(f"Could not parse branch length: {''.join(brlen_digits)}")
                        tree.nodes[this_node_id].branch_length = brlen
                        return idx
                    case ",":
                        # End of this node, next char is start of sibling node
                        try:
                            brlen = float("".join(brlen_digits))
                        except ValueError:
                            raise ValueError(f"Could not parse branch length: {''.join(brlen_digits)}")
                        tree.nodes[this_node_id].branch_length = brlen
                        # reset vars for next node
                        this_node_id = tree.add_node(parent_node)
                        node_chars = []
                        brlen_digits = []
                        in_node_name = True
                    case ":":
                        # End of node name. Next is brlen
                        in_node_name = False
                        tree.nodes[this_node_id].name = "".join(node_chars)
                        node_chars = []
                    case _:
                        if in_node_name:
                            node_chars.append(char)
                        else:
                            brlen_digits.append(char)

        tree = cls()
        idx = 0
        char = nwk_string[idx]
        if char != "(":
            raise ValueError("newick string format issue: tree should start with (")
        parse_newick(tree, nwk_string, idx, tree.root)
        return tree


    def to_newick(self, decimals: int=None) -> str:
        # Set global decimal precision for brlens in recursive function
        if decimals is not None:
            precision = f".{decimals}f"
        else:
            precision = ""

        def subtree_string(tree: "Tree", node_id: int) -> str:
            this_node = tree.nodes[node_id]
            if this_node.is_leaf:
                return f"{this_node.name}:{this_node.branch_length:{precision}}"

            # otherwise traverse subtree nodes to build the newick string
            subtrees = []
            for child in this_node.children:
                subtrees.append(subtree_string(tree, child))
            return f"({','.join(subtrees)}):{this_node.branch_length:{precision}}"

        root_id = self.root
        root =  self.nodes[root_id]
        if root.children == []:
            raise ValueError("Cannot produce newick for a tree whose root has no children.")

        subtrees = []
        for child in root.children:
            subtrees.append(subtree_string(self, child))
        return f"({','.join(subtrees)});"


    def to_grapetree(self, decimals: int=None):
        """Produce a newick representation of a minimum spanning tree that is compatible with grapetree

        Similar to `Tree.to_newick`, except all nodes in the tree are expected to be named leaves and
        each leaf is replaced by a new internal node from which the leaf is a 0-brlen descendent
        """
        # Set global decimal precision for brlens in recursive function
        if decimals is not None:
            precision = f".{decimals}f"
        else:
            precision = ""

        def subtree_string(tree: "Tree", node_id: int) -> str:
            this_node = tree.nodes[node_id]
            if this_node.is_leaf:
                return f"{this_node.name}:{this_node.branch_length:{precision}}"

            # otherwise traverse subtree nodes to build the newick string
            subtrees = [f"{this_node.name}:{0.0:{precision}}"]
            for child in this_node.children:
                subtrees.append(subtree_string(tree, child))
            return f"({','.join(subtrees)}):{this_node.branch_length:{precision}}"

        root_id = self.root
        root =  self.nodes[root_id]
        if root.children == []:
            raise ValueError("Cannot produce newick for a tree whose root has not children.")

        subtrees = [f"{root.name}:{0.0:{precision}}"]
        for child in root.children:
            subtrees.append(subtree_string(self, child))

        return f"({','.join(subtrees)});"


    def midpoint_root(self):
        def subtree_dist_scan(tree: "Tree", node_idx: int, max_dists: list[float]) -> float:
            """calculate maximum internode distance and distance to furthest leaf below this node

            Store maximum internode distance in shared `max_dists` variable.
            Return distance to furthest leaf to parent call
            """
            node = tree.nodes[node_idx]
            if node.is_leaf:
                # don't change max_dist for this node
                # send branch length up to parent call
                return node.branch_length

            # get distance to furthest leaf below each child
            child_dists = []
            for child_id in node.children:
                child_dists.append((child_id, subtree_dist_scan(tree, child_id, max_dists)))
            child_dists.sort(key=lambda x: x[1], reverse=True)

            # handle internal node with a single child
            if len(child_dists) == 1:
                child_id, dist = child_dists[0]
                max_dists[node_idx] = (node_idx, 0, (child_id,), (dist,))
                return dist + node.branch_length

            # furthest internode distance is the sum of the two largest returned vals
            (child_a, dist_a), (child_b, dist_b) = child_dists[:2]
            internode_dist = dist_a + dist_b
            # Update distances in shared list
            max_dists[node_idx] = (node_idx, internode_dist, (child_a, child_b), (dist_a, dist_b))
            # return furthest distance to a leaf below this node
            # plus this node's branch length
            # That's the distance from the parent to that furthest node
            return dist_a + node.branch_length

        def find_midpoint(tree: "Tree") -> tuple[int, float]:
            """Find midpoint in tree

            Returns index of node below midpoint and the branch distance above it to midpoint
            """
            # find largest tree dist between nodes subtended by this node
            max_dists = [(0, 0, (-1, -1), (0, 0)) for _ in self.nodes]
            subtree_dist_scan(self, self.root, max_dists)
            # find node with largest internode distance
            max_node_id, _, child_idxs, child_dists = max(max_dists, key=lambda x: x[1])
            # traverse to branch containing midpoint
            # and determine brlen either side of midpoint
            traverse_target_len = (child_dists[0] - child_dists[1]) / 2
            this_node_id = child_idxs[0]
            traversed = 0
            while True:
                this_node = tree.nodes[this_node_id]
                if (traversed + this_node.branch_length) >= traverse_target_len:
                    # midpoint is in the branch above this node
                    break
                traversed += this_node.branch_length
                _, _, child_idxs, child_dists = max_dists[this_node_id]
                (this_node_id, _) = max(zip(child_idxs, child_dists), key=lambda x: x[1])

            # return node and how far above it midpoint is located
            remaining_length = traverse_target_len - traversed
            return this_node_id, remaining_length


        midpoint_id, length_above = find_midpoint(self)
        midpoint_node = self.nodes[midpoint_id]
        if length_above == 0:
            self.reroot_on_node(midpoint_node)
        else:
            self.reroot_on_branch(midpoint_node, length_above)

    def reroot_on_node(self, node: Node):
        new_root_id = node.id
        new_brlen_above = 0
        old_parent_id = node.parent
        new_parent_id = None
        while True:
            # rotate this node
            if old_parent_id is not None:
                # don't add old root's null parent as a child
                node.children.append(old_parent_id)
            try:
                node.children.remove(new_parent_id)
            except ValueError:
                # valueerror thrown on first round
                pass
            node.parent = new_parent_id
            new_child_brlen = node.branch_length
            node.branch_length = new_brlen_above
            # move to next node and update variables
            new_brlen_above = new_child_brlen
            new_parent_id = node.id

            if old_parent_id is None:
                # just rotated the old root
                break
            node = self.nodes[old_parent_id]
            old_parent_id = node.parent

        self.root = new_root_id

    def reroot_on_branch(self, node: Node, height_above_node: float):
        # node's branch length becomes its distance from the midpoint
        # other branch out of root will be remaining branch length
        node.branch_length = node.branch_length - height_above_node
        # make the new root
        # We'll insert it with the parent of the other node and then
        # call reroot on node as the process is the same from that point
        new_root_id = len(self.nodes)
        new_root = Node.new(
            id=len(self.nodes),
            children=[node.id],
            branch_length=height_above_node
        )
        self.nodes.append(new_root)
        # reparent node
        old_parent_id = node.parent
        node.parent = new_root_id
        # insert the root into the branch we are splitting.
        new_root.parent = old_parent_id
        old_parent = self.nodes[old_parent_id]
        old_parent.children.remove(node.id)
        old_parent.children.append(new_root_id)

        self.reroot_on_node(new_root)

    def collapse_polytomies(self):
        """collapse chains of internal zero-length branches to polytomies"""
        for node in self.nodes:
            # only consider internal non-root nodes
            if node.is_leaf or node.id == self.root:
                continue
            # Is this node zero-length?
            if node.branch_length != 0:
                continue
            # Otherwise this node should be collapsed
            # Move all child nodes to the parent and unparent this node
            if node.parent is None:
                raise AttributeError("Unparented node found when collapsing polytomies")
            parent_node = self.nodes[node.parent]
            parent_node.children.remove(node.id)
            for child_id in node.children:
                child = self.nodes[child_id]
                child.parent = parent_node.id
                parent_node.children.append(child_id)
            node.children = []
            node.parent = None

    @classmethod
    def from_allele_codes(cls, codes: dict[str, list[str]], max_digits: int=6):

        # store index of nodes in tree.nodes for partial allele codes
        code_idxs = {}
        tree = cls()
        for sample, nums in codes.items():
            parent_idx = tree.root
            partial_code = []
            n = 0 # in case no code
            for n, digit in enumerate(nums):
                partial_code.append(digit)
                digit_id = ".".join(partial_code)
                node_idx = code_idxs.get(digit_id)
                if node_idx is None:
                    # Add allele code to tree as child of prior partial code
                    parent_idx = tree.add_node(
                        parent=parent_idx,
                        name=digit_id,
                        branch_length=1.
                    )
                    code_idxs[digit_id] = parent_idx
                else:
                    parent_idx = node_idx

            # Now add sample to tree as child of its final digit
            brlen = max_digits - (n + 1.)
            tree.add_node(
                parent=parent_idx,
                name=sample,
                branch_length=brlen
            )
        return tree
