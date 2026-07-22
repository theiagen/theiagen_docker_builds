import re
from dataclasses import dataclass, field
from math import sqrt
from abc import ABC

import numpy as np
import numpy.typing as npt

from tree_model import Tree


class CondensedMatrix(ABC):
    """Base class with methods for working with a condensed matrix"""
    nrows: int

    @staticmethod
    def condense_matrix(matrix: np.ndarray) -> npt.NDArray[np.float32]:
        """Take a square 2d matrix and convert it into a 1d matrix of the upper triangle"""
        # Init empty array
        nrows = matrix.shape[0]
        array_size = nrows * (nrows -1) // 2
        new_mat = np.empty(array_size, dtype=np.float32)
        # Populate array with top triangle
        new_idx = 0
        for i in range(nrows):
            for j in range(i+1, nrows):
                new_mat[new_idx] = matrix[i, j]
                new_idx += 1
        return new_mat

    def condensed_idx(self, i: int, j: int) -> int:
        """Convert a square matrix index into the equivalent condensed index"""
        if i == j:
            raise ValueError("diagonal indices are not represented in a condensed matrix")
        if i >= self.nrows or j >= self.nrows:
            raise ValueError("indices cannot be greater than dimensions of the square array")
        if i > j:
            # convert i to lower number to work on upper triangle
            i, j = j, i
        idx = (self.nrows * i) - (i * (i + 1) //2) - 1 - i + j
        return idx

    def square_idx(self, idx: int, nrows: int=None) -> tuple[int, int]:
        if nrows is None:
            nrows = self.nrows
        b = 1 - (2 * nrows)
        i = (-b - sqrt(b ** 2 - 8 * idx)) // 2
        j = idx + i * (b + i + 2) // 2 + 1
        return int(i), int(j)

    def condensed_offset(self, i: int) -> int:
        """get the offset of the beginning of a row"""
        off = (self.nrows * i) - (i * (i + 1) //2) - 1 - i
        return off

    def as_square(self, matrix: np.ndarray) -> npt.NDArray[np.float32]:
        """create squareform of provided condensed matrix"""
        nrows = int((1 + sqrt(1 + 8*matrix.shape[0])) / 2)
        square = np.zeros((nrows, nrows), dtype=np.float32)
        for cond_idx, val in enumerate(matrix):
            i, j = self.square_idx(cond_idx)
            if val is np.ma.masked:
                val = -1
            square[i,j] = square[j,i] = val
        return square


@dataclass
class DistanceMatrix(CondensedMatrix):
    # private attributes handled by class methods or properties
    _implemented_distance_metrics = None
    _distance_metric: str = None
    _implemented_algorithms = None
    _algorithm: callable = None
    _samples: np.array = None
    _redundants: dict[str, list[str]] = None # Mapping representatives to other samples with identical profiles
    _profiles: np.array = None
    _distmat: np.array = None # distmat with groups of identical samples collapsed to 1 representative
    _tree: str = None # newick format tree
    _leaf_order: list[str] = None
    _algorithm_settings: dict = field(default_factory=dict)

    def __post_init__(self):
        self._implemented_algorithms = {
            "upgma": self.upgma,
            "single_linkage": self.single_linkage,
            "complete_linkage": self.complete_linkage,
            "neighbor_joining": self.neighbor_joining,
            "minimum_spanning": self.minimum_spanning
        }
        self._implemented_distance_metrics = {
            "normalized_allele_differences": self.normalized_allele_differences,
            "absolute_allele_differences": self.absolute_allele_differences
        }

    @property
    def distance_metric(self) -> callable:
        if self._distance_metric is None:
            raise Exception("distance metric has not been set")
        return self._distance_metric


    @distance_metric.setter
    def distance_metric(self, metric: str):
        if metric not in self._implemented_distance_metrics:
            raise ValueError(
                f"Distance metric {metric} not implemented. "
                f"Distance metric must be one of {', '.join(self._implemented_distance_metrics)}"
                )
        else:
            self._distance_metric = self._implemented_distance_metrics[metric]


    @property
    def algorithm(self) -> callable:
        if self._algorithm is None:
            raise Exception("algorithm has not been set")
        return self._algorithm


    @algorithm.setter
    def algorithm(self, algorithm: str):
        if algorithm not in self._implemented_algorithms:
            raise ValueError(
                f"Algorithm {algorithm} not implemented. "
                f"Algorithm must be one of {', '.join(self._implemented_algorithms)}"
                )
        else:
            self._algorithm = self._implemented_algorithms[algorithm]


    @property
    def samples(self) -> list[str]:
        return self._samples


    @property
    def redundants(self) -> dict[str, list[str]]:
        return self._redundants


    @property
    def profiles(self) -> npt.NDArray[np.intp]:
        """Check profiles exist before returning
        """
        if self._profiles is not None:
            return self._profiles
        else:
            raise AttributeError("No profiles loaded")

    @property
    def distmat(self) -> npt.NDArray[np.float32]:
        """Return a copy of the distmat to ensure the original is unchanged
        This is the distance matrix used for tree inference. Removing redundant samples is faster
        """
        if self._distmat is not None:
            return self._distmat.copy()
        else:
            raise AttributeError("Distance matrix has not been generated")


    @property
    def tree(self) -> str:
        """return newick string if a tree has been inferred else infer one"""
        if self._tree is not None:
            return self._tree
        else:
            self.infer_tree()
            return self._tree

    @property
    def leaf_order(self) -> list[str]:
        if self._leaf_order is None:
            return []
        return self._leaf_order


    def read_profile_dict(self, profile_dict: dict[str, list[str]]):
        """Read sample names and their allele profiles into np.arrays
        profile_dict:
            dict of profile. {"Headers": [headers], "Samples_N": [profile, for, sample, N]}
        """
        names, profiles = [], []
        for k,v in profile_dict.items():
            if k == "Headers":
                continue
            names.append(k)
            profiles.append(np.array([int(i) if i not in ("", "N", "-", "0") else 0 for i in v], dtype=np.int64))
        self._profiles = np.array(profiles)
        self._samples = np.array([re.sub(r'[\(\)\ \,\"\';]', '_', n) for n in names])

    def read_alignment_dict(self, alignment_dict: dict[str, str]):
        """Read sample names and their allele profiles into np.arrays
        alignment_dict:
            dict of aligned seqs. {"Samples_N": "AT-NCG"...}
        """
        names, profiles = [], []
        for k,v in alignment_dict.items():
            if k == "Headers":
                continue
            names.append(k)
            profiles.append(np.array([i for i in v]))

        self._profiles = np.array(profiles)
        self._samples = np.array([re.sub(r'[\(\)\ \,\"\';]', '_', n) for n in names])

    def nonredundant(self):
        """Remove redundant and uninformative information to speed up processing
        Remove loci from profile if no samples differ, collapse identical samples to a single representative
        """
        samples = self.samples
        profiles = self.profiles
        encoded_profiles = np.array([np.unique(p, return_inverse=True)[1]+1 for p in profiles.T]).T
        encoded_profiles[(profiles == 0)] = 0

        samples = samples[np.lexsort(encoded_profiles.T)]
        profiles = encoded_profiles[np.lexsort(encoded_profiles.T)]
        presence = (np.sum(profiles > 0, 1) > 0)
        samples, profiles = samples[presence], profiles[presence]

        uniqueness = np.concatenate([[1], np.sum(np.diff(profiles, axis=0) != 0, 1) > 0])

        self._redundants = {samples[0]:[]}
        redundant_group = self.redundants[samples[0]]
        for n, u in zip(samples, uniqueness) :
            if u == 0 :
                redundant_group.append(n)
            else :
                self._redundants[n] = [n]
                redundant_group = self.redundants[n]
        self._samples = samples[uniqueness>0]
        self._profiles = profiles[uniqueness>0]


    def normalized_allele_differences(self):
        profiles = self.profiles
        self.nrows = profiles.shape[0]
        presences = (profiles > 0)
        array_size = self.nrows * (self.nrows -1) // 2
        distances = np.zeros(shape=array_size, dtype=np.float32)
        offset = 0
        for n, i in enumerate(np.arange(0, self.nrows-1)):
            profile, presence = profiles[i], presences[i]
            diffs = 100*(
                np.sum(
                    # from n onward as this is only half the triangle
                    ((profiles[n+1:,:] != profile) & (presences[n+1:,:] * presence)), axis=1)
                /np.sum(presences[n+1:,:] & presence, axis=1)
            )
            this_row_len = self.nrows - (n + 1)
            idxs = np.arange(offset, offset + this_row_len, dtype=np.intp)
            offset += this_row_len
            distances[idxs] = diffs
        self._distmat = distances

    def absolute_allele_differences(self):
        profiles = self.profiles
        self.nrows = profiles.shape[0]
        presences = (profiles > 0)
        array_size = self.nrows * (self.nrows -1) // 2
        distances = np.zeros(shape=array_size, dtype=np.float32)
        offset = 0
        for n, i in enumerate(np.arange(0, self.nrows-1)):
            profile, presence = profiles[i], presences[i]
            diffs = np.sum(
                # from n onward as this is only half the triangle
                ((profiles[n+1:,:] != profile) & (presences[n+1:,:] * presence)),
                axis=1,
            )
            this_row_len = self.nrows - (n + 1)
            idxs = np.arange(offset, offset + this_row_len, dtype=np.intp)
            offset += this_row_len
            distances[idxs] = diffs
        self._distmat = distances

    def symmetric_distance(self):
        self.distance_metric()
        del self._profiles


    def add_redundant_samples_to_distmat(self, distmat: np.ndarray):
        samples = self.samples
        for i, name in enumerate(samples[::-1]):
            idx = len(samples) - (i+1)
            if len(self.redundants[name]) > 1:
                num_reps = len(self.redundants[name])-1
                distmat = np.insert(distmat, [idx]*num_reps, distmat[idx], axis=0)
                distmat = np.insert(distmat, [idx]*num_reps, distmat[:,[idx]], axis=1)

        return distmat


    def reorder_distmat(self, distmat: np.ndarray) -> npt.NDArray[np.float32]:
        """Reorder distmat to match leaves and add in redundant samples
        """
        redundant_names = []
        for n in self.samples:
            redundant_names += self.redundants[n]

        loc_lookup = {i:n for n,i in enumerate(redundant_names)}
        new_order = [loc_lookup[n] for n in self.leaf_order]
        return distmat[np.ix_(new_order, new_order)]


    def upgma(self, **_):
        self.cluster("upgma", **self._algorithm_settings)

    def complete_linkage(self, **_ ):
        self.cluster("complete_linkage", **self._algorithm_settings)

    def single_linkage(self, **_):
        self.cluster("single_linkage", **self._algorithm_settings)

    def cluster( #del
        self,
        linkage_method: str,
        max_tree_height: int=0,
        collapse_zero_length_nodes: bool=True,
        brlen_rounding: int=None,
        **_
    ):
        # need a node lookup rather than directly indexing tree nodes
        # because identical samples will be added as a polytomy descended from one sample
        tree = Tree.empty()
        node_lookup = []
        for sample in self.samples:
            leaves = self.redundants[sample]
            # convert to str for leaf naming
            leaves = [str(name) for name in leaves]
            if len(leaves) == 1:
                this_id = tree.add_node(None, leaves[0], 0.)
            else:
                this_id = tree.add_leaf_polytomy(
                    leaf_names=leaves,
                )
            node_lookup.append(this_id)

        arr = self.distmat
        node_heights = {}
        # set initial number of subtended leaves of each node in the lookup
        clust_sizes = {node_id: len(reds) for node_id, reds in zip(node_lookup, self.redundants.values())}
        # track number of times to repeat until end.
        remaining = self.nrows - 1
        while remaining > 0:
            minloc = arr.argmin()
            min_idx, max_idx = self.square_idx(minloc, remaining+1)
            edge_len = arr[minloc]
            min_node_id = node_lookup[min_idx]
            max_node_id = node_lookup[max_idx]

            # If we are capping tree height, adjust the node height
            if max_tree_height:
                edge_len = min(edge_len, max_tree_height)

            # Adjust edge len to account for incorporated dists
            original_edge_len = edge_len
            min_edge = (edge_len/2) - node_heights.get(min_node_id, 0)
            max_edge = (edge_len/2) - node_heights.get(max_node_id, 0)

            # add to tree
            tree.nodes[min_node_id].branch_length = min_edge
            tree.nodes[max_node_id].branch_length = max_edge

            new_node = tree.add_node_with_children([min_node_id, max_node_id])
            node_lookup[max_idx] = new_node
            del node_lookup[min_idx]
            node_heights[new_node] = original_edge_len/2

            # update matrix with new values at max index row, column based on linkage method used
            # add the first value in the matrix just to keep the correspondence between the numbers
            # in the min and max dists lists.
            # The value won't be used - it's just for index correspondence
            min_idxs = self.get_condensed_idxs_of_sample(min_idx, remaining+1, to_skip=max_idx)
            min_dists = arr[min_idxs]
            max_idxs = self.get_condensed_idxs_of_sample(max_idx, remaining+1, to_skip=min_idx)
            max_dists = arr[max_idxs]

            if linkage_method == "upgma":
                # get cluster size of min and max to weight distances
                min_idx_weight = clust_sizes[min_node_id]
                max_idx_weight = clust_sizes[max_node_id]
                new_dists = (
                    (
                        (min_dists*min_idx_weight)
                        +(max_dists*max_idx_weight)
                    )
                    /(min_idx_weight+max_idx_weight)
                )
                # add sizes of joined clusters
                clust_sizes[new_node] = min_idx_weight + max_idx_weight
            elif linkage_method == "single_linkage":
                new_dists = np.ma.min([min_dists, max_dists], axis=0)
            elif linkage_method == "complete_linkage":
                new_dists = np.ma.max([min_dists, max_dists], axis=0)
            else:
                raise ValueError(f"Algorithm {self._algorithm_settings['algorithm']} not implemented as a method, "
                    "but is in the list of implemented algorithms. This is a bug.")
                break
            # Set new distances for max node
            arr[max_idxs] = new_dists
            to_del = self.get_condensed_idxs_of_sample(min_idx, remaining+1)
            arr = np.delete(arr, to_del)
            remaining -= 1

        # root on last added node
        tree.root = new_node
        if collapse_zero_length_nodes:
            tree.collapse_polytomies()
        self._tree = tree.to_newick(brlen_rounding)
        self._leaf_order = tree.leaf_names

    def get_condensed_idxs_of_sample(
            self,
            target_idx,
            nrows: int=None,
            to_skip=None
        ) -> npt.NDArray[np.uintp]:
        """For a row index, return the condensed indices of all cells in that row

        If to_skip is provided, skips the index corresponding to that value in the matrix
        """
        if nrows is None:
            nrows = self.nrows
        if to_skip is not None:
            idxs_count = nrows - 2
        else:
            idxs_count = nrows - 1

        idxs = np.zeros(idxs_count, dtype=np.uintp)
        i = 0 # counter to track where in the index array we are adding values
        total_offset = 0 # how far into the condensed matrix is the start of this row
        for rownum in range(target_idx+1):
            # how far into this row before the first value in the condensed matrix
            this_row_offset = rownum + 1
            # how many values from this row are in the condensed matrix
            this_row_len = nrows - (rownum + 1)
            if rownum == target_idx:
                idx_range = np.arange(total_offset, total_offset + this_row_len)
                if to_skip is not None and to_skip > target_idx:
                    # remove value to skip from range
                    idx_range = np.delete(idx_range, to_skip-this_row_offset)
                idxs[i:i+len(idx_range)] = idx_range
                i += len(idx_range)
            else:
               # Check if this row corresponds to the value to skip
               if rownum != to_skip:
                idxs[i]  = total_offset + target_idx - this_row_offset
                i += 1
            total_offset += this_row_len
        return idxs


    def neighbor_joining(self):
        brlen_rounding = self._algorithm_settings.get("brlen_rounding", None)
        tree = NJ.from_matrix(self)
        self._tree = tree.to_newick(brlen_rounding)
        self._leaf_order = tree.leaf_names

    def minimum_spanning(self, **_):
        """Prims algorithm"""
        # MS tree is a non-directional tree where every node is a leaf
        # need a node lookup rather than directly indexing tree nodes
        # because identical samples will be added as a polytomy descended from one sample
        network = Tree.empty()
        node_lookup = []
        for sample in self.samples:
            leaves = [str(l) for l in self.redundants[sample]]
            if len(leaves) == 1:
                this_id = network.add_node(None, leaves[0], 0.)
            else:
                this_id = network.add_leaf_polytomy(
                    internal_name=leaves[0], # use representative as name for internal node
                    leaf_names=leaves[1:],
                )
            node_lookup.append(this_id)
        # track which nodes have already been joined to something
        node_idxs = set()
        # pick a starting point
        start_idx = 0
        network.root = node_lookup[start_idx] # root arbitrarily on first sample
        node_idxs.add(start_idx)

        # Set up masks for use in identifying next node to add to network
        # Use two masks to keep things straight.
        # joined_mask will mask edges between already joined nodes in distmat
        # candidate_mask will mask edges that do not include a node currently in the matrix
        # joined_mask will mask more and more edges, while chandidate_mask will mask fewer
        # as the algorithm proceeds
        joined_mask = np.array([False] * self._distmat.shape[0]) # start with none masked
        candidate_mask = np.array([True] * self._distmat.shape[0]) # start with all masked
        # Choose candidate to start with and unmask its edges
        start_edges = [self.condensed_idx(i, start_idx) for i in range(self.nrows) if i != start_idx]
        candidate_mask[start_edges] = False

        while len(node_idxs) < len(self.samples):
            mask = np.ma.mask_or(joined_mask, candidate_mask)
            masked_array = np.ma.masked_array(self._distmat, mask=mask)
            # Choose node to add to
            min_dist_idx = np.argmin(masked_array)
            weight = self._distmat[min_dist_idx]
            idx_a, idx_b = self.square_idx(min_dist_idx)
            # determine which node is new
            if idx_a in node_idxs:
                existing_idx, partner_idx = idx_a, idx_b
            else:
                existing_idx, partner_idx = idx_b, idx_a

            # Add new node as child of node already in the tree.
            # Direction doesn't matter here, "child" is only meaningful for the other
            # use of the Tree class.
            # Adding new nodes as children of existing nodes means that the parent
            # attr is only set once (when the node is added) so parent-child
            # relationships can be used to easily traverse the tree structure
            # when producing a newick
            existing = network.nodes[node_lookup[existing_idx]]
            new = network.nodes[node_lookup[partner_idx]]
            existing.children.append(new.id)
            new.parent = existing.id
            new.branch_length = weight
            # mask edges between added node and nodes already in the network
            joined_mask[
                    [self.condensed_idx(existing_idx, partner_idx) for existing_idx in node_idxs]
                ] = True
            # unmask edges from added node to nodes not yet in the network
            mask_idxs = self.get_condensed_idxs_of_sample(partner_idx)
            candidate_mask[mask_idxs] = False

            node_idxs.add(partner_idx)
        brlen_rounding = self._algorithm_settings.get("brlen_rounding", None)
        self._tree = network.to_grapetree(brlen_rounding)

    def infer_tree(self):
        self.algorithm()

    def dmat_as_list(self, precision: int=None) -> list[list[float]]:
        distmat = self.as_square(self._distmat)
        distmat = self.add_redundant_samples_to_distmat(distmat)
        if self._tree is not None:
            distmat = self.reorder_distmat(distmat)

        return [[round(i, precision) for i in row] for row in distmat.tolist()]

    @classmethod
    def from_json(cls, json: dict, algorithm: str="upgma", distance: str="normalized_allele_differences", settings: dict=None):
        if settings is None:
            settings = {}
        dmat = cls()
        dmat.algorithm = algorithm
        dmat.distance_metric = distance

        dmat.read_profile_dict(json)
        dmat.nonredundant()
        dmat.symmetric_distance()
        dmat._algorithm_settings.update(settings)

        return dmat


    @classmethod
    def from_seqs(cls, seqs: dict[str, str], algorithm: str="upgma", distance: str="absolute_allele_differences", settings: dict=None):
        if settings is None:
            settings = {}
        dmat = cls()
        dmat.algorithm = algorithm
        dmat.distance_metric = distance

        dmat.read_alignment_dict(seqs)
        dmat.nonredundant()
        dmat.symmetric_distance()
        dmat._algorithm_settings.update(settings)

        return dmat


class NJ(CondensedMatrix):
    """Class to handle neighbor joining"""
    def __init__(self, dmat: np.ndarray, tree: Tree, node_ids: list[int]):
        self.tree = tree
        self.nodes = node_ids
        self.dmat = dmat
        self.nrows = int((1 + sqrt(1+ 8*dmat.shape[0])) / 2)
        self.init_qmat_dist_sums()
        self.joined = set()

    @classmethod
    def from_matrix(cls, distmat: DistanceMatrix, brlen_rounding = None) -> Tree:
        dmat = distmat.distmat
        tree = Tree.empty()
        node_ids = []
        for sample in distmat.samples:
            leaves = distmat.redundants[sample]
            # convert to str for use in tree
            leaves = [str(name) for name in leaves]
            if len(leaves) == 1:
                this_id = tree.add_node(None, leaves[0], 0.)
            else:
                this_id = tree.add_leaf_polytomy(leaves)
            node_ids.append(this_id)
        nj = cls(dmat, tree, node_ids)
        nj._run()
        return nj.tree

    def _run(self):
        self.join()
        for _ in range(self.nrows-4): # do n-4 joins (until three nodes left)
            self.find_best_pair()
            self.join()

        # join remaining 3 nodes as a polytomy
        # first get indices of remaining nodes - the two best matches (i,j) and the third (k)
        self.find_best_pair()
        i, j = self.join_idxs
        k = [k for k in range(self.nrows) if k not in self.joined and k not in [i, j]][0]
        # get distances for branched
        d_i, d_j = self.calculate_distance(i, j)
        d_k = 0.5 * (
            self.dmat[self.condensed_idx(k, i)]
            + self.dmat[self.condensed_idx(k, j)]
            - self.dmat[self.condensed_idx(i, j)]
        )
        # Add final branches to the tree and root it
        self.tree.nodes[self.nodes[i]].branch_length = d_i
        self.tree.nodes[self.nodes[j]].branch_length = d_j
        self.tree.nodes[self.nodes[k]].branch_length = d_k
        root = self.tree.add_node_with_children(
            children=[self.nodes[x] for x in [i, j, k]]
        )
        self.tree.root = root
        self.tree.midpoint_root()

    def join(self):
        # get distances
        x, y = self.join_idxs
        d_xu, d_yu = self.calculate_distance(x, y, disallow_negative=True)
        # combine nodes into new node
        self.tree.nodes[self.nodes[x]].branch_length = d_xu
        self.tree.nodes[self.nodes[y]].branch_length = d_yu
        self.nodes[x] = self.tree.add_node_with_children(
            children=[self.nodes[x], self.nodes[y]]
        )
        self.joined.add(y)
        self.update_dists(x, y)

    def init_qmat_dist_sums(self):
        array_size = self.nrows * (self.nrows -1) // 2
        self.qmat = np.empty(array_size, dtype=np.float32)
        self.q_min_idxs = np.empty(self.nrows, dtype=np.uintp)
        self.q_min_vals = np.empty(self.nrows, dtype=np.float32)
        self.dist_sums = np.zeros(self.nrows, dtype=np.float32)
        # init distance sums
        k = 0
        for i in range(self.nrows):
            for j in range(i+1, self.nrows):
                d = self.dmat[k]
                self.dist_sums[i] += d
                self.dist_sums[j] += d
                k += 1
        # init Q matrix
        k = 0
        lowest_q = np.finfo(np.float32).max
        for i in range(self.nrows):
            min_q = np.finfo(np.float32).max
            for j in range(i+1, self.nrows):
                # Q = (N-2)*d[i,j] - sum_dist[i] - sum_dist[j]
                q = (self.nrows-2)*self.dmat[k] - self.dist_sums[i] - self.dist_sums[j]
                self.qmat[k] = q
                if q < min_q:
                    min_q = q
                    self.q_min_idxs[i] = j
                    self.q_min_vals[i] = q
                if q < lowest_q:
                    lowest_q = q
                    self.join_idxs = (i, j)
                k += 1

    def update_dists(self, x: int, y: int):
        # Update dists and distance_sum for the new node
        self.dist_sums[x] = 0
        d_xy = self.dmat[self.condensed_idx(x, y)]
        for z in range(self.nrows):
            if z in self.joined or z == x or z == y:
                continue
            # Get indices of xz and yz and the correspnding distances
            xz = self.condensed_idx(x, z)
            yz = self.condensed_idx(y, z)
            d_xz = self.dmat[xz]
            d_yz = self.dmat[yz]
            # overwrite xz with new distance (i.e., row X becomes new node U)
            d_uz = 0.5 * (d_xz + d_yz - d_xy)
            self.dmat[xz] = d_uz
            # Subtract the old node distances from this row's sum and add in the new node distance
            self.dist_sums[z] += (d_uz - d_xz - d_yz)
            # add the new dist to the sum for node U
            self.dist_sums[x] += d_uz

    def calculate_distance(self, x, y, rounding: int=None, disallow_negative: bool=True):
        """Get the distance from X and Y to the new node, U"""
        # dAU = (0.5*dAB) + 1/(2*(num_rows-2))*(sum(d_mat[row])-sum(d_mat[col]))
        c_i = self.condensed_idx(x, y)
        d_xy = self.dmat[c_i]
        n_remaining = self.nrows - len(self.joined)
        d_xu = 0.5*(d_xy) + 1/(2*(n_remaining-2))*(self.dist_sums[x] - self.dist_sums[y])
        d_yu = d_xy - d_xu
        if disallow_negative:
            d_xu = max(0, d_xu)
            d_yu = max(0, d_yu)
        if rounding is not None:
            d_xu = round(d_xu, rounding)
            d_yu = round(d_yu, rounding)

        return d_xu, d_yu

    def find_best_pair(self):
        last_x, _ = self.join_idxs
        # Update the just-joined row
        best_q, y = self.check_row_q(last_x)
        # init with the previous row to be overwritten if a better Q found
        x = last_x
        # iterate over all rows of the distance matrix to find the new best Q
        for i in range(self.nrows):
            if i in self.joined:
                continue
            if i == last_x:
                # just checked this row above
                continue

            if i < last_x:
                # this row of the triangle includes a column for the
                # just joined node so we need to update that column
                sum_d_i = self.dist_sums[i]
                sum_d_lastx = self.dist_sums[last_x]
                d_i_lastx = self.dmat[self.condensed_idx(i, last_x)]
                n_remaining = self.nrows - len(self.joined)
                q = (n_remaining-2)*d_i_lastx - sum_d_i - sum_d_lastx
                if q < self.q_min_vals[i]:
                    # If the updated value is better than the previous
                    # best for this column then update it
                    self.q_min_vals[i] = q
                    self.q_min_idxs[i] = last_x

            if self.q_min_vals[i] > best_q:
                # this row's best Q is worse than the best we've seen so skip
                continue

            # Otherwise, this row could container a better Q. Search the whole row
            q, j = self.check_row_q(i)

            if q < best_q:
                best_q = q
                x = i
                y = j
        if x > y:
            x, y = y, x
        self.join_idxs = (x, y)

    def check_row_q(self, row):
        x = row
        y = x+1
        q_xy = np.finfo(np.float32).max
        sum_d_x = self.dist_sums[x]
        offset = self.condensed_offset(x) # calc offset and add j each iter rather than calculating every iter
        for j in range(x+1, self.nrows):
            if j in self.joined:
                # Not considering this column moving forward
                continue

            sum_d_j = self.dist_sums[j]
            c_idx = offset + j
            d_xj = self.dmat[c_idx]
            n_remaining = self.nrows - len(self.joined)
            q = (n_remaining-2)*d_xj - sum_d_x - sum_d_j

            if q < q_xy:
                q_xy = q
                y = j
        # store the best q value and its index
        self.q_min_idxs[x] = y
        self.q_min_vals[x] = q_xy
        return q_xy, y
