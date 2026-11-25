from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
import numpy as np


def _alignment_fitch_scorer(alignment: MultipleSeqAlignment):
    from Bio.Phylo.TreeConstruction import ParsimonyScorer

    if not isinstance(alignment, MultipleSeqAlignment):
        return ParsimonyScorer()
    width = alignment.get_alignment_length()
    if not 512 <= len(alignment) * width <= 2_000_000:
        return ParsimonyScorer()
    if not all(record.seq.defined for record in alignment):
        return ParsimonyScorer()

    residues = np.array([bytes(record.seq) for record in alignment], dtype=f"S{width}")
    residues = residues.view(np.uint8).reshape(len(alignment), width)
    alphabet = np.unique(residues)
    if len(alphabet) > 64 or alphabet[-1] > 127:
        return ParsimonyScorer()
    state_dtype = np.min_scalar_type((1 << len(alphabet)) - 1)
    lookup = np.zeros(128, dtype=state_dtype)
    lookup[alphabet] = np.left_shift(
        state_dtype.type(1), np.arange(len(alphabet), dtype=state_dtype)
    )
    leaf_states = dict(zip((record.id for record in alignment), lookup[residues]))

    class BatchedFitchScorer(ParsimonyScorer):
        def get_score(self, tree, alignment):
            if not tree.is_bifurcating():
                return super().get_score(tree, alignment)
            if not tree.rooted:
                tree.root_at_midpoint()
            terminals = sorted(tree.get_terminals(), key=lambda term: term.name)
            alignment.sort()
            if not all(term.name == record.id for term, record in zip(terminals, alignment)):
                raise ValueError(
                    "Taxon names of the input tree should be the same with the alignment."
                )
            states = {term: leaf_states[term.name] for term in terminals}
            score = 0
            for clade in tree.get_nonterminals(order="postorder"):
                left, right = clade.clades[:2]
                lhs, rhs = states.pop(left), states.pop(right)
                # Fitch uses the intersection when nonempty, and charges one
                # step at each site where it must use the union.
                intersection = lhs & rhs
                disjoint = intersection == 0
                score += np.count_nonzero(disjoint)
                np.bitwise_or(lhs, rhs, out=intersection, where=disjoint)
                states[clade] = intersection
            return score

    return BatchedFitchScorer()


def build_parsimony_guide_tree(alignment: MultipleSeqAlignment) -> Tree:
    """
    Build a parsimony guide tree for heterotachy-aware trimming.

    Uses an NJ tree from identity distances as a start tree, followed by
    NNI-based parsimony optimization.
    """
    from Bio.Phylo.TreeConstruction import (
        DistanceCalculator,
        DistanceTreeConstructor,
        NNITreeSearcher,
        ParsimonyTreeConstructor,
    )

    distance_calculator = DistanceCalculator("identity")
    distance_matrix = distance_calculator.get_distance(alignment)
    nj_tree = DistanceTreeConstructor().nj(distance_matrix)

    parsimony_scorer = _alignment_fitch_scorer(alignment)
    tree_searcher = NNITreeSearcher(parsimony_scorer)
    parsimony_constructor = ParsimonyTreeConstructor(tree_searcher, nj_tree)
    return parsimony_constructor.build_tree(alignment)
