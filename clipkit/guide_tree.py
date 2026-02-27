from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator,
    DistanceTreeConstructor,
    NNITreeSearcher,
    ParsimonyScorer,
    ParsimonyTreeConstructor,
)


def build_parsimony_guide_tree(alignment: MultipleSeqAlignment) -> Tree:
    """
    Build a parsimony guide tree for heterotachy-aware trimming.

    Uses an NJ tree from identity distances as a start tree, followed by
    NNI-based parsimony optimization.
    """
    distance_calculator = DistanceCalculator("identity")
    distance_matrix = distance_calculator.get_distance(alignment)
    nj_tree = DistanceTreeConstructor().nj(distance_matrix)

    parsimony_scorer = ParsimonyScorer()
    tree_searcher = NNITreeSearcher(parsimony_scorer)
    parsimony_constructor = ParsimonyTreeConstructor(tree_searcher, nj_tree)
    return parsimony_constructor.build_tree(alignment)
