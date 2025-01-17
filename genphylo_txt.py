
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

alignment_file = "test1.aln"
alignment = AlignIO.read(alignment_file, "clustal")

calculator = DistanceCalculator("identity")
distance_matrix = calculator.get_distance(alignment)

constructor = DistanceTreeConstructor(calculator, method="upgma")
tree = constructor.build_tree(alignment)

with open("phylogenetic_tree.txt", "w") as file:
    Phylo.draw_ascii(tree, file=file)

