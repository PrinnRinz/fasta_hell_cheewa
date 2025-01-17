from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

alignment_file = "test1.aln"
alignment = AlignIO.read(alignment_file, "clustal")

calculator = DistanceCalculator("identity")
distance_matrix = calculator.get_distance(alignment)

constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)

output_file = "phylogenetic_tree.xml"
Phylo.write(tree, output_file, "phyloxml")

print(f"Phylogenetic tree saved as {output_file}")

