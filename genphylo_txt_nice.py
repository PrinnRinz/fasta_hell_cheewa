from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
from Bio import SeqIO

def parse_species_names(fasta_file):
    species_mapping = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        identifier = record.id  
        description = record.description
        species_name = description.split(":")[-1].strip()
        species_mapping[identifier] = species_name
    return species_mapping

alignment_file = "test1.aln"
fasta_file = "test1.fasta"
alignment = AlignIO.read(alignment_file, "clustal")

calculator = DistanceCalculator("identity")
distance_matrix = calculator.get_distance(alignment)
constructor = DistanceTreeConstructor(calculator, method="upgma")
tree = constructor.build_tree(alignment)

species_mapping = parse_species_names(fasta_file)
for clade in tree.find_clades():
    if clade.name in species_mapping:
        clade.name = species_mapping[clade.name]

print("\nPhylogenetic Tree:")
Phylo.draw_ascii(tree)

with open("phylogenetic_tree_with_species.txt", "w") as file:
    Phylo.draw_ascii(tree, file=file)

