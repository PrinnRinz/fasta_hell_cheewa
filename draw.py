from Bio import Phylo
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

tree_file = "phylogenetic_tree.xml" 
logging.info(f"Loading the phylogenetic tree from '{tree_file}'.")

try:
    phylo_tree = Phylo.read(tree_file, "phyloxml")
    logging.info("Phylogenetic tree loaded successfully.")
except Exception as e:
    logging.exception("Error loading the phylogenetic tree.")
    raise

logging.info("Customizing tree labels to show only species names.")
for clade in phylo_tree.find_clades():
    if clade.name:
        clade.name = clade.name.split("_")[0]

logging.info("Visualizing the phylogenetic tree with simplified labels.")
Phylo.draw(phylo_tree)  

