import logging
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import os

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("debug.log"),
        logging.StreamHandler()
    ]
)

# Define paths
fasta_file = "test1.fasta"  # Replace with your FASTA file
clustalw_exe = "/usr/bin/clustalw2"  # Update this to your ClustalW installation path

logging.info("Starting phylogenetic tree construction script.")

if not os.path.exists(fasta_file):
    logging.error(f"File not found: {fasta_file}")
    raise FileNotFoundError(f"File not found: {fasta_file}")

# Step 1: Extract species names from the FASTA file
logging.info("Extracting species names from the FASTA file.")
species_names = []
for record in SeqIO.parse(fasta_file, "fasta"):
    # Extract the portion of the header after the colon and remove any leading/trailing spaces
    if ":" in record.description:
        species_name = record.description.split(":", 1)[1].strip()
    else:
        # Fallback in case there's no colon
        species_name = record.description.strip()
    species_names.append(species_name)
logging.debug(f"Extracted species names: {species_names}")

# Step 2: Align sequences using ClustalW
logging.info("Aligning sequences using ClustalW.")
aligned_file = "test1.aln"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile=fasta_file)
try:
    stdout, stderr = clustalw_cline()
    logging.debug(f"ClustalW stdout: {stdout}")
    logging.debug(f"ClustalW stderr: {stderr}")
except Exception as e:
    logging.exception("Error during sequence alignment.")
    raise

# Step 3: Read the alignment
logging.info("Reading the alignment file.")
try:
    alignment = AlignIO.read(aligned_file, "clustal")
    logging.debug(f"Alignment loaded: {alignment}")
except Exception as e:
    logging.exception("Error reading alignment file.")
    raise

# Step 4: Calculate distance matrix
logging.info("Calculating distance matrix.")
try:
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
    logging.debug(f"Distance matrix: {distance_matrix}")
except Exception as e:
    logging.exception("Error calculating distance matrix.")
    raise

# Step 5: Build the phylogenetic tree using the distance matrix
logging.info("Constructing the phylogenetic tree.")
try:
    constructor = DistanceTreeConstructor(calculator, method="nj")  # Use 'upgma' for UPGMA
    phylo_tree = constructor.build_tree(alignment)
    phylo_tree.ladderize()  # Arrange branches for better visualization

except Exception as e:
    logging.exception("Error constructing the phylogenetic tree.")
    raise

# Step 6: Save the tree
logging.info("Saving the phylogenetic tree.")
try:
    Phylo.write(phylo_tree, "phylogenetic_tree_with_species.xml", "phyloxml")
    logging.info("Phylogenetic tree saved as 'phylogenetic_tree_with_species.xml'.")
except Exception as e:
    logging.exception("Error saving the phylogenetic tree.")
    raise

