import logging
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
import os

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("alignment_debug.log"),
        logging.StreamHandler()
    ]
)

fasta_file = "test1.fasta" 
clustalw_exe = "/usr/bin/clustalw2"

logging.info("Starting alignment script.")

if not os.path.exists(fasta_file):
    logging.error(f"File not found: {fasta_file}")
    raise FileNotFoundError(f"File not found: {fasta_file}")

logging.info("Get Species from FASTA")
species_names = []
for record in SeqIO.parse(fasta_file, "fasta"):
    if ":" in record.description:
        species_name = record.description.split(":", 1)[1].strip()
    else:
        species_name = record.description.strip()
    species_names.append(species_name)
logging.debug(f"FASTA Species: {species_names}")

logging.info("Starting ClustalW.")
aligned_file = "test1.aln"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile=fasta_file)
try:
    stdout, stderr = clustalw_cline()
    logging.debug(f"ClustalW stdout: {stdout}")
    logging.debug(f"ClustalW stderr: {stderr}")
    logging.info(f"Alignment file created: {aligned_file}")
except Exception as e:
    logging.exception("Error during sequence alignment.")
    raise

