from Bio.PDB.PDBList import PDBList
from Bio.PDB.Structure import Structure
from Bio import PDB


def download_pdb(pdb_id: str, output_dir: str):
    pdbl = PDBList()
    pdbl.retrieve_pdb_file(pdb_id, pdir=output_dir, file_format="pdb")
    
def get_structure(pdb_id: str) -> Structure:
    download_pdb(pdb_id, "./tmp")
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(pdb_id, f"./tmp/pdb{pdb_id.lower()}.ent")
    return structure
