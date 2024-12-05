from pathlib import Path
from get_structure import get_structure
import argparse
from Bio.PDB import *


def save_residue_as_pdb(residue, output_path):
    structure = Structure.Structure('template')
    model = Model.Model(0)
    chain = Chain.Chain('A')
    
    new_residue = residue.copy()
    chain.add(new_residue)
    model.add(chain)
    structure.add(model)
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))

def create_rna_templates(structure_pdb_id, output_dir):
    structure = get_structure(structure_pdb_id)
    A = None
    U = None
    C = None
    G = None
    for model in structure:
        for chain in model:
            for res in chain:
                if not A and res.get_resname() == 'A':
                    atom_names = [atom.get_name() for atom in res.get_atoms()]
                    if 'P' in atom_names and 'C4\'' in atom_names:
                        A = res
                if not U and res.get_resname() == 'U':
                    atom_names = [atom.get_name() for atom in res.get_atoms()]
                    if 'P' in atom_names and 'C4\'' in atom_names:
                        U = res
                if not G and res.get_resname() == 'G':
                    atom_names = [atom.get_name() for atom in res.get_atoms()]
                    if 'P' in atom_names and 'C4\'' in atom_names:
                        G = res
                if not C and res.get_resname() == 'C':
                    atom_names = [atom.get_name() for atom in res.get_atoms()]
                    if 'P' in atom_names and 'C4\'' in atom_names:
                        C = res

    output_dir = Path(output_dir)
    if A:
        save_residue_as_pdb(A, output_dir / 'A.pdb')
    if U:
        save_residue_as_pdb(U, output_dir / 'U.pdb')
    if G:
        save_residue_as_pdb(G, output_dir / 'G.pdb')
    if C:
        save_residue_as_pdb(C, output_dir / 'C.pdb')
                        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb-id', required=True, help="PDB ID of RNA to use for template creation")
    parser.add_argument('--output-dir', required=True, help="Output directory for template PDB files")

    args = parser.parse_args()

    create_rna_templates(args.pdb_id, args.output_dir)

if __name__ == "__main__":
    main()