import argparse
from Bio import PDB
from pathlib import Path
import numpy as np
from consts import *

def load_template(template_dir, residue_name):
    pdb_parser = PDB.PDBParser()
    template_path = Path(template_dir)/f"{residue_name}.pdb"
    structure = pdb_parser.get_structure(residue_name, template_path)
    return structure

def superimpose(fixed_atoms, template_atoms):

    if len(fixed_atoms) != len(template_atoms):
        raise ValueError("Mismatched atom counts.")

    sup = PDB.Superimposer()
    sup.set_atoms(fixed_atoms, template_atoms)
    sup.apply(template_atoms)

    return sup.rotran


def rebuild_structure(input_path, output_path, template_dir):
    pdb_parser = PDB.PDBParser()
    io = PDB.PDBIO()

    structure = pdb_parser.get_structure("cg_structure", input_path)
    output_structure = PDB.Structure.Structure("restored_structure")

    processed_residues = 0

    for model in structure:
        output_model = PDB.Model.Model(model.id)
        for chain in model:
            output_chain = PDB.Chain.Chain(chain.id)
            for residue in chain:
                if residue.resname in [*PYRIMIDINES , *PURINES]:
                    template_structure = load_template(template_dir, residue.resname)
                    template_residue = next(template_structure.get_residues())

                    if residue.resname in PURINES:
                        coarse_atoms = PURINE_ATOMS
                    elif residue.resname in PYRIMIDINES:
                        coarse_atoms = PYRIMIDINE_ATOMS
                    else:
                        continue

                    coarse_atoms.extend(["P", "C4'", "O5'", "C3'", "O3'"])

                    fixed_atoms = []
                    for atom in residue.get_atoms():
                        if atom.get_name() in coarse_atoms:
                            fixed_atoms.append(atom)

                    template_atoms = []
                    for atom in template_residue.get_atoms():
                        if atom.get_name() in coarse_atoms:
                            template_atoms.append(atom)

                    if len(fixed_atoms) != len(template_atoms):
                        matched_atoms = min(len(fixed_atoms), len(template_atoms))
                        fixed_atoms = fixed_atoms[:matched_atoms]
                        template_atoms = template_atoms[:matched_atoms]

                    if len(fixed_atoms) < 3:
                        continue

                    rotran = superimpose(fixed_atoms, template_atoms)

                    new_residue = PDB.Residue.Residue(residue.id, residue.resname, residue.segid)
                    for atom in template_residue.get_atoms():
                        new_atom = atom.copy()
                        new_atom.transform(rotran[0], rotran[1])
                        new_residue.add(new_atom)

                    for atom in residue.get_atoms():
                        if atom.get_name() not in coarse_atoms:
                            new_residue.add(atom.copy())

                    output_chain.add(new_residue)
                    processed_residues += 1
                else:
                    continue

            output_model.add(output_chain)
        output_structure.add(output_model)

    io.set_structure(output_structure)
    io.save(output_path)

    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help="Coarse-grained PDB input file")
    parser.add_argument('--output', required=True, help="Restored PDB output file")

    args = parser.parse_args()

    rebuild_structure(args.input, args.output, './templates')

if __name__ == "__main__":
    main()