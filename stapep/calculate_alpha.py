import os
import sys

from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.PDB import DSSP
from Bio.PDB.Polypeptide import is_aa
import pymol
from pymol import cmd

def calculate_alpha(pdb_file_path):
    parser = PDBParser()
    ligand_structure = parser.get_structure("ligand", pdb_file_path)
    model = ligand_structure[0]
    dssp = DSSP(model, pdb_file_path)

    alpha_helix_count = sum(1 for res in dssp if res[2] == 'H')  # alpha-helix
    beta_sheet_count = sum(1 for res in dssp if res[2] == 'E')  # beta-sheet
    total_res = len(dssp)
    helix_ratio = alpha_helix_count / total_res
    return helix_ratio


def get_chain_atoms(structure, chain_id=None):
    """获取指定链的所有 alpha carbon 原子列表（用于对齐）"""
    atoms = []
    for model in structure:
        for chain in model:
            if chain_id is None or chain.id == chain_id:
                for residue in chain:
                    if is_aa(residue) and 'CA' in residue:
                        atoms.append(residue['CA'])
                if chain_id is not None:
                    return atoms
    return atoms


def pymol_align(protein_path, ligand_root_path, filenames, ouput_path):
    pymol.finish_launching(['pymol', '-cq'])
    cmd.load(protein_path, 'protein')

    for ligand_name in filenames:
        ligand_path = os.path.join(ligand_root_path, ligand_name)
        ligand_name = os.path.splitext(os.path.basename(ligand_name))[0]
        cmd.load(ligand_path, ligand_name)
        cmd.align(f"{ligand_name} and name CA", "protein and chain X and name CA")
        cmd.save(f"{ouput_path}/{ligand_name}.pdb", ligand_name)
    cmd.quit()



if __name__ == '__main__':
    ligand_name = "1gng_1"
    protein_path = "/home/d3008/Documents/zhr/pycharmprojects/RFdiffusion/examples/input_pdbs/1gng.pdb"
    ligand_root_path = f"example/data/gsk3beta/{ligand_name}/denove"
    filenames = sorted(f for f in os.listdir(ligand_root_path) if f.endswith(".pdb"))
    align_filenames = []
    for filename in filenames:
        ligand_path = os.path.join(ligand_root_path, filename)
        helix_ratio = calculate_alpha(ligand_path)
        if helix_ratio > 0.3:
            align_filenames.append(filename)

    ligand_output_path = f"example/data/gsk3beta/{ligand_name}/align_denove"
    if not os.path.exists(ligand_output_path):
        os.makedirs(ligand_output_path)
    pymol_align(protein_path, ligand_root_path, align_filenames, ligand_output_path)
