import os
import torch
from Bio.PDB import PDBParser
import numpy as np
import pandas as pd
import torch
import re
from Bio import PDB

def pdb_to_pt(input_pdb, output_pt, id, seq):
    """
    Convert a PDB file to ProteinMPNN-compatible .pt file.

    Args:
        input_pdb (str): Path to the input PDB file.
        output_pt (str): Path to save the output .pt file.
        id (str): Unique identifier for the structure.
        seq (str): Sequence of the protein from the PDB file.
    """
    # Load the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(id, input_pdb)

    # Extract atomic coordinates for chains
    chains = []
    coords = []
    for model in structure:
        for chain in model:
            # chain_id = chain.get_id()
            chain_id = 'A'
            chains.append(chain_id)
            chain_coords = []
            for residue in chain:
                if "CA" in residue:
                    atom = residue["CA"]
                    chain_coords.append(atom.coord)
            coords.append(torch.tensor(chain_coords, dtype=torch.float32))

    # Dummy transformation matrix (identity)
    asmb_xform = torch.eye(4).unsqueeze(0)
    asmb_method = ["?"]
    asmb_details = ["author_defined_assembly"]

    # Define the metadata for .pt file
    data = {
        "method": "X-RAY_DIFFRACTION",  # Modify if needed
        "date": "2004-12-30",  # Modify with actual date if known
        "resolution": 1.7,  # Update with actual resolution if available
        "chains": chains,
        "seq": [[seq, seq]],  # Ensure the sequence matches
        "id": id,
        "asmb_chains": [",".join(chains)],
        "asmb_details": asmb_details,
        "asmb_method": asmb_method,
        "asmb_ids": ["1"],  # Default assembly ID
        "asmb_xform0": asmb_xform,
        "tm": torch.ones((1, 1, 3)),  # Dummy translation matrix
    }

    # Save as .pt file
    torch.save(data, output_pt)
    print(f"Saved ProteinMPNN-compatible file: {output_pt}")


def pdb_to_pt_chain(pdb_file, output_file):
    # Initialize data structures
    seq = ""
    xyz_list = []
    bfac_list = []
    occ_list = []
    mask_list = []

    # Define a mapping from three-letter to one-letter amino acid codes
    three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',  # standard amino acids
                       'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',  # standard amino acids
                       'GLY': 'G', 'HIS': 'H', 'HIE': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',  # standard amino acids
                       'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',  # standard amino acids
                       'PS3': 'S3', 'PS5': 'S5', 'PS8': 'S8',  # stapled amino acids
                       'PR3': 'R3', 'PR5': 'R5', 'PR8': 'R8',  # stapled amino acids
                       'NLE': 'B'  # non-standard amino acids
                       }

    # Read the PDB file
    with open(pdb_file, "r") as file:
        for line in file:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                resseq = line[22:26].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occupancy = float(line[54:60].strip())
                bfactor = float(line[60:66].strip())

                # Initialize xyz entry with 4 slots for N, CA, C, O atoms
                if len(xyz_list) == 0 or resseq != xyz_list[-1][0]:
                    xyz_list.append([resseq, np.full((4, 3), np.nan)])
                    occ_list.append(np.full(4, 0.0))
                    bfac_list.append(np.full(4, 0.0))
                    mask_list.append(np.zeros(4))

                current_residue = xyz_list[-1]

                # Assign coordinates and other data for specific atom types
                atom_indices = {"N": 0, "CA": 1, "C": 2, "O": 3}
                if atom_name in atom_indices:
                    idx = atom_indices[atom_name]
                    current_residue[1][idx] = [x, y, z]
                    occ_list[-1][idx] = occupancy
                    bfac_list[-1][idx] = bfactor
                    mask_list[-1][idx] = 1.0

                # Convert three-letter code to one-letter code
                if atom_name == "CA":
                    if resname in three_to_one:
                        seq += three_to_one[resname]
                    else:
                        seq += "X"  # Unknown residues

    # Convert lists to tensors
    xyz_tensor = torch.tensor([residue[1] for residue in xyz_list], dtype=torch.float32)
    bfac_tensor = torch.tensor(bfac_list, dtype=torch.float32)
    occ_tensor = torch.tensor(occ_list, dtype=torch.float32)
    mask_tensor = torch.tensor(mask_list, dtype=torch.float32)

    # Create output dictionary
    chain_data = {
        "A": {  # Default chain ID is 'A'
            "seq": seq,
            "xyz": xyz_tensor,
            "mask": mask_tensor,
            "bfac": bfac_tensor,
            "occ": occ_tensor
        }
    }

    # Save to PT file
    torch.save(chain_data, output_file)
    print(f"Chain data saved to {output_file}")

# Example Usage
if __name__ == "__main__":
    # dataset_dir = "example/stapep_data/stapep292"
    # csv_file = "example/datasets/Filtered_peptides.csv"
    # dataset_dir = "example/stapep_data/amp"
    # csv_file = "example/datasets/filter_stapep_amp.csv"
    dataset_dir = "example/stapep_data/pdb_S5_sup_4986"
    csv_file = "example/datasets/list.csv"
    # 加载 CSV 文件
    csv_data = pd.read_csv(csv_file)

    # 获取并排序所有的 .pdb 文件
    pdb_files = sorted(
        [file_name for file_name in os.listdir(dataset_dir) if file_name.endswith(".pdb")]
    )
    # 检查是否存在足够的行数
    if len(csv_data) >= len(pdb_files):
        # 遍历文件列表，基于索引 i
        for i in pdb_files:
            input_pdb_path = os.path.join(dataset_dir, i)

            _i = i.replace("_", "")
            pt_file_name = _i.replace(".pdb", ".pt")
            output_pt_path = f"example/stapep_data/stapep_pt_S5/{pt_file_name}"
            protein_id = i
            # num = int(re.search(r"stapep_(\d+)\.pdb", i).group(1))
            # protein_seq = csv_data.loc[num, "SEQUENCE"]

            num = re.sub(r"(stapep\d+)\.pdb", r"\1_A", _i)
            protein_seq = csv_data[csv_data['CHAINID'] == num]['SEQUENCE'].values[0]
            pdb_to_pt(input_pdb_path, output_pt_path, protein_id, protein_seq)

            pt_file_name_A = _i.replace(".pdb", "_A.pt")
            output_pt_path_A = f"example/stapep_data/stapep_pt_S5_chain/{pt_file_name_A}"
            pdb_to_pt_chain(input_pdb_path, output_pt_path_A)
    else:
        print("Error: The number of rows in the CSV file is less than the number of .pdb files.")



    # pdb_to_pt_chain(input_pdb_path, output_pt_path)
