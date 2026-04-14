import os
import csv
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from rdkit import Chem

# =========================
# 路径配置
# =========================
BASE_DIR = "/home/d3008/Documents/zhr/PDL1/PDL1_500_2"
OUTPUT_CSV = "/home/d3008/Documents/zhr/PDL1/PDL1_500_2/sting_sequences.csv"

# =========================
# 提取序列
# =========================
def extract_sequence(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", pdb_file)
    
    seq = ""
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":
                    try:
                        seq += seq1(residue.resname)
                    except:
                        seq += "X"
    return seq


# =========================
# PDB → SMILES
# =========================
def pdb_to_smiles(pdb_file):
    try:
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol is None:
            return "ERROR"

        # sanitize（很重要）
        Chem.SanitizeMol(mol)

        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        return smiles

    except Exception as e:
        print(f"[SMILES ERROR] {pdb_file}: {e}")
        return "ERROR"


# =========================
# 主函数
# =========================
def main():
    results = []

    for i in range(500):
        folder = os.path.join(BASE_DIR, f"pdl1_{i}", "cealign_denovo")

        if not os.path.exists(folder):
            continue

        for file in os.listdir(folder):
            if file.endswith(".pdb"):
                pdb_path = os.path.join(folder, file)

                print(f"[PROCESS] {pdb_path}")

                seq = extract_sequence(pdb_path)
                smiles = pdb_to_smiles(pdb_path)

                results.append([
                    file,
                    pdb_path,
                    seq,
                    smiles
                ])

    # 写入 CSV
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file_name", "full_path", "sequence", "smiles"])
        writer.writerows(results)

    print(f"\n✅ DONE! Saved to: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()

