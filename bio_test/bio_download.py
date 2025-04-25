from Bio.PDB import PDBList
import torch
from Bio import PDB
import numpy as np
import os


def extract_sequence_and_coordinates(pdb_file):
    """
    提取PDB文件中的序列和坐标。
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("PDB_structure", pdb_file)

    chains = []
    sequences = []
    coordinates = []

    for model in structure:
        for chain in model:
            chain_id = chain.get_id()
            chains.append(chain_id)
            seq = []
            coords = []
            for residue in chain:
                # 获取氨基酸序列
                if PDB.is_aa(residue):
                    seq.append(residue.get_resname())
                    # 获取3D坐标
                    coords.append(residue['CA'].get_coord())  # 使用 Cα 原子的坐标
            sequences.append(seq)
            coordinates.append(np.array(coords))

    return chains, sequences, coordinates


def generate_pt_file(pdb_file, output_pt_file):
    """
    从PDB文件生成.pt文件
    """
    # 提取PDB中的序列和坐标
    chains, sequences, coordinates = extract_sequence_and_coordinates(pdb_file)

    # 假设只有一个组装体，且没有变换矩阵信息
    # 如果有多个组装体和变换矩阵，根据需求添加
    asmb_ids = ['1']  # 这里只用了一个组装体，若有多个可以扩展
    asmb_chains = [",".join(chains)]  # 假设所有链属于同一个组装体
    asmb_xform = np.zeros((1, 4, 4))  # 没有坐标变换矩阵，矩阵大小为1x4x4

    # 创建输出数据字典
    pt_data = {
        'seq': sequences[0],  # 假设第一个结构作为序列
        'xyz': torch.tensor(coordinates[0], dtype=torch.float32),  # 转为torch tensor
        'chains': chains,
        'asmb_ids': asmb_ids,
        'asmb_chains': asmb_chains,
        'asmb_xform': torch.tensor(asmb_xform, dtype=torch.float32),  # 默认没有变换矩阵
        'tm': np.array([[1.0]]),  # 假设序列相似度为1（实际情况下应根据同源性计算）
    }

    # 保存为.pt文件
    torch.save(pt_data, output_pt_file)
    print(f"PT file saved to {output_pt_file}")

# 下载pdb文件
# pdb_list = PDBList()
# pdb_list.retrieve_pdb_file("5H7K", pdir="pdb_files", file_format="pdb")

# 使用示例
# pdb_file = '../stapep/example/data_20/denovo_R8_S5_5_0.pdb'
pdb_file = "pdb_files/pdb5h7k.ent"  # 输入你的PDB文件路径
output_pt_file = "pt_files/pdb5h7k.pt"  # 输出生成的.pt文件路径
generate_pt_file(pdb_file, output_pt_file)
