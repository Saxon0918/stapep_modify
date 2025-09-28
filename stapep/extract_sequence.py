import os
import pandas as pd
from Bio import PDB
import params


def extract_sequence_from_pdb(pdb_path):
    """从不含chain ID的PDB文件中提取三字母残基名并转为一字母序列"""
    sequence = []
    seen_residues = set()

    try:
        with open(pdb_path, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    res_name = line[17:20].strip()
                    res_id = (line[22:26].strip(), line[21].strip())  # resSeq 和 chain ID

                    if res_id not in seen_residues:
                        seen_residues.add(res_id)
                        one_letter = params.amino_acid_dict.get(res_name, 'X')  # 未知残基记为 X
                        sequence.append(one_letter)
    except Exception as e:
        print(f"Error reading {pdb_path}: {e}")
        return ''

    return ''.join(sequence)


# 设置根目录路径和CSV文件路径
root_dir = '/home/d3008/Documents/zhr/gsk3beta/'  # 或者是你的实际根目录路径
csv_path = os.path.join(root_dir, 'result_final.csv')

# 读取CSV文件
df = pd.read_csv(csv_path, sep=',')  # 可能是tab分隔，也可能是逗号，按需修改

# 遍历每一行并提取序列
sequences = []
for index, row in df.iterrows():
    name = row['Name']
    folder_lvl1 = name.split('_')[0] + '_' + name.split('_')[1]  # e.g., 1gng_X
    subfolder_path = os.path.join(root_dir, folder_lvl1, 'openbpmd', name)
    pdb_file_path = os.path.join(subfolder_path, 'complex_B_NCACO.pdb')

    if os.path.exists(pdb_file_path):
        seq = extract_sequence_from_pdb(pdb_file_path)
    else:
        seq = ''  # 文件不存在时设为空

    sequences.append(seq)

# 更新DataFrame并保存
df['Sequence'] = sequences
df.to_csv('/home/d3008/Documents/zhr/gsk3beta/result_final_with_sequences.csv', index=False, sep=',')

