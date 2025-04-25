import pandas as pd

# # 读取原始文件
# input_file = "example/datasets/Stapled-peptide_permeability_filtered.csv"
# output_file = "example/datasets/Filtered_peptides.csv"
# target_elements = ['S3', 'S5', 'S8', 'R3', 'R5', 'R8']
#
# # 加载CSV文件
# df = pd.read_csv(input_file)
#
# # 确保第一列为 SEQUENCE，并提取包含目标元素的行
# filtered_df = df[df['SEQUENCE'].apply(lambda seq: any(element in seq for element in target_elements))]
#
# # 格式化SEQUENCE列，去除指定符号
# def format_sequence(seq):
#     # 去除指定符号
#     for char in ["Ac", "-", "NH2", "(EEA)"]:
#         seq = seq.replace(char, "")
#     return seq
#
# filtered_df['SEQUENCE'] = filtered_df['SEQUENCE'].apply(format_sequence)
#
# # 将过滤后的数据保存到新文件
# filtered_df.to_csv(output_file, index=False)
#
# print(f"已成功提取并保存至 {output_file}")


import pandas as pd

# 定义氨基酸字典
amino_acid_dict = {'CYS': 'C',  'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', # standard amino acids
                   'ILE': 'I',  'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', # standard amino acids
                   'GLY': 'G',  'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', # standard amino acids
                   'ALA': 'A',  'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', # standard amino acids
                   'PS3': 'S3', 'PS5': 'S5', 'PS8': 'S8', # stapled amino acids
                   'PR3': 'R3', 'PR5': 'R5', 'PR8': 'R8', # stapled amino acids
                   'ACE': 'Ac', 'NME': 'NH2',  # N/C terminal
                   'NLE': 'B', 'AIB': 'Aib', # non-standard amino acids
                   }

# 读取 CSV 文件
df = pd.read_csv('example/datasets/filter_stapep_amp.csv')
# 获取 SEQUENCE 列
sequences = df['SEQUENCE']

# 判断每个序列是否有不在氨基酸字典内的数据
def check_sequence_validity(sequence):
    # 列出所有在字典中的符号
    valid_symbols = ''.join(amino_acid_dict.values())

    # 判断序列中的每个字符是否在有效符号中
    for char in sequence:
        if char not in valid_symbols:
            print(f"Invalid character '{char}' found in sequence: {sequence}")
            return False  # 如果有无效字符则返回 False
    return True  # 如果所有字符都有效，返回 True
# 检查所有序列并打印无效序列
df['Valid'] = sequences.apply(check_sequence_validity)
# 输出验证结果
print(df[['SEQUENCE', 'Valid']])

