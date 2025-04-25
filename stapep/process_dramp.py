import csv

# 打开并读取dramp.csv文件
with open('example/datasets/dramp_data.csv', mode='r', encoding='utf-8') as file:
    reader = csv.reader(file)
    header = next(reader)  # 读取表头
    sequence_index = header.index("Hiden_Sequence")  # 找到"Hiden_Sequence"列的索引

    # 遍历每一行，判断是否包含两个X字母
    for row in reader:
        sequence = row[sequence_index]
        if sequence.count('X') % 2 != 0:
            print(sequence)
