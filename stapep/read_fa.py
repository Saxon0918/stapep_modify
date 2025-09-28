import os
import csv
import re

def parse_fa_files(folder_path, same_threshold, res_threshold, rfdiffusion_template):
    result = []
    names = []
    seqs = []
    folder_name = os.path.basename(os.path.normpath(folder_path))
    file_path = f"{folder_path}/epoch1000_step1000/seqs"
    filenames = sorted(f for f in os.listdir(file_path) if f.endswith(".fa"))
    for filename in filenames:
        filepath = os.path.join(file_path, filename)
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # 从第3行开始（索引2），每两行为一组：奇数行为信息，偶数行为序列
        for i in range(2, len(lines), 2):
            score_line = lines[i].strip()
            seq_line = lines[i+1].strip() if i+1 < len(lines) else ''

            reward_match = re.search(r'reward=([0-9.]+)', score_line)
            reward = float(reward_match.group(1)) if reward_match else 0

            sample_match = re.search(r'sample=(\d+)', score_line)
            sample = sample_match.group(1) if sample_match else 'unknown'

            if reward <= 0:
                continue

            tokens = re.findall(r'[SR]\d+|[A-Z]', seq_line)

            # for 2gv2 check F, W, L is in sequence or not
            required_2gv2 = {"F", "W", "L"}
            if not required_2gv2.issubset(tokens) and rfdiffusion_template.startswith("2gv2"):
                continue

            # check G number, for 1gng must minimize than 1/4 of whole list, for 2gv2 must minimize than 1/3 of whole list
            G_token_count = tokens.count('G')
            if G_token_count >= 3 and rfdiffusion_template.startswith("1gng"):
                continue
            elif G_token_count >= len(tokens)/4 and rfdiffusion_template.startswith("2gv2"):
                continue

            # check S5 gap
            S5_flag = detect_S5_distance(tokens, filepath)
            if not S5_flag:
                continue

            count = 1
            repeat_flag = False
            for i in range(1,len(tokens)):
                if tokens[i] == tokens[i-1]:
                    count += 1
                    if count >= same_threshold:  # same residues sequence length threshold
                        repeat_flag = True
                        break
                else:
                    count = 1
            if repeat_flag:
                continue
            unique_aas = set(tokens)

            # sequence length threshold
            if rfdiffusion_template.startswith("1gng"):  # for 1gng
                final_res_threshold = res_threshold
            elif rfdiffusion_template.startswith("2gv2"):
                if len(tokens) <= 10:
                    final_res_threshold = len(tokens) - 2
                elif len(tokens) <= 12:
                    final_res_threshold = len(tokens) - 3
                else:
                    final_res_threshold = min(len(tokens) - 4, 10)
            else:
                final_res_threshold = 1

            if len(unique_aas) >= final_res_threshold:  # residues classes threshold
                name = f"{os.path.splitext(filename)[0]}_{sample}"
                if rfdiffusion_template.startswith("1gng"):
                    name = name.replace("1gng_B", rfdiffusion_template)
                elif rfdiffusion_template.startswith("2gv2"):
                    name = name.replace("2gv2_B", rfdiffusion_template)
                names.append(name)
                seqs.append(seq_line)
                result.append([name, seq_line])
    names, seqs = remove_duplicate_seqs(names, seqs)
    return folder_name, names, seqs


def remove_duplicate_seqs(names, seqs):
    seen = set()
    unique_seqs = []
    unique_names = []

    for seq, name in zip(seqs, names):
        if seq not in seen:
            seen.add(seq)
            unique_seqs.append(seq)
            unique_names.append(name)

    return unique_names, unique_seqs


def detect_S5_distance(tokens, filepath):
    tokens_stapled = [i for i in tokens if re.search(r'\d', i)]
    if not all(i == 'S5' for i in tokens_stapled):
        return True
    S5_index = [i for i, token in enumerate(tokens) if token == "S5"]
    S5_flag = True
    for i in range(0, len(S5_index), 2):
        if S5_index[i + 1] - S5_index[i] != 4:
            S5_flag = False
            break
    return S5_flag


