import os
import csv
import re

def parse_fa_files(folder_path, same_threshold, res_threshold, rfdiffusion_template):
    result = []
    names = []
    seqs = []
    global_scores = []
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

            global_score_match = re.search(r'global_score=((?:0|[1-9]\d*)(?:\.\d+)?)', score_line)
            global_score = float(global_score_match.group(1)) if global_score_match else 0

            if reward <= 0:
                continue

            tokens = re.findall(r'[SR]\d+|[A-Z]', seq_line)

            # for 2gv2 check F, W, L is in sequence or not
            required_2gv2 = {"F", "W", "L"}
            if not required_2gv2.issubset(tokens) and rfdiffusion_template.startswith("2gv2"):
                continue

            required_pdl1 = {"F", "W"}
            if not required_pdl1.issubset(tokens) and rfdiffusion_template.startswith("pdl1"):
                continue

            required_clpp = {"F", "W", "Y"}
            if required_clpp.isdisjoint(tokens) and rfdiffusion_template.startswith("clpp"):
                continue

            required_sting = {"Y"}
            if not required_sting.issubset(tokens) and rfdiffusion_template.startswith("sting"):
                continue

            # check G number, for 1gng must minimize than 1/4 of whole list, for 2gv2 must minimize than 1/3 of whole list
            G_token_count = tokens.count('G')
            if G_token_count >= 3 and rfdiffusion_template.startswith("1gng"):
                continue
            elif G_token_count >= len(tokens)/4:
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
            else:
                if len(tokens) <= 10:
                    final_res_threshold = len(tokens) - 2
                elif len(tokens) <= 12:
                    final_res_threshold = len(tokens) - 3
                else:
                    final_res_threshold = min(len(tokens) - 4, 10)

            if len(unique_aas) >= final_res_threshold:  # residues classes threshold
                name = f"{os.path.splitext(filename)[0]}_{sample}"
                if rfdiffusion_template.startswith("1gng"):
                    name = name.replace("1gng_B", rfdiffusion_template)
                elif rfdiffusion_template.startswith("2gv2"):
                    name = name.replace("2gv2_B", rfdiffusion_template)
                elif rfdiffusion_template.startswith("pdl1"):
                    name = name.replace("pdl1_B", rfdiffusion_template)
                elif rfdiffusion_template.startswith("clpp"):
                    name = name.replace("clpp_B", rfdiffusion_template)
                elif rfdiffusion_template.startswith("sting"):
                    name = name.replace("sting_B", rfdiffusion_template)
                elif rfdiffusion_template.startswith("ipad"):
                    name = name.replace("ipad_B", rfdiffusion_template)
                elif rfdiffusion_template.startswith("1ou8"):
                    name = name.replace("1ou8_B", rfdiffusion_template)
                elif rfdiffusion_template.startswith("7aa4"):
                    name = name.replace("7aa4_B", rfdiffusion_template)
                elif rfdiffusion_template.startswith("8e4a"):
                    name = name.replace("8e4a_B", rfdiffusion_template)
                elif rfdiffusion_template.startswith("6ih0"):
                    name = name.replace("6ih0_B", rfdiffusion_template)
                elif rfdiffusion_template.startswith("9co2"):
                    name = name.replace("9co2_B", rfdiffusion_template)
                names.append(name)
                seqs.append(seq_line)
                global_scores.append(global_score)
                result.append([name, seq_line])
    names, seqs, global_scores = remove_duplicate_seqs(names, seqs, global_scores)

    # output_csv = "/home/d3008/Documents/zhr/9CO2/global_scores_result.csv"
    # assert len(names) == len(global_scores) == len(seqs), (
    #     f"Length mismatch: "
    #     f"names={len(names)}, "
    #     f"global_scores={len(global_scores)}, "
    #     f"seqs={len(seqs)}"
    # )
    # file_exists = os.path.exists(output_csv)
    # with open(output_csv, mode="a", newline="", encoding="utf-8") as f:
    #     writer = csv.writer(f)
    #     if not file_exists:
    #         writer.writerow(["name", "global_score", "sequence"])
    #     for name, score, seq in zip(names, global_scores, seqs):
    #         writer.writerow([name, score, seq])

    return folder_name, names, seqs


def remove_duplicate_seqs(names, seqs, global_scores):
    seen = set()
    unique_seqs = []
    unique_names = []
    unique_global_scores = []

    for seq, name, global_score in zip(seqs, names, global_scores):
        if seq not in seen:
            seen.add(seq)
            unique_seqs.append(seq)
            unique_names.append(name)
            unique_global_scores.append(global_score)

    return unique_names, unique_seqs, unique_global_scores


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


