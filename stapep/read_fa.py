import os
import re


# =========================
# 靶点筛选规则配置
# =========================

# 各靶点对生成序列的氨基酸硬性要求
# mode: "ALL" 表示必须全部包含, "ANY" 表示至少包含一个
TEMPLATE_RESIDUE_RULES = {
    "2gv2": ({"F", "W", "L"}, "ALL"),
    "pdl1": ({"F", "W"}, "ALL"),
    "clpp": ({"F", "W", "Y"}, "ANY"),
    "sting": ({"Y"}, "ALL"),
}

# 1gng 靶点对 G 数量的特殊限制
MAX_G_COUNT_1GNG = 3
# 其它靶点要求 G 占比低于 1/4
MAX_G_RATIO = 1 / 4


def parse_fa_files(folder_path, same_threshold, res_threshold, rfdiffusion_template,
                   template_prefix):
    """解析 RFdiffusion 输出的 .fa 文件并筛选候选序列。

    template_prefix: 模板名前缀（对应 CLI 的 --template-prefix），
                     用于把 .fa 文件名中的 "{prefix}_B" 替换为实际模板名。
    """
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

            # 各靶点要求的氨基酸组成检查
            residues_required, mode = _match_required_residues(rfdiffusion_template)
            if residues_required is not None:
                if mode == "ALL" and not residues_required.issubset(tokens):
                    continue
                if mode == "ANY" and residues_required.isdisjoint(tokens):
                    continue

            # G 数量检查: 1gng 最多 3 个 G; 所有靶点 G 占比需低于 1/4
            g_token_count = tokens.count('G')
            if g_token_count >= len(tokens) * MAX_G_RATIO:
                continue
            if rfdiffusion_template.startswith("1gng") and g_token_count >= MAX_G_COUNT_1GNG:
                continue

            # check S5 gap
            S5_flag = detect_S5_distance(tokens)
            if not S5_flag:
                continue

            # 连续相同残基不能超过 same_threshold
            count = 1
            repeat_flag = False
            for j in range(1, len(tokens)):
                if tokens[j] == tokens[j-1]:
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
                name = replace_template_name(name, rfdiffusion_template, template_prefix)
                names.append(name)
                seqs.append(seq_line)

    names, seqs = remove_duplicate_seqs(names, seqs)

    return folder_name, names, seqs


def _match_required_residues(rfdiffusion_template):
    """返回 (required_residues, mode); 若模板不在配置中则返回 (None, None)"""
    for prefix, (residues, mode) in TEMPLATE_RESIDUE_RULES.items():
        if rfdiffusion_template.startswith(prefix):
            return residues, mode
    return None, None


def replace_template_name(name, rfdiffusion_template, template_prefix):
    """把 name 中的 "{template_prefix}_B" 替换为实际模板名 rfdiffusion_template。"""
    if rfdiffusion_template.startswith(template_prefix):
        name = name.replace(f"{template_prefix}_B", rfdiffusion_template)
    return name


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


def detect_S5_distance(tokens):
    tokens_stapled = [i for i in tokens if re.search(r'\d', i)]
    if not all(i == 'S5' for i in tokens_stapled):
        return True
    S5_index = [i for i, token in enumerate(tokens) if token == "S5"]
    if len(S5_index) % 2 != 0:
        return False
    S5_flag = True
    for i in range(0, len(S5_index), 2):
        if S5_index[i + 1] - S5_index[i] != 4:
            S5_flag = False
            break
    return S5_flag