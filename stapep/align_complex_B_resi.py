import os


def replace_mol_residues(equil_path, complex_path, output_path):
    # 1. 读取 equil_system.pdb，记录 MOL 残基行和位置
    equil_lines = []
    mol_residues = []  # 存储 MOL 残基号
    if not equil_path or not os.path.exists(equil_path):
        return ""
    with open(equil_path, 'r') as f:
        for line in f:
            equil_lines.append(line)
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                if resname in ["HOH", "Cl-", "Na+"]:
                    continue
                resi = int(line[22:26].strip())
                if resname == "MOL":
                    if resi not in mol_residues:
                        mol_residues.append(resi)

    print(f"Detected MOL residues in equil_system: {mol_residues}")

    # 2. 读取 complex_B_NCACO.pdb 的氨基酸序列
    complex_residues = []
    with open(complex_path, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                resi = int(line[22:26].strip())
                if not complex_residues or complex_residues[-1][0] != resi:
                    complex_residues.append((resi, resname))

    print(f"Detected {len(complex_residues)} residues in complex_B_NCACO")

    # 校验数量是否一致
    if len(mol_residues) != len(complex_residues):
        raise ValueError(f"Number of MOL residues ({len(mol_residues)}) does not match "
                         f"complex residues ({len(complex_residues)})")

    # 3. 构建 residue mapping：MOL 残基号 → 正确氨基酸名
    mapping = {mol_residues[i]: complex_residues[i][1] for i in range(len(mol_residues))}
    print("Residue mapping (equil → complex):")
    for k, v in mapping.items():
        print(f"  {k} → {v}")

    # 4. 替换 equil_system.pdb 中对应残基名
    new_lines = []
    max_ket_flag = False
    for line in equil_lines:
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            if resname in ["HOH", "Cl-", "Na+"]:
                continue
            if not max_ket_flag:
                resi = int(line[22:26].strip())
                max_key = max(mapping.keys())
                if resi > max_key:
                    max_ket_flag = True
                if resi in mapping:
                    # 替换残基名（列 17–20）
                    new_line = line[:17] + f"{mapping[resi]:>3}" + line[20:]
                    new_lines.append(new_line)
                else:
                    new_lines.append(line)
            else:
                new_lines.append(line)
        else:
            new_lines.append(line)

    # 5. 保存结果
    with open(output_path, 'w') as f:
        f.writelines(new_lines)
    print(f"✅ New file saved as: {output_path}")


if __name__ == '__main__':
    base_dir = "/home/d3008/Documents/zhr/STING/STING_100_1"
    for i in range(60):
        folder_name = f"sting_{i}"
        source_folder = os.path.join(base_dir, folder_name)
        cealign_folder = os.path.join(source_folder, "cealign_denovo")
        if os.path.isdir(cealign_folder) and os.listdir(cealign_folder):
            pdb_files = [f for f in os.listdir(cealign_folder) if f.endswith(".pdb")]
            for pdb_file in pdb_files:
                name = os.path.splitext(pdb_file)[0]
                equil_path = os.path.join(source_folder, "openbpmd", name, "equil_system.pdb")
                complex_path = os.path.join(source_folder, "openbpmd", name, "complex_B_NCACO.pdb")
                output_path = os.path.join(source_folder, "openbpmd", name, "equil_system_corrected.pdb")
                replace_mol_residues(equil_path, complex_path, output_path)
