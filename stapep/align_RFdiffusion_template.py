import os
from pymol import cmd

# 设置参考结构路径和目标文件夹
ref_pdb = "/home/d3008/Documents/zhr/MDM2/2gv2.pdb"
pdb_folder = "/home/d3008/Documents/zhr/pycharmprojects/RFdiffusion/examples/example_outputs/2gv2_500_4/"
output_folder = "/home/d3008/Documents/zhr/MDM2_500_4/RFdiffusion/"

# 加载参考结构
cmd.reinitialize()
cmd.load(ref_pdb, "ref")

# 遍历 test 文件夹下所有 .pdb 文件
for filename in os.listdir(pdb_folder):
    if filename.endswith(".pdb"):
        filepath = os.path.join(pdb_folder, filename)
        name = os.path.splitext(filename)[0]

        # 加载并align
        cmd.load(filepath, name)
        cmd.align(name, "ref")

        # 保存对齐后的结构
        aligned_path = os.path.join(output_folder, f"{name}.pdb")
        cmd.save(aligned_path, name)
        print(f"Aligned and saved: {aligned_path}")

cmd.quit()

