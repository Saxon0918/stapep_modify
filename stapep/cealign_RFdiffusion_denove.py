import os
from pymol import cmd

def cealign_proteins(protein_path, ligand_path, align_output_pdb_path):
    """
    使用 PyMOL 的 cealign 将 ligand 结构对齐到参考 protein 的 A 链，并保存结果。

    Parameters
    ----------
    protein_path : str
        参考结构的 PDB 文件路径（用于提取 A 链）
    ligand_path : str
        待对齐的结构文件路径
    align_output_pdb_path : str
        输出对齐后结构的保存路径（PDB 文件）
    """

    # 1. 清空 PyMOL 当前环境
    cmd.reinitialize()
    cmd.load(protein_path, "ref_protein")
    cmd.load(ligand_path, "ligand")
    cmd.select("ref_chainA", "ref_protein and chain A")
    cmd.cealign("ref_chainA", "ligand")
    cmd.save(align_output_pdb_path, "ligand")
    cmd.delete("ref_chainA")
    print(f"[INFO] Alignment completed. Output saved to: {align_output_pdb_path}")


# 示例调用（请修改为你自己的路径）
if __name__ == "__main__":
    root_folder = 'STING/STING_100_1'
    for i in range(100):
        rfdiffusion_template = f'sting_{i}'
        align_output_folder = f'/home/d3008/Documents/zhr/{root_folder}/{rfdiffusion_template}/cealign_denovo'
        if not os.path.exists(align_output_folder):
            os.makedirs(align_output_folder)

        root_path = f'/home/d3008/Documents/zhr/{root_folder}/{rfdiffusion_template}/denovo'
        # 遍历 test 文件夹下所有 .pdb 文件
        for filename in os.listdir(root_path):
            if filename.endswith(".pdb"):
                try:
                    ligand_path = os.path.join(root_path, filename)
                    name = os.path.splitext(filename)[0]
                    rfdiffusion_name = "_".join(name.split("_")[:2])
                    protein_path = f'/home/d3008/Documents/zhr/{root_folder}/RFdiffusion/{rfdiffusion_name}.pdb'
                    align_output_pdb_path = f'{align_output_folder}/{filename}'
                    cealign_proteins(protein_path, ligand_path, align_output_pdb_path)
                except Exception as e:
                    print(f"error {root_path}")

    cmd.quit()

