"""使用 PyMOL cealign 将 denovo 肽段对齐到 RFdiffusion 模板蛋白 A 链（二次对齐）。

用法示例：
    python -m stapep.cealign_RFdiffusion_denove --base-dir /data \
        --root-folder 2KOY/2KOY_100_4 --template-prefix 2koy --start 53 --end 100
"""
import argparse
import os
from pathlib import Path

from pymol import cmd

DEFAULT_BASE_DIR = "/home/d3008/Documents"
DEFAULT_ROOT_FOLDER = "2KOY/2KOY_100_4"
DEFAULT_TEMPLATE_PREFIX = "2koy"
DEFAULT_START = 53
DEFAULT_END = 100


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


def run_cealign(base_dir, root_folder, template_prefix, start, end):
    """对 [start, end) 区间内每个模板的 align_denovo 结构执行 cealign 二次对齐。"""
    for i in range(start, end):
        rfdiffusion_template = f"{template_prefix}_{i}"
        template_dir = Path(base_dir) / root_folder / rfdiffusion_template
        align_output_folder = template_dir / "cealign_denovo"
        align_output_folder.mkdir(parents=True, exist_ok=True)

        ligand_root_path = template_dir / "align_denovo"
        protein_path = Path(base_dir) / root_folder / "RFdiffusion" / f"{rfdiffusion_template}.pdb"

        if not ligand_root_path.exists():
            print(f"[WARN] {ligand_root_path} 不存在，跳过")
            continue

        # 遍历 align_denovo 文件夹下所有 .pdb 文件
        for filename in sorted(f for f in os.listdir(ligand_root_path) if f.endswith(".pdb")):
            try:
                ligand_path = ligand_root_path / filename
                align_output_pdb_path = align_output_folder / filename
                cealign_proteins(str(protein_path), str(ligand_path), str(align_output_pdb_path))
            except Exception as e:
                print(f"[ERROR] {filename}: {e}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="用 PyMOL cealign 将 align_denovo 肽段二次对齐到模板蛋白 A 链",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--base-dir", default=DEFAULT_BASE_DIR, help="数据根目录")
    parser.add_argument("--root-folder", default=DEFAULT_ROOT_FOLDER, help="数据集子目录（相对 base-dir）")
    parser.add_argument("--template-prefix", default=DEFAULT_TEMPLATE_PREFIX,
                        help="模板名前缀，模板名 = 前缀_{i}")
    parser.add_argument("--start", type=int, default=DEFAULT_START, help="模板序号起始（含）")
    parser.add_argument("--end", type=int, default=DEFAULT_END, help="模板序号结束（不含）")
    args = parser.parse_args(argv)

    try:
        run_cealign(args.base_dir, args.root_folder, args.template_prefix,
                    args.start, args.end)
    finally:
        cmd.quit()
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())