"""RFdiffusion 订书肽结合物生成流水线（由原 structure.py 的 __main__ 拆分而来）。

支持两种建模方式：
- denovo   : 解析 RFdiffusion 候选序列 -> ESMFold 从头建模 + OpenMM 短时 MD 优化
             -> 按 α 螺旋占比过滤 -> PyMOL align 对齐到靶点蛋白 A 链（align_denovo/）
             -> PyMOL cealign 二次对齐（cealign_denovo/，可用 --skip-cealign 跳过）
- modeller : 在 RFdiffusion 骨架上遍历插入 S5 订书位置 -> Modeller 同源建模 -> 对齐

用法示例：
    python -m stapep.rfdiffusion_denovo --method denovo --base-dir /data \
        --root-folder 2KOY/2KOY_100_4 --template-prefix 2koy --start 53 --end 60
"""
import argparse
import os
import sys
from pathlib import Path

# 支持 `python stapep/rfdiffusion_denovo.py` 直接运行时导入 stapep 包
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from Bio.PDB import PDBIO, PDBParser, Select

from stapep.calculate_alpha import calculate_alpha
from stapep.cealign_RFdiffusion_denove import cealign_proteins
from stapep.read_fa import parse_fa_files
from stapep.structure import AlignStructure, Structure

HELIX_RATIO_THRESHOLD = 0.3


class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id


def random_insert_stapep_peptide(root_path, pdb_path):
    """
    Args:
        root_path: 文件根路径
        pdb_path: RFdiffusion 生成的多肽骨架的 PDB 文件

    Returns: 依次添加 S5 订书肽后采用 Modeller 生成的订书肽骨架
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))

    # 获取 A 链数据单独保存为一个 pdb 文件
    root_path = Path(root_path)
    chain_A_path = root_path / "chain_A.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(chain_A_path), select=ChainSelect("A"))

    # 获取 A 链的残基数量
    chain_A = structure[0]['A']
    residues = [res for res in chain_A if res.id[0] == ' ']
    residue_count = len(residues)
    peptide_sequence = list('G' * residue_count)

    noalign_root_path = root_path / "noalign"
    align_root_path = root_path / "align"
    noalign_root_path.mkdir(parents=True, exist_ok=True)
    align_root_path.mkdir(parents=True, exist_ok=True)

    # 向 RFdiffusion 骨架中依次添加 S5 订书肽，然后采用 Modeller 建模
    for i in range(1, residue_count - 5):
        peptide_sequence_copy = peptide_sequence.copy()
        peptide_sequence_copy[i] = 'S5'
        peptide_sequence_copy[i + 4] = 'S5'
        seq = ''.join(peptide_sequence_copy)
        st = Structure(verbose=True)
        output_pdb_path = noalign_root_path / f"homology_model_{i}.pdb"
        st.generate_3d_structure_from_template(seq=seq,
                                               output_pdb=str(output_pdb_path),
                                               template_pdb=str(chain_A_path))
        align_output_pdb_path = align_root_path / f"aligned_homology_model_{i}.pdb"
        AlignStructure.align(ref_pdb=str(chain_A_path),
                             pdb=str(output_pdb_path),
                             output_pdb=str(align_output_pdb_path))


def _get_thresholds(same_threshold, res_threshold):
    """未显式指定阈值时使用默认值。"""
    same_threshold = 4 if same_threshold is None else same_threshold
    res_threshold = 6 if res_threshold is None else res_threshold
    return same_threshold, res_threshold


def run_denovo_pipeline(base_dir, root_folder, template_prefix, start, end,
                        same_threshold=None, res_threshold=None,
                        helix_threshold=HELIX_RATIO_THRESHOLD,
                        skip_cealign=False):
    """解析 RFdiffusion 候选序列 -> ESMFold 建模 + MD 优化 -> 螺旋过滤
    -> PyMOL align 对齐到靶点蛋白 -> cealign 二次对齐。"""
    for i in range(start, end):
        rfdiffusion_template = f"{template_prefix}_{i}"
        root_path = Path(base_dir) / root_folder / rfdiffusion_template
        same_threshold, res_threshold = _get_thresholds(same_threshold, res_threshold)

        fa_file_path = root_path / "fa"
        de_novo_folder_path = root_path / "denovo"
        de_novo_folder_path.mkdir(parents=True, exist_ok=True)

        _, index_list, seq_list = parse_fa_files(
            str(fa_file_path), same_threshold, res_threshold,
            rfdiffusion_template, template_prefix)
        for seq, index in zip(seq_list, index_list):
            st = Structure(verbose=True, save_tmp_dir=False)
            de_novo_path = de_novo_folder_path / f"{index}.pdb"
            if de_novo_path.exists():
                print(f"skip existing file: {de_novo_path}")
                continue
            print(f"NO file: {de_novo_path}")
            output_de_novo = st.de_novo_3d_structure(seq=seq, output_pdb=str(de_novo_path))
            if not output_de_novo:
                print(f"[ERROR] Iteration {index} - error")

        filenames = sorted(f for f in os.listdir(de_novo_folder_path) if f.endswith(".pdb"))
        align_filenames = []
        for filename in filenames:
            ligand_path = de_novo_folder_path / filename
            helix_ratio = calculate_alpha(str(ligand_path))
            if helix_ratio > helix_threshold:
                align_filenames.append(filename)

        protein_path = Path(base_dir) / root_folder / "RFdiffusion" / f"{rfdiffusion_template}.pdb"
        align_denovo_path = root_path / "align_denovo"
        align_denovo_path.mkdir(parents=True, exist_ok=True)
        for ligand_name in align_filenames:
            ligand_path = de_novo_folder_path / ligand_name
            ligand_stem = os.path.splitext(ligand_name)[0]
            align_output_pdb_path = align_denovo_path / f"{ligand_stem}.pdb"
            AlignStructure.align_denovo(ref_pdb=str(protein_path),
                                        pdb=str(ligand_path),
                                        output_pdb=str(align_output_pdb_path),
                                        rfdiffusion_template=rfdiffusion_template)

        if skip_cealign:
            continue

        # cealign 二次对齐（原独立 cealign 脚本的逻辑，现并入 denovo 流程）
        cealign_denovo_path = root_path / "cealign_denovo"
        cealign_denovo_path.mkdir(parents=True, exist_ok=True)
        for filename in sorted(f for f in os.listdir(align_denovo_path)
                               if f.endswith(".pdb")):
            try:
                cealign_proteins(str(protein_path),
                                 str(align_denovo_path / filename),
                                 str(cealign_denovo_path / filename))
            except Exception as e:
                print(f"[ERROR] cealign {filename}: {e}")


def run_modeller_pipeline(base_dir, root_folder, template_prefix, start, end):
    """遍历 RFdiffusion 骨架上的 S5 插入位置，用 Modeller 同源建模。"""
    for i in range(start, end):
        rfdiffusion_template = f"{template_prefix}_{i}"
        root_path = Path(base_dir) / root_folder / rfdiffusion_template
        root_path.mkdir(parents=True, exist_ok=True)
        pdb_path = Path(base_dir) / root_folder / "RFdiffusion" / f"{rfdiffusion_template}.pdb"
        random_insert_stapep_peptide(root_path, pdb_path)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="RFdiffusion 订书肽结合物流水线 "
                    "(denovo: ESMFold+MD+align+cealign; modeller: S5 遍历+同源建模)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--method", choices=["denovo", "modeller"], required=True,
                        help="建模方式（必填）")
    parser.add_argument("--base-dir", required=True,
                        help="数据根目录（必填）")
    parser.add_argument("--root-folder", required=True,
                        help="数据集子目录，相对 base-dir（必填）")
    parser.add_argument("--template-prefix", required=True,
                        help="模板名前缀，模板名 = 前缀_{i}（必填）")
    parser.add_argument("--start", type=int, required=True,
                        help="模板序号起始（含，必填）")
    parser.add_argument("--end", type=int, required=True,
                        help="模板序号结束（不含，必填）")
    parser.add_argument("--same-threshold", type=int, default=None,
                        help="连续相同残基数阈值（默认 4）")
    parser.add_argument("--res-threshold", type=int, default=None,
                        help="唯一残基种类数阈值（默认 6）")
    parser.add_argument("--helix-threshold", type=float, default=HELIX_RATIO_THRESHOLD,
                        help="denovo 流程保留 PDB 的 α 螺旋占比下限")
    parser.add_argument("--skip-cealign", action="store_true",
                        help="默认执行 cealign 二次对齐，加此参数则跳过")
    args = parser.parse_args(argv)

    if args.method == "denovo":
        run_denovo_pipeline(args.base_dir, args.root_folder, args.template_prefix,
                            args.start, args.end, args.same_threshold,
                            args.res_threshold, args.helix_threshold,
                            skip_cealign=args.skip_cealign)
    else:
        run_modeller_pipeline(args.base_dir, args.root_folder, args.template_prefix,
                              args.start, args.end)
    return 0


if __name__ == "__main__":
    sys.exit(main())