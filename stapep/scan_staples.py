"""扫描不同订书类型/插入位置对肽段理化性质的影响（由原 structure.py 的
calculate_all_kinds_stapep 拆分而来）。

用法示例：
    python -m stapep.scan_staples --seq ACDEFGHIKLMNPQRSTVWY \
        --out-file example/insulin/insulin.txt --pdb-dir example/insulin
"""
import argparse
import os
import sys
from pathlib import Path

# 支持 `python stapep/scan_staples.py` 直接运行时导入 stapep 包
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from stapep.structure import Structure
from stapep.utils import PhysicochemicalPredictor

STAPEP_WAYS = ['S5_S5_3', 'R8_S5_6', 'R5_S8_6', 'R5_S5_2', 'S5_S5_4',
               'R8_S5_5', 'R8_S5_3', 'S5_R8_6', 'R5_S8_7']
COLUMNS = [
    "Iteration", "seq", "helix percent", "sheet percent", "loop percent",
    "mean bfactor", "mol surf", "mean gyrate", "psa",
    "total number of hydrogen bonds", "smiles"
]
DEFAULT_OUT_FILE = "example/insulin/insulin.txt"
DEFAULT_PDB_DIR = "example/insulin"


def calculate_all_kinds_stapep(seq, st, out_file=DEFAULT_OUT_FILE, pdb_dir=DEFAULT_PDB_DIR):
    """遍历订书类型（首尾残基及间隔），对每个插入位置生成结构并计算理化性质。

    Args:
        seq: 待扫描的多肽序列
        st: Structure 实例（建议 save_tmp_dir=True，缓存 MD 轨迹供特征计算）
        out_file: 结果记录文件路径
        pdb_dir: denovo 结构与平均结构 PDB 的输出目录
    """
    out_file = Path(out_file)
    pdb_dir = Path(pdb_dir)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    pdb_dir.mkdir(parents=True, exist_ok=True)

    if not out_file.exists():
        with open(out_file, 'w') as file:
            file.write('\t'.join(COLUMNS) + '\n')  # 写入列名

    # 使用 for 循环遍历 seq，生成每一轮的 new_seq
    for stapep_way in STAPEP_WAYS:
        stapep_first, stapep_second, interval = stapep_way.split('_')
        interval = int(interval)

        for i in range(0, len(seq)):
            new_seq = seq[:i]
            # 插入 stapep_first 和 stapep_second 后，跳过 interval 个字符
            if i == 0:
                new_seq = new_seq + stapep_first + '-' + seq[i:i + interval]
            else:
                new_seq = new_seq + '-' + stapep_first + '-' + seq[i:i + interval]

            back_i = i + interval
            if back_i < len(seq):  # 确保不会越界
                new_seq = new_seq + '-' + stapep_second + '-' + seq[back_i:]
            elif back_i == len(seq):
                new_seq = new_seq + '-' + stapep_second
            else:
                break

            # 去模板预测，采用 ESMFold
            de_novo_path = pdb_dir / f"denovo_{stapep_first}_{stapep_second}_{interval}_{i}.pdb"
            output_de_novo = st.de_novo_3d_structure(seq=new_seq, output_pdb=str(de_novo_path))
            if not output_de_novo:
                record = f"Iteration {stapep_first}_{stapep_second}_{interval}_{i} - error"
                with open(out_file, 'a') as file:
                    file.write(record + '\n')
                continue

            pathname = st.tmp_dir  # 拓扑文件和轨迹文件所在路径
            pcp = PhysicochemicalPredictor(sequence=new_seq,
                                           topology_file=os.path.join(pathname, 'pep_vac.prmtop'),
                                           trajectory_file=os.path.join(pathname, 'traj.dcd'),
                                           start_frame=0)

            # Get the features
            helix_percent = pcp.calc_helix_percent()
            sheet_percent = pcp.calc_extend_percent()
            loop_percent = pcp.calc_loop_percent()

            print('helix percent: ', helix_percent)  # 计算α螺旋占比
            print('sheet percent: ', sheet_percent)  # 计算β折叠占比
            print('loop percent: ', loop_percent)  # 计算环结构（loop）占比
            # save the mean structure of the trajectory
            mean_structure_path = pdb_dir / f"mean_structure_{stapep_first}_{stapep_second}_{interval}_{i}.pdb"
            pcp._save_mean_structure(str(mean_structure_path))  # 保存分子的平均结构的 pdb

            # calculate the Mean B-factor, Molecular Surface, Mean Gyration Radius and 3D-PSA
            mean_bfactor = pcp.calc_mean_bfactor()
            mol_surf = pcp.calc_mean_molsurf()
            mean_gyrate = pcp.calc_mean_gyrate()
            psa = pcp.calc_psa(str(mean_structure_path))
            total_number_hydrogen_bonds = pcp.calc_n_hbonds()
            print('mean bfactor: ', mean_bfactor)  # 平均B因子
            print('mol surf: ', mol_surf)  # 分子表面积
            print('mean gyrate: ', mean_gyrate)  # 旋转半径
            print('psa: ', psa)  # 极性表面积
            print('total number of hydrogen bonds: ', total_number_hydrogen_bonds)  # 氢键总数

            # extract 2D structure of the peptide
            smiles = pcp.extract_2d_structure(str(mean_structure_path))
            print(smiles)
            # 构建记录内容
            record = f"{stapep_first}_{stapep_second}_{interval}_{i}\t{new_seq}\t{helix_percent}\t{sheet_percent}\t{loop_percent}\t" \
                     f"{mean_bfactor}\t{mol_surf}\t{mean_gyrate}\t{psa}\t" \
                     f"{total_number_hydrogen_bonds}\t{smiles}"

            # 追加记录到文件
            with open(out_file, 'a') as file:
                file.write(record + '\n')


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="扫描不同订书类型/插入位置，生成结构并计算理化性质",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--seq", required=True, help="待扫描的多肽序列")
    parser.add_argument("--out-file", default=DEFAULT_OUT_FILE, help="结果记录文件路径")
    parser.add_argument("--pdb-dir", default=DEFAULT_PDB_DIR, help="PDB 输出目录")
    args = parser.parse_args(argv)

    st = Structure(verbose=True, save_tmp_dir=True)
    calculate_all_kinds_stapep(args.seq, st, out_file=args.out_file, pdb_dir=args.pdb_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())