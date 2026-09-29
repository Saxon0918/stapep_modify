"""批量从数据集 CSV 生成订书肽 PDB（由原 structure.py 的 generate_pdb 拆分而来）。

用法示例：
    python -m stapep.batch_generate_dataset_pdb --csv example/datasets/0001-3906_sup.csv \
        --out-dir example/stapep_data/pdb_sup --start 101 --end 200
"""
import argparse
import csv
import gc
import os
import sys
from pathlib import Path

os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'
import torch  # noqa: E402

# 支持 `python stapep/batch_generate_dataset_pdb.py` 直接运行时导入 stapep 包
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from stapep.structure import Structure  # noqa: E402

DEFAULT_CSV = "example/datasets/0001-3906_sup.csv"
DEFAULT_OUT_DIR = "example/stapep_data/pdb_sup"
DEFAULT_START = 101
DEFAULT_END = 19993


def generate_pdb(st, csv_file=DEFAULT_CSV, out_dir=DEFAULT_OUT_DIR,
                 start=DEFAULT_START, end=DEFAULT_END):
    """
    Args:
        st: 传入的 Structure(verbose=True) 实例

    Returns: 采用 CSV 数据集中 Sequence 生成的订书肽 PDB 文件
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sequence_column = []
    index_column = []

    with open(csv_file, mode='r', encoding='utf-8') as file:
        csv_reader = csv.reader(file)
        headers = next(csv_reader)
        sequence_index = headers.index('seq')
        stapep_id = headers.index('stapep_id')
        for row in csv_reader:
            sequence_column.append(row[sequence_index])
            index_column.append(row[stapep_id])

    end = min(end, len(sequence_column))  # 防止序号越界
    for i in range(start, end):
        seq = sequence_column[i]
        torch.cuda.empty_cache()
        gc.collect()
        de_novo_path = out_dir / f"stapep_{index_column[i]}.pdb"
        output_de_novo = st.de_novo_3d_structure(seq=seq, output_pdb=str(de_novo_path))
        if not output_de_novo:
            record = f"Iteration {i} - error"
            print(record)
    print("FINISH")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="批量从数据集 CSV 生成订书肽 PDB",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--csv", default=DEFAULT_CSV, help="数据集 CSV 路径（需含 seq、stapep_id 列）")
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR, help="PDB 输出目录")
    parser.add_argument("--start", type=int, default=DEFAULT_START, help="起始行号（含）")
    parser.add_argument("--end", type=int, default=DEFAULT_END, help="结束行号（不含，超长自动截断）")
    parser.add_argument("--quiet", action="store_true", help="关闭 Structure 详细日志")
    args = parser.parse_args(argv)

    st = Structure(verbose=not args.quiet)
    generate_pdb(st, csv_file=args.csv, out_dir=args.out_dir,
                 start=args.start, end=args.end)
    return 0


if __name__ == "__main__":
    sys.exit(main())