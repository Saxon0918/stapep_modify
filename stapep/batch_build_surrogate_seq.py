"""批量构建 surrogate 序列结构（为 AF3 准备输入），由原 batch_build_surrogate_seq.py 整理而来。

把 cealign_denovo 中的订书肽（PS5/S5）替换为代理残基 GLY（仅保留主链原子），
输出 surrogate PDB + mmCIF + meta JSON，并生成 AF3 批量运行所需的 sequences.csv。

用法示例：
    python -m stapep.batch_build_surrogate_seq --base-folder 2KOY --protein-name 2koy

    # 完全自定义路径：
    python -m stapep.batch_build_surrogate_seq \
        --source-root /data/2KOY/2KOY_100_4 \
        --target-root /data/2KOY/2KOY_100_4_cif \
        --protein-pdb /data/2KOY/2KOY_100_4/RFdiffusion/2koy_0.pdb \
        --protein-chain B
"""
import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Dict, List

from Bio.PDB import PDBParser
from Bio.PDB import MMCIFIO
from Bio.SeqUtils import seq1

# =========================
# 默认路径配置（与旧实现一致）
# =========================
DEFAULT_BASE_DIR = "/home/d3008/Documents"
DEFAULT_BASE_FOLDER = "2KOY"
DEFAULT_PROTEIN_NAME = "2koy"
SOURCE_FOLDER_SUFFIX = "_100_4"
TARGET_FOLDER_SUFFIX = "_100_4_cif"
DEFAULT_REVISION_DATE = "2024-01-01"
DEFAULT_PROTEIN_CHAIN_ID = "B"

# 当前方案：PS5 -> GLY
DEFAULT_STAPLED_RESNAMES = {"PS5", "S5"}
DEFAULT_SURROGATE_RESNAME = "GLY"

STANDARD_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIE": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# 代理残基只保留这些主链原子
SURROGATE_ALLOWED_ATOMS = {
    "GLY": {"N", "CA", "C", "O", "OXT"},
}

# manifest 记录字段
MANIFEST_FIELDS = [
    "source_pdb",
    "surrogate_pdb",
    "surrogate_cif",
    "meta_json",
    "status",
    "message",
    "length",
    "sequence",
    "converted_stapled_count",
]


# =========================
# 通用工具
# =========================
def clean_resname(resname: str) -> str:
    return resname.strip().upper()


def is_hydrogen(atom) -> bool:
    element = (atom.element or "").strip().upper()
    name = atom.get_name().strip().upper()
    return element == "H" or name.startswith("H")


def is_supported_residue(residue, stapled_resnames) -> bool:
    resname = clean_resname(residue.get_resname())
    return (resname in STANDARD_3TO1) or (resname in stapled_resnames)


def get_first_chain(structure):
    model = next(structure.get_models())
    chains = list(model.get_chains())
    if not chains:
        return None
    return chains[0]


def format_atom_name_for_pdb(atom_name: str) -> str:
    atom_name = atom_name.strip()
    if len(atom_name) >= 4:
        return atom_name[:4]
    return atom_name.rjust(4)


def pdb_atom_line(
    serial: int,
    atom_name: str,
    resname: str,
    chain_id: str,
    resseq: int,
    x: float,
    y: float,
    z: float,
    occupancy: float,
    bfactor: float,
    element: str,
) -> str:
    atom_name_fmt = format_atom_name_for_pdb(atom_name)
    element = (element or atom_name[0]).strip().upper()[:2]
    return (
        f"{'ATOM':<6}{serial:>5} {atom_name_fmt:<4} "
        f"{resname:>3} {chain_id:1}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"{occupancy:>6.2f}{bfactor:>6.2f}"
        f"          {element:>2}\n"
    )


def _manifest_row(source_pdb: str = "", status: str = "skip", message: str = "") -> Dict:
    row = {field: "" for field in MANIFEST_FIELDS}
    row["source_pdb"] = source_pdb
    row["status"] = status
    row["message"] = message
    return row


def extract_protein_sequence_from_pdb(pdb_path: Path, chain_id: str) -> str:
    """
    从固定蛋白 PDB 中提取指定链的氨基酸序列。
    只读取标准残基（hetflag == " "）。
    非标准残基如果 Biopython 无法识别，则记为 X。
    """
    if not pdb_path.exists():
        raise FileNotFoundError(f"蛋白 PDB 不存在: {pdb_path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))

    model = next(structure.get_models())

    target_chain = None
    for chain in model:
        if str(chain.id).strip() == str(chain_id).strip():
            target_chain = chain
            break

    if target_chain is None:
        raise ValueError(f"在 {pdb_path} 中未找到链 {chain_id}")

    seq_chars = []
    for residue in target_chain:
        hetflag = residue.id[0]
        if hetflag != " ":
            continue

        resname = clean_resname(residue.get_resname())
        try:
            seq_chars.append(seq1(resname))
        except Exception:
            seq_chars.append("X")

    sequence = "".join(seq_chars)
    if not sequence:
        raise ValueError(f"从 {pdb_path} 的链 {chain_id} 未提取到有效序列")

    return sequence


# =========================
# surrogate 结构生成
# =========================
def build_surrogate_records(chain, stapled_resnames, surrogate_resname) -> (List[Dict], Dict):
    records = []
    residue_meta = []
    sequence_chars = []

    atom_serial = 1
    new_resseq = 1
    stapled_count = 0

    for residue in chain:
        if not is_supported_residue(residue, stapled_resnames):
            continue

        original_resname = clean_resname(residue.get_resname())
        is_stapled = original_resname in stapled_resnames

        if is_stapled:
            surrogate_resname_for_res = surrogate_resname
            allowed_atoms = SURROGATE_ALLOWED_ATOMS[surrogate_resname]
            stapled_count += 1
        else:
            surrogate_resname_for_res = original_resname
            allowed_atoms = None

        residue_atoms = []
        seen_atom_names = set()

        for atom in residue.get_atoms():
            atom_name = atom.get_name().strip().upper()

            if is_hydrogen(atom):
                continue

            if allowed_atoms is not None and atom_name not in allowed_atoms:
                continue

            if atom_name in seen_atom_names:
                continue
            seen_atom_names.add(atom_name)

            coord = atom.get_coord()
            occupancy = 1.0 if atom.get_occupancy() is None else float(atom.get_occupancy())
            bfactor = 0.0 if atom.get_bfactor() is None else float(atom.get_bfactor())
            element = (atom.element or atom_name[0]).strip().upper()

            residue_atoms.append({
                "serial": atom_serial,
                "atom_name": atom_name,
                "resname": surrogate_resname_for_res,
                "chain_id": "A",
                "resseq": new_resseq,
                "x": float(coord[0]),
                "y": float(coord[1]),
                "z": float(coord[2]),
                "occupancy": occupancy,
                "bfactor": bfactor,
                "element": element,
            })
            atom_serial += 1

        atom_names = {a["atom_name"] for a in residue_atoms}
        if not {"N", "CA", "C", "O"}.issubset(atom_names):
            continue

        records.extend(residue_atoms)
        sequence_chars.append(STANDARD_3TO1[surrogate_resname_for_res])

        hetflag, orig_resseq, orig_icode = residue.id
        residue_meta.append({
            "new_chain_id": "A",
            "new_resseq": new_resseq,
            "source_chain_id": str(chain.id).strip(),
            "source_hetflag": str(hetflag).strip(),
            "source_resseq": int(orig_resseq),
            "source_icode": str(orig_icode).strip(),
            "source_resname": original_resname,
            "surrogate_resname": surrogate_resname_for_res,
            "is_stapled_source": is_stapled,
        })

        new_resseq += 1

    meta = {
        "source_chain_id": str(chain.id).strip(),
        "surrogate_chain_id": "A",
        "surrogate_resname_for_stapled": surrogate_resname,
        "stapled_source_resnames": sorted(list(stapled_resnames)),
        "sequence": "".join(sequence_chars),
        "length": len(sequence_chars),
        "n_stapled_positions_converted": stapled_count,
        "residue_map": residue_meta,
    }

    return records, meta


def write_surrogate_pdb(records: List[Dict], out_pdb: Path):
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    with open(out_pdb, "w", encoding="utf-8") as f:
        for rec in records:
            f.write(
                pdb_atom_line(
                    serial=rec["serial"],
                    atom_name=rec["atom_name"],
                    resname=rec["resname"],
                    chain_id=rec["chain_id"],
                    resseq=rec["resseq"],
                    x=rec["x"],
                    y=rec["y"],
                    z=rec["z"],
                    occupancy=rec["occupancy"],
                    bfactor=rec["bfactor"],
                    element=rec["element"],
                )
            )
        f.write("TER\nEND\n")


def ensure_release_date_in_cif(cif_path: Path, revision_date: str = DEFAULT_REVISION_DATE):
    content = cif_path.read_text(encoding="utf-8")
    if "_pdbx_audit_revision_history.revision_date" in content:
        return

    with open(cif_path, "a", encoding="utf-8") as f:
        f.write("\n#\n")
        f.write("loop_\n")
        f.write("_pdbx_audit_revision_history.ordinal\n")
        f.write("_pdbx_audit_revision_history.data_content_type\n")
        f.write("_pdbx_audit_revision_history.major_revision\n")
        f.write("_pdbx_audit_revision_history.minor_revision\n")
        f.write("_pdbx_audit_revision_history.revision_date\n")
        f.write(f"1 'Structure model' 1 0 {revision_date}\n")


def convert_pdb_to_cif(pdb_path: Path, cif_path: Path, revision_date: str = DEFAULT_REVISION_DATE):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("surrogate", str(pdb_path))

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(cif_path))

    ensure_release_date_in_cif(cif_path, revision_date=revision_date)


def process_single_pdb(src_pdb: Path, dst_dir: Path,
                       stapled_resnames, surrogate_resname,
                       revision_date: str = DEFAULT_REVISION_DATE) -> Dict:
    result = _manifest_row(source_pdb=str(src_pdb), status="error")

    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(src_pdb.stem, str(src_pdb))
    except Exception as e:
        result["message"] = f"PDB 解析失败: {e}"
        return result

    chain = get_first_chain(structure)
    if chain is None:
        result["message"] = "未找到链。"
        return result

    records, meta = build_surrogate_records(chain, stapled_resnames, surrogate_resname)
    if not records or meta["length"] == 0:
        result["message"] = "未生成有效 surrogate 结构。"
        return result

    stem = src_pdb.stem + "_surrogate"
    out_pdb = dst_dir / f"{stem}.pdb"
    out_cif = dst_dir / f"{stem}.cif"
    out_meta = dst_dir / f"{stem}_meta.json"

    try:
        write_surrogate_pdb(records, out_pdb)
        convert_pdb_to_cif(out_pdb, out_cif, revision_date=revision_date)

        with open(out_meta, "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2, ensure_ascii=False)

        result.update({
            "surrogate_pdb": str(out_pdb),
            "surrogate_cif": str(out_cif),
            "meta_json": str(out_meta),
            "status": "ok",
            "message": "success",
            "length": meta["length"],
            "sequence": meta["sequence"],
            "converted_stapled_count": meta["n_stapled_positions_converted"],
        })
        return result

    except Exception as e:
        result["message"] = f"写出 surrogate 文件失败: {e}"
        return result


# =========================
# sequences.csv 生成
# =========================
def write_sequences_csv(
    output_csv: Path,
    manifest_rows: List[Dict],
    protein_sequence: str,
):
    """
    生成后续 AF3 运行所需的 sequences.csv
    字段包含：
    - file_name
    - full_path
    - protein_sequence
    - peptide_sequence
    """
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "file_name",
        "full_path",
        "protein_sequence",
        "peptide_sequence",
    ]

    rows = []
    for row in manifest_rows:
        if row.get("status") != "ok":
            continue

        surrogate_pdb = row.get("surrogate_pdb", "").strip()
        peptide_sequence = row.get("sequence", "").strip()

        if not surrogate_pdb or not peptide_sequence:
            continue

        rows.append({
            "file_name": Path(row["source_pdb"]).name,
            "full_path": surrogate_pdb,
            "protein_sequence": protein_sequence,
            "peptide_sequence": peptide_sequence,
        })

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


# =========================
# 主流程
# =========================
def run_batch(source_root: Path, target_root: Path, protein_pdb: Path,
              protein_chain_id: str, base_protein_name: str,
              stapled_resnames, surrogate_resname, revision_date: str) -> None:
    target_root.mkdir(parents=True, exist_ok=True)

    manifest_rows = []

    if not source_root.exists():
        raise FileNotFoundError(f"源目录不存在: {source_root}")

    protein_sequence = extract_protein_sequence_from_pdb(protein_pdb, protein_chain_id)
    print(f"[INFO] Protein sequence length ({protein_chain_id} chain): {len(protein_sequence)}")

    template_dirs = sorted([
        p for p in source_root.iterdir()
        if p.is_dir() and p.name.startswith(f"{base_protein_name}_")
    ])

    for template_dir in template_dirs:
        src_cealign = template_dir / "cealign_denovo"
        dst_cealign = target_root / template_dir.name / "cealign_denovo"
        dst_cealign.mkdir(parents=True, exist_ok=True)

        if not src_cealign.exists():
            manifest_rows.append(
                _manifest_row(message=f"{src_cealign} 不存在")
            )
            continue

        pdb_files = sorted(src_cealign.glob("*.pdb"))
        if len(pdb_files) == 0:
            manifest_rows.append(
                _manifest_row(status="empty", message=f"{src_cealign} 为空")
            )
            continue

        for pdb_file in pdb_files:
            row = process_single_pdb(pdb_file, dst_cealign,
                                     stapled_resnames, surrogate_resname,
                                     revision_date=revision_date)
            manifest_rows.append(row)
            print(f"[{row['status']}] {pdb_file}")

    manifest_path = target_root / "surrogate_manifest.csv"
    with open(manifest_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"\nDone. Manifest saved to: {manifest_path}")

    sequences_csv_path = target_root / "sequences.csv"
    write_sequences_csv(
        output_csv=sequences_csv_path,
        manifest_rows=manifest_rows,
        protein_sequence=protein_sequence,
    )
    print(f"Sequences CSV saved to: {sequences_csv_path}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="把订书肽替换为 GLY 代理残基，生成 AF3 输入（PDB/CIF/meta/sequences.csv）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--base-dir", default=DEFAULT_BASE_DIR, help="数据根目录")
    parser.add_argument("--base-folder", default=DEFAULT_BASE_FOLDER, help="数据集文件夹名（如 2KOY）")
    parser.add_argument("--protein-name", default=DEFAULT_PROTEIN_NAME, help="蛋白模板名前缀（如 2koy）")
    parser.add_argument("--source-root", default=None,
                        help="源目录；默认 = base-dir/base-folder/base-folder_100_4")
    parser.add_argument("--target-root", default=None,
                        help="输出目录；默认 = 源目录 + '_cif'")
    parser.add_argument("--protein-pdb", default=None,
                        help="固定蛋白结构 PDB；默认 = source-root/RFdiffusion/{protein-name}_0.pdb")
    parser.add_argument("--protein-chain", default=DEFAULT_PROTEIN_CHAIN_ID, help="蛋白 PDB 中的链 ID")
    parser.add_argument("--stapled-resnames", default=",".join(sorted(DEFAULT_STAPLED_RESNAMES)),
                        help="订书残基名（逗号分隔，将被替换为代理残基）")
    parser.add_argument("--surrogate-resname", default=DEFAULT_SURROGATE_RESNAME, help="代理残基名")
    parser.add_argument("--revision-date", default=DEFAULT_REVISION_DATE,
                        help="写入 CIF 的 revision_date（AF3 兼容）")
    args = parser.parse_args(argv)

    source_root = (Path(args.source_root) if args.source_root
                   else Path(args.base_dir) / args.base_folder / f"{args.base_folder}{SOURCE_FOLDER_SUFFIX}")
    target_root = (Path(args.target_root) if args.target_root
                   else Path(args.base_dir) / args.base_folder / f"{args.base_folder}{TARGET_FOLDER_SUFFIX}")
    protein_pdb = (Path(args.protein_pdb) if args.protein_pdb
                   else source_root / "RFdiffusion" / f"{args.protein_name}_0.pdb")
    stapled_resnames = {name.strip().upper() for name in args.stapled_resnames.split(",") if name.strip()}
    surrogate_resname = args.surrogate_resname.strip().upper()

    if surrogate_resname not in SURROGATE_ALLOWED_ATOMS:
        raise ValueError(f"不支持的代理残基 {surrogate_resname}，可选: {list(SURROGATE_ALLOWED_ATOMS)}")

    run_batch(source_root=source_root,
              target_root=target_root,
              protein_pdb=protein_pdb,
              protein_chain_id=args.protein_chain,
              base_protein_name=args.protein_name,
              stapled_resnames=stapled_resnames,
              surrogate_resname=surrogate_resname,
              revision_date=args.revision_date)
    return 0


if __name__ == "__main__":
    sys.exit(main())