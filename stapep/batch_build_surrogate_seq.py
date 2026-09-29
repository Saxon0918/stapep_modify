import os
import csv
import json
from pathlib import Path
from typing import Dict, List, Tuple

from Bio.PDB import PDBParser
from Bio.PDB import MMCIFIO
from Bio.SeqUtils import seq1


# =========================
# 路径配置
# =========================
BASE_FOLDER = "2KOY"
BASE_PROTEINNAME = "2koy"
SOURCE_ROOT = Path(f"/home/d3008/Documents/{BASE_FOLDER}/{BASE_FOLDER}_100_4")
TARGET_ROOT = Path(f"/home/d3008/Documents/{BASE_FOLDER}/{BASE_FOLDER}_100_4_cif")

# 固定的 protein 结构文件
PROTEIN_PDB = f"/home/d3008/Documents/{BASE_FOLDER}/{BASE_FOLDER}_100_4/RFdiffusion/{BASE_PROTEINNAME}_0.pdb"
PROTEIN_CHAIN_ID = "B"

# # =========================
# # 路径配置
# # =========================
# BASE_FOLDER = "MDM2"
# BASE_PROTEINNAME = "mdm2"
# SOURCE_ROOT = Path(f"/home/d3008/Documents/MDM2/MDM2_after_sample/mdm2_0")
# TARGET_ROOT = Path(f"/home/d3008/Documents/MDM2/MDM2_after_sample_cif/mdm2_0")
#
# # 固定的 protein 结构文件
# PROTEIN_PDB = f"/home/d3008/Documents/MDM2/MDM2_after_sample/mdm2_0/mdm2_0_complex.pdb"
# PROTEIN_CHAIN_ID = "B"

# 当前方案：PS5 -> GLY
STAPLED_RESNAMES = {"PS5", "S5"}
SURROGATE_RESNAME = "GLY"

STANDARD_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIE": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

SURROGATE_ALLOWED_ATOMS = {
    "GLY": {"N", "CA", "C", "O", "OXT"},
}


# =========================
# 通用工具
# =========================
def clean_resname(resname: str) -> str:
    return resname.strip().upper()


def is_hydrogen(atom) -> bool:
    element = (atom.element or "").strip().upper()
    name = atom.get_name().strip().upper()
    return element == "H" or name.startswith("H")


def is_supported_residue(residue) -> bool:
    resname = clean_resname(residue.get_resname())
    return (resname in STANDARD_3TO1) or (resname in STAPLED_RESNAMES)


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
def build_surrogate_records(chain) -> Tuple[List[Dict], Dict]:
    records = []
    residue_meta = []
    sequence_chars = []

    atom_serial = 1
    new_resseq = 1
    stapled_count = 0

    for residue in chain:
        if not is_supported_residue(residue):
            continue

        original_resname = clean_resname(residue.get_resname())
        is_stapled = original_resname in STAPLED_RESNAMES

        if is_stapled:
            surrogate_resname = SURROGATE_RESNAME
            allowed_atoms = SURROGATE_ALLOWED_ATOMS[SURROGATE_RESNAME]
            stapled_count += 1
        else:
            surrogate_resname = original_resname
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
                "resname": surrogate_resname,
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
        sequence_chars.append(STANDARD_3TO1[surrogate_resname])

        hetflag, orig_resseq, orig_icode = residue.id
        residue_meta.append({
            "new_chain_id": "A",
            "new_resseq": new_resseq,
            "source_chain_id": str(chain.id).strip(),
            "source_hetflag": str(hetflag).strip(),
            "source_resseq": int(orig_resseq),
            "source_icode": str(orig_icode).strip(),
            "source_resname": original_resname,
            "surrogate_resname": surrogate_resname,
            "is_stapled_source": is_stapled,
        })

        new_resseq += 1

    meta = {
        "source_chain_id": str(chain.id).strip(),
        "surrogate_chain_id": "A",
        "surrogate_resname_for_stapled": SURROGATE_RESNAME,
        "stapled_source_resnames": sorted(list(STAPLED_RESNAMES)),
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


def ensure_release_date_in_cif(cif_path: Path, revision_date: str = "2024-01-01"):
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


def convert_pdb_to_cif(pdb_path: Path, cif_path: Path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("surrogate", str(pdb_path))

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(cif_path))

    ensure_release_date_in_cif(cif_path, revision_date="2024-01-01")


def process_single_pdb(src_pdb: Path, dst_dir: Path) -> Dict:
    result = {
        "source_pdb": str(src_pdb),
        "surrogate_pdb": "",
        "surrogate_cif": "",
        "meta_json": "",
        "status": "error",
        "message": "",
        "length": 0,
        "sequence": "",
        "converted_stapled_count": 0,
    }

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

    records, meta = build_surrogate_records(chain)
    if not records or meta["length"] == 0:
        result["message"] = "未生成有效 surrogate 结构。"
        return result

    stem = src_pdb.stem + "_surrogate"
    out_pdb = dst_dir / f"{stem}.pdb"
    out_cif = dst_dir / f"{stem}.cif"
    out_meta = dst_dir / f"{stem}_meta.json"

    try:
        write_surrogate_pdb(records, out_pdb)
        convert_pdb_to_cif(out_pdb, out_cif)

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
def main():
    TARGET_ROOT.mkdir(parents=True, exist_ok=True)

    manifest_rows = []

    if not SOURCE_ROOT.exists():
        raise FileNotFoundError(f"源目录不存在: {SOURCE_ROOT}")

    protein_sequence = extract_protein_sequence_from_pdb(
        Path(PROTEIN_PDB),
        PROTEIN_CHAIN_ID,
    )
    print(f"[INFO] Protein sequence length ({PROTEIN_CHAIN_ID} chain): {len(protein_sequence)}")

    mdm2_dirs = sorted([
        p for p in SOURCE_ROOT.iterdir()
        if p.is_dir() and p.name.startswith(f"{BASE_PROTEINNAME}_")
    ])

    # mdm2_dirs = [1]
    for mdm2_dir in mdm2_dirs:
        src_cealign = mdm2_dir / "cealign_denovo"
        # src_cealign = Path("/home/d3008/Documents/MDM2/MDM2_after_sample/mdm2_0/cealign_denovo")
        dst_cealign = TARGET_ROOT / mdm2_dir.name / "cealign_denovo"
        # dst_cealign = Path("/home/d3008/Documents/MDM2/MDM2_after_sample_cif/mdm2_0/cealign_denovo")
        dst_cealign.mkdir(parents=True, exist_ok=True)

        if not src_cealign.exists():
            manifest_rows.append({
                "source_pdb": "",
                "surrogate_pdb": "",
                "surrogate_cif": "",
                "meta_json": "",
                "status": "skip",
                "message": f"{src_cealign} 不存在",
                "length": 0,
                "sequence": "",
                "converted_stapled_count": 0,
            })
            continue

        pdb_files = sorted(src_cealign.glob("*.pdb"))
        if len(pdb_files) == 0:
            manifest_rows.append({
                "source_pdb": "",
                "surrogate_pdb": "",
                "surrogate_cif": "",
                "meta_json": "",
                "status": "empty",
                "message": f"{src_cealign} 为空",
                "length": 0,
                "sequence": "",
                "converted_stapled_count": 0,
            })
            continue

        for pdb_file in pdb_files:
            row = process_single_pdb(pdb_file, dst_cealign)
            manifest_rows.append(row)
            print(f"[{row['status']}] {pdb_file}")

    manifest_path = TARGET_ROOT / "surrogate_manifest.csv"
    fieldnames = [
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

    with open(manifest_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"\nDone. Manifest saved to: {manifest_path}")

    sequences_csv_path = TARGET_ROOT / "sequences.csv"
    write_sequences_csv(
        output_csv=sequences_csv_path,
        manifest_rows=manifest_rows,
        protein_sequence=protein_sequence,
    )
    print(f"Sequences CSV saved to: {sequences_csv_path}")


if __name__ == "__main__":
    main()

