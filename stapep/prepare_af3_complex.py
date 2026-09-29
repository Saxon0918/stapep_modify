import os
from Bio.PDB import MMCIFParser, PDBIO, Select


# =========================================================
# config
# =========================================================

ROOT_PATH = "/home/d3008/Documents/IPAD/IPAD_TOP"

AF3_ROOT = os.path.join(
    ROOT_PATH,
    "AF3_result"
)

PROTEIN_CHAIN = "A"
PEPTIDE_CHAIN = "B"


# =========================================================
# chain selector
# =========================================================

class ChainSelect(Select):

    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id


# =========================================================
# sample name
# =========================================================

def extract_sample_name(filename):
    """
    Example:
    ipad_21_9_28_ipad_21_9_28_model.cif
        ->
    ipad_21_9_28
    """

    base = filename.replace(".cif", "")

    if "_model" in base:
        base = base.split("_model")[0]

    parts = base.split("_")

    half = len(parts) // 2

    return "_".join(parts[:half])


# =========================================================
# cif -> pdb
# =========================================================

def prepare_complex_from_cif(
        cif_path,
        output_dir,
        protein_chain='A',
        peptide_chain='B'
):

    structure_name = os.path.basename(cif_path)

    parser = MMCIFParser(QUIET=True)

    structure = parser.get_structure(
        structure_name,
        cif_path
    )

    os.makedirs(output_dir, exist_ok=True)

    protein_pdb = os.path.join(
        output_dir,
        "protein.pdb"
    )

    peptide_pdb = os.path.join(
        output_dir,
        "peptide.pdb"
    )

    complex_pdb = os.path.join(
        output_dir,
        "complex_B_NCACO.pdb"
    )

    io = PDBIO()

    # ---------------- protein ----------------
    io.set_structure(structure)

    io.save(
        protein_pdb,
        ChainSelect(protein_chain)
    )

    # ---------------- peptide ----------------
    io.set_structure(structure)

    io.save(
        peptide_pdb,
        ChainSelect(peptide_chain)
    )

    # ---------------- merge ----------------
    with open(complex_pdb, 'w') as outfile:

        with open(protein_pdb, 'r') as infile:

            for line in infile:

                if not line.startswith("END"):
                    outfile.write(line)

        with open(peptide_pdb, 'r') as infile:

            for line in infile:

                if not line.startswith("END"):
                    outfile.write(line)

        outfile.write("END\n")

    print(f"Generated: {complex_pdb}")


# =========================================================
# main
# =========================================================

if __name__ == "__main__":

    cif_files = sorted([
        f for f in os.listdir(AF3_ROOT)
        if f.endswith(".cif")
    ])

    print(f"Total cif files: {len(cif_files)}")

    for cif_file in cif_files:

        try:

            sample_name = extract_sample_name(cif_file)

            print(f"\nProcessing: {sample_name}")

            cif_path = os.path.join(
                AF3_ROOT,
                cif_file
            )

            # output folder
            output_dir = os.path.join(
                ROOT_PATH,
                sample_name,
                "openbpmd",
                sample_name
            )

            protein_pdb = os.path.join(
                output_dir,
                "protein.pdb"
            )

            peptide_pdb = os.path.join(
                output_dir,
                "peptide.pdb"
            )

            complex_pdb = os.path.join(
                output_dir,
                "complex_B_NCACO.pdb"
            )

            # =====================================================
            # skip already processed
            # =====================================================

            if (
                    os.path.exists(protein_pdb)
                    and os.path.exists(peptide_pdb)
                    and os.path.exists(complex_pdb)
            ):
                print(f"{sample_name} already processed, skip")

                continue

            # =====================================================
            # prepare output dir
            # =====================================================

            os.makedirs(
                output_dir,
                exist_ok=True
            )

            # =====================================================
            # preprocess
            # =====================================================

            prepare_complex_from_cif(
                cif_path,
                output_dir,
                protein_chain=PROTEIN_CHAIN,
                peptide_chain=PEPTIDE_CHAIN
            )

        except Exception as e:

            print(f"Error processing {cif_file}")
            print(str(e))