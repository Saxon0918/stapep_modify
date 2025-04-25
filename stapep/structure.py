import os
import uuid
import shutil
import logging
import csv
import gc
import torch

import pandas as pd
from Bio import AlignIO
from Bio.PDB import PDBParser, Superimposer, PDBIO, Select
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio import pairwise2
# from Bio.SubsMat import MatrixInfo as matlist
import Bio.Align.substitution_matrices as matlist
from collections import defaultdict
from io import StringIO

try:
    from Bio.Data.PDBData import protein_letters_3to1 as aa3to1
except ImportError:
    print('Warning: Bio.Data.PDBData is deprecated, please use Bio.Data.SCOPData instead')
    from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1

from stapep.molecular_dynamics import PrepareProt, Simulation
from stapep.utils import PhysicochemicalPredictor, SeqPreProcessing
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'

class Structure(object):
    """
    The Structure class represents a structure and provides methods for generating 3D structures of peptides.

    Args:
        solvent (str, optional): The solvent to use for the simulation. Defaults to 'water'.
        save_tmp_dir (bool, optional): Whether to save the temporary directory used for the simulation. Defaults to False.
        verbose (bool, optional): Whether to enable verbose logging. Defaults to False.

    Examples:
        ```python
        # Create an instance of the Structure class
        structure = Structure()

        # Show the available solvent options
        structure.show_solvent_options()

        # Generate a 3D structure from a template
        structure.generate_3d_structure_from_template('Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ', 'output.pdb', 'template.pdb')

        # Generate a de novo 3D structure
        structure.de_novo_3d_structure('Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ', 'output.pdb')

        # Generate a 3D structure from a sequence
        structure.generate_3d_structure_from_sequence('ACDEFG', 'output.pdb')
        ```
    """

    def __init__(self, solvent: str = 'water', save_tmp_dir: bool = False, verbose: bool = False):
        self.solvent = solvent  # 溶剂，默认值是水溶剂
        self.tmp_dir = os.path.join('/tmp', str(uuid.uuid4()))  # 临时目录路径
        self.save_tmp_dir = save_tmp_dir  # 临时目录是否保存
        self.verbose = verbose  # 是否启动详细的日志记录

        # 如果临时路径不存在文件夹，则新建文件夹
        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)
        # 如果启动详细日志，则记录INFO级别日志
        if verbose:
            logging.basicConfig(level=logging.INFO)
            if self.save_tmp_dir:
                logging.info(f'Temporarily directory: {self.tmp_dir}')

    def show_solvent_options(self):
        print('You can choose the following solvent options:')
        print('- water: default')
        print('- chloroform')
        print('- DMF: dimethylformamide')
        print('- DMSO: dimethyl sulfoxide')
        print('- ethanol')
        print('- acetone')

    def _check_solvent(self, solvent: str):
        if solvent not in ['water', 'chloroform', 'DMF', 'DMSO', 'ethanol', 'acetone']:
            raise ValueError(
                f'{solvent} is not a valid solvent option, please choose from the following options: water, chloroform, DMF, DMSO, ethanol, acetone')

    # 采用分子动力学默认进行100000步
    def _short_time_simulation(self, nsteps: int = 100000):
        sim = Simulation(self.tmp_dir)
        sim.setup(type='implicit',  # 隐式溶剂模型
                  solvent=self.solvent,  # 溶剂类型
                  temperature=300,  # 温度为300K
                  friction=1,  # 摩擦系数控制粒子运动的阻力大小
                  timestep=2,  # 时间步长为 2 fs（飞秒），分子动力学模拟中每一步的时间间隔
                  interval=10,  # 记录间隔为 10 步
                  nsteps=nsteps)  # 模拟步数
        # 是否启动详细日志
        if self.verbose:
            logging.info(f'Running short time simulation for {nsteps} steps')

        try:
            sim.minimize()  # 能量最小化
            sim.run()  # 开始运行模拟
            return True
        except ValueError as e:
            logging.error(f"Simulation error: {e}. Skipping this simulation.")
            return False

    def _get_opt_structure(self, seq, pdb):
        # PhysicochemicalPredictor从分子动力学（MD）轨迹中预测蛋白质的物理化学性质
        pcp = PhysicochemicalPredictor(sequence=seq,
                                       topology_file=os.path.join(self.tmp_dir, 'pep_vac.prmtop'),
                                       # topology file　(default: pep_vac.prmtop in the data folder)
                                       trajectory_file=os.path.join(self.tmp_dir, 'traj.dcd'),
                                       # trajectory file (default: traj.dcd in the data folder)
                                       start_frame=0)  # start frame (default: 500)
        pcp._save_mean_structure(pdb)  # 计算平均结构

    def _del_tmp_dir(self):
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)

    def generate_3d_structure_from_template(self,
                                            seq: str,
                                            output_pdb: str,
                                            template_pdb: str,
                                            additional_residues: dict = None):
        '''
            Generate a 3D structure of a peptide from a template using Modeller.

            Args:
                seq (str): The sequence of the peptide.
                output_pdb (str): The path to save the output PDB file.
                template_pdb (str): The path to the template PDB file.
                additional_residues: additional residues to be added to the system (default: None)
                    For example, {'AIB': ('/path/to/AIB.prepin', '/path/to/frcmod.AIB')
                                  'NLE': ('/path/to/NLE.prepin', '/path/to/frcmod.NLE')}

            Returns:
                str: The path to the generated PDB file.
        '''
        # 如果输入有非标准残基，那就需要额外逻辑处理
        spp = SeqPreProcessing(additional_residues=additional_residues)  # 预处理序列，处理非标准残基
        spp.check_seq_validation(seq)  # 检查输入氨基酸的合法性
        # get absolute path of template pdb file
        template_pdb = os.path.abspath(template_pdb)  # 转化为绝对路径
        # 进行预测，传入序列、临时路径目录、方法、模板pdb文件路径
        pp = PrepareProt(seq, self.tmp_dir, method='modeller', template_pdb_file_path=template_pdb)
        pp._gen_prmtop_and_inpcrd_file()

        self._short_time_simulation()  # 优化结构
        self._get_opt_structure(seq, output_pdb)
        if not self.save_tmp_dir:
            self._del_tmp_dir()
        return output_pdb

    def de_novo_3d_structure(self,
                             seq: str,
                             output_pdb: str,
                             additional_residues: dict = None,
                             proxy=None):
        '''
            Generate a de novo 3D structure of a peptide using ESMFold.

            Args:
                seq (str): The sequence of the peptide.
                output_pdb (str): The path to save the output PDB file.
                additional_residues: additional residues to be added to the system (default: None)
                    For example, {'AIB': ('/path/to/AIB.prepin', '/path/to/frcmod.AIB')
                                  'NLE': ('/path/to/NLE.prepin', '/path/to/frcmod.NLE')}

            Returns:
                str: The path to the generated PDB file.
        '''
        spp = SeqPreProcessing(additional_residues=additional_residues)
        spp.check_seq_validation(seq)
        pp = PrepareProt(seq, self.tmp_dir, method='alphafold', additional_residues=additional_residues)
        pp._gen_prmtop_and_inpcrd_file()
        # self._short_time_simulation()
        simulation_result = self._short_time_simulation()
        if not simulation_result:
            print("Short-time simulation failed. Aborting de novo 3D structure generation.")
            return False  # 返回 None 表示中止

        self._get_opt_structure(seq, output_pdb)
        if not self.save_tmp_dir:
            self._del_tmp_dir()
        return output_pdb

    def generate_3d_structure_from_sequence(self,
                                            seq: str,
                                            output_pdb: str,
                                            additional_residues: dict = None):
        '''
            Generate a 3D structure of a peptide using Ambertools.

            Args:
                seq (str): The sequence of the peptide.
                output_pdb (str): The path to save the output PDB file.
                additional_residues: additional residues to be added to the system (default: None)
                    For example, {'AIB': ('/path/to/AIB.prepin', '/path/to/frcmod.AIB')
                                  'NLE': ('/path/to/NLE.prepin', '/path/to/frcmod.NLE')}

            Returns:
                str: The path to the generated PDB file.

            Note:
                This method is not recommended as the generated structure is not stable.
        '''
        spp = SeqPreProcessing(additional_residues=additional_residues)
        spp.check_seq_validation(seq)
        pp = PrepareProt(seq, self.tmp_dir, method=None)
        pp._gen_prmtop_and_inpcrd_file()
        self._short_time_simulation()
        self._get_opt_structure(seq, output_pdb)
        if not self.save_tmp_dir:
            self._del_tmp_dir()
        return output_pdb


class AlignStructure(object):

    @staticmethod
    def convert_pdb_to_seq(id: str, pdb_file: str) -> list[str]:
        '''
            Convert a pdb file to a fasta file.
        '''
        res_list = AlignStructure._get_pdb_sequence(id, pdb_file)
        seq = [res[1] for res in res_list]
        seq = ''.join(seq)
        return seq

    @staticmethod
    def _get_pdb_sequence(id, pdb_file) -> list[tuple]:
        '''
            Return a list of tuples (idx, sequence).
            eg:[(6, 'P'),
                (7, 'D'),
                (8, 'I'),
                (9, 'F'),]
        '''
        parser = PDBParser()
        structure = parser.get_structure(id, pdb_file)
        _aainfo = lambda r: (r.id[1], aa3to1.get(r.resname, 'X'))
        return [_aainfo(r) for r in structure.get_residues() if is_aa(r)]

    @staticmethod
    def align(ref_pdb: str, pdb: str, output_pdb: str):
        '''
            Align structures using BioPython.

            Args:
                ref_pdb (str): The path to the reference PDB file.
                pdb (str): The path to the PDB file to align.
                output_pdb (str): The path to save the output PDB file.

            Returns:
                str: The path to the generated PDB file.
        '''

        try:
            import pymol
            from pymol import cmd
            cmd.reinitialize()
        except Exception as e:
            raise ImportError('Please install PyMOL to use this method. mamba install -c conda-forge pymol-open-source')

        # Load the PDB files
        pymol.cmd.load(ref_pdb, 'ref')  # 加载模板文件
        pymol.cmd.load(pdb, 'denovo')  # 加载预测pdb文件

        # Perform alignment on alpha carbons (CA atoms)
        out = pymol.cmd.align('denovo and name CA', 'ref and name CA')  # 基于α-碳原子（CA 原子）进行对齐
        # RMSD、对齐原子数量、迭代次数、对齐前RMSD、对齐前对齐的原子数量、对齐评分、对齐涉及的氨基酸残基数量
        rmsd, n_atoms, n_cycles, n_rmsd_pre, n_atom_pre, score, n_res = out

        # Make sure to update coordinates of the denovo structure
        pymol.cmd.alter_state(1, 'denovo', 'x, y, z = x, y, z')  # 确保目标结构denovo中的坐标在对齐后得到更新
        # Apply the transformation matrix after alignment
        pymol.cmd.matrix_copy('denovo', 'ref')  # 确保目标结构与参考结构完全对齐
        # Save the aligned denovo structure to a new PDB file
        pymol.cmd.save(output_pdb, 'denovo')  # 保存新的目标结构
        return rmsd  # 返回RMSD

    # @staticmethod
    # def align(ref_pdb: str, pdb: str, output_pdb: str):
    #     '''
    #         Align structures using BioPython.

    #         Args:
    #             ref_pdb (str): The path to the reference PDB file.
    #             pdb (str): The path to the PDB file to align.
    #             output_pdb (str): The path to save the output PDB file.

    #         Returns:
    #             str: The path to the generated PDB file.
    #     '''

    #     ref_id = os.path.basename(ref_pdb).split('.')[0]
    #     pdb_id = os.path.basename(pdb).split('.')[0]

    #     parser = PDBParser()
    #     ref_structure = parser.get_structure(ref_id, ref_pdb)
    #     pdb_structure = parser.get_structure(pdb_id, pdb)
    #     ref_model = list(ref_structure.get_models())[0]
    #     pdb_model = list(pdb_structure.get_models())[0]

    #     resseq_A = AlignStructure._get_pdb_sequence(ref_id, ref_pdb)
    #     resseq_B = AlignStructure._get_pdb_sequence(pdb_id, pdb)
    #     sequence_A = AlignStructure.convert_pdb_to_seq(ref_id, ref_pdb)
    #     sequence_B = AlignStructure.convert_pdb_to_seq(pdb_id, pdb)

    #     alns = pairwise2.align.globalds(sequence_A, sequence_B, matlist.load("BLOSUM62"), -10.0, -0.5,
    #                                         penalize_end_gaps=(False, False) )
    #     best_aln = alns[0]
    #     aligned_A, aligned_B, score, begin, end = best_aln
    #     mapping = {}
    #     aa_i_A, aa_i_B = 0, 0
    #     for aa_aln_A, aa_aln_B in zip(aligned_A, aligned_B):
    #         if aa_aln_A == '-':
    #             if aa_aln_B != '-':
    #                 aa_i_B += 1
    #         elif aa_aln_B == '-':
    #             if aa_aln_A != '-':
    #                 aa_i_A += 1
    #         else:
    #             assert resseq_A[aa_i_A][1] == aa_aln_A
    #             assert resseq_B[aa_i_B][1] == aa_aln_B
    #             mapping[resseq_A[aa_i_A][0]] = resseq_B[aa_i_B][0]
    #             aa_i_A += 1
    #             aa_i_B += 1

    #     # Extract CA atoms using the helper function
    #     refe_ca_list = AlignStructure.get_CA_atoms_from_model(ref_model, list(mapping.keys()))
    #     mobi_ca_list = AlignStructure.get_CA_atoms_from_model(pdb_model, list(mapping.values()))

    #     # Superimpose matching residues
    #     try:
    #         si = Superimposer()
    #         si.set_atoms(refe_ca_list, mobi_ca_list)
    #         si.apply(pdb_structure.get_atoms())

    #         io = PDBIO()
    #         io.set_structure(pdb_structure)
    #         io.save(output_pdb)
    #         return output_pdb
    #     except Exception as e:
    #         print(e)

    @staticmethod
    def rmsd(ref_pdb: str, pdb: str):
        try:
            import pymol
            from pymol import cmd
            cmd.reinitialize()
        except Exception as e:
            print(e)
            return None

        pymol.cmd.load(ref_pdb, 'ref')
        pymol.cmd.load(pdb, 'denovo')
        out = pymol.cmd.align('ref and name CA', 'denovo and name CA')
        rmsd, n_atoms, n_cyles, n_rmsd_pre, n_atom_pre, score, n_res = out
        return rmsd

    @staticmethod
    def get_CA_atoms_from_model(model, residue_numbers):
        """
        Get CA atoms from a given model based on specified residue numbers.

        Args:
        - model: BioPython Model object
        - residue_numbers: List of residue numbers to extract CA atoms for

        Returns:
        - List of CA atoms
        """
        ca_atoms = []

        for chain in model:
            ca_atoms.extend(
                res['CA']
                for res in chain
                if res.id[1] in residue_numbers and 'CA' in res
            )
        return ca_atoms


def calculate_all_kinds_stapep(seq, st):
    stapep_ways = ['S5_S5_3', 'R8_S5_6', 'R5_S8_6', 'R5_S5_2', 'S5_S5_4', 'R8_S5_5', 'R8_S5_3', 'S5_R8_6', 'R5_S8_7']
    columns = [
        "Iteration", "seq", "helix percent", "sheet percent", "loop percent",
        "mean bfactor", "mol surf", "mean gyrate", "psa",
        "total number of hydrogen bonds", "smiles"
    ]
    file_path = 'example/insulin/insulin.txt'
    if not os.path.exists(file_path):
        with open(file_path, 'w') as file:
            file.write('\t'.join(columns) + '\n')  # 写入列名

    # 使用 for 循环遍历 seq，生成每一轮的 new_seq
    for stapep_way in stapep_ways:
        stapep_first, stapep_second, interval = stapep_way.split('_')
        interval = int(interval)

        for i in range(0, len(seq)):
            new_seq = seq[:i]
            # 插入stapep_first和stapep_second后，跳过interval个字符
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

            # 同源建模，通过相似蛋白质序列的结构预测三维结构
            # st.generate_3d_structure_from_template(seq=seq, output_pdb='example/data/homology_model_test0.pdb',
            #                                        template_pdb='example/data/template.pdb')
            # 将预测出的pdb中的C原子与模板中对齐，三维坐标对齐
            # AlignStructure.align(ref_pdb='example/data/template.pdb', pdb='example/data/homology_model.pdb',
            #                      output_pdb='example/data/aligned.pdb')
            #
            # # 采用AmberTools工具建模
            # st.generate_3d_structure_from_sequence(seq=seq, output_pdb='example/data/sequence.pdb')

            # 去模板预测，采用ESMFold
            de_novo_path = f"example/insulin/denovo_{stapep_first}_{stapep_second}_{interval}_{i}.pdb"
            output_de_novo = st.de_novo_3d_structure(seq=new_seq, output_pdb=de_novo_path)
            if not output_de_novo:
                record = f"Iteration {stapep_first}_{stapep_second}_{interval}_{i} - error"
                # 追加记录到文件
                with open(file_path, 'a') as file:  # 使用追加模式
                    file.write(record + '\n')
                continue
            pathname = st.tmp_dir  # The path where the topology file and trajectory file are located
            # pathname = '/tmp/9fbdc077-e700-4bb2-a959-9bd3f8ab8192'
            pcp = PhysicochemicalPredictor(sequence=new_seq,
                                           topology_file=os.path.join(pathname, 'pep_vac.prmtop'),
                                           # topology file　(default: pep_vac.prmtop in the data folder)
                                           trajectory_file=os.path.join(pathname, 'traj.dcd'),
                                           # trajectory file (default: traj.dcd in the data folder)
                                           start_frame=0)  # start frame (default: 500)

            # Get the features
            helix_percent = pcp.calc_helix_percent()
            sheet_percent = pcp.calc_extend_percent()
            loop_percent = pcp.calc_loop_percent()

            print('helix percent: ', helix_percent)  # 计算α螺旋占比
            print('sheet percent: ', sheet_percent)  # 计算β折叠占比
            print('loop percent: ', loop_percent)  # 计算环结构（loop）占比
            # save the mean structure of the trajectory
            mean_structure_path = f"example/insulin/mean_structure_{stapep_first}_{stapep_second}_{interval}_{i}.pdb"
            pcp._save_mean_structure(mean_structure_path)  # 保存分子的平均结构的pdb

            # calculate the Mean B-factor, Molecular Surface, Mean Gyration Radius, Hydrophobic Index, and 3D-PSA
            mean_bfactor = pcp.calc_mean_bfactor()
            mol_surf = pcp.calc_mean_molsurf()
            mean_gyrate = pcp.calc_mean_gyrate()
            psa = pcp.calc_psa(mean_structure_path)
            total_number_hydrogen_bonds = pcp.calc_n_hbonds()
            print('mean bfactor: ', mean_bfactor)  # 平均B因子，较低的B因子表明原子较稳定，柔性较小
            print('mol surf: ', mol_surf)  # 分子表面积
            print('mean gyrate: ', mean_gyrate)  # 分子的旋转半径，分子质量分布的中心到外部的平均距离，评估分子的形状和紧密性
            # print('hydrophobic index: ',
            #       pcp.calc_hydrophobic_index(os.path.join(pathname, 'mean_structure.pdb')))  # TODO 疏水指数环境没配好
            print('psa: ', psa)  # 极性表面积，评价分子吸收、渗透和跨膜运输能力
            print('total number of hydrogen bonds: ', total_number_hydrogen_bonds)  # 氢键的总数，决定了分子的稳定性和结合能力

            # extract 2D structure of the peptide
            smiles = pcp.extract_2d_structure(mean_structure_path)
            print(smiles)
            # 构建记录内容
            record = f"{stapep_first}_{stapep_second}_{interval}_{i}\t{new_seq}\t{helix_percent}\t{sheet_percent}\t{loop_percent}\t" \
                     f"{mean_bfactor}\t{mol_surf}\t{mean_gyrate}\t{psa}\t" \
                     f"{total_number_hydrogen_bonds}\t{smiles}"

            # 追加记录到文件
            with open(file_path, 'a') as file:  # 使用追加模式
                file.write(record + '\n')

def generate_pdb(st):
    """

    Args:
        st: 传入的Structure(verbose=True)样本

    Returns: 采用PDB数据集中2176条数据的Sequence生成的订书肽PDB文件

    """

    # file_path = 'example/datasets/Filtered_peptides.csv'
    # file_path = 'example/datasets/filter_stapep_amp.csv'
    file_path = 'example/datasets/0001-3906_sup.csv'
    sequence_column = []
    index_column = []

    with open(file_path, mode='r', encoding='utf-8') as file:
        csv_reader = csv.reader(file)
        headers = next(csv_reader)
        sequence_index = headers.index('seq')
        stapep_id = headers.index('stapep_id')
        for row in csv_reader:
            sequence_column.append(row[sequence_index])
            index_column.append(row[stapep_id])

    # for i in range(4, len(sequence_column)):
    for i in range(101, 19993):
        seq = sequence_column[i]
        torch.cuda.empty_cache()
        gc.collect()
        de_novo_path = f"example/stapep_data/pdb_sup/stapep_{index_column[i]}.pdb"
        output_de_novo = st.de_novo_3d_structure(seq=seq, output_pdb=de_novo_path)
        if not output_de_novo:
            record = f"Iteration {i} - error"
            print(record)
    print("FINISH")


class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id


def random_insert_stapep_peptide(pdb_file):
    """

    Args:
        pdb_file: RFdiffusion生成的多肽骨架

    Returns: 依次添加S5订书肽后采用Modeller生成的订书肽骨架

    """

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    # 获取A链数据单独保存为一个pdb文件
    io = PDBIO()
    io.set_structure(structure)
    io.save("example/data/gsk3beta/1gng_5/chain_A.pdb", select=ChainSelect("A"))

    # 获取 A 链的残基数量
    chain_A = structure[0]['A']
    residues = [res for res in chain_A if res.id[0] == ' ']
    residue_count = len(residues)
    peptide_sequence = list('G' * residue_count)

    # 向RFDiffusion中依次添加S5订书肽，然后采用Modeller肽链
    for i in range(residue_count - 4):
        peptide_sequence_copy = peptide_sequence.copy()
        peptide_sequence_copy[i] = 'S5'
        peptide_sequence_copy[i+4] = 'S5'
        seq = ''.join(peptide_sequence_copy)
        st = Structure(verbose=True)
        output_pdb_path = fr'example/data/gsk3beta/1gng_5/homology_model_{i}.pdb'
        template_pdb_path = 'example/data/gsk3beta/1gng_5/chain_A.pdb'
        st.generate_3d_structure_from_template(seq=seq,
                                               output_pdb=output_pdb_path,
                                               template_pdb=template_pdb_path)

        align_output_pdb_path = fr'example/data/gsk3beta/1gng_5/align/aligned_homology_model_{i}.pdb'
        AlignStructure.align(ref_pdb=template_pdb_path, pdb=output_pdb_path, output_pdb=align_output_pdb_path)

if __name__ == '__main__':
    # 1. 初始转换（Initial Conversion）
    # 目标：将非标准氨基酸转换为 ALA（丙氨酸），简化建模初期的处理。
    # 原因：非标准氨基酸可能缺乏标准的三维结构模板，因此暂时用ALA替代，有助于初步建模。
    #
    # 2. 肽链结构建模（Peptide Structure Modeling）
    # 包含三种方法：
    # （1）同源建模（Homology Modeling）：使用 Modeller 工具。假设相似序列具有类似的三维结构，选择已知结构作为模板，通过比对序列生成目标肽链的模型。
    # （2）去模板预测（De Novo Prediction）：使用 ESMFold，根据肽链序列直接预测三维结构，不依赖模板。
    # （3）AmberTools：使用其 tleap 模块从肽链序列生成初步三维结构。但是准确性相对较低，这种方法通常不被推荐。
    #
    # 3. 氨基酸重转换（Reconversion of Amino Acids）：在建模完成后，将所有临时替换为ALA的非标准氨基酸恢复到其原始命名。
    # 4. 结构完善（Structural Completion）：使用 AmberTools 的 tleap 模块，调整三维结构，确保非标准氨基酸正确定位并被包含在模型中。
    # 5. 动力学优化（Dynamics Optimization）：使用 OpenMM 工具对模型进行分子动力学（MD）模拟。模拟短时间（如100皮秒，ps）的动力学行为，优化模型的稳定性和精确度。

    random_insert_stapep_peptide('example/data/gsk3beta/1gng_5/1gng_5.pdb')

    # --------------------- Modeller ---------------------
    # seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ'
    # st = Structure(verbose=True)
    # st.generate_3d_structure_from_template(seq=seq,
    #                                        output_pdb='example/data/homology_model.pdb',
    #                                        template_pdb='example/data/template.pdb')
    #
    # AlignStructure.align(ref_pdb='example/data/template.pdb',
    #                      pdb='example/data/homology_model.pdb',
    #                      output_pdb='example/data/aligned.pdb')

    # --------------------- Denove ---------------------
    # seq = 'S5KKIS5KKIKKKLK'
    # st = Structure(verbose=True, save_tmp_dir=True)
    # # 去模板预测，采用ESMFold
    # de_novo_path = f"example/stapep_data/test/pred_stapep_58075.pdb"
    # output_de_novo = st.de_novo_3d_structure(seq=seq, output_pdb=de_novo_path)
    # if not output_de_novo:
    #     record = f"Iteration {385} - error"
    #     print(record)

    # generate_pdb(st)



