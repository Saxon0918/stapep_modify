import os
import sys
import uuid
import shutil
import logging
from pathlib import Path

# 支持 `python stapep/structure.py` 直接运行时导入 stapep 包
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from parmed.exceptions import OpenMMError

from stapep.molecular_dynamics import PrepareProt, Simulation
from stapep.utils import PhysicochemicalPredictor, SeqPreProcessing

try:
    from Bio.Data.PDBData import protein_letters_3to1 as aa3to1
except ImportError:
    print('Warning: Bio.Data.PDBData is deprecated, please use Bio.Data.SCOPData instead')
    from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa


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
                  timestep=1,  # 时间步长为 1 fs（飞秒），分子动力学模拟中每一步的时间间隔
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
        except OpenMMError as e:
            logging.error(f"OpenMM error: {e}. Skipping this simulation.")
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
        try:
            spp = SeqPreProcessing(additional_residues=additional_residues)
            spp.check_seq_validation(seq)
        except Exception as e:
            print(f"[ERROR] Sequence validation failed for seq={seq}: {e}")
            return False

        try:
            pp = PrepareProt(
                seq,
                self.tmp_dir,
                method='alphafold',
                additional_residues=additional_residues
            )
            pp._gen_prmtop_and_inpcrd_file()
        except Exception as e:
            print(f"[ERROR] Failed to prepare topology/coordinate files for seq={seq}: {e}")
            return False

        try:
            simulation_result = self._short_time_simulation()
            if not simulation_result:
                print(f"[ERROR] Short-time simulation failed for seq={seq}")
                return False
        except Exception as e:
            print(f"[ERROR] Exception during short-time simulation for seq={seq}: {e}")
            return False

        try:
            self._get_opt_structure(seq, output_pdb)
        except Exception as e:
            print(f"[ERROR] Failed to generate optimized structure for seq={seq}, output={output_pdb}: {e}")
            return False

        try:
            if not self.save_tmp_dir:
                self._del_tmp_dir()
        except Exception as e:
            print(f"[WARNING] Failed to delete tmp dir for seq={seq}: {e}")

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
    def convert_pdb_to_seq(id: str, pdb_file: str) -> str:
        '''
            Convert a pdb file to a sequence string.
        '''
        res_list = AlignStructure._get_pdb_sequence(id, pdb_file)
        seq = [res[1] for res in res_list]
        seq = ''.join(seq)
        return seq

    @staticmethod
    def _get_pdb_sequence(id, pdb_file) -> list:
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
            Align structures using PyMOL.

            Args:
                ref_pdb (str): The path to the reference PDB file.
                pdb (str): The path to the PDB file to align.
                output_pdb (str): The path to save the output PDB file.

            Returns:
                float: The RMSD of the alignment.
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

    @staticmethod
    def align_denovo(ref_pdb: str, pdb: str, output_pdb: str, rfdiffusion_template: str):
        '''
            Align structures using PyMOL.

            Args:
                ref_pdb (str): The path to the reference PDB file.
                pdb (str): The path to the PDB file to align.
                output_pdb (str): The path to save the output PDB file.
                rfdiffusion_template (str): The RFdiffusion template name,
                    used to decide which chain of the reference to align to.

            Returns:
                float: The RMSD of the alignment.
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
        # for 1gng chain X; for 2gv2 chain B
        if rfdiffusion_template.startswith("1gng"):
            out = pymol.cmd.align('denovo', 'ref and chain A')  # 基于α-碳原子（CA 原子）进行对齐
        else:
            out = pymol.cmd.align('denovo and name CA', 'ref and chain A and name CA')  # 基于α-碳原子（CA 原子）进行对齐
        # RMSD、对齐原子数量、迭代次数、对齐前RMSD、对齐前对齐的原子数量、对齐评分、对齐涉及的氨基酸残基数量
        rmsd, n_atoms, n_cycles, n_rmsd_pre, n_atom_pre, score, n_res = out

        # Make sure to update coordinates of the denovo structure
        pymol.cmd.alter_state(1, 'denovo', 'x, y, z = x, y, z')  # 确保目标结构denovo中的坐标在对齐后得到更新
        # Apply the transformation matrix after alignment
        pymol.cmd.matrix_copy('denovo', 'ref')  # 确保目标结构与参考结构完全对齐
        # Save the aligned denovo structure to a new PDB file
        pymol.cmd.save(output_pdb, 'denovo')  # 保存新的目标结构
        return rmsd  # 返回RMSD

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


if __name__ == '__main__':
    # 保持旧入口兼容：原 structure.py 的 __main__（RFdiffusion 订书肽流水线）
    # 已拆分到 stapep.rfdiffusion_denovo，此处直接转发。
    from stapep.rfdiffusion_denovo import main

    sys.exit(main())