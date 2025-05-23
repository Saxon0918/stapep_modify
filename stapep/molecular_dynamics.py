import os
import subprocess
import time
import shutil
import requests
# OpenMM Imports
import openmm as mm
import openmm.app as app

# ParmEd Imports
import parmed as pmd
from parmed import unit as u
from utils import SeqPreProcessing

class PrepareProt(object):
    '''
        Prepare protein structure

        Args:
            seq (str): protein sequence
            output (str): output file name
            method (str): if method is None, the protein structure is from sequence using AmberTools. 
                          if method is 'alphafold', the protein structure is from ESMFold.
                          if method is 'modeller', the protein structure is from homology modeling using Modeller.
            additional_residues: additional residues to be added to the system (default: None)
                For example, {'AIB': ('/path/to/AIB.prepin', '/path/to/frcmod.AIB')
                                'NLE': ('/path/to/NLE.prepin', '/path/to/frcmod.NLE')}

    '''

    def __init__(self, seq: str, output: str, method: str = None, template_pdb_file_path: str = None,
                 additional_residues: dict = None) -> None:
        self.seq = seq  # 序列
        self.output = output  # 输出文件的路径
        if not os.path.exists(self.output):
            os.makedirs(self.output)
        self.additional_residues = additional_residues  # 非标准残基
        self.seqpp = SeqPreProcessing(self.additional_residues)  # 处理非标准残基
        self.method = method  # Modeller、ESMFold、AmberTools
        self.template_pdb_file_path = template_pdb_file_path  # 依据模板预测的话就是模板路径，例如Modeller方法
        self._check_programs_installed()  # AmberTools方法需要确认是否有tleap

    def _check_programs_installed(self) -> None:
        os.environ['PATH'] += os.pathsep + '/data/zxd/miniconda3/envs/stap/bin'
        if shutil.which('tleap') is None:
            raise Exception('tleap not found. Please install AmberTools.')

    @property
    def _has_cap(self) -> bool:
        return self.seqpp._seq_to_list(self.seq)[0] in ['Ac', 'ACE']

    @property
    def _has_ter(self) -> bool:
        return self.seqpp._seq_to_list(self.seq)[-1] in ['NH2', 'NME']

    def _calculate_atom_coord(self, lines, step, orientation, offset):
        CA_coord = [float(lines[step][30:38]), float(lines[step][38:46]), float(lines[step][46:54])]
        atom_coord = [CA_coord[0] + offset * orientation[0], CA_coord[1] + offset * orientation[1],
                      CA_coord[2] + offset * orientation[2]]
        return atom_coord

    def _insert_ace_and_nme(self, lines: str) -> str:
        '''
            inserts ACE and NME residues into the protein's PDB file.
            This function calculates the coordinates for the ACE and NME atoms based on the coordinates of the first and last residues in the protein, respectively, 
            and then inserts the ACE and NME atoms into the appropriate positions in the PDB file.
        '''
        # add ACE and NME to the PDB file
        # calculate atom coordinates of ACE according to the first residue
        if self._has_cap:
            for step, line in enumerate(lines):
                if line.startswith('ATOM') and line[13:15] == 'CA':
                    try:
                        orientation = [float(lines[step + 1][30:38]) - float(lines[step - 1][30:38]),
                                       float(lines[step + 1][38:46]) - float(lines[step - 1][38:46]),
                                       float(lines[step + 1][46:54]) - float(lines[step - 1][46:54])]
                        ace_atom_coord = self._calculate_atom_coord(lines, step, orientation, -1)
                        ace_atom_line = f'ATOM  {len(lines) + 1:5d}  C   ACE A   0    {ace_atom_coord[0]:8.3f}{ace_atom_coord[1]:8.3f}{ace_atom_coord[2]:8.3f}  1.00  0.00           C  '
                        # insert ACE to the first line of ATOM
                        lines.insert(step - 1, ace_atom_line)
                        break
                    except IndexError:
                        print("Error: Failed to calculate ACE atom coordinates")
                        break
        # calculate atom coordinates of NME according to the last residue
        if self._has_ter:
            for step, line in enumerate(reversed(lines)):
                if line.startswith('ATOM') and line[13:15] == 'CA':
                    try:
                        orientation = [float(lines[-step - 2][30:38]) - float(lines[-step][30:38]),
                                       float(lines[-step - 2][38:46]) - float(lines[-step][38:46]),
                                       float(lines[-step - 2][46:54]) - float(lines[-step][46:54])]
                        nme_atom_coord = self._calculate_atom_coord(lines, -step - 1, orientation, 1)
                        nme_atom_line = f'ATOM  {len(lines) + 1:5d}  N   NME A{len(self.seqpp._seq_to_list(self.seq)) + 1:4d}    {nme_atom_coord[0]:8.3f}{nme_atom_coord[1]:8.3f}{nme_atom_coord[2]:8.3f}  1.00  0.00           N  '
                        # insert NME to the last line of ATOM
                        lines.append(nme_atom_line)
                        break
                    except IndexError:
                        print("Error: Failed to calculate NME atom coordinates")
                        break
        return lines

    # def _insert_ace_and_nme(self, lines: str) -> str:
    #     '''
    #         inserts ACE and NME residues into the protein's PDB file.
    #         This function calculates the coordinates for the ACE and NME atoms based on the coordinates of the first and last residues in the protein, respectively, 
    #         and then inserts the ACE and NME atoms into the appropriate positions in the PDB file.
    #     '''
    #     # add ACE and NME to the PDB file
    #     # calculate atom coordinates of ACE according to the first residue
    #     if self._has_cap:
    #         for step, line in enumerate(lines):
    #             if line.startswith('ATOM') and line[13:15] == 'CA':
    #                 CA_coord = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    #                 oritentation = [float(lines[step+1][30:38]) - float(lines[step-1][30:38]),
    #                                 float(lines[step+1][38:46]) - float(lines[step-1][38:46]),
    #                                 float(lines[step+1][46:54]) - float(lines[step-1][46:54])]
    #                 ace_atom_coord = [CA_coord[0] - oritentation[0], CA_coord[1] - oritentation[1], CA_coord[2] - oritentation[2]]
    #                 ace_atom_line = f'ATOM  {len(lines)+1:5d}  C   ACE A   0    {ace_atom_coord[0]:8.3f}{ace_atom_coord[1]:8.3f}{ace_atom_coord[2]:8.3f}  1.00  0.00           C  '
    #                 print(ace_atom_line)
    #                 # insert ACE to the first line of ATOM
    #                 lines.insert(step-1, ace_atom_line)
    #                 break
    #     # calculate atom coordinates of NME according to the last residue
    #     if self._has_ter:
    #         for step, line in enumerate(lines[::-1]):
    #             if line.startswith('ATOM') and line[13:15] == 'CA':
    #                 CA_coord = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    #                 oritentation = [float(lines[step+1][30:38]) - float(lines[step-1][30:38]),
    #                                 float(lines[step+1][38:46]) - float(lines[step-1][38:46]),
    #                                 float(lines[step+1][46:54]) - float(lines[step-1][46:54])]
    #                 nme_atom_coord = [CA_coord[0] + oritentation[0], CA_coord[1] + oritentation[1], CA_coord[2] + oritentation[2]]
    #                 nme_atom_line = f'ATOM  {len(lines)+1:5d}  N   NME A{len(self.seqpp._seq_to_list(self.seq))+1:4d}    {nme_atom_coord[0]:8.3f}{nme_atom_coord[1]:8.3f}{nme_atom_coord[2]:8.3f}  1.00  0.00           N  '
    #                 # insert NME to the last line of ATOM
    #                 print(nme_atom_line)
    #                 lines.append(nme_atom_line)
    #                 break
    #     return lines

    def _pdb_to_seq(self, pdb_file_path):
        """
        Extract the amino acid sequence from a PDB file. If the amino acid is not standard, return 'X'.
        
        Args:
        - pdb_file_path (str): Path to the PDB file.
        
        Returns:
        - str: Amino acid sequence.
        """
        # Define the mapping of three-letter code to one-letter code for amino acids.
        amino_acid_map = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
            'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
            'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
            'TYR': 'Y', 'VAL': 'V'
        }

        sequence = ""
        last_residue_number = None

        with open(pdb_file_path, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    residue_name = line[17:20].strip()
                    residue_number = line[22:26].strip()

                    if residue_number != last_residue_number:
                        sequence += amino_acid_map.get(residue_name, 'X')
                        last_residue_number = residue_number

        return sequence

    def homology_modeling(self, seq):
        if self.template_pdb_file_path is None:
            raise ValueError('Please provide template pdb file path')

        try:
            from modeller import Environ, Model, Alignment
            from modeller.automodel import automodel
        except ImportError:
            raise ImportError('Please install modeller first')

        # self.seq -> sequence.fasta, fasta文件仅包含序列和描述信息，不包含额外的元数据
        with open(f'{self.output}/sequence.fasta', 'w') as f:
            f.write(f'>target\n{seq}')

        try:
            pwd = os.getcwd()
            shutil.copy(self.template_pdb_file_path, f'{self.output}/template.pdb')  # 复制template文件到输出目录
            os.chdir(self.output)  # 调整工作文件夹到输出文件夹

            # generate alignment file (alignment.ali)
            e = Environ()  # 创建Modeller环境e
            m = Model(e, file='template')  # 加载模板模型template.pdb并生成Model对象
            aln = Alignment(e)  # 对齐模板模型和目标序列
            aln.append_model(m, align_codes='template', atom_files='template.pdb')
            aln.append(file='sequence.fasta', align_codes='target', alignment_format='FASTA')
            aln.align2d()  # 二维序列对齐，以便将模板结构与目标序列对齐
            aln.write(file='alignment.ali', alignment_format='PIR')  # 生成对齐文件

            # homology modeling
            # 创建automodel对象a，并指定对齐文件alignment.ali、已知模板template和目标序列target。
            a = automodel(e, alnfile='alignment.ali', knowns='template', sequence='target')
            a.starting_model = 1  # 只生成一个模型
            a.ending_model = 1  # 只生成一个模型

            a.make()  # 执行建模过程
        finally:
            os.chdir(pwd)  # 执行完建模过程后，恢复到原来的工作目录

        return os.path.join(pwd, self.output, "target.B99990001.pdb")  # 返回生成的PDB文件路径

    def _seq_to_pdb(self,
                    method: str = 'alphafold',
                    max_retries: int = 3,
                    local: bool = True) -> str:
        '''
            Generate PDB file from sequence using ESMFold or homology modeling using Modeller
            DOI: 10.1101/2022.07.20.500902
            DOI: 10.1007/978-1-0716-0892-0_14
        '''

        seq_list = self.seqpp._seq_to_list(self.seq)

        # replace Non-standard amino acids with Alanine
        non_std_aa = ['B', 'Aib', 'X', 'S3', 'S5', 'S8', 'R3', 'R5', 'R8']
        if self.additional_residues is not None:
            non_std_aa.extend(list(self.additional_residues.keys()))

        # std_seq_list标准化后的氨基酸序列，其中非标准氨基酸被替换为A
        std_seq_list = [seq if seq not in non_std_aa else 'A' for seq in seq_list]
        # std_seq_list修饰基团如ACE、NME、Ac、NH2
        std_seq_list = [seq for seq in std_seq_list if seq not in ['ACE', 'NME', 'Ac', 'NH2']]

        if method == 'alphafold':
            # 如果local等于True，就用本地的ESMFold模型预测
            if local:
                from esmfold import predict_pdb
                device = 'cuda:0'
                lines = predict_pdb(''.join(std_seq_list))  # 要传一个gpu
                lines = lines.split('\n')
            else:
                url = 'https://api.esmatlas.com/foldSequence/v1/pdb/'
                for i in range(max_retries):
                    try:
                        r = requests.post(url, data=''.join(std_seq_list))
                        r.raise_for_status()  # This will raise an HTTPError if the HTTP request returned an unsuccessful status code
                        lines = r.text.split('\n')
                        break
                    except requests.exceptions.RequestException as e:
                        if i < max_retries - 1:  # If it's not the last retry
                            time.sleep(2)  # Wait for 2 seconds before retrying
                        else:  # If it's the last retry
                            raise e  # Propagate the error up, so it can be caught and handled outside of this function
        elif method == 'modeller':
            hm_file = self.homology_modeling(''.join(std_seq_list))  # 建模，返回的是pdb文件的路径
            with open(hm_file, 'r') as f:
                lines = f.readlines()  # 一个氨基酸有多个原子构成，这个lines就是原肽链里每个氨基酸原子的序列
        else:
            raise ValueError(
                'Please choose a method from alphafold or modeller, and provide template pdb if you choose modeller')

        stapled_idx_list = []
        for step, seq in enumerate(seq_list):  # 提取肽链中非标准氨基酸的位置信息
            if seq in non_std_aa:
                if self._has_cap:
                    stapled_idx_list.append(step)
                else:
                    stapled_idx_list.append(step + 1)

        seq_3_letter = self.seqpp._one_to_three(self.seq).split(' ')  # 肽链转化为三字表示
        # if additional_residues is not None:
        #     seq_3_letter = [additional_residues[seq][0] if seq in additional_residues else seq for seq in seq_3_letter]
        stapled_aa_type_list = [seq_3_letter[step] for step, seq in enumerate(seq_list) if seq in non_std_aa]  # 提取肽链中非标准氨基酸的三字表达

        new_pdb = []
        for line in lines:
            if line.startswith('ATOM'):
                if int(line[22:26]) in stapled_idx_list and line[13:15] != 'CB':  # 查看该原子对应的氨基酸编号是否在非标准氨基酸位置序列中
                    idx = int(line[22:26])
                    if idx in stapled_idx_list:
                        aa_type = stapled_aa_type_list[stapled_idx_list.index(idx)]
                        line = line[:17] + aa_type + line[20:]  # 因为之前用A代替非标准氨基酸，所以现在要还回去
                        new_pdb.append(line)
                elif int(line[22:26]) not in stapled_idx_list:
                    new_pdb.append(line)
            else:
                new_pdb.append(line)

        new_pdb = self._insert_ace_and_nme(new_pdb)  # 添加ACE和NME
        if not os.path.exists(self.output):
            os.makedirs(self.output)

        with open(f'{self.output}/model.pdb', 'w') as f:
            f.write('\n'.join(new_pdb))  # 保存
        return f'{self.output}/model.pdb'

    @property
    def _stapled_idx(self) -> list[str]:
        seq_list = self.seqpp._seq_to_list(self.seq)
        return [step + 1 for step, seq in enumerate(seq_list) if seq in ['S3', 'S5', 'S8', 'R3', 'R5', 'R8']]

    def _gen_tleap_file(self,
                        prmtop_vac: str = 'pep_vac.prmtop',
                        inpcrd_vac: str = 'pep_vac.inpcrd',
                        prmtop_sol: str = 'pep.prmtop',
                        inpcrd_sol: str = 'pep.inpcrd',
                        additional_residues: dict = None) -> str:
        '''
            Generate tleap file for AMBER

            Args:
                prmtop_vac: prmtop file name of vacuum system (default: pep_vac.prmtop)
                inpcrd_vac: inpcrd file name of vacuum system (default: pep_vac.inpcrd)
                prmtop_sol: prmtop file name of solvated system (default: pep.prmtop)
                inpcrd_sol: inpcrd file name of solvated system (default: pep.inpcrd)
                additional_residues: additional residues to be added to the system (default: None)
                    For example, {'AIB': ('/path/to/AIB.prepin', '/path/to/frcmod.AIB')
                                  'NLE': ('/path/to/NLE.prepin', '/path/to/frcmod.NLE')}

            Returns:
                tleap file name
        '''
        # 设置模拟的力场
        # ff14SB是蛋白质模型力场
        # tip3p是模型水环境下的生物分子行为
        # mbondi3是一组默认的原子半径，在模拟中更适合处理蛋白质的静电溶剂化能
        lines = ['source leaprc.protein.ff14SB',
                 'source leaprc.water.tip3p',
                 'set default PBRadii mbondi3']

        # 非标准残基
        residues = {'S3': 'PS3',
                    'S5': 'PS5',
                    'S8': 'PS8',
                    'R3': 'PR3',
                    'R5': 'PR5',
                    'R8': 'PR8',
                    'B': 'NLE',
                    'Aib': 'AIB'}
        # 加载residues中残基的AMBER力场参数文件
        for residue, prefix in residues.items():
            if residue in self.seqpp._seq_to_list(self.seq):
                prepin_file = f"{os.path.dirname(os.path.realpath(__file__))}/templates/{prefix}/{prefix.lower()}.prepin"
                frcmod_file = f"{os.path.dirname(os.path.realpath(__file__))}/templates/{prefix}/frcmod.{prefix.lower()}"
                lines.extend((f'loadAmberPrep {prepin_file}', f'loadAmberParams {frcmod_file}'))

        # 如果存在额外残基，就加载额外残基的AMBER力场参数文件
        if additional_residues is not None:
            for residue, prefix in additional_residues.items():
                residues[residue] = prefix
                additional_prepin = additional_residues[residue][0]
                additional_prepin = os.path.abspath(additional_prepin)
                additional_frcmod = additional_residues[residue][1]
                additional_frcmod = os.path.abspath(additional_frcmod)
                lines.extend((f'loadAmberPrep {additional_prepin}', f'loadAmberParams {additional_frcmod}'))

        # 判断方法，Modeller、ESMFold、AmberTools
        if self.method is None:
            lines.append('pep = sequence { ' + self.seqpp._one_to_three(self.seq) + ' }')  # 将单字氨基酸名称转化为三字
        elif self.method == 'alphafold':
            pdb_file = self._seq_to_pdb(method='alphafold', local=True)
            base_pdb_file = os.path.basename(pdb_file)
            lines.append(f'pep = loadpdb {base_pdb_file}')
        elif self.method == 'modeller':
            pdb_file = self._seq_to_pdb(method='modeller')
            base_pdb_file = os.path.basename(pdb_file)
            lines.append(f'pep = loadpdb {base_pdb_file}')
        # 判断_stapled_idx 列表的长度是否不为零，并且长度是否为偶数
        if len(self._stapled_idx) != 0 and len(self._stapled_idx) % 2 == 0:
            if self._has_cap and self.method is not None:
                # 在两个位置self._stapled_idx[i] 和 self._stapled_idx[i + 1]中创建一个C2原子的化学键
                lines.extend(f'bond pep.{self._stapled_idx[i] - 1}.C2 pep.{self._stapled_idx[i + 1] - 1}.C2' for i in
                             range(0, len(self._stapled_idx), 2))
            else:
                lines.extend(f'bond pep.{self._stapled_idx[i]}.C2 pep.{self._stapled_idx[i + 1]}.C2' for i in
                             range(0, len(self._stapled_idx), 2))

        # 创建准备分子动力学模拟所需的tleap输入文件
        lines.extend([f'saveAmberParm pep {prmtop_vac} {inpcrd_vac}',
                      'solvatebox pep TIP3PBOX 10.0',
                      'addions pep Na+ 0',
                      'addions pep Cl- 0',
                      f'saveAmberParm pep {prmtop_sol} {inpcrd_sol}',
                      'savepdb pep pep_solv.pdb',
                      'quit'])

        output_file = f'{self.output}/tleap.in'
        with open(output_file, 'w') as f:
            f.write('\n'.join(lines))
        return output_file

    def _gen_prmtop_and_inpcrd_file(self, prmtop_vac: str = 'pep_vac.prmtop', inpcrd_vac: str = 'pep_vac.inpcrd',
                                    prmtop_sol: str = 'pep.prmtop', inpcrd_sol: str = 'pep.inpcrd') -> None:
        '''
            Generate prmtop and inpcrd file

            Args:
                prmtop_vac: prmtop file name for vacuum (default: pep_vac.prmtop)
                inpcrd_vac: inpcrd file name for vacuum (default: pep_vac.inpcrd)
                prmtop_sol: prmtop file name for solvated (default: pep.prmtop)
                inpcrd_sol: inpcrd file name for solvated (default: pep.inpcrd)

            Returns:
                None
        '''
        # 获得肽链的拓扑文件prmtop和初始坐标文件inpcrd
        # prmtop_vac和inpcrd_vac表示真空环境
        # prmtop_sol和inpcrd_sol表示溶剂环境
        tleap_file = self._gen_tleap_file(
            prmtop_vac=prmtop_vac,
            inpcrd_vac=inpcrd_vac,
            prmtop_sol=prmtop_sol,
            inpcrd_sol=inpcrd_sol,
            additional_residues=self.additional_residues,

        )
        base_tleap_file = os.path.basename(tleap_file)
        cmd = ['tleap', '-f', base_tleap_file]
        subprocess.run(cmd, cwd=self.output)  # 采用tleap执行刚刚撰写的tleap.in文件


class Simulation(object):
    '''
        Molecular dynamics using OpenMM
    '''

    def __init__(self, work_path) -> None:
        self.work_path = work_path

    def setup(self, type: str = 'implicit',  # explicit or implicit
              solvent: str = 'water',  # water or membrane
              temperature: float = 300,  # Kelvin
              friction: float = 1.0,  # ps^-1
              timestep: float = 2,  # fs
              interval: int = 1000,
              nsteps: int = 5000000) -> None:
        '''
            The function has several parameters, some of which have default values. The function has the following parameters:

            args:
                type: the type of solvent model to use, either "explicit" or "implicit".
                solvent: the type of solvent to use, either "water" or "membrane".
                temperature: the temperature of the simulation, in Kelvins.
                friction: the friction coefficient, in ps^(-1).
                timestep: the time interval between iterations, in fs.
                interval: the interval at which system information is output, in ps.
                nsteps: the total time duration of the simulation or the number of iterations.
        '''
        if interval > nsteps or nsteps / interval < 100:
            raise ValueError('interval should be smaller than nsteps and nsteps / interval should be larger than 100')
        self.nsteps = nsteps

        if type == 'explicit':
            self.prmtop = os.path.join(self.work_path, 'pep.prmtop')
            self.inpcrd = os.path.join(self.work_path, 'pep.inpcrd')
            self.pep_solv = self._load_file
            self.system = self.pep_solv.createSystem(nonbondedMethod=app.PME,
                                                     nonbondedCutoff=8.0 * u.angstroms,
                                                     constraints=app.HBonds)
        elif type == 'implicit':
            self.prmtop = os.path.join(self.work_path, 'pep_vac.prmtop')
            self.inpcrd = os.path.join(self.work_path, 'pep_vac.inpcrd')
            self.pep_solv = self._load_file
            if solvent == 'water':
                self.system = self.pep_solv.createSystem(implicitSolvent=app.OBC2,
                                                         soluteDielectric=1.0,
                                                         solventDielectric=78.5)
            elif solvent == 'chloroform':  # Membrane solvent
                self.system = self.pep_solv.createSystem(implicitSolvent=app.OBC2,
                                                         soluteDielectric=1.0,
                                                         solventDielectric=4.7)
            elif solvent == 'DMF':
                self.system = self.pep_solv.createSystem(implicitSolvent=app.OBC2,
                                                         soluteDielectric=1.0,
                                                         solventDielectric=36.7)
            elif solvent == 'DMSO':
                self.system = self.pep_solv.createSystem(implicitSolvent=app.OBC2,
                                                         soluteDielectric=1.0,
                                                         solventDielectric=48.9)
            elif solvent == 'ethanol':
                self.system = self.pep_solv.createSystem(implicitSolvent=app.OBC2,
                                                         soluteDielectric=1.0,
                                                         solventDielectric=24.4)
            elif solvent == 'acetone':
                self.system = self.pep_solv.createSystem(implicitSolvent=app.OBC2,
                                                         soluteDielectric=1.0,
                                                         solventDielectric=20.7)
            else:
                raise ValueError(
                    'Solvent not supported, please choose from water, chloroform, DMF, DMSO, ethanol, acetone')
        else:
            raise ValueError('Type of simulation must be explicit or implicit')

        self.integrator = mm.LangevinIntegrator(
            temperature * u.kelvin,  # Temperature of heat bath
            friction / u.picoseconds,  # Friction coefficient
            timestep * u.femtoseconds)  # Time step
        self.platform = mm.Platform.getPlatformByName('CUDA')
        prop = dict(CudaPrecision='mixed')  # Use mixed single/double precision
        self.sim = app.Simulation(self.pep_solv.topology, self.system, self.integrator, self.platform, prop)
        self.sim.context.setPositions(self.pep_solv.positions)
        self.sim.reporters.append(
            mm.app.DCDReporter(f'{self.work_path}/traj.dcd', interval)
        )
        self.sim.reporters.append(
            pmd.openmm.StateDataReporter(f'{self.work_path}/traj.log', reportInterval=interval,
                                         volume=True, density=True, separator='\t')
        )
        self.sim.reporters.append(
            pmd.openmm.ProgressReporter(f'{self.work_path}/traj.log.info', interval * 5, self.nsteps)
        )

    def minimize(self):
        print('Minimizing...')
        self.sim.minimizeEnergy()

    def run(self):
        start_time = time.time()
        print('Running...')
        self.sim.step(self.nsteps)
        end_time = time.time()
        print(f'Done! Time elapsed: {end_time - start_time} s')

    @property
    def _load_file(self) -> pmd.Structure:
        return pmd.load_file(self.prmtop, self.inpcrd, structure=True)


if __name__ == '__main__':
    # Initialize a sequence and output directory
    # 序列规则:
    #   使用氨基酸单字母缩写用
    #   单字母氨基酸之间可以连续写，也可以用'-'分隔，非单字母氨基酸之间必须用'-'分隔。例如: Ac-BATPR8RRR-Aib-LBRR3FKRLQ-NH2
    #   N端可以使用'Ac-'开头, 添加Acetyl group
    #   C端可以使用'-NH2'结尾, 添加Amide group
    #   Aib是2-Aminoisobutyric acid
    #   B是Nle的缩写，Nle是Norleucine
    seq = 'Ac-BATPR8RRR-Aib-B-LBRR3FKRLQ'
    output = 'data'

    # Prepare the protein for molecular dynamics simulation
    pp = PrepareProt(seq, output, method='alphafold')
    pp._gen_prmtop_and_inpcrd_file()
    sim = Simulation(output)
    sim.setup(
        type='implicit',  # 'explicit' or 'implicit'
        solvent='water',  # 'water' or 'chloroform'
        temperature=300,  # Kelvin
        friction=1,  # ps^-1
        timestep=2,  # fs
        interval=100,  # ps
        nsteps=500000  # 1 ns
    )
    sim.minimize()

    # seq = 'Ac-RRRRR-NH2'
    # output = 'data'
    # pp = PrepareProt(seq, output, method='alphafold')
    # pp._gen_prmtop_and_inpcrd_file()
    #
    # sim = Simulation(output)
    # sim.setup(interval=100, nsteps=50000)
    # sim.minimize()
    # sim.run()
