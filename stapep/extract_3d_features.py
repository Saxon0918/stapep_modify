from stapep.utils import PhysicochemicalPredictor
import os

seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ'  # Define the peptide sequence
pathname = '/tmp/60f9a3d2-a7ce-401f-a590-0996b40a6540'  # The path where the topology file and trajectory file are located
pcp = PhysicochemicalPredictor(sequence=seq,
                               topology_file=os.path.join(pathname, 'pep_vac.prmtop'),
                               # topology file　(default: pep_vac.prmtop in the data folder)
                               trajectory_file=os.path.join(pathname, 'traj.dcd'),
                               # trajectory file (default: traj.dcd in the data folder)
                               start_frame=0)  # start frame (default: 500)

# Get the features
print('helix percent: ', pcp.calc_helix_percent())  # 计算α螺旋占比
print('sheet percent: ', pcp.calc_extend_percent())  # 计算β折叠占比
print('loop percent: ', pcp.calc_loop_percent())  # 计算环结构（loop）占比

# save the mean structure of the trajectory
pcp._save_mean_structure(os.path.join(pathname, 'mean_structure.pdb'))  # 保存分子的平均结构的pdb

# calculate the Mean B-factor, Molecular Surface, Mean Gyration Radius, Hydrophobic Index, and 3D-PSA
print('mean bfactor: ', pcp.calc_mean_bfactor())  # 平均B因子，较低的B因子表明原子较稳定，柔性较小
print('mol surf: ', pcp.calc_mean_molsurf())  # 分子表面积
print('mean gyrate: ', pcp.calc_mean_gyrate())  # 分子的旋转半径，分子质量分布的中心到外部的平均距离，评估分子的形状和紧密性
# print('hydrophobic index: ', pcp.calc_hydrophobic_index(os.path.join(pathname, 'mean_structure.pdb')))  # TODO 疏水指数环境没配好
print('psa: ', pcp.calc_psa(os.path.join(pathname, 'mean_structure.pdb')))  # 极性表面积，评价分子吸收、渗透和跨膜运输能力
print('total number of hydrogen bonds: ', pcp.calc_n_hbonds())  # 氢键的总数，决定了分子的稳定性和结合能力

# extract 2D structure of the peptide
smiles = pcp.extract_2d_structure(os.path.join(pathname, 'mean_structure.pdb'))
print(smiles)
