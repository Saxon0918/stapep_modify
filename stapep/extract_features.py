from stapep.utils import ProtParamsSeq
import os

# seq = 'Ac-BATP-R8-RRR-Aib-BLBR-R3-FKRLQ'  # Define the peptide sequence
seq = 'EEEKEKK-R8-KKEEEEKK-R3-QKEKT'
pps = ProtParamsSeq(seq)  # Initialize the class

# Get the features
print('length: ', pps.seq_length)  # 序列长度（不包含N端和C端）
print('weight: ', pps.weight)  # 分子量，穿膜（越小越好）
print('hydrophobicity index: ', pps.hydrophobicity_index)  # TODO 疏水性指数（不准确）
print('charge', pps.calc_charge(pH=7.0))  # pH=7.0时的净电荷
print('charge_density', pps.calc_charge_density(pH=7.0))  # pH=7.0时的电荷密度
print('aromaticity', pps.aromaticity)  # 香族氨基酸比例
print('fraction_arginine', pps.fraction_arginine)  # 精氨酸占比
print('fraction_lysine', pps.fraction_lysine)  # 赖氨酸占比
print('lyticity index: ', pps.lyticity_index)  # 溶血性指数
print('isoelectric_point: ', pps.isoelectric_point)  # 等电点

pathname = 'data'
pps.plot_lyticity_index(os.path.join(pathname, 'lyticity_index.svg'))  # 间隔3个或4个位点的溶解性指数图像，左上角溶血性结果
