import numpy as np
import matplotlib.pyplot as plt
import noiselib
from QPTunneling import *
from ChargeOffset import *


base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Q3Q4Corr\General\Parameter\\'
CO = ChargeOffset()
CO.add_dataset(base_path + 'cvc0232imo_correlation.hdf5').limit_dataset('Q2')
CO.add_dataset(base_path + 'cvd0430bxy_correlation.hdf5').limit_dataset('Q4')
CO.add_dataset(base_path + 'cvd0802ltu_correlation.hdf5').limit_dataset('Q4')
CO.add_dataset(base_path + 'cvd1745asc_correlation.hdf5')
CO.add_dataset(base_path + 'cve0104ppt_correlation.hdf5')
CO.add_dataset(base_path + 'cvf0542fei_correlation.hdf5')
base_path = 'fluxNoise2\DR1 - 2019-12-17\CorrFar\Q1Q2Q3Q4Corr\General\Parameter\\'
CO.add_dataset(base_path + 'cvj0203nun_correlation.hdf5')
CO.add_dataset(base_path + 'cvj1806ywn_correlation.hdf5')
CO.add_dataset(base_path + 'cvk0412mns_correlation.hdf5')
CO.plot_charge_offset()
CO.plot_jump_sizes()
CO.get_jump_sizes(plot=True)
CO.plot_charge_correlation('Q1','Q2')
CO.plot_charge_correlation('Q3','Q4')
CO.plot_charge_correlation('Q1','Q3')
CO.plot_time_steps()
plt.show()

"""######### Interleaved Data ##########"""

# base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q3Q4Corr\General\Parameter\\'
# CO = ChargeOffset()
# CO.add_dataset(base_path + 'csn0100hzy_correlation.hdf5')
# CO.add_dataset(base_path + 'csn1803vvg_correlation.hdf5')
# CO.add_dataset(base_path + 'cso0030jjg_correlation.hdf5')
# CO.add_dataset(base_path + 'csq0128zge_correlation.hdf5')
# # CO.add_dataset(base_path + 'csy2025oeg_correlation.hdf5')
# # CO.add_dataset(base_path + 'csz0145znf_correlation.hdf5')
# CO.plot_charge_offset()
# CO.plot_jump_sizes()
# CO.plot_charge_correlation('Q3','Q4')
# CO.plot_time_steps()
# plt.show()

# base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Corr\General\Parameter\\'
# CO = ChargeOffset()
# CO.add_dataset(base_path + 'csu0052mrt_correlation.hdf5')
# CO.add_dataset(base_path + 'csv0034shs_correlation.hdf5')
# CO.add_dataset(base_path + 'csv0824hxv_correlation.hdf5')
# # CO.add_dataset(base_path + 'cvh0508hpk_correlation.hdf5') # dif timing, etc?  good data though
# CO.plot_charge_offset()
# CO.plot_jump_sizes()
# CO.plot_charge_correlation('Q1','Q2')
# CO.plot_time_steps()
# plt.show()

# base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q3Corr\General\Parameter\\'
# CO = ChargeOffset()
# CO.add_dataset(base_path + 'csw0353opr_correlation.hdf5')
# CO.add_dataset(base_path + 'csx0545cpf_correlation.hdf5')
# CO.add_dataset(base_path + 'csx1710hcz_correlation.hdf5')
# CO.add_dataset(base_path + 'csy0726yvn_correlation.hdf5')
# # base_path = 'fluxNoise\DR1 - 2019-12-17\CorrFar\Q4\General\Parameter\\'
# # CO.add_dataset(base_path + 'crc2345hmc_parameters.hdf5') # might say Q3,Q4 instead?
# CO.plot_charge_offset()
# CO.plot_jump_sizes()
# CO.plot_charge_correlation('Q1','Q3')
# CO.plot_time_steps()
# plt.show()