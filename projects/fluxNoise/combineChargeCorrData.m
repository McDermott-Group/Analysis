% base_path = 'Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\';
% files = { ...
%     'Q3Q4Corr\General\Parameter\csn0100hzy_correlation.hdf5', ...
%     'Q3Q4Corr\General\Parameter\csn1803vvg_correlation.hdf5', ...
%     'Q3Q4Corr\General\Parameter\cso0030jjg_correlation.hdf5', ...
%     'Q3Q4Corr\General\Parameter\csq0128zge_correlation.hdf5', ...
%     }
% qA = 'Q3'; qB = 'Q4'; fignum = 1111; 

% base_path = 'Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\';
% files = { ...
%     'Q1Q2Corr\General\Parameter\csu0052mrt_correlation.hdf5', ...
%     'Q1Q2Corr\General\Parameter\csv0034shs_correlation.hdf5', ...
%     }
% qA = 'Q1'; qB = 'Q2'; fignum = 1112; 

%%%% Actually Qubits 1,3
base_path = 'Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-11-12\CorrFar\';
files = { ...
    'Q4\General\Parameter\crc0922kvk_parameters.hdf5', ...
    'Q4\General\Parameter\crc2345hmc_parameters.hdf5', ...
    }
qA = 'Q3'; qB = 'Q4'; fignum = 1113; 

N_all = [];
nA_all = [];
nB_all = [];
nCorrAB_all = [];
nCorrBA_all = [];
nCorrExpected_all = [];
for f = files
    [N, nA, nB, nCorrAB, nCorrBA, nCorrExpected] = AnalyzeChargeJumpCorrelations([base_path f{1}], qA, qB, fignum);
    N_all = [N_all N];
    nA_all = [nA_all nA];
    nB_all = [nB_all nB];
    nCorrAB_all = [nCorrAB_all nCorrAB];
    nCorrBA_all = [nCorrBA_all nCorrBA];
    nCorrExpected_all = [nCorrExpected_all nCorrExpected];
end
N_tot = sum(N_all)
nA_tot = sum(nA_all)
nB_tot = sum(nB_all)
nCorrAB_tot = sum(nCorrAB_all)
nCorrBA_tot = sum(nCorrBA_all)
nCorrExpected_tot = sum(nCorrExpected_all)