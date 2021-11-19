"""
PSD for Q1, Q2, Q4 at different J2 Bias
Data:
Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJRadiatorQPT_2021Aug19
Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJRadiatorQPT_2021Aug19_HighDensity
Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJJRadiatorQPT_HighFreq_2021Aug27
Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJJRadiatorQPT_LowFreq_2021Aug27
Fitting method Chris' no white noise verison

"""
import noiselib
import matplotlib.pyplot as plt
import numpy as np

Q1= np.array([[0.00000000e+00, 2.06462838e+03],
       [2.50000000e+00, 2.14159974e+03],
       [5.00000000e+00, 2.11753776e+03],
       [7.50000000e+00, 2.10100522e+03],
       [1.00000000e+01, 1.98788372e+03],
       [1.25000000e+01, 1.80400137e+03],
       [1.50000000e+01, 1.99009825e+03],
       [1.58000000e+01, 2.15122391e+03],
       [1.63000000e+01, 2.21216568e+03],
       [1.68000000e+01, 2.21065742e+03],
       [1.73000000e+01, 2.50545713e+03],
       [1.78000000e+01, 2.47681769e+03],
       [1.83000000e+01, 2.46035263e+03],
       [1.88000000e+01, 2.54705682e+03],
       [1.93000000e+01, 2.24051649e+03],
       [1.98000000e+01, 2.29891743e+03],
       [2.03000000e+01, 2.37333672e+03],
       [2.08000000e+01, 2.40258296e+03],
       [2.13000000e+01, 2.40969300e+03],
       [2.18000000e+01, 2.43144508e+03],
       [2.23000000e+01, 2.45468133e+03],
       [2.28000000e+01, 3.34289962e+03],
       [2.33000000e+01, 3.80406426e+03],
       [2.38000000e+01, 4.92192346e+03],
       [2.43000000e+01, 4.85653587e+03],
       [2.48000000e+01, 6.05877427e+03],
       [2.53000000e+01, 5.15048320e+03],
       [2.58000000e+01, 5.16333980e+03],
       [2.63000000e+01, 4.49718396e+03],
       [2.68000000e+01, 3.84549917e+04],
       [2.73000000e+01, 7.51188202e+04],
       [2.78000000e+01, 7.56242052e+04],
       [2.83000000e+01, 6.08098331e+04],
       [2.88000000e+01, 5.24566938e+04],
       [2.93000000e+01, 9.82375524e+03],
       [2.98000000e+01, 5.27141580e+03],
       [3.03000000e+01, 1.09139357e+04],
       [3.08000000e+01, 9.00209034e+03],
       [3.13000000e+01, 1.22792841e+04],
       [3.18000000e+01, 7.59313029e+03],
       [3.23000000e+01, 1.03385242e+04],
       [3.28000000e+01, 5.99419972e+03],
       [3.33000000e+01, 5.92398810e+03],
       [3.38000000e+01, 5.23191368e+03],
       [3.43000000e+01, 4.95831460e+03],
       [3.48000000e+01, 5.07252232e+03],
       [3.53000000e+01, 5.59550751e+03],
       [3.58000000e+01, 5.76881276e+03],
       [3.63000000e+01, 5.60744165e+03],
       [3.68000000e+01, 4.87476449e+03],
       [3.73000000e+01, 4.87167986e+03],
       [3.78000000e+01, 4.43962268e+03],
       [3.83000000e+01, 4.39982964e+03],
       [3.88000000e+01, 4.84167021e+03],
       [3.93000000e+01, 5.45804007e+03],
       [3.98000000e+01, 5.27754405e+03],
       [4.03000000e+01, 6.46861330e+03],
       [4.08000000e+01, 6.25228974e+03],
       [4.13000000e+01, 4.83929048e+03],
       [4.18000000e+01, 4.50698046e+03],
       [4.23000000e+01, 5.29617862e+03],
       [4.28000000e+01, 5.65108251e+03],
       [4.33000000e+01, 5.64339140e+03],
       [4.38000000e+01, 5.63787385e+03],
       [4.43000000e+01, 5.93678317e+03],
       [4.48000000e+01, 6.31534429e+03],
       [4.53000000e+01, 6.37649422e+03],
       [4.58000000e+01, 7.19105199e+03],
       [4.63000000e+01, 6.77868091e+03],
       [4.68000000e+01, 6.30307753e+03],
       [4.73000000e+01, 6.62688704e+03],
       [4.78000000e+01, 5.41197339e+03],
       [4.83000000e+01, 4.75156843e+03],
       [4.88000000e+01, 5.70873142e+03],
       [4.93000000e+01, 5.88193636e+03],
       [4.98000000e+01, 5.65812037e+03],
       [5.03000000e+01, 5.33482048e+03],
       [5.08000000e+01, 5.47988814e+03],
       [5.13000000e+01, 5.57163521e+03],
       [5.18000000e+01, 5.41144626e+03],
       [5.23000000e+01, 7.58924423e+03],
       [5.28000000e+01, 7.45164140e+03],
       [5.33000000e+01, 6.62803166e+03],
       [5.38000000e+01, 6.23525335e+03],
       [5.43000000e+01, 6.60989820e+03],
       [5.48000000e+01, 6.19348417e+03],
       [5.53000000e+01, 5.55917175e+03],
       [5.58000000e+01, 4.77236675e+03],
       [5.63000000e+01, 5.69461469e+03],
       [5.68000000e+01, 5.66342058e+03],
       [5.73000000e+01, 5.81542345e+03],
       [5.78000000e+01, 5.10132664e+03],
       [5.83000000e+01, 4.69881189e+03],
       [5.88000000e+01, 4.76791522e+03],
       [5.93000000e+01, 4.93757574e+03],
       [5.98000000e+01, 5.06820973e+03],
       [6.03000000e+01, 4.85394983e+03],
       [6.08000000e+01, 4.63741567e+03],
       [6.50000000e+01, 3.94578719e+03],
       [7.00000000e+01, 3.77087635e+03],
       [7.50000000e+01, 3.74657434e+03],
       [8.00000000e+01, 3.87538282e+03],
       [8.50000000e+01, 3.67523525e+03],
       [9.00000000e+01, 3.93150625e+03],
       [9.50000000e+01, 4.06390394e+03],
       [1.00000000e+02, 7.47805672e+03],
       [1.05000000e+02, 1.75656775e+04],
       [1.10000000e+02, 4.15690840e+04],
       [1.15000000e+02, 6.14916579e+04],
       [1.20000000e+02, 6.04498428e+04],
       [1.25000000e+02, 6.11725234e+04],
       [1.30000000e+02, 7.43910557e+04],
       [1.35000000e+02, 4.53429740e+04],
       [1.40000000e+02, 6.90157039e+03],
       [1.45000000e+02, 6.05090167e+03]])

Q2=np.array([[0.00000000e+00, 5.25033283e+01],
       [2.50000000e+00, 5.32239100e+01],
       [5.00000000e+00, 5.23701869e+01],
       [7.50000000e+00, 5.15186194e+01],
       [1.00000000e+01, 5.21718672e+01],
       [1.25000000e+01, 5.32478620e+01],
       [1.50000000e+01, 5.27486062e+01],
       [1.58000000e+01, 5.34666665e+01],
       [1.63000000e+01, 5.23874397e+01],
       [1.68000000e+01, 5.16505189e+01],
       [1.73000000e+01, 5.18825975e+01],
       [1.78000000e+01, 5.41630186e+01],
       [1.83000000e+01, 5.33239319e+01],
       [1.88000000e+01, 5.53976725e+01],
       [1.93000000e+01, 5.31118443e+01],
       [1.98000000e+01, 5.32724007e+01],
       [2.03000000e+01, 5.17795943e+01],
       [2.08000000e+01, 5.39159321e+01],
       [2.13000000e+01, 5.15293093e+01],
       [2.18000000e+01, 5.29990198e+01],
       [2.23000000e+01, 5.46373044e+01],
       [2.28000000e+01, 5.27560980e+01],
       [2.33000000e+01, 5.62433048e+01],
       [2.38000000e+01, 5.73668340e+01],
       [2.43000000e+01, 5.99317727e+01],
       [2.48000000e+01, 6.93308686e+01],
       [2.53000000e+01, 8.59803249e+01],
       [2.58000000e+01, 1.17759073e+02],
       [2.63000000e+01, 1.55629110e+02],
       [2.68000000e+01, 1.89997122e+02],
       [2.73000000e+01, 2.32974995e+02],
       [2.78000000e+01, 2.79884282e+02],
       [2.83000000e+01, 3.06409045e+02],
       [2.88000000e+01, 3.05147956e+02],
       [2.93000000e+01, 3.12950593e+02],
       [2.98000000e+01, 3.00964079e+02],
       [3.03000000e+01, 3.03698234e+02],
       [3.08000000e+01, 2.85805324e+02],
       [3.13000000e+01, 2.75768551e+02],
       [3.18000000e+01, 2.84041993e+02],
       [3.23000000e+01, 2.89325260e+02],
       [3.28000000e+01, 2.70370588e+02],
       [3.33000000e+01, 2.31346977e+02],
       [3.38000000e+01, 2.05289456e+02],
       [3.43000000e+01, 2.02864368e+02],
       [3.48000000e+01, 2.13193661e+02],
       [3.53000000e+01, 2.22571263e+02],
       [3.58000000e+01, 2.21849115e+02],
       [3.63000000e+01, 2.20752430e+02],
       [3.68000000e+01, 2.30130299e+02],
       [3.73000000e+01, 2.49345238e+02],
       [3.78000000e+01, 2.86931400e+02],
       [3.83000000e+01, 3.18232441e+02],
       [3.88000000e+01, 3.45638789e+02],
       [3.93000000e+01, 3.88986018e+02],
       [3.98000000e+01, 4.63047101e+02],
       [4.03000000e+01, 4.74786316e+02],
       [4.08000000e+01, 4.89603938e+02],
       [4.13000000e+01, 4.52175916e+02],
       [4.18000000e+01, 4.47308320e+02],
       [4.23000000e+01, 4.69386903e+02],
       [4.28000000e+01, 4.96506036e+02],
       [4.33000000e+01, 4.88940486e+02],
       [4.38000000e+01, 4.66073918e+02],
       [4.43000000e+01, 4.54863684e+02],
       [4.48000000e+01, 5.01847251e+02],
       [4.53000000e+01, 5.58350250e+02],
       [4.58000000e+01, 5.79308331e+02],
       [4.63000000e+01, 6.03753769e+02],
       [4.68000000e+01, 6.14033809e+02],
       [4.73000000e+01, 2.64011418e+02],
       [4.78000000e+01, 2.64083147e+02],
       [4.83000000e+01, 4.42874293e+02],
       [4.88000000e+01, 4.64215581e+02],
       [4.93000000e+01, 4.48171584e+02],
       [4.98000000e+01, 4.54481595e+02],
       [5.03000000e+01, 5.06066945e+02],
       [5.08000000e+01, 5.61552412e+02],
       [5.13000000e+01, 5.60861110e+02],
       [5.18000000e+01, 5.44650493e+02],
       [5.23000000e+01, 5.38270395e+02],
       [5.28000000e+01, 5.31621733e+02],
       [5.33000000e+01, 5.19806371e+02],
       [5.38000000e+01, 4.72130348e+02],
       [5.43000000e+01, 4.49895726e+02],
       [5.48000000e+01, 4.12777162e+02],
       [5.53000000e+01, 4.00954311e+02],
       [5.58000000e+01, 4.09294659e+02],
       [5.63000000e+01, 4.05860796e+02],
       [5.68000000e+01, 3.74934484e+02],
       [5.73000000e+01, 3.44337119e+02],
       [5.78000000e+01, 3.33962979e+02],
       [5.83000000e+01, 3.82462274e+02],
       [5.88000000e+01, 4.20585455e+02],
       [5.93000000e+01, 4.42554059e+02],
       [5.98000000e+01, 4.17756476e+02],
       [6.03000000e+01, 3.70763802e+02],
       [6.08000000e+01, 3.46505067e+02],
       [6.58000000e+01, 3.19662191e+02],
       [6.78000000e+01, 3.54386958e+02],
       [6.98000000e+01, 5.00855011e+02],
       [7.18000000e+01, 5.48816725e+02],
       [7.38000000e+01, 5.51485770e+02],
       [7.58000000e+01, 7.14622934e+02],
       [7.78000000e+01, 7.87194752e+02],
       [7.98000000e+01, 9.34588071e+02],
       [8.18000000e+01, 1.05074044e+03],
       [8.38000000e+01, 1.10293067e+03],
       [8.58000000e+01, 9.55914721e+02],
       [8.78000000e+01, 9.01683697e+02],
       [8.98000000e+01, 7.42965694e+02],
       [9.18000000e+01, 6.08392470e+02],
       [9.38000000e+01, 4.16777440e+02],
       [9.58000000e+01, 3.64175788e+02],
       [9.78000000e+01, 3.47333610e+02],
       [9.98000000e+01, 3.56323278e+02],
       [1.01800000e+02, 2.98776259e+02],
       [1.03800000e+02, 3.18277636e+02],
       [1.05800000e+02, 3.21073982e+02],
       [1.07800000e+02, 3.16479103e+02],
       [1.09800000e+02, 3.28225885e+02],
       [1.11800000e+02, 3.37329038e+02],
       [1.13800000e+02, 3.55749416e+02],
       [1.15800000e+02, 3.75278583e+02],
       [1.17800000e+02, 3.86217527e+02],
       [1.19800000e+02, 4.20168101e+02],
       [1.21800000e+02, 4.35482869e+02],
       [1.23800000e+02, 4.53672837e+02],
       [1.25800000e+02, 4.89414783e+02],
       [1.27800000e+02, 5.22940061e+02],
       [1.29800000e+02, 5.61401473e+02],
       [1.31800000e+02, 5.98717570e+02],
       [1.33800000e+02, 6.36781083e+02],
       [1.35800000e+02, 6.80247511e+02],
       [1.37800000e+02, 7.16320910e+02],
       [1.39800000e+02, 7.66307504e+02],
       [1.41800000e+02, 7.77022385e+02],
       [1.43800000e+02, 8.24935997e+02],
       [1.45800000e+02, 8.71800120e+02],
       [1.47800000e+02, 9.01106905e+02],
       [1.49800000e+02, 9.58130610e+02],
       [1.51800000e+02, 9.77608362e+02],
       [1.53800000e+02, 1.03460750e+03],
       [1.55800000e+02, 1.05634739e+03],
       [1.57800000e+02, 1.12596475e+03],
       [1.59800000e+02, 1.15106914e+03],
       [1.61800000e+02, 1.18782283e+03],
       [1.63800000e+02, 1.23711790e+03],
       [1.65800000e+02, 1.22224732e+03],
       [1.67800000e+02, 1.28852980e+03],
       [1.69800000e+02, 1.30047586e+03],
       [1.71800000e+02, 1.32774017e+03],
       [1.73800000e+02, 1.33118019e+03],
       [1.75800000e+02, 1.37803139e+03],
       [1.77800000e+02, 1.36191417e+03],
       [1.79800000e+02, 1.44001442e+03],
       [1.81800000e+02, 1.42041189e+03],
       [1.83800000e+02, 1.46665350e+03],
       [1.85800000e+02, 1.47481982e+03],
       [1.87800000e+02, 1.49619041e+03],
       [1.89800000e+02, 1.46985936e+03],
       [1.91800000e+02, 1.50474805e+03],
       [1.93800000e+02, 1.48747104e+03],
       [1.95800000e+02, 1.54787253e+03],
       [1.97800000e+02, 1.58129608e+03],
       [1.99800000e+02, 1.66252235e+03],
       [2.01800000e+02, 1.55616960e+03],
       [2.03800000e+02, 1.57564537e+03],
       [2.05800000e+02, 1.49296415e+03],
       [2.07800000e+02, 1.81836675e+03],
       [2.09800000e+02, 1.56200964e+03],
       [2.11800000e+02, 1.49696424e+03],
       [2.13800000e+02, 1.18047469e+04],
       [2.15800000e+02, 1.74242504e+04],
       [2.17800000e+02, 1.24470309e+04],
       [2.19800000e+02, 2.17047325e+04],
       [2.21800000e+02, 2.37847874e+04],
       [2.23800000e+02, 2.91236915e+04],
       [2.25800000e+02, 3.87700186e+04],
       [2.27800000e+02, 4.44672146e+04],
       [2.29800000e+02, 5.02109214e+04]])

Q4=np.array([[0.00000000e+00, 6.19621125e+02],
       [2.50000000e+00, 6.18274241e+02],
       [5.00000000e+00, 6.05067302e+02],
       [7.50000000e+00, 6.17025068e+02],
       [1.00000000e+01, 6.10066088e+02],
       [1.25000000e+01, 5.99012014e+02],
       [1.50000000e+01, 6.07573224e+02],
       [1.58000000e+01, 6.01630220e+02],
       [1.63000000e+01, 6.05795325e+02],
       [1.68000000e+01, 6.04401373e+02],
       [1.73000000e+01, 6.07696451e+02],
       [1.78000000e+01, 6.16369237e+02],
       [1.83000000e+01, 6.10566213e+02],
       [1.88000000e+01, 6.01113099e+02],
       [1.93000000e+01, 6.09382604e+02],
       [1.98000000e+01, 6.03682193e+02],
       [2.03000000e+01, 6.15969087e+02],
       [2.08000000e+01, 6.24182031e+02],
       [2.13000000e+01, 6.20333024e+02],
       [2.18000000e+01, 6.35223408e+02],
       [2.23000000e+01, 6.18014933e+02],
       [2.28000000e+01, 6.06272577e+02],
       [2.33000000e+01, 6.09520583e+02],
       [2.38000000e+01, 6.24782523e+02],
       [2.43000000e+01, 6.41923912e+02],
       [2.48000000e+01, 7.57073173e+02],
       [2.53000000e+01, 8.78882433e+02],
       [2.58000000e+01, 1.10999430e+03],
       [2.63000000e+01, 1.31140624e+03],
       [2.68000000e+01, 1.41728093e+03],
       [2.73000000e+01, 1.55280768e+03],
       [2.78000000e+01, 1.60209589e+03],
       [2.83000000e+01, 1.66978276e+03],
       [2.88000000e+01, 1.64776202e+03],
       [2.93000000e+01, 1.64229384e+03],
       [2.98000000e+01, 1.60428512e+03],
       [3.03000000e+01, 1.66748388e+03],
       [3.08000000e+01, 2.02131630e+03],
       [3.13000000e+01, 2.19074465e+03],
       [3.18000000e+01, 2.58224472e+03],
       [3.23000000e+01, 2.58961292e+03],
       [3.28000000e+01, 2.22349999e+03],
       [3.33000000e+01, 2.08560597e+03],
       [3.38000000e+01, 2.05107088e+03],
       [3.43000000e+01, 2.35452697e+03],
       [3.48000000e+01, 2.59895003e+03],
       [3.53000000e+01, 3.14188196e+03],
       [3.58000000e+01, 3.26540746e+03],
       [3.63000000e+01, 3.19510338e+03],
       [3.68000000e+01, 2.91563318e+03],
       [3.73000000e+01, 2.61729702e+03],
       [3.78000000e+01, 2.57465240e+03],
       [3.83000000e+01, 2.78141413e+03],
       [3.88000000e+01, 3.24949159e+03],
       [3.93000000e+01, 3.53303992e+03],
       [3.98000000e+01, 4.41068788e+03],
       [4.03000000e+01, 4.39874811e+03],
       [4.08000000e+01, 4.94909481e+03],
       [4.13000000e+01, 5.06813390e+03],
       [4.18000000e+01, 5.06665587e+03],
       [4.23000000e+01, 5.31059320e+03],
       [4.28000000e+01, 5.72189726e+03],
       [4.33000000e+01, 5.92181548e+03],
       [4.38000000e+01, 6.10365351e+03],
       [4.43000000e+01, 6.20216977e+03],
       [4.48000000e+01, 6.60518781e+03],
       [4.53000000e+01, 6.23485601e+03],
       [4.58000000e+01, 5.96213192e+03],
       [4.63000000e+01, 5.95363344e+03],
       [4.68000000e+01, 6.16846053e+03],
       [4.73000000e+01, 6.31893439e+03],
       [4.78000000e+01, 6.51092962e+03],
       [4.83000000e+01, 5.88985595e+03],
       [4.88000000e+01, 6.01646570e+03],
       [4.93000000e+01, 6.29014582e+03],
       [4.98000000e+01, 6.26660384e+03],
       [5.03000000e+01, 6.37042678e+03],
       [5.08000000e+01, 6.09907466e+03],
       [5.13000000e+01, 6.17898532e+03],
       [5.18000000e+01, 6.55254613e+03],
       [5.23000000e+01, 7.55831677e+03],
       [5.28000000e+01, 8.10495215e+03],
       [5.33000000e+01, 8.76637128e+03],
       [5.38000000e+01, 7.72228775e+03],
       [5.43000000e+01, 7.70992601e+03],
       [5.48000000e+01, 7.98077428e+03],
       [5.53000000e+01, 7.85255996e+03],
       [5.58000000e+01, 8.30783558e+03],
       [5.63000000e+01, 7.41139531e+03],
       [5.68000000e+01, 6.20519526e+03],
       [5.73000000e+01, 6.64101397e+03],
       [5.78000000e+01, 5.86206074e+03],
       [5.83000000e+01, 5.66778844e+03],
       [5.88000000e+01, 5.61098192e+03],
       [5.93000000e+01, 5.49788359e+03],
       [5.98000000e+01, 5.49560704e+03],
       [6.03000000e+01, 5.36958182e+03],
       [6.08000000e+01, 5.03707967e+03],
       [6.50000000e+01, 3.66047715e+03],
       [7.00000000e+01, 2.75304075e+03],
       [7.50000000e+01, 2.61791197e+03],
       [8.00000000e+01, 2.08231875e+03],
       [8.50000000e+01, 2.04904789e+03],
       [9.00000000e+01, 2.15759340e+03],
       [9.50000000e+01, 2.26049619e+03],
       [1.00000000e+02, 2.46073525e+03],
       [1.05000000e+02, 2.80912303e+03],
       [1.10000000e+02, 3.07963685e+03],
       [1.15000000e+02, 3.41498104e+03],
       [1.20000000e+02, 3.74304398e+03],
       [1.25000000e+02, 4.10846496e+03],
       [1.30000000e+02, 4.71665608e+03],
       [1.35000000e+02, 5.46833319e+03],
       [1.40000000e+02, 6.00446939e+03],
       [1.45000000e+02, 6.17067560e+03]])



# offset = 19.2
# print('Q1[:, 0]=', Q1[:, 0])
# print('Q2[:, 0]=', Q2[:, 0])
# print('Q4[:, 0]=', Q4[:, 0])
# plt.plot(Q1[:, 0], Q1[:, 1], color='b', label='Q1')
# plt.plot(Q2[:, 0], Q2[:, 1], color='r', label='Q2')
# plt.plot(Q4[:, 0], Q4[:, 1], color='y', label='Q4')
# plt.xlabel('Radiator Josephson Frequency (mDAC)')
# plt.ylabel('PSD')
# plt.yscale('log')
# plt.grid()
# plt.legend(loc=1)
# plt.show()

f = 4.604
plt.plot(Q1[:, 0]*f, Q1[:, 1], color='b', label='Q1')
plt.plot(Q2[:, 0]*f, Q2[:, 1], color='r', label='Q2')
plt.plot(Q4[:, 0]*f, Q4[:, 1], color='y', label='Q4')
plt.xlabel('Radiator Josephson Frequency (GHz)')
plt.ylabel('PSD')
plt.yscale('log')
plt.xscale('log')
plt.grid(True, which="both")
plt.legend(loc=2)
plt.show()



