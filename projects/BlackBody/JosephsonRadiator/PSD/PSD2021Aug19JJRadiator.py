"""
PSD for Q1, Q2, Q4 at different J2 Bias
Z:\mcdermott-group\data\BlackBody\Circmon\LIU\CW20180514A_Ox2\JJRadiatorQPT_2021Aug19
"""
import numpy as np
import matplotlib.pyplot as plt

Q1 = np.array([
    [40, 2235.5179104863346], [42, 2981.0616813309734], [44, 5950.149229218026],
    [46, 6189.936076280245], [48, 9589.316698210283], [50, 6043.995172261502],
    [55, 6414.099382926821], [57, 5473.670509160685], [59, 6702.015229610355],
    [61, 6225.628969187118], [63, 5557.32577070727], [65, 6685.567248816772],
    [67, 5325.968673918348], [69, 5913.947627420805], [71, 5934.69302284728],
    [73, 5822.622146360961], [75, 5293.434508792141], [77, 4968.368740278515],
    [79, 5137.8393876761675], [81, 4449.584295033822], [83, 5000.0],
    [85, 4314.081798237319]])

Q2 = np.array([
    [40, 53.618002587631956], [42, 51.34845389221705], [44, 70.11655240131854],
    [46, 187.41218622320065], [48, 304.7456358290479], [50, 286.64682064343367],
    [55, 220.7810843830252], [57, 295.70896403010227], [59, 453.753251555533],
    [61, 445.77794284731294], [63, 466.4695896159714], [65, 572.1022829482743],
    [67, 435.33431250626023], [69, 508.51697311146785], [71, 526.7609333389622],
    [73, 430.61188978895336], [75, 393.4080389484947], [77, 341.06249080707795],
    [79, 419.1239775341623], [81, 336.30317863717437], [83, 320.9756931918318],
    [85, 318.75037467125316]])

Q2_more = np.array([
    [85, 319.66219092265396], [87, 354.38695756161417], [89, 500.8550106864243],
    [91, 548.8167252522882], [93, 551.4857701064179], [95, 714.622933629825],
    [97, 787.1947517927705], [99, 934.5880705643397], [101, 1050.7404397665698],
    [103, 1102.9306737771517], [105, 955.9147207975321], [107, 901.6836970623358],
    [109, 742.9656940712134], [111, 608.3924703037844], [113, 416.7774404494524],
    [115, 364.175788200114], [117, 347.3336099153992], [119, 356.3232782844768],
    [121, 298.77625905787977], [123, 318.27763598172726], [125, 321.07398246377346],
    [127, 316.4791026113248], [129, 328.22588509113723], [131, 337.32903784969346],
    [133, 355.74941634706363], [135, 375.2785830561038], [137, 386.21752658077435],
    [139,420.1681011213989], [141, 435.4828686138838], [143, 453.67283711361677],
    [145, 489.4147828788428], [147, 522.9400605994541], [149, 561.4014726543671],
    [151, 598.7175698822276], [153, 636.7810831784419], [155, 680.24751109491],
    [157, 716.3209097497648], [159, 766.3075035087892], [161, 777.0223847876968],
    [163, 824.935997130911], [165, 871.8001201342591], [167, 901.1069046579375],
    [169, 958.1306100222354], [171, 977.6083620671726], [173, 1034.6074968034125],
    [175, 1056.3473944010823], [177, 1125.9647537591509], [179, 1151.069140547619],
    [181, 1187.822827772019], [183, 1237.1178974655957], [185, 1222.2473197309628],
    [187, 1288.5298041741714], [189, 1300.4758594363464], [191, 1327.7401693180398],
    [193, 1331.1801888076182], [195, 1378.0313870805105], [197, 1361.9141666743951],
    [199, 1440.0144248408415], [201, 1420.4118917638598], [203, 1466.6534963513768],
    [205, 1474.819821402907], [207, 1496.1904116992791], [209, 1469.8593571964734],
    [211, 1504.7480532000936], [213, 1487.4710353334751], [215, 1547.8725311991693],
    [217, 1581.296082576847], [219, 1662.5223526974764], [221, 1556.1696047230405],
    [223, 1575.6453674711904], [225, 1492.96415275399], [227, 1818.3667465684657],
    [229, 1562.0096411556385], [231, 1496.9642352259275], [233, 11804.746858281269],
    [235, 17424.25044684786], [237, 12447.030865277324], [239, 21704.73252636716],
    [241, 23784.787355666474], [243, 29123.69150162886], [245, 38770.01858035095],
    [247, 44467.21459416427], [249, 50210.921388859846]])

Q2_all = np.concatenate((Q2, Q2_more))

Q4 = np.array([
    [40, 617.22564405253], [42, 597.65754027423], [44, 729.493218794977],
    [46, 1359.1337736283635], [48, 1642.9805441333394], [50, 1956.1886984156724],
    [55, 3291.911904007495], [57, 2626.321419158653], [59, 4176.030291861123],
    [61, 5151.927893012501], [63, 6480.543800633735], [65, 6082.222543926706],
    [67, 5939.985537958111], [69, 6039.458870205478], [71, 8466.05597435597],
    [73, 7674.386283931807], [75, 7180.200478210239], [77, 5873.526322158399],
    [79, 5513.411701619655], [81, 4285.692033387549], [83, 4076.5726427237823],
    [85, 3585.3093928941466]])

### J2 Bias offset
J2_Offset = 19.2
Q1[:, 0] = Q1[:, 0]-J2_Offset
Q2_all[:, 0] = Q2_all[:, 0]-J2_Offset
Q4[:, 0] = Q4[:, 0]-J2_Offset

plt.plot(Q1[:, 0], Q1[:, 1], label='Q1')
plt.plot(Q2_all[:, 0], Q2_all[:, 1], label='Q2')
plt.plot(Q4[:, 0], Q4[:, 1], label='Q4')
plt.xlabel('J2 Bias (mDAC)')
plt.ylabel('QPT (Hz)')
plt.yscale('log')
plt.grid()
plt.legend()
plt.draw()
plt.show()

# plt.plot(Q1[:, 0]*4.8, Q1[:, 1], 'b', label='Q1')
# plt.plot(Q2_all[:, 0]*4.8, Q2_all[:, 1], 'r', label='Q2')
# plt.plot(Q4[:, 0]*4.8, Q4[:, 1], 'y', label='Q4')
# plt.xlabel('J2 Bias (GHz)')
# plt.ylabel('QPT (Hz)')
# plt.yscale('log')
# plt.grid()
# plt.legend()
# plt.draw()
# plt.show()