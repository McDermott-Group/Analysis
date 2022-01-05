import matplotlib.pyplot as plt

Q1_J6Bias= [110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259]
Q1_J6ParityRate= [1650.3861137388612, 2075.826701605176, 2441.64598668285, 3630.923667298121, 2097.1636440049197, 1910.705694424216, 1721.535423472802, 2049.744223835382, 2666.5530360953408, 2970.354693964023, 3652.399175319042, 5185.025766031896, 6840.8579567272645, 5737.8544266407425, 10079.137305125034, 10347.860191106469, 10563.978329569823, 8042.711961566178, 8255.218969531894, 6921.373412093129, 7293.207633698482, 6033.9782744272425, 4547.199830677469, 4459.568137761193, 4933.359506818861, 5712.7737501434685, 6793.968518447016, 7158.239768523455, 6569.617123133405, 6656.2072612376805, 5683.104603695633, 4850.913245902843, 4325.927041458489, 3760.858897958199, 3654.2540280914563, 3633.961550846235, 3602.291173856141, 3708.5464363157507, 3027.6225374992773, 2910.4900724014806, 2914.632652510983, 3216.9854302777676, 3753.5457891590777, 3494.5060093954216, 3627.8279347144107, 3868.743463665192, 5046.477367975774, 3457.7093334320725, 2593.238676860829, 3265.471462539448, 3909.857888475861, 3933.2136616051002, 3676.156178514995, 4129.150707522869, 3754.5221087255873, 3607.031766281465, 3865.9035713787507, 1856.2064825920336, 3938.8108090231567, 3221.775875958141, 4462.296115075852, 3762.7927481773277, 3193.8603985494024, 4081.9285522094865, 3902.313322127078, 4011.352231918253, 3106.285817252664, 3328.5723257163972, 3888.708033527462, 3804.5847834467663, 1468.2690928732868, 1381.4898396229676, 1552.1133924556375, 3157.685497027096, 3419.4643276546, 2741.488641593783, 2869.0775305334005, 2731.1400552733867, 2605.120526670422, 3186.858292800519, 3568.7415273883594, 3222.6676016780143, 2340.3020057246367, 2235.4386812326456, 2297.835723906047, 3427.967235395041, 3510.4104130894393, 3657.8878363437657, 1572.8242081075123, 2218.0385654946595, 2635.2197528679844, 2622.5837363828227, 3078.8450251363934, 2447.899092757938, 2281.177279306325, 2551.0923240348425, 2876.793412800907, 2327.5031745240476, 2594.5014921403417, 2634.3948235403213, 2823.5117276146198, 2232.168479806405, 2174.838386927689, 2175.1649675402127, 3288.155066065317, 3401.8814772216892, 3254.7338250875855, 2799.0480041796754, 2934.5743227243247]
Q1_J6ParityUncertainty= [37.63155830195622, 60.452798186411506, 22.904316819606098, 41.581435030476044, 20.33504986726265, 36.22727079114777, 44.573241091587796, 146.15005397230516, 36.38526070340112, 51.938063804615375, 103.76394250676437, 214.16392541121536, 210.7844079439541, 1274.3634027808118, 121.6479107806503, 147.20528913742285, 132.28364253144323, 283.97155311366777, 562.1231237416846, 283.1702171864233, 60.58361475813755, 182.65880401284682, 98.99460113539402, 81.24705741882399, 175.54276469895365, 119.59299924536049, 220.8854059931617, 108.98382977926093, 229.35184004223817, 258.8096588339142, 394.25693932395296, 57.30011733818069, 92.62431926814547, 57.29117292720282, 39.339647324537985, 46.146589444544595, 80.23379772622454, 41.31053424024163, 41.92110540542216, 44.113066899324174, 46.793686700127104, 33.39846922836419, 102.00306219135553, 85.74697128547764, 279.20388702218025, 39.11167153629804, 208.66738115397825, 69.92571494557154, 460.86034432911663, 134.549166586642, 63.7585509253093, 113.78971000754804, 66.74773428725423, 73.8947968764186, 53.26230630055568, 198.08803814101475, 55.081442751049494, 408.20321161059667, 86.88882425330276, 240.53626960004294, 282.1147644437374, 67.1475153078914, 82.29496351408332, 47.80511358107216, 38.99069846216343, 407.50381053831916, 82.48754313714346, 30.062973935060125, 128.18453072561377, 55.83522468436939, 41.606669223773146, 51.7161633676027, 35.19020354098975, 58.95155760225498, 71.184180623516, 38.54635518627277, 129.50163652222278, 50.71953604266633, 64.90018515998702, 153.34237132518174, 80.51859389895448, 66.86906749602544, 87.35113739373953, 44.19527816735505, 269.1296817335914, 179.68304746910815, 90.73096321004286, 78.1112462170985, 316.04159880622143, 53.140226973300045, 72.18348703195528, 69.37050017426586, 275.9464331310495, 76.7571815800974, 45.76030663521067, 39.05286563022819, 83.05240435303034, 124.66264543995001, 89.3985397730081, 59.44364488483309, 208.6062625654445, 46.06302678718216, 61.844522154588596, 42.17768395623416, 101.68735412490106, 156.01667224907544, 126.58525845104398, 145.94225918729217, 139.3454133924734]
Q1_J6Fidelity= [0.47151700439186023, 0.48440036142453924, 0.4877052629085045, 0.500632569082717, 0.6517056456003621, 0.670393510739365, 0.681206412315865, 0.6408553214521778, 0.6060786603802861, 0.5838263387763838, 0.5731662092894502, 0.565796200885404, 0.5605503589156675, 0.6387385499481768, 0.5484131939380528, 0.5531390721471404, 0.5499766752611291, 0.5621253410605704, 0.5528513455204307, 0.5579370277996911, 0.5538978417814033, 0.5612243655825935, 0.5733866746770189, 0.5810505494607195, 0.5743224636522936, 0.564007657865819, 0.5500573391024872, 0.5421485423409463, 0.5471750622998784, 0.5439577466240411, 0.5226562976314273, 0.557696838571175, 0.5612563118972075, 0.5443004201342939, 0.5383008326639627, 0.5132442985178689, 0.5166179567283983, 0.5160295578936336, 0.5151251855614811, 0.5308989218499748, 0.4899039772313516, 0.483438028850721, 0.5207810431954316, 0.5153481363039245, 0.5469846428167146, 0.5137462140003989, 0.5079504460707284, 0.5514894931736567, 0.635144861016634, 0.5641290649227958, 0.5143994340907105, 0.5132149169989098, 0.5097444962502554, 0.519687912193904, 0.5574820329537198, 0.5992876287783215, 0.5716923370706243, 0.6344902572803038, 0.5933668317739668, 0.592022568832531, 0.5436267064933681, 0.5538235914689869, 0.6128103460805208, 0.5522044519092538, 0.5741561534151324, 0.5942874928213121, 0.6341284068481283, 0.5400198028728558, 0.4822148453231408, 0.38723403485647245, 0.5576339861313936, 0.5037947846187731, 0.5597114301144807, 0.3952745857631298, 0.3840434908885506, 0.3847809101684608, 0.3689676332135475, 0.38245000119743866, 0.3906979695790298, 0.38072517250780524, 0.402724428237025, 0.44790148528048845, 0.5731760208149428, 0.5585776676564598, 0.5526823861591571, 0.4447711635212601, 0.4821954295400498, 0.47818697640308294, 0.6350063916395576, 0.5577087881235895, 0.45012853193394065, 0.3701400899685299, 0.37857151808459694, 0.3895505185374856, 0.4090321095144311, 0.4480987763085994, 0.42803571152048353, 0.45350251989573703, 0.452840087490081, 0.4556262702879654, 0.43462572690903223, 0.4858744351493466, 0.4656019728334699, 0.5013605504439279, 0.47229681128687, 0.46172346888487475, 0.4722897894356389, 0.46642610930692635, 0.4553030184595185]
Q2_J6Bias= [190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259]
Q2_J6ParityRate= [72.37353720178726, 77.34452123022453, 82.24303877398074, 90.18211267410676, 89.8466342847496, 88.01528088371833, 86.90656289617746, 81.59032740748003, 84.43472494758934, 94.5186472245878, 101.65566536990664, 109.44512520829184, 119.1430017719559, 129.5864386856868, 135.87320251875425, 131.9272744010228, 131.95947607350655, 140.6930484430729, 143.40524714951795, 155.05266003275335, 152.83349814932745, 147.02485595536217, 127.63809608765291, 114.33781778917407, 95.95802503595374, 95.12915189814427, 101.97289170363084, 112.75036911155999, 122.11587275321173, 131.23976510534868, 131.1967177581841, 130.83227216672964, 130.13389340447577, 128.21947431648795, 130.50282148419623, 131.46489951767825, 118.94746158311948, 108.55056793276572, 99.33693789623567, 89.75939052038929, 90.31576231255272, 94.3175971826566, 93.63266555092149, 105.33688237276904, 111.87460078741742, 115.15355420846326, 115.81795023734885, 115.49067780123173, 113.31151317953608, 109.98942114554734, 97.41624735509347, 89.51083170548519, 83.77507436481365, 85.50702742890708, 92.08798304485511, 94.74394574367318, 97.27112764847476, 95.43466957352095, 87.16424922254893, 80.16134632447233, 82.29639338449235, 78.68153357912236, 74.52207769732085, 73.72787253814961, 67.99999755267785, 67.94156388580897, 72.82827037266631, 78.22909985985214, 87.50058345873133, 94.00622949757856]
Q2_J6ParityUncertainty= [0.6987945685279944, 1.1023295586215882, 1.4910693065040592, 0.8957756778514722, 0.7375341504768008, 1.4512403439693449, 0.7308119852778634, 0.44268972917701627, 2.504238358249534, 2.1022766058396756, 3.577482095332835, 2.426396836824336, 0.8596554969739912, 1.6855525197916983, 1.99516968890546, 1.1099715044928706, 1.7456142100102576, 1.3485232313347486, 1.0669627479049486, 0.765149268446978, 1.6354629930731746, 1.8065204937056711, 1.9675485792973697, 1.2846040394527944, 1.8242271637614311, 1.5426906923852532, 2.0971321598775066, 1.544597593554614, 2.50374407113618, 1.0450773630828585, 3.7753690259122465, 1.9061793659430761, 1.8201010576680177, 0.6465350125421534, 2.071650136578643, 1.8472444652399422, 0.9645737516020297, 1.0854327513748474, 1.5587483396590982, 1.2598078965370825, 1.0498314590146203, 0.7681629912367117, 0.7939723601629621, 1.0037172673614148, 1.1951535823780222, 2.5566503041585817, 1.9133818978052108, 1.357038149122378, 1.5834487467151181, 3.195416198879965, 1.5925187628294948, 1.8199283788810363, 1.6381356202591628, 1.227551949845834, 0.646689970385639, 1.9626922435769942, 1.3221234849527335, 1.5889473701503283, 1.1077285589524306, 1.6814270479548405, 1.3139095696461969, 1.9435717756067963, 2.2684241191807124, 1.5809690642370506, 1.2414254267575604, 0.6711218011534098, 0.660099217510758, 2.8762931286658038, 1.2573422196769755, 0.9444232087509912]
Q2_J6Fidelity= [0.6901339450068663, 0.688598775427782, 0.6902607750121869, 0.6360811209629883, 0.6770126870172641, 0.656088227322227, 0.6835163707291084, 0.6842876693008636, 0.6893449444126114, 0.6588085880240617, 0.6887822618449995, 0.6856304514486495, 0.6901486957283448, 0.6675655595131981, 0.6567463906415903, 0.6907905512057124, 0.6941368172558835, 0.6821672400130461, 0.6953119037899201, 0.6942577057980224, 0.649286255457092, 0.6900053873847615, 0.6836467543817467, 0.6480767510326113, 0.672105075686922, 0.6904960481496463, 0.6958418012039058, 0.677653580397265, 0.6946176119115272, 0.682500287501383, 0.687532979165988, 0.6919986802121129, 0.689526623397797, 0.6933803181843412, 0.6938769556223834, 0.6939321556967682, 0.6881381662111924, 0.6909515499473761, 0.6982094930578377, 0.695196508718239, 0.6960589410838762, 0.6630206113182965, 0.6721650574500049, 0.6858772630534694, 0.6931723009822739, 0.6932587667258207, 0.6925070682016586, 0.6909280955648025, 0.6750516779724605, 0.6895287083273792, 0.6798508272592533, 0.6931353979235906, 0.6507310629381309, 0.689046302624718, 0.6759380234486987, 0.6608806833332311, 0.6952538289442728, 0.6661392646429489, 0.6897525678151906, 0.6700979342223624, 0.6892099754663525, 0.6776321942823106, 0.6803071133895628, 0.6823456814986095, 0.6886562291500766, 0.6374395754641531, 0.6970086443081798, 0.6744246292193774, 0.691819522869958, 0.6899908002831816]


### J6 Bias offset
J6_Offset = 0
Q1_J6Bias[:] = [bias-J6_Offset for bias in Q1_J6Bias]
Q2_J6Bias[:] = [bias-J6_Offset for bias in Q2_J6Bias]

plt.errorbar(Q1_J6Bias, Q1_J6ParityRate,yerr=Q1_J6ParityUncertainty, label='Q1', ecolor='k', capthick=4,color='r')
plt.errorbar(Q2_J6Bias, Q2_J6ParityRate,yerr=Q2_J6ParityUncertainty, label='Q2', ecolor='k', capthick=4,color='g')

plt.xlabel('J6 Bias (mV)')
plt.ylabel('QPT (Hz)')
plt.yscale('log')
plt.grid()
plt.legend()
plt.draw()
plt.show()