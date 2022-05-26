import matplotlib.pyplot as plt

Q1_J6Bias= [0, 10, 20, 30, 40, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410, 415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570, 575, 580]
Q1_J6ParityRate= [1027.6999811130868, 985.7414490947174, 1094.717952229256, 1070.1088029421448, 1192.2222700421287, 1201.451390009225, 1028.8806996512008, 985.1913052499343, 1171.8131635569275, 1050.9410452095024, 969.750686464911, 1162.763100024833, 1033.4948631731081, 1030.7265933956137, 959.2472461310151, 1246.6703420656922, 3950.4339602641194, 6558.676587037766, 7346.150787339924, 5874.855746406484, 12251.642973164588, 4606.533713930108, 3060.256566531238, 2750.168235782164, 2959.2785571912755, 2480.5762211971614, 3341.149467550169, 3235.816396059684, 4100.334953888143, 2351.2576603955126, 2067.01713748069, 2817.9742360399614, 1493.5871871818035, 3672.0734471274905, 3325.418148859013, 1349.8409086740116, 3373.079391924842, 2705.3251192429357, 2342.4351741802056, 2974.3574446985185, 2418.2269975597333, 2153.633764068821, 2277.7495882460344, 1760.4657698060885, 1886.1602565628637, 1473.6343773817587, 1616.9493786118023, 1551.168898267046, 1705.774596542017, 1507.7761615912532, 1575.038684659755, 1636.5656171611076, 1902.5686635799225, 1755.420394244821, 1690.5141595214081, 1651.179839686399, 1806.5953770880399, 1620.724928211042, 1628.1990373142633, 1807.5749280810032, 1950.7346777007376, 1741.9081734039912, 1483.0074736007168, 1397.183530843422, 1573.2853603013054, 1518.967077083475, 1379.8528212791346, 1261.9204505464518, 1145.6496638731421, 1334.9511170553592, 1259.7118463043796, 1447.8061079469333, 1448.7582689638966, 1473.9543653506473, 1466.4308766829781, 1598.2849747160021, 1269.3635860701072, 1257.63570556575, 1422.0874297514756, 1773.0714901145225, 1978.8931971585594, 496.3411240642167, 661.8597495140036, 742.5528770988632, 6566.086746442973, 804.968026805961, 6518.3276740523725, 12722.708920003617, 14942.397799980052, 12277.166685910383, 15341.195444461466, 18287.746197605138, 14716.706756926504, 14162.131226766027, 14692.820368589248, 15364.42053085245, 15394.158495141946, 13936.147131660524, 13717.812413203988, 14282.923219479766, 14564.384589905052, 15184.629439095863, 7088.80690711874, 15919.44319687589, 11467.028957196384, 17525.42816548337, 15143.577485258209, 16766.675683565434, 15753.08048506938, 17139.18165024498]
Q1_J6ParityUncertainty= [23.233504471320902, 26.87162643467588, 36.54083359779698, 18.258263718319576, 24.262708024985024, 28.455658579124194, 30.682438033936247, 24.752006070034664, 31.92028658133762, 32.080545244911185, 19.642093139135987, 34.32714967251894, 28.66487447591004, 20.574470209382717, 23.51766468777442, 24.9424506633306, 711.8349898089006, 133.36520579260923, 69.1800998691987, 87.42878903338026, 912.806258608145, 1219.0131292982542, 103.56976929302164, 31.45744642581063, 105.06063234819067, 66.46948913900965, 42.37381475313023, 48.91740115222207, 1093.9355455078976, 15.519502305853628, 41.73903316995414, 74.37723935763833, 61.82720569717619, 120.73460566061607, 71.86780523613038, 434.2730301570759, 117.74297112728281, 46.7590563156324, 52.421808866055486, 71.93115028742643, 41.33506225904131, 32.52832358651263, 87.16115188437628, 24.791626587427455, 49.731324147112524, 227.10870789834874, 29.935456810816373, 24.559981852621878, 34.514828307194435, 20.10299023016404, 37.767454027710166, 20.9131509866929, 33.51281905406265, 43.84558335666095, 42.76761351601355, 25.840957601674116, 34.81912295605596, 22.479692366851467, 28.048736065747732, 80.20265450143577, 33.58850002359413, 82.55999426284383, 26.277131839950968, 34.589796847099684, 41.15562479775006, 30.16242248071916, 74.97568626171869, 17.98689030596812, 228.48235447780644, 29.857078350120542, 30.71937554778986, 39.46527947339495, 43.41073030366604, 103.8293305116793, 34.095829570287066, 25.386357586785564, 191.1932632745876, 28.777869080588236, 37.38403390664095, 67.72419984526582, 62.24114513630783, 76.94908083136905, 18.31362519248018, 16.29306082417282, 2740.6685519166717, 20.04145346853793, 2734.5608396954963, 2190.7422987651526, 308.29129551582315, 2090.786493526533, 1829.16769450456, 261.0219318976312, 391.79160905054715, 196.75196194075306, 304.817210329897, 229.59569447147206, 217.66321144116017, 319.60851500018856, 182.25889432349766, 322.7544788238045, 243.92011216505813, 262.39414051342555, 2633.176390197768, 300.33932432243495, 1995.3795538602774, 181.5245117980447, 223.50756416330208, 356.81122316185605, 1080.9083225308839, 391.2449333657705]
Q1_J6Fidelity= [0.548836631872006, 0.5947870336208269, 0.4153049938920921, 0.5218549997777842, 0.48581070286855843, 0.4818871275974933, 0.5567906718076406, 0.5558386324469835, 0.5203018818955221, 0.5211577530566173, 0.5598659993919746, 0.49997318954881675, 0.5436404550278275, 0.5595180999584379, 0.6155575155290045, 0.5624852968853582, 0.601086420760224, 0.5852892778879317, 0.49962691936286646, 0.6326216338783478, 0.5357659911712976, 0.44571143696321414, 0.4614958369450018, 0.6005071759038763, 0.41105477214085695, 0.47126809046796636, 0.5894179845509465, 0.5765731775859338, 0.4135322571615922, 0.5805114484525523, 0.44325149483965903, 0.4412398822182973, 0.4553881132854561, 0.425065607454563, 0.4895926972334871, 0.5093772825468207, 0.5274458631106987, 0.5444640420617799, 0.4950107132326249, 0.5011089409504733, 0.5157084213301683, 0.5022692230707466, 0.4666128958095407, 0.5020998147042438, 0.4802197628965706, 0.4615163771335985, 0.5290649610007716, 0.5527575606580076, 0.47719727176725923, 0.49204859553627667, 0.4744713510694325, 0.46184963716026023, 0.42399429647054426, 0.41533282106645225, 0.4028697322527129, 0.4295557893732248, 0.4508302681711333, 0.48561125991416326, 0.5093033995577725, 0.4959889175623813, 0.46868327774642415, 0.44052065417212777, 0.49883995691466004, 0.5341289139555434, 0.48046125023912295, 0.46452918951570954, 0.5181182573844005, 0.5342384139693453, 0.373336516619549, 0.42208799135960173, 0.44892747281152257, 0.41540772843695145, 0.39596244228351185, 0.2685928969145388, 0.4051368706463359, 0.3818893018267315, 0.33725917940107625, 0.403876692314255, 0.3763958365304816, 0.2849485519828875, 0.35027498203951946, 0.36325304946436426, 0.5151143424422251, 0.42964001380687955, 0.46481684256101663, 0.4446281064144278, 0.46863578012736595, 0.5099318604880906, 0.5067128463335357, 0.5179911394754118, 0.5055190214570674, 0.5448544606584691, 0.5231256189748246, 0.5303515499343918, 0.5828192250781634, 0.5980715720419519, 0.6307368571918611, 0.5454671253036981, 0.6534302450593964, 0.5540390943556375, 0.5529814628570238, 0.5454388475674462, 0.6997374898560672, 0.5568848726472938, 0.4233388543887601, 0.5870726492531415, 0.5781238568712481, 0.49950248221860133, 0.43859820279049444, 0.4934778134990394]
Q2_J6Bias= [0, 10, 20, 30, 40, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410, 415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570, 575, 580]
Q2_J6ParityRate= [12.41254627890859, 12.519418873460548, 12.95430006436121, 11.565078753588141, 12.744246337774495, 12.77042396106378, 12.475466397666608, 12.08662883192739, 13.648269694917712, 16.72963055783721, 10.95939140291254, 12.503219314009417, 13.451516512032132, 13.606313890056228, 13.478676794096462, 13.321928113473799, 12.728649802097351, 13.951429714022755, 45.44740439184656, 69.06834360304424, 60.72467402803805, 48.972626294711915, 46.29022726705678, 33.88778647699016, 31.80739502615266, 30.38734743057439, 39.08265578693789, 48.691401169470126, 74.27115828965631, 68.03747298353515, 53.980171066351645, 76.3861634223426, 87.03380065066078, 102.71753479517093, 139.12776148335655, 135.7204340867366, 102.8325247399761, 132.24303664883487, 111.38242931534856, 102.19092891113618, 114.16723668311764, 99.17083433022658, 93.17007959453593, 79.28718549891886, 79.12437563835914, 93.29547457460572, 73.30156567588101, 75.86862411436564, 89.51441398832539, 82.07188261971837, 90.67073811984855, 79.70812013344428, 114.4479561816222, 91.68290417377995, 104.47730420139317, 138.3231511845986, 246.45041704699878, 286.8036376697918, 320.4233801749724, 413.76378490635545, 474.5125881344328, 418.874630294809, 339.7308499251966, 343.6341347442207, 394.7217604115802, 417.4770905756551, 427.5882828055328, 438.6807363797246, 534.5981284825494, 653.4772111295017, 585.4869539598839, 657.0544118431079, 638.9636526353089, 658.7458526195489, 646.5803634148963, 596.3142804868318, 494.6993898309914, 355.1689290500558, 336.72952342992625, 343.3295539455752, 630.9564631025339, 276.4363459326377, 238.94138814382688, 281.50644084555904, 232.4328305368075, 190.9126417391317, 181.45303367070068, 177.57725449654043, 144.55268914677555, 143.97731651393218, 144.1502091255896, 139.79176083416274, 154.49163433471378, 166.5579220784, 162.68970136868654, 157.93410166143065, 141.15171171795978, 124.8663389665879, 125.38965738078849, 127.53096557019691, 123.92536638225471, 128.60566233478912, 120.88316689598311, 125.10780040900548, 128.2549576091543, 135.31508433518366, 145.26541720883154, 152.09842260944052, 154.22147528068805, 155.21296609567813]
Q2_J6ParityUncertainty= [0.20044632397252116, 0.42701869178093155, 1.1630354716485733, 0.24664432263583969, 0.04794652810305333, 0.1544971280388605, 1.1248064915877456, 0.22111192686390524, 0.18785437992760112, 2.546990730246226, 0.3097075477612128, 0.5189775319552048, 0.22623340315946106, 0.39347936759706137, 0.4681818708492509, 0.418667165067542, 0.7067059771472861, 0.29235192736289584, 0.47334814858893104, 0.3807509007853369, 0.5287357001563997, 0.43812274906627735, 1.8402173851963912, 0.02769021563249474, 1.1200017373648627, 0.2592728370967432, 0.9251625474152494, 1.7495382868757225, 1.4430188373704098, 2.6795603167772235, 0.8870195575426187, 2.8203212761807315, 2.1180018601116153, 1.999623989525027, 1.145788334608973, 2.5957716067689276, 0.5544678883887357, 1.4915059312549275, 1.3430333065021713, 1.1901801328835049, 2.007187783111867, 0.8274020495872279, 0.462058230294609, 0.4388116040512635, 0.7552352689823145, 2.1300720082196953, 0.7633742125471876, 1.0112673171326527, 0.10508832523825262, 0.45435590509403667, 1.1485579238560604, 0.3198802970937358, 2.784724436566364, 1.4148960644553306, 0.02881558485705682, 1.5602469505343421, 1.7287563361592788, 9.489403370086634, 3.512539339160668, 7.765912177115678, 9.675335906393665, 5.741815991075276, 2.406320544297813, 1.8148270275027159, 5.825396643596406, 3.5511719905980286, 4.916162057547364, 8.676239931943615, 9.523823563340066, 9.412239940007055, 10.301959546638876, 10.753701769786353, 7.1522155918830626, 9.714675360006655, 7.391314342851995, 8.85958476729866, 8.172980255541885, 1.4585640180571744, 5.403053505555505, 2.4209151527608412, 168.3090758630442, 1.1580263447063288, 1.4124709613150515, 7.439820929003858, 2.7651675915273444, 2.264757827135611, 4.526799412900672, 2.991998277538549, 2.7621711837553207, 2.181672220513704, 3.5974554754610253, 3.2328517279925677, 0.7043953436900878, 0.026486745882152718, 0.6357567448067414, 2.171082937554182, 3.692328789128389, 0.1755750061786614, 1.5264146449107088, 0.5605609543325158, 2.682812587319262, 4.652653980624116, 0.8528338070464301, 1.303394212102539, 0.942093481927202, 3.2455817674984075, 1.0487109969040418, 4.307378087828901, 1.992622879616249, 3.3049761714541006]
Q2_J6Fidelity= [0.6947895432046824, 0.6954416083088915, 0.7003706295571948, 0.6930192820635661, 0.6978770514763946, 0.4003028094079798, 0.6922248772822026, 0.5623861947844012, 0.6471417352290698, 0.8549387468594075, 0.6823739685056867, 0.6946781777014072, 0.6989400337388789, 0.70679135578448, 0.5333553884187727, 0.6999653024953617, 0.6939167384506555, 0.6754737543480153, 0.5309783768067114, 0.647901589837003, 0.6956840626659645, 0.6669013987936823, 0.6930967408233653, 0.6930687618157106, 0.6752653569449537, 0.6960388584102518, 0.6930651516242563, 0.6967908323028746, 0.6093302932307141, 0.6864606815862091, 0.6948977512680493, 0.6919961051117333, 0.6984026733091361, 0.6947310993358029, 0.6959426781581801, 0.6965427727083542, 0.6903936535583133, 0.695768330545268, 0.6926419911016819, 0.7022900714745215, 0.6938978620807417, 0.6900708187788889, 0.68624745751119, 0.6869231064137181, 0.641420014736984, 0.6950282964659238, 0.6879319731650564, 0.6216587867908474, 0.6824957755792862, 0.6929718835577288, 0.6260882231563962, 0.6981028137132796, 0.6920751906096736, 0.6863812287116475, 0.6940487931572928, 0.6827143997661913, 0.6860880356993744, 0.4654322798008318, 0.6760874878806427, 0.6561002543576844, 0.6787730082791392, 0.6750335318806225, 0.6814813758478412, 0.6751206669568484, 0.6718869031787329, 0.6873943172392625, 0.6901991373822678, 0.6803734901212002, 0.672985494460547, 0.6731575145916185, 0.6776868194384333, 0.6827583923608375, 0.6729787667141631, 0.6870856141865149, 0.6799028122440077, 0.6872682992988003, 0.6859995915231388, 0.5965138684588935, 0.6954638972372776, 0.6259079913295262, 0.4562363635593371, 0.6953334608095243, 0.6972350483585501, 0.693111726506376, 0.7050203680025452, 0.6999445117620373, 0.7036034897066199, 0.7032942321869091, 0.700551947850498, 0.7040248182884785, 0.7038322409990347, 0.700405796403744, 0.7017768385785987, 0.6841127528328398, 0.6798134854460056, 0.7085917947534337, 0.7035560193860368, 0.7069680083716633, 0.6996178908543633, 0.7033814430967117, 0.7027019354707611, 0.6998402971543463, 0.7032670022321754, 0.6969587052890482, 0.7056640902765359, 0.6744120428406004, 0.692841363806268, 0.6965978278353272, 0.7015145247731847, 0.6896619023872106]
Q4_J6Bias= [0, 10, 20, 30, 40, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410, 415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500, 505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570, 575, 580]
Q4_J6ParityRate= [191.44356342562892, 188.90025823771003, 189.5452214794004, 190.0261344054667, 193.24667717278902, 187.87128771841844, 204.8028361725847, 188.21462951319955, 188.9646566588445, 193.12598189025084, 183.59587412419484, 190.40137130660838, 194.50043211824888, 192.73170253051774, 195.3437703459071, 169.67009592484044, 198.54590759686607, 235.8228428893632, 295.22480212817305, 330.64663852735634, 447.3042479799527, 453.0287771991792, 603.2743851596797, 552.7244787968829, 347.5922283454224, 337.58580456163696, 498.9374460457916, 663.3112531461004, 902.0137760356786, 815.4956852025746, 753.2062781854531, 1349.58507476166, 2312.677012397331, 2951.910061634431, 2658.131395774906, 3984.8317774609177, 3440.896802368123, 2962.0358959355785, 4703.1165533473, 3496.127438971389, 3611.80206172426, 2973.8285837527965, 2770.0160346855046, 2383.8270051831496, 2485.731697000874, 2609.349459015681, 2588.2821075515044, 2274.5676314213633, 2673.646733027699, 2761.6118266947456, 2009.3956860972091, 1791.1016475511199, 1662.9131051738761, 1211.0712495398602, 1271.4929131910735, 1286.0129422575658, 1365.6047642286371, 1849.461315881978, 1898.0963397496275, 1398.4899035543735, 1230.5404403028379, 1524.918556460855, 1282.6062244636598, 1148.0052833474908, 849.0458948960475, 720.9611391126106, 761.2911339006979, 671.2223800506986, 768.0590347116101, 821.1358257109865, 693.9799160285622, 630.7069136644304, 677.260962254159, 658.3455045630053, 660.5310466942336, 600.1736297793932, 610.5123229706235, 669.6864223557526, 721.5708259524592, 687.5035861349783, 890.3809932054616, 860.0746186474544, 867.6538043288466, 925.1120998872702, 926.2480669026572, 969.734565247892, 469.9129608381548, 965.5070774364913, 891.1092372752897, 894.6138996547149, 942.3269709653765, 935.8814549411505, 957.9060036070277, 971.6131091500434, 1032.5557496447107, 1019.4791095911472, 1001.9080875078387, 1062.4784000034529, 1109.622020140117, 1038.1728366680322, 1058.8671893304129, 1092.3173523613202, 1071.301555344408, 1072.5727995100674, 1157.5889698748263, 1235.2574810071092, 1287.4987267378715, 1300.012869858657, 1358.0179322817482, 1413.0592425188993]
Q4_J6ParityUncertainty= [2.3299587775045154, 6.502786635713775, 1.8840454275150051, 3.584929478876205, 3.162854376668676, 3.0544352329620494, 4.031691877538179, 3.377199932084693, 5.850350496000459, 1.7482945184152072, 10.833543545131718, 0.17910735800556665, 1.6745934554209612, 6.329910363726199, 2.0617627785744617, 3.0224644026528082, 5.03686512706615, 5.437103545688771, 7.950834701413783, 4.705381454421982, 11.495363071771571, 8.11975796326452, 18.824307297906714, 12.28648679891808, 3.849645280561262, 1.1389235617113656, 10.059019058055872, 18.101363709385375, 17.689803291015735, 13.06568216113432, 10.972916966130501, 16.69204222752532, 44.680560181081546, 38.88406117433233, 42.101718798286896, 63.52572522853174, 39.10207298670139, 54.47193665731168, 61.0863200792955, 74.89064548239412, 33.07237858492213, 45.132420768615404, 42.86552848360237, 32.66538130790975, 38.43620214156327, 28.971490979357686, 47.65349738772886, 44.73990741452617, 36.03658425821525, 42.55724102724202, 26.849891754710885, 30.031095281005673, 27.89033490831196, 16.019810733893166, 25.797159556408964, 27.308848476280943, 36.55849820233722, 20.443400935033125, 28.323946566929557, 14.46491410397254, 20.47508993006125, 27.33259971356928, 23.207087046071095, 16.238184680754294, 15.464309706647054, 12.491389419413926, 13.796504971511936, 15.63692730181666, 16.323965852981893, 21.97975381760305, 13.621854220597424, 14.383545343904299, 11.966117905321692, 13.243057823150053, 16.942293515681047, 14.645947018634843, 11.386920230473343, 10.927063518922331, 14.13871893902998, 11.745859311294653, 19.588335098221577, 17.221667531371132, 23.004552650185023, 14.98775673892792, 18.250471082871876, 19.781728823537737, 156.76577137969247, 20.80863761696834, 16.493753626711857, 15.106490992921207, 24.408641505168312, 29.050178578837127, 19.593729770015667, 13.631627619124163, 16.016305519393786, 15.210088891564993, 21.79789274834903, 194.12964500719292, 25.676941294563676, 16.052553235426526, 12.127673825391554, 19.375859462467798, 14.819947891864617, 22.090356148627496, 17.79976729658752, 19.62312567259699, 16.238232144937886, 14.244839944330186, 10.38320494908792, 19.767422719476425]
Q4_J6Fidelity= [0.5462251946450336, 0.5974410707677162, 0.5852634411690458, 0.4670968976346529, 0.5584400615628741, 0.5721722143864361, 0.5795547593659326, 0.5719260856515813, 0.5762790614568322, 0.5790138051749064, 0.5751803714542147, 0.5253335424462451, 0.5811712916858753, 0.44386701406422235, 0.5789780871910244, 0.5740158219143906, 0.5807400841832316, 0.5807423092815107, 0.48889501707361604, 0.47696616785479906, 0.5851074383185039, 0.5834683559237274, 0.5878930941640814, 0.585837588852408, 0.5716531582041496, 0.5840226571969078, 0.5657184502017213, 0.5786247500039206, 0.5685494261702712, 0.576674113564781, 0.582059199700822, 0.584287563987514, 0.5991314581125469, 0.5772926135393328, 0.5868360944385235, 0.5718280886592916, 0.5796893331885168, 0.5720741076348009, 0.5671846482341931, 0.5672669706666199, 0.5660625181375517, 0.5710969273275703, 0.5688576366457242, 0.554967978270455, 0.5689774640892106, 0.5795057555112202, 0.5759483425907097, 0.5803464681304328, 0.5767005320924581, 0.5879895096129013, 0.5875611595687811, 0.5832633230280982, 0.578746761832717, 0.5893949200414823, 0.590725605850231, 0.501495685263935, 0.5803401259763461, 0.5870693941673392, 0.5761350916607053, 0.5849034165639594, 0.5936951070201351, 0.5821483176795803, 0.5792555676228069, 0.5835202784514595, 0.589870151209977, 0.5973756767591086, 0.5544584758736689, 0.5825521136187746, 0.4244577068682799, 0.5921969875182856, 0.591927090806896, 0.5891516785285157, 0.5906630279470502, 0.5868880047632683, 0.5908289303237201, 0.5969775847287306, 0.5986334305371359, 0.5942604450590169, 0.586427773763269, 0.5903941186711145, 0.5961234454750037, 0.587882463494441, 0.5124250382888644, 0.594173180672017, 0.5959378081040868, 0.5957027577078751, 0.7763977169104906, 0.5915154384586501, 0.5921062906852219, 0.5845097151110554, 0.597803641960906, 0.595542819532545, 0.5939243286940961, 0.5912738585600411, 0.589917818105507, 0.5940554696801689, 0.6007710929700573, 0.37772483636191445, 0.5972643342022594, 0.6032161595383321, 0.592170078959817, 0.5618087474048183, 0.5844965439039858, 0.5931100683588675, 0.5987099304363432, 0.591551075300659, 0.5956548268498663, 0.59434568092095, 0.6004342583964613, 0.5931607987963315]


### J6 Bias offset
J6_Offset = 0
Q1_J6Bias[:] = [bias-J6_Offset for bias in Q1_J6Bias]
Q2_J6Bias[:] = [bias-J6_Offset for bias in Q2_J6Bias]
Q4_J6Bias[:] = [bias-J6_Offset for bias in Q4_J6Bias]

plt.errorbar(Q1_J6Bias, Q1_J6ParityRate,yerr=Q1_J6ParityUncertainty, label='Q1', ecolor='k', capthick=4,color='r')
plt.errorbar(Q2_J6Bias, Q2_J6ParityRate,yerr=Q2_J6ParityUncertainty, label='Q2', ecolor='k', capthick=4,color='g')
plt.errorbar(Q4_J6Bias, Q4_J6ParityRate,yerr=Q4_J6ParityUncertainty, label='Q4', ecolor='k', capthick=4,color='b')
plt.xlabel('J6 Bias (mV)')
plt.ylabel('QPT (Hz)')
plt.yscale('log')
plt.grid()
plt.legend()
plt.draw()
plt.show()