import urllib

import numpy as np
import pandas as pd
import peptacular.sequence
import streamlit as st
from peptacular.fragment import build_fragments
from peptacular.mass import calculate_mass, calculate_mz
from peptacular.score import compute_fragment_matches, hyper_score, binomial_score

import plotly.express as px
import matplotlib as mpl
import ms_deisotope
import plotly.graph_objs as go

from plot_util import generate_annonated_spectra_plotly, coverage_string

VALID_MASS_TYPES = {'monoisotopic', 'average'}

base_color_dict = {'+i': 'magenta', '++i': 'magenta', '+++i': 'magenta',
                   '+a': 'brown', '++a': 'chocolate', '+++a': 'sandybrown',
                   '+b': 'blue', '++b': 'steelblue', '+++b': 'skyblue',
                   '+c': 'green', '++c': 'limegreen', '+++c': 'lawngreen',
                   '+x': 'orange', '++x': 'tomato', '+++x': 'coral',
                   '+y': 'red', '++y': 'crimson', '+++y': 'salmon',
                   '+z': 'purple', '++z': 'darkorchid', '+++z': 'fuchsia',
                   'unassigned': 'grey'}

VALID_MASS_TOLERANCE_TYPES = {'ppm', 'Da'}
VALID_PEAK_ASSIGNMENTS = {'largest', 'closest'}  # Add more valid options as needed
VALID_Y_AXIS_SCALES = {'linear', 'log'}  # Add more valid scales as needed
MAX_CHARGE = 3

H = 1.008  # Atomic mass of Hydrogen
O = 15.999  # Atomic mass of Oxygen
N = 14.007  # Atomic mass of Nitrogen
P = 30.974  # Atomic mass of Phosphorus

NEUTRAL_LOSSES = {
    'H2O': -1 * (2 * H + O),
    'NH3': -1 * (N + 3 * H),
    'H3PO4': -1 * (3 * H + P + 4 * O)
}

IONS = 'abcxyz'
BASE_URL = 'http://localhost:8503/'

DEFAULT_SPECTRA = '283.751526:6.493506;287.601379:11.096813;295.031097:2.801403;305.472137:2.626226;307.404083:3.930447;308.360321:9.893921;310.225128:3.961838;311.703918:3.872004;313.032806:3.668097;316.054626:5.160672;320.380554:5.133078;321.298553:13.925406;326.093964:7.292274;332.248291:13.145130;333.495422:4.659065;335.068054:18.960356;336.147919:7.259778;337.641632:13.158745;338.382996:9.423209;341.457825:10.963329;343.543579:1.641261;350.107208:14.771483;353.221008:13.144381;355.704102:12.644948;358.252045:52.690765;359.490479:14.436722;362.997833:23.139444;364.979797:1.644492;365.830322:4.412334;368.139893:1.788786;371.229553:7.809958;372.367676:3.643448;377.157013:1.648062;378.191528:3.817382;380.201416:25.570644;381.213379:9.733955;383.437012:9.783660;385.548309:6.963592;389.269531:17.547304;390.042816:13.935720;391.471161:3.108196;394.526886:5.130955;395.299805:3.567818;398.154816:27.128242;399.385803:2.952545;404.270020:13.601854;407.019470:8.618811;408.281067:18.146912;409.774048:1.355445;413.205719:11.547805;415.212311:73.431404;415.839935:0.893067;417.280792:11.975701;418.890106:8.434253;421.298279:10.868628;422.462738:7.692440;426.243835:31.430904;427.224670:2.512977;429.384583:9.423120;430.648376:4.450206;432.221619:1.788626;433.866150:30.993015;434.596252:1.646225;440.453857:4.504568;441.353943:2.224833;444.435242:3.818611;450.200714:18.508629;450.991760:4.519960;452.172302:22.264191;453.226776:11.091013;455.147949:6.710865;457.207214:1.500907;458.410553:15.481137;463.061707:2.223360;466.992645:4.993165;467.639313:4.073272;469.348083:9.050010;470.393494:6.463762;473.183380:5.893355;474.337402:16.059982;476.350952:5.299544;477.385620:2.647703;478.320892:7.154815;479.264374:14.371003;480.104492:7.295769;485.130615:52.425201;486.415497:10.451238;487.197937:1.931733;490.431274:5.386128;495.331726:35.877655;497.432495:9.165106;498.240784:205.154526;499.368774:32.136749;502.205872:43.880852;503.751038:16.804586;504.372833:1.217326;508.100189:14.077838;509.434875:8.070131;513.244873:57.381676;513.892517:4.781098;517.095825:10.307590;518.391357:7.687596;519.081177:4.015530;520.527100:9.813676;521.991821:26.945614;523.370972:12.161747;524.278442:3.381213;525.198730:3.853895;526.070801:3.820176;527.147827:25.483242;528.342773:11.368365;529.232910:8.339369;530.080994:7.815114;532.809937:6.275185;534.071167:8.151394;535.136719:8.521399;536.179688:19.452126;537.128174:10.514593;538.413452:18.023069;544.333984:9.910419;545.273132:11.143578;551.329102:8.597209;552.420044:22.845009;553.281860:9.001399;555.260742:15.997867;559.417725:4.437509;562.404053:405.051514;563.498840:48.098305;565.600952:3.387759;566.737915:5.851104;569.473877:12.515773;570.589050:20.032322;571.230835:6.776375;572.278748:19.095318;573.149292:17.406981;576.461426:14.950231;579.109985:10.360884;581.514038:14.560740;584.206909:30.957813;585.478271:3.959915;588.279968:2.945799;589.511108:10.556177;590.388306:26.146528;591.337036:14.897098;592.107544:4.106139;593.386719:3.813607;594.614868:11.580721;596.236755:24.995085;597.317261:22.972145;597.943237:9.280321;599.349854:14.812905;603.377197:5.381608;605.454956:9.477353;606.883667:21.293558;608.134216:24.021933;609.302917:3.112401;611.359680:70.052658;612.350830:18.329905;613.035706:5.502901;613.796448:26.041834;614.528076:9.602185;615.393188:5.696084;617.266235:25.979660;619.470947:5.926153;620.351318:7.735986;621.099426:17.342321;623.369934:4.083014;624.489319:28.252970;630.051208:5.977171;631.418579:6.180241;632.366882:59.349007;633.249756:6.176957;634.552917:47.269478;635.357300:9.258817;636.320435:4.388174;637.270996:18.594582;638.201416:17.334301;642.301575:63.441235;643.334106:6.157311;645.329346:2.262579;647.864929:11.704320;648.759888:17.945156;651.322754:21.489046;652.223389:7.210659;653.374146:3.931837;654.224792:18.874992;655.326111:2.228123;657.757080:19.246632;658.415100:12.592502;660.134033:101.692833;661.193115:10.691292;661.904053:6.576667;662.516602:4.279579;664.242065:18.312006;665.157166:20.826456;666.211670:2.076535;667.019531:4.251490;668.473267:3.958891;669.413818:28.253675;670.335449:2.223964;671.225647:19.438091;672.297607:10.943571;675.311035:367.105042;676.390076:56.040432;677.274902:2.779641;678.241089:3.061245;679.293030:16.921597;681.385132:22.106785;682.576660:71.931320;683.337036:75.866539;684.422241:9.401711;685.715088:11.327285;686.376404:8.983388;687.370178:12.646134;690.408325:3.990864;694.330078:22.521143;695.348999:6.941419;698.137329:3.668891;699.333252:27.351761;700.453369:12.757339;704.297119:13.704563;704.986328:11.427057;707.710693:1.798411;709.259766:8.880344;710.178589:33.539955;711.190735:62.974567;712.227539:16.642796;713.228394:22.969948;714.115906:11.897521;714.968384:2.656370;716.742188:20.264675;717.348999:39.501438;720.382446:3.775619;721.317688:20.167229;722.506958:9.227962;723.299988:6.964044;724.248169:6.849056;725.511230:7.901319;728.155029:16.628355;729.426697:10.212529;731.317505:3.374717;732.255615:3.239853;734.284668:24.261179;735.587402:5.172235;737.013184:6.724487;738.274902:7.320567;740.276611:8.154497;741.148560:13.379368;743.352173:3.285241;745.010193:15.715961;746.427368:449.830200;747.415588:71.976730;749.643005:15.057947;751.233032:27.113340;752.309631:7.464421;754.591858:23.379366;755.771851:18.915462;756.640747:24.991356;757.367065:3.327466;758.285400:12.644419;762.493042:9.303966;764.344238:20.904579;765.182495:36.106560;766.047363:44.043243;767.351746:93.445885;768.647461:45.062241;769.492249:49.004608;770.466187:22.710712;773.349915:20.466806;775.872864:24.173382;776.571167:12.895178;779.926025:71.074783;781.227539:20.878002;782.137817:24.390005;782.912109:15.161629;784.051025:18.374445;785.101196:12.601723;786.147583:18.534374;787.149414:6.842381;788.437866:21.625475;790.202271:54.209568;791.045532:30.317631;792.775513:7.952184;794.207642:32.385384;794.885620:12.549341;796.421143:76.460754;798.049927:75.988861;798.663208:71.781517;799.443359:209.214066;800.319946:38.709007;801.893677:1.503646;803.322021:6.021701;804.709351:38.214836;805.433594:2.278663;806.286987:2.336695;808.165039:11.217024;809.360596:2.344270;810.098267:1.211102;811.546631:2.510352;812.395142:12.395452;813.171631:6.356281;813.881592:113.044014;814.634033:90.342583;815.365234:66.921944;816.088135:4.924015;817.227539:16.973310;820.229858:13.633173;821.318726:13.969715;821.975830:3.959853;822.611328:25.306316;824.213867:11.723189;825.117798:54.699104;830.105957:17.464411;830.850586:18.916851;832.385010:311.003601;833.427979:979.646790;834.436646:215.636597;835.075684:16.184921;836.556274:19.470211;837.765625:19.225105;839.373291:46.029331;840.378540:108.001488;841.547119:24.809801;844.396729:15.088540;845.276855:3.935811;848.456909:94.455627;850.608643:7.264608;852.140991:23.430210;854.625244:34.734818;855.423950:27.207273;856.728638:183.597626;857.678955:48.737762;858.507690:20.543966;859.484009:21.632872;862.401855:44.788620;864.676636:3.573800;865.399658:7.493271;866.201416:14.760412;867.356445:25.953390;869.311401:21.269566;870.990967:3.855475;872.290771:4.370858;875.432617:4.977350;878.199829:1.646701;879.992310:11.020376;880.596191:30.899002;882.348999:20.763721;883.343384:29.157408;884.627563:19.478388;886.764282:34.861149;887.680298:16.320459;888.419067:28.837633;889.872925:25.790997;890.497437:4.758450;891.749268:5.014816;892.388184:16.043514;893.140137:12.282856;894.220581:7.399072;895.659668:5.926276;897.411621:11.247326;898.477173:8.560768;899.306030:36.849598;900.273926:24.534241;901.123657:31.042982;902.442749:29.487898;904.416382:497.802307;905.352661:107.720070;906.463257:31.142368;909.061035:64.796158;910.182495:14.582764;911.039429:12.176151;915.368774:12.317439;917.316162:55.859604;919.448730:4.802111;920.208740:16.427258;921.361694:21.958057;921.979614:9.060344;923.254761:11.032974;924.682739:8.303555;925.586792:7.690586;927.187744:304.723755;927.882568:40.639816;928.501953:76.819817;929.349243:14.867746;930.167725:30.868053;931.146484:13.427282;932.044556:3.526408;933.240601:6.819841;936.336914:42.253407;937.202271:6.540798;938.217896:7.398048;939.460083:16.356503;940.879639:13.613162;941.551147:7.549289;943.104980:5.383054;944.357300:26.002014;945.336060:514.112061;946.364258:158.641571;947.341187:58.288780;948.350342:52.267628;949.202393:9.984763;950.454346:18.264801;953.539062:36.800938;954.298706:3.925758;954.984863:7.846375;955.801147:6.416565;956.738037:29.660776;958.598511:30.174816;960.617920:13.462435;961.954956:26.801933;964.685791:9.572931;966.092285:13.461919;967.562866:18.168354;968.440796:7.562728;969.729614:58.977028;970.730591:172.426971;971.613281:27.080954;972.427002:7.571497;973.345093:25.057213;974.245361:17.889034;974.929199:5.122307;978.763306:174.032135;979.492798:621.792542;980.229980:821.247681;981.645508:63.345230;982.493042:6.777298;983.343140:4.108294;984.822266:53.706779;985.691650:31.978601;987.477905:52.893188;988.312622:17.456316;989.402588:16.742617;991.427368:1042.701172;992.476807:298.718201;994.179443:396.866638;994.853882:30.691385;996.377441:59.738918;997.650269:18.692841;998.586182:2.225106;1025.497681:29.486681;1028.456055:6.753218;1030.479248:2.228055;1032.459839:2.480675;1033.503662:18.229073;1036.303467:12.162532;1038.201660:16.165134;1039.307617:7.980670;1040.379150:5.440311;1040.979492:7.187830;1042.299683:69.587288;1043.451904:27.655910;1044.548462:1.932230;1049.283813:5.820963;1053.765869:7.851082;1054.454468:3.043788;1055.359863:9.157644;1057.414795:47.249596;1059.063477:11.750919;1060.464355:110.334991;1061.240967:19.877192;1061.850098:48.165707;1062.459106:2.022357;1065.437012:12.424438;1066.324951:5.455384;1067.549805:23.094603;1068.597778:5.413030;1071.799561:3.667963;1076.606079:10.756640;1077.902588:113.072197;1078.509521:1543.968018;1079.529907:403.849792;1082.739014:11.438251;1083.652100:60.678062;1084.357910:39.330742;1085.079468:32.408241;1086.149902:20.945984;1086.806396:9.531338;1090.305664:11.417344;1092.486572:3.693444;1096.873047:9.814255;1101.339966:126.459076;1102.520264:36.289654;1103.616211:12.251237;1107.628418:6.094751;1110.905762:5.358802;1115.560059:6.108119;1119.317017:121.802902;1120.509277:69.353142;1122.307373:2.656713;1123.615845:6.962000;1124.326172:2.971273;1125.467041:6.552139;1126.704468:22.162647;1127.611572:61.035873;1128.631348:17.246387;1129.519897:15.741142;1133.684082:3.931945;1136.337402:9.022869;1137.696533:13.108788;1138.659668:12.340773;1139.562134:6.178851;1144.098267:22.034798;1146.544312:11.868593;1149.281738:19.254013;1151.326294:31.879799;1152.234619:5.220381;1154.266357:102.956070;1155.712524:51.562366;1156.606812:33.257694;1158.216553:7.032396;1162.519897:6.382627;1164.563110:1.094425;1165.363037:3.940457;1166.629883:29.609102;1167.378296:35.449104;1169.220703:41.219803;1170.325684:22.420496;1172.399536:293.539917;1173.380127:101.494263;1174.279053:9.931173;1178.696289:8.723358;1181.421021:4.238530;1190.368042:301.993439;1191.473511:757.853394;1192.454590:274.651184;1193.155273:14.667293;1198.364380:1.936964;1199.149658:5.023824;1200.965576:1.784602;1205.643188:52.082390;1208.918091:18.258869;1213.761963:4.253032;1214.572510:16.779848;1215.402710:3.864663;1217.238525:12.336935;1223.812012:3.236175;1224.736084:4.446353;1226.616211:19.375689;1228.294678:3.084185;1229.135132:5.374313;1230.744751:5.920449;1231.450439:6.109221;1233.541870:3.928753;1237.142700:3.695366;1239.778320:15.212591;1241.428589:51.327785;1242.603394:20.564140;1243.826172:5.750034;1249.720947:5.952960;1250.680908:4.951524;1253.561890:3.407656;1256.649536:1.979903;1259.197266:145.335358;1260.031128:17.104628;1260.650635:11.347934;1262.692383:10.603176;1263.818359:9.955406;1267.717285:2.071668;1269.407715:2.660241;1270.406250:6.261440;1277.253662:166.012421;1278.536865:88.384178;1279.569092:4.182921;1280.639893:12.693683;1281.757080:32.080242;1285.239624:18.666552;1285.853149:1.168572;1287.993896:15.721422;1288.884521:7.222568;1289.642212:19.486639;1291.822266:5.868732;1294.781372:33.283390;1295.757324:6.549904;1296.592773:10.694402;1298.491089:45.541962;1299.353760:29.731981;1300.371948:12.320871;1302.936523:15.345958;1303.938965:9.683271;1306.554565:161.443283;1307.521729:122.924423;1308.180908:2.345035;1311.209106:2.917070;1312.297729:54.099438;1313.388550:65.652885;1315.641113:8.838790;1316.688354:6.263451;1318.843506:17.979279;1323.372437:5.236959;1327.708984:29.739098;1328.394653:3.262248;1330.326782:82.029991;1331.365112:45.539597;1340.049927:5.806835;1340.807983:8.932220;1342.391724:3.411675;1345.632690:55.687668;1346.715454:39.691803;1348.319336:73.651413;1349.068848:43.869858;1349.803223:20.080372;1351.747070:6.190201;1354.340942:6.573697;1355.222778:11.552413;1356.462036:18.910576;1361.678223:7.807085;1363.541748:952.613892;1364.678345:377.923950;1370.736572:11.532435;1371.587769:15.245440;1372.715698:2.217935;1373.702515:18.248135;1376.031372:14.299359;1381.615112:4.656536;1382.986816:2.685962;1390.837769:15.506762;1391.541138:9.292782;1395.145874:10.974713;1408.527100:23.743408;1409.434448:12.054468;1416.421875:9.006732;1421.586548:8.945581;1425.598877:61.589577;1426.549194:27.886532;1428.401855:6.388432;1429.901733:3.320348;1433.128418:30.681986;1433.891968:9.279342;1437.509644:14.370709;1439.898682:7.308615;1443.460571:150.656891;1444.644653:44.598263;1448.410034:3.875106;1450.821899:2.399039;1452.293091:8.425796;1455.716553:5.403679;1461.479980:224.135635;1462.352173:103.046066;1463.460938:7.004594;1468.606323:5.683524;1469.336670:37.897030;1470.626343:10.847969;1474.402466:3.413245;1477.389526:3.725972;1478.073853:8.795669;1484.711426:17.292351;1485.752319:14.534218;1486.898193:16.326588;1492.260376:7.281396;1492.932861:5.523360;1494.065552:4.399163;1495.265503:16.447903;1503.872314:2.786598;1508.825806:2.225055;1510.614136:240.094635;1511.693115:79.985481;1516.638794:9.632721;1526.601929:12.610979;1527.404541:7.291560;1533.964233:7.974834;1536.188843:6.263948;1540.243164:3.733452;1541.497681:4.747337;1544.046875:18.565285;1544.941895:2.343512;1545.768677:6.547971;1547.140991:13.864877;1559.415894:18.763100;1568.263794:7.275522;1569.483887:5.987229;1571.599609:9.165901;1572.616577:49.792252;1573.633301:6.704363;1577.103882:7.855889;1579.654663:15.046736;1586.627319:4.794915;1590.526001:65.342926;1591.485962:44.237915;1597.729370:81.960213;1598.824097:69.959351;1608.538086:184.265015;1609.462158:41.348549;1610.105835:34.686352;1613.638916:11.629824;1616.565063:2.133877;1620.355469:10.171824;1625.926392:8.419968;1628.513062:10.672646;1630.957642:10.739971;1633.547607:4.763734;1634.626099:2.366058;1646.982666:21.267290;1648.979858:26.819963;1661.571411:3.035370;1665.685425:52.008629;1666.737915:30.760693;1667.640747:14.584148;1680.840088:11.928331;1683.041626:7.551940;1683.700439:8.729392;1690.709229:10.300920;1692.552856:7.314527;1704.977173:1.493919;1712.177002:9.987165;1712.842407:44.945171;1713.655884:31.200359;1714.962280:5.667592;1715.935669:34.665314;1725.400513:3.321785;1726.815674:3.146703;1735.107300:5.831262;1743.909424:61.078598;1744.900024:32.709404;1745.692993:4.817210;1760.944092:19.819504;1761.642456:143.824478;1762.603638:95.878731;1763.625732:31.149204;1779.663330:366.637024;1780.679077:305.051819;1798.638428:6.996554;1804.001221:26.200655;1810.885376:3.268579;1812.401245:3.987965;1813.833862:5.700371;1823.732666:11.294738;1829.043701:8.697560;1842.471436:10.691337;1858.751587:2.804528;1870.911499:1.789929;1875.026001:10.042640;1875.768311:1.793097;1895.849243:7.419658'

st.set_page_config(page_title="Spectra Viewer", page_icon=":glasses:", layout="wide")

# 1) Initialize Query Parameters
params = st.experimental_get_query_params()
query_peptide_sequence = params.get('sequence', ['[164.0700]FDSFGDLSSASAIM(16)GNPK'])[0]
query_mass_type = params.get('mass_type', ['monoisotopic'])[0]
query_fragment_types = params.get('fragment_types', ['1b;1y'])[0]
query_fragment_types = set('+' * int(f[0]) + f[1] for f in query_fragment_types.split(';'))
query_mass_tolerance_type = params.get('mass_tolerance_type', ['ppm'])[0]
query_mass_tolerance = params.get('mass_tolerance', [50])[0]
query_peak_assignment = params.get('peak_assignment', ['largest'])[0]
query_internal_fragments = params.get('internal_fragments', ['False'])[0]
query_min_intensity = params.get('min_intensity', [0.0])[0]
query_y_axis_scale = params.get('y_axis_scale', ['linear'])[0]
query_neutral_losses = params.get('neutral_losses', [''])[0]
query_neutral_losses = set(nl for nl in query_neutral_losses.split(';') if nl)
query_hide_unassigned_peaks = params.get('hide_unassigned_peaks', ['False'])[0]
query_spectra = params.get('spectra', [DEFAULT_SPECTRA])[0]
query_peak_picker = params.get('peak_picker', ['False'])[0]
query_peak_picker_min_intensity = params.get('peak_picker_min_intensity', [1.0])[0]
query_peak_picker_mass_tolerance = params.get('peak_picker_mass_tolerance', [0.02])[0]
query_isotopes = params.get('isotopes', [0])[0]
query_filter_missing_mono = params.get('filter_missing_mono', ['False'])[0]
query_filter_interrupted_iso = params.get('filter_interrupted_iso', ['False'])[0]

# 2) Validate Query Parameters

# Mass Type Validation
if query_mass_type not in VALID_MASS_TYPES:
    st.error("Invalid mass type.")
    st.stop()

# Fragment Types Validation
if not query_fragment_types.issubset(base_color_dict.keys()):
    st.error(f"Invalid fragment types: {query_fragment_types}")
    st.stop()

# Mass Tolerance Type Validation
if query_mass_tolerance_type not in VALID_MASS_TOLERANCE_TYPES:
    st.error("Invalid mass tolerance type.")
    st.stop()

# Mass Tolerance Validation
try:
    query_mass_tolerance = float(query_mass_tolerance)
    if query_mass_tolerance <= 0:
        st.error("Mass tolerance must be a positive number.")
        st.stop()
except ValueError:
    st.error("Mass tolerance must be a number.")
    st.stop()

# Peak Assignment Validation
if query_peak_assignment not in VALID_PEAK_ASSIGNMENTS:
    st.error("Invalid peak assignment.")
    st.stop()

# Internal Fragments Validation
if query_internal_fragments not in {'True', 'False'}:
    st.error("Internal fragments must be a boolean value.")
    st.stop()
query_internal_fragments = query_internal_fragments == 'True'

# Minimum Intensity Validation
try:
    query_min_intensity = float(query_min_intensity)
    if query_min_intensity < 0:
        st.error("Minimum intensity must be a non-negative number.")
        st.stop()
except ValueError:
    st.error("Minimum intensity must be a number.")
    st.stop()

# Y-axis Scale Validation
if query_y_axis_scale not in VALID_Y_AXIS_SCALES:
    st.error("Invalid y-axis scale.")
    st.stop()

# Hide Unassigned Peaks Validation
if query_hide_unassigned_peaks not in {'True', 'False'}:
    st.error("Hide unassigned peaks must be a boolean value.")
    st.stop()
query_hide_unassigned_peaks = query_hide_unassigned_peaks == 'True'

spectra_pairs = query_spectra.split(';')
for pair in spectra_pairs:
    if ':' not in pair:
        st.error(f"Invalid format for pair '{pair}' in spectra. Expected format: 'mass:intensity'")
        st.stop()
        continue

    mass, intensity = pair.split(':')

    try:
        # Validate mass
        mass = float(mass)
        if mass < 0:
            st.error(f"Mass '{mass}' must be a non-negative number.")
            st.stop()
    except ValueError:
        st.error(f"Mass '{mass}' must be a number.")
        st.stop()

    try:
        # Validate intensity
        intensity = float(intensity)
        if intensity < 0:
            st.error(f"Intensity '{intensity}' must be a non-negative number.")
            st.stop()
    except ValueError:
        st.error(f"Intensity '{intensity}' must be a number.")
        st.stop()

query_spectra = '\n'.join([f'{pair.split(":")[0]} {pair.split(":")[1]}' for pair in query_spectra.split(';')])

# Peak Picker Validation
if query_peak_picker not in {'True', 'False'}:
    st.error("Peak picker must be a boolean value.")
    st.stop()
query_peak_picker = query_peak_picker == 'True'

# Peak Picker Minimum Intensity Validation
try:
    query_peak_picker_min_intensity = float(query_peak_picker_min_intensity)
    if query_peak_picker_min_intensity < 0:
        st.error("Peak picker minimum intensity must be a non-negative number.")
        st.stop()
except ValueError:
    st.error("Peak picker minimum intensity must be a number.")
    st.stop()

# Peak Picker Mass Tolerance Validation
try:
    query_peak_picker_mass_tolerance = float(query_peak_picker_mass_tolerance)
    if query_peak_picker_mass_tolerance <= 0:
        st.error("Peak picker mass tolerance must be a positive number.")
        st.stop()
except ValueError:
    st.error("Peak picker mass tolerance must be a number.")
    st.stop()

# Neutral losses validation
if query_neutral_losses:
    for nl in query_neutral_losses:
        if nl not in NEUTRAL_LOSSES:
            st.error(f"Invalid neutral loss: {nl}")
            st.stop()

# Isotopes validation
try:
    query_isotopes = int(query_isotopes)
    if query_isotopes < 0:
        st.error("Isotopes must be a non-negative integer.")
        st.stop()
except ValueError:
    st.error("Isotopes must be an integer.")
    st.stop()

# Filter interrupted isotopes validation
if query_filter_missing_mono not in {'True', 'False'}:
    st.error("Filter interrupted isotopes must be a boolean value.")
    st.stop()
query_filter_missing_mono = query_filter_missing_mono == 'True'

if query_filter_interrupted_iso not in {'True', 'False'}:
    st.error("Filter interrupted isotopes must be a boolean value.")
    st.stop()
query_filter_interrupted_iso = query_filter_interrupted_iso == 'True'

# 3) Get User inputs
with st.sidebar:
    st.title('Spectra Viewer')

    sequence = st.text_input(label='Sequence', value=query_peptide_sequence)
    unmodified_sequence = peptacular.sequence.strip_modifications(sequence)

    st.write('Ions')
    for ion in IONS:
        cols = st.columns(MAX_CHARGE)
        for i, col in enumerate(cols, 1):
            label = '+' * i + ion
            col.checkbox(label, value=label in query_fragment_types, key=label)

    internal_fragments = st.checkbox(label='Internal Fragments', value=query_internal_fragments)
    num_isotopes = st.number_input(label='Isotopes', value=query_isotopes, min_value=0, max_value=5)
    filter_missing_mono = st.checkbox(label='Filter missing mono peaks',
                                      value=query_filter_missing_mono)
    filter_interrupted_iso = st.checkbox(label='Filter interrupted isotopes',
                                         value=query_filter_interrupted_iso)

    st.write('Neutral Losses')
    cols = st.columns(len(NEUTRAL_LOSSES))
    losses = [0.0]
    nls = []
    for i, (nl, mass) in enumerate(NEUTRAL_LOSSES.items()):
        l = cols[i].checkbox(nl, value=nl in query_neutral_losses, key=nl)
        if l:
            losses.append(mass)
            nls.append(nl)

    mass_type = st.radio('Mass Type', ['monoisotopic', 'average'], horizontal=True,
                         index=0 if query_mass_type == 'monoisotopic' else 1)
    mass_tolerance_type = st.radio('Mass Tolerance Type', ['ppm', 'th'], horizontal=True,
                                   index=0 if query_mass_tolerance_type == 'ppm' else 1)
    mass_tolerance = st.number_input('Mass Tolerance', value=query_mass_tolerance)
    peak_assignment = st.radio('Peak Assignment', ['largest', 'closest'], horizontal=True,
                               index=0 if query_peak_assignment == 'largest' else 1)

    min_intensity = st.number_input('Min Intensity', value=query_min_intensity)
    y_axis_scale = st.radio('Y Axis Scale', ['linear', 'log'], horizontal=True,
                            index=0 if query_y_axis_scale == 'linear' else 1)
    hide_unassigned_peaks = st.checkbox('Hide Unassigned Peaks', value=query_hide_unassigned_peaks)
    peak_picker = st.checkbox('Peak Picker', value=query_peak_picker)

    peak_picker_min_intensity = query_peak_picker_min_intensity
    peak_picker_mass_tolerance = query_peak_picker_mass_tolerance
    if peak_picker:
        c1, c2 = st.columns(2)
        peak_picker_min_intensity = c1.number_input('Peak Picker Min Intensity', value=query_peak_picker_min_intensity)
        peak_picker_mass_tolerance = c2.number_input('Peak Picker Mass Tolerance',
                                                     value=query_peak_picker_mass_tolerance)

    spectra = st.text_area(label='Spectra', value=query_spectra)

if mass_tolerance_type == 'th' and mass_tolerance > 1:
    st.error("Mass tolerance in th must be less than or equal to 1, when type is th.")
    st.stop()

# retrieve ion selections
ion_types = []
charges = []
for ion in IONS:
    for i in range(1, MAX_CHARGE + 1):
        label = '+' * i + ion
        if st.session_state.get(label, False):
            ion_types.append(ion)
            charges.append(i)


def generate_app_url():
    base_url = BASE_URL
    params = {
        'sequence': urllib.parse.quote(sequence),
        'mass_type': mass_type,
        'fragment_types': ';'.join([str(charge) + ion for ion, charge in zip(ion_types, charges)]),
        'mass_tolerance_type': str(mass_tolerance_type),
        'peak_assignment': str(peak_assignment),
        'internal_fragments': internal_fragments,
        'min_intensity': str(min_intensity),
        'y_axis_scale': str(y_axis_scale),
        'neutral_losses': ';'.join(nls),
        'hide_unassigned_peaks': hide_unassigned_peaks,
        'spectra': ';'.join([f'{line.split(" ")[0]}:{line.split(" ")[1]}' for line in spectra.split('\n')]),
        'peak_picker': peak_picker,
        'peak_picker_min_intensity': str(peak_picker_min_intensity),
        'peak_picker_mass_tolerance': str(peak_picker_mass_tolerance),
        'isotopes': str(num_isotopes),
        'filter_missing_mono': str(filter_missing_mono),
        'filter_interrupted_iso': str(filter_interrupted_iso),
    }
    query_string = '&'.join([f'{key}={value}' for key, value in params.items() if value is not None])
    return f'{base_url}?{query_string}'


# Show Analysis URL
url = generate_app_url()
st.write(f'##### [Analysis URL]({url}) (copy me and send to your friends!)')

# Show Sequence Info
sequence_charge = 2
st.header(sequence)
c1, c2, c3, c4 = st.columns(4)
c1.metric('Mass', round(calculate_mass(sequence), 4))
c2.metric('m/z', round(calculate_mz(sequence, sequence_charge), 4))
c3.metric('Charge', sequence_charge)
c4.metric('Length', len(unmodified_sequence))

fragments = []
for ion, charge in zip(ion_types, charges):
    fragments.extend(build_fragments(sequence=sequence,
                                     ion_types=ion,
                                     charges=charge,
                                     monoisotopic=(mass_type == 'monoisotopic'),
                                     internal=internal_fragments,
                                     isotopes=list(range(num_isotopes + 1)),
                                     losses=losses))

frag_df = pd.DataFrame([fragment.to_dict() for fragment in fragments])

#  for all ends that are negative add seq_len
frag_df.loc[frag_df['end'] < 0, 'end'] += len(unmodified_sequence)
frag_df.loc[frag_df['start'] < 0, 'start'] += len(unmodified_sequence)

# where end is None, set to seq_len
frag_df.loc[frag_df['end'].isna(), 'end'] = len(unmodified_sequence)

if spectra:

    mzs, ints = [], []

    for line in spectra.split('\n'):
        mz, intensity = line.split(' ')

        mz = float(mz)
        intensity = float(intensity)

        if intensity <= min_intensity:
            continue

        mzs.append(mz)
        ints.append(intensity)

    if peak_picker:
        peaks = [(mz, intensity) for mz, intensity in zip(mzs, ints)]
        deconvoluted_peaks, _ = ms_deisotope.deconvolute_peaks(peaks, averagine=ms_deisotope.peptide,
                                                               scorer=ms_deisotope.MSDeconVFitter(
                                                                   query_peak_picker_min_intensity,
                                                                   query_peak_picker_mass_tolerance))
        mzs = [peak.mz for peak in deconvoluted_peaks]
        ints = [peak.intensity for peak in deconvoluted_peaks]

    max_spectra_mz = max(mzs)

    fragment_matches = compute_fragment_matches(fragments, mzs, ints, mass_tolerance, mass_tolerance_type)
    fragment_matches.sort(key=lambda x: abs(x.error), reverse=True)
    fragment_matches = {fm.mz: fm for fm in fragment_matches}  # keep the best error for each fragment

    match_data = {'sequence': [], 'charge': [], 'ion_type': [], 'number': [], 'internal': [], 'parent_number': [],
                  'monoisotopic': [], 'mz': [], 'intensity': [], 'error': [], 'abs_error': [], 'theo_mz': [],
                  'label': [],
                  'isotope': [], 'loss': []}
    data = {'sequence': [], 'charge': [], 'ion_type': [], 'number': [], 'internal': [], 'parent_number': [],
            'monoisotopic': [], 'mz': [], 'intensity': [], 'error': [], 'abs_error': [], 'theo_mz': [], 'label': [],
            'isotope': [], 'loss': []}

    for mz, i in zip(mzs, ints):
        fm = fragment_matches.get(mz, None)

        if fm:
            match_data['sequence'].append(fm.fragment.sequence)
            match_data['charge'].append(fm.fragment.charge)
            match_data['ion_type'].append(fm.fragment.ion_type)
            match_data['number'].append(fm.fragment.number)
            match_data['internal'].append(fm.fragment.internal)
            match_data['parent_number'].append(fm.fragment.parent_number)
            match_data['monoisotopic'].append(fm.fragment.monoisotopic)
            match_data['mz'].append(mz)
            match_data['intensity'].append(i)
            match_data['error'].append(fm.error)
            match_data['abs_error'].append(abs(fm.error))
            match_data['theo_mz'].append(fm.fragment.mz)
            match_data['label'].append(fm.fragment.label)
            match_data['isotope'].append(fm.fragment.isotope)
            match_data['loss'].append(fm.fragment.loss)

        data['sequence'].append('')
        data['charge'].append(0)
        data['ion_type'].append('')
        data['number'].append(0)
        data['internal'].append(False)
        data['parent_number'].append(0)
        data['monoisotopic'].append(True)
        data['mz'].append(mz)
        data['intensity'].append(i)
        data['error'].append(0)
        data['abs_error'].append(0)
        data['theo_mz'].append(0)
        data['label'].append('')
        data['isotope'].append(0)
        data['loss'].append(0.0)

    spectra_df = pd.DataFrame(data)
    spectra_df['matched'] = False
    match_df = pd.DataFrame(match_data)
    match_df['matched'] = True

    # for keep only the lowest abs_error for ion_type, charge, num
    if peak_assignment == 'most intense':
        match_df.sort_values(by='intensity', inplace=True, ascending=False)
        match_df.drop_duplicates(subset=['theo_mz'], inplace=True)

    else:
        match_df.sort_values(by='abs_error', inplace=True)
        match_df.drop_duplicates(subset=['theo_mz'], inplace=True)

    if filter_missing_mono:
        # remove peaks that skip isotopes
        mono_labels = set(match_df[match_df['isotope'] == 0]['label'].unique())
        labels_to_remove = set()
        for label in match_df[match_df['isotope'] != 0]['label'].unique():
            mono_label = label.replace('*', '')
            if mono_label not in mono_labels:
                labels_to_remove.add(label)

        match_df = match_df[~match_df['label'].isin(labels_to_remove)]

    if filter_interrupted_iso:
        # remove peaks any peaks after a missing isotope
        labels_to_remove = set()

        # get min isotope label for each mono label
        labels = set(match_df['label'].unique())
        mono_iso_labels = set([label.replace('*', '') for label in labels])
        min_iso_labels = set()
        for mono_label in mono_iso_labels:
            for i in range(1, num_isotopes+1):
                label = mono_label + '*' * i
                if label in labels:
                    min_iso_labels.add(label)
                    break

        for min_label in min_iso_labels:
            break_flag = False
            iso_cnt = min_label.count('*')
            for i in range(iso_cnt, num_isotopes+1):
                label = min_label + '*' * i
                if label not in labels:
                    break_flag = True
                if break_flag:
                    labels_to_remove.add(label)
        match_df = match_df[~match_df['label'].isin(labels_to_remove)]





    spectra_df = spectra_df[~spectra_df['mz'].isin(match_df['mz'])]
    spectra_df = pd.concat([spectra_df, match_df])

    if hide_unassigned_peaks:
        spectra_df = spectra_df[spectra_df['matched']]

    spectra_df['ion_color_type'] = spectra_df['ion_type']
    spectra_df.loc[spectra_df['internal'], 'ion_color_type'] = 'i'


    def create_labels(row):
        if row['ion_type'] != '':
            charge_str = '+' * int(row['charge'])
            ion_type_str = row['ion_type']
            parent_number_str = str(int(row['parent_number']))
            internal_str = 'i' if row['internal'] else ''

            color_label = f"{charge_str}{ion_type_str}"
            ion_label = f"{charge_str}{ion_type_str}{parent_number_str}{internal_str}"
        else:
            color_label = 'unassigned'
            ion_label = 'unassigned'

        return ion_label, color_label


    # List comprehension to create ion and color labels
    labels = [create_labels(row) for _, row in spectra_df.iterrows()]
    spectra_df['ion_label'], spectra_df['color_label'] = zip(*labels)

    # Assigning colors based on color labels
    spectra_df['color'] = [base_color_dict[label] for label in spectra_df['color_label']]

    cmap = mpl.colormaps.get_cmap('Blues')

    st.caption('Sequence coverage')
    for ion, charge in zip(ion_types, charges):
        cov_arr = [0] * len(unmodified_sequence)
        tmp_df = spectra_df[(spectra_df['ion_type'] == ion) & (spectra_df['charge'] == charge)]
        nums = tmp_df['parent_number'].unique()

        if ion in 'abc':
            for num in nums:
                cov_arr[num - 1] = 1
        else:
            for num in nums:
                cov_arr[len(unmodified_sequence) - (num - 1) - 1] = 1

        if len(cov_arr) > 0:
            c = base_color_dict['+' * charge + ion]
            s = coverage_string(cov_arr, unmodified_sequence, c)

            # center text
            ion_span = f'<span style="color:{c}">{ion}<sup>+{charge}</sup></span>'
            st.markdown(f'{ion_span} {s}', unsafe_allow_html=True)

    fig = generate_annonated_spectra_plotly(spectra_df, scale=y_axis_scale)
    st.plotly_chart(fig, use_container_width=True)

    dfs = []
    combined_data = {'AA': list(unmodified_sequence)}
    for ion, charge in zip(ion_types, charges):
        data = {'AA': list(unmodified_sequence)}
        ion_df = frag_df[
            (frag_df['ion_type'] == ion) & (frag_df['charge'] == charge) &
            (frag_df['internal'] == False)]
        ion_df.sort_values(by=['number'], inplace=True)

        # keep only a single number
        ion_df.drop_duplicates(subset=['number'], inplace=True)

        frags = ion_df['mz'].tolist()

        if ion in 'xyz':
            frags = frags[::-1]

        data[ion] = frags

        combined_data['+' * charge + ion] = frags

        # Displaying the table
        df = pd.DataFrame(data)
        df['# (abc)'] = list(range(1, len(df) + 1))
        df['# (xyz)'] = list(range(1, len(df) + 1))[::-1]

        # reorder columns so that # is first # +1 is last and AA is in the middle
        df = df[
            ['AA'] + ['# (abc)'] + [col for col in df.columns if col not in ['AA', '# (abc)', '# (xyz)']] + ['# (xyz)']]
        dfs.append(df)

    combined_df = pd.DataFrame(combined_data)
    # sort columns based on alphabetical order
    combined_df = combined_df.reindex(sorted(combined_df.columns), axis=1)

    styled_dfs = []


    def color_by_ion_type(col):
        ion_type = col.name[-1]
        color = base_color_dict.get(ion_type, 'grey')  # get color or default to grey if not found
        return ['color: %s' % color] * len(col)


    for df in dfs:
        styled_df = df.style.apply(color_by_ion_type)

        # Set table styles with increased horizontal padding for more space between columns,
        # centered text, and no borders
        styles = [
            dict(selector="td", props=[("padding", "2px 2px"), ("text-align", "center"), ("border", "none")]),
            dict(selector="th", props=[("padding", "2px 2px"), ("text-align", "center"), ("border", "none")])
        ]
        styled_df = styled_df.set_table_styles(styles)
        styled_dfs.append(styled_df)


    def highlight_cells(data):
        # Initialize empty DataFrame with same index and columns as original
        styled = pd.DataFrame('', index=data.index, columns=data.columns)

        # Iterate over cells and update `styled` based on cell position
        for row in data.index:
            for col in data.columns:
                if col == 'AA' or col == '# (abc)' or col == '# (xyz)':
                    styled.loc[
                        row, col] = f'background-color: gainsboro; color: black; text-align: center; font-weight: bold;'
                    continue

                ion = col[-1]
                if ion in 'abc':
                    ion_number = row + 1
                else:
                    ion_number = len(unmodified_sequence) - row
                label = col + str(ion_number)
                if label in accepted_normal_ions:
                    styled.loc[
                        row, col] = f'background-color: {base_color_dict[col]}; color: white; text-align: center;'
                elif label in accepted_internal_ions:
                    styled.loc[
                        row, col] = f'background-color: {base_color_dict[col]}; color: magenta; text-align: center; font-style: italic;'
                else:
                    styled.loc[row, col] = f'background-color: white; color: black; text-align: center;'

        return styled


    matched_ions = spectra_df[spectra_df['ion_type'] != '']
    accepted_normal_ions = matched_ions[matched_ions['internal'] == False]['ion_label'].tolist()
    accepted_internal_ions = matched_ions[matched_ions['internal'] == True]['ion_label'].tolist()
    accepted_internal_ions = [ion[:-1] for ion in accepted_internal_ions]

    combined_df['# (abc)'] = list(range(1, len(unmodified_sequence) + 1))
    combined_df['# (xyz)'] = list(range(1, len(unmodified_sequence) + 1))[::-1]

    # reorder columns so that # is first # +1 is last and AA is in the middle
    combined_cols = combined_df.columns.tolist()
    combined_cols.remove('# (abc)')
    combined_cols.remove('# (xyz)')
    combined_cols.remove('AA')
    forward_cols = [col for col in combined_cols if 'a' in col or 'b' in col or 'c' in col]
    reverse_cols = [col for col in combined_cols if 'x' in col or 'y' in col or 'z' in col]

    # sort
    forward_cols.sort()
    reverse_cols.sort(reverse=True)

    new_cols = ['# (abc)'] + forward_cols + ['AA'] + reverse_cols + ['# (xyz)']
    combined_df = combined_df[new_cols]

    combined_df = combined_df.style.format(precision=4).apply(highlight_cells, axis=None)

    st.markdown('---')

    st.subheader('Fragment Ions')

    st.dataframe(combined_df, height=int(35.2 * (len(df) + 1)), hide_index=True)

    st.markdown('---')

    st.subheader('Stats')

    # calculate percerntage if intensity accounted for
    total_intensity = spectra_df['intensity'].sum()
    c1, c2, c3, c4 = st.columns(4)
    c1.metric(label='Total Intensity', value=round(total_intensity, 1))
    c2.metric(label='Matched Intensity', value=round(match_df['intensity'].sum(), 1))
    c3.metric(label='Unmatched Intensity', value=round(spectra_df['intensity'].sum() - match_df['intensity'].sum(), 1))
    c4.metric(label='Matched Intensity %', value=round(match_df['intensity'].sum() / total_intensity * 100, 2))

    c1, c2 = st.columns(2)
    hs = hyper_score(fragments, spectra_df['mz'].tolist(), spectra_df['intensity'].tolist(), mass_tolerance,
                     mass_tolerance_type)
    c1.metric(f'Hyperscore', round(hs, 5))
    bs = binomial_score(fragments, spectra_df['mz'].tolist(), spectra_df['intensity'].tolist(), mass_tolerance,
                        mass_tolerance_type)
    c2.metric(f'Binomial Score', round(bs, 5))

    # isotopes

    st.markdown('---')

    st.subheader('Isotopes')

    iso_data = []
    for i, row in spectra_df.iterrows():
        mz = row['mz']
        intensity = row['intensity']
        charge = row['charge']

        if charge == 0:
            continue

        isotopic_pattern = ms_deisotope.peptide.isotopic_cluster(mz, charge)
        for peak in isotopic_pattern:
            iso_data.append([i, peak.mz, peak.intensity])

with st.expander('Fragments'):
    st.dataframe(frag_df)

with st.expander('Peaks'):
    st.dataframe(spectra_df)
