import numpy as np
CMAES=[[0.24711362652718327, 0.013408866996237508, 101.33707596282798, 806.5254169541688, 7.743919229112687e-05, 0.0018886700374539794, -0.00033587180721899625, 7.896536259502126e-11, 8.940643685592182, 4.431191900797479, 5.118761570372981, 0.5816941444272093, 0.13997198513370698, 1.0144430613416717] ,
[0.2437136033540455, 0.041916959869871935, 136.91071718538336, 535.9775743351688, 7.719208829248059e-05, 0.002419672072987142, -0.00039798693040548304, 7.389022703687543e-11, 8.940524733361801, 4.357762581087894, 4.974566671527844, 0.6065682978303066, 0.17229158014845414, 0.8388595014110765] ,
[0.2441617016505956, 0.039364964752397265, 129.63299534053405, 566.3832308179987, 7.675869934460854e-05, 0.002110040474056385, -0.00037716035453121306, 7.418176499131204e-11, 8.940535789828703, 4.36718783811148, 4.994014066060285, 0.6004786421049727, 0.170416862420205, 0.8556047227891193] ,
[0.24371075391162972, 0.03620334850289272, 116.21904685871607, 593.2387950043394, 7.651157477868371e-05, 0.0019226672079532131, -0.0003703299642881259, 7.450073156164885e-11, 8.940537383480539, 4.374221715529037, 5.016020419361728, 0.5885550477158792, 0.1668535738891655, 0.8426632305034059] ,
[0.24256288257254335, 0.04189344127166237, 132.7412818851098, 535.0738541808678, 7.646463290727133e-05, 0.001833433242548442, -0.000378944972896339, 7.347777107474285e-11, 8.940517774295671, 4.357590583645216, 4.976638544716193, 0.5983044930766688, 0.1746954789216934, 0.7998841372117569] ,
[0.24154243644795073, 0.044598443408388325, 140.0629741223467, 502.223965116431, 7.647752034272181e-05, 0.0017553595146381565, -0.00038470953946489507, 7.295719995663411e-11, 8.940504028716738, 4.34871686729649, 4.958004115599524, 0.6021762793401005, 0.1779650906761758, 0.7678043847822007] ,
[0.2412122267816212, 0.04458757671090765, 141.4012320566926, 501.7238889858575, 7.640603385071265e-05, 0.0015654836888094673, -0.0003793449051687544, 7.27039628798611e-11, 8.940489722797832, 4.34916320517287, 4.9575883041467455, 0.601251596869569, 0.17836430477179255, 0.7567570859273287] ,
[0.24098843491379984, 0.04487340955980655, 138.41393346805972, 497.01721973880177, 7.615863211092515e-05, 0.0016498905042449871, -0.00038402749348127006, 7.233525893017939e-11, 8.940503541631484, 4.346649914048442, 4.955558867141757, 0.5992196347402237, 0.1784266288141846, 0.7529279195536173] ,
[0.24049901473118884, 0.04458198556041415, 134.5525370361804, 499.088790616767, 7.606658289450992e-05, 0.0014113104215398783, -0.00037680776786745696, 7.232621690026429e-11, 8.940505024458282, 4.346287189609495, 4.958385759285501, 0.5955919554239513, 0.17945675354461255, 0.7394138501735974] ,
[0.24116928285797873, 0.044189985533550226, 135.53440184454033, 510.6019261842285, 7.57155256924869e-05, 0.001459719577206274, -0.00037125091413634306, 7.200978526571967e-11, 8.940522370065947, 4.349034891035952, 4.961815285416956, 0.5973834903322666, 0.17881384642004294, 0.7614168539552315]]
MCMC_stds=[[0.0001253679955162606, 0.002895613641839063, 3.1252180810254915, 18.849961717990514, 1.1022391327601339e-07, 0.0003696160080219368, 1.7544696385022283e-05, 6.284549666150088e-13, 2.8571153334294526e-05, 0.005369092322632974, 0.0068750799126064575, 0.008852803296053835, 0.0097099597803178] ,
[0.00012506228879330568, 0.0010276392858232367, 3.3943252838317504, 14.377334012599432, 7.078816122315669e-08, 0.0002752662906474063, 1.3325585826216273e-05, 4.1096820214702057e-13, 2.5671795439058407e-05, 0.0041409393083647475, 0.006995619004994169, 0.007825119971865364, 0.008051230260239576] ,
[0.0001277966457738442, 0.0011789372602185792, 3.4907074983816617, 16.114786931325927, 7.660293559816026e-08, 0.00029518925264785756, 1.4256252769809632e-05, 4.552353445751688e-13, 2.5901236481383688e-05, 0.004601369053550331, 0.007633827472365744, 0.007900005995236796, 0.00823502345589884] ,
[0.00012177165361576796, 0.0013035369251618006, 3.162264122480119, 17.0858239827902, 7.734136760485048e-08, 0.00030233056231061253, 1.4681624985082344e-05, 4.715177523499399e-13, 2.5567829876817843e-05, 0.0048517958879027425, 0.007792526361308455, 0.0073837478201809, 0.007782089407920716] ,
[0.00011794189454181035, 0.0010171160644777064, 3.5377134757928794, 14.147085494651037, 6.752443527565798e-08, 0.00026874390657436256, 1.3038695116427224e-05, 3.9889027950859027e-13, 2.4703934789704675e-05, 0.004030153492466415, 0.006884737788910197, 0.0074037022490257, 0.007310121253715233] ,
[0.0001169269792943754, 0.0008305894462864009, 3.498925512817276, 11.756176892467957, 6.147869780604016e-08, 0.00023581717944887856, 1.152051571538602e-05, 3.304251109582605e-13, 2.3829953668351576e-05, 0.0033652974648172244, 0.005933245757118431, 0.00730104551513345, 0.0072864141519479375] ,
[0.00011423531675517217, 0.0008293563114805797, 3.702132264532235, 11.746678209357857, 6.167899181528225e-08, 0.0002413333228332799, 1.178377354651386e-05, 3.3514455663647447e-13, 2.2288711888353114e-05, 0.003330585722229247, 0.005992564695067563, 0.007369641616946625, 0.007469450401646866] ,
[0.0001168396518608515, 0.0008096553804910855, 3.48603732823611, 11.600422482145332, 6.160591524935496e-08, 0.00023454767024456153, 1.139517656602675e-05, 3.2591578422565465e-13, 2.3562777929477305e-05, 0.003279634460831279, 0.005872824190067626, 0.007174445434441141, 0.007220778220217468] ,
[0.00011317078192679746, 0.0008244718426788368, 3.400899469677774, 11.730921478807657, 5.926114780019426e-08, 0.00023264328010241846, 1.1367330968500481e-05, 3.315119467205256e-13, 2.268304504198771e-05, 0.0033155430932922613, 0.005963311325516575, 0.006593645285944483, 0.0072115475654292325] ,
[0.00011587817787598058, 0.0008593300890708407, 3.6964352783051333, 12.19608990289096, 6.099952331761285e-08, 0.00023640075188081493, 1.152025125816077e-05, 3.3850644186196096e-13, 2.332493512773059e-05, 0.003428292234150651, 0.006154262515036874, 0.007187996019639223, 0.007206067242433421]]
MCMC_means=[[0.2471220201635999, 0.012941773545923779, 101.68033456963573, 808.0261899777262, 7.745155375130099e-05, 0.0018551113728669923, -0.0003343653396312254, 7.901419820933974e-11, 8.940643865168989, 4.431632701288219, 5.119021325341933, 0.14029961547418038, 1.0158090565691607] ,
[0.24371349835072367, 0.041936858006944175, 136.94984672076333, 535.7422688370601, 7.718805149832768e-05, 0.0024279435479043158, -0.00039823680382825065, 7.389241287105244e-11, 8.940525430018987, 4.357671632392124, 4.9744574866828914, 0.1719878792743887, 0.8398172908978135] ,
[0.2441499221987774, 0.039369959271601536, 129.77358193064788, 566.5295782931634, 7.675421174480399e-05, 0.002127678396729643, -0.0003779498859312667, 7.418307979000883e-11, 8.940536126112725, 4.367209118752778, 4.994051645495269, 0.16904023339328733, 0.8564031119763135] ,
[0.24371695030159093, 0.03619432973609918, 116.48110064000034, 593.821847467774, 7.651154072307574e-05, 0.0019124162060485325, -0.000369709520767548, 7.452386397100383e-11, 8.940538588093382, 4.374369577143438, 5.016052992819522, 0.16628604694152554, 0.8436276972333014] ,
[0.2425603899035272, 0.041868128343087674, 133.6390372449569, 536.3073477053446, 7.646335295865882e-05, 0.001825409568032108, -0.00037844786085061607, 7.351683588062141e-11, 8.940517797767612, 4.357921017979627, 4.9766052169346935, 0.1749078120698646, 0.8013745538764748] ,
[0.24153910732168946, 0.044603891113748496, 140.47824310452017, 502.33792593039504, 7.647385675474214e-05, 0.0017712378385968135, -0.00038532095477947643, 7.294913544203554e-11, 8.940503412695545, 4.3487637277733215, 4.957850585551708, 0.17757928994786973, 0.7684043096153139] ,
[0.24121750374766188, 0.044614142937894535, 141.9583866287712, 501.5623795195233, 7.640462015541115e-05, 0.001573448821389732, -0.00037962537521388234, 7.269210703911166e-11, 8.940489587237728, 4.349140090733003, 4.957190783420997, 0.1783266391512056, 0.7578279816855213] ,
[0.24097546764255365, 0.04481883685433181, 138.24756157958754, 497.36547270896585, 7.616304688764342e-05, 0.0016539056339240717, -0.00038428564081840794, 7.231888565513295e-11, 8.940504264967933, 4.346756890589177, 4.955871426388552, 0.177330114686983, 0.753551341034815] ,
[0.24049810162497975, 0.04465310233293925, 134.771285037017, 498.27830472943623, 7.606134642330584e-05, 0.0014339297596182342, -0.00037788902935737126, 7.230414492381035e-11, 8.940506811448596, 4.346040609502586, 4.957886372731034, 0.17878627096840077, 0.7401488125358374] ,
[0.24117167531813616, 0.044130560446392623, 135.64446572670033, 511.7165123899569, 7.57152316894211e-05, 0.0014490868743300348, -0.0003706986852740553, 7.20449700317367e-11, 8.94052256122038, 4.3493190681991125, 4.962286820156788, 0.17813247128394036, 0.7624356927435779]]
i_0=6.038052132625424e-06
desired_MCMC_params=list(range(0, len(MCMC_means[0])))
desired_CMAES_params=list(range(0, len(CMAES[0])))
desired_CMAES_params=(np.delete(desired_CMAES_params, 11))
diffs=np.zeros((len(CMAES), len(desired_CMAES_params)))
def list_indexer(list, index):
    return [list[x] for x in index]
for i in range(0, len(MCMC_means)):
    diffs[i, :]=np.divide(abs(np.subtract(list_indexer(CMAES[i],desired_CMAES_params), list_indexer(MCMC_means[i],desired_MCMC_params))), MCMC_stds[i])
coeffs=np.divide(MCMC_stds, MCMC_means)
print(np.mean(np.mean(coeffs, axis=0))*100)

#print(np.mean(np.mean(diffs, axis=0)))
