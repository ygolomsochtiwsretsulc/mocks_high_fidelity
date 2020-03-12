#!/bin/bash 

cd /data40s/erosim/eRASS/eSASS-v11-erass8-jc
pyCONDA
source /home/erosita/sw/sass-setup.sh eSASSusers_190925

copy healpix_radius_all.dat into healpix_radius.datb and select the fields of interest !

nohup sh cal_evt.sh > cal_evt.log &

cd 000
../eRASS_det_pipe_V11.py evt_000.fits --bands=c 
cd ../001
../eRASS_det_pipe_V11.py evt_001.fits --bands=c 
cd ../002
../eRASS_det_pipe_V11.py evt_002.fits --bands=c 

nohup sh all_sky_simulate.sh > all_sky_simulate.log & 

cd 103
../eRASS_det_pipe_V11.py evt_103.fits --bands=c

cd ../104
../eRASS_det_pipe_V11.py evt_104.fits --bands=c

cd ../105
../eRASS_det_pipe_V11.py evt_105.fits --bands=c

cd ../106
../eRASS_det_pipe_V11.py evt_106.fits --bands=c

cd ../107
../eRASS_det_pipe_V11.py evt_107.fits --bands=c

cd ../108
../eRASS_det_pipe_V11.py evt_108.fits --bands=c

cd ../109
../eRASS_det_pipe_V11.py evt_109.fits --bands=c

cd ../110
../eRASS_det_pipe_V11.py evt_110.fits --bands=c


cd 201
../eRASS_det_pipe_V11.py evt_201.fits --bands=c

cd ../202
../eRASS_det_pipe_V11.py evt_202.fits --bands=c

cd ../203
../eRASS_det_pipe_V11.py evt_203.fits --bands=c

cd ../204
../eRASS_det_pipe_V11.py evt_204.fits --bands=c

cd ../205
../eRASS_det_pipe_V11.py evt_205.fits --bands=c

cd ../206
../eRASS_det_pipe_V11.py evt_206.fits --bands=c

cd ../207
../eRASS_det_pipe_V11.py evt_207.fits --bands=c

cd ../208
../eRASS_det_pipe_V11.py evt_208.fits --bands=c

cd ../209
../eRASS_det_pipe_V11.py evt_209.fits --bands=c

cd ../210
../eRASS_det_pipe_V11.py evt_210.fits --bands=c

cd 211
../eRASS_det_pipe_V11.py evt_211.fits --bands=c

cd ../212
../eRASS_det_pipe_V11.py evt_212.fits --bands=c

cd 301
../eRASS_det_pipe_V11.py evt_301.fits --bands=c

cd ../302
../eRASS_det_pipe_V11.py evt_302.fits --bands=c

cd ../303
../eRASS_det_pipe_V11.py evt_303.fits --bands=c

cd ../304
../eRASS_det_pipe_V11.py evt_304.fits --bands=c

cd ../305
../eRASS_det_pipe_V11.py evt_305.fits --bands=c

cd ../306
../eRASS_det_pipe_V11.py evt_306.fits --bands=c

cd ../307
../eRASS_det_pipe_V11.py evt_307.fits --bands=c

cd ../308
../eRASS_det_pipe_V11.py evt_308.fits --bands=c

cd ../380
../eRASS_det_pipe_V11.py evt_380.fits --bands=c


# nohup sh sim_000.sh > sim_000.log & 
# nohup sh sim_001.sh > sim_001.log & 
# nohup sh sim_002.sh > sim_002.log & 
# nohup sh sim_003.sh > sim_003.log & 
# nohup sh sim_004.sh > sim_004.log & 
# nohup sh sim_005.sh > sim_005.log & 
# nohup sh sim_006.sh > sim_006.log & 
# nohup sh sim_007.sh > sim_007.log & 
# nohup sh sim_008.sh > sim_008.log & 
# nohup sh sim_009.sh > sim_009.log & 
# nohup sh sim_010.sh > sim_010.log & 
# nohup sh sim_011.sh > sim_011.log & 
# nohup sh sim_012.sh > sim_012.log & 
# nohup sh sim_013.sh > sim_013.log & 
# nohup sh sim_014.sh > sim_014.log & 
# nohup sh sim_015.sh > sim_015.log & 
# nohup sh sim_016.sh > sim_016.log & 
# nohup sh sim_017.sh > sim_017.log & 
# nohup sh sim_018.sh > sim_018.log & 
# nohup sh sim_019.sh > sim_019.log & 
# nohup sh sim_020.sh > sim_020.log & 
# nohup sh sim_021.sh > sim_021.log & 
# nohup sh sim_022.sh > sim_022.log & 
# nohup sh sim_023.sh > sim_023.log & 
# nohup sh sim_024.sh > sim_024.log & 
# nohup sh sim_025.sh > sim_025.log & 
# nohup sh sim_026.sh > sim_026.log & 
# nohup sh sim_027.sh > sim_027.log & 
# nohup sh sim_028.sh > sim_028.log & 
# nohup sh sim_029.sh > sim_029.log & 
# nohup sh sim_030.sh > sim_030.log & 
# nohup sh sim_031.sh > sim_031.log & 
# nohup sh sim_032.sh > sim_032.log & 
# nohup sh sim_033.sh > sim_033.log & 
# nohup sh sim_034.sh > sim_034.log & 
# nohup sh sim_035.sh > sim_035.log & 
# nohup sh sim_036.sh > sim_036.log & 
# nohup sh sim_037.sh > sim_037.log & 
# nohup sh sim_038.sh > sim_038.log & 
# nohup sh sim_039.sh > sim_039.log & 
# nohup sh sim_040.sh > sim_040.log & 
# nohup sh sim_041.sh > sim_041.log & 
# nohup sh sim_042.sh > sim_042.log & 
# nohup sh sim_043.sh > sim_043.log & 
# nohup sh sim_044.sh > sim_044.log & 
# nohup sh sim_045.sh > sim_045.log & 
# nohup sh sim_046.sh > sim_046.log & 
# nohup sh sim_047.sh > sim_047.log & 
# nohup sh sim_048.sh > sim_048.log & 
# nohup sh sim_049.sh > sim_049.log & 
# nohup sh sim_050.sh > sim_050.log & 
# nohup sh sim_051.sh > sim_051.log & 
# nohup sh sim_052.sh > sim_052.log & 
# nohup sh sim_053.sh > sim_053.log & 
# nohup sh sim_054.sh > sim_054.log & 
# nohup sh sim_055.sh > sim_055.log & 
# nohup sh sim_056.sh > sim_056.log & 
# nohup sh sim_057.sh > sim_057.log & 
# nohup sh sim_058.sh > sim_058.log & 
# nohup sh sim_059.sh > sim_059.log & 
# nohup sh sim_060.sh > sim_060.log & 
nohup sh sim_061.sh > sim_061.log & 
nohup sh sim_062.sh > sim_062.log & 
nohup sh sim_063.sh > sim_063.log & 
nohup sh sim_064.sh > sim_064.log & 
nohup sh sim_065.sh > sim_065.log & 
nohup sh sim_066.sh > sim_066.log & 
nohup sh sim_067.sh > sim_067.log & 
nohup sh sim_068.sh > sim_068.log & 
nohup sh sim_069.sh > sim_069.log & 
nohup sh sim_070.sh > sim_070.log & 
nohup sh sim_071.sh > sim_071.log & 
nohup sh sim_072.sh > sim_072.log & 
nohup sh sim_073.sh > sim_073.log & 
nohup sh sim_074.sh > sim_074.log & 
nohup sh sim_075.sh > sim_075.log & 
nohup sh sim_076.sh > sim_076.log & 
nohup sh sim_077.sh > sim_077.log & 
nohup sh sim_078.sh > sim_078.log & 
nohup sh sim_079.sh > sim_079.log & 
nohup sh sim_080.sh > sim_080.log & 
nohup sh sim_081.sh > sim_081.log & 
nohup sh sim_082.sh > sim_082.log & 
nohup sh sim_083.sh > sim_083.log & 
nohup sh sim_084.sh > sim_084.log & 
nohup sh sim_085.sh > sim_085.log & 
nohup sh sim_086.sh > sim_086.log & 
nohup sh sim_087.sh > sim_087.log & 
nohup sh sim_088.sh > sim_088.log & 
nohup sh sim_089.sh > sim_089.log & 
nohup sh sim_090.sh > sim_090.log & 
nohup sh sim_091.sh > sim_091.log & 
nohup sh sim_092.sh > sim_092.log & 
nohup sh sim_093.sh > sim_093.log & 
nohup sh sim_094.sh > sim_094.log & 
nohup sh sim_095.sh > sim_095.log & 
nohup sh sim_096.sh > sim_096.log & 
nohup sh sim_097.sh > sim_097.log & 
nohup sh sim_098.sh > sim_098.log & 
nohup sh sim_099.sh > sim_099.log & 
nohup sh sim_100.sh > sim_100.log & 
nohup sh sim_101.sh > sim_101.log & 
nohup sh sim_102.sh > sim_102.log & 
nohup sh sim_103.sh > sim_103.log & 
nohup sh sim_104.sh > sim_104.log & 
nohup sh sim_105.sh > sim_105.log & 
nohup sh sim_106.sh > sim_106.log & 
nohup sh sim_107.sh > sim_107.log & 
nohup sh sim_108.sh > sim_108.log & 
nohup sh sim_109.sh > sim_109.log & 
nohup sh sim_110.sh > sim_110.log & 
nohup sh sim_111.sh > sim_111.log & 
nohup sh sim_112.sh > sim_112.log & 
nohup sh sim_113.sh > sim_113.log & 
nohup sh sim_114.sh > sim_114.log & 
nohup sh sim_115.sh > sim_115.log & 
nohup sh sim_116.sh > sim_116.log & 
nohup sh sim_117.sh > sim_117.log & 
nohup sh sim_118.sh > sim_118.log & 
nohup sh sim_119.sh > sim_119.log & 
nohup sh sim_120.sh > sim_120.log & 
nohup sh sim_121.sh > sim_121.log & 
nohup sh sim_122.sh > sim_122.log & 
nohup sh sim_123.sh > sim_123.log & 
nohup sh sim_124.sh > sim_124.log & 
nohup sh sim_125.sh > sim_125.log & 
nohup sh sim_126.sh > sim_126.log & 
nohup sh sim_127.sh > sim_127.log & 
nohup sh sim_128.sh > sim_128.log & 
nohup sh sim_129.sh > sim_129.log & 
nohup sh sim_130.sh > sim_130.log & 
nohup sh sim_131.sh > sim_131.log & 
nohup sh sim_132.sh > sim_132.log & 
nohup sh sim_133.sh > sim_133.log & 
nohup sh sim_134.sh > sim_134.log & 
nohup sh sim_135.sh > sim_135.log & 
nohup sh sim_136.sh > sim_136.log & 
nohup sh sim_137.sh > sim_137.log & 
nohup sh sim_138.sh > sim_138.log & 
nohup sh sim_139.sh > sim_139.log & 
nohup sh sim_140.sh > sim_140.log & 
nohup sh sim_141.sh > sim_141.log & 
nohup sh sim_142.sh > sim_142.log & 
nohup sh sim_143.sh > sim_143.log & 
nohup sh sim_144.sh > sim_144.log & 
nohup sh sim_145.sh > sim_145.log & 
nohup sh sim_146.sh > sim_146.log & 
nohup sh sim_147.sh > sim_147.log & 
nohup sh sim_148.sh > sim_148.log & 
nohup sh sim_149.sh > sim_149.log & 
nohup sh sim_150.sh > sim_150.log & 
nohup sh sim_151.sh > sim_151.log & 
nohup sh sim_152.sh > sim_152.log & 
nohup sh sim_153.sh > sim_153.log & 
nohup sh sim_154.sh > sim_154.log & 
nohup sh sim_155.sh > sim_155.log & 
nohup sh sim_156.sh > sim_156.log & 
nohup sh sim_157.sh > sim_157.log & 
nohup sh sim_158.sh > sim_158.log & 
nohup sh sim_159.sh > sim_159.log & 
nohup sh sim_160.sh > sim_160.log & 
nohup sh sim_161.sh > sim_161.log & 
nohup sh sim_162.sh > sim_162.log & 
nohup sh sim_163.sh > sim_163.log & 
nohup sh sim_164.sh > sim_164.log & 
nohup sh sim_165.sh > sim_165.log & 
nohup sh sim_166.sh > sim_166.log & 
nohup sh sim_167.sh > sim_167.log & 
nohup sh sim_168.sh > sim_168.log & 
nohup sh sim_169.sh > sim_169.log & 
nohup sh sim_170.sh > sim_170.log & 
nohup sh sim_171.sh > sim_171.log & 
nohup sh sim_172.sh > sim_172.log & 
nohup sh sim_173.sh > sim_173.log & 
nohup sh sim_174.sh > sim_174.log & 
nohup sh sim_175.sh > sim_175.log & 
nohup sh sim_176.sh > sim_176.log & 
nohup sh sim_177.sh > sim_177.log & 
nohup sh sim_178.sh > sim_178.log & 
nohup sh sim_179.sh > sim_179.log & 
nohup sh sim_180.sh > sim_180.log & 
nohup sh sim_181.sh > sim_181.log & 
nohup sh sim_182.sh > sim_182.log & 
nohup sh sim_183.sh > sim_183.log & 
nohup sh sim_184.sh > sim_184.log & 
nohup sh sim_185.sh > sim_185.log & 
nohup sh sim_186.sh > sim_186.log & 
nohup sh sim_187.sh > sim_187.log & 
nohup sh sim_188.sh > sim_188.log & 
nohup sh sim_189.sh > sim_189.log & 
nohup sh sim_190.sh > sim_190.log & 
nohup sh sim_191.sh > sim_191.log & 
nohup sh sim_192.sh > sim_192.log & 
nohup sh sim_193.sh > sim_193.log & 
nohup sh sim_194.sh > sim_194.log & 
nohup sh sim_195.sh > sim_195.log & 
nohup sh sim_196.sh > sim_196.log & 
nohup sh sim_197.sh > sim_197.log & 
nohup sh sim_198.sh > sim_198.log & 
nohup sh sim_199.sh > sim_199.log & 
nohup sh sim_200.sh > sim_200.log & 
nohup sh sim_201.sh > sim_201.log & 
nohup sh sim_202.sh > sim_202.log & 
nohup sh sim_203.sh > sim_203.log & 
nohup sh sim_204.sh > sim_204.log & 
nohup sh sim_205.sh > sim_205.log & 
nohup sh sim_206.sh > sim_206.log & 
nohup sh sim_207.sh > sim_207.log & 
nohup sh sim_208.sh > sim_208.log & 
nohup sh sim_209.sh > sim_209.log & 
nohup sh sim_210.sh > sim_210.log & 
nohup sh sim_211.sh > sim_211.log & 
nohup sh sim_212.sh > sim_212.log & 
nohup sh sim_213.sh > sim_213.log & 
nohup sh sim_214.sh > sim_214.log & 
nohup sh sim_215.sh > sim_215.log & 
nohup sh sim_216.sh > sim_216.log & 
nohup sh sim_217.sh > sim_217.log & 
nohup sh sim_218.sh > sim_218.log & 
nohup sh sim_219.sh > sim_219.log & 
nohup sh sim_220.sh > sim_220.log & 
nohup sh sim_221.sh > sim_221.log & 
nohup sh sim_222.sh > sim_222.log & 
nohup sh sim_223.sh > sim_223.log & 
nohup sh sim_224.sh > sim_224.log & 
nohup sh sim_225.sh > sim_225.log & 
nohup sh sim_226.sh > sim_226.log & 
nohup sh sim_227.sh > sim_227.log & 
nohup sh sim_228.sh > sim_228.log & 
nohup sh sim_229.sh > sim_229.log & 
nohup sh sim_230.sh > sim_230.log & 
nohup sh sim_231.sh > sim_231.log & 
nohup sh sim_232.sh > sim_232.log & 
nohup sh sim_233.sh > sim_233.log & 
nohup sh sim_234.sh > sim_234.log & 
nohup sh sim_235.sh > sim_235.log & 
nohup sh sim_236.sh > sim_236.log & 
nohup sh sim_237.sh > sim_237.log & 
nohup sh sim_238.sh > sim_238.log & 
nohup sh sim_239.sh > sim_239.log & 
nohup sh sim_240.sh > sim_240.log & 
nohup sh sim_241.sh > sim_241.log & 
nohup sh sim_242.sh > sim_242.log & 
nohup sh sim_243.sh > sim_243.log & 
nohup sh sim_244.sh > sim_244.log & 
nohup sh sim_245.sh > sim_245.log & 
nohup sh sim_246.sh > sim_246.log & 
nohup sh sim_247.sh > sim_247.log & 
nohup sh sim_248.sh > sim_248.log & 
nohup sh sim_249.sh > sim_249.log & 
nohup sh sim_250.sh > sim_250.log & 
nohup sh sim_251.sh > sim_251.log & 
nohup sh sim_252.sh > sim_252.log & 
nohup sh sim_253.sh > sim_253.log & 
nohup sh sim_254.sh > sim_254.log & 
nohup sh sim_255.sh > sim_255.log & 
nohup sh sim_256.sh > sim_256.log & 
nohup sh sim_257.sh > sim_257.log & 
nohup sh sim_258.sh > sim_258.log & 
nohup sh sim_259.sh > sim_259.log & 
nohup sh sim_260.sh > sim_260.log & 
nohup sh sim_261.sh > sim_261.log & 
nohup sh sim_262.sh > sim_262.log & 
nohup sh sim_263.sh > sim_263.log & 
nohup sh sim_264.sh > sim_264.log & 
nohup sh sim_265.sh > sim_265.log & 
nohup sh sim_266.sh > sim_266.log & 
nohup sh sim_267.sh > sim_267.log & 
nohup sh sim_268.sh > sim_268.log & 
nohup sh sim_269.sh > sim_269.log & 
nohup sh sim_270.sh > sim_270.log & 
nohup sh sim_271.sh > sim_271.log & 
nohup sh sim_272.sh > sim_272.log & 
nohup sh sim_273.sh > sim_273.log & 
nohup sh sim_274.sh > sim_274.log & 
nohup sh sim_275.sh > sim_275.log & 
nohup sh sim_276.sh > sim_276.log & 
nohup sh sim_277.sh > sim_277.log & 
nohup sh sim_278.sh > sim_278.log & 
nohup sh sim_279.sh > sim_279.log & 
nohup sh sim_280.sh > sim_280.log & 
nohup sh sim_281.sh > sim_281.log & 
nohup sh sim_282.sh > sim_282.log & 
nohup sh sim_283.sh > sim_283.log & 
nohup sh sim_284.sh > sim_284.log & 
nohup sh sim_285.sh > sim_285.log & 
nohup sh sim_286.sh > sim_286.log & 
nohup sh sim_287.sh > sim_287.log & 
nohup sh sim_288.sh > sim_288.log & 
nohup sh sim_289.sh > sim_289.log & 
nohup sh sim_290.sh > sim_290.log & 
nohup sh sim_291.sh > sim_291.log & 
nohup sh sim_292.sh > sim_292.log & 
nohup sh sim_293.sh > sim_293.log & 
nohup sh sim_294.sh > sim_294.log & 
nohup sh sim_295.sh > sim_295.log & 
nohup sh sim_296.sh > sim_296.log & 
nohup sh sim_297.sh > sim_297.log & 
nohup sh sim_298.sh > sim_298.log & 
nohup sh sim_299.sh > sim_299.log & 
nohup sh sim_300.sh > sim_300.log & 
nohup sh sim_301.sh > sim_301.log & 
nohup sh sim_302.sh > sim_302.log & 
nohup sh sim_303.sh > sim_303.log & 
nohup sh sim_304.sh > sim_304.log & 
nohup sh sim_305.sh > sim_305.log & 
nohup sh sim_306.sh > sim_306.log & 
nohup sh sim_307.sh > sim_307.log & 
nohup sh sim_308.sh > sim_308.log & 
nohup sh sim_309.sh > sim_309.log & 
nohup sh sim_310.sh > sim_310.log & 
nohup sh sim_311.sh > sim_311.log & 
nohup sh sim_312.sh > sim_312.log & 
nohup sh sim_313.sh > sim_313.log & 
nohup sh sim_314.sh > sim_314.log & 
nohup sh sim_315.sh > sim_315.log & 
nohup sh sim_316.sh > sim_316.log & 
nohup sh sim_317.sh > sim_317.log & 
nohup sh sim_318.sh > sim_318.log & 
nohup sh sim_319.sh > sim_319.log & 
nohup sh sim_320.sh > sim_320.log & 
nohup sh sim_321.sh > sim_321.log & 
nohup sh sim_322.sh > sim_322.log & 
nohup sh sim_323.sh > sim_323.log & 
nohup sh sim_324.sh > sim_324.log & 
nohup sh sim_325.sh > sim_325.log & 
nohup sh sim_326.sh > sim_326.log & 
nohup sh sim_327.sh > sim_327.log & 
nohup sh sim_328.sh > sim_328.log & 
nohup sh sim_329.sh > sim_329.log & 
nohup sh sim_330.sh > sim_330.log & 
nohup sh sim_331.sh > sim_331.log & 
nohup sh sim_332.sh > sim_332.log & 
nohup sh sim_333.sh > sim_333.log & 
nohup sh sim_334.sh > sim_334.log & 
nohup sh sim_335.sh > sim_335.log & 
nohup sh sim_336.sh > sim_336.log & 
nohup sh sim_337.sh > sim_337.log & 
nohup sh sim_338.sh > sim_338.log & 
nohup sh sim_339.sh > sim_339.log & 
nohup sh sim_340.sh > sim_340.log & 
nohup sh sim_341.sh > sim_341.log & 
nohup sh sim_342.sh > sim_342.log & 
nohup sh sim_343.sh > sim_343.log & 
nohup sh sim_344.sh > sim_344.log & 
nohup sh sim_345.sh > sim_345.log & 
nohup sh sim_346.sh > sim_346.log & 
nohup sh sim_347.sh > sim_347.log & 
nohup sh sim_348.sh > sim_348.log & 
nohup sh sim_349.sh > sim_349.log & 
# nohup sh sim_350.sh > sim_350.log & 
# nohup sh sim_351.sh > sim_351.log & 
# nohup sh sim_352.sh > sim_352.log & 
# nohup sh sim_353.sh > sim_353.log & 
# nohup sh sim_354.sh > sim_354.log & 
# nohup sh sim_355.sh > sim_355.log & 
# nohup sh sim_356.sh > sim_356.log & 
# nohup sh sim_357.sh > sim_357.log & 
# nohup sh sim_358.sh > sim_358.log & 
# nohup sh sim_359.sh > sim_359.log & 
nohup sh sim_360.sh > sim_360.log & 
nohup sh sim_361.sh > sim_361.log & 
nohup sh sim_362.sh > sim_362.log & 
nohup sh sim_363.sh > sim_363.log & 
nohup sh sim_364.sh > sim_364.log & 
nohup sh sim_365.sh > sim_365.log & 
nohup sh sim_366.sh > sim_366.log & 
nohup sh sim_367.sh > sim_367.log & 
nohup sh sim_368.sh > sim_368.log & 
nohup sh sim_369.sh > sim_369.log & 
nohup sh sim_370.sh > sim_370.log & 
nohup sh sim_371.sh > sim_371.log & 
nohup sh sim_372.sh > sim_372.log & 
nohup sh sim_373.sh > sim_373.log & 
nohup sh sim_374.sh > sim_374.log & 
nohup sh sim_375.sh > sim_375.log & 
nohup sh sim_376.sh > sim_376.log & 
nohup sh sim_377.sh > sim_377.log & 
nohup sh sim_378.sh > sim_378.log & 
nohup sh sim_379.sh > sim_379.log & 
nohup sh sim_380.sh > sim_380.log & 
nohup sh sim_381.sh > sim_381.log & 
nohup sh sim_382.sh > sim_382.log & 
nohup sh sim_383.sh > sim_383.log & 
nohup sh sim_384.sh > sim_384.log & 
nohup sh sim_385.sh > sim_385.log & 
nohup sh sim_386.sh > sim_386.log & 
nohup sh sim_387.sh > sim_387.log & 
nohup sh sim_388.sh > sim_388.log & 
nohup sh sim_389.sh > sim_389.log & 
nohup sh sim_390.sh > sim_390.log & 
nohup sh sim_391.sh > sim_391.log & 
nohup sh sim_392.sh > sim_392.log & 
nohup sh sim_393.sh > sim_393.log & 
nohup sh sim_394.sh > sim_394.log & 
nohup sh sim_395.sh > sim_395.log & 
nohup sh sim_396.sh > sim_396.log & 
nohup sh sim_397.sh > sim_397.log & 
nohup sh sim_398.sh > sim_398.log & 
nohup sh sim_399.sh > sim_399.log & 
nohup sh sim_400.sh > sim_400.log & 
nohup sh sim_401.sh > sim_401.log & 
nohup sh sim_402.sh > sim_402.log & 
nohup sh sim_403.sh > sim_403.log & 
nohup sh sim_404.sh > sim_404.log & 
nohup sh sim_405.sh > sim_405.log & 
nohup sh sim_406.sh > sim_406.log & 
nohup sh sim_407.sh > sim_407.log & 
nohup sh sim_408.sh > sim_408.log & 
nohup sh sim_409.sh > sim_409.log & 
nohup sh sim_410.sh > sim_410.log & 
nohup sh sim_411.sh > sim_411.log & 
nohup sh sim_412.sh > sim_412.log & 
nohup sh sim_413.sh > sim_413.log & 
nohup sh sim_414.sh > sim_414.log & 
nohup sh sim_415.sh > sim_415.log & 
nohup sh sim_416.sh > sim_416.log & 
nohup sh sim_417.sh > sim_417.log & 
nohup sh sim_418.sh > sim_418.log & 
nohup sh sim_419.sh > sim_419.log & 
nohup sh sim_420.sh > sim_420.log & 
nohup sh sim_421.sh > sim_421.log & 
nohup sh sim_422.sh > sim_422.log & 
nohup sh sim_423.sh > sim_423.log & 
nohup sh sim_424.sh > sim_424.log & 
nohup sh sim_425.sh > sim_425.log & 
nohup sh sim_426.sh > sim_426.log & 
nohup sh sim_427.sh > sim_427.log & 
nohup sh sim_428.sh > sim_428.log & 
nohup sh sim_429.sh > sim_429.log & 
nohup sh sim_430.sh > sim_430.log & 
nohup sh sim_431.sh > sim_431.log & 
nohup sh sim_432.sh > sim_432.log & 
nohup sh sim_433.sh > sim_433.log & 
nohup sh sim_434.sh > sim_434.log & 
nohup sh sim_435.sh > sim_435.log & 
nohup sh sim_436.sh > sim_436.log & 
nohup sh sim_437.sh > sim_437.log & 
nohup sh sim_438.sh > sim_438.log & 
nohup sh sim_439.sh > sim_439.log & 
nohup sh sim_440.sh > sim_440.log & 
nohup sh sim_441.sh > sim_441.log & 
nohup sh sim_442.sh > sim_442.log & 
nohup sh sim_443.sh > sim_443.log & 
nohup sh sim_444.sh > sim_444.log & 
nohup sh sim_445.sh > sim_445.log & 
nohup sh sim_446.sh > sim_446.log & 
nohup sh sim_447.sh > sim_447.log & 
nohup sh sim_448.sh > sim_448.log & 
nohup sh sim_449.sh > sim_449.log & 
nohup sh sim_450.sh > sim_450.log & 
nohup sh sim_451.sh > sim_451.log & 
nohup sh sim_452.sh > sim_452.log & 
nohup sh sim_453.sh > sim_453.log & 
nohup sh sim_454.sh > sim_454.log & 
nohup sh sim_455.sh > sim_455.log & 
nohup sh sim_456.sh > sim_456.log & 
nohup sh sim_457.sh > sim_457.log & 
nohup sh sim_458.sh > sim_458.log & 
nohup sh sim_459.sh > sim_459.log & 
nohup sh sim_460.sh > sim_460.log & 
nohup sh sim_461.sh > sim_461.log & 
nohup sh sim_462.sh > sim_462.log & 
nohup sh sim_463.sh > sim_463.log & 
nohup sh sim_464.sh > sim_464.log & 
nohup sh sim_465.sh > sim_465.log & 
nohup sh sim_466.sh > sim_466.log & 
nohup sh sim_467.sh > sim_467.log & 
nohup sh sim_468.sh > sim_468.log & 
nohup sh sim_469.sh > sim_469.log & 
nohup sh sim_470.sh > sim_470.log & 
nohup sh sim_471.sh > sim_471.log & 
nohup sh sim_472.sh > sim_472.log & 
nohup sh sim_473.sh > sim_473.log & 
nohup sh sim_474.sh > sim_474.log & 
nohup sh sim_475.sh > sim_475.log & 
nohup sh sim_476.sh > sim_476.log & 
nohup sh sim_477.sh > sim_477.log & 
nohup sh sim_478.sh > sim_478.log & 
nohup sh sim_479.sh > sim_479.log & 
nohup sh sim_480.sh > sim_480.log & 
nohup sh sim_481.sh > sim_481.log & 
nohup sh sim_482.sh > sim_482.log & 
nohup sh sim_483.sh > sim_483.log & 
nohup sh sim_484.sh > sim_484.log & 
nohup sh sim_485.sh > sim_485.log & 
nohup sh sim_486.sh > sim_486.log & 
nohup sh sim_487.sh > sim_487.log & 
nohup sh sim_488.sh > sim_488.log & 
nohup sh sim_489.sh > sim_489.log & 
nohup sh sim_490.sh > sim_490.log & 
nohup sh sim_491.sh > sim_491.log & 
nohup sh sim_492.sh > sim_492.log & 
nohup sh sim_493.sh > sim_493.log & 
nohup sh sim_494.sh > sim_494.log & 
nohup sh sim_495.sh > sim_495.log & 
nohup sh sim_496.sh > sim_496.log & 
nohup sh sim_497.sh > sim_497.log & 
nohup sh sim_498.sh > sim_498.log & 
nohup sh sim_499.sh > sim_499.log & 
nohup sh sim_500.sh > sim_500.log & 
nohup sh sim_501.sh > sim_501.log & 
nohup sh sim_502.sh > sim_502.log & 
nohup sh sim_503.sh > sim_503.log & 
nohup sh sim_504.sh > sim_504.log & 
nohup sh sim_505.sh > sim_505.log & 
nohup sh sim_506.sh > sim_506.log & 
nohup sh sim_507.sh > sim_507.log & 
nohup sh sim_508.sh > sim_508.log & 
nohup sh sim_509.sh > sim_509.log & 
nohup sh sim_510.sh > sim_510.log & 
nohup sh sim_511.sh > sim_511.log & 
nohup sh sim_512.sh > sim_512.log & 
nohup sh sim_513.sh > sim_513.log & 
nohup sh sim_514.sh > sim_514.log & 
nohup sh sim_515.sh > sim_515.log & 
nohup sh sim_516.sh > sim_516.log & 
nohup sh sim_517.sh > sim_517.log & 
nohup sh sim_518.sh > sim_518.log & 
nohup sh sim_519.sh > sim_519.log & 
nohup sh sim_520.sh > sim_520.log & 
nohup sh sim_521.sh > sim_521.log & 
nohup sh sim_522.sh > sim_522.log & 
nohup sh sim_523.sh > sim_523.log & 
nohup sh sim_524.sh > sim_524.log & 
nohup sh sim_525.sh > sim_525.log & 
nohup sh sim_526.sh > sim_526.log & 
nohup sh sim_527.sh > sim_527.log & 
nohup sh sim_528.sh > sim_528.log & 
nohup sh sim_529.sh > sim_529.log & 
nohup sh sim_530.sh > sim_530.log & 
nohup sh sim_531.sh > sim_531.log & 
nohup sh sim_532.sh > sim_532.log & 
nohup sh sim_533.sh > sim_533.log & 
nohup sh sim_534.sh > sim_534.log & 
nohup sh sim_535.sh > sim_535.log & 
nohup sh sim_536.sh > sim_536.log & 
nohup sh sim_537.sh > sim_537.log & 
nohup sh sim_538.sh > sim_538.log & 
nohup sh sim_539.sh > sim_539.log & 
nohup sh sim_540.sh > sim_540.log & 
nohup sh sim_541.sh > sim_541.log & 
nohup sh sim_542.sh > sim_542.log & 
nohup sh sim_543.sh > sim_543.log & 
nohup sh sim_544.sh > sim_544.log & 
nohup sh sim_545.sh > sim_545.log & 
nohup sh sim_546.sh > sim_546.log & 
nohup sh sim_547.sh > sim_547.log & 
nohup sh sim_548.sh > sim_548.log & 
nohup sh sim_549.sh > sim_549.log & 
nohup sh sim_550.sh > sim_550.log & 
nohup sh sim_551.sh > sim_551.log & 
nohup sh sim_552.sh > sim_552.log & 
nohup sh sim_553.sh > sim_553.log & 
nohup sh sim_554.sh > sim_554.log & 
nohup sh sim_555.sh > sim_555.log & 
nohup sh sim_556.sh > sim_556.log & 
nohup sh sim_557.sh > sim_557.log & 
nohup sh sim_558.sh > sim_558.log & 
nohup sh sim_559.sh > sim_559.log & 
nohup sh sim_560.sh > sim_560.log & 
nohup sh sim_561.sh > sim_561.log & 
nohup sh sim_562.sh > sim_562.log & 
nohup sh sim_563.sh > sim_563.log & 
nohup sh sim_564.sh > sim_564.log & 
nohup sh sim_565.sh > sim_565.log & 
nohup sh sim_566.sh > sim_566.log & 
nohup sh sim_567.sh > sim_567.log & 
nohup sh sim_568.sh > sim_568.log & 
nohup sh sim_569.sh > sim_569.log & 
nohup sh sim_570.sh > sim_570.log & 
nohup sh sim_571.sh > sim_571.log & 
nohup sh sim_572.sh > sim_572.log & 
nohup sh sim_573.sh > sim_573.log & 
nohup sh sim_574.sh > sim_574.log & 
nohup sh sim_575.sh > sim_575.log & 
nohup sh sim_576.sh > sim_576.log & 
nohup sh sim_577.sh > sim_577.log & 
nohup sh sim_578.sh > sim_578.log & 
nohup sh sim_579.sh > sim_579.log & 
nohup sh sim_580.sh > sim_580.log & 
nohup sh sim_581.sh > sim_581.log & 
nohup sh sim_582.sh > sim_582.log & 
nohup sh sim_583.sh > sim_583.log & 
nohup sh sim_584.sh > sim_584.log & 
nohup sh sim_585.sh > sim_585.log & 
nohup sh sim_586.sh > sim_586.log & 
nohup sh sim_587.sh > sim_587.log & 
nohup sh sim_588.sh > sim_588.log & 
nohup sh sim_589.sh > sim_589.log & 
nohup sh sim_590.sh > sim_590.log & 
nohup sh sim_591.sh > sim_591.log & 
nohup sh sim_592.sh > sim_592.log & 
nohup sh sim_593.sh > sim_593.log & 
nohup sh sim_594.sh > sim_594.log & 
nohup sh sim_595.sh > sim_595.log & 
nohup sh sim_596.sh > sim_596.log & 
nohup sh sim_597.sh > sim_597.log & 
nohup sh sim_598.sh > sim_598.log & 
nohup sh sim_599.sh > sim_599.log & 
nohup sh sim_600.sh > sim_600.log & 
nohup sh sim_601.sh > sim_601.log & 
nohup sh sim_602.sh > sim_602.log & 
nohup sh sim_603.sh > sim_603.log & 
nohup sh sim_604.sh > sim_604.log & 
nohup sh sim_605.sh > sim_605.log & 
nohup sh sim_606.sh > sim_606.log & 
nohup sh sim_607.sh > sim_607.log & 
nohup sh sim_608.sh > sim_608.log & 
nohup sh sim_609.sh > sim_609.log & 
nohup sh sim_610.sh > sim_610.log & 
nohup sh sim_611.sh > sim_611.log & 
nohup sh sim_612.sh > sim_612.log & 
nohup sh sim_613.sh > sim_613.log & 
nohup sh sim_614.sh > sim_614.log & 
nohup sh sim_615.sh > sim_615.log & 
nohup sh sim_616.sh > sim_616.log & 
nohup sh sim_617.sh > sim_617.log & 
nohup sh sim_618.sh > sim_618.log & 
nohup sh sim_619.sh > sim_619.log & 
nohup sh sim_620.sh > sim_620.log & 
nohup sh sim_621.sh > sim_621.log & 
nohup sh sim_622.sh > sim_622.log & 
nohup sh sim_623.sh > sim_623.log & 
nohup sh sim_624.sh > sim_624.log & 
nohup sh sim_625.sh > sim_625.log & 
nohup sh sim_626.sh > sim_626.log & 
nohup sh sim_627.sh > sim_627.log & 
nohup sh sim_628.sh > sim_628.log & 
nohup sh sim_629.sh > sim_629.log & 
nohup sh sim_630.sh > sim_630.log & 
nohup sh sim_631.sh > sim_631.log & 
nohup sh sim_632.sh > sim_632.log & 
nohup sh sim_633.sh > sim_633.log & 
nohup sh sim_634.sh > sim_634.log & 
nohup sh sim_635.sh > sim_635.log & 
nohup sh sim_636.sh > sim_636.log & 
nohup sh sim_637.sh > sim_637.log & 
nohup sh sim_638.sh > sim_638.log & 
nohup sh sim_639.sh > sim_639.log & 
nohup sh sim_640.sh > sim_640.log & 
nohup sh sim_641.sh > sim_641.log & 
nohup sh sim_642.sh > sim_642.log & 
nohup sh sim_643.sh > sim_643.log & 
nohup sh sim_644.sh > sim_644.log & 
nohup sh sim_645.sh > sim_645.log & 
nohup sh sim_646.sh > sim_646.log & 
nohup sh sim_647.sh > sim_647.log & 
nohup sh sim_648.sh > sim_648.log & 
nohup sh sim_649.sh > sim_649.log & 
nohup sh sim_650.sh > sim_650.log & 
nohup sh sim_651.sh > sim_651.log & 
nohup sh sim_652.sh > sim_652.log & 
nohup sh sim_653.sh > sim_653.log & 
nohup sh sim_654.sh > sim_654.log & 
nohup sh sim_655.sh > sim_655.log & 
nohup sh sim_656.sh > sim_656.log & 
nohup sh sim_657.sh > sim_657.log & 
nohup sh sim_658.sh > sim_658.log & 
nohup sh sim_659.sh > sim_659.log & 
nohup sh sim_660.sh > sim_660.log & 
nohup sh sim_661.sh > sim_661.log & 
nohup sh sim_662.sh > sim_662.log & 
nohup sh sim_663.sh > sim_663.log & 
nohup sh sim_664.sh > sim_664.log & 
nohup sh sim_665.sh > sim_665.log & 
nohup sh sim_666.sh > sim_666.log & 
nohup sh sim_667.sh > sim_667.log & 
nohup sh sim_668.sh > sim_668.log & 
nohup sh sim_669.sh > sim_669.log & 
nohup sh sim_670.sh > sim_670.log & 
nohup sh sim_671.sh > sim_671.log & 
nohup sh sim_672.sh > sim_672.log & 
nohup sh sim_673.sh > sim_673.log & 
nohup sh sim_674.sh > sim_674.log & 
nohup sh sim_675.sh > sim_675.log & 
nohup sh sim_676.sh > sim_676.log & 
nohup sh sim_677.sh > sim_677.log & 
nohup sh sim_678.sh > sim_678.log & 
nohup sh sim_679.sh > sim_679.log & 
nohup sh sim_680.sh > sim_680.log & 
nohup sh sim_681.sh > sim_681.log & 
nohup sh sim_682.sh > sim_682.log & 
nohup sh sim_683.sh > sim_683.log & 
nohup sh sim_684.sh > sim_684.log & 
nohup sh sim_685.sh > sim_685.log & 
nohup sh sim_686.sh > sim_686.log & 
nohup sh sim_687.sh > sim_687.log & 
nohup sh sim_688.sh > sim_688.log & 
nohup sh sim_689.sh > sim_689.log & 
nohup sh sim_690.sh > sim_690.log & 
nohup sh sim_691.sh > sim_691.log & 
nohup sh sim_692.sh > sim_692.log & 
nohup sh sim_693.sh > sim_693.log & 
nohup sh sim_694.sh > sim_694.log & 
nohup sh sim_695.sh > sim_695.log & 
nohup sh sim_696.sh > sim_696.log & 
nohup sh sim_697.sh > sim_697.log & 
nohup sh sim_698.sh > sim_698.log & 
nohup sh sim_699.sh > sim_699.log & 
nohup sh sim_700.sh > sim_700.log & 
nohup sh sim_701.sh > sim_701.log & 
nohup sh sim_702.sh > sim_702.log & 
nohup sh sim_703.sh > sim_703.log & 
nohup sh sim_704.sh > sim_704.log & 
nohup sh sim_705.sh > sim_705.log & 
nohup sh sim_706.sh > sim_706.log & 
nohup sh sim_707.sh > sim_707.log & 
nohup sh sim_708.sh > sim_708.log & 
nohup sh sim_709.sh > sim_709.log & 
nohup sh sim_710.sh > sim_710.log & 
nohup sh sim_711.sh > sim_711.log & 
nohup sh sim_712.sh > sim_712.log & 
nohup sh sim_713.sh > sim_713.log & 
nohup sh sim_714.sh > sim_714.log & 
nohup sh sim_715.sh > sim_715.log & 
nohup sh sim_716.sh > sim_716.log & 
nohup sh sim_717.sh > sim_717.log & 
nohup sh sim_718.sh > sim_718.log & 
nohup sh sim_719.sh > sim_719.log & 
nohup sh sim_720.sh > sim_720.log & 
nohup sh sim_721.sh > sim_721.log & 
nohup sh sim_722.sh > sim_722.log & 
nohup sh sim_723.sh > sim_723.log & 
nohup sh sim_724.sh > sim_724.log & 
nohup sh sim_725.sh > sim_725.log & 
nohup sh sim_726.sh > sim_726.log & 
nohup sh sim_727.sh > sim_727.log & 
nohup sh sim_728.sh > sim_728.log & 
nohup sh sim_729.sh > sim_729.log & 
nohup sh sim_730.sh > sim_730.log & 
nohup sh sim_731.sh > sim_731.log & 
nohup sh sim_732.sh > sim_732.log & 
nohup sh sim_733.sh > sim_733.log & 
nohup sh sim_734.sh > sim_734.log & 
nohup sh sim_735.sh > sim_735.log & 
nohup sh sim_736.sh > sim_736.log & 
nohup sh sim_737.sh > sim_737.log & 
nohup sh sim_738.sh > sim_738.log & 
nohup sh sim_739.sh > sim_739.log & 
nohup sh sim_740.sh > sim_740.log & 
nohup sh sim_741.sh > sim_741.log & 
nohup sh sim_742.sh > sim_742.log & 
nohup sh sim_743.sh > sim_743.log & 
nohup sh sim_744.sh > sim_744.log & 
nohup sh sim_745.sh > sim_745.log & 
nohup sh sim_746.sh > sim_746.log & 
nohup sh sim_747.sh > sim_747.log & 
nohup sh sim_748.sh > sim_748.log & 
nohup sh sim_749.sh > sim_749.log & 
nohup sh sim_750.sh > sim_750.log & 
nohup sh sim_751.sh > sim_751.log & 
nohup sh sim_752.sh > sim_752.log & 
nohup sh sim_753.sh > sim_753.log & 
nohup sh sim_754.sh > sim_754.log & 
nohup sh sim_755.sh > sim_755.log & 
nohup sh sim_756.sh > sim_756.log & 
nohup sh sim_757.sh > sim_757.log & 
nohup sh sim_758.sh > sim_758.log & 
nohup sh sim_759.sh > sim_759.log & 
nohup sh sim_760.sh > sim_760.log & 
nohup sh sim_761.sh > sim_761.log & 
nohup sh sim_762.sh > sim_762.log & 
nohup sh sim_763.sh > sim_763.log & 
nohup sh sim_764.sh > sim_764.log & 
nohup sh sim_765.sh > sim_765.log & 
nohup sh sim_766.sh > sim_766.log & 
nohup sh sim_767.sh > sim_767.log & 
