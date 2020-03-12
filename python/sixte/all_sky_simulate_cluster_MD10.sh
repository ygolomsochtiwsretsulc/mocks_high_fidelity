
ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd1_evt.fits > list_ccd1.lis
stilts tcat in=@list_ccd1.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd1.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd2_evt.fits > list_ccd2.lis
stilts tcat in=@list_ccd2.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd2.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd3_evt.fits > list_ccd3.lis
stilts tcat in=@list_ccd3.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd3.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd4_evt.fits > list_ccd4.lis
stilts tcat in=@list_ccd4.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd4.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd5_evt.fits > list_ccd5.lis
stilts tcat in=@list_ccd5.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd5.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd6_evt.fits > list_ccd6.lis
stilts tcat in=@list_ccd6.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd6.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/???/erass_ccd7_evt.fits > list_ccd7.lis
stilts tcat in=@list_ccd7.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd7.fits ofmt=fits

ls /data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccd?.fits > list_ccdA.lis
stilts tcat in=@list_ccdA.lis ifmt=fits icmd=explodeall icmd='keepcols "RA DEC SRC_ID_1 TIME SIGNAL"' icmd='select "SRC_ID_1>0"' omode=out out=/data40s/erosim/eRASS/eRASS8_cluster_MD10/simulated_photons_ccdA.fits ofmt=fits


pyCONDA
cd $GIT_AGN_MOCK/python/sixte
# source /home/erosita/sw/sass-setup.sh eSASSusers_190925
python simulate_cluster_only.py 000 MD10                        # 000
python simulate_cluster_only.py 001 MD10                        # 001
python simulate_cluster_only.py 002 MD10                        # 002
python simulate_cluster_only.py 003 MD10                        # 003
python simulate_cluster_only.py 004 MD10                        # 004
python simulate_cluster_only.py 005 MD10                        # 005
python simulate_cluster_only.py 006 MD10                        # 006
python simulate_cluster_only.py 007 MD10                        # 007
python simulate_cluster_only.py 008 MD10                        # 008
python simulate_cluster_only.py 009 MD10                        # 009
python simulate_cluster_only.py 010 MD10                        # 010
python simulate_cluster_only.py 011 MD10                        # 011
python simulate_cluster_only.py 012 MD10                        # 012
python simulate_cluster_only.py 013 MD10                        # 013
python simulate_cluster_only.py 014 MD10                        # 014
python simulate_cluster_only.py 015 MD10                        # 015
python simulate_cluster_only.py 016 MD10                        # 016
python simulate_cluster_only.py 017 MD10                        # 017
python simulate_cluster_only.py 018 MD10                        # 018
python simulate_cluster_only.py 019 MD10                        # 019
python simulate_cluster_only.py 020 MD10                        # 020
python simulate_cluster_only.py 021 MD10                        # 021
python simulate_cluster_only.py 022 MD10                        # 022
python simulate_cluster_only.py 023 MD10                        # 023
python simulate_cluster_only.py 024 MD10                        # 024
python simulate_cluster_only.py 025 MD10                        # 025
python simulate_cluster_only.py 026 MD10                        # 026
python simulate_cluster_only.py 027 MD10                        # 027
python simulate_cluster_only.py 028 MD10                        # 028
python simulate_cluster_only.py 029 MD10                        # 029
python simulate_cluster_only.py 030 MD10                        # 030
python simulate_cluster_only.py 031 MD10                        # 031
python simulate_cluster_only.py 032 MD10                        # 032
python simulate_cluster_only.py 033 MD10                        # 033
python simulate_cluster_only.py 034 MD10                        # 034
python simulate_cluster_only.py 035 MD10                        # 035
python simulate_cluster_only.py 036 MD10                        # 036
python simulate_cluster_only.py 037 MD10                        # 037
python simulate_cluster_only.py 038 MD10                        # 038
python simulate_cluster_only.py 039 MD10                        # 039
python simulate_cluster_only.py 040 MD10                        # 040
python simulate_cluster_only.py 041 MD10                        # 041
python simulate_cluster_only.py 042 MD10                        # 042
python simulate_cluster_only.py 043 MD10                        # 043
python simulate_cluster_only.py 044 MD10                        # 044
python simulate_cluster_only.py 045 MD10                        # 045
python simulate_cluster_only.py 046 MD10                        # 046
python simulate_cluster_only.py 047 MD10                        # 047
python simulate_cluster_only.py 048 MD10                        # 048
python simulate_cluster_only.py 049 MD10                        # 049
python simulate_cluster_only.py 050 MD10                        # 050
python simulate_cluster_only.py 051 MD10                        # 051
python simulate_cluster_only.py 052 MD10                        # 052
python simulate_cluster_only.py 053 MD10                        # 053
python simulate_cluster_only.py 054 MD10                        # 054
python simulate_cluster_only.py 055 MD10                        # 055
python simulate_cluster_only.py 056 MD10                        # 056
python simulate_cluster_only.py 057 MD10                        # 057
python simulate_cluster_only.py 058 MD10                        # 058
python simulate_cluster_only.py 059 MD10                        
python simulate_cluster_only.py 060 MD10                        
python simulate_cluster_only.py 061 MD10                        
python simulate_cluster_only.py 062 MD10                        
python simulate_cluster_only.py 063 MD10                        
python simulate_cluster_only.py 064 MD10                        
python simulate_cluster_only.py 065 MD10                        
python simulate_cluster_only.py 066 MD10                        
python simulate_cluster_only.py 067 MD10                        
python simulate_cluster_only.py 068 MD10                        
python simulate_cluster_only.py 069 MD10                        
python simulate_cluster_only.py 070 MD10                        
python simulate_cluster_only.py 071 MD10                        
python simulate_cluster_only.py 072 MD10                        
python simulate_cluster_only.py 073 MD10                        
python simulate_cluster_only.py 074 MD10                        
python simulate_cluster_only.py 075 MD10                        
python simulate_cluster_only.py 076 MD10                        
python simulate_cluster_only.py 077 MD10                        
python simulate_cluster_only.py 078 MD10                        
python simulate_cluster_only.py 079 MD10                        
python simulate_cluster_only.py 080 MD10                        
python simulate_cluster_only.py 081 MD10                        
python simulate_cluster_only.py 082 MD10                        
python simulate_cluster_only.py 083 MD10                        
python simulate_cluster_only.py 084 MD10                        
python simulate_cluster_only.py 085 MD10                        
python simulate_cluster_only.py 086 MD10                        
python simulate_cluster_only.py 087 MD10                        
python simulate_cluster_only.py 088 MD10                        
python simulate_cluster_only.py 089 MD10                        
python simulate_cluster_only.py 090 MD10                        
python simulate_cluster_only.py 091 MD10                        
python simulate_cluster_only.py 092 MD10                        
python simulate_cluster_only.py 093 MD10                        
python simulate_cluster_only.py 094 MD10                        
python simulate_cluster_only.py 095 MD10                        
python simulate_cluster_only.py 096 MD10                        
python simulate_cluster_only.py 097 MD10                        
python simulate_cluster_only.py 098 MD10                        
python simulate_cluster_only.py 099 MD10                        

# #!/bin/bash                                                   
pyCONDA                                                         
cd $GIT_AGN_MOCK/python/sixte                                   
source /home/erosita/sw/sass-setup.sh eSASSusers_190520         
python simulate_cluster_only.py 100 MD10                        # 100
python simulate_cluster_only.py 101 MD10                        # 101
python simulate_cluster_only.py 102 MD10                        # 102
python simulate_cluster_only.py 103 MD10                        # 103
python simulate_cluster_only.py 104 MD10                        # 104
python simulate_cluster_only.py 105 MD10                        # 105
python simulate_cluster_only.py 106 MD10                        # 106
python simulate_cluster_only.py 107 MD10                        # 107
python simulate_cluster_only.py 108 MD10                        # 108
python simulate_cluster_only.py 109 MD10                        # 109
python simulate_cluster_only.py 110 MD10                        
python simulate_cluster_only.py 111 MD10                        
python simulate_cluster_only.py 112 MD10                        
python simulate_cluster_only.py 113 MD10                        
python simulate_cluster_only.py 114 MD10                        
python simulate_cluster_only.py 115 MD10                        
python simulate_cluster_only.py 116 MD10                        
python simulate_cluster_only.py 117 MD10                        
python simulate_cluster_only.py 118 MD10                        
python simulate_cluster_only.py 119 MD10                        
python simulate_cluster_only.py 120 MD10                        
python simulate_cluster_only.py 121 MD10                        
python simulate_cluster_only.py 122 MD10                        
python simulate_cluster_only.py 123 MD10                        
python simulate_cluster_only.py 124 MD10                        
python simulate_cluster_only.py 125 MD10                        
python simulate_cluster_only.py 126 MD10                        
python simulate_cluster_only.py 127 MD10                        
python simulate_cluster_only.py 128 MD10                        
python simulate_cluster_only.py 129 MD10                        
python simulate_cluster_only.py 130 MD10                        
python simulate_cluster_only.py 131 MD10                        
python simulate_cluster_only.py 132 MD10                        
python simulate_cluster_only.py 133 MD10                        
python simulate_cluster_only.py 134 MD10                        
python simulate_cluster_only.py 135 MD10                        
python simulate_cluster_only.py 136 MD10                        
python simulate_cluster_only.py 137 MD10                        
python simulate_cluster_only.py 138 MD10                        
python simulate_cluster_only.py 139 MD10                        
python simulate_cluster_only.py 140 MD10                        
python simulate_cluster_only.py 141 MD10                        
python simulate_cluster_only.py 142 MD10                        
python simulate_cluster_only.py 143 MD10                        
python simulate_cluster_only.py 144 MD10                        
python simulate_cluster_only.py 145 MD10                        
python simulate_cluster_only.py 146 MD10                        
python simulate_cluster_only.py 147 MD10                        
python simulate_cluster_only.py 148 MD10                        
python simulate_cluster_only.py 149 MD10                        
python simulate_cluster_only.py 150 MD10                        
python simulate_cluster_only.py 151 MD10                        
python simulate_cluster_only.py 152 MD10                        
python simulate_cluster_only.py 153 MD10                        
python simulate_cluster_only.py 154 MD10                        
python simulate_cluster_only.py 155 MD10                        
python simulate_cluster_only.py 156 MD10                        
python simulate_cluster_only.py 157 MD10                        
python simulate_cluster_only.py 158 MD10                        
python simulate_cluster_only.py 159 MD10                        
python simulate_cluster_only.py 160 MD10                        
python simulate_cluster_only.py 161 MD10                        
python simulate_cluster_only.py 162 MD10                        
python simulate_cluster_only.py 163 MD10                        
python simulate_cluster_only.py 164 MD10                        
python simulate_cluster_only.py 165 MD10                        
python simulate_cluster_only.py 166 MD10                        
python simulate_cluster_only.py 167 MD10                        
python simulate_cluster_only.py 168 MD10                        
python simulate_cluster_only.py 169 MD10                        
python simulate_cluster_only.py 170 MD10                        
python simulate_cluster_only.py 171 MD10                        
python simulate_cluster_only.py 172 MD10                        
python simulate_cluster_only.py 173 MD10                        
python simulate_cluster_only.py 174 MD10                        
python simulate_cluster_only.py 175 MD10                        
python simulate_cluster_only.py 176 MD10                        
python simulate_cluster_only.py 177 MD10                        
python simulate_cluster_only.py 178 MD10                        
python simulate_cluster_only.py 179 MD10                        
python simulate_cluster_only.py 180 MD10                        
python simulate_cluster_only.py 181 MD10                        
python simulate_cluster_only.py 182 MD10                        
python simulate_cluster_only.py 183 MD10                        
python simulate_cluster_only.py 184 MD10                        
python simulate_cluster_only.py 185 MD10                        
python simulate_cluster_only.py 186 MD10                        
python simulate_cluster_only.py 187 MD10                        
python simulate_cluster_only.py 188 MD10                        
python simulate_cluster_only.py 189 MD10                        
python simulate_cluster_only.py 190 MD10                        
python simulate_cluster_only.py 191 MD10                        
python simulate_cluster_only.py 192 MD10                        
python simulate_cluster_only.py 193 MD10                        
python simulate_cluster_only.py 194 MD10                        
python simulate_cluster_only.py 195 MD10                        
python simulate_cluster_only.py 196 MD10                        
python simulate_cluster_only.py 197 MD10                        
python simulate_cluster_only.py 198 MD10                        
python simulate_cluster_only.py 199 MD10                        

#!/bin/bash                                                     
pyCONDA                                                         
cd $GIT_AGN_MOCK/python/sixte                                   
# source /home/erosita/sw/sass-setup.sh eSASSusers_190520         
python simulate_cluster_only.py 200 MD10                        # 200
python simulate_cluster_only.py 201 MD10                        # 201
python simulate_cluster_only.py 202 MD10                        # 202
python simulate_cluster_only.py 203 MD10                        # 203
python simulate_cluster_only.py 204 MD10                        # 204
python simulate_cluster_only.py 205 MD10                        # 205
python simulate_cluster_only.py 206 MD10                        # 206
python simulate_cluster_only.py 207 MD10                        # 207
python simulate_cluster_only.py 208 MD10                        # 208
python simulate_cluster_only.py 209 MD10                        # 209
python simulate_cluster_only.py 210 MD10                        # 210
python simulate_cluster_only.py 211 MD10                        # 211
python simulate_cluster_only.py 212 MD10                        # 212
python simulate_cluster_only.py 213 MD10                        # 213
python simulate_cluster_only.py 214 MD10                        # 214
python simulate_cluster_only.py 215 MD10                        # 215
python simulate_cluster_only.py 216 MD10                        # 216
python simulate_cluster_only.py 217 MD10                        # 217
python simulate_cluster_only.py 218 MD10                        # 218
python simulate_cluster_only.py 219 MD10                        # 219
python simulate_cluster_only.py 220 MD10                        # 220
python simulate_cluster_only.py 221 MD10                        # 221
python simulate_cluster_only.py 222 MD10                        # 222
python simulate_cluster_only.py 223 MD10                        # 223
python simulate_cluster_only.py 224 MD10                        # 224
python simulate_cluster_only.py 225 MD10                        # 225
python simulate_cluster_only.py 226 MD10                        # 226
python simulate_cluster_only.py 227 MD10                        # 227
python simulate_cluster_only.py 228 MD10                        # 228
python simulate_cluster_only.py 229 MD10                        # 229
python simulate_cluster_only.py 230 MD10                        
python simulate_cluster_only.py 231 MD10                        
python simulate_cluster_only.py 232 MD10                        
python simulate_cluster_only.py 233 MD10                        
python simulate_cluster_only.py 234 MD10                        
python simulate_cluster_only.py 235 MD10                        
python simulate_cluster_only.py 236 MD10                        
python simulate_cluster_only.py 237 MD10                        
python simulate_cluster_only.py 238 MD10                        
python simulate_cluster_only.py 239 MD10                        
python simulate_cluster_only.py 240 MD10                        
python simulate_cluster_only.py 241 MD10                        
python simulate_cluster_only.py 242 MD10                        
python simulate_cluster_only.py 243 MD10                        
python simulate_cluster_only.py 244 MD10                        
python simulate_cluster_only.py 245 MD10                        
python simulate_cluster_only.py 246 MD10                        
python simulate_cluster_only.py 247 MD10                        
python simulate_cluster_only.py 248 MD10                        
python simulate_cluster_only.py 249 MD10                        
python simulate_cluster_only.py 250 MD10                        
python simulate_cluster_only.py 251 MD10                        
python simulate_cluster_only.py 252 MD10                        
python simulate_cluster_only.py 253 MD10                        
python simulate_cluster_only.py 254 MD10                        
python simulate_cluster_only.py 255 MD10                        
python simulate_cluster_only.py 256 MD10                        
python simulate_cluster_only.py 257 MD10                        
python simulate_cluster_only.py 258 MD10                        
python simulate_cluster_only.py 259 MD10                        
python simulate_cluster_only.py 260 MD10                        
python simulate_cluster_only.py 261 MD10                        
python simulate_cluster_only.py 262 MD10                        
python simulate_cluster_only.py 263 MD10                        
python simulate_cluster_only.py 264 MD10                        
python simulate_cluster_only.py 265 MD10                        
python simulate_cluster_only.py 266 MD10                        
python simulate_cluster_only.py 267 MD10                        
python simulate_cluster_only.py 268 MD10                        
python simulate_cluster_only.py 269 MD10                        
python simulate_cluster_only.py 270 MD10                        
python simulate_cluster_only.py 271 MD10                        
python simulate_cluster_only.py 272 MD10                        
python simulate_cluster_only.py 273 MD10                        
python simulate_cluster_only.py 274 MD10                        
python simulate_cluster_only.py 275 MD10                        
python simulate_cluster_only.py 276 MD10                        
python simulate_cluster_only.py 277 MD10                        
python simulate_cluster_only.py 278 MD10                        
python simulate_cluster_only.py 279 MD10                        
python simulate_cluster_only.py 280 MD10                        
python simulate_cluster_only.py 281 MD10                        
python simulate_cluster_only.py 282 MD10                        
python simulate_cluster_only.py 283 MD10                        
python simulate_cluster_only.py 284 MD10                        
python simulate_cluster_only.py 285 MD10                        
python simulate_cluster_only.py 286 MD10                        
python simulate_cluster_only.py 287 MD10                        
python simulate_cluster_only.py 288 MD10                        
python simulate_cluster_only.py 289 MD10                        
python simulate_cluster_only.py 290 MD10                        
python simulate_cluster_only.py 291 MD10                        
python simulate_cluster_only.py 292 MD10                        
python simulate_cluster_only.py 293 MD10                        
python simulate_cluster_only.py 294 MD10                        
python simulate_cluster_only.py 295 MD10                        
python simulate_cluster_only.py 296 MD10                        
python simulate_cluster_only.py 297 MD10                        
python simulate_cluster_only.py 298 MD10                        
python simulate_cluster_only.py 299 MD10                        
                                                                 
#!/bin/bash                                                      
pyCONDA                                                          
cd $GIT_AGN_MOCK/python/sixte                                    
# source /home/erosita/sw/sass-setup.sh eSASSusers_190520          
python simulate_cluster_only.py 300 MD10                        # 300 
python simulate_cluster_only.py 301 MD10                        # 301 
python simulate_cluster_only.py 302 MD10                        # 302 
python simulate_cluster_only.py 303 MD10                        # 303 
python simulate_cluster_only.py 304 MD10                        # 304 
python simulate_cluster_only.py 305 MD10                        # 305 
python simulate_cluster_only.py 306 MD10                        # 306 
python simulate_cluster_only.py 307 MD10                        # 307 
python simulate_cluster_only.py 308 MD10                        # 308 
python simulate_cluster_only.py 309 MD10                        # 309 
python simulate_cluster_only.py 310 MD10                        # 310 
python simulate_cluster_only.py 311 MD10                        # 311 
python simulate_cluster_only.py 312 MD10                        # 312 
python simulate_cluster_only.py 313 MD10                        # 313 
python simulate_cluster_only.py 314 MD10                        # 314 
python simulate_cluster_only.py 315 MD10                        # 315 
python simulate_cluster_only.py 316 MD10                        # 316 
python simulate_cluster_only.py 317 MD10                        # 317 
python simulate_cluster_only.py 318 MD10                        # 318 
python simulate_cluster_only.py 319 MD10                        # 319 
python simulate_cluster_only.py 320 MD10                        # 320 
python simulate_cluster_only.py 321 MD10                        # 321#  
python simulate_cluster_only.py 322 MD10                        # 322#  
python simulate_cluster_only.py 323 MD10                        # 323#  
python simulate_cluster_only.py 324 MD10                        # 324#  
python simulate_cluster_only.py 325 MD10                        # 325#  
python simulate_cluster_only.py 326 MD10                        # 326#  
python simulate_cluster_only.py 327 MD10                        # 327#  
python simulate_cluster_only.py 328 MD10                        # 328#  
python simulate_cluster_only.py 329 MD10                        # 329#  
python simulate_cluster_only.py 330 MD10                        
python simulate_cluster_only.py 331 MD10                        
python simulate_cluster_only.py 332 MD10                        
python simulate_cluster_only.py 333 MD10                        
python simulate_cluster_only.py 334 MD10                        
python simulate_cluster_only.py 335 MD10                        
python simulate_cluster_only.py 336 MD10                        
python simulate_cluster_only.py 337 MD10                        
python simulate_cluster_only.py 338 MD10                        
python simulate_cluster_only.py 339 MD10                        
python simulate_cluster_only.py 340 MD10                        
python simulate_cluster_only.py 341 MD10                        
python simulate_cluster_only.py 342 MD10                        
python simulate_cluster_only.py 343 MD10                        
python simulate_cluster_only.py 344 MD10                        
python simulate_cluster_only.py 345 MD10                        
python simulate_cluster_only.py 346 MD10                        
python simulate_cluster_only.py 347 MD10                        
python simulate_cluster_only.py 348 MD10                        
python simulate_cluster_only.py 349 MD10                        
python simulate_cluster_only.py 350 MD10                        
python simulate_cluster_only.py 351 MD10                        

pyCONDA                                                          
cd $GIT_AGN_MOCK/python/sixte                                    

python simulate_cluster_only.py 352 MD10                        
python simulate_cluster_only.py 353 MD10                        
python simulate_cluster_only.py 354 MD10                        
python simulate_cluster_only.py 355 MD10                        
python simulate_cluster_only.py 356 MD10                        
python simulate_cluster_only.py 357 MD10                        
python simulate_cluster_only.py 358 MD10                        
python simulate_cluster_only.py 359 MD10                        
python simulate_cluster_only.py 360 MD10                        
python simulate_cluster_only.py 361 MD10                        
python simulate_cluster_only.py 362 MD10                        
python simulate_cluster_only.py 363 MD10                        
python simulate_cluster_only.py 364 MD10                        
python simulate_cluster_only.py 365 MD10                        
python simulate_cluster_only.py 366 MD10                        
python simulate_cluster_only.py 367 MD10                        
python simulate_cluster_only.py 368 MD10                        
python simulate_cluster_only.py 369 MD10                        
python simulate_cluster_only.py 370 MD10                        
python simulate_cluster_only.py 371 MD10                        
python simulate_cluster_only.py 372 MD10                        
python simulate_cluster_only.py 373 MD10                        
python simulate_cluster_only.py 374 MD10                        
python simulate_cluster_only.py 375 MD10                        
python simulate_cluster_only.py 376 MD10                        
python simulate_cluster_only.py 377 MD10                        
python simulate_cluster_only.py 378 MD10                        
python simulate_cluster_only.py 379 MD10                        
python simulate_cluster_only.py 380 MD10                        
python simulate_cluster_only.py 381 MD10                        
python simulate_cluster_only.py 382 MD10                        
python simulate_cluster_only.py 383 MD10                        
python simulate_cluster_only.py 384 MD10                        
python simulate_cluster_only.py 385 MD10                        
python simulate_cluster_only.py 386 MD10                        
python simulate_cluster_only.py 387 MD10                        
python simulate_cluster_only.py 388 MD10                        
python simulate_cluster_only.py 389 MD10                        
python simulate_cluster_only.py 390 MD10                        
python simulate_cluster_only.py 391 MD10                        
python simulate_cluster_only.py 392 MD10                        
python simulate_cluster_only.py 393 MD10                        
python simulate_cluster_only.py 394 MD10                        
python simulate_cluster_only.py 395 MD10                        
python simulate_cluster_only.py 396 MD10                        
python simulate_cluster_only.py 397 MD10                        
python simulate_cluster_only.py 398 MD10                        
python simulate_cluster_only.py 399 MD10                        
                                                                
#!/bin/bash                                                     
pyCONDA                                                         
cd $GIT_AGN_MOCK/python/sixte                                   
# source /home/erosita/sw/sass-setup.sh eSASSusers_190925         
python simulate_cluster_only.py 400 MD10                        # 400#  
python simulate_cluster_only.py 401 MD10                        # 401#  
python simulate_cluster_only.py 402 MD10                        # 402#  
python simulate_cluster_only.py 403 MD10                        # 403#  
python simulate_cluster_only.py 404 MD10                        # 404#  
python simulate_cluster_only.py 405 MD10                        # 405#  
python simulate_cluster_only.py 406 MD10                        # 406#  
python simulate_cluster_only.py 407 MD10                        # 407#  
python simulate_cluster_only.py 408 MD10                        # 408#  
python simulate_cluster_only.py 409 MD10                        # 409#  
python simulate_cluster_only.py 410 MD10                        # 410#  # 
python simulate_cluster_only.py 411 MD10                        # 411# # # 
python simulate_cluster_only.py 412 MD10                        # 412# # # 
python simulate_cluster_only.py 413 MD10                        # 413# # # 
python simulate_cluster_only.py 414 MD10                        # 414# # # 
python simulate_cluster_only.py 415 MD10                        # 415# # # 
python simulate_cluster_only.py 416 MD10                        # 416# # # 
python simulate_cluster_only.py 417 MD10                        # 417# # # 
python simulate_cluster_only.py 418 MD10                        # 418# # # 
python simulate_cluster_only.py 419 MD10                        # 419# # # 
python simulate_cluster_only.py 420 MD10                        # 420# # # 
python simulate_cluster_only.py 421 MD10                        # 421# # # 
python simulate_cluster_only.py 422 MD10                        # 422# # # 
python simulate_cluster_only.py 423 MD10                        # 423# # # 
python simulate_cluster_only.py 424 MD10                        # 424# # # 
python simulate_cluster_only.py 425 MD10                        # 425# # # 
python simulate_cluster_only.py 426 MD10                        # 426# # # 
python simulate_cluster_only.py 427 MD10                        # 427# # # 
python simulate_cluster_only.py 428 MD10                        # 428# # # 
python simulate_cluster_only.py 429 MD10                        # 429# # # 
python simulate_cluster_only.py 430 MD10                        # 430# # # 
python simulate_cluster_only.py 431 MD10                        # 431# # # 
python simulate_cluster_only.py 432 MD10                        # 432# # # 
python simulate_cluster_only.py 433 MD10                        # 433# # # 
python simulate_cluster_only.py 434 MD10                        # 434# # # 
python simulate_cluster_only.py 435 MD10                        # 435# # # 
python simulate_cluster_only.py 436 MD10                        # 436# # # 
python simulate_cluster_only.py 437 MD10                        # 437# # # 
python simulate_cluster_only.py 438 MD10                        # 438# # # 
python simulate_cluster_only.py 439 MD10                        # 439# # # 
python simulate_cluster_only.py 440 MD10                        # 440# ## # 
python simulate_cluster_only.py 441 MD10                        # 441# ## # 
python simulate_cluster_only.py 442 MD10                        # 442# ## # 
python simulate_cluster_only.py 443 MD10                        # 443# ## # 
python simulate_cluster_only.py 444 MD10                        # 444# ## # 
python simulate_cluster_only.py 445 MD10                        # 445# ## # 
python simulate_cluster_only.py 446 MD10                        # 446# ## # 
python simulate_cluster_only.py 447 MD10                        # 447# ## # 
python simulate_cluster_only.py 448 MD10                        # 448# ## # 
python simulate_cluster_only.py 449 MD10                        
python simulate_cluster_only.py 450 MD10                        
                                                                
pyCONDA                                                         
cd $GIT_AGN_MOCK/python/sixte                                   
python simulate_cluster_only.py 451 MD10                        
python simulate_cluster_only.py 452 MD10                        
python simulate_cluster_only.py 453 MD10                        
python simulate_cluster_only.py 454 MD10                        
python simulate_cluster_only.py 455 MD10                        
python simulate_cluster_only.py 456 MD10                        
python simulate_cluster_only.py 457 MD10                        
python simulate_cluster_only.py 458 MD10                        
python simulate_cluster_only.py 459 MD10                        
python simulate_cluster_only.py 460 MD10                        
python simulate_cluster_only.py 461 MD10                        
python simulate_cluster_only.py 462 MD10                        
python simulate_cluster_only.py 463 MD10                        
python simulate_cluster_only.py 464 MD10                        
python simulate_cluster_only.py 465 MD10                        
python simulate_cluster_only.py 466 MD10                        
python simulate_cluster_only.py 467 MD10                        
python simulate_cluster_only.py 468 MD10                        
python simulate_cluster_only.py 469 MD10                        
python simulate_cluster_only.py 470 MD10                        
python simulate_cluster_only.py 471 MD10                        
python simulate_cluster_only.py 472 MD10                        
python simulate_cluster_only.py 473 MD10                        
python simulate_cluster_only.py 474 MD10                        
python simulate_cluster_only.py 475 MD10                        
python simulate_cluster_only.py 476 MD10                        
python simulate_cluster_only.py 477 MD10                        
python simulate_cluster_only.py 478 MD10                        
python simulate_cluster_only.py 479 MD10                        
python simulate_cluster_only.py 480 MD10                        
python simulate_cluster_only.py 481 MD10                        
python simulate_cluster_only.py 482 MD10                        
python simulate_cluster_only.py 483 MD10                        
python simulate_cluster_only.py 484 MD10                        
python simulate_cluster_only.py 485 MD10                        # 
python simulate_cluster_only.py 486 MD10                        # 
python simulate_cluster_only.py 487 MD10                        # 
python simulate_cluster_only.py 488 MD10                        # 
python simulate_cluster_only.py 489 MD10                        # 
python simulate_cluster_only.py 490 MD10                        # 
python simulate_cluster_only.py 491 MD10                        # 
python simulate_cluster_only.py 492 MD10                        # 
python simulate_cluster_only.py 493 MD10                        # 
python simulate_cluster_only.py 494 MD10                        # 
python simulate_cluster_only.py 495 MD10                        # 
python simulate_cluster_only.py 496 MD10                        # 
python simulate_cluster_only.py 497 MD10                        # 
python simulate_cluster_only.py 498 MD10                        # 
python simulate_cluster_only.py 499 MD10                        # 
                                                                # 
#!/bin/bash                                                     # 
pyCONDA                                                         # 
cd $GIT_AGN_MOCK/python/sixte                                   # 
# source /home/erosita/sw/sass-setup.sh eSASSusers_190520         # 
python simulate_cluster_only.py 500 MD10                        # 500# ## # # 
python simulate_cluster_only.py 501 MD10                        # 501# ## # # 
python simulate_cluster_only.py 502 MD10                        # 502# ## # # 
python simulate_cluster_only.py 503 MD10                        # 503# ## # # 
python simulate_cluster_only.py 504 MD10                        # 504# ## # # 
python simulate_cluster_only.py 505 MD10                        # 505# ## # # 
python simulate_cluster_only.py 506 MD10                        # 506# ## # # 
python simulate_cluster_only.py 507 MD10                        # 507# ## # # 
python simulate_cluster_only.py 508 MD10                        # 508# ## # # 
python simulate_cluster_only.py 509 MD10                        # 509# ## # # 
python simulate_cluster_only.py 510 MD10                        # 510# ## # # 
python simulate_cluster_only.py 511 MD10                        # 511# ## # # 
python simulate_cluster_only.py 512 MD10                        # 512# ## # # 
python simulate_cluster_only.py 513 MD10                       
python simulate_cluster_only.py 514 MD10                       
python simulate_cluster_only.py 515 MD10                       
python simulate_cluster_only.py 516 MD10                       
python simulate_cluster_only.py 517 MD10                       
python simulate_cluster_only.py 518 MD10                       
python simulate_cluster_only.py 519 MD10                       
python simulate_cluster_only.py 520 MD10                       
python simulate_cluster_only.py 521 MD10                       
python simulate_cluster_only.py 522 MD10                       
python simulate_cluster_only.py 523 MD10                       
python simulate_cluster_only.py 524 MD10                       
python simulate_cluster_only.py 525 MD10                       
python simulate_cluster_only.py 526 MD10                       
python simulate_cluster_only.py 527 MD10                       
python simulate_cluster_only.py 528 MD10                       
python simulate_cluster_only.py 529 MD10                       
python simulate_cluster_only.py 530 MD10                       
python simulate_cluster_only.py 531 MD10                       
python simulate_cluster_only.py 532 MD10                       
python simulate_cluster_only.py 533 MD10                       
python simulate_cluster_only.py 534 MD10                       
python simulate_cluster_only.py 535 MD10                       
python simulate_cluster_only.py 536 MD10                       
python simulate_cluster_only.py 537 MD10                        # 
python simulate_cluster_only.py 538 MD10                        # 
python simulate_cluster_only.py 539 MD10                        # 
python simulate_cluster_only.py 540 MD10                        # 
python simulate_cluster_only.py 541 MD10                        # 
python simulate_cluster_only.py 542 MD10                        # 
python simulate_cluster_only.py 543 MD10                        # 
python simulate_cluster_only.py 544 MD10                        # 
python simulate_cluster_only.py 545 MD10                        # 
python simulate_cluster_only.py 546 MD10                        # 
python simulate_cluster_only.py 547 MD10                        # 
python simulate_cluster_only.py 548 MD10                        # 
python simulate_cluster_only.py 549 MD10                        # 
python simulate_cluster_only.py 550 MD10                        # 
                                                                # 
pyCONDA                                                         # 
cd $GIT_AGN_MOCK/python/sixte                                   # 
python simulate_cluster_only.py 551 MD10                        # 
python simulate_cluster_only.py 552 MD10                        # 
python simulate_cluster_only.py 553 MD10                        # 
python simulate_cluster_only.py 554 MD10                        # 
python simulate_cluster_only.py 555 MD10                        # 
python simulate_cluster_only.py 556 MD10                        # 
python simulate_cluster_only.py 557 MD10                        # 
python simulate_cluster_only.py 558 MD10                        # 
python simulate_cluster_only.py 559 MD10                        # 
python simulate_cluster_only.py 560 MD10                        # 
python simulate_cluster_only.py 561 MD10                        # 
python simulate_cluster_only.py 562 MD10                        # 
python simulate_cluster_only.py 563 MD10                        # 
python simulate_cluster_only.py 564 MD10                        # 
python simulate_cluster_only.py 565 MD10                        # 
python simulate_cluster_only.py 566 MD10                        # 
python simulate_cluster_only.py 567 MD10                        # 
python simulate_cluster_only.py 568 MD10                        # 
python simulate_cluster_only.py 569 MD10                        # 
python simulate_cluster_only.py 570 MD10                        # 
python simulate_cluster_only.py 571 MD10                        # 
python simulate_cluster_only.py 572 MD10                        # 
python simulate_cluster_only.py 573 MD10                        # 
python simulate_cluster_only.py 574 MD10                        # 
python simulate_cluster_only.py 575 MD10                        # 
python simulate_cluster_only.py 576 MD10                        # 
python simulate_cluster_only.py 577 MD10                        # 
python simulate_cluster_only.py 578 MD10                        # 
python simulate_cluster_only.py 579 MD10                        # 
python simulate_cluster_only.py 580 MD10                        # 
python simulate_cluster_only.py 581 MD10                        # 
python simulate_cluster_only.py 582 MD10                        # 
python simulate_cluster_only.py 583 MD10                        # 
python simulate_cluster_only.py 584 MD10                        # 
python simulate_cluster_only.py 585 MD10                        # 
python simulate_cluster_only.py 586 MD10                        # 
python simulate_cluster_only.py 587 MD10                        # 
python simulate_cluster_only.py 588 MD10                        # 
python simulate_cluster_only.py 589 MD10                        # 
python simulate_cluster_only.py 590 MD10                        # 
python simulate_cluster_only.py 591 MD10                        # 
python simulate_cluster_only.py 592 MD10                        # 
python simulate_cluster_only.py 593 MD10                        # 
python simulate_cluster_only.py 594 MD10                        # 
python simulate_cluster_only.py 595 MD10                        # 
python simulate_cluster_only.py 596 MD10                        # 
python simulate_cluster_only.py 597 MD10                        # 
python simulate_cluster_only.py 598 MD10                        # 
python simulate_cluster_only.py 599 MD10                        # 
                                                                # 
#!/bin/bash                                                     # 
pyCONDA                                                         # 
cd $GIT_AGN_MOCK/python/sixte                                   # 
# source /home/erosita/sw/sass-setup.sh eSASSusers_190520         # 
python simulate_cluster_only.py 600 MD10                         # 600# ## # # # 
python simulate_cluster_only.py 601 MD10                         # 601# ## # # # 
python simulate_cluster_only.py 602 MD10                         # 602# ## # # # 
python simulate_cluster_only.py 603 MD10                         # 603# ## # # # 
python simulate_cluster_only.py 604 MD10                         # 604# ## # # # 
python simulate_cluster_only.py 605 MD10                         # 605# ## # # # 
python simulate_cluster_only.py 606 MD10                         # 606# ## # # # 
python simulate_cluster_only.py 607 MD10                         # 607# ## # # # 
python simulate_cluster_only.py 608 MD10                         # 609# ## # # # 
python simulate_cluster_only.py 609 MD10                         
python simulate_cluster_only.py 610 MD10                         
python simulate_cluster_only.py 611 MD10                         # 611# ## # # # 
python simulate_cluster_only.py 612 MD10                         # 612# ## # # # 
python simulate_cluster_only.py 613 MD10                         # 613# ## # # # 
python simulate_cluster_only.py 614 MD10                         # 614# ## # # # 
python simulate_cluster_only.py 615 MD10                         # 615# ## # # # 
python simulate_cluster_only.py 616 MD10                         
python simulate_cluster_only.py 617 MD10                         # 617# ## # # # 
python simulate_cluster_only.py 618 MD10                         # 618# ## # # # 
python simulate_cluster_only.py 619 MD10                         
python simulate_cluster_only.py 620 MD10                         
python simulate_cluster_only.py 621 MD10                         
python simulate_cluster_only.py 622 MD10                         
python simulate_cluster_only.py 623 MD10                         
python simulate_cluster_only.py 624 MD10                        # 
python simulate_cluster_only.py 625 MD10                        # 
python simulate_cluster_only.py 626 MD10                        # 
python simulate_cluster_only.py 627 MD10                        # 
python simulate_cluster_only.py 628 MD10                        # 
python simulate_cluster_only.py 629 MD10                        # 
python simulate_cluster_only.py 630 MD10                        # 
python simulate_cluster_only.py 631 MD10                        # 
python simulate_cluster_only.py 632 MD10                        # 
python simulate_cluster_only.py 633 MD10                        # 
python simulate_cluster_only.py 634 MD10                        # 
python simulate_cluster_only.py 635 MD10                        # 
python simulate_cluster_only.py 636 MD10                        # 
python simulate_cluster_only.py 637 MD10                        # 
python simulate_cluster_only.py 638 MD10                        # 
python simulate_cluster_only.py 639 MD10                        # 
python simulate_cluster_only.py 640 MD10                        # 
python simulate_cluster_only.py 641 MD10                        # 
python simulate_cluster_only.py 642 MD10                        # 
python simulate_cluster_only.py 643 MD10                        # 
python simulate_cluster_only.py 644 MD10                        # 
python simulate_cluster_only.py 645 MD10                        # 
python simulate_cluster_only.py 646 MD10                        # 
python simulate_cluster_only.py 647 MD10                        # 
python simulate_cluster_only.py 648 MD10                        # 
python simulate_cluster_only.py 649 MD10                        # 
python simulate_cluster_only.py 650 MD10                        # 
                                                                # 
pyCONDA                                                         # 
cd $GIT_AGN_MOCK/python/sixte                                   # 
python simulate_cluster_only.py 651 MD10                        # 
python simulate_cluster_only.py 652 MD10                        # 
python simulate_cluster_only.py 653 MD10                        # 
python simulate_cluster_only.py 654 MD10                        # 
python simulate_cluster_only.py 655 MD10                        # 
python simulate_cluster_only.py 656 MD10                        # 
python simulate_cluster_only.py 657 MD10                        # 
python simulate_cluster_only.py 658 MD10                        # 
python simulate_cluster_only.py 659 MD10                        # 
python simulate_cluster_only.py 660 MD10                        # 
python simulate_cluster_only.py 661 MD10                        # 
python simulate_cluster_only.py 662 MD10                        # 
python simulate_cluster_only.py 663 MD10                        # 
python simulate_cluster_only.py 664 MD10                        # 
python simulate_cluster_only.py 665 MD10                        # 
python simulate_cluster_only.py 666 MD10                        # 
python simulate_cluster_only.py 667 MD10                        # 
python simulate_cluster_only.py 668 MD10                        # 
python simulate_cluster_only.py 669 MD10                        # 
python simulate_cluster_only.py 670 MD10                        # 
python simulate_cluster_only.py 671 MD10                        # 
python simulate_cluster_only.py 672 MD10                        # 
python simulate_cluster_only.py 673 MD10                        # 
python simulate_cluster_only.py 674 MD10                        # 
python simulate_cluster_only.py 675 MD10                        # 
python simulate_cluster_only.py 676 MD10                        # 
python simulate_cluster_only.py 677 MD10                        # 
python simulate_cluster_only.py 678 MD10                        # 
python simulate_cluster_only.py 679 MD10                        # 
python simulate_cluster_only.py 680 MD10                        # 
python simulate_cluster_only.py 681 MD10                        # 
python simulate_cluster_only.py 682 MD10                        # 
python simulate_cluster_only.py 683 MD10                        # 
python simulate_cluster_only.py 684 MD10                        # 
python simulate_cluster_only.py 685 MD10                        # 
python simulate_cluster_only.py 686 MD10                        # 
python simulate_cluster_only.py 687 MD10                        # 
python simulate_cluster_only.py 688 MD10                        # 
python simulate_cluster_only.py 689 MD10                        # 
python simulate_cluster_only.py 690 MD10                        # 
python simulate_cluster_only.py 691 MD10                        # 
python simulate_cluster_only.py 692 MD10                        # 
python simulate_cluster_only.py 693 MD10                        # 
python simulate_cluster_only.py 694 MD10                        # 
python simulate_cluster_only.py 695 MD10                        # 
python simulate_cluster_only.py 696 MD10                        # 
python simulate_cluster_only.py 697 MD10                        # 
python simulate_cluster_only.py 698 MD10                        # 
python simulate_cluster_only.py 699 MD10                        # 
                                                                # 
#!/bin/bash                                                     # 
pyCONDA                                                         # 
cd $GIT_AGN_MOCK/python/sixte                                   # 
# source /home/erosita/sw/sass-setup.sh eSASSusers_190520         # 
python simulate_cluster_only.py 700 MD10                        # 700# ## # # # # 
python simulate_cluster_only.py 701 MD10                        # 701# ## # # # # 
python simulate_cluster_only.py 702 MD10                        # 702# ## # # # # 
python simulate_cluster_only.py 703 MD10                        # 703# ## # # # # 
python simulate_cluster_only.py 704 MD10                        # 704# ## # # # # 
python simulate_cluster_only.py 705 MD10                        # 705# ## # # # # 
python simulate_cluster_only.py 706 MD10                        # 706# ## # # # # 
python simulate_cluster_only.py 707 MD10                        # 707# ## # # # # 
python simulate_cluster_only.py 708 MD10                        # 
python simulate_cluster_only.py 709 MD10                        # 
python simulate_cluster_only.py 710 MD10                        # 
python simulate_cluster_only.py 711 MD10                        # 
python simulate_cluster_only.py 712 MD10                        # 
python simulate_cluster_only.py 713 MD10                        # 
python simulate_cluster_only.py 714 MD10                        # 
python simulate_cluster_only.py 715 MD10                        # 
python simulate_cluster_only.py 716 MD10                        # 
python simulate_cluster_only.py 717 MD10                        # 
python simulate_cluster_only.py 718 MD10                        # 
python simulate_cluster_only.py 719 MD10                        # 
python simulate_cluster_only.py 720 MD10                        # 
python simulate_cluster_only.py 721 MD10                        # 
python simulate_cluster_only.py 722 MD10                        # 
python simulate_cluster_only.py 723 MD10                        # 
python simulate_cluster_only.py 724 MD10                        # 
python simulate_cluster_only.py 725 MD10                        # 
python simulate_cluster_only.py 726 MD10                        # 
python simulate_cluster_only.py 727 MD10                        # 
python simulate_cluster_only.py 728 MD10                        # 
python simulate_cluster_only.py 729 MD10                        # 
python simulate_cluster_only.py 730 MD10                        # 
python simulate_cluster_only.py 731 MD10                        # 
python simulate_cluster_only.py 732 MD10                        # 
python simulate_cluster_only.py 733 MD10                        # 
python simulate_cluster_only.py 734 MD10                        # 
                                                                # 
python simulate_cluster_only.py 735 MD10                        # 
python simulate_cluster_only.py 736 MD10                        # 
python simulate_cluster_only.py 737 MD10                        # 
python simulate_cluster_only.py 738 MD10                        # 
python simulate_cluster_only.py 739 MD10                        # 
python simulate_cluster_only.py 740 MD10                        # 
python simulate_cluster_only.py 741 MD10                        # 
python simulate_cluster_only.py 742 MD10                        # 
python simulate_cluster_only.py 743 MD10                        # 
python simulate_cluster_only.py 744 MD10                        # 
python simulate_cluster_only.py 745 MD10                        # 
python simulate_cluster_only.py 746 MD10                        # 
python simulate_cluster_only.py 747 MD10                        # 
python simulate_cluster_only.py 748 MD10                        # 
python simulate_cluster_only.py 749 MD10                        # 
python simulate_cluster_only.py 750 MD10                        # 
python simulate_cluster_only.py 751 MD10                        # 
python simulate_cluster_only.py 752 MD10                        # 
python simulate_cluster_only.py 753 MD10                        # 
python simulate_cluster_only.py 754 MD10                        # 
python simulate_cluster_only.py 755 MD10                        # 
python simulate_cluster_only.py 756 MD10                        # 
python simulate_cluster_only.py 757 MD10                        # 
python simulate_cluster_only.py 758 MD10                        # 
python simulate_cluster_only.py 759 MD10                        # 
python simulate_cluster_only.py 760 MD10                        # 
python simulate_cluster_only.py 761 MD10                        # 
python simulate_cluster_only.py 762 MD10                        # 
python simulate_cluster_only.py 763 MD10                        # 
python simulate_cluster_only.py 764 MD10                        # 
python simulate_cluster_only.py 765 MD10                        # 
python simulate_cluster_only.py 766 MD10                        # 
python simulate_cluster_only.py 767 MD10                        # 

