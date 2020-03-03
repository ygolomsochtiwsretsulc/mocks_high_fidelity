# on ds52
import os
import glob
import sys
import numpy as n
lll = sorted(n.array(glob.glob('all_*.fits')))
for ll in lll:
    print(ll,)
    out = os.system("fkeyprint " + ll + " NAXIS2 | tail --lines=1")

lll = n.array(glob.glob('replicated_*/all_*.list'))
lll.sort()
for ll in lll:
    print(ll,)
    out = os.system(" wc -l " + ll)
