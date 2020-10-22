#!/usr/bin/env python3

import mms
import unittest
from mooseutils import fuzzyEqual, fuzzyAbsoluteEqual

labels = ['L2u', 'L2v', 'L2p']
df1 = mms.run_spatial('2d-rc.i', [5,6,7,8,9], y_pp=labels, mpi=32, console=True)

fig = mms.ConvergencePlot(xlabel='Element Size ($h$)', ylabel='$L_2$ Error')
fig.plot(df1, label=labels, marker='o', markersize=8, num_fitted_points=3, slope_precision=1)
fig.save('2d-rc.png')
for key,value in fig.label_to_slope.items():
    print("%s, %f" % (key, value))
