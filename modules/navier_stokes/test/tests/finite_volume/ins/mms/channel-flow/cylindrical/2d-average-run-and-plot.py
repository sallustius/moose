#!/usr/bin/env python3

import mms

labels = ['L2u', 'L2v', 'L2p']
df1 = mms.run_spatial('2d-average.i', 6, y_pp=labels)

fig = mms.ConvergencePlot(xlabel='Element Size ($h$)', ylabel='$L_2$ Error')
fig.plot(df1, label=labels, marker='o', markersize=8, num_fitted_points=3, slope_precision=1)
fig.save('2d-average.png')
