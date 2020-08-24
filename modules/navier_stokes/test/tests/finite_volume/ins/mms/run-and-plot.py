#!/usr/bin/env python3

import mms
df1 = mms.run_spatial('lid-driven.i', 9, y_pp=['L2u'])

fig = mms.ConvergencePlot(xlabel='Element Size ($h$)', ylabel='$L_2$ Error')
fig.plot(df1, label=['L2u'], marker='o', markersize=8, num_fitted_points=3, slope_precision=1)
fig.save('lid-driven.png')
