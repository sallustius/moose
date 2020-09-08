#!/usr/bin/env python3

import mms
df1 = mms.run_spatial('rc.i', 8, y_pp=['L2u', 'L2v', 'L2p'], mpi=8)

fig = mms.ConvergencePlot(xlabel='Element Size ($h$)', ylabel='$L_2$ Error')
fig.plot(df1, label=['L2u', 'L2v', 'L2p'], marker='o', markersize=8, num_fitted_points=3, slope_precision=1)
fig.save('rc.png')
