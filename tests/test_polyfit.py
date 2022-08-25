from pyspextool.fit.poly_fit_1d import poly_fit_1d
import numpy as np
from pyspextool.utils.for_print import for_print
import matplotlib.pyplot as pl

#x = np.array([10. ,20. ,30. ,40. ,50., 60. ,70. ,80. ,90.])
#y = np.array([0.37, 0.58, 0.83, 1.15, 1.36, 1.62, 1.90, 2.18, 2.45])
#yerr = np.array([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
#
#goodbad = np.ones(9)
##goodbad[1] = 0
##goodbad[2] = 2
#result = poly_fit_1d(x,y,1,yunc=yerr,goodbad=goodbad, silent=False,
#                     robust={'thresh':4, 'eps':0.1})






N = 100
idx = np.arange(N)
x = np.linspace(0, 4, N)
y = x**3 - 6*x**2 + 12*x - 9 
y += np.random.normal(0, 2, size=len(y))
sigma = 1.5
yerr = np.ones(N)*sigma 
y[50] = 30
y[30] = np.nan

#goodbad = np.ones(N)
#goodbad[30] = 2

result = poly_fit_1d(x, y, 3, yunc=yerr, robust={'thresh':4, 'eps':0.1},
                    silent=False)


yfit = result['yfit']
pl.plot(x,y,'o')
pl.plot(x,yfit,'r-')
z = np.where(result['goodbad'] == 0)
pl.plot(x[z],y[z],'ro')
pl.show()

