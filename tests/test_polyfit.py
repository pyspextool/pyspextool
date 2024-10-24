from pyspextool.fit.polyfit import polyfit_1d
import numpy as np
# import matplotlib.pyplot as pl


N = 100
idx = np.arange(N)
x = np.linspace(0, 4, N)
y = x**3 - 6*x**2 + 12*x - 9
y += np.random.normal(0, 2, size=len(y))
sigma = 1.5
yerr = np.ones(N)*sigma
y[50] = 30
y[30] = np.nan

def test_polyfit_1d():

    result = polyfit_1d(x, y, 3, yunc=yerr, robust={'thresh':4, 'eps':0.1},
                    silent=False)

    fit = result['yfit']
    # TODO: Test the results, which don't seem to be reproducible

    z = np.where(result['goodbad'] == 0)
    assert x[z] == 2.0202020202020203
    assert y[z][0] == 30.0

    # Plot the fit
    # pl.plot(x,y,'o')
    # pl.plot(x,yfit,'r-')
    # pl.plot(x[z],y[z],'ro')
    # pl.show()

