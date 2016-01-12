#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt

l2OmegaOld = -1
l2OmegaOld = -1
h1OmegaOld = -1
h1GammaOld = -1

hOld = -1
tauOld = -1

def eoc( a, b ):
    return np.log( a / b )

for i in range( 5 ):
    try:
        filename = 'poisson--{0}-.txt'.format( i )
        my_data = np.genfromtxt( filename , delimiter='  ')

        h = -1
        tau = -1
        file = open( filename )
        lines = file.readlines()
        for line in lines:
            if '# h:' in line:
                data = line.split()
                h = float( data[-1] )
            if '# tau:' in line:
                data = line.split()
                tau = float( data[-1] )
        file.close()
    except:
        print 'unable to open', filename
        continue

    # extract data
    time    = [ d[0] for d in my_data ]
    l2Omega = [ d[1] for d in my_data ]
    h1Omega = [ d[2] for d in my_data ]
    l2Gamma = [ d[3] for d in my_data ]
    h1Gamma = [ d[4] for d in my_data ]

    # plot
    plt.plot( time, l2Omega, 'b.-' )
    plt.semilogy( time, l2Gamma, 'r.-' )

    # print results
    if i == 0:
        print '    h         tau      L2 Omega  (eoc h)  L2 Gamma   (eoc h)'

    if i > 0 and h >= 0:
        eocOmegah = eoc( l2Omega[-1], l2OmegaOld ) / eoc( h, hOld )
        eocGammah = eoc( l2Gamma[-1], l2GammaOld ) / eoc( h, hOld )
        eocOmegaTau = eoc( l2Omega[-1], l2OmegaOld ) / eoc( tau, tauOld )
        eocGammaTau = eoc( l2Gamma[-1], l2GammaOld ) / eoc( tau, tauOld )

        print '{0:7.4e} {1:7.4e} {2:7.4e} {3:7.4f} {4:7.4f} {5:7.4e} {6:7.4f} {7:7.4f}'.format( h, tau,
                                                                                       l2Omega[-1], eocOmegah, eocOmegaTau, \
                                                                                       l2Gamma[-1], eocGammah, eocGammaTau )
    else:
        print '{0:7.4e} {1:7.4e} {2:7.4e} {3} {4} {5:7.4e} {6} {7}'.format( h, tau, l2Omega[-1], '  ---  ', '  ---  ', \
                                                                    l2Gamma[-1], '  ---  ', '  ---  ' )

    l2OmegaOld = l2Omega[-1]
    l2GammaOld = l2Gamma[-1]
    h1OmegaOld = h1Omega[-1]
    h1GammaOld = h1Gamma[-1]
    hOld = h
    tauOld = tau

plt.show()
