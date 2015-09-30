#!/usr/bin/python
#-*-coding:utf8-*-

from scipy.stats import linregress as lr
from scipy import zeros, random, stats, std, sqrt, log, arange

def vectoriza( seq ):
    """
    Vectorizacion del cogido genetico segun Yu and Chen (2000)
    """
    N = len( seq )
    vec = zeros( N, int )
    a = (-2, -1, 1, 2 )
    b = 'ACGT'
    for i in range( N ):
        idx = b.find( seq[i] )
        vec[i] = a[idx]
    return vec

def Hurst( seq , vectoriza=None ):
    """
    Obtiene el exponente de Hurst de alguna secuencia 
    """
    N = len( seq )
    if vectoriza==None: 
        X = vectoriza( seq )
    else:
        X = seq
    mn = X.mean()
    Y = X-mn
    Z = Y.cumsum()
    R,S,E = zeros(Z.size), zeros(Z.size), zeros(Z.size) 
    for i in range(1,len(Z)+1):
        R[i] = Z[:i].max() - Z[:i].min()
        S[i] = seq[:i].std()
        E[i] = R[i]/S[i]
    for i in range(1,len(E)):
        Em[i] = E[:i].mean()
    lnEm = log(Em[1:])
    lnX = log(range(1,len(Em)))
    m,b,r,p,s = linregress(lnX,lnEm)
    return m,r,p
