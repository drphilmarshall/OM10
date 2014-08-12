from numpy import *
'''
To convert Adri's mathematica version to Python

Kai 
'''

def bs(n):
    return (2.*n-1.)/3.

def Sersic(R, n):
    return exp(-bs(n)*R**(1./n))

def flase(x, y, flat, pa, n):
    return Sersic((flat*(x*cos(pa)+y*sin(pa))**2.+flat**(-1.0)*(-sin(pa)*x+cos(pa)*y)**2.)**0.5, n)

def G(x, dx):
    return exp(-0.5*x**2./dx**2.)/((2*pi)**0.5*dx)

def GG(x, dx):
    return G(abs(x)**0.5, dx)

def Gint(x, y, dx, dy):
    return (9./16)*G(x,dx)*G(y,dy)+(3./32)*(G(x+1.,dx)*G(y,dy)+G(x-1.,dx)*G(y,dy)+\
           G(x,dx)*G(y+1.,dy)+G(x,dx)*G(y-1.,dy))+(1/64.)*(G(x-1.,dx)*G(y-1.,dy)+\
           G(x-1.,dx)G(y+1.,dy)+G(x+1.,dx)*G(y-1.,dt)+G(x+1.,dx)*G(y+1.,dy))
#Gint is useful to interpolate the Gaussian psf on 3*3 grid.
#SDSS
pixscale = 0.4
meanIQ = 1.4/2
meanIQ = meanIQ/(log(2.)*2.**0.5)    # make sure the log is log_e!
meandepth = 20.8 #magnitudes per arcsecond
errdepth = 0.3
#more specific: band fluxes and fluctuations
gmean = 21.9
egd = 0.3
gsky = pixscale**2.*10.**(9.-0.4*gmean)
rmean = 20.9
erd = 0.3
rsky = pixscale**2.*10.**(9.-0.4*rmean)
imean = 20.2
eid = 0.4
isky = pixscale**2.*10.**(9.-0.4*imean)
zmean = 18.9
ezd = 0.5
zsky = pixscale**2.*10.**(9.-0.4*zmean)
#psf width distributions
mgIQ = 1.65/(2.*2.**0.5*log(2.))
dgIQ = 0.4/(2.*2.**0.5*log(2.))
moIQ = 1.4/(2.*2.**0.5*log(2.))
doIQ = 0.3/(2.*2.**0.5*log(2.))  #psf width in the other bands

dr = pixscale**2.*10.**(9.-0.4*meandepth)/5. #five sigma detection of deepest source
expo = (log(10.)*erd/(2.5*(2*pi)**0.5))/dr**2.
dg = (log(10.)*egd/(2.5*(2*pi)**0.5))**0.5/expo**0.5
di = (log(10.)*eid/(2.5*(2*pi)**0.5))**0.5/expo**0.5
dz = (log(10.)*ezd/(2.5*(2*pi)**0.5))**0.5/expo**0.5








