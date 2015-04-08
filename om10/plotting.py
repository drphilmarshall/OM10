# ======================================================================

# Globally useful modules, imported here and then accessible by all
# functions in this file:

import matplotlib
# Force matplotlib to not use any Xwindows backend:
try: matplotlib.use('Agg')
except: pass

# Fonts, latex:
matplotlib.rc('font',**{'family':'serif', 'serif':['TimesNewRoman']})
matplotlib.rc('text', usetex=True)

import pylab,sys,numpy as np

import om10

# ======================================================================

def plot_lens(lens):

    USAGE = """
    NAME
      plot_lens

    PURPOSE
      Given an OM10 lens, compute some basic quantities,
      then plot them on the sky.

    COMMENTS

    AUTHORS
      This file is part of the OM10 project, distributed under the
      GPL v2 by Phil Marshall (KIPAC).
      Please cite: Oguri & Marshall (2010), MNRAS, 405, 2579.

    HISTORY
      2010-06-13 started as standalone script Marshall (KIPAC)
      2013-09-27 adapted for OM10 project Marshall (KIPAC)
    """

    # --------------------------------------------------------------------

    # Pull out data for ease of use:
    id = lens.LENSID[0]
    xi = lens.XIMG[0]
    yi = lens.YIMG[0]
    mui = lens.MAG[0]
    nim = lens.NIMG[0]
    md = lens.APMAG_I[0]
    ms = lens.MAGI_IN[0]
    xs = lens.XSRC[0]
    ys = lens.YSRC[0]
    xd = 0.0
    yd = 0.0
    zd = lens.ZLENS[0]
    zs = lens.ZSRC[0]
    q = 1.0 - lens.ELLIP[0]
    phi = lens.PHIE[0]

    print "om10.plot_lens: plotting image configuration of lens ID ",id

    # Compute image magnitudes:
    mi = np.zeros(nim)
    lfi = np.zeros(nim)
    for i in range(nim):
      mi[i] = ms - 2.5*np.log10(np.abs(mui[i]))
      lfi[i] = 0.4*(24-mi[i])
    print "om10.plot_lens: lens, image magnitudes:",md,mi
    lfd = 0.4*(24-md)
    # print "om10.plot_lens: lens, image log fluxes:",lfd,lfi

    # ------------------------------------------------------------------
    # Compute caustics and critical curves:




    # ------------------------------------------------------------------

    # Start figure:
    fig = pylab.figure(figsize=(8,8))
    # ,aspect='equal')

    # Axes limits, useful sizes:
    xmax = 1.99
    dm = 1.0/10

    # Plot command sets its own axes. 'bp' = blue pentagons
    # pylab.plot(xi, yi, 'bp')
    pylab.plot(xi, yi, color='blue', \
               marker='+', markersize=10, markeredgewidth=2, \
               linestyle='')
    pylab.plot(xs, ys, color='lightblue', \
               marker='+', markersize=10, markeredgewidth=2, \
               linestyle='')
    pylab.plot(xd, yd, color='orange', \
               marker='+', markersize=10, markeredgewidth=2, \
               linestyle='')

    # Ellipse to represent lens brightness:
    ell = matplotlib.patches.Ellipse((xd,yd), width=2*dm*lfd, height=2*q*dm*lfd, angle=phi, alpha=0.2, fc='orange')
    pylab.gca().add_patch(ell)

    # Circles to represent image brightness:
    for i in range(nim):
      cir = pylab.Circle((xi[i],yi[i]), radius=dm*lfi[i], alpha=0.2, fc='blue')
      pylab.gca().add_patch(cir)

    # Circle to represent 0.7" seeing:
    cir = pylab.Circle((1.5,-1.5), radius=0.7/2.0, alpha=0.1, fc='grey')
    pylab.gca().add_patch(cir)
    text = '0.7" seeing'
    pylab.annotate(text, (370,5), xytext=None, fontsize=14, \
                     xycoords='axes points',textcoords='axes points')

    # Legend giving lens, source redshift:
    text1 = "$z_d$ = %5.2f" % zd
    text2 = "$z_s$ = %5.2f" % zs
    pylab.annotate(text1, (10,430), xytext=None, fontsize=14, \
                     xycoords='axes points',textcoords='axes points')
    pylab.annotate(text2, (10,410), xytext=None, fontsize=14, \
                     xycoords='axes points',textcoords='axes points')

    # Plot title:
    title = "OM10 lensed QSO, ID="+str(id)
    pylab.title(title,fontsize=20)
    # Set axes labels:
    pylab.xlabel("x / arcsec",fontsize=20)
    pylab.ylabel("y / arcsec",fontsize=20)
    # Set axis limits:
    pylab.axis([-xmax,xmax,-xmax,xmax])
    # Add a grid:
    pylab.grid(color='grey', linestyle='--', linewidth=0.5)

    # Plot graph to file:
    pngfile = "om10_qso_ID="+str(id)+".png"
    pylab.savefig(pngfile)
    print "om10.plot_lens: figure saved to file:",pngfile


# ======================================================================

if __name__ == '__main__':

    db = om10.DB(catalog="data/qso_mock.fits")

    # Pull out a specific lens and plot it:
    id = 7176527
    lens = db.get_lens(id)
    om10.plot_lens(lens)

    # Plot 3 random lenses and plot them:
    lenses = db.select_random(maglim=21.4,area=30000.0,IQ=1.0,Nlens=3)
    if lenses is not None:
        for id in lenses.LENSID:
            lens = db.get_lens(id)
            om10.plot_lens(lens)


# ======================================================================
