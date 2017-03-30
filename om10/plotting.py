# ======================================================================

# Globally useful modules, imported here and then accessible by all
# functions in this file:

# Fonts, latex:
import matplotlib
matplotlib.rc('font',**{'family':'serif', 'serif':['TimesNewRoman']})
matplotlib.rc('text', usetex=True)
import corner

import pylab, sys, numpy as np


# ======================================================================

def plot_sample(sample, saveImg=False, fig=None, color='black',
                parameters=('MAGI','IMSEP','VELDISP','ZLENS','ZSRC')):
    """
    Given an OM10 sample, make a corner plot of the required quantities.

    Parameters
    ----------
    parameters : str, tuple
        Names of the lens parameters to plot
    saveImg : bool
        If true, save image with standardized name.
    IQ : float
        Image quality, for reference.
    fig : matplotlib figure object
        Overlay plot on an existing figure

    Returns
    -------
    fig : matplotlib figure object
        New or updated figure
    """

    features, labels = extract_features(sample, parameters)

    if fig is None:
        fig = corner.corner(features, labels=labels, color=color, smooth=1.0)
    else:
        _   = corner.corner(features, labels=labels, color=color, smooth=1.0, fig=fig)

    if saveImg:
        pngfile = "om10_sample.png"
        pylab.savefig(pngfile)
        print "OM10: Sample plot saved to file:", pngfile

    return fig

# ======================================================================

def extract_features(x, names):
    """
    Given an OM10 table of lenses, extract the required parameters and
    provide labels for them.

    Parameters
    ----------
    x : Table
        OM10 lens sample.
    names : str, tuple
        Names of features required.

    Returns
    -------
    features : float, ndarray
        Values of requested features, for each lens in the Table
    labels : str, list
        Corresponding axis labels
    """

    features = np.array([])
    labels = []

    p = len(names)
    n = len(x)

    for name in names:
        features = np.append(features, x[name])
        labels.append(axis_labels[name])

    return features.reshape(p,n).transpose(), labels

# ======================================================================

def plot_lens(lens, saveImg=False, IQ=0.7):
    """
    Given an OM10 lens, compute some basic quantities
    and use them to plot a cartoon visualization of the lens.

    Parameters
    ----------
    saveImg : bool
        If true, save image with standardized name.
    IQ : float
        Image quality, for reference.
    """

    # # Force matplotlib to not use any Xwindows backend:
    # if saveImg:
    #     try: matplotlib.use('Agg')
    #     except: pass
    # else:
    #     try: matplotlib.use('TkAgg')
    #     except: pass

    # Pull out data for ease of use:
    id = lens['LENSID'][0]
    xi = lens['XIMG'][0]
    yi = lens['YIMG'][0]
    nim = lens['NIMG'][0]
    mui = lens['MAG'][0]
    md = lens['APMAG_I'][0]
    ms = lens['MAGI_IN'][0]
    xs = lens['XSRC'][0]
    ys = lens['YSRC'][0]
    xd = 0.0
    yd = 0.0
    zd = lens['ZLENS'][0]
    zs = lens['ZSRC'][0]
    q = 1.0 - lens['ELLIP'][0]
    phi = lens['PHIE'][0]

    print "OM10: Plotting image configuration of lens ID ",id

    # Compute image magnitudes:
    mi = np.zeros(nim)
    lfi = np.zeros(nim)
    for i in range(nim):
      mi[i] = ms - 2.5*np.log10(np.abs(mui[i]))
      lfi[i] = 0.4*(24-mi[i])
    print "OM10: lens, image magnitudes:",md,mi
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

    # Circle to represent seeing:
    cir = pylab.Circle((1.5,-1.5), radius=IQ/2.0, alpha=0.1, fc='grey')
    pylab.gca().add_patch(cir)
    text = '{:3.1f}" seeing'.format(IQ)
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
    if saveImg:
        pngfile = "om10_qso_ID="+str(id)+".png"
        pylab.savefig(pngfile)
        print "OM10: Lens plot saved to file:",pngfile

# ======================================================================

axis_labels = {}
axis_labels['ZLENS'] = '$z_d$'
axis_labels['VELDISP'] = '$\sigma_d$ / km/s'
axis_labels['ELLIP'] = '$\epsilon_d$'
axis_labels['PHIE'] = '$\phi_d$ / km/s'
axis_labels['GAMMA'] = '$\gamma$'
axis_labels['PHIG'] = '$\phi_{\gamma}$'
axis_labels['ZSRC'] = '$z_s$'
axis_labels['MAGI'] = '$i_3$'
axis_labels['MAGI_IN'] = '$i_s$'
axis_labels['IMSEP'] = '$\Delta \\theta$ / arcsec'
axis_labels['i_SDSS'] = '$i_{\\rm SDSS}$ (AB mag)'
axis_labels['ug'] = '$u-g$ color'
axis_labels['gr'] = '$g-r$ color'
axis_labels['ri'] = '$r-i$ color'
axis_labels['iz'] = '$i-z$ color'
axis_labels['ug'] = '$u-g$ color'
