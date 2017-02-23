import os
from math import pi
import numpy as np
from auxiliaryfunctions import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from matplotlib.widgets import Slider, Button, RadioButtons
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from subprocess import call
from matplotlib.ticker import ScalarFormatter

class Index(object):

    def next(self, event):
        global i
        if i<len(id):
            i += 1
        refresh_fig(i)


    def prev(self, event):
        global i
        if i>0:
            i -= 1
        refresh_fig(i)

    def change_units(self, event):
        global units
        units = not units
        refresh_fig(i)

    def change_ext_law(self, event):
        global extlawN
        extlawN+=1
        if extlawN==len(ext_law): extlawN=0
        refresh_fig(i)


def onpick(event):
    global i
    
    if event.artist!=line: return True
    
    if len(event.ind)==1:
        i = int(event.ind)
        refresh_fig(i)
    else:
        print 'failed to choose one point'

    #fig.canvas.mpl_connect('button_press_event', zoom2dspec)
    return True

def get_spectral_position(event):
    
    if event.inaxes == axFullSpec1D:
        wave_selected=event.xdata
        plot_specflux(i, wave_selected)
        fig.canvas.draw_idle()
    
    if event.inaxes == axZoomSpec1D:
        global pos1, pos2
        if pos1==None: pos1 = [event.xdata, event.ydata]
        else:
            pos2 = [event.xdata, event.ydata]
            EW, wl, fl = ew(pos1, pos2, wavelength_zoom, spectra1Dzoom)
            print 'Jen', EW
            imax.cla()
            #imax.plot(x, f)
            imax.plot(wl, fl)
            imax.set_ylim([-1e-18, +5e-18])
            fig.canvas.draw_idle()
            pos1, pos2 = None, None

    return None

def plot_sfr(ax1, ax2, i):
    
    ax1.cla()
    ax2.cla()
    global line
    
    ax1.set_xlim([8.5,13.5])
    ax2.set_xlim([8.5,13.5])
    ax1.set_ylim([0.5,4.5])
    ax2.set_ylim([-10.0,-7.0])
    ax2.set_xlabel('log($M_*/M_{\odot}$)', size=12)
    ax1.set_ylabel('log(SFR [$M_{\odot}$/year])', labelpad = 15)
    ax2.set_ylabel('log(sSFR [year$^{-1}$])', labelpad = 0)
    ax1.xaxis.set_ticks(np.arange(9.0, 12.5, 1.0))
    ax2.xaxis.set_ticks(np.arange(9.0, 12.5, 1.0))
    ax1.yaxis.set_ticks(np.arange(1.0, 4.5, 1.0))
    ax2.yaxis.set_ticks(np.arange(-10, -7, 1.0))
    
    ax1.errorbar(mass, sfr, xerr = [masslow, masshigh], yerr = [sfrlow, sfrhigh], marker = 'o', ms = 5, linewidth=0, elinewidth=1, color = 'k')
    line, = ax1.plot(massplot, sfrplot, marker = 'o', ms = 5, linewidth=0, color = 'k', picker = 5.0)
    ax2.errorbar(mass, ssfr, xerr = [masslow, masshigh], yerr = [ssfrlow, ssfrhigh], marker = 'o', ms = 5, linewidth=0, color = 'k', elinewidth=1)
    
    if i>=0:
        ax1.plot(mass[i], sfr[i], marker = 'o', ms = 5, color = 'red', linewidth = 0)
        ax2.plot(mass[i], ssfr[i], marker = 'o', ms = 5, color = 'red', linewidth = 0)

def plot_specflux(i, wave_zoom):
    global spectra1Dzoom, wavelength_zoom
    
    axFullSpec1D.cla()
    axFullSpec2D.cla()
    axZoomSpec1D.cla()
    axZoomSpec2D.cla()
    delete_ticks(axFullSpec2D)
    delete_ticks(axZoomSpec2D)
    
    wavelength, spectra1D, cdelt = getspectrum(fn[i])
    zoom_selection = [k for k in range(len(wavelength)) if (wavelength[k]>(wave_zoom-wave_window) and wavelength[k]<(wave_zoom+wave_window))]
    wavelength_zoom = wavelength[zoom_selection[0]:zoom_selection[-1]]
    spectra1Dzoom = spectra1D[zoom_selection[0]:zoom_selection[-1]]
    
    
    axZoomSpec1D.annotate('$\lambda$, [$\AA$]', xy=(0.98, 0), ha='left', va='top', xycoords='axes fraction', fontsize=12)
    axZoomSpec1D.set_ylabel('Flux, [$10^{-18}erg/s/cm^2/\AA$]', size = 12)
    axFullSpec1D.annotate('$\lambda$, [$\AA$]', xy=(0.98, 0), ha='left', va='top', xycoords='axes fraction', fontsize=12)
    axFullSpec1D.set_ylabel('Flux, [$10^{-18}erg/s/cm^2/\AA$]', size = 12)
    axFullSpec1D.xaxis.set_tick_params(width = 0.5, length = 5.5, labelsize = 12)
    axFullSpec1D.yaxis.set_tick_params(width = 0.5, length = 5.5, labelsize = 12)
    axZoomSpec1D.xaxis.set_tick_params(width = 0.5, length = 5.5, labelsize = 12)
    axZoomSpec1D.yaxis.set_tick_params(width = 0.5, length = 5.5, labelsize = 12)
    axFullSpec1D.plot(wavelength, spectra1D)
    axFullSpec1D.fill_betweenx([0.95*np.amin(spectra1D),1.05*np.amax(spectra1D)],wave_zoom-wave_window,wave_zoom+wave_window,color='LimeGreen',alpha=0.15,zorder=-1)
    
    imax.cla()
    imax.plot(wavelength_zoom, spectra1Dzoom)
    speczoom,=axZoomSpec1D.plot(wavelength_zoom,spectra1Dzoom,'-',color='DodgerBlue',lw=1)
    
    imax.set_ylim([-1e-18, +5e-18])
    imax.set_ylabel('Flux, [$10^{-18}erg/s/cm^2/\AA$]', size = 12)
    imax.get_yaxis().get_offset_text().set_visible(False)
    
    imax.set_xlim([wavelength_zoom[0], wavelength_zoom[-1]])
    
    """
    fmt = ScalarFormatter()
    axZoomSpec1D.yaxis.set_major_formatter(fmt)
    axFullSpec1D.yaxis.set_major_formatter(fmt)
    imax.yaxis.set_major_formatter(fmt)
    """
    axZoomSpec1D.get_yaxis().get_offset_text().set_visible(False)
    axFullSpec1D.get_yaxis().get_offset_text().set_visible(False)
    
    
    if int(model[i])>0 and float(age[i])>0:
        tempname = '/Users/yanahusanova/progs/phd/templates/template'+id[i]+'.dat'
        if os.path.isfile(tempname):
            tempx, tempy = np.loadtxt(('./templates/template'+id[i]+'.dat'), unpack=True)
            for k in range(len(tempx)):
                tempx[k]=(1.0+float(zs[i]))*tempx[k]
                tempy[k]=tempy[k]
            axFullSpec1D.plot(tempx,tempy, linewidth = 2.0)
    
    #***CONVERT MAGS IN FLUXES***
    magx, magy, errl, erru = [], [], [], []
    
    for dp in mags[i]:
        if dp[1]>-9.0:
            k = len(magy)
            magx.append(dp[0])
            magy.append(flux(dp[1], dp[0]))
            errl.append(magy[k]*(1-np.power(10,-0.4*dp[2])))
            erru.append(magy[k]*(np.power(10,0.4*dp[2])-1))

    axFullSpec1D.errorbar(magx, magy, yerr=[errl,erru], color = 'red', marker = 'o', ms = 5, elinewidth=2.0, linewidth=0, capthick = 2.0, capsize = 5)
    axZoomSpec1D.errorbar(magx, magy, yerr=[errl,erru], color = 'red', marker = 'o', ms = 5, elinewidth=2.0, linewidth=0, capthick = 2.0, capsize = 5)
    
    """
    #plot spectrum lines
    for line in lines:
        l = lines[line][1]*(1.0+float(zs[i]))
        if (l>(x[0]+100) and l<x[len(x)-1]-100):
            ax.text(l, 3.5e-18, lines[line][0], fontsize = 12, ha = 'center', va = 'bottom', rotation=90)
            ax.axvline(x=l, ymin=0.35, ymax = 0.65, ls='-', linewidth=0.8, color='k', alpha=0.5)
    """

    spectra2D, header = fits.getdata(imfile[i], header=True)
    spectra2Dzoom = spectra2D[:, zoom_selection[0]:zoom_selection[-1]]
    axFullSpec2D.imshow(spectra2D,cmap='hot',vmin=-1,vmax=10,extent=(wavelength[0],wavelength[-1],0,spectra2D.shape[0]),aspect='auto')
    speczoomimage=axZoomSpec2D.imshow(spectra2Dzoom,cmap='hot',vmin=-1,vmax=10,extent=(wavelength_zoom[0],wavelength_zoom[-1],0,spectra2Dzoom.shape[0]),aspect='auto')

    ## LIMITS FULL SPECTRA
    axFullSpec1D.set_ylim(0.95*np.amin(spectra1D),1.05*np.amax(spectra1D))
    axFullSpec1D.set_xlim(wavelength[0],wavelength[-1])

    ## LIMITS ZOOM SPECTRA
    axZoomSpec1D.set_ylim(0.95*np.amin(spectra1Dzoom),1.05*np.amax(spectra1Dzoom))
    axZoomSpec1D.set_xlim(wavelength_zoom[0],wavelength_zoom[-1])

    fig.canvas.mpl_connect('button_press_event',get_spectral_position)

def plot_SED(ax, i, units):
    
    ax.cla()
    
    x, y, cdelt = getspectrum(fn[i])
    
    ax.set_xlabel('$\lambda$, [$\AA$]', size = 12)
    if units: ax.set_ylabel('Mag AB', size = 12)
    else: ax.set_ylabel('Flux, [$10^{-18}erg/s/cm^2/\AA$]', size = 12)
    title = str(id[i]) + ", flag " + str(zf[i])
    ax.set_title(title, size = 14)
    ax.xaxis.set_tick_params(width = 0.5, length = 5.5, labelsize = 12)
    ax.yaxis.set_tick_params(width = 0.5, length = 5.5, labelsize = 12)
    ax.set_xscale('log')
    if units: ax.set_ylim([36,18])
    else: ax.set_ylim([-1.5e-18,6.1e-18])
    ax.set_xlim([1500,50000])
    ax.get_yaxis().get_offset_text().set_visible(False)
    
    if units:
        for k in range(len(x)):
            AB=mAB(x[k],y[k])
            if AB!= None: y[k]=AB
            else: y[k] = 36.0
    
    ax.plot(x,y)

    hastemplate = False

    if int(model[i])>0 and float(age[i])>0:
        tempname = '/Users/yanahusanova/progs/phd/templates/template'+id[i]+'.dat'
        if os.path.isfile(tempname):
            hastemplate = True
            tempx, tempy = np.loadtxt((tempname), unpack=True)
            calzl, calzk = np.loadtxt(ext_law[extlawN], unpack=True)
            
            for k in range(len(tempx)):
                tempx[k]=(1.0+float(zs[i]))*tempx[k]
                calzl[k]=(1.0+float(zs[i]))*calzl[k]
            
            tempk = np.interp(tempx, calzl, calzk)
            for k in range(len(tempx)):
                tempy[k]=tempy[k]*np.power(10, -0.4*float(ebv[i])*tempk[k])
            conv = fnu(float(magi[i]))/specflux(tempx,tempy)
            
            for k in range(len(tempx)):
                tempy[k]=tempy[k]*conv
                if units: tempy[k]=mAB(tempx[k],tempy[k])
            
            ax.plot(tempx,tempy, linewidth = 2.0)

    #***Plot mags***
    if units: err = []
    else: errl, erru = [], []
    magx, magy = [], []
    
    for dp in mags[i]:
        if dp[1]>-9.0:
            magx.append(dp[0])
            if units:
                magy.append(dp[1])
                err.append(dp[2])
                #if hastemplate:
                    #magfit.append(specmag(tempx, tempy, filterpath = '', center = ''))
            else:
                magy.append(flux(dp[1], dp[0]))
                errl.append(magy[-1]*(1-np.power(10,-0.4*dp[2])))
                erru.append(magy[-1]*(np.power(10,0.4*dp[2])-1))

    if units: p3, t1, t2 = ax.errorbar(magx, magy, yerr=err, color = 'red', marker = 'o', ms = 5, elinewidth=2.0, linewidth=0, capthick = 2.0, capsize = 5)
    else:
        p3, t1, t2 = ax.errorbar(magx, magy, yerr=[errl,erru], color = 'red', marker = 'o', ms = 5, elinewidth=2.0, linewidth=0, capthick = 2.0, capsize = 5)

    p1 = plt.axhline(y=0.6,xmin=0,xmax=1,color='blue', label='Spectrum')
    p2 = plt.axhline(y=0.6,xmin=0,xmax=1,color='green', label='Best fit model')
    extra1 = mpatches.Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0, label=("$z=%6.4f$" % float(zs[i])))
    
    ax.legend([p1, p2, p3, extra1],[p1.get_label(), p2.get_label(), 'Photometry', extra1.get_label()], loc='upper left', prop={'size':12})

    global text
    text.remove()
    text = axext.text(-2.5, 0.0, ext_law[extlawN])

def create_templates():
    
    f = open('BC03COMB_MOD.list', 'r')
    modelfiles = f.readlines()
    f.close()
    
    for k in range(len(model)):
        if int(model[k])>0 and float(age[k])>0:
            tempname = '/Users/yanahusanova/progs/phd/templates/template'+id[k]+'.dat'
            if not os.path.isfile(tempname):
                f = open('create_template.sh', 'w')
                s = r"""#!/bin/sh

cd ../../bc03/src/
echo ""
echo ""
echo ""
echo ""
echo ""
echo ""
echo ""
echo ""
echo ""
echo ""
echo ""
echo "%8.6e"
echo "100, 15000, -240, %.9e, %s"
echo "%s"
./galaxevpl ./%s
cd ../../progs/phd/""" % (age[k], fluxi(float(magi[k])), zs[k], tempname, modelfiles[int(model[k])-1].split()[0])
                f.write(s.encode('utf8'))
                f.close()
                call(['./create_template.sh'])

def change_vmin(val):
    global VMIN
    data1 = axFullSpec2D.get_images()[0]
    data1.set_clim(vmin=VMIN*val,vmax=VMAX)
    data2 = axZoomSpec2D.get_images()[0]
    data2.set_clim(vmin=VMIN*val,vmax=VMAX)
    fig.canvas.draw_idle()
    return None

def change_vmax(val):
    global VMAX
    data1 = axFullSpec2D.get_images()[0]
    data1.set_clim(vmin=VMIN,vmax=VMAX*val)
    data2 = axZoomSpec2D.get_images()[0]
    data2.set_clim(vmin=VMIN,vmax=VMAX*val)
    fig.canvas.draw_idle()
    return None

def add_sliders(z_lim):
    vmin.cla()
    vmax.cla()
    zchange.cla()
    slider_vmin = Slider(vmin, 'Vmin', 0.0, 1.0, valinit=1.0, color='Gold', closedmin=False)
    slider_vmax = Slider(vmax, 'Vmax', 0.0, 1.0, valinit=1.0, color='Gold', closedmin=False)
    slider_zchange = Slider(zchange, 'z_lim', 5.0, 7.0, valinit=z_lim, color='green', closedmin=False)
    slider_vmin.on_changed(change_vmin)
    slider_vmax.on_changed(change_vmax)
    slider_zchange.on_changed(change_z_lim)
    return slider_vmin, slider_vmax, slider_zchange

def refresh_fig(i):
    global Vminslider, Vmaxslider, zchangeslider
    
    plot_SED(specax, i, units)
    plot_sfr(sfrmass, ssfrmass, i)
    fig.canvas.mpl_connect('pick_event', onpick)
    plot_specflux(i, wave_zoom)
    fig.canvas.draw_idle()
    Vminslider, Vmaxslider, zchangeslider = add_sliders(z_lim)


def change_z_lim(val):
    
    global i, z_lim
    global id, fn, imfile, zs, zf, magi, mags, model, mass, masslow, masshigh, sfr, sfrlow, sfrhigh, ssfr, ssfrlow, ssfrhigh, age, agelow, agehigh, massplot, sfrplot
    i, z_lim = -1, val
    id, fn, imfile, zs, zf, magi, ebv, mags, model, mass, masslow, masshigh, sfr, sfrlow, sfrhigh, ssfr, ssfrlow, ssfrhigh, age, agelow, agehigh, massplot, sfrplot = data_at_z(val)
    plot_sfr(sfrmass, ssfrmass, i)
    fig.canvas.draw_idle()
    #refresh_fig(i)


if __name__ == "__main__":
    lines = define_lines()
    z_lim = 2.0
    pos1, pos2 = None, None
    units = True #if units = True, spec+SED will be plotted in mags
    ext_law = ['calzetti.dat','SMC_prevot.dat']
    #ext_law = ['calzetti_mod.dat', 'calzetti_modified.dat', 'calzetti.dat', 'extinc_ctio.dat', 'extinc_eso.dat', 'LMC_Fitzpatrick_noBump.dat', 'LMC_Fitzpatrick.dat', 'MW_Allen.dat', 'MW_seaton.dat', 'SB_calzetti_2200bump.dat', 'SB_calzetti_bump.dat', 'SB_calzetti_bump1.dat', 'SB_calzetti_bump2.dat', 'SB_calzetti_bump3.dat', 'SB_calzetti_mod.dat', 'SB_calzetti.dat', 'SMC_prevot.dat]
    extlawN=0
    
    id, fn, imfile, zs, zf, magi, ebv, mags, model, mass, masslow, masshigh, sfr, sfrlow, sfrhigh, ssfr, ssfrlow, ssfrhigh, age, agelow, agehigh, massplot, sfrplot = data_at_z(z_lim)

    create_templates()

    i = 0
    VMIN = -34
    VMAX = +34
    wave_zoom = 5000
    wave_window = 300
    spectra1Dzoom, wavelength_zoom = 1, 1
    
    fig = plt.figure(figsize=(16,9.0), facecolor='white')
    sfrmass = fig.add_axes([0.05,0.73,0.25,0.25])
    ssfrmass = fig.add_axes([0.05,0.48,0.25,0.25],sharex=sfrmass)
    specax = fig.add_axes([0.05,0.07,0.65,0.32])
    axFullSpec2D = fig.add_axes([0.35,0.93,0.63,0.05])
    axFullSpec1D = fig.add_axes([0.35,0.78,0.63,0.15],sharex=axFullSpec2D)
    axZoomSpec2D = fig.add_axes([0.35,0.65,0.35,0.13])
    axZoomSpec1D = fig.add_axes([0.35,0.45,0.35,0.2],sharex=axZoomSpec2D)
    axprev = plt.axes([0.05, 0.03, 0.05, 0.035])
    axnext = plt.axes([0.11, 0.03, 0.05, 0.035])
    axunits = plt.axes([0.65, 0.03, 0.05, 0.035])
    axext = plt.axes([0.57, 0.03, 0.05, 0.035])
    vmin = plt.axes([0.77, 0.71, 0.18, 0.025])
    vmax = plt.axes([0.77, 0.67, 0.18, 0.025])
    zchange = plt.axes([0.77, 0.63, 0.18, 0.025])
    imax = plt.axes([0.77, 0.07, 0.18, 0.32])
    text = axext.text(-2.5, 0.0, ext_law[extlawN])

    line = 1
    callback = Index()
    bnext = Button(axnext, 'Next')
    bprev = Button(axprev, 'Previous')
    bunits = Button(axunits, 'Flux/Mag')
    bext = Button(axext, 'Ext_law')
    bnext.on_clicked(callback.next)
    bprev.on_clicked(callback.prev)
    bunits.on_clicked(callback.change_units)
    bext.on_clicked(callback.change_ext_law)
    Vminslider, Vmaxslider, zchangeslider = 0, 0, 0
    
    refresh_fig(i)

    plt.show()
