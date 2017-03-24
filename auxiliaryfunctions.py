from __future__ import division # confidence high
import os
import glob
import numpy as np
from scipy.integrate import simps
#import matplotlib.pyplot as plt
#from matplotlib._png import read_png
from astropy.io import fits
import collections
from scipy.integrate import simps
from scipy.integrate import quad
import math
c_AAs     = 2.99792458e18
c_m = c_AAs*1e-10
sigmal = 0.76
sigmam = 0.24
#from matplotlib.offsetbox import OffsetImage, AnnotationBbox


#zscale algorith from http://stsdas.stsci.edu/stsci_python_epydoc/numdisplay/numdisplay.zscale-pysrc.html#zscale
MAX_REJECT = 0.5
MIN_NPIXELS = 5
GOOD_PIXEL = 0
BAD_PIXEL = 1
KREJ = 2.5
MAX_ITERATIONS = 5

def zscale(image, nsamples=1000, contrast=0.25, bpmask=None, zmask=None):
    """Implement IRAF zscale algorithm nsamples=1000 and contrast=0.25 are the IRAF display task defaults bpmask and zmask not implemented yet image is a 2-d numpy array returns (z1, z2) """

    # Sample the image
    samples = zsc_sample (image, nsamples, bpmask, zmask)
    npix = len(samples)
    samples.sort()
    zmin = samples[0]
    zmax = samples[-1]
    # For a zero-indexed array
    center_pixel = (npix - 1) // 2
    if npix%2 == 1:
        median = samples[center_pixel]
    else:
        median = 0.5 * (samples[center_pixel] + samples[center_pixel + 1])

    #
    # Fit a line to the sorted array of samples
    minpix = max(MIN_NPIXELS, int(npix * MAX_REJECT))
    ngrow = max (1, int (npix * 0.01))
    ngoodpix, zstart, zslope = zsc_fit_line (samples, npix, KREJ, ngrow, MAX_ITERATIONS)

    if ngoodpix < minpix:
        z1 = zmin
        z2 = zmax
    else:
        if contrast > 0: zslope = zslope / contrast
        z1 = max (zmin, median - (center_pixel - 1) * zslope)
        z2 = min (zmax, median + (npix - center_pixel) * zslope)
    return z1, z2

def zsc_sample(image, maxpix, bpmask=None, zmask=None):
    # Figure out which pixels to use for the zscale algorithm
    # Returns the 1-d array samples
    # Don't worry about the bad pixel mask or zmask for the moment
    # Sample in a square grid, and return the first maxpix in the sample
    nc = image.shape[0]
    nl = image.shape[1]
    stride = max (1.0, math.sqrt((nc - 1) * (nl - 1) / float(maxpix)))
    stride = int (stride)
    samples = image[::stride,::stride].flatten()
    samples = samples[~np.isnan(samples)]
    return samples[:maxpix]

def zsc_fit_line (samples, npix, krej, ngrow, maxiter): 
    #
    # First re-map indices from -1.0 to 1.0
    xscale = 2.0 / (npix - 1)
    xnorm = np.arange(npix)
    xnorm = xnorm * xscale - 1.0
    
    ngoodpix = npix
    minpix = max (MIN_NPIXELS, int (npix*MAX_REJECT))
    last_ngoodpix = npix + 1

    # This is the mask used in k-sigma clipping.  0 is good, 1 is bad
    badpix = np.zeros(npix, dtype="int32")

    #
    #  Iterate

    for niter in range(maxiter):

        if (ngoodpix >= last_ngoodpix) or (ngoodpix < minpix):
            break

        # Accumulate sums to calculate straight line fit
        goodpixels = np.where(badpix == GOOD_PIXEL)
        sumx = xnorm[goodpixels].sum()
        sumxx = (xnorm[goodpixels]*xnorm[goodpixels]).sum()
        sumxy = (xnorm[goodpixels]*samples[goodpixels]).sum()
        sumy = samples[goodpixels].sum()
        sum = len(goodpixels[0])

        delta = sum * sumxx - sumx * sumx
        # Slope and intercept
        intercept = (sumxx * sumy - sumx * sumxy) / delta
        slope = (sum * sumxy - sumx * sumy) / delta
        
        # Subtract fitted line from the data array
        fitted = xnorm*slope + intercept
        flat = samples - fitted

        # Compute the k-sigma rejection threshold
        ngoodpix, mean, sigma = zsc_compute_sigma (flat, badpix, npix)
    
        threshold = sigma * krej
    
        # Detect and reject pixels further than k*sigma from the fitted line
        lcut = -threshold
        hcut = threshold
        below = np.where(flat < lcut)
        above = np.where(flat > hcut)

        badpix[below] = BAD_PIXEL
        badpix[above] = BAD_PIXEL

        # Convolve with a kernel of length ngrow
        kernel = np.ones(ngrow,dtype="int32")
        badpix = np.convolve(badpix, kernel, mode='same')

        ngoodpix = len(np.where(badpix == GOOD_PIXEL)[0])

        niter += 1
    
    # Transform the line coefficients back to the X range [0:npix-1]
    zstart = intercept - slope
    zslope = slope * xscale
    
    return ngoodpix, zstart, zslope

def zsc_compute_sigma (flat, badpix, npix):

    # Compute the rms deviation from the mean of a flattened array.
    # Ignore rejected pixels

    # Accumulate sum and sum of squares
    goodpixels = np.where(badpix == GOOD_PIXEL)
    sumz = flat[goodpixels].sum()
    sumsq = (flat[goodpixels]*flat[goodpixels]).sum()
    ngoodpix = len(goodpixels[0])
    if ngoodpix == 0:
        mean = None
        sigma = None
    elif ngoodpix == 1:
        mean = sumz
        sigma = None
    else:
        mean = sumz / ngoodpix
        temp = sumsq / (ngoodpix - 1) - sumz*sumz / (ngoodpix * (ngoodpix - 1))
        if temp < 0:
            sigma = 0.0
        else:
            sigma = math.sqrt (temp)

    return ngoodpix, mean, sigma

cwldict = {'u_au': [3823.29], 'g_au': [4720.0], 'r_au': [6404.0], 'i_au': [7764.0], 'z_au': [9256.0], 'j_au': [12520.0], 'h_au': [16310.0], 'k_au': [21460.0], 'FUV': [1539.78], 'NUV': [2313.89], 'u': [3823.29], 'B': [4458.34], 'V': [5477.83], 'g': [4777.07], 'r': [6288.71], 'i': [7683.88], 'z': [9036.88], 'J': [12534.65], 'K': [21539.88], 'ks_ISAAC': [21539.88], 'Ks_HAWKI': [21539.88], 'irac1_servs': [35634.28], 'irac2_servs': [45110.13], 'I1': [35634.28], 'I2': [45110.13], 'I3': [57593.39], 'I4': [79594.95], 'U_MUSYC': [3650], 'U38_MUSYC': [3800], 'B_MUSYC': [4450], 'V_MUSYC': [5510], 'R_MUSYC': [6580], 'I_MUSYC': [8060], 'z_MUSYC': [9000], 'J_MUSYC': [12200], 'H_MUSYC': [16300], 'K_MUSYC': [21900], 'IA427_MUSYC': [4263.45], 'IA464_MUSYC': [4635.13], 'IA445_MUSYC': [4450], 'IA484_MUSYC': [4849.2], 'IA505_MUSYC': [5050], 'IA527_MUSYC': [5270], 'IA550_MUSYC': [5500], 'IA574_MUSYC': [5740], 'IA598_MUSYC': [5980], 'IA624_MUSYC': [6240], 'IA651_MUSYC': [6510], 'IA679_MUSYC': [6790], 'IA709_MUSYC': [7090], 'IA738_MUSYC': [7380], 'IA767_MUSYC': [7970], 'IA797_MUSYC': [7970], 'IA827_MUSYC': [8270], 'IA856_MUSYC': [8560], 'I1_SIMPLE': [35634.28], 'I2_SIMPLE': [45110.13], 'I3_SIMPLE': [57593.39], 'I4_SIMPLE': [79594.95], 'U_CTIO': [3650], 'U_VIMOS': [3650], 'F435w': [4350], 'F606w': [6060], 'F775w': [7750], 'F814w': [8140], 'F850lp': [8500], 'F098m': [980], 'F105w': [1050], 'F125w': [125], 'F160w': [160]}

magerrdict = {'eI2_SEDS': 'I2_SEDS', 'eR_MUSYC': 'R_MUSYC', 'eV': 'V', 'eIA738_MUSYC': 'IA738_MUSYC', 'eZ_MUSYC': 'z_MUSYC', 'eg': 'g', 'k_er_au': 'k_au', 'eF105w': 'F105w', 'eIA527_MUSYC': 'IA527_MUSYC', 'eI3_SIMPLE': 'I3_SIMPLE', 'j_er_au': 'j_au', 'eIA598_MUSYC': 'IA598_MUSYC', 'eF775w': 'F775w', 'eF850lp': 'F850lp', 'eFUV': 'FUV', 'ej_TENIS': 'j_TENIS', 'eKs_ISAAC': 'ks_ISAAC', 'eI4_GOODS': 'I4_GOODS', 'er': 'r', 'eB_MUSYC': 'B_MUSYC', 'eF125w': 'F125w', 'eIA574_MUSYC': 'IA574_MUSYC', 'eI_MUSYC': 'I_MUSYC', 'eIA484_MUSYC': 'IA484_MUSYC', 'eI1_SEDS': 'I1_SEDS', 'eu': 'u', 'eIA679_MUSYC': 'IA679_MUSYC', 'eIA827_MUSYC': 'IA827_MUSYC', 'eU_VIMOS': 'U_VIMOS', 'eB': 'B', 'eNUV': 'NUV', 'eIA624_MUSYC': 'IA624_MUSYC', 'ez': 'z', 'eIA767_MUSYC': 'IA767_MUSYC', 'eH_MUSYC': 'H_MUSYC', 'eI1_SIMPLE': 'I1_SIMPLE', 'eIA464_MUSYC': 'IA464_MUSYC', 'eF435w': 'F435w', 'ev_TENIS': 'v_TENIS', 'eI1': 'I1', 'eIA651_MUSYC': 'IA651_MUSYC', 'eI3': 'I3', 'eI2': 'I2', 'eI4': 'I4', 'ez_TENIS': 'z_TENIS', 'eIA797_MUSYC': 'IA797_MUSYC', 'h_er_au': 'h_au', 'eK': 'K', 'eIA427_MUSYC': 'IA427_MUSYC', 'eIA505_MUSYC': 'IA505_MUSYC', 'ei': 'i', 'i_er_au': 'i_au', 'g_er_au': 'g_au', 'eJ_MUSYC': 'J_MUSYC', 'u_er_au': 'u_au', 'eKs_HAWKI': 'Ks_HAWKI', 'eks_TENIS': 'ks_TENIS', 'eF814w': 'F814w', 'eU_MUSYC': 'U_MUSYC', 'eF098m': 'F098m', 'eIA550_MUSYC': 'IA550_MUSYC', 'eU38_MUSYC': 'U38_MUSYC', 'eV_MUSYC': 'V_MUSYC', 'eI2_SIMPLE': 'I2_SIMPLE', 'eK_MUSYC': 'K_MUSYC', 'eIA709_MUSYC': 'IA709_MUSYC', 'eJ': 'J', 'r_er_au': 'r_au', 'eF606w': 'F606w', 'ei3_GOODS': 'I3_GOODS', 'eIA445_MUSYC': 'IA445_MUSYC', 'z_er_au': 'z_au', 'eF160w': 'F160w', 'eU_CTIO': 'U_CTIO', 'eIA856_MUSYC': 'IA856_MUSYC', 'eI4_SIMPLE': 'I4_SIMPLE', 'irac1_err_servs': 'irac1_servs', 'irac2_err_servs': 'irac2_servs'}

def SFR_f(M, z):
    t = age_at_z(z)
    return (0.8-0.022*t)*M - (6.09 - 0.07*t)

def find_nearest(array,value):
    return (np.abs(array-value)).argmin()

def readcat(cat):
    f = open(cat, 'r')
    entries = f.readlines()
    f.close()
    return entries

def data_at_z(z_lim, vvdslist, cosmoslist, ecdfslist):

    vvdsdata = readdata(vvdslist)
    cosmosdata = readdata(cosmoslist)
    ecdfsdata = readdata(ecdfslist)
    vvdsdata['specdir'], cosmosdata['specdir'], ecdfsdata['specdir'] = (vvdslist[:-4]+'/VVDS2h'), (cosmoslist[:-4]+'/COSMOS'), (ecdfslist[:-4]+'/ECDFS')
    
    mass, masslow, masshigh = (vvdsdata['mass']+cosmosdata['mass']+ecdfsdata['mass']), (vvdsdata['mass_low68']+cosmosdata['mass_low68']+ecdfsdata['mass_low68']), (vvdsdata['mass_high68']+cosmosdata['mass_high68']+ecdfsdata['mass_high68'])
    sfr, sfrlow, sfrhigh = (vvdsdata['sfr']+cosmosdata['sfr']+ecdfsdata['sfr']), (vvdsdata['sfr_low68']+cosmosdata['sfr_low68']+ecdfsdata['sfr_low68']), (vvdsdata['sfr_high68']+cosmosdata['sfr_high68']+ecdfsdata['sfr_high68'])
    ssfr, ssfrlow, ssfrhigh = (vvdsdata['ssfr']+cosmosdata['ssfr']+ecdfsdata['ssfr']), (vvdsdata['ssfr_low68']+cosmosdata['ssfr_low68']+ecdfsdata['ssfr_low68']), (vvdsdata['ssfr_high68']+cosmosdata['ssfr_high68']+ecdfsdata['ssfr_high68'])
    age, agelow, agehigh = (vvdsdata['Age']+cosmosdata['Age']+ecdfsdata['Age']), (vvdsdata['age_low68']+cosmosdata['age_low68']+ecdfsdata['age_low68']), (vvdsdata['age_high68']+cosmosdata['age_high68']+ecdfsdata['age_high68'])
    zs = (vvdsdata['z_spec']+cosmosdata['z_spec']+ecdfsdata['z_spec'])
    zphot = (vvdsdata['photo_z']+cosmosdata['photo_z']+ecdfsdata['zphot_MUSYC'])
    mags = (vvdsdata['mags']+cosmosdata['mags']+ecdfsdata['mags'])
    id = (vvdsdata['#ident']+cosmosdata['#ident']+ecdfsdata['#ident'])
    zf = (vvdsdata['zflags']+cosmosdata['zflags']+ecdfsdata['zflags'])
    model = (vvdsdata['Model']+cosmosdata['Model']+ecdfsdata['Model'])
    magi = (vvdsdata['magi']+cosmosdata['magi']+ecdfsdata['magi'])
    ebv = (vvdsdata['ebv']+cosmosdata['ebv']+ecdfsdata['ebv'])
    
    
    f = open('VVDS2h_obj.list', 'w')
    for ident in vvdsdata['#ident']:
        f.write(ident+'\n')
    f.close()

    f = open('COSMOS_obj.list', 'w')
    for ident in cosmosdata['#ident']:
        f.write(ident+'\n')
    f.close()

    f = open('ECDFS_obj.list', 'w')
    for ident in ecdfsdata['#ident']:
        f.write(ident+'\n')
    f.close()
    
    fn = []
    for k in range(len(vvdsdata['spec1d'])):
        fn.append(os.path.join(vvdsdata['specdir'],vvdsdata['spec1d'][k]))
    for k in range(len(cosmosdata['spec1d'])):
        fn.append(os.path.join(cosmosdata['specdir'],cosmosdata['spec1d'][k]))
    for k in range(len(ecdfsdata['spec1d'])):
        fn.append(os.path.join(ecdfsdata['specdir'],ecdfsdata['spec1d'][k]))
    
    imfile = []
    for k in range(len(fn)):
        temp = fn[k].replace('atm_clean', 'pix*2D').replace('Spec1D', 'spec2d')
        imfile.append(glob.glob(temp)[0])
    
    k = 0
    l = 0
    sfrplot = sfr
    massplot = mass

    while (k<len(mass)):
        if float(zs[k])>=z_lim:
            if float(mass[k])<0.0:
                massplot[k]=13.1
                sfrplot[k]=1.0+l*0.2
                l+=1
            mass[k]=float(mass[k])
            masslow[k]=mass[k]-float(masslow[k])
            masshigh[k]=float(masshigh[k])-mass[k]
            sfr[k]=float(sfr[k])
            sfrlow[k]=sfr[k]-float(sfrlow[k])
            sfrhigh[k]=float(sfrhigh[k])-sfr[k]
            ssfr[k]=float(ssfr[k])
            ssfrlow[k]=ssfr[k]-float(ssfrlow[k])
            ssfrhigh[k]=float(ssfrhigh[k])-ssfr[k]
            age[k]=float(age[k])*1e-9
            agelow[k]=age[k]-float(agelow[k])
            agehigh[k]=float(agehigh[k])-age[k]
            k+=1
        else:
            del id[k], fn[k], imfile[k], zs[k], zf[k], magi[k], mags[k], model[k], mass[k], masslow[k], masshigh[k], sfr[k], sfrlow[k], sfrhigh[k], ssfr[k], ssfrlow[k], ssfrhigh[k], age[k], agelow[k], agehigh[k]

    #print np.sort(age)


    return id, fn, imfile, zs, zphot, zf, magi, ebv, mags, model, mass, masslow, masshigh, sfr, sfrlow, sfrhigh, ssfr, ssfrlow, ssfrhigh, age, agelow, agehigh, massplot, sfrplot

def readdata(catname, showpar = False):
    cat = readcat(catname)

    x = cat[0].split()
    Ncol, parameters = len(x), {}
    
    allinfo = {}
    
    for i in range(len(x)):
        parameters[i]=x[i]
        if x[i]=='offset_servs':
            offk = i
            #print offk
        allinfo[x[i]]=[]
    allinfo['mags']=[]
    
    for i in range(1, len(cat)):
        x = cat[i].split()
        magnitudes, magerrors = {}, {}
        cwl = []
        for k in range(Ncol):
            
            if not cwldict.has_key(parameters[k]):
                if not magerrdict.has_key(parameters[k]):
                    allinfo[parameters[k]].append(x[k])
                else:
                    magerrors[magerrdict[parameters[k]]] = float(x[k])#np.sqrt(float(x[k])**2+(cwldict[magerrdict[parameters[k]]][1])**2)
            else:
                if parameters[k] in ['irac1_servs', 'irac2_servs']:
                    #print x[k]
                    if float(x[offk])>0 and float(x[k])>0: magnitudes[parameters[k]] = float(x[k])# + float(x[offk])
                    else: magnitudes[parameters[k]] = float(x[k])
                else: magnitudes[parameters[k]] = float(x[k])
        mags = []
        for key in magnitudes.keys():
            a = magerrors[key]
            mags.append([cwldict[key][0], magnitudes[key], magerrors[key], key])
        allinfo['mags'].append(mags)

    if showpar: print parameters

    return allinfo

def mstar(z):
    return -18.56 - 1.37*z + 0.18*z*z

def Ez(z):
    return 1.0/np.sqrt(sigmal+sigmam*((1+z)**3))

def dl(z):
    #luminosity distance at z in pc
    return (1+z)*13.39*quad(Ez, 0.0, z)[0]

def agesbb(x):
    return 13.39/(x*np.sqrt(sigmal+sigmam/(x*x*x)))

def age_at_z(z):
    return quad(agesbb, 0.0, 1.0/(1.0+z))[0]


def define_lines():
    alpha = r'$\alpha$'
    beta  = r'$\beta$'
    delta = r'$\delta$'
    gamma = r'$\gamma$'
    eps   = r'$\epsilon$'

    return collections.OrderedDict([
        ('Ly-limit',        ['Ly-limit', 911.74]),
        ('Ly-gamma',        ['Ly'+gamma, 949.73]),
        #('Ly-delta',        ['Ly'+delta, 972.53]),
        ('Ly-beta',        ['Ly'+beta ,1025.72]),
        #('OVI 1033',        ['OVI'     ,1033.83]),
        ('Si-II',      ['Si-IV',   1396.76]),
        ('Ly-alpha 1215',   ['Ly'+alpha,1215.67]),
        #('N-V 1240',        ['N-V',     1240.14]),
        ('OI 1304',         ['OI',      1304.35]),
        ('Si-IV 1396',      ['Si-IV',   1396.76]),
        ('C-IV 1549',       ['C-IV',    1549.06]),
        ('C-III 1908',      ['C-III',   1908.73]),
        ('Mg-II 2798',      ['Mg-II',   2798.75]),
        ('[OII] 3728',      ['[OII]',   3728.48]),
        ('[NeIII] 3869',    ['[NeIII]', 3869.85]),
        ('H-8 3890',        ['H-8',     3890.15]),
        ('H-eps 3971',      ['H'+eps,   3971.20]),
        ('H-delta 4102',    ['H'+delta, 4102.89]),
        ('H-gamma 4341',    ['H'+gamma, 4341.68]),
        ('H-beta 4862',     ['H'+beta,  4862.68]),
        ('[OIII] 4960',     ['[OIII]',  4960.30]),
        ('[OIII] 5008',     ['[OIII]',  5008.24]),
        ('[NI] 5199',       ['[NI]',    5199   ]),
        ('HeI 5877',        ['HeI',     5877.29]),
        ('[OI] 6302',       ['[OI]',    6302.05]),
        ('[NII] 6549',      ['[NII]',   6549.85]),
        ('H-alpha 6564',    ['H'+alpha, 6564.61]),
        ('[NII] 6585',      ['[NII]',   6585.28]),
        ('SII 6718',        ['SII',     6718.29]),
        ('SII 6732',        ['SII',     6732.67]),
        ('A:Ca(H) 3934',    ['A:Ca(H)', 3934.78]),
        ('A:Ca(K) 3969',    ['A:Ca(K)', 3969.59]),
        ('A:G-band 4300',   ['A:G-band',4300.4 ]),
        ('A:Mg-1 5167',     ['A:Mg-1',  5167.3222]),
        ('A:Mg-2 5172',     ['A:Mg-2',  5172.6847]),
        ('A:Mg-3 5183',     ['A:Mg-3',  5183.6046]),
        ('A:Na 5894',       ['A:Na',    5894.57  ]),
        ('A:H-delta 4102',  ['A:H'+delta,4102.89]),
        ('A:H-gamma 4341',  ['A:H'+gamma,4341.68]),
        ('A:H-beta 4862',   ['A:H'+beta, 4862.68]),
        ('A:H-alpha 6564',  ['A:H'+alpha,6564.61])
    ])

def getspectrum(filename):
    hdu = fits.open(filename)
    y = hdu[0].data
    if len(y)==1:
        y = hdu[0].data[0]

    header = hdu[0].header
    crval = header['CRVAL1']
    cdelt = header['CDELT1']
    x=[]
    for i in range(len(y)):
        x.append(crval+i*cdelt)
    y = y.tolist()
    return x, y, cdelt

def flux(mag, cw = 7.6905E+03, c = 1.0):
    return c*c_AAs*np.power(10,(-0.4*(mag+48.6)))/cw/cw

def fnu(mag):
    return np.power(10,(-0.4*(mag+48.6)))

def fluxi(mag):
    cw = 7.6905E+03
    return c_AAs*np.power(10,(-0.4*(mag+48.6)))/cw/cw

def mAB(wl, flux):
    if flux>0: return -2.5*np.log10(wl*wl*flux/c_AAs)-48.5
    else: return None

def mAB_fnu(flux):
    if flux>0: return -2.5*np.log10(flux)-48.5
    else: return None

def specmag(x, y, filterpath = './filters/ip.pb', center = 7.6905E+03):
    xrf, yrf = [], []
    xrf, yrf = np.loadtxt(filterpath, unpack=True)

    filt_int  = np.interp(x,xrf,yrf)

    I1        = simps(y*filt_int*x,x)                     #Denominator
    I2        = simps(  filt_int/x,x)                     #Numerator
    fnu       = I1/I2 / c_AAs                          #Average flux density
    #fl = I1/I2/center/center #flambda

    return mAB_fnu(fnu)

def specflux(x, y, filterpath = './filters/ip.pb', center = 7.6905E+03):
    xrf, yrf = [], []
    xrf, yrf = np.loadtxt(filterpath, unpack=True)

    filt_int  = np.interp(x,xrf,yrf)

    I1        = simps(y*filt_int*x,x)                     #Denominator
    I2        = simps(  filt_int/x,x)                     #Numerator

    fnu       = I1/I2 / c_AAs                          #Average flux density
    return fnu

def specfnu(x, y, filterpath = './filters/ip.pb', center = 7.6905E+03):
    #this is all wrong!!!!!
    xrf, yrf = [], []
    xrf, yrf = np.loadtxt(filterpath, unpack=True)

    dnu = c_AAs / x

    filt_int  = np.interp(x,xrf,yrf)

    fnu = simps((y*filt_int),dnu)


    return fnu

def ew(pos1, pos2, wl, fl):
    #print fl
    zoom = [k for k in range(len(wl)) if (wl[k]>pos1[0] and wl[k]<pos2[0])]
    x = wl[zoom[0]:zoom[-1]]
    #x = np.arange(pos1[0], pos2[0], 1.0)
    f = np.interp(x, wl, fl)
    c = np.interp(x, [pos1[0], pos2[0]], [pos1[1], pos2[1]])
    s = 1.0 - np.asarray(f)/c
    return np.sum(s), wl, fl

def delete_ticks(ax):
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    #ax.axis('off')

#*****************************CREATE_TEMPLATES***************************************************
def create_templates():
    
    f = open('BC03COMB_MOD.list', 'r')
    modelfiles = f.readlines()
    f.close()
    
    modeldict = { 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: []}
    
    for k in range(len(model)):
        if int(model[k])>0 and float(age[k])>0: modeldict[int(model[k])].append(str(age[k]))

    for k in range(len(modeldict)):
        if len(modeldict[k+1])>0:
            #print k+1
            tempname = '/Users/yanahusanova/progs/phd/templates/template'+str(k+1)+'nu.dat'
            f = open('create_template.sh', 'w')
            galages = str(modeldict[k+1])
            galages = galages.translate(None, " ' [ ]")
            #print k+1, galages
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
echo "%s"
echo "-100, 15000"
echo "%s"
./galaxevpl /Users/yanahusanova/progs/phd/models/%s
cd ../../progs/phd/""" % (galages, tempname, modelfiles[k].split()[0])
            f.write(s.encode('utf8'))
            f.close()
            #call(['./create_template.sh'])
            tempname = '/Users/yanahusanova/progs/phd/templates/template'+str(k+1)+'l.dat'
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
echo "%s"
echo "100, 15000"
echo "%s"
./galaxevpl ./%s
cd ../../progs/phd/""" % (galages, tempname, modelfiles[k].split()[0])
            f.write(s.encode('utf8'))
            f.close()
    #call(['./create_template.sh'])
    return modeldict


def blabla():
    
    f = open('temp.txt', 'w')

    for i in range(1, len(cat)):
        
        x = cat[i].split()
        
        id = int(x[0])

        s = r'''%d
''' % (id)

        f.write(s.encode('utf8'))

    f.close()

    flags = []
    
    for i in range(1, len(cat)):
        
        x = cat[i].split()

        if not int(x[7]) in flags:
            flags.append(int(x[7]))

    print sorted(flags)

    inputdir = os.path.join(dir,'test3')
    fn = os.path.join(inputdir,'iselspecallcosmos.txt')
    cat = readcat(fn)
    
    f = open(os.path.join(inputdir,'izselcosmos.txt'), 'w')
    f0 = open('temp.txt', 'w')
    
    for i in range(1, len(cat)):
        
        x = cat[i].split()

        flag = int(x[7])

        if int(x[7])>0 and int(x[7])!=20:
            z = float(x[6])
            if not (int(x[7]) in [1, 11, 21, 31, 41, 241]):
                if float(x[12])>-21.74 and float(x[12])<-20.25 and z >= 3.8 and z<=5.0:
                    s = r'''%s
''' % x[0]
                    f0.write(s.encode('utf8'))
        else:
            z = float(x[9])

        if z>=3.0 and z<=5.0:
            s = cat[i]
            f.write(s.encode('utf8'))

    f.close()
    f0.close()

    dir = "./templates"
    fn = os.path.join(dir,'1.txt')
    cat = readcat(fn)
    
    f = open(os.path.join(dir,'template13.dat'), 'w')
    
    for i in range(len(cat)):
        
        x = cat[i].split()
        
        wl = float(x[0])
        fl = float(x[1])
        
        if wl<10000:
            fl*=0.13
        
        fl*=10e-18
            
        s = r'''%s  %.4e
''' % (x[0], fl)
        f.write(s.encode('utf8'))

    f.close()
