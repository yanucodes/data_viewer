import os
import glob
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt
from matplotlib._png import read_png
from astropy.io import fits
import collections
from scipy.integrate import simps
from scipy.integrate import quad
c_AAs     = 2.99792458e18
#from matplotlib.offsetbox import OffsetImage, AnnotationBbox

cwldict = {'u_au': [3823.29], 'g_au': [4720.0], 'r_au': [6404.0], 'i_au': [7764.0], 'z_au': [9256.0], 'j_au': [12520.0], 'h_au': [16310.0], 'k_au': [21460.0], 'FUV': [1539.78], 'NUV': [2313.89], 'u': [3823.29], 'B': [4458.34], 'V': [5477.83], 'g': [4777.07], 'r': [6288.71], 'i': [7683.88], 'z': [9036.88], 'J': [12534.65], 'K': [21539.88], 'I1': [35634.28], 'I2': [45110.13], 'I3': [57593.39], 'I4': [79594.95], 'U_MUSYC': [3650], 'B_MUSYC': [4450], 'V_MUSYC': [5510], 'R_MUSYC': [6580], 'I_MUSYC': [8060], 'z_MUSYC': [9000], 'J_MUSYC': [12200], 'H_MUSYC': [16300], 'K_MUSYC': [21900], 'IA427_MUSYC': [4263.45], 'IA464_MUSYC': [4635.13], 'IA445_MUSYC': [4450], 'IA484_MUSYC': [4849.2], 'IA505_MUSYC': [5050], 'IA527_MUSYC': [5270], 'IA550_MUSYC': [5500], 'IA574_MUSYC': [5740], 'IA598_MUSYC': [5980], 'IA624_MUSYC': [6240], 'IA651_MUSYC': [6510], 'IA679_MUSYC': [6790], 'IA709_MUSYC': [7090], 'IA738_MUSYC': [7380], 'IA767_MUSYC': [7970], 'IA797_MUSYC': [7970], 'IA827_MUSYC': [8270], 'IA856_MUSYC': [8560]}

magerrdict = {'eI2_SEDS': 'I2_SEDS', 'eR_MUSYC': 'R_MUSYC', 'eV': 'V', 'eIA738_MUSYC': 'IA738_MUSYC', 'eZ_MUSYC': 'z_MUSYC', 'eg': 'g', 'k_er_au': 'k_au', 'eF105w': 'F105w', 'eIA527_MUSYC': 'IA527_MUSYC', 'eI3_SIMPLE': 'I3_SIMPLE', 'j_er_au': 'j_au', 'eIA598_MUSYC': 'IA598_MUSYC', 'eF775w': 'F775w', 'eF850lp': 'F850lp', 'eFUV': 'FUV', 'ej_TENIS': 'j_TENIS', 'eKs_ISAAC': 'ks_ISAAC', 'eI4_GOODS': 'I4_GOODS', 'er': 'r', 'eB_MUSYC': 'B_MUSYC', 'eF125w': 'F125w', 'eIA574_MUSYC': 'IA574_MUSYC', 'eI_MUSYC': 'I_MUSYC', 'eIA484_MUSYC': 'IA484_MUSYC', 'eI1_SEDS': 'I1_SEDS', 'eu': 'u', 'eIA679_MUSYC': 'IA679_MUSYC', 'eIA827_MUSYC': 'IA827_MUSYC', 'eU_VIMOS': 'U_VIMOS', 'eB': 'B', 'eNUV': 'NUV', 'eIA624_MUSYC': 'IA624_MUSYC', 'ez': 'z', 'eIA767_MUSYC': 'IA767_MUSYC', 'eH_MUSYC': 'H_MUSYC', 'eI1_SIMPLE': 'I1_SIMPLE', 'eIA464_MUSYC': 'IA464_MUSYC', 'eF435w': 'F435w', 'ev_TENIS': 'v_TENIS', 'eI1': 'I1', 'eIA651_MUSYC': 'IA651_MUSYC', 'eI3': 'I3', 'eI2': 'I2', 'eI4': 'I4', 'ez_TENIS': 'z_TENIS', 'eIA797_MUSYC': 'IA797_MUSYC', 'h_er_au': 'h_au', 'eK': 'K', 'eIA427_MUSYC': 'IA427_MUSYC', 'eIA505_MUSYC': 'IA505_MUSYC', 'ei': 'i', 'i_er_au': 'i_au', 'g_er_au': 'g_au', 'eJ_MUSYC': 'J_MUSYC', 'u_er_au': 'u_au', 'eKs_HAWKI': 'Ks_HAWKI', 'eks_TENIS': 'ks_TENIS', 'eF814w': 'F814w', 'eU_MUSYC': 'U_MUSYC', 'eF098m': 'F098m', 'eIA550_MUSYC': 'IA550_MUSYC', 'eU38_MUSYC': 'U38_MUSYC', 'eV_MUSYC': 'V_MUSYC', 'eI2_SIMPLE': 'I2_SIMPLE', 'eK_MUSYC': 'K_MUSYC', 'eIA709_MUSYC': 'IA709_MUSYC', 'eJ': 'J', 'r_er_au': 'r_au', 'eF606w': 'F606w', 'ei3_GOODS': 'I3_GOODS', 'eIA445_MUSYC': 'IA445_MUSYC', 'z_er_au': 'z_au', 'eF160w': 'F160w', 'eU_CTIO': 'U_CTIO', 'eIA856_MUSYC': 'IA856_MUSYC', 'eI4_SIMPLE': 'I4_SIMPLE'}

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

def data_at_z(z_lim):

    vvdsdata = readdata('./data/VVDS2h.txt')
    cosmosdata = readdata('./data/cosmos.txt')
    ecdfsdata = readdata('./data/ECDF-2017-01-16.txt')
    vvdsdata['specdir'], cosmosdata['specdir'], ecdfsdata['specdir'] = './data/VVDS2h', './data/COSMOS', './data/ECDFS'
    
    mass, masslow, masshigh = (vvdsdata['mass']+cosmosdata['mass']+ecdfsdata['mass']), (vvdsdata['mass_low68']+cosmosdata['mass_low68']+ecdfsdata['mass_low68']), (vvdsdata['mass_high68']+cosmosdata['mass_high68']+ecdfsdata['mass_high68'])
    sfr, sfrlow, sfrhigh = (vvdsdata['sfr']+cosmosdata['sfr']+ecdfsdata['sfr']), (vvdsdata['sfr_low68']+cosmosdata['sfr_low68']+ecdfsdata['sfr_low68']), (vvdsdata['sfr_high68']+cosmosdata['sfr_high68']+ecdfsdata['sfr_high68'])
    ssfr, ssfrlow, ssfrhigh = (vvdsdata['ssfr']+cosmosdata['ssfr']+ecdfsdata['ssfr']), (vvdsdata['ssfr_low68']+cosmosdata['ssfr_low68']+ecdfsdata['ssfr_low68']), (vvdsdata['ssfr_high68']+cosmosdata['ssfr_high68']+ecdfsdata['ssfr_high68'])
    age, agelow, agehigh = (vvdsdata['Age']+cosmosdata['Age']+ecdfsdata['Age']), (vvdsdata['age_low68']+cosmosdata['age_low68']+ecdfsdata['age_low68']), (vvdsdata['age_high68']+cosmosdata['age_high68']+ecdfsdata['age_high68'])
    zs = (vvdsdata['z_spec']+cosmosdata['z_spec']+ecdfsdata['z_spec'])
    mags = (vvdsdata['mags']+cosmosdata['mags']+ecdfsdata['mags'])
    id = (vvdsdata['#ident']+cosmosdata['#ident']+ecdfsdata['#ident'])
    zf = (vvdsdata['zflags']+cosmosdata['zflags']+ecdfsdata['zflags'])
    model = (vvdsdata['Model']+cosmosdata['Model']+ecdfsdata['Model'])
    magi = (vvdsdata['magi']+cosmosdata['magi']+ecdfsdata['magi'])
    
    fn = []
    for k in range(len(vvdsdata['spec1d'])):
        fn.append(os.path.join(vvdsdata['specdir'],vvdsdata['spec1d'][k]))
    for k in range(len(cosmosdata['spec1d'])):
        fn.append(os.path.join(cosmosdata['specdir'],cosmosdata['spec1d'][k]))
    for k in range(len(ecdfsdata['spec1d'])):
        fn.append(os.path.join(ecdfsdata['specdir'],ecdfsdata['spec1d'][k]))
    
    imfile = []
    for k in range(len(fn)):
        temp = fn[k].replace('atm_clean', 'pix*2D').replace('pec1', 'pec2')
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

    return id, fn, imfile, zs, zf, magi, mags, model, mass, masslow, masshigh, sfr, sfrlow, sfrhigh, ssfr, ssfrlow, ssfrhigh, age, agelow, agehigh, massplot, sfrplot

def readdata(catname):
    cat = readcat(catname)

    x = cat[0].split()
    Ncol, parameters = len(x), {}
    
    allinfo = {}
    
    for i in range(len(x)):
        parameters[i]=x[i]
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
                    magerrors[magerrdict[parameters[k]]] = float(x[k])
            else:
                magnitudes[parameters[k]] = float(x[k])
        mags = []
        for key in magnitudes.keys():
            mags.append([cwldict[key][0], magnitudes[key], magerrors[key], key])
        allinfo['mags'].append(mags)

    return allinfo

def mstar(z):
    return -18.56 - 1.37*z + 0.18*z*z

def agesbb(x):
    return 13.39/(x*np.sqrt(0.76+0.24/(x*x*x)))

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

def flux(mag, cw, c = 1.0):
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

def specmag(x, y, filterpath = './filters/i_SDSS.res', center = 7683.88):
    xrf, yrf = [], []
    xrf, yrf = np.loadtxt(filterpath, unpack=True)

    filt_int  = np.interp(x,xrf,yrf)

    I1        = simps(y*filt_int*x,x)                     #Denominator
    I2        = simps(  filt_int/x,x)                     #Numerator
    fnu       = I1/I2 / c_AAs                          #Average flux density
    #fl = I1/I2/center/center #flambda

    return mAB_fnu(fnu)

def delete_ticks(ax):
    ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    ax.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    #ax.axis('off')
