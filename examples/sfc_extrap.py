import numpy as np
#***********************************************************************
#          Shift variables from HT = HQ to wind height HU = 10m
#          compute potential tempertature at HU,  air density
#                                 Bill Large Feb 1996
#          Note: it is assumed that HT=HQ; and only HT is used in this
#          routine.  Pressure is now in Pa (Jan Morzel, Feb 1997)
#          and u-,v-winds are inputs instead of wind speed.
#
#
#          Assume:
#           1) A NEUTRAL 10M DRAG COEFFICIENT = CDN = 
#                    .0027/U10 + .000142 + .0000764 U10
#           2) A NEUTRAL 10M STANTON NUMBER  CTN= .0327 SQRT(CDN), UNSTABLE
#                                         CTN= .0180 SQRT(CDN), STABLE
#           3) A NEUTRAL 10M DALTON  NUMBER  CEN= .0346 SQRT(CDN)
#           4) THE SATURATION HUMIDITY OF AIR AT T(K)  = QSAT(T)  (kg/kg)
#  NOTE  1) HERE, TSTAR = <WT>/U*, AND QSTAR = <WQ>/
#        2) WIND SPEEDS SHOULD ALL BE ABOVE A MINIMUM SPEED SAY 0.5 m/
#        3) WITH OPTIONAL ITERATION LOOP , NITER=3, SHOULD 
# ****  THIS VERSION IS FOR ANALYSES INPUTS WITH HU = 10m and HT = HQ *****
#      QSAT(TKELVIN) = 640380. / EXP(5107.4/TKELVIN)          ! kg / m3                                                                                                                              
#      CDN(Umps) = 0.0027 / Umps + .000142 + .0000764 * Umps                   
#***********************************************************************

ht = 2.0
z_ref = 10.0 # reference height (m)
f1 = 0.608
c1 = 1.0
c2 = 2.0
umin = 0.5 # minimum wind speed
zolmin  = -100. # minimum stability parameter
vonk = 0.4 # von Karman's constant
gamma = 0.0098  # adiabatic lapse rate (deg/m)
grav = 9.80616  # gravity (m/s2)


czol = z_ref*vonk*grav

def qsat(tk):

    qsat = 640380.0/np.exp(5107.4/tk)

    return qsat

def cdn(u):

    cdn = 0.0027/u + 0.000142 + 0.0000764*u

    return cdn
    
def extrap_sfc(uref,vref,tsfc,qsfc,tskin,bp):

    aln = np.log(ht/z_ref)
    sh = uref*uref + vref*vref
    sh = np.maximum(np.sqrt(sh),umin)
    t0 = tsfc*(c1 + f1*qsfc)  # initial virtual temperature guess, z/L=0;HT=HU-Hq=Z
    rho_air = 1.29*(273.16/t0)*(bp/101000.) # ambient air density (kg/m3)
    ssq = 0.98*(qsat(tskin)/rho_air) # sea surface humidity (kg/kg)
    delq = qsfc - ssq
    delp = tsfc + gamma*ht - tskin
    stable = 0.5 + 0.5*(np.sign(delp))
    rdn = np.sqrt(cdn(sh))
    rhn = (c1-stable)*0.0327 + stable*0.0180
    ren = 0.0346
    ustar = rdn*sh
    tstar = rhn*delp
    qstar = ren*delq

    for i in np.arange(0,5):
        huol = czol * (tstar/t0 + qstar/(c1/f1 + qsfc)) / ustar**2.0
        huol = np.maximum(huol,zolmin)
        stable = 0.5+ 0.5*np.sign(huol)
        htol = huol * ht/z_ref
        xsq = np.maximum(np.sqrt(np.abs(c1-16.*huol)),c1)
        x = np.sqrt(xsq)
        psimh = -5.0 * huol*stable + (c1-stable)*np.log((c1+x*(c2+x))*(c1+x*x)/8.0)-c2*np.arctan(x) + 1.571
        psixz = -5.0 * huol*stable+(c1-stable)*c2*np.log((c1+xsq)/c2)
        xsq = np.maximum(np.sqrt(np.abs(c1-16.0*htol)),c1)
        psixh = -5.0 * htol * stable + (c1-stable)*c2*np.log((c1+xsq)/c2)
        rd = rdn / (c1-rdn/vonk*psimh)
        uzn = np.maximum(sh*rd/rdn,umin)
        rdn=np.sqrt(cdn(uzn))
        ren = 0.0346
        rhn = (c1-stable)*0.0327 + stable*0.0180
        rd = rdn / (c1-rdn/vonk*psimh)
        rh = rhn / (c1+rhn/vonk*(aln - psixh))
        re = ren / (c1+ren/vonk*(aln - psixh))
        ustar = rd*sh
        qstar = re*delq
        tstar = rh*delp

    t10 = delp + tskin + tstar/vonk*(psixh - psixz - aln)
    q10 = qsfc + qstar/vonk * (psixh - psixz - aln)
    
    return t10,q10

    
    
    

  
