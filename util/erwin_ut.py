class Halo() :
    def __init__(self,id,aexp) :
        self.id = id
        self.aexp = aexp
    def __repr__(self):
        return 'Halo(id=%d, aexp=%6.4f)' % (self.id, self.aexp)
    def getId(self) :
        return self.id
    def getAexpn(self) :
        return self.aexp

def integrate_profile ( pro, rbins, rcore=0.0 ) :
    from numpy import zeros
    num_bins = len(pro)
    sum_profile = zeros(num_bins)
    sum = 0.0
    for i in range(num_bins) :
        if ( rbins[i] >= rcore ):
            sum += pro[i]
            sum_profile[i] = sum

    return sum_profile


def average_integrate_profile ( pro, rbins, weight, rcore=0.0 ) :
    from numpy import zeros
    num_bins = len(pro)
    sum_profile = zeros(num_bins)
    sum = 0.0
    sum_weight = 0.0
    for i in range(num_bins) :
        if ( rbins[i] >= rcore ):
            sum += pro[i]*weight[i]
            sum_weight += weight[i]
            sum_profile[i] = sum/sum_weight

    return sum_profile


def find_profile_max_smooth ( pro, rbins, s=0, delta=0.1) :

    from numpy import log10, amin, amax, argmin, argmax, linspace
    from scipy.interpolate import UnivariateSpline
    minx = amin(rbins)
    maxx = amax(rbins)
    x = linspace(minx,maxx, 1000)

    yspline = UnivariateSpline(rbins, pro, s = s )
    ys = yspline(x)

    p = amax(ys)
    ir = argmax(ys)
    r = x[ir]
    i = argmax(pro)
    #rmin, rmax  = (rbins[i]-rbins[i-1])/2.0,(rbins[i+1]-rbins[i])/2.0

    return (p, r)

def find_profile_min_smooth ( pro, rbins, s=0, delta=0.1 ) :

    from numpy import log10, amin, amax, argmin, argmax, linspace
    from scipy.interpolate import UnivariateSpline, InterpolatedUnivariateSpline
    from scipy.signal import savgol_filter
    lrbins = log10(rbins)
    minx = amin(lrbins)
    maxx = amax(lrbins)
    x = linspace(minx,maxx, 1000)
    yspline = UnivariateSpline(lrbins, pro)
    yspline.set_smoothing_factor(0.1)
    ys = yspline(x)

    maxtab, mintab = peakdet(ys, delta, x)
    imin =  argmin(mintab[:,1])
    rmin = 10**(mintab[:,0][imin])
    pmin = mintab[:,1][imin]

    return (pmin, rmin)


def find_profile_max ( pro, rbins, delta=0.1 ) :
 
    from numpy import log10, argmax
    x = log10(rbins)
    y = pro

    maxtab, mintab = peakdet(y, delta, x)

    imax =  argmax(maxtab[:,1])
    rmax = 10**(maxtab[:,0][imax])
    pmax = maxtab[:,1][imax]

    return (pmax, rmax)


def find_profile_min ( pro, rbins, delta=0.1 ) :
 
    from numpy import log10, argmin
    x = log10(rbins)
    y = pro

    maxtab, mintab = peakdet(y, delta, x)   
    #print mintab[:,1]
    imin =  argmin(mintab[:,1])
    rmin = 10**(mintab[:,0][imin])
    pmin = mintab[:,1][imin]

    return (pmin, rmin)

def compute_density_profile ( rmid, rr, mass ) :
    '''
    Compute density profile (Msun h^2/kpc^3) for a given mass profile in Msun/h. Takes in radial bin in kpc/h (rmid,rr)
    '''
    from numpy import array, pi
    num_bins = len(rmid)

    rmid0 = rmid[0]
    rr0 = rr[0]
    mass0 = mass[0]
    vol = 4./3.*pi*rr0**3
    lbinr = [rmid0]
    lrho = [mass0/vol]

    for i in xrange(1,num_bins) :
        # vol in (kpc/h)^3
        vol = 4./3.*pi*(rr[i]**3-rr0**3) 
			
        lbinr.append( rmid[i] )	
        lrho.append( (mass[i]-mass0)/vol )
		
        mass0 = mass[i]
        rr0 = rr[i]

    rho = array( lrho )
	
    return ( rho )

def compute_density_profiles ( rmid, rr, Mdm, Mgas, Mtot ) :
    '''
    Compute density profile (Msun h^2/kpc^3) for DM, gas and total given their mass profiles in Msun/h. Takes in radial bin in kpc/h (rmid,rr)
    '''
    from numpy import array, pi
    num_bins = len(rmid)

    rmid0 = rmid[0]
    rr0 = rr[0]
    Mdm0 = Mdm[0]
    Mgas0 = Mgas[0]
    Mtot0 = Mtot[0]
    vol = 4./3.*pi*rr0**3
    lbinr = [rmid0]
    lrhodm = [Mdm0/vol]
    lrhoicm = [Mgas0/vol]
    lrhotot = [Mtot0/vol]

    for i in range(1,num_bins) :
        # vol in (kpc/h)^3
        vol = 4./3.*pi*(rr[i]**3-rr0**3) 
			
        lbinr.append( rmid[i] )	
        lrhodm.append( (Mdm[i]-Mdm0)/vol )
        lrhoicm.append( (Mgas[i]-Mgas0)/vol )
        lrhotot.append( (Mtot[i]-Mtot0)/vol )
		
        rr0 = rr[i]
        Mdm0 = Mdm[i]
        Mgas0 = Mgas[i]
        Mtot0 = Mtot[i]

    rhodm = array(lrhodm)
    rhoicm = array(lrhoicm)
    rhotot = array(lrhotot)
	
    return ( rhodm, rhoicm, rhotot )

def compute_diffmass_profiles ( rmid, rr, Mdm, Mgas, Mtot ) :
    '''
    Compute differential mass profile (Msun/h) for DM, gas and total given their mass profiles in Msun/h. Takes in radial bin in kpc/h (rmid,rr)
    '''
    from numpy import array, pi
    num_bins = len(rmid)

    rmid0 = rmid[0]
    rr0 = rr[0]
    Mdm0 = Mdm[0]
    Mgas0 = Mgas[0]
    Mtot0 = Mtot[0]
    #vol = 4./3.*pi*rr0**3
    lbinr = [rmid0]
    lrhodm = [Mdm0]
    lrhoicm = [Mgas0]
    lrhotot = [Mtot0]

    for i in range(1,num_bins) :
        # vol in (kpc/h)^3
        #vol = 4./3.*pi*(rr[i]**3-rr0**3) 
			
        lbinr.append( rmid[i] )	
        lrhodm.append( (Mdm[i]-Mdm0) )
        lrhoicm.append( (Mgas[i]-Mgas0) )
        lrhotot.append( (Mtot[i]-Mtot0) )
		
        rr0 = rr[i]
        Mdm0 = Mdm[i]
        Mgas0 = Mgas[i]
        Mtot0 = Mtot[i]

    rhodm = array(lrhodm)
    rhoicm = array(lrhoicm)
    rhotot = array(lrhotot)
	
    return ( rhodm, rhoicm, rhotot )
        
def calc_log_slope (x,y, m = 5) :
    '''
    Compute log slope. Assumes that radius are spaced uniformly in log. 
    Uses SG filter with default window size m = 5; polynomial or order 2
    '''

    #from sg_filter import *
    import scipy.signal
    import numpy as np

    num_bins = len(x)

    dlogy  = scipy.signal.savgol_filter( np.log10(y), m, 2, deriv=1 )
    dlogx  = scipy.signal.savgol_filter( np.log10(x), m, 2, deriv=1 )
    #dlogx = np.zeros(num_bins)
    #for i in np.arange(len(dlogx)) :
    #    dlogx[i] = np.log10(x[2]) - np.log10(x[1])

    dlydlx = dlogy/dlogx

    return ( dlydlx )

 
def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %        PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """

    from numpy import NaN, Inf, arange, isscalar, asarray, array
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
 
    return array(maxtab), array(mintab)
 

def find_overdensity(rhalo, rr, mtot, aexpn, crit=True, omega_m=0.27, omega_l=0.73 ) :
    '''
    find overdensity for given halo radius and mass profile
    '''
    from numpy import sqrt, array, pi, log10

    Msun = 1.989e33
    kpc = 3.08567758e21
    H_0 = 100.0*1.e5/(1.e3*kpc)
    grav_c = 6.673e-8

    rho_crit = (3.0/(8.0*pi*grav_c))*H_0*H_0 / Msun * (kpc**3.0) # h^2Msun/kpc^3

    Ez = sqrt(omega_m/aexpn**3.0+omega_l)
    z = 1.0/aexpn-1.0

    if crit :
        rho_background = rho_crit*Ez**2 # h^2Msun/kpc^3

    else:
        # this should be the mean density, not critical
        rho_background = omega_m* rho_crit /aexpn**3.0 # h^2Msun/kpc^3

    # this is the average density enclosed in radius rr in kpc/h
    sphvol = 4.0*pi/3.0*rr**3.0
    rho_enc = log10(mtot/sphvol)
    lrr = log10(rr)

    i = 0
    while log10(rhalo) > lrr[i] :
        i += 1

    # interpolate to find rhovir
    # rvir = 10.0**(lrr[i-1] + (lrr[i]-lrr[i-1])/(rho_enc[i]-rho_enc[i-1]) * (log10(rho_vir) - rho_enc[i-1]))
    rho_vir = 10.0**(rho_enc[i-1] + (rho_enc[i]-rho_enc[i-1])/(lrr[i]-lrr[i-1]) * (log10(rhalo)-lrr[i-1]))
    overdensity = rho_vir/rho_background

    return overdensity, rho_vir

def find_rvir_mvir(rr, mtot, aexpn, overdensity = 200.0, crit = True, vir=False, omega_m=0.27, omega_l=0.73):
    '''
    find rvir, mvir, and rho_vir given mass profile mtot, radial bins rr, expansion factor and overdensity with respect to critical (default) or mean (critical = False). Radial bins should be in [kpc/h] and mass bins in [Msun/h]. 
    '''
    from numpy import sqrt, array, pi, log10

    Msun = 1.989e33
    kpc = 3.08567758e21
    H_0 = 100.0*1.e5/(1.e3*kpc)
    grav_c = 6.673e-8

    rho_crit = (3.0/(8.0*pi*grav_c))*H_0*H_0 / Msun * (kpc**3.0) # h^2Msun/kpc^3

    Ez = sqrt(omega_m/aexpn**3.0+omega_l)
    z = 1.0/aexpn-1.0

    if crit or vir:
        rho_background = rho_crit*Ez**2 # h^2Msun/kpc^3

    else:
        # this should be the mean density, not critical
        rho_background = omega_m* rho_crit /aexpn**3.0 # h^2Msun/kpc^3

    if vir :
        ## Use Bryan & Norman (1998) formula for virial overdensty wrt critical
        omega_z = omega_m/aexpn**3./Ez**2.
        overdensity = 18.0*pi*pi+82.0*(omega_z-1)-39.0*(omega_z-1)*(omega_z-1)

    rho_vir = overdensity * rho_background

    # this is the average density enclosed in radius rr in kpc/h
    sphvol = 4.0*pi/3.0*rr**3.0
    rho_enc = log10(mtot/sphvol)
    lrr = log10(rr)


    i = 0
    while rho_enc[i] > log10(rho_vir) :
        i += 1
        if i >= len(rho_enc) :
            break

    # interpolate to find rvir
    if i <  len(rho_enc) :
        rvir = 10.0**(lrr[i-1] + (lrr[i]-lrr[i-1])/(rho_enc[i]-rho_enc[i-1]) * (log10(rho_vir) - rho_enc[i-1]))
    else :
        rvir = 10.0**(lrr[-1])

    mvir =  4.0*pi/3.0*rvir**3.0*rho_vir

    return rvir, mvir, rho_vir

def getSelfSimilarValues (mvir, delta, crit=True, aexp=1.0, omega_m=0.27, omega_l=0.73,	omega_b = 0.0469, hubble=0.7 ) :
    ''' 
    Return self-similar quantites (e.g. T500c) for temperature [keV] pressure [erg/cm^3], and entropy [keV cm^2] given halo mass mvir [Msun/h], its overdensity delta, with respect to critical (crit=True) or mean (crit=False) at expansion factor aexp (= 1.0 by default). Assume WMAP5 cosmological parameters (true for Bolshoi and L500 simulations).
    '''
    from numpy import sqrt

    fb = omega_b/omega_m
    Ez = sqrt(omega_m/aexp**3.0+omega_l)
    mu = 0.59
    mue = 1.14

    mpivot = 1e15 # in Msun/h

    delta = float(delta)
	
    if crit :
        Tfit_norm = 11.055*(mu/0.59)*(hubble/0.7)**(2./3.)
        Tfit_alpha = 2./3.
        Tfit_beta = 2./3.
        Pfit_norm = 1.4458e-11*(fb/0.1737)*(hubble/0.7)**2
        Pfit_alpha = 2./3.
        Pfit_beta = 8./3.
        Kfit_alpha = 2./3.
        Kfit_beta = -2./3.
        #Kfit_norm = 1963.6*(mu/0.59)*(mue/1.14)**(2./3.)*(fb/0.1737)**(-2./3.)*(hubble/0.7)**(-2./3.)
        Kfit_norm = 1265.7*(mu/0.59)**(5./3.)*(fb/0.1737)**(-2./3.)*(hubble/0.7)**(-2./3.)

        Tvir = (delta/500.0)**(1./3.)*Tfit_norm*(mvir/mpivot)**(Tfit_alpha)*Ez**(Tfit_beta)
        Kvir = (delta/500.0)**(-1./3.)*Kfit_norm*(mvir/mpivot)**(Kfit_alpha)*Ez**(Kfit_beta) 
        Pvir = (delta/500.0)**(4./3.)*Pfit_norm*(mvir/mpivot)**(Pfit_alpha)*Ez**(Pfit_beta)
    else :
        Tfit_norm = 5.2647*(mu/0.59)*(hubble/0.7)**(2./3.)*(omega_m/0.27)**(1./3.)
        Tfit_alpha = 2./3.
        Tfit_beta = -1.0
        Pfit_norm = 7.4359e-13*(fb/0.1737)*(hubble/0.7)**2*(omega_m/0.27)**(4./3.)
        Pfit_alpha = 2./3.
        Pfit_beta = 4.0
        Kfit_alpha = 2./3.
        Kfit_beta = 1.0
        #Kfit_norm = 4123.3*(mu/0.59)*(mue/1.14)**(2./3.)*(fb/0.1737)**(-2./3.)*(hubble/0.7)**(-2./3.)*(omega_m/0.27)**(-1./3.)
        Kfit_norm = 2799.75*(mu/0.59)**(5./3.)*(fb/0.1737)**(-2./3.)*(hubble/0.7)**(-2./3.)*(omega_m/0.27)**(-1./3.)
        Tvir = (delta/200.0)**(1./3.)*Tfit_norm*(mvir/mpivot)**(Tfit_alpha)*(aexp)**(Tfit_beta)
        Kvir = (delta/200.0)**(-1./3.)*Kfit_norm*(mvir/mpivot)**(Kfit_alpha)*(aexp)**(Kfit_beta)
        Pvir = (delta/200.0)**(4./3.)*Pfit_norm*(mvir/mpivot)**(Pfit_alpha)/(aexp)**(Pfit_beta)

    return (Tvir, Pvir, Kvir)

def getMnl(z) :
    """
    Return characteristic collapse mass Mnl at a given redshift z.
    """

    #from numpy import *
    #from scipy import *
    import numpy as np
    import cosmolopy.perturbation as cp

    #Cosmological Parameters for L500 and Bolshoi
    cosmo = {'omega_M_0' : 0.27, 'omega_lambda_0' : 0.73, 'omega_b_0' : 0.0469388, 'omega_n_0' : 0.0, 'N_nu' : 0,'h' : 0.70, 'n' : 0.96, 'sigma_8' : 0.82, 'baryonic_effects': True}
    
    delta_c = 1.686

    # mass bins in Msun/h
    mbins = np.array(10**(arange(1,16,0.01)))
    nu = np.zeros (len(mbins))

    # Build nu array
    sigma = cp.sigma_r(cp.mass_to_radius(mbins/cosmo['h'], **cosmo), z, **cosmo)[0]
    nu = delta_c/sigma

    # Compute Mnl which equals M when nu = 1
    closest = where(nu < 1.)[0][-1] #nu increases with M

    # Interplolate to get Mnl
    Mnl = np.log10(mbins[closest]) + np.log10(mbins[closest+1]/mbins[closest])/np.log10(nu[closest+1]/nu[closest])*np.fabs(np.log10(nu[closest]))
    Mnl = 10**Mnl

    return Mnl

def getNu(mvir,z) :
    """
    Return peak height nu for a given halo mass Mvir and redshift z. Mvir is in Msun/h

    """
    #from numpy import *
    #from scipy import *
    import cosmolopy.perturbation as cp

    #Cosmological Parameters for L500 and Bolshoi
    cosmo = {'omega_M_0' : 0.27, 'omega_lambda_0' : 0.73, 'omega_b_0' : 0.0469388, 'omega_n_0' : 0.0, 'N_nu' : 0,'h' : 0.70, 'n' : 0.96, 'sigma_8' : 0.820, 'baryonic_effects': True}

    delta_c = 1.686

    sigma = cp.sigma_r(cp.mass_to_radius(mvir/cosmo['h'], **cosmo), z, **cosmo)[0]
    nu = delta_c/sigma

    return nu

