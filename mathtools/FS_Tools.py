from scipy import *
from pylab import *
from mathtools.utils import map_to_interval


# tool for making a multi-dimensional basis
def cartesian_product_basis(BX,BY):
    Bamp = array([x*y for x in BX.transpose() for y in BY.transpose()]).transpose()

    return Bamp


# tools for doing Fourier Series analysis on data

# construct a fourier series whose phase is 2*pi*alpha*R2
def FS_basis(R2,alpha,NFterms):
    N = len(R2)
    basis = zeros((N,2*NFterms+1),dtype='double')
    basis[:,0] = 1.0
    for k in r_[0:NFterms]:
        basis[:,k+1] = 2.0*cos(2.0*pi*(k+1)*alpha*R2)
        basis[:,k+1+NFterms] = 2.0*sin(2.0*pi*(k+1)*alpha*R2)

    return basis

# derivative of that basis w/r to phase = 2*pi*[alpha*R^2]
def FS_dbasis(R2,alpha,NFterms):
    N = len(R2)
    dbasis = zeros((N,2*NFterms+1),dtype='double')
    dbasis[:,0] = 0.0
    for k in r_[0:NFterms]:
        dbasis[:,k+1] = -2.0*sin(2.0*pi*(k+1)*alpha*R2)*(k+1)
        dbasis[:,k+1+NFterms] = 2.0*cos(2.0*pi*(k+1)*alpha*R2)*(k+1)

    return dbasis

def FS_d2basis(R2,alpha,NFterms):
    N = len(R2)
    d2basis = zeros((N,2*NFterms+1),dtype='double')
    d2basis[:,0] = 0.0
    for k in r_[0:NFterms]:
        d2basis[:,k+1] = -2.0*cos(2.0*pi*(k+1)*alpha*R2)*((k+1)**2)
        d2basis[:,k+1+NFterms] = -2.0*sin(2.0*pi*(k+1)*alpha*R2)*((k+1)**2)

    return d2basis

# Nterm fourier series coefs for unit rectangle shifted delta orders to the right
def FS_rect_coef(sigma,Nterms,delta):
    zz = (sigma*pi*r_[1:Nterms+1])
    ak = r_[1.0,sin(zz)/zz]
    ck = zeros(1+2*Nterms,dtype='double')
    ck[0:Nterms+1] = ak*cos(2.0*pi*r_[0:Nterms+1]*delta)
    ck[Nterms+1:] = ak[1:]*sin(2.0*pi*r_[1:Nterms+1]*delta)
    return ck

# Nterm fourier series coefs for gaussian shifted delta orders to the right
def FS_gaussian_coef(sigma,Nterms,delta):
    ak = exp(-0.5*(2*pi*sigma*r_[0:Nterms+1])**2)
    ck = zeros(1+2*Nterms,dtype='double')
    ck[0:Nterms+1] = ak*cos(2.0*pi*r_[0:Nterms+1]*delta)
    ck[Nterms+1:] = ak[1:]*sin(2.0*pi*r_[1:Nterms+1]*delta)
    return ck

# compute FS of an ideal Fabry-Perot response
def FS_FP_coef(R,NFterms,delta):
    ak = ((1.0-R)/(1.0+R))*R**r_[0:NFterms+1]
    ck = zeros(1+2*NFterms,dtype='double')
    ck[0:NFterms+1] = ak*cos(2.0*pi*r_[0:NFterms+1]*delta)
    ck[NFterms+1:] = ak[1:]*sin(2.0*pi*r_[1:NFterms+1]*delta)
    return ck

# compute FS of a unit function
def FS_unit(NFterms):
    ck = zeros(1+2*NFterms,dtype='double')
    ck[0] = 1.0
    return ck

# compute the 'continuous' inner product using FS coefs
def FS_inner_product(cc,dd):
    inner_product = cc[0]*dd[0] + 2*(cc[1:]*dd[1:]).sum()
    return inner_product

# shift a signal to the right delta orders - using its fourier series
def FS_shift(cc,delta):
    NFterms = (len(cc)-1)/2
    ck = 0*cc
    ck[0] = cc[0]

    kset = r_[1:NFterms+1]
    cset = r_[1:NFterms+1]
    sset = r_[NFterms+1:2*NFterms+1]

    ck[cset] = cc[cset]*cos(kset*2*pi*delta)-cc[sset]*sin(kset*2*pi*delta)
    ck[sset] = cc[cset]*sin(kset*2*pi*delta)+cc[sset]*cos(kset*2*pi*delta)
    return ck

# do the FS shift on a whole matrix of coefs using a vector of deltas
# here the FS coefs are stored vertically (each column in a new set of FS coefs)
def FS_big_shift(cc,delta):
    
    NFterms = (cc.shape[0]-1)/2
    Ncols = cc.shape[1]

    ck = 0*cc
    ck[0,:] = cc[0,:]

    kset = r_[1:NFterms+1]
    cset = r_[1:NFterms+1]
    sset = r_[NFterms+1:2*NFterms+1]
    KD = outer(kset,delta)

    ck[cset,:] = cc[cset,:]*cos(2*pi*KD)-cc[sset,:]*sin(2*pi*KD)

    ck[sset,:] = cc[cset,:]*sin(2*pi*KD)+cc[sset,:]*cos(2*pi*KD)


    return ck



# cross correlate a reference with a signal using their Fourier Series representations
def FS_xcorr(Fsig,Fref):
    Fh = zeros(Fsig.shape,dtype='double')
    NFterms = (len(Fh)-1)/2
    Fh[0] = Fsig[0]*Fref[0]

    cset = r_[1:NFterms+1]
    sset = r_[NFterms+1:2*NFterms+1]

    Fh[cset] = Fref[cset]*Fsig[cset] + Fref[sset]*Fsig[sset]
    Fh[sset] = Fref[cset]*Fsig[sset] - Fref[sset]*Fsig[cset]

    return Fh

# convolve a reference from a signal using their Fourier Series representations
def FS_conv(Fsig,Fref):
    Fh = zeros(Fsig.shape,dtype='double')
    NFterms = (len(Fh)-1)/2
    Fh[0] = Fsig[0]*Fref[0]

    cset = r_[1:NFterms+1]
    sset = r_[NFterms+1:2*NFterms+1]

    Fh[cset] = Fref[cset]*Fsig[cset] - Fref[sset]*Fsig[sset]
    Fh[sset] = Fref[cset]*Fsig[sset] + Fref[sset]*Fsig[cset]

    return Fh


# deconvolve a reference from a signal using their Fourier Series representation
# sigma2 is regularizer
def FS_deconv(Fsig,Fref,sigma2=0.0,beta1=0.0,beta2=0.0):
    Fh = zeros(Fsig.shape,dtype='double')
    NFterms = (len(Fh)-1)/2
    kvec = r_[1:NFterms+1]

    Fh[0] = Fsig[0]/Fref[0]

    cset = r_[1:NFterms+1]
    sset = r_[NFterms+1:2*NFterms+1]
    den_k = Fref[cset]**2 + Fref[sset]**2 + sigma2 + beta1*kvec + beta2*kvec**2
    Fh[cset] = (Fref[cset]*Fsig[cset] + Fref[sset]*Fsig[sset])/den_k
    Fh[sset] = (Fref[cset]*Fsig[sset] - Fref[sset]*Fsig[cset])/den_k

    return Fh

# estimate the width of the nearest gaussian to a signal using its Fourier Series representation
# use frequencies from NqI:NqF (to avoid aerosol portion)
def FS_gauss_width(Fsig,NqI,NqF,betaR=1.0e-24,iweight=0):
    NFterms = (len(Fsig)-1)/2    
    kvec = r_[0:NFterms+1]
    EB = c_[ones(NFterms+1),kvec,kvec**2] 

    ccexp = sqrt(Fsig[0:NFterms+1]**2 + r_[0.0,Fsig[NFterms+1:2*NFterms + 1]]**2) # spectrum shifted back to origin

    if iweight == 0:
        fit_ecc,res,rnk,sv = lstsq(EB[NqI:NqF,:],log(abs(ccexp))[NqI:NqF])
    else:
        weight = 1.0/((kvec)+r_[1.0,zeros(NFterms)])
        fit_ecc,res,rnk,sv = lstsq(multiply(weight,EB.transpose()).transpose()[NqI:NqF,:],(weight*log(abs(ccexp)+betaR))[NqI:NqF])

    if fit_ecc[2] >= 0:
        width = 0.0
        avg = Fsig[0]
    else:
        width = sqrt(2)*sqrt(-fit_ecc[2])/(2*pi)
        avg = exp(fit_ecc[0]) # this is the approx of the 0th FS coef == avg

    return width,avg

# more of a blunt way to do this - what fraction of data are above the mean?
def FS_width(Fsig,pBasis):
    N = pBasis.shape[0]
    recon = dot(pBasis,Fsig)
    pavg = recon.mean()
    pwidth = (recon > pavg).sum()/double(N)

    return pwidth,pavg

# yields truncated fourier series approximation of pixel space boxcar blur of
# function
def FS_blur_basis(B,delta,pix,x0,alpha):
    NFterms = (B.shape[1]-1)/2
    blrB = 1.0*B
    for k in r_[1:NFterms+1]:
        defect_k = sinc(2*alpha*delta*(pix-x0)) # variable defect coefficient
        blrB[:,k] = defect_k*B[:,k] # set the cos terms
        blrB[:,k+NFterms] = defect_k*B[:,k+NFterms] # set the sin terms

    return blrB


# Basis lof Legendre Polynomials
# either x needs to be defined on [-1,1] or we need to center and normalize it
def Legendre_basis(Nterms,abscissa,icenter=1,inorm=1):
    N = len(abscissa)

    if icenter != 0:
        x0 = abscissa.mean()
    else:
        x0 = 0.0
        
    if inorm != 0:
        L = 0.5*((abscissa-x0).max() - (abscissa-x0).min())
    else:
        L = 1.0

    x = (abscissa-x0)/L
    x = map_to_interval(abscissa, [-1, 1])

    Lbasis = zeros((N,Nterms),dtype='double')
    Pkm1 = ones(N,dtype='double')
    Pk = x
    Lbasis[:,0] = Pkm1
    Lbasis[:,1] = Pk
    for k in r_[2:Nterms]:
        m = k-1
        c = (2.0*m+1.0)
        Pkp1 = (c*x*Pk - m*Pkm1)/double(m+1)
        Lbasis[:,k] = Pkp1
        Pkm1 = Pk
        Pk = Pkp1
        

    return Lbasis

# derivative Basis lof Legendre Polynomials 
# either x needs to be defined on [-1,1] or we need to center and normalize it
# tjis will be the variable we take derivative with respect to
def dLegendre_basis(Nterms,abscissa,icenter=1,inorm=1):
    N = len(abscissa)

    if icenter != 0:
        x0 = abscissa.mean()
    else:
        x0 = 0.0
        
    if inorm != 0:
        L = 0.5*((abscissa-x0).max() - (abscissa-x0).min())
    else:
        L = 1.0

    x = (abscissa-x0)/L
    x = map_to_interval(abscissa, [-1, 1])

    dLbasis = zeros((N,Nterms),dtype='double')
    Pkm1 = ones(N,dtype='double')
    Pk = x
    dPkm1 = zeros(N,dtype='double')
    dPk = ones(N,dtype='double')

    dLbasis[:,0] = dPkm1
    dLbasis[:,1] = dPk
    for k in r_[2:Nterms]:
        m = k-1
        c = (2.0*m+1.0)
        Pkp1 = (c*x*Pk - m*Pkm1)/double(m+1)
        dPkp1 = (c*(Pk+x*dPk) - m*dPkm1)/double(m+1)
        dLbasis[:,k] = dPkp1
        Pkm1 = Pk
        Pk = Pkp1
        dPkm1 = dPk
        dPk = dPkp1        

    return dLbasis

# derivative Basis lof Legendre Polynomials 
# either x needs to be defined on [-1,1] or we need to center and normalize it
# tjis will be the variable we take derivative with respect to
def d2Legendre_basis(Nterms,abscissa,icenter=1,inorm=1):
    N = len(abscissa)

    if icenter != 0:
        x0 = abscissa.mean()
    else:
        x0 = 0.0
        
    if inorm != 0:
        L = 0.5*((abscissa-x0).max() - (abscissa-x0).min())
    else:
        L = 1.0

    x = (abscissa-x0)/L
    x = map_to_interval(abscissa, [-1, 1])

    d2Lbasis = zeros((N,Nterms),dtype='double')
    Pkm1 = ones(N,dtype='double')
    Pk = x
    dPkm1 = zeros(N,dtype='double')
    dPk = ones(N,dtype='double')

    d2Pkm1 = zeros(N,dtype='double')
    d2Pk = zeros(N,dtype='double')

    
    d2Lbasis[:,1] = d2Pkm1
    d2Lbasis[:,2] = d2Pk
    for k in r_[2:Nterms]:
        m = k-1
        c = (2.0*m+1.0)
        Pkp1 = (c*x*Pk - m*Pkm1)/double(m+1)
        dPkp1 = (c*(Pk+x*dPk) - m*dPkm1)/double(m+1)
        d2Pkp1 = (c*(dPk + x*d2Pk + dPk) - m*d2Pkm1)/double(m+1)

        d2Lbasis[:,k] = d2Pkp1

        Pkm1 = Pk
        Pk = Pkp1
        dPkm1 = dPk
        dPk = dPkp1        
        d2Pkm1 = d2Pk
        d2Pk = d2Pkp1        


    return d2Lbasis


# given the fourier series of a signal and a the fourier series of a positive kernel
# produce an estimate of the deconvolution of the reference from the signal
# with the restriction that the result be positive
# Basis is the FS basis used to construct the FS coefs
# cx0 is the fourier series of the initial guess
# cx0 = zeros(2NFterms+1) cx0[0] = cc_sig[0]/cc_ref[0] is a good starting point
#
def FS_pos_deconv(cc_sig,cc_ref,cx0,Basis,Niter,thresh=1.0e-2,areg=1.0e-6):
    NFterms = (len(cc_sig)-1)/2
    # will need 
    # FS_conv(Fsig,Fref)
    # FS_xcorr(Fsig,Fref)
    # and maybe: FS_shift(cc,delta)
    b = dot(Basis,cc_sig)
    #cx0 = zeros(2*NFterms+1)
    
    #cx0[0] = cc_sig[0]/cc_ref[0]
    x0 = abs(dot(Basis,cx0))

    for k in r_[0:Niter]:
        cy0 =  FS_conv(cx0,cc_ref)
        y0 = dot(Basis,cy0)
        cc_byk,res,rnk,sv = lstsq(Basis,b/y0) 
        ccAsbyk = FS_xcorr(cc_byk,cc_ref)
        pk = (dot(Basis,ccAsbyk) - areg*x0)/cc_ref[0]

        x1 = x0*pk
        x0 = x1
        cx0,res,rnk,sv = lstsq(Basis,x0)

        check_cvg = abs(log(abs(pk))).max()
        if (check_cvg < thresh):
            break

    return x0,cx0,k


# return the kernel and background from deconvolving ref from sig shifted to origin
# for example cc_ref = (1.0/Tfac)*cFP0
# shift to origin
# trick - 1st step: use a standard (oscillatory deconv) - grab center of reconstruction to estimate bck
#       - 2nd step: use a positive deconv with DC offset by bck
# cc_k2,background,sig_k2,approx_sig,iter2 = deconv_two_step_w_shift(cc_sig, cc_ref, Basis, beta2=1.0e-2, Niter_pos=30, thresh_pos=1.0e-2,areg=1.0e-10)
def deconv_two_step_w_shift(cc_sig, cc_ref, Basis, beta2=1.0e-2, Niter_pos=30, thresh_pos=1.0e-2,areg=1.0e-10):
    NFterms = (len(cc_sig)-1)/2
    Nsamples = Basis.shape[0]
    delta = arctan2(cc_sig[NFterms+1],cc_sig[1])/(2.0*pi)
    cS = FS_shift(cc_sig,-delta)

    # 1st estimate of the kernel (with background)
    cc_k1 = FS_deconv(cS,cc_ref,beta2=beta2) 
    sig_k1 = dot(Basis,cc_k1)
    # estimate background from assumption that kernel should have compact support - so any offset is bck
    background = sig_k1[Nsamples/2]

    # now try a positive deconv
    cx0 = 0*cc_k1 # initial guess
    cx0[0] = cc_sig[0]

    dc_adjust = 0.0*cc_k1
    dc_adjust[0] = background
    sig_k2,cc_k2,iter2 = FS_pos_deconv(cS-dc_adjust,cc_ref,cx0-dc_adjust,Basis,Niter_pos,thresh=thresh_pos,areg=areg)

    # a quick reconstruction of the signal as a convolution of the positive kernel (+ bck) with the reference
    # to assisit in error analysis
    cc_approx = FS_conv(cc_k2+dc_adjust,cc_ref)
    approx_sig = dot(Basis,cc_approx)
    
    return cc_k2,background,delta,sig_k2,approx_sig,iter2  


# build a weighted projection for L2 penalizing fit and smoothness, with segmented interest 
def weighted_projection(Basis,dBasis,d2Basis,data,idx_fit,idx_exclude,b0,b1fit,b1exc,b2fit,b2exc):
    Nsamples,Ncoef = Basis.shape
    B = r_[Basis[idx_fit,:],b0*eye(Ncoef),\
               b1fit*dBasis[idx_fit,:],b2fit*d2Basis[idx_fit,:],\
               b1exc*dBasis[idx_exclude,:],b2exc*d2Basis[idx_exclude,:]]

    f = r_[data[idx_fit],zeros(Ncoef),0.0*data[idx_fit],0.0*data[idx_fit],0.0*data[idx_exclude],0.0*data[idx_exclude]]
    cc,res,rnk,sv = lstsq(B,f)

    return cc


# Basis lof Legendre Polynomials
# either x needs to be defined on [-1,1] or we need to center and normalize it
def Poly_basis(Nterms,abscissa,icenter=1,inorm=1):
    N = len(abscissa)

    if icenter != 0:
        x0 = abscissa.mean()
    else:
        x0 = 0.0
        
    if inorm != 0:
        L = 0.5*((abscissa-x0).max() - (abscissa-x0).min())
    else:
        L = 1.0

    x = (abscissa-x0)/L

    Pbasis = zeros((N,Nterms),dtype='double')
    
    for k in r_[0:Nterms]:
        Pbasis[:,k]=x**k
        
    return Pbasis

def dPoly_basis(Nterms,abscissa,icenter=1,inorm=1):
    N = len(abscissa)

    if icenter != 0:
        x0 = abscissa.mean()
    else:
        x0 = 0.0
        
    if inorm != 0:
        L = 0.5*((abscissa-x0).max() - (abscissa-x0).min())
    else:
        L = 1.0

    x = (abscissa-x0)/L

    dPbasis = zeros((N,Nterms),dtype='double')
    
    for k in r_[1:Nterms]:
        dPbasis[:,k]=k*(x**(k-1))
        
    return dPbasis

def d2Poly_basis(Nterms,abscissa,icenter=1,inorm=1):
    N = len(abscissa)

    if icenter != 0:
        x0 = abscissa.mean()
    else:
        x0 = 0.0
        
    if inorm != 0:
        L = 0.5*((abscissa-x0).max() - (abscissa-x0).min())
    else:
        L = 1.0

    x = (abscissa-x0)/L

    d2Pbasis = zeros((N,Nterms),dtype='double')
    
    for k in r_[2:Nterms]:
        d2Pbasis[:,k]=(k-1)*k*(x**(k-2))
        
    return d2Pbasis


# quick and dirty cubic spline: 
# return the four basis vectors values over the unit interval t in [0:1]
def cubic_spline_interval(t):
    Nsamples = len(t)
    u = 1.0-t
    v = t
    N = zeros((Nsamples,4),dtype='double')
    c = 1.0/6.0
    q = 1.0+u*v
    N[:,0] = c*(u**3)
    N[:,1] = c + 0.5*u*q
    N[:,2] = c + 0.5*v*q
    N[:,3] = c*(v**3)

    return N

def dcubic_spline_interval(t):
    Nsamples = len(t)
    u = 1.0-t
    v = t
    dudt = -1.0
    dvdt = 1.0

    dN = zeros((Nsamples,4),dtype='double')

    q = 1.0+u*v
    dqdu = v
    dqdv = u

    dN[:,0] = 0.5*(u**2)*dudt
    dN[:,1] = 0.5*((u*dqdu + q)*dudt + (u*dqdv)*dvdt)
    dN[:,2] = 0.5*((v*dqdu)*dudt + (v*dqdv + q)*dvdt)
    dN[:,3] = 0.5*(v**2)*dvdt

    return dN

def d2cubic_spline_interval(t):
    Nsamples = len(t)
    u = 1.0-t
    v = t
    dudt = -1.0
    dvdt = 1.0

    d2N = zeros((Nsamples,4),dtype='double')

    q = 1.0+u*v
    dqdt = dudt*v + u*dvdt
    dqdu = v
    dqdv = u
    d2qdudt = dvdt
    d2qdvdt = dudt


    d2N[:,0] = u*(dudt**2)
    #d2N[:,1] = 0.5*((u*dqdu + q)*dudt + (u*dqdv)*dvdt)
    d2N[:,1] = 0.5*((dudt*dqdu + u*d2qdudt + dqdt)*dudt + (dudt*dqdv+u*d2qdvdt)*dvdt)
    #d2N[:,2] = 0.5*((v*dqdu)*dudt + (v*dqdv + q)*dvdt)
    d2N[:,2] = 0.5*((dvdt*dqdu+v*d2qdudt)*dudt + (dvdt*dqdv + v*d2qdvdt + dqdt)*dvdt)

    d2N[:,3] = v*(dvdt**2)

    return d2N


# knots = x.min() + (x.max()-x.min())*r_[-1:Nknots]/double(Nknots-2) # unifom knots
def get_cubic_spline_basis(x,knots):
    Nsamples = len(x)
    Nknots = len(knots)

    Basis = zeros((Nsamples,Nknots),dtype='double')

    for k in r_[1:Nknots-2]:
        idx_k = nonzero((x >= knots[k]) & (x <= knots[k+1]))
        t = (x[idx_k] - knots[k])/(knots[k+1]-knots[k])
        
        local_basis = cubic_spline_interval(t)
        
        Basis[idx_k,k-1] = local_basis[:,0] 
        Basis[idx_k,k] = local_basis[:,1] 
        Basis[idx_k,k+1] = local_basis[:,2] 
        Basis[idx_k,k+2] = local_basis[:,3] 


    return Basis


def get_dcubic_spline_basis(x,knots):
    Nsamples = len(x)
    Nknots = len(knots)

    dBasis = zeros((Nsamples,Nknots),dtype='double')

    for k in r_[1:Nknots-2]:
        idx_k = nonzero((x >= knots[k]) & (x <= knots[k+1]))
        h = knots[k+1]-knots[k]
        t = (x[idx_k] - knots[k])/h
        
        dlocal_basis = dcubic_spline_interval(t)
        
        dBasis[idx_k,k-1] = dlocal_basis[:,0]/h
        dBasis[idx_k,k] = dlocal_basis[:,1]/h
        dBasis[idx_k,k+1] = dlocal_basis[:,2]/h 
        dBasis[idx_k,k+2] = dlocal_basis[:,3]/h 

    return dBasis
    
def get_d2cubic_spline_basis(x,knots):
    Nsamples = len(x)
    Nknots = len(knots)

    d2Basis = zeros((Nsamples,Nknots),dtype='double')

    for k in r_[1:Nknots-2]:
        idx_k = nonzero((x >= knots[k]) & (x <= knots[k+1]))
        h = knots[k+1]-knots[k]
        t = (x[idx_k] - knots[k])/h
        
        d2local_basis = d2cubic_spline_interval(t)

        h2 = h**2
        
        d2Basis[idx_k,k-1] = d2local_basis[:,0]/h2
        d2Basis[idx_k,k] = d2local_basis[:,1]/h2
        d2Basis[idx_k,k+1] = d2local_basis[:,2]/h2 
        d2Basis[idx_k,k+2] = d2local_basis[:,3]/h2

    return d2Basis

def get_uniform_knots(x,Nknots):
    knots = x.min() + (x.max()-x.min())*r_[-1:Nknots-1]/double(Nknots-3)
    return knots


# regularize smoothnes with weighted penalty
# d is data ordinates (array of Npoints)
# xp is abscissa (array of Npoints)
# weight array of N points to weight the impact of the smoothness penalty (larger is more)
# areg1, areg2 1st and second derivative wieghting constants
def weighted_spline_fit(d,xp,weight,Nknots,areg1,areg2):
    smooth_knots = get_uniform_knots(xp,Nknots)
    Basis = get_cubic_spline_basis(xp,smooth_knots)
    dBasis = get_dcubic_spline_basis(xp,smooth_knots) 
    d2Basis = get_d2cubic_spline_basis(xp,smooth_knots)

    B = r_[Basis,areg1*multiply(weight,dBasis.transpose()).transpose(),areg2*multiply(weight,d2Basis.transpose()).transpose()]
    RHS = r_[d,0*d,0*d]
    
    cc,res,rnk,sv = lstsq(B,RHS)
    fit = dot(Basis,cc)
    dfit = dot(dBasis,cc)
    d2fit = dot(d2Basis,cc)

    return fit,dfit,d2fit

# compute smooth spline and derivative
def pos_spline_fit(data,xp,Nknots,areg1,areg2,ewL=1.0,ewR=100.0,rhl=0.95):
    d0 = data.min() - 2.0*abs(data.min())

    fixd = data - d0

    eweight = ewL*(1.0-exp(-xp/xp.max())) + ewR*exp(xp-rhl*xp.max())


    efit,defit,d2efit = weighted_spline_fit(log(fixd),xp,eweight,Nknots,areg1,areg2)

    fit = exp(efit) + d0
    dfit = exp(efit)*defit

    return fit,dfit,eweight

# build the set of basis tools needed for spline fitting with uniform knots and
# weighted regularization of 1st and 2nd derivatives
def build_weighted_spline_tools(xp,weight,Nknots,areg1,areg2):
    smooth_knots = get_uniform_knots(xp,Nknots)
    Basis = get_cubic_spline_basis(xp,smooth_knots)
    dBasis = get_dcubic_spline_basis(xp,smooth_knots) 
    d2Basis = get_d2cubic_spline_basis(xp,smooth_knots)

    B = r_[Basis,areg1*multiply(weight,dBasis.transpose()).transpose(),areg2*multiply(weight,d2Basis.transpose()).transpose()]    
    U,S,VH = svd(B,full_matrices=0)
    invB = dot(VH.transpose(),multiply(1.0/S,U).transpose())
    
    return Basis,dBasis,d2Basis,B,invB

def fit_weighted_lstsq_with_bases(data,Basis,dBasis,d2Basis,B,invB):
    rhs = r_[data,0*data,0*data]
    cc = dot(invB,rhs)
    fit = dot(Basis,cc)
    dfit = dot(dBasis,cc)
    
    residual = fit-data
    return fit,dfit,residual
    


# how about a couple of unit tests?
def FS_Test():
    N = 1000
    NFterms = 20
    alpha = 1.0
    
    x = r_[0:N]/double(N)
    dx = x[2]-x[1]

    Basis = FS_basis(x,alpha,NFterms)


    dshift = 0.3
    R = 0.76
    cc_FP = FS_FP_coef(R,NFterms,dshift)    
    cc_FPs = FS_shift(cc_FP,-dshift) # FP shift back to origin

    sigma = 0.1
    cc_gauss = FS_gaussian_coef(sigma,NFterms,0.0) # centered gaussian - unit integral area

    cc_conv = FS_conv(cc_FP,cc_gauss) # convolution
    
    # project coef onto basis for reconstruction
    yFP = dot(Basis,cc_FP) 
    yFPs = dot(Basis,cc_FPs)

    ygauss = dot(Basis,cc_gauss)
    yconv = dot(Basis,cc_conv)

    cc_unit = FS_unit(NFterms)
    yunit = dot(Basis,cc_unit) # should end up being all ones...

    figure(9);clf();hold(True)
    plot(x,yFP,'-b.',mec='b',label='FP')
    plot(x,yFPs,'-r.',mec='r',label='FPshift')
    hold(False)
    grid(True)
    legend(loc=0)
    title('FP signal at $\delta$ = {0} - and shift back'.format(dshift))
    
    figure(10);clf();hold(True)
    plot(x,yFP,'-b.',mec='b',label='FP')
    plot(x,ygauss,'-k',mec='k',label='G')
    plot(x,yconv,'-r.',mec='r',label='FP*G')
    hold(False)
    grid(True)
    legend(loc=0)
    title('FP signal convolved with Gaussian')

    FGA = FS_inner_product(cc_gauss,cc_unit)
    FGE = FS_inner_product(cc_gauss,cc_gauss)
    NGA = dx*(ygauss.sum())
    NGE = dx*((ygauss**2).sum())
    print 'Area of Gaussian FS/Numerical: {0:g} / {1:g},\nEnergy of Gaussian FS/Numerical: {2:g} / {3:g}'.format(FGA,NGA,FGE,NGE)
    print 'These are identical - so - Awesome!'

    # example of deconvolution...
    noise = (1.0e-6)*randn(N)
    ydata = yconv + noise
    cc_data,res,rnk,sv = lstsq(Basis,ydata) # first get coefs from data
    #beta2 = 1.0e-26 # regularization - for near perfect data
    #beta2 = 3.0e-8 # regularization - for data with noise
    beta2 = 1.0e-5
    #beta2 = dx*((noise*noise).sum())/FS_inner_product(cc_conv,cc_conv) # ideal regularization

    cc_dcv = FS_deconv(cc_data,cc_gauss,sigma2=1.0e-18,beta1=1.0e-14,beta2=0*1.0e-13)
    ydcv = dot(Basis,cc_dcv)

    figure(11);clf();hold(True)
    plot(x,ydata,'k.',mec='k',label='FP*G')
    #plot(x,yconv,'-k',mec='k',label='FP*G')
    plot(x,ydcv,'-r',mec='r',label='dcv')
    plot(x,yFP,'-b',mec='b',label='FP')
    hold(False)
    title('deconvolution example with $\gamma^2$ = {0:g}'.format(beta2))
    legend(loc=0)
    grid(True)

    figure(44);clf();hold(True)
    plot(cc_FP,'-b.',mec='b',label='Truth')
    plot(cc_dcv,'ro',mec='r',label='DCV')
    hold(False)
    grid(True)
    title('coefficients')
    legend(loc=0)
    
    # New example: find the width of a gaussian from fourier coefs
    ccexp = sqrt(cc_gauss[0:NFterms+1]**2 + r_[0.0,cc_gauss[NFterms+1:2*NFterms + 1]]**2)

    kvec = r_[0:NFterms+1]
    EB = c_[ones(NFterms+1),kvec,kvec**2]
    fit_ecc,res,rnk,sv = lstsq(EB,log(ccexp))
    synth_ec = dot(EB,fit_ecc)

    figure(20);clf();hold(True)
    plot(log(ccexp[0:NFterms+1]),'-b.',label='cc_exp');
    plot(synth_ec,'-r',label='synth from quadratic');
    grid(True)
    title('log of fourier coefs vs projection onto quadratic freqs')
    hold(False)
    print 'sigma/est {0:g} / {1:g}'.format(sigma,sqrt(2)*sqrt(-fit_ecc[2])/(2*pi)), 'Awesome Again !!'
    legend(loc=0)

    # now let's do it the other way around - yank the gaussian out of the data - and estimate its width...
    weight = 1.0/((kvec)+r_[1.0,zeros(NFterms)])
    
    cc_dcv_g = FS_deconv(cc_data,cc_FP,sigma2=1.0e-12,beta1=0.0,beta2=1.0e-6)
    synth_gg = dot(Basis,cc_dcv_g)

    # because gaussian is even and coefs are positive - can collapse everthing to first NFterms+1 terms
    ccgg = sqrt(cc_dcv_g[0:NFterms+1]**2 +r_[0,cc_dcv_g[NFterms+1:2*NFterms+1]]**2)
    #betaR = 1.0e-24
    betaR = 1.0e-24
    Nqr = 5

    #fit_ecc0,res,rnk,sv = lstsq(EB[0:Nqr,:],log(abs(cc_gauss[0:NFterms+1])+betaR)[0:Nqr])
    #fit_ecc1,res,rnk,sv = lstsq(EB[0:Nqr,:],log(abs(ccgg)+betaR)[0:Nqr])
    fit_ecc2,res,rnk,sv = lstsq(multiply(weight,EB.transpose()).transpose()[0:Nqr,:],(weight*log(abs(ccgg)+betaR))[0:Nqr])

    figure(50);clf();hold(True);
    plot(log(abs(cc_dcv_g)),'-b.',label='data_prj');
    grid(True);
    plot(log(abs(cc_gauss)+betaR),'r.',mec='r',label='gaussian');
    hold(False)
    title('log of coefs')
    legend(loc=0)

    figure(22);clf();hold(True)
    plot(cc_gauss[0:NFterms+1],'-b.',mec='b',label='ideal');grid(True)
    #plot(cc_dcv_g[0:NFterms+1],'-r');
    plot(ccgg,'ro',mec='r',label='data');
    hold(False)
    title('coefs of gaussian')
    legend(loc=0)

    figure(21);clf();hold(True)
    plot(x,ygauss,'-b.',label='gaussian')
    plot(x,synth_gg,'r.',label='recon')
    grid(True)
    hold(False)
    title('gaussian reconstruction')
    legend(loc=0)

    print 'sigma/est2 {0:g} / {1:g}, using 0:{2:d}'.format(sigma,sqrt(2)*sqrt(-fit_ecc2[2])/(2*pi),Nqr), 'Awesome Yet Again (with regularization)!!!'
    #print 'sigma/est1 {0:g} / {1:g}, using 0:{2:d}'.format(sigma,sqrt(2)*sqrt(-fit_ecc1[2])/(2*pi),Nqr), 'Awesome Yet Again (w/o regularization)!!!'

    # Now do it with our canned width calculator...
    Nqi = 1
    #width1,avg1 = FS_gauss_width(cc_dcv_g,Nqi,Nqr,betaR=betaR,iweight=0) # regularization is better
    width2,avg2 = FS_gauss_width(cc_dcv_g,Nqi,Nqr,betaR=betaR,iweight=1)

    #print 'function: sigma/est1 {0:g} / {1:g}, using {2:d}:{3:d}'.format(sigma,width1,Nqi,Nqr) # regularization is better
    print 'function: sigma/est2 {0:g} / {1:g}, using {2:d}:{3:d}'.format(sigma,width2,Nqi,Nqr)
 
    # Basic width process is...
    #cc_dcv_g = FS_deconv(cc_data,cc_FP,sigma2=1.0e-12,beta1=0.0,beta2=1.0e-6)
    #width2,avg2 = FS_gauss_width(cc_dcv_g,Nqi,Nqr,betaR=betaR,iweight=1)

    # what if there is some multiplier - will we still get the right width and average? - YEP!
    cc_g3 = 15.0*cc_gauss
    width3,avg3 = FS_gauss_width(cc_g3,Nqi,Nqr,betaR=betaR,iweight=1)
    print 'width3, avg3 {0:g}, {1:g}'.format(width3,avg3)
        
    Bprd = np.load('basis.npy')
    zprd = np.load('dmn.npy')
    cdata = np.load('dcv_data.npy')
    cusky1 = cdata[:,0]
    curef = cdata[:,1]
    cc_dcv_S1 = FS_deconv(curef,curef,sigma2=1.0e-12,beta1=1.0e-6,beta2=1.0e-4) # deconv ref1

    fringe_dcv = dot(Bprd,cc_dcv_S1)
    figure(2032);clf();plot(zprd,fringe_dcv,'-b.',mec='b');grid(True);
    title('{0}, {1}, {2}'.format(FS_width(cusky1,Bprd)[0],FS_width(curef,Bprd)[0],FS_width(cc_dcv_S1,Bprd)[0]))


if __name__ == '__main__':
    ion()
    N = 1000
    x = 4*r_[0:N]/double(N)
    y = exp(-x)*cos(2.0*pi*x)

    Nknots = 31

    knots = get_uniform_knots(x,Nknots)
    B = get_cubic_spline_basis(x,knots)
    dB =get_dcubic_spline_basis(x,knots)

    cc,res,rnk,sv = lstsq(B,y)
    g = dot(B,cc)

    dy = diff(y)/(x[2]-x[1])
    dg = dot(dB,cc)


    figure(11);clf();hold(True)
    plot(x,y,'-bo',mec='b')
    plot(x,g,'r.',mec='r')
    grid(True)
    hold(False)
    title('Spline Fit')

    figure(12);clf();hold(True)
    plot(x[1:],dy,'-bo',mec='b')
    plot(x,dg,'r.',mec='r')
    grid(True)
    hold(False)
    title('Spline Derivative')
    
    # NOW - Lets do it again - but with NON-Uniform knots
    kn2 = r_[-0.1, 0.0, 1.0, 2.0, 2.5, 3.0, 4.5, 5.321,5.67]
    B2 = get_cubic_spline_basis(x,kn2)
    dB2 =get_dcubic_spline_basis(x,kn2)
    
    figure(888);clf();plot(x,B2);grid(True)
