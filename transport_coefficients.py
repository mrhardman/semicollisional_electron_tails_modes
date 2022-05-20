# -*- coding: utf-8 -*-
import numpy as np
from collision_matrices import KKe_matrix, QQe_matrix, arhs_vector, crhs_vector
from numpy import matmul as cdot 

# Flux surface average of integrand
# wgts_theta = d theta /  (a B . grad theta / B0) 
# < Integrand >
def flux_surface_average(integrand, wgts_theta):
    average = np.sum(integrand*wgts_theta)/np.sum(wgts_theta)
    return average

def spitzer_harm_transport_coefficients(Zi,l):
    
    KKe = KKe_matrix(l,Zi)
    KKeInv = np.linalg.inv(KKe)
    
    arhs = arhs_vector(l)
    crhs = crhs_vector(l)
    
    aSH = cdot(KKeInv,arhs)
    cSH = cdot(KKeInv,crhs)
    
    return arhs, crhs, aSH, cSH

# Document this
def neoclassical_transport_coefficients(fc,Bmag,kpar,wgts_theta,Zi,l):
    #  < a b.grad theta >
    av_kpar = flux_surface_average(kpar, wgts_theta)
    #  < a B.grad theta / B0 >
    av_Bkpar = flux_surface_average(Bmag*kpar, wgts_theta)
    #  < a b.grad theta B0 / B >
    av_kparBm = flux_surface_average(kpar/Bmag, wgts_theta)
    #  < (B0 / B)^2 >    
    av_Bm2 = flux_surface_average(1./(Bmag**2), wgts_theta) 
    #  < (B / B0)^2 >
    av_B2 = flux_surface_average(Bmag**2, wgts_theta)
    #  < (B / B0)^2 > * fraction of circulating particles
    fc_eff = fc*av_B2

    arhs, crhs, aSH, cSH = spitzer_harm_transport_coefficients(Zi,l)

    KKe = KKe_matrix(l,Zi)
    QQe = QQe_matrix(l,Zi,fc_eff)
    
    QQeInv = np.linalg.inv(QQe)
    #KKeInv = np.linalg.inv(KKe)
    
    #arhs = arhs_vector(l)
    #crhs = crhs_vector(l)
    
    #aSH = cdot(KKeInv,arhs)
    #cSH = cdot(KKeInv,crhs)
    
    
    L00 = av_Bm2*KKe[0,0] - fc*cdot( KKe[:,0], cdot( QQeInv, KKe[:,0] ) ) 
    
    L01 = av_Bm2*(KKe[0,0] - KKe[0,1]) - fc*cdot( KKe[:,0], cdot( QQeInv, (KKe[:,0] - KKe[:,1]) ) ) 
    
    L02 = -av_kparBm*arhs[0] + av_Bkpar*fc*cdot( cdot(KKe[:,0], QQeInv), arhs) 
    
    L03 = -av_kparBm*crhs[0] + av_Bkpar*fc*cdot( cdot(KKe[:,0], QQeInv), crhs) 
    
    L11 = av_Bm2*(KKe[0,0] - KKe[0,1] - KKe[1,0] + KKe[1,1]) - fc*cdot( (KKe[:,0] - KKe[:,1]), cdot( QQeInv, (KKe[:,0] - KKe[:,1]) ) )
    
    L12 = -av_kparBm*(arhs[0] - arhs[1]) + av_Bkpar*fc*cdot( cdot( (KKe[:,0] - KKe[:,1]), QQeInv), arhs)
    
    L13 = -av_kparBm*(crhs[0] - crhs[1]) + av_Bkpar*fc*cdot( cdot( (KKe[:,0] - KKe[:,1]), QQeInv), crhs)
    
    L22 = ((av_kpar)**2)*cdot(aSH,arhs) - fc*((av_Bkpar)**2)*cdot( arhs, cdot(QQeInv, arhs))
    
    L23 = ((av_kpar)**2)*cdot(cSH,arhs) - fc*((av_Bkpar)**2)*cdot( crhs, cdot(QQeInv, arhs))
    
    L33 = ((av_kpar)**2)*cdot(cSH,crhs) - fc*((av_Bkpar)**2)*cdot( crhs, cdot(QQeInv, crhs))
    
    LNeo = np.array([[L00,L01,L02,L03],
                     [L01,L11,L12,L13],
                     [L02,L12,L22,L23],
                     [L03,L13,L23,L33]])

    return LNeo
    
def classical_transport_coefficients(Bmag,gradrho,wgts_theta,Zi):
    #gradrho = a grad rho = grad r 
    #Bmag = B / B0
    #wgts_theta = d theta /  (a B . grad theta / B0)
    #Zi ion charge number 
    
    L00 = Zi*flux_surface_average((gradrho/Bmag)**2,wgts_theta)
    L01 = -0.5*Zi*flux_surface_average((gradrho/Bmag)**2,wgts_theta)
    L11 = ((5./4.)*Zi + np.sqrt(2))*flux_surface_average((gradrho/Bmag)**2,wgts_theta)
    
    LCla = np.array([[L00,L01],
                     [L01,L11]])

    return LCla
    
# matching constant from 
# Integral [ (B / kperp2) (d theta/b.grad theta) ]
# assumes input theta grid ranges from -Pi to Pi    
def Omicron_function(theta,Bmag,kpar,geo2,geo1,geo0,wgts_theta,theta0,shat,kxfac): #gds2,gds21,gds22,
    
    # number of theta in 2 Pi segment
    ntheta = Bmag.size
    kperp2norm = np.zeros(ntheta)
    integrand = np.zeros(ntheta)
    
    sum = 0.
    nseg = 1001 # should be a large number 
    for iseg in range(0,nseg):
        thetalocal = theta + 2.*np.pi*(iseg-nseg//2)
        for itheta in range(0,ntheta):
            kperp2norm[itheta]= ( geo0[itheta]  + 2.0*(theta0-thetalocal[itheta])*geo1[itheta]  + ((theta0-thetalocal[itheta])**2)*geo2[itheta] )
            integrand[itheta] = ((Bmag[itheta])**2)/kperp2norm[itheta]
        sum = sum + np.sum(integrand*wgts_theta) 
        #print("kperp2norm",kperp2norm)    
        # computes line integral
        # Integral [ (B^2 / kperp2) (d theta/B.grad theta) ]
        # Noting that wgts_theta contains 1/B.grad theta
    # divide by  < a B.grad theta / B0 >
    av_Bkpar = flux_surface_average(Bmag*kpar, wgts_theta)
    sum = sum/av_Bkpar
    # divide by Pi/(shat*kxfac)
    # makes Omicron = 1 for a/R << 1 and circular surfaces with q R = 1
    sum = sum/(np.pi/shat*kxfac)
    print("Om test",sum*av_Bkpar)
    return sum
    
def trapped_fraction_integrand(pitch,Bmag,wgts_theta):
    I = np.sqrt(1. - pitch*Bmag)
    integrand = np.sum(wgts_theta)*pitch/np.sum(I*wgts_theta)
    return integrand

def trapped_fraction(Bmag,wgts_theta):
    from scipy.special import roots_legendre
    Bmagmax= Bmag[-1]
    pitchmax = 1/Bmagmax
    npitch = 1000
    x, w = roots_legendre(npitch)
    pitch  = pitchmax*(x + 1.)/2.
    wgts_pitch = pitchmax*w/2.
    
    integrand = np.zeros(npitch)
    for ipitch in range(0,npitch):
        integrand[ipitch] = trapped_fraction_integrand(pitch[ipitch],Bmag,wgts_theta)
    
    trapped_fraction = 1. - (3./4.)*np.sum(integrand*wgts_pitch)
    return trapped_fraction

