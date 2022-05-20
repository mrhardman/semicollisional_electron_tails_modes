# -*- coding: utf-8 -*-
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg') # this line allows plots to be made without using a display environment variable
from os import getcwd
from matplotlib.backends.backend_pdf import PdfPages
from utils import plot_1d_list_pdf,plot_1d_loglog_list_pdf 
import numpy as np
from scipy.special import roots_legendre
from scipy.optimize import curve_fit    
from utils import lin_func, round_sig, power_func, quad_func, linsqrt_func
from utils import lin_func2, quad_func2

from transport_coefficients import neoclassical_transport_coefficients
from transport_coefficients import classical_transport_coefficients
from transport_coefficients import Omicron_function
from transport_coefficients import trapped_fraction
def Bmag_func(B0,rminor,Rmajor,theta):
    #large R/r approximation of B
    Bmag = B0*( 1. - (rminor/Rmajor)*np.cos(theta))
    return Bmag 
def jacob_func(q,Rmajor,B0,theta):
    #large R/r approximation of 1/B.grad theta
    jacob = q*Rmajor/B0
    return jacob

def plot_test_with_sqrt_epsilon(pdf,title,test_value,slope,sqrteps,trapped_fraction=False):

    ylab_list = ["Numerical",]
    field_list = [test_value]
    a_fit = slope 
    if(trapped_fraction):
        fit = lin_func2(sqrteps,a_fit)
        ylab_list.append("Fit $= f^0_{\\rm t} =$ "+ str("%.4g" % round_sig(a_fit,sig=4))+"$\\sqrt{\\epsilon}$")
    else:
        fit = lin_func2(sqrteps,a_fit*1.462) #1.462 is prefactor of trapped fraction
        #ylab_list.append("Fit $=$ "+ str("%.6g" % round_sig(a_fit,sig=6))+"$\\sqrt{\\epsilon}$")
        ylab_list.append("Fit $=$ "+ str("%.3g" % round_sig(a_fit,sig=3))+"$f^0_{\\rm t}$")
    field_list.append(fit)
    xlist = [sqrteps,sqrteps,sqrteps]
    marker_list = ["bx","k","g","k","y"]
    xlab = "$\\sqrt{\\epsilon}$"
    legend_title = ""
    plot_1d_loglog_list_pdf  (xlist,field_list,marker_list,xlab, pdf,
      title=title,ylab='',xlims=None,ylims=None,aspx=9,aspy=6, xticks = None, yticks = None,
      markersize=10, legend_title=legend_title, use_legend=True,loc_opt='lower right', ylab_list = ylab_list,
      bbox_to_anchor_opt=(1.0, 0.0), legend_fontsize=25, ncol_opt=1)

ntheta = 10
x, w = roots_legendre(ntheta)
#print(x,w)
theta = x*np.pi

shat = 0.8    
theta0 = np.pi/2. #0.0   
Rmajor = 1.0
B0 = 1.0
rminor_list = np.logspace(-5,-0.1,num=20)
print(rminor_list)
neps = rminor_list.size
sqrteps = np.sqrt(rminor_list/Rmajor)
ftrapped = np.zeros(neps)
LNeo = np.zeros([neps,4,4])
LCla = np.zeros([neps,2,2])
for iel in range(0,neps):
    Bmag = np.zeros(ntheta)
    for itheta in range(ntheta):
        Bmag[itheta] = Bmag_func(B0,rminor_list[iel],Rmajor,theta[itheta])
    Zi = 1.
    l = 5 # number of Sonine polys
    kpar = np.ones(ntheta)
    wgts_theta = w*np.pi/(kpar*Bmag) # d theta / a B . grad theta / B0 
    ftrapped[iel] = trapped_fraction(Bmag,wgts_theta)
    fc = 1. - ftrapped[iel]
    
    LNeo[iel,:,:] = neoclassical_transport_coefficients(fc,Bmag,kpar,wgts_theta,Zi,l)
    #print(LNeo[iel,:,:])
    grho = np.ones(ntheta)
    LCla[iel,:,:] = classical_transport_coefficients(Bmag,grho,wgts_theta,Zi)
    #print(LCla[iel,:,:])
    
    #geo0 = np.ones(ntheta)
    #geo1 = np.zeros(ntheta)
    #geo2 = (shat**2)*np.ones(ntheta)
    #kxfac = 1.0
    #Omicron = Omicron_function(theta,Bmag,kpar,geo2,geo1,geo0,wgts_theta,theta0,shat,kxfac)
    #print(Omicron)
    
#print(ftrapped)    
file= getcwd() + "/tests/ftrapped.LNeo.smalleps.test.pdf"
print(file)
pdf = PdfPages(file)
plot_test_with_sqrt_epsilon(pdf, "$f_{\\rm t} = f_{\\rm trapped}$", ftrapped,1.462,sqrteps,trapped_fraction=True)      
plot_test_with_sqrt_epsilon(pdf,"$L^N_{00}$", LNeo[:,0,0],1.53,sqrteps)      
plot_test_with_sqrt_epsilon(pdf,"$-L^N_{01}$", -LNeo[:,0,1],0.59,sqrteps)      
plot_test_with_sqrt_epsilon(pdf,"$ L^N_{02} $", LNeo[:,0,2],1.66,sqrteps)      
plot_test_with_sqrt_epsilon(pdf,"$ L^N_{03} $", LNeo[:,0,3],1.75,sqrteps)      
plot_test_with_sqrt_epsilon(pdf,"$ L^N_{11}$", LNeo[:,1,1],(2.51-0.59),sqrteps)      
plot_test_with_sqrt_epsilon(pdf,"$ L^N_{12} $", LNeo[:,1,2],(1.66-1.19),sqrteps)      
plot_test_with_sqrt_epsilon(pdf,"$ L^N_{13} $", LNeo[:,1,3],(0.11+1.75),sqrteps)      
plot_test_with_sqrt_epsilon(pdf,"$L^N_{22}  $", LNeo[:,2,2],2.55,sqrteps)      
plot_test_with_sqrt_epsilon(pdf,"$ L^N_{23}  $", LNeo[:,2,3],3.51,sqrteps)      
plot_test_with_sqrt_epsilon(pdf,"$ L^N_{33}  $", LNeo[:,3,3],(3.50547754567+2.98999156016),sqrteps)      
pdf.close()     


