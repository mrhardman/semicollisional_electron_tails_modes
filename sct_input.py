# -*- coding: utf-8 -*-
import numpy as np
#from scipy.io import netcdf
from netCDF4 import Dataset
import f90nml

class sctphysicsinput:    

    def __init__(self,wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,
 wdrift,b0,b1,b2,b3,m0,m1,m2,m3):
 
        self.wstarn = wstarn
        self.wstart = wstart 
        self.tite = tite 
        self.shat = shat 
        self.beta = beta 
        self.a0 = a0 
        self.a1 = a1 
        self.a2 = a2 
        self.d0 = d0 
        self.d1 = d1 
        self.d2 = d2 
        self.f0 = f0 
        self.wdrift = wdrift 
        self.b0 = b0 
        self.b1 = b1 
        self.b2 = b2 
        self.b3 = b3
        self.m0 = m0 
        self.m1 = m1 
        self.m2 = m2 
        self.m3 = m3 
    
        print("Running with \n ")
        print("wstarn = "+str(wstarn)+ "\n")
        print("wstart = "+str(wstart)+ "\n")
        print("tite = "+str(tite)+ "\n")
        print("shat = "+str(shat)+ "\n")
        print("beta = "+str(beta)+ "\n")
        print("a0 = "+str(a0)+ "\n")
        print("a1 = "+str(a1)+ "\n")
        print("a2 = "+str(a2)+ "\n")
        print("d0 = "+str(d0)+ "\n")
        print("d1 = "+str(d1)+ "\n")
        print("d2 = "+str(d2)+ "\n")
        print("f0 = "+str(f0)+ "\n")
        print("wdrift = "+str(wdrift)+ "\n")
        print("b0 = "+str(b0)+ "\n")
        print("b1 = "+str(b1)+ "\n")
        print("b2 = "+str(b2)+ "\n")
        print("b3 = "+str(b3)+ "\n")
        print("m0 = "+str(m0)+ "\n")
        print("m1 = "+str(m1)+ "\n")
        print("m2 = "+str(m2)+ "\n")
        print("m3 = "+str(m3)+ "\n")
        
        self.list =  [self.wstarn, self.wstart, self.tite, self.shat, self.beta, self.a0, 
                                self.a1, self.a2, self.d0, self.d1, self.d2, self.f0,
                                self.wdrift, self.b0, self.b1, self.b2, self.b3,
                                self.m0, self.m1, self.m2, self.m3]
                                
class geodata: #the class of objects which contain semi-collisional mode data
    #Variables global to the class here

    total_instance_number = 0
    total_instance_in_memory = 0
    #object functions here
    def __init__(self, work_dir, run_name, output_dir):
    
        
        self.work_dir = work_dir
        self.run_name = run_name
        self.output_dir = output_dir
        geodata.total_instance_number +=1 # increment the total number of instances 
        geodata.total_instance_in_memory += 1 #Increment the total number of instances in memory
        
        with open(work_dir+"/"+run_name+ ".in") as nml_file:
            namelist = f90nml.read(nml_file)
            nml_file.close()
        q_safety_factor = namelist["theta_grid_parameters"]["qinp"]
        Rgeo = namelist["theta_grid_parameters"]["R_geo"]
        Rmaj = namelist["theta_grid_parameters"]["Rmaj"]
        rhoc = namelist["theta_grid_parameters"]["rhoc"]
        
        infile = work_dir+"/"+run_name+".out.nc" 
        #ncfile = netcdf.netcdf_file(infile,'r')
        ncfile = Dataset(infile,'r')
        #print(ncfile.variables)
        self.theta = np.copy(ncfile.variables['theta'][:]) # theta
        self.ntheta = self.theta.size
        self.gbdrift = np.copy(ncfile.variables['gbdrift'][:]) # theta
        self.gbdrift0 = np.copy(ncfile.variables['gbdrift0'][:]) # theta
        self.cvdrift = np.copy(ncfile.variables['cvdrift'][:]) # theta
        self.cvdrift0 = np.copy(ncfile.variables['cvdrift0'][:]) # theta
        self.bmag = np.copy(ncfile.variables['bmag'][:]) # theta
        self.gds2 = np.copy(ncfile.variables['gds2'][:]) # theta
        self.gds21 = np.copy(ncfile.variables['gds21'][:]) # theta
        self.gds22 = np.copy(ncfile.variables['gds22'][:]) # theta
        self.jacob = np.copy(ncfile.variables['jacob'][:]) # theta
        self.gradpar = np.copy(ncfile.variables['gradpar'][:]) # theta
        self.grho = np.copy(ncfile.variables['grho'][:]) # theta
        self.drhodpsi = np.copy(ncfile.variables['drhodpsi'][...])
        self.kxfac = np.abs((q_safety_factor/rhoc)/self.drhodpsi)#np.copy(ncfile.variables['kxfac'][...])
        #need to revisit this variable print(self.kxfac)
        self.shat = np.copy(ncfile.variables['shat'][...])
        self.Rgeo = Rgeo
        self.Rmaj = Rmaj
        self.q_safety_factor = q_safety_factor
        
        
        dtheta  = np.copy(self.theta)
        dtheta[0:self.ntheta-1] = np.subtract(self.theta[1:self.ntheta],self.theta[0:self.ntheta-1])
        dtheta[self.ntheta-1] = 0.
        Bdotgradtheta = self.bmag*self.gradpar
        wgt = np.divide(dtheta,Bdotgradtheta)
        self.wgts_theta = wgt
        #print("wgt",wgt)
        #print("sum wgt",np.sum(wgt))
        
        
def finite_aspect_ratio_coeffs(workdir, file, outputdir, beta_effective, tite, Zi, fprim, tprim, theta0, lpoly):
    from transport_coefficients import trapped_fraction
    from transport_coefficients import spitzer_harm_transport_coefficients
    from transport_coefficients import neoclassical_transport_coefficients
    from transport_coefficients import classical_transport_coefficients
    from transport_coefficients import Omicron_function
    from transport_coefficients import flux_surface_average
    
    geodata0 = geodata(workdir,file,outputdir)
    bmag = geodata0.bmag
    wgts_theta = geodata0.wgts_theta
    gradpar = geodata0.gradpar
    grho = geodata0.grho
    theta = geodata0.theta
    gds22 = geodata0.gds22
    gds21 = geodata0.gds21
    gds2 = geodata0.gds2
    gbdrift = geodata0.gbdrift
    gbdrift0 = geodata0.gbdrift0
    cvdrift = geodata0.cvdrift
    cvdrift0 = geodata0.cvdrift0
    shat = geodata0.shat 
    kxfac = geodata0.kxfac 
    q_safety_factor = geodata0.q_safety_factor
    Rmaj = geodata0.Rmaj
    Rgeo = geodata0.Rgeo
    drhodpsi = geodata0.drhodpsi
    
    ftrapped = trapped_fraction(bmag,wgts_theta)
    fc = 1. - ftrapped
    kpar = gradpar*q_safety_factor*Rmaj # change normalisation of z 
    
    arhs, crhs, aSH, cSH = spitzer_harm_transport_coefficients(Zi,lpoly)
    LNeo = neoclassical_transport_coefficients(fc,bmag,kpar,wgts_theta,Zi,lpoly)
    LCla = classical_transport_coefficients(bmag,grho,wgts_theta,Zi)
    
    #print(LNeo)
    #print(kpar**2)
    av_kpar2 = flux_surface_average(kpar**2,wgts_theta)
    #print(av_kpar2)
    
    a0 = -LNeo[2,2] - aSH[0]*av_kpar2
    a1 = -LNeo[2,3] - cSH[0]*av_kpar2
    a2 = -LNeo[3,3] - (cSH[0]- (5./2.)*cSH[1])*av_kpar2
    
    neofactor = Rgeo*drhodpsi/(q_safety_factor*Rmaj)
    clafactor = 1.0/(q_safety_factor*Rmaj)
    
    testbootstrap = False#True
    testcrossfield = False#True
    if(testbootstrap):
        m0 = 0.
        m1 = 0.
        m2 = 0.
        m3 = 0.
        b0 = 0.
        b1 = 0.
        b2 = 0.
        b3 = 0.
    else:
        m0 = 0. #Onsager Symmetry
        m1 = (LNeo[0,3] - LNeo[2,1])*neofactor
        m2 = (LNeo[1,2] - LNeo[3,0])*neofactor
        m3 = 0. #Onsager Symmetry
        
        b0 = LNeo[2,0]*neofactor
        b1 = LNeo[2,1]*neofactor
        b2 = LNeo[3,0]*neofactor
        b3 = LNeo[3,1]*neofactor
       
    if(testcrossfield):
        d0 = 0.
        d1 = 0.
        d2 = 0.
    else:    
        d0 = LCla[0,0]*(clafactor)**2 + LNeo[0,0]*(neofactor)**2
        d1 = LCla[0,1]*(clafactor)**2 + LNeo[0,1]*(neofactor)**2
        d2 = LCla[1,1]*(clafactor)**2 + LNeo[1,1]*(neofactor)**2
    
    f0 = 1. #redundant variable
    
    kperp2norm_test = gds2 + 2.0*theta0*gds21 + gds22*theta0**2
    #print("kperp2norm_test",kperp2norm_test)
    kperpgeo2 = gds22
    kperpgeo1 = gds21 + theta*gds22
    kperpgeo0 = gds2 - (theta**2)*gds22 + 2.*theta*(gds21 + theta*gds22)
    #print(kperpgeo0)
    #print(kperpgeo1)
    #print(kperpgeo2)
    Omicron = Omicron_function(theta,bmag,gradpar,kperpgeo2,kperpgeo1,kperpgeo0,wgts_theta,theta0,shat,kxfac)
    print(Omicron)
    # include Omicron with beta for now 
    #beta_effective = beta_effective*Omicron
    # need more tests of Omicron fn
    
    # normalised drive frequencies
    drift = -(1./4.)*(gbdrift + theta*gbdrift0 + cvdrift + theta*cvdrift0)
    wdrift = flux_surface_average(drift,wgts_theta)
    wstarn = -fprim/2.
    wstart = -tprim/2.
    
    # include kxfac with shat for now 
    shat_code = shat*kxfac
    
    physics_input = sctphysicsinput(wstarn,wstart,tite,shat_code,beta_effective,a0,a1,a2,d0,d1,d2,f0,wdrift,b0,b1,b2,b3,m0,m1,m2,m3)
    
    return physics_input
    