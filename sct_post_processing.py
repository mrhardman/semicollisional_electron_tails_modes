# -*- coding: utf-8 -*-
from scipy.io import netcdf
import numpy as np


class sctdata: #the class of objects which contain semi-collisional mode data
    #Variables global to the class here

    total_instance_number = 0
    total_instance_in_memory = 0
    #object functions here
    def __init__(self, work_dir, run_name, output_dir):
    
        
        self.work_dir = work_dir
        self.run_name = run_name
        self.output_dir = output_dir
        sctdata.total_instance_number +=1 # increment the total number of instances 
        sctdata.total_instance_in_memory += 1 #Increment the total number of instances in memory
        
        
        infile = work_dir+"/"+run_name+".nc" 
        ncfile = netcdf.netcdf_file(infile,'r')

        self.zed = np.copy(ncfile.variables['zed'][:]) # nzed
        self.densr = np.copy(ncfile.variables['densr'][:,:]) # nzed,nmodes
        self.densi = np.copy(ncfile.variables['densi'][:,:]) # nzed,nmodes
        self.tempr = np.copy(ncfile.variables['tempr'][:,:]) # nzed,nmodes
        self.tempi = np.copy(ncfile.variables['tempi'][:,:]) # nzed,nmodes
        self.currentr = np.copy(ncfile.variables['currentr'][:]) # nmodes
        self.currenti = np.copy(ncfile.variables['currenti'][:]) # nmodes
        self.omegar = np.copy(ncfile.variables['omegar'][:]) # nmodes
        self.omegai = np.copy(ncfile.variables['omegai'][:]) # nmodes
        
        self.dens = self.densr + 1j*self.densi
        self.temp = self.tempr + 1j*self.tempi
        self.current = self.currentr + 1j*self.currenti
        self.omega = self.omegar + 1j*self.omegai
        
        self.wstarn = np.copy(ncfile.variables['wstarn'][...])
        self.wstart = np.copy(ncfile.variables['wstart'][...])
        self.tite = np.copy(ncfile.variables['tite'][...])
        self.shat = np.copy(ncfile.variables['shat'][...])
        self.beta = np.copy(ncfile.variables['beta'][...])
        self.a0 = np.copy(ncfile.variables['a0'][...])
        self.a1 = np.copy(ncfile.variables['a1'][...])
        self.a2 = np.copy(ncfile.variables['a2'][...])
        self.d0 = np.copy(ncfile.variables['d0'][...])
        self.d1 = np.copy(ncfile.variables['d1'][...])
        self.d2 = np.copy(ncfile.variables['d2'][...])
        self.f0 = np.copy(ncfile.variables['f0'][...])
        self.wdrift = np.copy(ncfile.variables['wdrift'][...])
        self.b0 = np.copy(ncfile.variables['b0'][...])
        self.b1 = np.copy(ncfile.variables['b1'][...])
        self.b2 = np.copy(ncfile.variables['b2'][...])
        self.b3 = np.copy(ncfile.variables['b3'][...])
        self.m0 = np.copy(ncfile.variables['m0'][...])
        self.m1 = np.copy(ncfile.variables['m1'][...])
        self.m2 = np.copy(ncfile.variables['m2'][...])
        self.m3 = np.copy(ncfile.variables['m3'][...])

        self.physics_input = [self.wstarn, self.wstart, self.tite, self.shat, self.beta, self.a0, 
                                self.a1, self.a2, self.d0, self.d1, self.d2, self.f0,
                                self.wdrift, self.b0, self.b1, self.b2, self.b3,
                                self.m0, self.m1, self.m2, self.m3]
    def __del__(self):
        sctdata.total_instance_in_memory += -1 #Increment the total number of instances in memory
    