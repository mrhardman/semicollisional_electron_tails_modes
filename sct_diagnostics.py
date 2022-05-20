# -*- coding: utf-8 -*-
import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg') # this line allows plots to be made without using a display environment variable
from matplotlib.backends.backend_pdf import PdfPages
from utils import plot_1d_list_pdf, round_sig
from scipy.io import netcdf

    
def save_diagnostic_data(work_dir,run_name,data):
    file_name  = work_dir + run_name + ".nc"
    print(file_name)
    [zed,dens,temp,current,omega,pin] = data
    [wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,wdrift,b0,b1,b2,b3,m0,m1,m2,m3] = pin
    
    [nzed,nmodes] = np.shape(dens)
    
    f = netcdf.netcdf_file(file_name, 'w') #We open a new netcdf file
    f.history = "A (semi)-collisional tails mode output file \n Created by sct_main.py \n Author: Michael R. Hardman michael.hardman@physics.ox.ac.uk"
    
    #Create the dimensions
    f.createDimension('nzed', nzed)
    f.createDimension('nmodes', nmodes)
    
    #Create the variables
    zed_handle = f.createVariable('zed', zed.dtype, ('nzed',))
    zed_handle[:] = zed
    
    densr = dens.real 
    densr_handle=f.createVariable('densr', densr.dtype, ('nzed','nmodes',))
    densr_handle[:,:] = densr
    
    densi = dens.imag 
    densi_handle=f.createVariable('densi', densi.dtype, ('nzed','nmodes',))
    densi_handle[:,:] = densi
    
    tempr = temp.real 
    tempr_handle=f.createVariable('tempr', tempr.dtype, ('nzed','nmodes',))
    tempr_handle[:,:] = tempr
    
    tempi = temp.imag 
    tempi_handle=f.createVariable('tempi', tempi.dtype, ('nzed','nmodes',))
    tempi_handle[:,:] = tempi
    
    currentr = current.real 
    currentr_handle=f.createVariable('currentr', currentr.dtype, ('nmodes',))
    currentr_handle[:] = currentr
    
    currenti = current.imag 
    currenti_handle=f.createVariable('currenti', currenti.dtype, ('nmodes',))
    currenti_handle[:] = currenti
    
    omegar = omega.real 
    omegar_handle=f.createVariable('omegar', omegar.dtype, ('nmodes',))
    omegar_handle[:] = omegar
    
    omegai = omega.imag 
    omegai_handle=f.createVariable('omegai', omegai.dtype, ('nmodes',))
    omegai_handle[:] = omegai
    
    wstarn_handle=f.createVariable('wstarn', float,())
    wstarn_handle[...] = wstarn
    wstart_handle=f.createVariable('wstart', float,())
    wstart_handle[...] = wstart
    tite_handle=f.createVariable('tite', float,())
    tite_handle[...] = tite
    shat_handle=f.createVariable('shat', float,())
    shat_handle[...] = shat
    beta_handle=f.createVariable('beta', float,())
    beta_handle[...] = beta
    a0_handle=f.createVariable('a0', float,())
    a0_handle[...] = a0
    a1_handle=f.createVariable('a1', float,())
    a1_handle[...] = a1
    a2_handle=f.createVariable('a2',float,())
    a2_handle[...] = a2
    d0_handle=f.createVariable('d0', float,())
    d0_handle[...] = d0
    d1_handle=f.createVariable('d1', float,())
    d1_handle[...] = d1
    d2_handle=f.createVariable('d2', float,())
    d2_handle[...] = d2
    f0_handle=f.createVariable('f0', float,())
    f0_handle[...] = f0
    wdrift_handle=f.createVariable('wdrift', float,())
    wdrift_handle[...] = wdrift
    b0_handle=f.createVariable('b0', float,())
    b0_handle[...] = b0
    b1_handle=f.createVariable('b1', float,())
    b1_handle[...] = b1
    b2_handle=f.createVariable('b2', float,())
    b2_handle[...] = b2
    b3_handle=f.createVariable('b3', float,())
    b3_handle[...] = b3
    m0_handle=f.createVariable('m0', float,())
    m0_handle[...] = m0
    m1_handle=f.createVariable('m1', float,())
    m1_handle[...] = m1
    m2_handle=f.createVariable('m2', float,())
    m2_handle[...] = m2
    m3_handle=f.createVariable('m3', float,())
    m3_handle[...] = m3
    
    f.close()

def plot_eigfns(data,work_dir,run_name):

    [zed,dens_all,temp_all,current_all,omega_all,pin] = data
    
    [nzed,nmodes]= np.shape(dens_all)
    
    file = work_dir+"/"+run_name + ".pdf"
    pdf = PdfPages(file)
    for imode in range(0,nmodes):
        dens = dens_all[:,imode]
        temp = temp_all[:,imode]
        current = current_all[imode]
        omega = omega_all[imode]
        
        field_list = [dens.real,dens.imag,
         temp.real,temp.imag]
        xlist = [zed for ifield in field_list]

        marker_list = ["b","r","g","k","b-.","r-.","g-.","k-."]
        xlab = "$z$"
        title = "$ \\omega$ = " +str(round_sig(omega.real)) + "+" +str(round_sig(omega.imag)) + "i" 
        legend_title = "$J = $" + str(round_sig(current.real)) + "+" +str(round_sig(current.imag)) + "i" 
        title = title + " \n " + legend_title 
        ylab_list = ["$n_{r}$","$n_{i}$","$T_{r}$","$T_{i}$"]
        plot_1d_list_pdf (xlist,field_list,marker_list,xlab, pdf,
          title=title,ylab='',xlims=None,ylims=None,aspx=9,aspy=6, xticks = None, yticks = None,
          markersize=10, legend_title="", use_legend=True,loc_opt='upper right', ylab_list = ylab_list,
          bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=25, ncol_opt=2)
        plot_1d_list_pdf (xlist,field_list,marker_list,xlab, pdf,
          title=title,ylab='',xlims=[-0.2,0.2],ylims=None,aspx=9,aspy=6, xticks = None, yticks = None,
          markersize=10, legend_title="", use_legend=True,loc_opt='upper right', ylab_list = ylab_list,
          bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=25, ncol_opt=2)
    
    pdf.close()     
    print(file)
    