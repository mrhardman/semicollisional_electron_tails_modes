# -*- coding: utf-8 -*-

from sct_advance_matrix import time_advance_matrix_electromagnetic
from sct_advance_matrix import time_advance_matrix_electrostatic
from sct_eigsolve import eigsolve
#from sct_eigsolve import get_electromagnetic_eigfns
#from sct_eigsolve import get_electrostatic_eigfns
from sct_diagnostics import save_diagnostic_data
from sct_diagnostics import plot_eigfns

def evaluate_semicollisional_tails(physics_input,nzh,Lz,work_dir,run_name,electromagnetic=False):

    pin = physics_input.list
    beta = physics_input.beta
    #[wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,wdrift,b0,b1,b2,b3,m0,m1,m2,m3] = pin

    if(beta > 0):
        electromagnetic = True
    
    if(electromagnetic):
        MM, zed = time_advance_matrix_electromagnetic(pin,nzh,Lz)
    else:
        MM, zed = time_advance_matrix_electrostatic(pin,nzh,Lz)
    
    omega, dens, temp, current, zed = eigsolve(MM,zed,nzh,electromagnetic,pin)
    
    data = [zed,dens,temp,current,omega,pin] 
    save_diagnostic_data(work_dir,run_name,data)
    plot_eigfns(data,work_dir,run_name)
    
