# -*- coding: utf-8 -*-
import numpy as np


def eigsolve(MM,zed,nzh,electromagnetic,pin):

    soln = np.linalg.eig(MM)
    
    E0 = soln[0] # the eigenvalues (E0 = -i omega = d/dt)
    indices = np.argsort(E0.real) # sort by growth rate 
    
    E = soln[0][indices] #get the eigenvalues into order
    eigfns = soln[1][:,indices] #get the eigenvectors into order
    
    # determine number of unstable modes to write output
    # or write out the 5 slowest decaying modes 
    nmodes = E.size
    nunstable = 0
    for imode in range(0,nmodes):
        if (E[imode].real > 0):
            omega = 1j*E[imode]
            #print("mode ",imode, "-i omega = ",E[imode],"\n")
            #print("mode ",imode, " omega = ",omega,"\n")
            nunstable = nunstable + 1
    
    if(nunstable == 0):
        ndiag = 5
    else:
        ndiag = nunstable
        
    omega = np.zeros(ndiag,dtype=complex)
    for imode in range(0,ndiag):
        omega[imode] = 1j*E[imode+nmodes-ndiag]
        print("mode ",imode, " omega = ",omega[imode],"\n")
    
    if(electromagnetic):
        dens, temp, current, zed_out = get_electromagnetic_eigfns(eigfns[:,nmodes-ndiag:],zed,nzh,pin,omega)
    else:
        dens, temp, current, zed_out = get_electrostatic_eigfns(eigfns[:,nmodes-ndiag:],zed,nzh,pin)
    
    return omega, dens, temp, current, zed_out 
    
def get_electromagnetic_eigfns(eigfns,zed,nzh,pin,omega,test=False):
    nzed = zed.size 
    [wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,wdrift,b0,b1,b2,b3,m0,m1,m2,m3] = pin
    [nmatrix,nmodes] = np.shape(eigfns)
    ized_current = nmatrix - 1
    ized_dzp = nzh
    ized_dzm = nzh - 1
    
    Dzed = zed[ized_dzp]
    #print(Dzed)
    zed_full = np.zeros([nzed+2])
    zed_full[:nzh] = zed[:nzh]
    zed_full[nzh+2:] = zed[nzh:]
    #print("zed_full",zed_full)
    
    dens_full = np.zeros( [nzed+2,nmodes],dtype=complex)
    temp_full = np.zeros( [nzed+2,nmodes],dtype=complex)
    current = np.zeros(nmodes,dtype=complex)
    for imode in range(0,nmodes):
        dens = eigfns[0:nzed,imode] 
        imax = np.argmax(np.abs(dens))
        norm = dens[imax]
        dens = dens/norm
        temp = eigfns[nzed:2*nzed,imode]/norm
        current[imode] = eigfns[ized_current,imode]/norm
        
        # postprocessing using imposed bc.
        dens0p = ( -(2./3.)*(Dzed/a0)*current[imode] + (4./3.)*dens[ized_dzp] - (1./3.)*dens[ized_dzp+1] 
        + (a1/a0)*( (2./3.)*(temp[ized_dzp]-temp[ized_dzm]) + (1./6.)*(temp[ized_dzm-1]-temp[ized_dzp+1]) )  - (a1/a0)*( 1j*np.pi*beta*wstart*current[imode]/(2.*shat) ) )
        dens0m = -dens0p + (4./3.)*(dens[ized_dzp] + dens[ized_dzm]) - (1./3.)*(dens[ized_dzp+1] + dens[ized_dzm-1])
        temp0m = (2./3.)*(temp[ized_dzm] + temp[ized_dzp]) - (1./6.)*(temp[ized_dzm-1] + temp[ized_dzp+1]) - 1j*np.pi*beta*wstart*current[imode]/(2.*shat) 
        temp0p = (2./3.)*(temp[ized_dzm] + temp[ized_dzp]) - (1./6.)*(temp[ized_dzm-1] + temp[ized_dzp+1]) + 1j*np.pi*beta*wstart*current[imode]/(2.*shat) 
        
        if(test):
            #temp0prim = (4.*temp[ized_dzp] - 3.*temp0p - temp[ized_dzp+1])/(2.*Dzed)
            temp0prim = (4.*(temp[ized_dzp]-temp[ized_dzm]) - 3.*(temp0p-temp0m) - temp[ized_dzp+1] + temp[ized_dzm-1])/(4.*Dzed)
            #dens0prim = (4.*dens[ized_dzp] - 3.*dens0p - dens[ized_dzp+1])/(2.*Dzed)
            dens0prim = (4.*(dens[ized_dzp]-dens[ized_dzm]) - 3.*(dens0p-dens0m) - dens[ized_dzp+1] + dens[ized_dzm-1])/(4.*Dzed)
            current_post = a0*dens0prim + a1*temp0prim
            print("sanity check comparisons of neighbouring grid points -- not second order accurate! \n ")
            print("dens0prim",dens0prim,"(dens[ized_dzp+1]-dens[ized_dzp])/Dzed,",(dens[ized_dzp+1]-dens[ized_dzp])/Dzed,
             "(dens[ized_dzm]-dens[ized_dzm-1])/Dzed,",(dens[ized_dzm]-dens[ized_dzm-1])/Dzed)
            print("temp0prim",temp0prim,"(temp[ized_dzp+1]-temp[ized_dzp])/Dzed,",(temp[ized_dzp+1]-temp[ized_dzp])/Dzed,
             "(temp[ized_dzm]-temp[ized_dzm-1])/Dzed,",(temp[ized_dzm]-temp[ized_dzm-1])/Dzed)
            print("n(Dz)-n(-Dz)",dens[ized_dzp]-dens[ized_dzm],"n(0+)-n(0-)",dens0p-dens0m) 
            print("T(Dz)-T(-Dz)",temp[ized_dzp]-temp[ized_dzm],"T(0+)-T(0-)",temp0p-temp0m)
            print("\n")
            print("Check b.c. using second order accurate formulas \n")
            print("J",current[imode],"a0*dens0prim + a1*temp0prim",current_post)
            print("n(0+)-n(0-)",dens0p-dens0m, "1j*np.pi*beta*(wstarn-omega)*current/shat",1j*np.pi*beta*(wstarn-omega[imode])*current[imode]/shat)
            print("T(0+)-T(0-)",temp0p-temp0m, "1j*np.pi*beta*wstart*current/shat",1j*np.pi*beta*wstart*current[imode]/shat)
            print("\n")
        
        
        dens_full[:nzh,imode] = dens[:nzh] 
        dens_full[nzh,imode] = dens0m 
        dens_full[nzh+1,imode] = dens0p 
        dens_full[nzh+2:,imode] =  dens[nzh:]

        temp_full[:nzh,imode] = temp[:nzh] 
        temp_full[nzh,imode] = temp0m 
        temp_full[nzh+1,imode] = temp0p 
        temp_full[nzh+2:,imode] =  temp[nzh:]

    return dens_full, temp_full, current, zed_full
    
def get_electrostatic_eigfns(eigfns,zed,nzh,pin):
    nzed = zed.size 
    [wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,wdrift,b0,b1,b2,b3,m0,m1,m2,m3] = pin
    [nmatrix,nmodes] = np.shape(eigfns)
    
    
    ized_dzp = nzh + 1
    ized_dzm = nzh - 1
    #print(zed[ized_dzm])
    #print(zed[ized_dzp])
    #print(zed)
    Dzed = zed[ized_dzp]
    dens_full = np.zeros( [nzed,nmodes],dtype=complex)
    temp_full = np.zeros( [nzed,nmodes],dtype=complex)
    current = np.zeros(nmodes,dtype=complex)
    for imode in range(0,nmodes):
        
        dens = eigfns[0:nzed,imode] 
        imax = np.argmax(np.abs(dens))
        norm = dens[imax]
        dens_full[:,imode] = dens/norm
        temp_full[:,imode] = eigfns[nzed:2*nzed,imode]/norm
        current[imode] = (a0*(dens_full[ized_dzp,imode]-dens_full[ized_dzm,imode])+ a1*(temp_full[ized_dzp,imode]-temp_full[ized_dzm,imode]) )/(2.*Dzed)
    
    return dens_full, temp_full, current, zed
    
