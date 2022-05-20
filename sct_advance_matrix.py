# -*- coding: utf-8 -*-
import numpy as np

# write matrix to screen for debugging

def write_matrix(MM):
    [nzx,nzy] = np.shape(MM)
    for izedx in range(0,nzx):
        string = ""
        for izedy in range(0,nzy):
            string = string + str(round(MM[izedx,izedy].real,2))+"+"+ str(round(MM[izedx,izedy].imag,2))+"j, " 
        print(string)
        print("\n")

# diagonal matrix elements 
# from terms without derivatives 

def Qnn(pin,zed):
    
    [wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,wdrift,b0,b1,b2,b3,m0,m1,m2,m3] = pin
    element = (    -0.5*f0*d0*(1.+tite)*(shat*zed)**2 
                   + 1j*wstarn*(1. + tite)/(1. + 1./tite) 
                   - 1j*wdrift*(1. + tite) + 1j*shat*b0*(1. + tite)/2. )
    return element
    
def QnT(pin,zed):
    
    [wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,wdrift,b0,b1,b2,b3,m0,m1,m2,m3] = pin
    element = (    -0.5*f0*d1*(1.+tite)*(shat*zed)**2 
                   - 1j*wdrift*(1. + tite) + 1j*shat*b1*(1. + tite)/2. )
    return element
    
def QTn(pin,zed):
    
    [wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,wdrift,b0,b1,b2,b3,m0,m1,m2,m3] = pin
    element = (    -(1./3.)*f0*d1*(shat*zed)**2 
                   - 1j*(2./3.)*wdrift + 1j*wstart/(1. + 1./tite) + 1j*shat*b2/3. )
    return element
    
def QTT(pin,zed):
    
    [wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,wdrift,b0,b1,b2,b3,m0,m1,m2,m3] = pin
    element = (    -(1./3.)*f0*d2*(shat*zed)**2 
                   - 1j*(7./3.)*wdrift  + 1j*shat*b3/3. )
    return element

    # time advance matrix for electromagnetic case.
    # fluxinspired_densityeqn_midpoints does not support b0-m3 coefficients
def time_advance_matrix_electromagnetic(physics_input,nzh,Lz,fluxinspired_densityeqn_midpoints=False):
    # get physics parameters
    pin = physics_input
    [wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,wdrift,b0,b1,b2,b3,m0,m1,m2,m3] = pin
    # get zed grid
    nzed = 2*nzh 
    Dzed = Lz/nzh
    zzl = [Dzed * ( j  - nzh ) for j in range(0,nzh+1)]
    zzu = [Dzed * ( j  - nzh ) for j in range(nzh,2*nzh+1)]

    zed_incl0pm = np.array(zzl+zzu)
    zed = np.array(zzl[:-1]+zzu[1:])
    nmatrix = 2*nzed + 1 # no. dens + no. temp + 1 ( no. of J)

    ized_dzm = nzh -1        # index of zed = -Dzed
    ized_dzp = nzh           # index of zed = Dzed
    ized_current = nmatrix-1 # index of J

    MM = np.zeros([nmatrix,nmatrix],dtype=complex) 
    #print(np.shape(MM))  
      # note for interpretation of below assignments  
      # MM[i,j] is the element on the ith row in the jth column
      # equations are listed by row, and couple to fields by the column
      # so MM[i,j] is the coefficient of the ith eqn for the jth field. 
    
    #coefficients appearing on derivative terms 

    nndzsq = 0.5*a0*(1. + tite)/(Dzed**2) 
    nTdzsq = 0.5*a1*(1. + tite)/(Dzed**2)    
        
    #nndz = -1j*(1+tite)*shat*m0/(4.*Dzed)    
    #nTdz = -1j*(1+tite)*shat*m1/(4.*Dzed)
    nndzfac = -1j*(1+tite)*shat*m0/4.    
    nTdzfac = -1j*(1+tite)*shat*m1/4.
    nndz = nndzfac/Dzed
    nTdz = nTdzfac/Dzed

    Tndzsq = (1./3.)*a1/(Dzed**2) 
    TTdzsq = (1./3.)*a2/(Dzed**2)   
        
    #Tndz = -1j*shat*m2/(6.*Dzed)    
    #TTdz = -1j*shat*m3/(6.*Dzed)
    Tndzfac = -1j*shat*m2/6.    
    TTdzfac = -1j*shat*m3/6.
    Tndz = Tndzfac/Dzed    
    TTdz = TTdzfac/Dzed    
    
    # first density equation terms  
    #density eqn density coupling 
    ized = 0
    MM[ized,ized] = Qnn(pin,zed[ized]) - 2.*nndzsq
    MM[ized,ized+1] =                       nndzsq + nndz*zed[ized]

    #density eqn temperature coupling
    ized_temp = ized+nzed # temp column for z(ized)
    MM[ized,ized_temp] = QnT(pin,zed[ized]) - 2.*nTdzsq
    MM[ized,ized_temp+1] =                       nTdzsq + nTdz*zed[ized]

    for ized in range(1,nzed-1): 
        #density eqn density coupling 
        MM[ized,ized] = Qnn(pin,zed[ized]) - 2.*nndzsq
        MM[ized,ized+1] =                       nndzsq + nndz*zed[ized]
        MM[ized,ized-1] =                       nndzsq - nndz*zed[ized]
        
        #density eqn temperature coupling
        ized_temp = ized+nzed # temp column for z(ized)
        MM[ized,ized_temp] = QnT(pin,zed[ized]) - 2.*nTdzsq
        MM[ized,ized_temp+1] =                       nTdzsq + nTdz*zed[ized]
        MM[ized,ized_temp-1] =                       nTdzsq - nTdz*zed[ized]

    ized = nzed-1
    #density eqn density coupling 
    MM[ized,ized] = Qnn(pin,zed[ized]) - 2.*nndzsq
    MM[ized,ized-1] =                       nndzsq - nndz*zed[ized]
     #density eqn temperature coupling
    ized_temp = ized + nzed
    MM[ized,ized_temp] = QnT(pin,zed[ized]) - 2.*nTdzsq
    MM[ized,ized_temp-1] =                       nTdzsq - nTdz*zed[ized]

    #Correct midpoints
    if(fluxinspired_densityeqn_midpoints):
        # z = -D z
        ized = ized_dzm
        MM[ized,:] = 0.
        #density eqn density coupling 
        MM[ized,ized] = Qnn(pin,zed[ized]) - nndzsq/4.
        MM[ized,ized-2] = nndzsq/4.
        #density eqn temperature coupling
        ized_temp = ized+nzed # temp column for z(ized)
        MM[ized,ized_temp] = QnT(pin,zed[ized]) - nTdzsq/4.
        MM[ized,ized_temp-2] = nTdzsq/4.
        #density eqn current coupling
        MM[ized,ized_current] =  (1. + tite)/(4.*Dzed)


        # z = Dz
        ized = ized_dzp
        MM[ized,:] = 0.
        #density eqn density coupling 
        MM[ized,ized] = Qnn(pin,zed[ized]) - nndzsq/4.
        MM[ized,ized+2] = nndzsq/4.
        #density eqn temperature coupling
        ized_temp = ized+nzed # temp column for z(ized)
        MM[ized,ized_temp] = QnT(pin,zed[ized]) - nTdzsq/4.
        MM[ized,ized_temp+2] = nTdzsq/4.
        #density eqn current coupling
        MM[ized,ized_current] = - (1. + tite)/(4.*Dzed)
    else:
        # z = -D z
        ized = ized_dzm
        MM[ized,:] = 0.
        #density eqn density coupling 
        MM[ized,ized] = Qnn(pin,zed[ized]) - nndzsq*(2./3.) - (4./3.)*nndzfac  #-Dz
        MM[ized,ized-1] =                    nndzsq*(2./3.) + (4./3.)*nndzfac  #-2Dz
        #density eqn temperature coupling
        ized_temp = ized+nzed # temp column for z(ized)
        MM[ized,ized_temp-1] =                    nTdzsq*(2./3) + nndzfac*(1./6.)*(a1/a0) + nTdzfac*(1. + (1./6.))#-2Dz
        MM[ized,ized_temp] = QnT(pin,zed[ized]) - nTdzsq*(2./3) - nndzfac*(2./3.)*(a1/a0) - nTdzfac*(2./3.)       #-Dz
        MM[ized,ized_temp+1] =                                    nndzfac*(2./3.)*(a1/a0) - nTdzfac*(2./3.)       # Dz
        MM[ized,ized_temp+2] =                                   -nndzfac*(1./6.)*(a1/a0) + nTdzfac*(1./6.)       # 2Dz
        #density eqn current coupling
        MM[ized,ized_current] =  ( (1. + tite)/(3.*Dzed)
         - nndzfac*( (2./3.)*(Dzed/a0) + (1j*np.pi*beta*wstart/(2.*shat))*(a1/a0))
         + nTdzfac*(1j*np.pi*beta*wstart/(2.*shat)) )


        # z = Dz
        ized = ized_dzp
        MM[ized,:] = 0.
        #density eqn density coupling 
        MM[ized,ized] = Qnn(pin,zed[ized]) - nndzsq*(2./3.) - nndzfac*(4./3.)     # Dz
        MM[ized,ized+1] =                    nndzsq*(2./3.) + nndzfac*(4./3.)     # 2Dz
        #density eqn temperature coupling
        ized_temp = ized+nzed # temp column for z(ized)
        MM[ized,ized_temp-2] =                                   - nndzfac*(1./6.)*(a1/a0) + nTdzfac*(1./6.)  #-2Dz
        MM[ized,ized_temp-1] =                                     nndzfac*(2./3.)*(a1/a0) - nTdzfac*(2./3.)  #-Dz
        MM[ized,ized_temp] = QnT(pin,zed[ized]) - nTdzsq*(2./3.) - nndzfac*(2./3.)*(a1/a0) - nTdzfac*(2./3.)  # Dz
        MM[ized,ized_temp+1] =                    nTdzsq*(2./3.) + nndzfac*(1./6.)*(a1/a0) + nTdzfac*(1. +(1./6.))  # 2Dz
        #density eqn current coupling
        MM[ized,ized_current] = ( - (1. + tite)/(3.*Dzed)
         + nndzfac*( (2./3.)*(Dzed/a0) + (1j*np.pi*beta*wstart/(2.*shat))*(a1/a0) )
         - nTdzfac*(1j*np.pi*beta*wstart/(2.*shat))   )


    # second temperature equation terms    
    ized = 0
    #temperature eqn density coupling 
    MM[ized+nzed,ized] = QTn(pin,zed[ized]) - 2.*Tndzsq
    MM[ized+nzed,ized+1] =                       Tndzsq + Tndz*zed[ized]

    #temperature eqn temperature coupling
    MM[ized+nzed,ized+nzed] = QTT(pin,zed[ized]) - 2.*TTdzsq
    MM[ized+nzed,ized+nzed+1] =                       TTdzsq + TTdz*zed[ized]
    #(1./3.)*a2/(Dzed**2)
        
    for ized in range(1,nzed-1):
        #temperature eqn density coupling 
        MM[ized+nzed,ized] = QTn(pin,zed[ized]) - 2.*Tndzsq
        MM[ized+nzed,ized+1] =                       Tndzsq + Tndz*zed[ized]
        MM[ized+nzed,ized-1] =                       Tndzsq - Tndz*zed[ized]
        
        #temperature eqn temperature coupling
        MM[ized+nzed,ized+nzed] = QTT(pin,zed[ized]) - 2.*TTdzsq
        MM[ized+nzed,ized+nzed+1] =                       TTdzsq + TTdz*zed[ized]
        MM[ized+nzed,ized+nzed-1] =                       TTdzsq - TTdz*zed[ized]


    ized = nzed-1
    #temperature eqn density coupling 
    MM[ized+nzed,ized] = QTn(pin,zed[ized]) - 2.*Tndzsq
    MM[ized+nzed,ized-1] =                       Tndzsq - Tndz*zed[ized]

    #temperature eqn temperature coupling
    MM[ized+nzed,ized+nzed] =  QTT(pin,zed[ized]) - 2.*TTdzsq
    MM[ized+nzed,ized+nzed-1] =                        TTdzsq - TTdz*zed[ized]

    #correct midpoints
    ized = ized_dzm
    MM[ized+nzed,:] = 0.
    #temperature eqn density coupling 
    MM[ized+nzed,ized_dzm] = QTn(pin,zed[ized]) - (2./3.)*Tndzsq  - Tndzfac*(4./3.)
    MM[ized+nzed,ized_dzm-1] =                    (2./3.)*Tndzsq  + Tndzfac*(4./3.)

    #temperature eqn temperature coupling
    MM[ized+nzed,ized_dzm+nzed-1] =                     TTdzsq*(1. - (1./6.)) - Tndzsq*(1./6.)*(a1/a0) + Tndzfac*(1./6.)*(a1/a0) + TTdzfac*(1. + (1./6.))
    MM[ized+nzed,ized_dzm+nzed] =  QTT(pin,zed[ized]) + TTdzsq*((2./3.)- 2.)  + Tndzsq*(2./3.)*(a1/a0) - Tndzfac*(2./3.)*(a1/a0) - TTdzfac*(2./3.)
    MM[ized+nzed,ized_dzp+nzed] =                       TTdzsq*(2./3.)        - Tndzsq*(2./3.)*(a1/a0) + Tndzfac*(2./3.)*(a1/a0) - TTdzfac*(2./3.) 
    MM[ized+nzed,ized_dzp+nzed+1] =                   - TTdzsq*(1./6.)        + Tndzsq*(1./6.)*(a1/a0) - Tndzfac*(1./6.)*(a1/a0) + TTdzfac*(1./6.)

    #temperature eqn current coupling 
    MM[ized+nzed,ized_current] = ( -TTdzsq*1j*np.pi*beta*wstart/(2.*shat) 
                                  + Tndzsq*( (1j*np.pi*beta*wstart/(2.*shat))*(a1/a0) + (2./3.)*(Dzed/a0) ) 
                                  - Tndzfac*( (1j*np.pi*beta*wstart/(2.*shat))*(a1/a0) + (2./3.)*(Dzed/a0) )
                                  + TTdzfac*(1j*np.pi*beta*wstart/(2.*shat)) ) 

    ized = ized_dzp
    MM[ized+nzed,:] = 0.

    #temperature eqn density coupling 
    MM[ized+nzed,ized_dzp] = QTn(pin,zed[ized]) - Tndzsq*(2./3.) - Tndzfac*(4./3.)
    MM[ized+nzed,ized_dzp+1] =                    Tndzsq*(2./3.) + Tndzfac*(4./3.)

    #temperature eqn temperature coupling
    MM[ized+nzed,ized_dzp+nzed+1] =                     TTdzsq*(1. - (1./6.)) - Tndzsq*(1./6.)*(a1/a0) + Tndzfac*(1./6.)*(a1/a0) + TTdzfac*(1. + (1./6.))
    MM[ized+nzed,ized_dzp+nzed] =  QTT(pin,zed[ized]) + TTdzsq*((2./3.)- 2.)  + Tndzsq*(2./3.)*(a1/a0) - Tndzfac*(2./3.)*(a1/a0) - TTdzfac*(2./3.)
    MM[ized+nzed,ized_dzm+nzed] =                       TTdzsq*(2./3.)        - Tndzsq*(2./3.)*(a1/a0) + Tndzfac*(2./3.)*(a1/a0) - TTdzfac*(2./3.)
    MM[ized+nzed,ized_dzm+nzed-1] =                   - TTdzsq*(1./6.)        + Tndzsq*(1./6.)*(a1/a0) - Tndzfac*(1./6.)*(a1/a0) + TTdzfac*(1./6.)

    #temperature eqn current coupling 
    MM[ized+nzed,ized_current] = ( TTdzsq*1j*np.pi*beta*wstart/(2.*shat)
                                 - Tndzsq*( (1j*np.pi*beta*wstart/(2.*shat))*(a1/a0) + (2./3.)*(Dzed/a0) )
                                 + Tndzfac*( (1j*np.pi*beta*wstart/(2.*shat))*(a1/a0) + (2./3.)*(Dzed/a0) )
                                 - TTdzfac*(1j*np.pi*beta*wstart/(2.*shat)) )
        
    # current equation 

    # current eqn current coupling
    MM[ized_current,ized_current] = -1j*( wstarn + wstart*(a1/a0) )  - (4./3.)*(shat*Dzed)/(np.pi*beta*a0)
    # current eqn density coupling
    MM[ized_current,ized_dzp+1] = -(1./3.)*shat/(np.pi*beta)
    MM[ized_current,ized_dzp] = (4./3.)*shat/(np.pi*beta)
    MM[ized_current,ized_dzm] = -(4./3.)*shat/(np.pi*beta)
    MM[ized_current,ized_dzm-1] = (1./3.)*shat/(np.pi*beta)
    # current eqn temperature coupling
    MM[ized_current,nzed+ized_dzp+1] = -(1./3.)*shat*a1/(a0*np.pi*beta)
    MM[ized_current,nzed+ized_dzp] = (4./3.)*shat*a1/(a0*np.pi*beta)
    MM[ized_current,nzed+ized_dzm] = -(4./3.)*shat*a1/(a0*np.pi*beta)
    MM[ized_current,nzed+ized_dzm-1] = (1./3.)*shat*a1/(a0*np.pi*beta)

    if(nzed < 10):
        write_matrix(MM) #write out matrix for visualising

    return MM, zed

    # time advance matrix for electrostatic case.    
def time_advance_matrix_electrostatic(physics_input,nzh,Lz):
    
    # get the physics parameters
    pin = physics_input
    [wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,wdrift,b0,b1,b2,b3,m0,m1,m2,m3] = pin
    
    # the zed coordinate 
    nzed = 2*nzh + 1
    Dzed = Lz/nzh
    zz = [Dzed * ( j  - nzh ) for j in range(0,nzed)]
    zed = np.array(zz)

    #nmatrix = no. dens + no. temp = 2*nzed points

    MM = np.zeros([2*nzed,2*nzed],dtype=complex) 
    #print(np.shape(MM))  
      # note for interpretation of below assignments  
      # MM[i,j] is the element on the ith row in the jth column
      # equations are listed by row, and couple to fields by the column
      # so MM[i,j] is the coefficient of the ith eqn for the jth field. 
        
    #coefficients appearing on derivative terms 

    nndzsq = 0.5*a0*(1. + tite)/(Dzed**2) 
    nTdzsq = 0.5*a1*(1. + tite)/(Dzed**2)    
        
    nndz = -1j*(1+tite)*shat*m0/(4.*Dzed)    
    nTdz = -1j*(1+tite)*shat*m1/(4.*Dzed)

    Tndzsq = (1./3.)*a1/(Dzed**2) 
    TTdzsq = (1./3.)*a2/(Dzed**2)   
        
    Tndz = -1j*shat*m2/(6.*Dzed)    
    TTdz = -1j*shat*m3/(6.*Dzed)
        
    # first density equation terms  
    #density eqn density coupling 
    ized = 0
    MM[ized,ized] = Qnn(pin,zed[ized]) - 2.*nndzsq
    MM[ized,ized+1] = nndzsq + nndz*zed[ized]

    #density eqn temperature coupling
    MM[ized,ized+nzed] = QnT(pin,zed[ized]) - 2.*nTdzsq
    MM[ized,ized+nzed+1] = nTdzsq + nTdz*zed[ized]

    for ized in range(1,nzed-1): 
        #density eqn density coupling 
        MM[ized,ized] = Qnn(pin,zed[ized]) - 2.*nndzsq
        MM[ized,ized+1] = nndzsq + nndz*zed[ized]
        MM[ized,ized-1] = nndzsq - nndz*zed[ized]
        
        #density eqn temperature coupling
        MM[ized,ized+nzed] = QnT(pin,zed[ized]) - 2.*nTdzsq
        MM[ized,ized+nzed+1] = nTdzsq + nTdz*zed[ized]
        MM[ized,ized+nzed-1] = nTdzsq - nTdz*zed[ized]

    ized = nzed-1
    #density eqn density coupling 
    MM[ized,ized] = Qnn(pin,zed[ized]) - 2.*nndzsq
    MM[ized,ized-1] = nndzsq - nndz*zed[ized]
     #density eqn temperature coupling
    MM[ized,ized+nzed] = QnT(pin,zed[ized]) - 2.*nTdzsq
    MM[ized,ized+nzed-1] =  nTdzsq - nTdz*zed[ized]

        
    # second temperature equation terms    
    ized = 0
    #temperature eqn density coupling 
    MM[ized+nzed,ized] = QTn(pin,zed[ized]) - 2.*Tndzsq
    MM[ized+nzed,ized+1] = Tndzsq + Tndz*zed[ized]

    #temperature eqn temperature coupling
    MM[ized+nzed,ized+nzed] = QTT(pin,zed[ized]) - 2.*TTdzsq
    MM[ized+nzed,ized+nzed+1] = TTdzsq + TTdz*zed[ized]
    #(1./3.)*a2/(Dzed**2)
        
    for ized in range(1,nzed-1):
        #temperature eqn density coupling 
        MM[ized+nzed,ized] = QTn(pin,zed[ized]) - 2.*Tndzsq
        MM[ized+nzed,ized+1] = Tndzsq + Tndz*zed[ized]
        MM[ized+nzed,ized-1] = Tndzsq - Tndz*zed[ized]
        
        #temperature eqn temperature coupling
        MM[ized+nzed,ized+nzed] = QTT(pin,zed[ized]) - 2.*TTdzsq
        MM[ized+nzed,ized+nzed+1] = TTdzsq + TTdz*zed[ized]
        MM[ized+nzed,ized+nzed-1] = TTdzsq - TTdz*zed[ized]


    ized = nzed-1
    #temperature eqn density coupling 
    MM[ized+nzed,ized] = QTn(pin,zed[ized]) - 2.*Tndzsq
    MM[ized+nzed,ized-1] = Tndzsq - Tndz*zed[ized]

    #temperature eqn temperature coupling
    MM[ized+nzed,ized+nzed] =  QTT(pin,zed[ized]) - 2.*TTdzsq
    MM[ized+nzed,ized+nzed-1] = TTdzsq - TTdz*zed[ized]

    if(nzed < 10):
        write_matrix(MM) #write out matrix for visualising

    return MM, zed
