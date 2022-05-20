# -*- coding: utf-8 -*-
from sct_main import evaluate_semicollisional_tails
from sct_input import sctphysicsinput
from sct_input import finite_aspect_ratio_coeffs
from os import getcwd
#input options
specification = "geofile" # "manual"
geowork = getcwd() #""
geofile = "geofiles/cbc" #cbc3 geo step_geo
geoout = "" #for diagnostic plots

#numerics input 
nzh = 50 # nzed/2 
Lz = 4. # length of half zed domain 

if(specification == "manual"):
    #physics input
    # normalised equilibrium quantities
    wstarn = -1.0
    wstart = -10.0
    tite = 1.0
    shat = 1.0
    beta = 1.0 # should be O(1) but large here for testing
    wdrift = -1.0 #should be 0 in model but included to get an ES instability
    # collisional coefficients
    a0 = 1.97
    a1 = 3.37
    a2 = 8.95
    d0 = 1.53
    d1 = -0.59 #-0.59 physical -- 0.59 unphysical sign for analytical eta >> 1 instability
    d2 = 1.92
    f0 = 1.462
    # aspect ratio small coeffs
    epsilon =  1./3.
    wdrift = wdrift*epsilon
    b0 = (epsilon**(1./4.))*2.43
    b1 = (epsilon**(1./4.))*0.69
    b2 = (epsilon**(1./4.))*0.69
    b3 = (epsilon**(1./4.))*4.54
    m0 = (epsilon**(1./4.))*0.01
    m1 = (epsilon**(1./4.))*1.87
    m2 = (epsilon**(1./4.))*0.01
    m3 = -(epsilon**(1./4.))*1.82
    physics_input = sctphysicsinput(wstarn,wstart,tite,shat,beta,a0,a1,a2,d0,d1,d2,f0,
 wdrift,b0,b1,b2,b3,m0,m1,m2,m3)

elif(specification == "geofile"):
    Zi = 1
    lpoly = 5 
    fprim = 0.773 #0.43
    tprim = 2.3 #2.77 #20.77
    theta0 = 0. #3.1415
    beta_effective = 100.0
    tite = 1.0
    physics_input = finite_aspect_ratio_coeffs(geowork,geofile,geoout,beta_effective,tite,Zi,fprim,tprim,theta0,lpoly)

work_dir = getcwd() + "/runs/"
run_name = "sct.test.cbc"

evaluate_semicollisional_tails(physics_input,nzh,Lz,work_dir,run_name)
