# -*- coding: utf-8 -*-

from sct_post_processing import sctdata
from sct_diagnostics import plot_eigfns
from os import getcwd 
work_dir = getcwd() + "/runs/" # path where simulation files are placed
run_name = "sct.test.cbc"
out_dir = work_dir # path where output files will be placed
data0 = sctdata(work_dir,run_name,out_dir)
#print(data0.zed)

zed = data0.zed
dens = data0.dens
temp = data0.temp
current = data0.current
omega = data0.omega
physics_input = data0.physics_input

data = [zed,dens,temp,current,omega,physics_input]
plot_eigfns(data,work_dir,run_name)