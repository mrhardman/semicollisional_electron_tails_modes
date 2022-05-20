# -*- coding: utf-8 -*-
import numpy as np
#import f90nml
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as colors
# setup some plot defaults
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=30)
rcParams.update({'text.latex.preamble' : r'\usepackage{bm}'})
rcParams.update({'figure.autolayout': True})

def round_sig(x, sig=3):
    import numpy as np
    from math import log10, floor
    if(np.isnan(x)): # this assumes non-array input just a single value
        return np.nan
    elif(abs(x)<1.0e-10 ):
        return 0.0
    else:
        return np.round(x, sig-int(floor(log10(abs(x))))-1)

def lin_func(x,a,b):
    return np.add(np.multiply(x,a),b)
def lin_func2(x,a):
    return np.multiply(x,a)
def quad_func(x,a,b,c):
    y = np.add(np.multiply(x,b),c)
    z = np.add(np.multiply(x,x)*a,y)
    return z
def quad_func2(x,a,b):
    y = np.multiply(x,b)
    z = np.add(np.multiply(x,x)*a,y)
    return z
def linsqrt_func(x,a,b,c):
    x = np.sqrt(x)
    y = np.add(np.multiply(x,b),c)
    z = np.add(np.multiply(x,x)*a,y)
    return z

def power_func(x,a,b,p):
    return np.add(np.multiply(np.power(x,p),a),b)

def lincos_func(x,a,b):
    return np.add(np.multiply(x,a),np.multiply(np.cos(x),b))
def lincossin_func(x,a,b,c):
    z1  = np.add(np.multiply(x,a),np.multiply(np.cos(x),b))
    z2 = np.add(z1,np.multiply(np.sin(x),c))
    return np.add(np.multiply(x,a),np.multiply(np.cos(x),b))
def lincosshift_func(x,a,b,c,d):
    z  = np.add(np.multiply(x,a),np.multiply(np.cos(np.add(x,c)),b))
    z = np.add(z,d)
    return z
def lincosshift2_func(x,a,b,c,d,e,f):
    z  = np.add(np.multiply(x,a),np.multiply(np.cos(np.add(x,c)),b))
    z  = np.add(z,np.multiply(np.cos(np.add(2.*x,e)),d))
    z = np.add(z,f)
    return z   
def log_func(x,m,n,c,shift):
    z = np.multiply(m,np.log(x))
    xprim = np.subtract(x,shift)
    z = np.add(z,np.multiply(n,np.log(xprim)))
    z = np.add(z,c)
    return z 

def plot_1d_pdf (x,y,xlab, pdf, title='',ylab='',end=False,zeros=False,marker='k',loglog=False,
    vlines = None,use_legend=False,legend_title="", loc_opt='upper right',
  bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=30, ncol_opt=1,
  legend_shadow=False,legend_frame=False):

    fig=plt.figure(figsize=(12,8))
    if(loglog):
        plt.loglog(x,y,marker,label=None)
    else:
        plt.plot(x,y,marker,label=None)
    if end:
        plt.plot(x,y[-1]*np.ones(len(x)),'g--')
    if zeros:
        plt.plot(x,np.zeros(len(x)),'r--')
    if(not vlines == None):
        for xin,xlabel,xcolor,xlinestyle in vlines:
            plt.axvline(x=xin, label=xlabel, color=xcolor,linestyle=xlinestyle,linewidth=3)
    plt.xlabel(xlab)
    if len(ylab) > 0:
        plt.ylabel(ylab)
    if len(title) > 0:
        plt.title(title)
    if(use_legend):
        plt.legend(title=legend_title,loc=loc_opt, bbox_to_anchor=bbox_to_anchor_opt,
        fontsize=legend_fontsize, frameon=legend_frame, handlelength=1, labelspacing=0.5,
        ncol=ncol_opt, columnspacing = 0.5 , handletextpad = 0.5, shadow=legend_shadow)
    
    pdf.savefig(fig)# pdf is the object of the current open PDF file to which the figures are appended
    plt.close (fig)
    return

def plot_1d_list_pdf (xlist,ylist,marker_list,xlab, pdf,
  title='',ylab='',xlims=None,ylims=None,aspx=12,aspy=8, xticks = None, yticks = None,
  markersize=5, legend_title="", use_legend=False,loc_opt='upper right', ylab_list = None,
  bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=10, ncol_opt=1,
  legend_shadow=False,legend_frame=False, vlines = None,marker_fill_style = None):

    fig=plt.figure(figsize=(aspx,aspy))
    nlist = len(ylist)
    if(ylab_list is None):
        ylab_list = [None for i in range(0,nlist)]
    for iy in range(0,nlist):
        plt.plot(xlist[iy],ylist[iy],marker_list[iy],markersize=markersize,label=ylab_list[iy],
        fillstyle = marker_fill_style)
    plt.xlabel(xlab)
    if len(ylab) > 0:
        plt.ylabel(ylab)
    if len(title) > 0:
        plt.title(title)
    if(not xlims is None):
        plt.xlim(xlims[0],xlims[1])
    if(not ylims is None):
        plt.ylim(ylims[0],ylims[1])
    if(not vlines == None):
        for xin,xlabel,xcolor,xlinestyle in vlines:
            plt.axvline(x=xin, label=xlabel, color=xcolor,linestyle=xlinestyle,linewidth=3)   
    if(use_legend):
        plt.legend(title=legend_title,loc=loc_opt, bbox_to_anchor=bbox_to_anchor_opt,
        fontsize=legend_fontsize, frameon=legend_frame, handlelength=1, labelspacing=0.5,
        ncol=ncol_opt, columnspacing = 0.5 , handletextpad = 0.5, shadow=legend_shadow)
    if(not xticks is None):
        plt.xticks(xticks)
    if(not yticks is None):
        plt.yticks(yticks)    
    pdf.savefig(fig)# pdf is the object of the current open PDF file to which the figures are appended
    plt.close (fig)
    return

def plot_1d_loglog_list_pdf (xlist,ylist,marker_list,xlab, pdf,
  title='',ylab='',xlims=None,ylims=None,aspx=12,aspy=8, xticks = None, yticks = None,
  markersize=5, legend_title="", use_legend=False,loc_opt='upper right', ylab_list = None,
  bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=10, ncol_opt=1,
  legend_shadow=False,legend_frame=False):

    fig=plt.figure(figsize=(aspx,aspy))
    nlist = len(ylist)
    if(ylab_list is None):
        ylab_list = [None for i in range(0,nlist)]
    for iy in range(0,nlist):
        plt.loglog(xlist[iy],ylist[iy],marker_list[iy],markersize=markersize,label=ylab_list[iy])
    plt.xlabel(xlab)
    if len(ylab) > 0:
        plt.ylabel(ylab)
    if len(title) > 0:
        plt.title(title)
    if(not xlims is None):
        plt.xlim(xlims[0],xlims[1])
    if(not ylims is None):
        plt.ylim(ylims[0],ylims[1])
    if(use_legend):
        plt.legend(title=legend_title,loc=loc_opt, bbox_to_anchor=bbox_to_anchor_opt,
        fontsize=legend_fontsize, frameon=legend_frame, handlelength=1, labelspacing=0.5,
        ncol=ncol_opt, columnspacing = 0.5 , handletextpad = 0.5, shadow=legend_shadow)
    if(not xticks is None):
        plt.xticks(xticks)
    if(not yticks is None):
        plt.yticks(yticks)    
    pdf.savefig(fig)# pdf is the object of the current open PDF file to which the figures are appended
    plt.close (fig)
    return

def plot_1d_semilog_list_pdf (xlist,ylist,marker_list,xlab, pdf,
  title='',ylab='',xlims=None,ylims=None,aspx=12,aspy=8, xticks = None, yticks = None,
  markersize=5, legend_title="", use_legend=False,loc_opt='upper right', ylab_list = None,
  bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=10, ncol_opt=1,
  legend_shadow=False,legend_frame=False):

    fig=plt.figure(figsize=(aspx,aspy))
    nlist = len(ylist)
    if(ylab_list is None):
        ylab_list = [None for i in range(0,nlist)]
    for iy in range(0,nlist):
        plt.semilogy(xlist[iy],ylist[iy],marker_list[iy],markersize=markersize,label=ylab_list[iy])
    plt.xlabel(xlab)
    if len(ylab) > 0:
        plt.ylabel(ylab)
    if len(title) > 0:
        plt.title(title)
    if(not xlims is None):
        plt.xlim(xlims[0],xlims[1])
    if(not ylims is None):
        plt.ylim(ylims[0],ylims[1])
    if(use_legend):
        plt.legend(title=legend_title,loc=loc_opt, bbox_to_anchor=bbox_to_anchor_opt,
        fontsize=legend_fontsize, frameon=legend_frame, handlelength=1, labelspacing=0.5,
        ncol=ncol_opt, columnspacing = 0.5 , handletextpad = 0.5, shadow=legend_shadow)
    if(not xticks is None):
        plt.xticks(xticks)
    if(not yticks is None):
        plt.yticks(yticks)    
    pdf.savefig(fig)# pdf is the object of the current open PDF file to which the figures are appended
    plt.close (fig)
    return

def plot_1d_semilogx_list_pdf (xlist,ylist,marker_list,xlab, pdf,
  title='',ylab='',xlims=None,ylims=None,aspx=12,aspy=8, xticks = None, yticks = None,
  markersize=5, legend_title="", use_legend=False,loc_opt='upper right', ylab_list = None,
  bbox_to_anchor_opt=(0.95, 0.95), legend_fontsize=10, ncol_opt=1,
  legend_shadow=False,legend_frame=False,xvals_vlines=None,
  columnspacing_opt = 2.0):

    fig=plt.figure(figsize=(aspx,aspy))
    nlist = len(ylist)
    if(ylab_list is None):
        ylab_list = [None for i in range(0,nlist)]
    for iy in range(0,nlist):
        plt.semilogx(xlist[iy],ylist[iy],marker_list[iy],markersize=markersize,label=ylab_list[iy])
    plt.xlabel(xlab)
    if len(ylab) > 0:
        plt.ylabel(ylab)
    if len(title) > 0:
        plt.title(title)
    if(not xlims is None):
        plt.xlim(xlims[0],xlims[1])
    if(not ylims is None):
        plt.ylim(ylims[0],ylims[1])
    if(not xvals_vlines is None):
        for xvals in xvals_vlines:
            xval = xvals[0]
            xval_lab = xvals[1]
            xval_col = xvals[2]
            xval_linestyle = xvals[3]
            plt.axvline(x=xval, label=xval_lab, color=xval_col,linestyle=xval_linestyle)
    if(use_legend):
        plt.legend(title=legend_title,loc=loc_opt, bbox_to_anchor=bbox_to_anchor_opt,
        fontsize=legend_fontsize, frameon=legend_frame, handlelength=1, labelspacing=0.5,
        ncol=ncol_opt, handletextpad = 0.5, shadow=legend_shadow,
        columnspacing = columnspacing_opt)
    if(not xticks is None):
        plt.xticks(xticks)
    if(not yticks is None):
        plt.yticks(yticks)    
    pdf.savefig(fig)# pdf is the object of the current open PDF file to which the figures are appended
    plt.close (fig)
    return
