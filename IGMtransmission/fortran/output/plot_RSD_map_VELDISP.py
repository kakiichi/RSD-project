from pylab import *
from scipy.interpolate import interp2d
from matplotlib import gridspec


# read from files
nfiles=1
s_perp=range(nfiles)
s_para=range(nfiles)
xi_s=range(nfiles)

nfiles=3
title=range(nfiles)
title[0]='REF'
title[1]='HIGH DISP'
title[2]='LOW DISP'

vmin=0.3
vmax=0.90

gs = gridspec.GridSpec(1, nfiles, wspace=0.1, top=0.95, bottom=0.15,left=0.1, right = 0.85)
fig=plt.figure(figsize=(10.5,4))
plt.rc('font',**{'family':'serif', 'size':16})

for n in range(nfiles):
    if n==0:
        s_para,s_perp,tau_eff=genfromtxt('RSD_01.dat',usecols=(4,5,6),unpack=True,comments='#')
    if n==1:
        s_para,s_perp,tau_eff=genfromtxt('RSD_11.dat',usecols=(4,5,6),unpack=True,comments='#')
    if n==2:
        s_para,s_perp,tau_eff=genfromtxt('RSD_12.dat',usecols=(4,5,6),unpack=True,comments='#')
    N=100

    F_values=exp(-tau_eff)
    Ng=500
    s_para_mesh=linspace(s_para.min(),s_para.max(),Ng)
    s_perp_mesh=linspace(s_perp.min(),s_perp.max(),Ng)
    F_mesh=griddata(s_para, s_perp, F_values, s_para_mesh, s_perp_mesh, interp='linear')


    F=exp(-tau_eff).reshape(N,N) # F[s_para,s_perp]
    #F=1.-F
    s_perp=s_perp.reshape(N,N)
    s_para=s_para.reshape(N,N)

    # plot
    ax=subplot(gs[n])
    extent=[0.08, 2.0, 0.08, 2.8]
    #im=ax.pcolormesh(s_perp,s_para,1-F,cmap='RdYlBu_r',vmin=vmin,vmax=vmax)
    im=ax.pcolormesh(s_perp_mesh,s_para_mesh,1-F_mesh.T,cmap='RdYlBu_r',vmin=vmin,vmax=vmax)

    if n==0:
        ylabel('$s_{\\parallel}$ $[\\rm pMpc]$',fontsize=20)
    xlabel('$s_{\\perp}$ $[\\rm pMpc]$',fontsize=20)
    plt.xscale('log')
    plt.yscale('log')
    plt.title(title[n])
    plt.axis(extent)
    plt.gca().set_aspect('equal')
    gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    if n==1 or n==2: ax.yaxis.set_ticklabels([])
    levels=arange(vmin,vmax,0.025)
    cs=ax.contour(s_perp,s_para,1-F,levels,colors='black',alpha=0.5)
    cax = fig.add_axes([0.87, 0.2, 0.02, 0.7])
    cbar=fig.colorbar(im,cax=cax)
    cbar.add_lines(cs)
    cbar.set_label('$1-\\langle F_{\\alpha}(s_{\\parallel},s_{\\perp})\\rangle$',fontsize=20)

