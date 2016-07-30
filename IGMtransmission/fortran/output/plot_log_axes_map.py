from pylab import *
from scipy.interpolate import interp2d

# read from files
nfiles=1
s_perp=range(nfiles)
s_para=range(nfiles)
xi_s=range(nfiles)

s_para,s_perp,tau_eff=genfromtxt('RSD_01.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_02.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_16.dat',usecols=(4,5,6),unpack=True,comments='#')

N=100

F=exp(-tau_eff).reshape(N,N) # F[s_para,s_perp]
#F=1.-F
s_perp=s_perp.reshape(N,N)
s_para=s_para.reshape(N,N)

Func=interp2d(s_para,s_perp,F,kind='linear')

Ng=100
x=np.logspace(log10(0.08), log10(3.0), Ng)
y=np.logspace(log10(0.08), log10(3.0), Ng)
s_para_mesh,s_perp_mesh=np.meshgrid(x,y)
F_mesh=zeros([Ng,Ng])
for i in range(Ng):
    for j in range(Ng):
        F_mesh[i,j]=Func(s_para_mesh[i,j],s_perp_mesh[i,j])

F4=zeros([2*N,2*N])
F4[N:2*N,N:2*N]=F[:,:]
F4[0:N,N:2*N]=F[::-1,:]
F4[N:2*N,0:N]=F[:,::-1]
F4[0:N,0:N]=F[::-1,::-1]

# plot
fig=plt.figure(figsize=(6.5,5))
plt.rc('font',**{'family':'serif', 'size':16})
#plt.subplots_adjust(left=0.0, right=0.92, top=0.95, bottom=0.1)
plt.subplots_adjust(left=0.0, right=0.92, top=0.95, bottom=0.15)

extent=[0.08, 2.0, 0.08, 3.0]
#im=plt.pcolormesh(s_perp,s_para,1-F,cmap='RdYlBu_r',vmin=0.3,vmax=0.75)
im=plt.pcolormesh(s_para_mesh,s_perp_mesh,1.0-F_mesh.T,cmap='RdYlBu_r',vmin=0.3,vmax=0.75)
ylabel('$s_{\\parallel}$ $[\\rm pMpc]$',fontsize=20)
xlabel('$s_{\\perp}$ $[\\rm pMpc]$',fontsize=20)
plt.xscale('log')
plt.yscale('log')
plt.axis(extent)
plt.gca().set_aspect('equal')
gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
levels=arange(0.3,0.75,0.025)
cs=plt.contour(s_para_mesh,s_perp_mesh,1-F_mesh.T,levels,colors='black',alpha=0.5)
#colorbar()
cax = fig.add_axes([0.8, 0.2, 0.03, 0.7])
cbar=fig.colorbar(im,cax=cax)
cbar.add_lines(cs)
cbar.set_label('$1-\\langle F_{\\alpha}(s_{\\parallel},s_{\\perp})\\rangle$',fontsize=20)
#plt.tight_layout()


# plot in linear scale
fig=plt.figure(figsize=(6.5,5))
plt.rc('font',**{'family':'serif', 'size':16})
#plt.subplots_adjust(left=0.0, right=0.92, top=0.95, bottom=0.1)
plt.subplots_adjust(left=0.0, right=0.92, top=0.95, bottom=0.15)

extent=[-s_perp.max(),s_perp.max(),-s_para.max(),s_para.max()]
##extent=[-2.5,2.5,-2.5,2.5]

#extent=[0,2.5,0,2.5]

im=imshow(1.0-F4,cmap='RdYlBu_r',extent=extent,vmin=0.3,vmax=0.75)
ylabel('$s_{\\parallel}$ $[\\rm pMpc]$',fontsize=20)
xlabel('$s_{\\perp}$ $[\\rm pMpc]$',fontsize=20)
plt.gca().set_aspect('equal')
xlim(-2.5,2.5)
ylim(-2.5,2.5)
levels=arange(0.3,0.75,0.025)
cs=plt.contour(1-F4,levels,extent=extent,colors='black',linestyles='solid',alpha=0.5)
#cs=plt.contour(1-F4,20,extent=extent,colors='black',linestyles='solid',alpha=0.5)
cax = fig.add_axes([0.8, 0.2, 0.03, 0.7])
cbar=fig.colorbar(im,cax=cax)
cbar.add_lines(cs)
cbar.set_label('$1-\\langle F_{\\alpha}(s_{\\parallel},s_{\\perp})\\rangle$',fontsize=20)

# plot as countour only map
fig=plt.figure(figsize=(6.5,5))
plt.rc('font',**{'family':'serif', 'size':16})
#plt.subplots_adjust(left=0.0, right=0.92, top=0.95, bottom=0.1)
plt.subplots_adjust(left=0.0, right=0.92, top=0.95, bottom=0.15)

extent=[-s_perp.max(),s_perp.max(),-s_para.max(),s_para.max()]
##extent=[-2.5,2.5,-2.5,2.5]

#extent=[0,2.5,0,2.5]

#im=imshow(1.0-F4,cmap='RdYlBu_r',extent=extent,vmin=0.3,vmax=0.90)
ylabel('$s_{\\parallel}$ $[\\rm pMpc]$',fontsize=20)
xlabel('$s_{\\perp}$ $[\\rm pMpc]$',fontsize=20)
plt.gca().set_aspect('equal')
xlim(-2.5,2.5)
ylim(-2.5,2.5)
levels=arange(0.3,0.95,0.025)
cs=plt.contour(1-F4,levels,extent=extent,linestyles='solid',alpha=0.5)
#cs=plt.contour(1-F4,20,extent=extent,colors='black',linestyles='solid',alpha=0.5)
cax = fig.add_axes([0.8, 0.2, 0.03, 0.7])
cbar=fig.colorbar(cs,cax=cax)
cbar.add_lines(cs)
cbar.set_label('$1-\\langle F_{\\alpha}(s_{\\parallel},s_{\\perp})\\rangle$',fontsize=20)

