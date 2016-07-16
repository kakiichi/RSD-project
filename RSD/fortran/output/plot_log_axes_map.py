from pylab import *
from scipy.interpolate import interp2d

# read from files
nfiles=1
s_perp=range(nfiles)
s_para=range(nfiles)
xi_s=range(nfiles)

#s_para,s_perp,tau_eff=genfromtxt('RSD_16.dat',usecols=(4,5,6),unpack=True,comments='#')
s_para,s_perp,tau_eff=genfromtxt('RSD_16.dat',usecols=(4,5,6),unpack=True,comments='#')

N=100

F=exp(-tau_eff).reshape(N,N) # F[s_para,s_perp]
#F=1.-F
s_perp=s_perp.reshape(N,N)
s_para=s_para.reshape(N,N)

Func=interp2d(s_para,s_perp,F)

x=np.logspace(log10(0.08), log10(3.0), N)
y=np.logspace(log10(0.08), log10(3.0), N)
s_para_mesh,s_perp_mesh=np.meshgrid(x,y)
F_mesh=zeros([N,N])
for i in range(N):
    for j in range(N):
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
#plt.pcolormesh(s_perp,s_para,F,cmap='RdYlBu',vmin=0.2)
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


figure()
extent=[-s_perp.max(),s_perp.max(),-s_para.max(),s_para.max()]
optdpt_eff=-log10(F4)
imshow(log10(optdpt_eff),cmap='RdYlBu_r',extent=extent)
colorbar()
levels=arange(0.1,1.0,0.05).tolist()
contour(F4,extent=extent,colors='black',alpha=0.5,levels=levels)
#imshow(F4,cmap='Spectral',extent=extent)
#contour(F4,extent=extent,colors='black',alpha=0.5,levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.675])
