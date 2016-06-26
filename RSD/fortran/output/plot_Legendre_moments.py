from pylab import *
from scipy.interpolate import interp2d

# read from files
nfiles=1
s_perp=range(nfiles)
s_para=range(nfiles)
xi_s=range(nfiles)

#s_para,s_perp,tau_eff=genfromtxt('RSD_00.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_02.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_03.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_04.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_05.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_06.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_07.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_08.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_09.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_10.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_11.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_12.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_13.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_14.dat',usecols=(4,5,6),unpack=True,comments='#')
#s_para,s_perp,tau_eff=genfromtxt('RSD_15.dat',usecols=(4,5,6),unpack=True,comments='#')
s_para,s_perp,tau_eff=genfromtxt('RSD_16.dat',usecols=(4,5,6),unpack=True,comments='#')
N=100

F=exp(-tau_eff).reshape(N,N) # F[s_para,s_perp]
#F=1.-F
s_perp=s_perp.reshape(N,N)
s_para=s_para.reshape(N,N)

F4=zeros([2*N,2*N])
F4[N:2*N,N:2*N]=F[:,:]
F4[0:N,N:2*N]=F[::-1,:]
F4[N:2*N,0:N]=F[:,::-1]
F4[0:N,0:N]=F[::-1,::-1]
figure()
extent=[-s_perp.max(),s_perp.max(),-s_para.max(),s_para.max()]
imshow(F4,cmap='RdYlBu',extent=extent)
#imshow(F4,cmap='Spectral',extent=extent)
colorbar()
contour(F4,extent=extent,colors='black',alpha=0.5,levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.675])

figure()
extent=[-s_perp.max(),s_perp.max(),-s_para.max(),s_para.max()]
imshow(log10(1-F4/F4.max()),cmap='RdYlBu_r',extent=extent,vmin=-2.0,vmax=0.0)
colorbar()
levels=linspace(-0.40,0.45,10)
contour(log10(1-F4/F4.max()),20,extent=extent,colors='black',linestyles='solid',alpha=0.5)



s_perp_grid=s_perp[0,:]
s_para_grid=s_para[:,0]

# Legendre moment of redshift space correlation function
def xi_moment(s,l):   

    def find_nearest(array,value):
        idx = (abs(array-value)).argmin()
        return array[idx],idx

    def Legendre(x,l):
        if l == 0:
            Legendre=1.
        elif l == 2:
            Legendre=1./2.*(3.*x*x-1.)
        elif l == 4:
            Legendre=1./8.*(35.*x**4-30.*x**2+3.)
        else:
            print 'no Legendre moment implemented yet'
        return Legendre

    Ng=500
    mu,dmu=linspace(0,1,Ng,retstep=True)

    xi_moment=0.
    for i in range(Ng):
        s_para_i=s*mu[i]    
        s_perp_i=s*sqrt(1.-mu[i]**2)
        s_para_idx=find_nearest(s_para_grid,s_para_i)[1]
        s_perp_idx=find_nearest(s_perp_grid,s_perp_i)[1]
        xi_moment+=F[s_para_idx,s_perp_idx]*Legendre(mu[i],l)*dmu
    xi_moment=(2.*l+1.)/2. * (xi_moment * 2.) 
    # x2 because we only integrated over the first quandrant mu=(0,1) 

    return xi_moment

# main
Nbins=30
s_min=0.01 # pMpc
s_max=3.0 # pMpc
s=linspace(s_min,s_max,Nbins)

xi_0=zeros(Nbins)
xi_2=zeros(Nbins)
xi_4=zeros(Nbins)
for i in range(Nbins):
    xi_0[i]=xi_moment(s[i],0)
    xi_2[i]=xi_moment(s[i],2)
    xi_4[i]=xi_moment(s[i],4)

def tau_eff_B13(z):
    # Becker+13 formula
    tau_eff_B13=0.751*((1.+z)/4.5)**2.90 - 0.132
    return tau_eff_B13


# load Adelberger+05
A05=genfromtxt('A05.dat',comments='#')
F_05=0.765

VLRS=genfromtxt('VLRS.dat',comments='#')


# plot
fig=plt.figure(figsize=(15,6))

figure(figsize=(15,6))
plt.rc('font',**{'family':'serif', 'size':22})
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('legend', fontsize=23)

subplot(1,3,1)
ylabel('$\\xi_0(s)$',fontsize=22)
xlabel('$s$ $[\\rm pMpc]$',fontsize=22)
plot(s,xi_0,'r-', linewidth=2,label='fiducial')
#errorbar(A05[:,0]/(1.+3.)/0.7,A05[:,1],yerr=A05[:,2]-A05[:,1],fmt='x',color='gray',label='Adelberger+05')
errorbar(VLRS[:,0]/(1.+3.)/0.7,VLRS[:,1],yerr=VLRS[:,2]-VLRS[:,1],fmt='o',color='black',label='VLRS+KBSS')
#hlines(F_05,0,s_max,linestyles='dotted')
hlines(exp(-tau_eff_B13(3.0)),0,s_max,linestyles='dotted')
xlim(0,s_max)
ylim(0,1.0)
lg=plt.legend(loc='upper right',fontsize=14)
lg.draw_frame(False)

subplot(1,3,2)
#ylabel('$s^2\\xi_2(s)$',fontsize=22)
ylabel('$\\xi_2(s),$ $\\xi_4(s)$',fontsize=22)
xlabel('$s$ $[\\rm pMpc]$',fontsize=22)
#plot(s,s**2*xi_2,'b-', linewidth=2,label='fiducial')
#plot(s,s**2*xi_4,'g-', linewidth=2,label='fiducial')
plot(s[1:],xi_2[1:],'b-', linewidth=2,label='$\\xi_2$')
plot(s[1:],-xi_2[1:]/(1-xi_0[1:]),'g-', linewidth=2,label='$-\\xi_2/(1-\\xi_0)$')
#plot(s,xi_4,'g-', linewidth=2,label='$\\xi_4$')
lg=plt.legend(loc='upper right',fontsize=14)
lg.draw_frame(False)

subplot(1,3,3)
extent=[-s_perp.max(),s_perp.max(),-s_para.max(),s_para.max()]
im=imshow(log10(1-F4/F4.max()),cmap='RdYlBu_r',extent=extent,vmin=-2.0,vmax=0.0)
#colorbar()
levels=linspace(-2.0,0.0,11)
contour(log10(1-F4/F4.max()),20,extent=extent,colors='black',linestyles='solid',alpha=0.5,levels=levels)
#cax = fig.add_axes([0.8, 0.2, 0.03, 0.7])
cbar=fig.colorbar(im)#, cax=cax)
cbar.set_label('$\\log_{10} 1-\\langle F_{\\alpha}\\rangle/\\bar{F}$',fontsize=20)

fig=plt.figure(figsize=(8,5))
extent=[-s_perp.max(),s_perp.max(),-s_para.max(),s_para.max()]
im=imshow(log10(1-F4/F4.max()),cmap='RdYlBu_r',extent=extent,vmin=-2.0,vmax=0.0)
#colorbar()
levels=linspace(-2.0,0.0,11)
contour(log10(1-F4/F4.max()),20,extent=extent,colors='black',linestyles='solid',alpha=0.5,levels=levels)
cax = fig.add_axes([0.8, 0.2, 0.03, 0.7])
cbar=fig.colorbar(im, cax=cax)
cbar.set_label('$\\log_{10} 1-\\langle F_{\\alpha}\\rangle/\\bar{F}$',fontsize=20)

#imshow(F4,cmap='RdYlBu')
#colorbar()
#contour(F4,20,colors='black',alpha=0.5)

#ylabel('$\\xi_2/\\xi_0$',fontsize=22)
#xlabel('$s$ $[\\rm pMpc]$',fontsize=22)
#plot(s,xi_2/xi_0,'g-', linewidth=2,label='fiducial')
#lg=plt.legend(loc='upper right',fontsize=14)
#lg.draw_frame(False)

plt.tight_layout()

#figure()
#imshow(xi_s_map)
