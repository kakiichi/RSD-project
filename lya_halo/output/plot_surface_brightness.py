from scipy.integrate import quad,romb,romberg
from scipy.interpolate import interp1d

# parameters
z=3.0

Mpc2cm=3.086e+24              # 1 Mpc=3.086e24 cm
rad2arcsec=180./pi*60*60      # 1 radian to 1 arcsec

# load data
r_proper,emissivity=genfromtxt('lya_halo_16.dat',comments='#',usecols=(1,2),unpack=True)
emission=interp1d(r_proper,emissivity)

# Lya surface brightness [erg/s/cm^2/arcsec]
# b impact parameter [pMpc]
def SB_lya(b):
    b=b*Mpc2cm # cm
    r_max=max(r_proper)*Mpc2cm
    integrand=lambda r: emission(r/Mpc2cm)*r/sqrt(r**2-b**2)
    # erg/s/cm2/sr
    SB_lya=1./(2.*pi*(1.+z)**4)*quad(integrand,b,r_max)[0] 
    SB_lya=SB_lya / (rad2arcsec**2)  # erg/s/cm^2/arcsec
    return SB_lya

b=r_proper # pMpc
SB=zeros(b.size)
for i in range(b.size):
    SB[i]=SB_lya(b[i])
    print b[i],SB[i]

Momose=array([[1.924, 1.760e-17],
              [7.225, 6.351e-18],
              [10.69, 2.230e-18],
              [15.77, 8.267e-19],
              [21.51, 4.121e-19],
              [29.33, 1.450e-19],
              [33.90, 2.624e-19],
              [41.06, 3.149e-20],
              [47.90, 1.003e-19],
              [54.53, 1.427e-19],
              [61.20, 6.057e-20],
              [74.75, 1.677e-20]])

fig=plt.figure(figsize=(6,5))
plt.rc('font',**{'family':'serif', 'size':20})

b=b*1000.0
semilogy(b,SB,'ro-',label='theory, z=3')
semilogy(b,SB*3,'m-',label='theory, z=3')
C1=2.4e-18
b1=25.2
semilogy(b,C1*exp(-b/b1),'b:',label='Steidel+2011')
semilogy(Momose[:,0],Momose[:,1],'gx-',label='Momose+2014')
xlabel('$b$ $[\\rm pkpc]$')
ylabel('$\\rm SB_\\alpha$ $\\rm [erg/s/cm^2/arcsec^2]$')
xlim(0,300)
ylim(1e-22,1e-17)
lg = legend(loc='upper right',fontsize=15)
lg.draw_frame(False)

plt.tight_layout()

figure()
loglog(r_proper,emissivity,'ro-')
xlabel('$r$ $[\\rm pMpc]$')
ylabel('$\\epsilon_\\alpha(r)$ $\\rm [erg/s/cm^3 (proper)]$')
