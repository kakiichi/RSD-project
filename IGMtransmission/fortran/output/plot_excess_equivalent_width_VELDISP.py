def tau_eff_FG08(z):
    # FG08 formula
    tau_eff=0.0018*(1.+z)**3.92
    return tau_eff

def tau_eff_B13(z):
    # Becker+13 formula
    tau_eff_B13=0.751*((1.+z)/4.5)**2.90 - 0.132
    return tau_eff_B13

# Steidel+2010 Table 4 data
impact_parameter=array([31,63,103])
EW_steidel=array([2.01,1.23,0.92])
xerr1=array([16,16,20])
xerr2=array([16,20,20])
err=array([0.15,0.20,0.12])
#impact_parameter2=array([170,240])
#EW_steidel2=array([0.78,0.31])
#err2=array([0.08,0.05])
b_Rakic=array([0.0900,0.152,0.215,0.303,0.426,0.601,0.854,1.20,1.71])*1000.0 # pMpc to pkpc
xerr1_Rakic=array([  31. ,   31. ,   31.5,   44. ,   61.5,   87.5,  126.5,  173. ,255. ])
xerr2_Rakic=array([31.0,31.5,44.,61.5,87.5,126.5,173., 255.,255. ])

EW_Rakic=array([1.08,0.414,0.653,0.448,0.425,0.395,0.452,0.406,0.308])
EW2=array([1.25,0.639,0.985,0.607,0.536,0.487,0.515,0.452,0.355])
err_Rakic=array([0.17,0.225,0.332,0.159,0.111,0.092,0.063,0.046,0.047])

nfiles=3
title=range(nfiles)
title[0]='REF'
title[1]='HIGH DISP'
title[2]='LOW DISP'

fig=plt.figure(figsize=(6,5))
plt.rc('font',**{'family':'serif', 'size':16})
for n in range(nfiles):
    # load data
    if n==0:
        v_para,r_perp,s_para,s_perp,tau_eff=genfromtxt('RSD_01.dat',usecols=(2,3,4,5,6),unpack=True,comments='#')
    if n==1:
        v_para,r_perp,s_para,s_perp,tau_eff=genfromtxt('RSD_11.dat',usecols=(2,3,4,5,6),unpack=True,comments='#')
    if n==2:
        v_para,r_perp,s_para,s_perp,tau_eff=genfromtxt('RSD_12.dat',usecols=(2,3,4,5,6),unpack=True,comments='#')

    N=100
    tau_eff=tau_eff.reshape((N,N))
    v_para=v_para.reshape((N,N))
    s_perp=s_perp.reshape((N,N))
    r=s_perp[0,:]

    F=exp(-tau_eff)

    #tau_eff_mean=tau_eff.min()
    #tau_eff_mean=tau_eff_FG08(3.0)
    tau_eff_mean=tau_eff_B13(3.0)
    F_mean=exp(-tau_eff_mean)

    #dF=F_mean-F
    dF=1.-F/F_mean
    wavelength_lya=1215.67      # [Angstrom]
    c=2.998e5                   # [km/s]
    dv= v_para[1,0]-v_para[0,0] # [km/s]
    dEW=sum(dF[0:55,:],axis=0)*wavelength_lya*(dv/c) # integrate upto 500km/s (Rakic+2012 consistent)
    #dEW=sum(dF[0:90,:],axis=0)*wavelength_lya*(dv/c) # integrate upto 500km/s (Rakic+2012 consistent)

    if n==0:
        errorbar(impact_parameter,EW_steidel,yerr=err,xerr=[xerr1,xerr2],fmt='o',label='Steidel+2010')
        errorbar(b_Rakic,EW_Rakic,yerr=err_Rakic,fmt='s',xerr=[xerr1_Rakic,xerr1_Rakic],label='Rakic+2012')

    if n==0:
        semilogx(1000*r,2.*dEW,'k-',lw=2,label=title[n])
    if n==1:
        semilogx(1000*r,2.*dEW,'r--',lw=2,label=title[n])
    if n==2:
        semilogx(1000*r,2.*dEW,'b:',lw=2,label=title[n])
    #ylim(0.1,10.0)
    ylim(0,3.5)
    xlim(10,4000)
    ylabel('$\\langle \\rm EW\\rangle$ $[\\AA]$',fontsize=20)
    xlabel('$r_{\\perp}$ $[\\rm pkpc]$',fontsize=20)
    lg = legend(loc='upper right',fontsize=15)
    lg.draw_frame(False)

plt.tight_layout()
