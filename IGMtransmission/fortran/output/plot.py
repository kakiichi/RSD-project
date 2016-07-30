v_para,r_perp,s_para,s_perp,tau_eff=genfromtxt('RSD_04.dat',usecols=(2,3,4,5,6),unpack=True,comments='#')

tau_eff=tau_eff.reshape((100,100))
F=exp(-tau_eff)


fig=plt.figure(figsize=(8,5))
plt.rc('font',**{'family':'serif', 'size':16})

extent=[s_perp.min(),s_perp.max(),s_para.min(),s_para.max()]

im=imshow(F,cmap='RdYlBu',extent=extent,origin='lower')#,vmax=0.8,vmin=0.0)
ylabel('$s_{\\parallel}$ $[\\rm pMpc]$',fontsize=20)
xlabel('$s_{\\perp}$ $[\\rm pMpc]$',fontsize=20)
contour(F,15,extent=extent,colors='black',alpha=0.5)
cax = fig.add_axes([0.8, 0.2, 0.03, 0.7])
cbar=fig.colorbar(im, cax=cax)
cbar.set_label('$\\langle F_{\\alpha}(s_{\\parallel},s_{\\perp})\\rangle$',fontsize=20)
plt.tight_layout()



fig=plt.figure(figsize=(8,5))
plt.rc('font',**{'family':'serif', 'size':16})

extent=[s_perp.min(),s_perp.max(),s_para.min(),s_para.max()]
im=imshow(log10(-log(F)),cmap='RdYlBu_r',extent=extent)#,vmax=0.8,vmin=0.0)
ylabel('$\\pi$ $[\\rm pMpc]$',fontsize=20)
xlabel('$\\sigma$ $[\\rm pMpc]$',fontsize=20)
contour(log10(-log(F)),15,extent=extent,colors='black',alpha=0.5)
cax = fig.add_axes([0.8, 0.2, 0.03, 0.7])
cbar=fig.colorbar(im, cax=cax)
cbar.set_label('$\\log_{10}$ $\\tau_{\\rm eff}(\\pi,\\sigma)$',fontsize=20)
plt.tight_layout()

# fig=plt.figure(figsize=(8,5))
# plt.rc('font',**{'family':'serif', 'size':16})

# extent=[s_perp.min(),s_perp.max(),s_para.min(),s_para.max()]
# im=imshow(log10(xi),aspect=1,extent=extent,cmap='RdYlBu_r')
# contour(log10(xi),15,extent=extent)
# ylabel('$s_{\\parallel}$ $[h^{-1}\\rm cMpc]$',fontsize=20)
# xlabel('$s_{\\perp}$ $[h^{-1}\\rm cMpc]$',fontsize=20)
# cax = fig.add_axes([0.8, 0.2, 0.03, 0.7])
# cbar=fig.colorbar(im, cax=cax)
# cbar.set_label('$\\log_{10}$ $\\xi_s$',fontsize=20)

# plt.tight_layout()
