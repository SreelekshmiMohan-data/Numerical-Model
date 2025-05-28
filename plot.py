import numpy as np
import matplotlib.pyplot as plt

nu = np.logspace(np.log10(0.01*1e9),np.log10(300e9),50)
total_flux0 = np.loadtxt("total_flux_qxprime0.txt", delimiter=",")

peak_index = np.where(total_flux0 == np.amax(total_flux0))
peak_nu = nu[peak_index]

plt.text(2.32e7,0.487,r"f${_{\nu}}$ ${\propto}$ ${\nu}^{+2.5}$")
plt.text(7.118e9,2.12,r"f${_{\nu}}$ ${\propto}$ ${\nu}^{-0.64}$")
plt.loglog(nu,total_flux0,color='k',linewidth = 2,label = r'$q_x^\prime$ = 0',zorder = 1)
ax = plt.axes()
plt.xlabel('Frequency (Hz)',fontsize = 17)
plt.ylabel('Flux density (mJy)',fontsize = 17)
ax.yaxis.set_label_coords(-0.14,0.5)
ax.xaxis.set_label_coords(0.5,-0.1)
plt.tick_params(axis='x', direction = 'in',length=8, labelsize=17,pad = 10)#pad = -15,
plt.tick_params(axis='y',direction = 'in',length=8, labelsize=17,pad = 10)
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e7,3e11)
#plt.ylim(0.00267,0.4)
#plt.minorticks_off()
#plt.legend(fontsize = 10,loc = 1)
ax = plt.axes()
#ax.yaxis.get_ticklocs(minor=True)
ax.tick_params(axis='y', which='minor', direction='in')
ax.tick_params(axis='x', which='minor', direction='in')
#plt.tight_layout()
plt.savefig('TotalFlux_x_const.eps', bbox_inches='tight')
plt.show()

frequency = np.ones(14)
fluxes = np.ones(14)
for k in range(14):
	frequency[k] = np.log10(nu[k+0]) 
	fluxes[k] = np.log10(total_flux0[k+0])
slope1, intercept = np.polyfit(frequency, fluxes, 1)
print('spectral index1',slope1)

frequency = np.ones(15)
fluxes = np.ones(15)
for k in range(15):
	frequency[k] = np.log10(nu[k+24]) 
	fluxes[k] = np.log10(total_flux0[k+24])
slope2, intercept = np.polyfit(frequency, fluxes, 1)
print('spectral index2',slope2)

