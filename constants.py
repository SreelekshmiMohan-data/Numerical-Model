import numpy as np
	
#general parms
au = 1.496e13
inc = 60.*(np.pi/180.)

r0 = 5000.*au #0.3*pc
theta_out = 30.*(np.pi/180.)
pc = 3.086e18
#y1max = 4500.*au #0.5*pc
y1max = 10000.*au
kpc = 3.086e21
d = 1.*kpc

y0 = r0*np.sin(inc)
x0 = 0.2 

#freefree parms
aj = 6.5e-38
ak = 0.212 

#synchrotron parms
q = 4.803*(10**-10)  #esu
m = 9.109e-28 #gms
c = 3e10 #cm/s
c1 = (3.*q)/(4.*np.pi*(m**3)*(c**5))
nu1 = 0.01*1e9 
nu2 = 100e9	
p = 2.3

w0_out = r0*(np.sin(theta_out/2.)/np.cos(theta_out/2.))
