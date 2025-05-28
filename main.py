import numpy as np
import matplotlib.pyplot as plt
from pylab import * 
from scipy import integrate
from scipy.special import gamma, factorial
from scipy.optimize import fsolve
import math
import angle
import constants
import thermal_functions
import nonthermal_functions
import time
import parameters_data
import jet_integration_summation
import jet_integration_summation2


t0 = time.time()
qx_prime = 0
relativistic_fraction = 1e-5
B0 = 0.3e-3
T0 = 1.e4
delta_theta = 0.5
epsilon = 1
qv = 0
qT = 0
qx = -0.5
qn = -qv - (2.*epsilon)
N0 = 500

nu = np.logspace(np.log10(0.01*1e9),np.log10(300e9),50)

#------------------------------------------------------------------------------------
def top_length(theta_in,epsilon):
	#inner cone
	theta_in_y1max = 2.*math.atan((np.sin(theta_in/2.)/np.cos(theta_in/2.))*((constants.y1max/constants.y0)**(epsilon-1.)) )
	if epsilon==1:
		theta_in_y2max_in = theta_in
	else:
		data = (epsilon,theta_in_y1max,constants.inc)
		theta_in_y2max_in = fsolve(angle.th1_solver,theta_in_y1max,args = data)
	
	s1_y1max_in =  (constants.y1max*np.sin(theta_in_y2max_in/2.))/(np.sin(constants.inc-theta_in_y2max_in/2.)*np.sin(constants.inc))
	y2max_in = constants.y1max + (s1_y1max_in*np.cos(constants.inc)*np.sin(constants.inc))
	w_y2max_in = (y2max_in/np.sin(constants.inc))*(np.sin(theta_in_y2max_in/2.)/np.cos(theta_in_y2max_in/2.))
	y3max_in = y2max_in +  (w_y2max_in*np.cos(constants.inc))
	#outer cone
	theta_out_y1max = 2.*math.atan((np.sin(constants.theta_out/2.)/np.cos(constants.theta_out/2.))*((constants.y1max/constants.y0)**(epsilon-1.)) )
	if epsilon==1:
		theta_out_y2max_out = constants.theta_out
	else:
		data = (epsilon,theta_out_y1max,constants.inc)
		theta_out_y2max_out = fsolve(angle.th1_solver,theta_out_y1max,args = data)
	s1_y1max_out =  (constants.y1max*np.sin(theta_out_y2max_out/2.))/(np.sin(constants.inc-theta_out_y2max_out/2.)*np.sin(constants.inc))
	y2max_out = constants.y1max + (s1_y1max_out*np.cos(constants.inc)*np.sin(constants.inc))
	w_y2max_out = (y2max_out/np.sin(constants.inc))*(np.sin(theta_out_y2max_out/2.)/np.cos(theta_out_y2max_out/2.))
	y3max_out = y2max_out +  (w_y2max_out*np.cos(constants.inc))

	return y2max_in,y3max_in,y2max_out,y3max_out

n0_ff = N0*(1.-relativistic_fraction)
n0_nt = N0*(relativistic_fraction)
theta_in = ((constants.theta_out*(180./np.pi)) - (2.*delta_theta))*np.pi/180.
 
free_parms = {'qx_prime':qx_prime,'N0':N0,'relativistic_fraction':relativistic_fraction,'n0_ff':n0_ff,'n0_nt':n0_nt,'B0':B0,'T0':T0,
'theta_in':theta_in,'epsilon':epsilon,'qv':qv,'qn':qn,'qT':qT,'qx':qx}


y2max_in,y3max_in,y2max_out,y3max_out = top_length(theta_in,epsilon)


ff_flux= np.ones(len(nu))
total_flux = np.ones(len(nu))
nt_flux = np.ones(len(nu))

for i in range(len(nu)):
			NU = nu[i]
			if NU<constants.nu2:
				total_flux[i] = jet_integration_summation.flux_final(constants.y1max,y2max_in,y2max_out,y3max_in,y3max_out,NU,free_parms)
			else:
				total_flux[i] = jet_integration_summation2.flux_final(constants.y1max,y2max_in,y2max_out,y3max_in,y3max_out,NU,free_parms)			
			print(total_flux[i],i)


np.savetxt('total_flux_r010au_theta040.txt',total_flux,delimiter=',')


t_end = time.time()
print('time(sec)',t_end-t0)




