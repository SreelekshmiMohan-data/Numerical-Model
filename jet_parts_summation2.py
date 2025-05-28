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


divisions = 5.

def tau_b(s,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2):
	delta_s = (upper2-lower2)/divisions
	n0_ff = free_parms['n0_ff']
	tau = combined_alpha_nu_b(s,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms)*delta_s
	return tau

def sigma_tau_b(s,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2,n_th):
	delta_s = (upper2-lower2)/divisions
	n0_ff = free_parms['n0_ff']
	sigma_tau=0.
	for i2 in range(int(n_th)):
		s_i = lower2 + 0.5*(delta_s*( 2.*(n_th-1.) + 1.))
		sigma_tau += tau_b(s_i,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2)
	return sigma_tau


def dI_b(s,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2,n_th):
	n0_ff = free_parms['n0_ff']
	s_nu2 = combined_j_nu_b(s,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms)/combined_alpha_nu_b(s,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms)
	tau_term1 = (1. - np.exp(-tau_b(s,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2)  ))
	tau_term2 = np.exp(-(sigma_tau_b(s,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2,n_th) - tau_b(s,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2)))
	return s_nu2*tau_term1*tau_term2


def combined_j_nu_b(s,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms):
	n0_ff = free_parms['n0_ff']

	j_nu_ff = thermal_functions.j_nu2(s,vprime,y,s2_vprime_out,theta_vprime_out,n0_ff,NU,w0_out,free_parms)
	j_nu_synch = nonthermal_functions.j_nu2(s,vprime,theta_vprime_out,y,s2_vprime_out,w0_out,free_parms,constants.nu2)*((NU/constants.nu2)**(-constants.p/2.))
	
	j_nu = j_nu_ff + j_nu_synch
	return j_nu

def combined_alpha_nu_b(s,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms):
	n0_ff = free_parms['n0_ff']

	alpha_nu_ff = thermal_functions.alpha_nu2(s,vprime,y,s2_vprime_out,theta_vprime_out,n0_ff,NU,w0_out,free_parms)
	alpha_nu_synch = nonthermal_functions.alpha_nu2(s,vprime,theta_vprime_out,y,s2_vprime_out,w0_out,free_parms,constants.nu2)*((NU/constants.nu2)**(-(constants.p+5.)/2.))
	
	alpha_nu = alpha_nu_ff + alpha_nu_synch
	return alpha_nu


def flux_b(vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2):
	delta_s = (upper2-lower2)/divisions
	#print(divisions,delta_s,y)
	intensity =0.
	for i1 in range(int(divisions)):
		s_i = lower2 + 0.5*(delta_s*( (2.*i1) + 1.))
		n_th = i1+1.
		intensity += dI_b(s_i,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2,n_th)
	return intensity*1.e26	

def free_tau_inner1(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out):

		#innner cone free-free tau
		N0 = free_parms['N0']

		tau_ff_1_in = integrate.quad(thermal_functions.alpha_nu1,0.,s1_vprime_in,args=(vprime,y,s1_vprime_in,theta_vprime_in,N0,NU,w0_in,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
		tau_ff_2_in = integrate.quad(thermal_functions.alpha_nu2,0.,s2_vprime_in,args=(vprime,y,s2_vprime_in,theta_vprime_in,N0,NU,w0_in,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
		tau_ff_inner = tau_ff_1_in + tau_ff_2_in
		return tau_ff_inner	

def free_inner1(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out):

		#innner cone free-free flux
		tau_ff_inner = free_tau_inner1(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)
		
		const = (constants.aj*free_parms['T0']*(NU**2))/(constants.ak) #2.*w**(constants.d**2)
		f_ff_inner = const*(1.- np.exp(-tau_ff_inner))*1.e26
		return f_ff_inner


def tau_f(s,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1):
	delta_s = (upper1-lower1)/divisions
	n0_ff = free_parms['n0_ff']
	tau = combined_alpha_nu_f(s,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms)*delta_s
	return tau

def sigma_tau_f(s,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1,n_th):
	delta_s = (upper1-lower1)/divisions
	n0_ff = free_parms['n0_ff']
	sigma_tau=0.
	for i2 in range(int(n_th)):
		s_i = lower1 + 0.5*(delta_s*( 2.*(n_th-1.) + 1.))
		sigma_tau += tau_f(s_i,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1)
	return sigma_tau


def dI_f(s,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1,n_th):
	n0_ff = free_parms['n0_ff']
	s_nu1 = combined_j_nu_f(s,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms)/combined_alpha_nu_f(s,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms)
	tau_term1 = (1. - np.exp(-tau_f(s,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1)  ))
	tau_term2 = np.exp(-(sigma_tau_f(s,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1,n_th) - tau_f(s,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1)))
	return s_nu1*tau_term1*tau_term2


def combined_j_nu_f(s,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms):
	n0_ff = free_parms['n0_ff']

	j_nu_ff = thermal_functions.j_nu1(s,vprime,y,s1_vprime_out,theta_vprime_out,n0_ff,NU,w0_out,free_parms)
	j_nu_synch = nonthermal_functions.j_nu1(s,vprime,theta_vprime_out,y,s1_vprime_out,w0_out,free_parms,constants.nu2)*((NU/constants.nu2)**(-constants.p/2.))
	
	j_nu = j_nu_ff + j_nu_synch
	return j_nu

def combined_alpha_nu_f(s,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms):
	n0_ff = free_parms['n0_ff']

	alpha_nu_ff = thermal_functions.alpha_nu1(s,vprime,y,s1_vprime_out,theta_vprime_out,n0_ff,NU,w0_out,free_parms)
	alpha_nu_synch = nonthermal_functions.alpha_nu1(s,vprime,theta_vprime_out,y,s1_vprime_out,w0_out,free_parms,constants.nu2)*((NU/constants.nu2)**(-(constants.p+5.)/2.))
	
	alpha_nu = alpha_nu_ff + alpha_nu_synch
	return alpha_nu


def flux_f(vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1):
	delta_s = (upper1-lower1)/divisions
	#print(divisions,delta_s,y)
	intensity =0.
	for i1 in range(int(divisions)):
		s_i = lower1 + 0.5*(delta_s*( (2.*i1) + 1.))
		n_th = i1+1.
		intensity += dI_f(s_i,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1,n_th)
	return intensity*1.e26	

############################################################################################################################################
#yamax to ybmax_in

def free_tau_inner2(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,WE_in,WE_out,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out):

		N0 = free_parms['N0']

		tau_ff_1_in = integrate.quad(thermal_functions.alpha_nu1,WE_in,s1_vprime_in,args=(vprime,y,s1_vprime_in,theta_vprime_in,N0,NU,w0_in,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
		tau_ff_2_in = integrate.quad(thermal_functions.alpha_nu2,0.,s2_vprime_in,args=(vprime,y,s2_vprime_in,theta_vprime_in,N0,NU,w0_in,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
		tau_ff_inner = tau_ff_1_in + tau_ff_2_in
		return tau_ff_inner


def free_inner2(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,WE_in,WE_out,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out):
		tau_ff_inner = free_tau_inner2(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,WE_in,WE_out,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)
		const = (constants.aj*free_parms['T0']*(NU**2))/(constants.ak) #2.*w**(constants.d**2)
		f_ff_inner = const*(1.- np.exp(-tau_ff_inner))*1.e26	
		return f_ff_inner



#############################################################################################################################################
#ybmax_in to ybmax_out


def free_tau_inner3(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,WE_in,WE_out,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out):
		#innner cone free-free
		N0 = free_parms['N0']
	
		tau_ff_2_in = integrate.quad(thermal_functions.alpha_nu2,GH_in,s2_vprime_in,args=(vprime,y,s2_vprime_in,theta_vprime_in,N0,NU,w0_in,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
		tau_ff_inner = tau_ff_2_in
		return tau_ff_inner


def free_inner3(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,WE_in,WE_out,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out):
		tau_ff_inner = free_tau_inner3(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,WE_in,WE_out,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)			
		const = (constants.aj*free_parms['T0']*(NU**2))/(constants.ak) #2.*w**(constants.d**2)
		f_ff_inner = const*(1.- np.exp(-tau_ff_inner))*1.e26	
		return f_ff_inner



################################################################################################################################################
#ybmax_out to ycmax_in

def free_tau_inner4(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out):
		#innner cone free-free
		N0 = free_parms['N0']

		tau_ff_2_in = integrate.quad(thermal_functions.alpha_nu2,GH_in,s2_vprime_in,args=(vprime,y,s2_vprime_in,theta_vprime_in,N0,NU,w0_in,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
		tau_ff_inner = tau_ff_2_in
		return tau_ff_inner 	

def free_inner4(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out):
		tau_ff_inner = free_tau_inner4(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)
		const = (constants.aj*free_parms['T0']*(NU**2))/(constants.ak) #2.*w**(constants.d**2)
		f_ff_inner = const*(1.- np.exp(-tau_ff_inner))*1.e26	
		return f_ff_inner





