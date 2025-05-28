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
import jet_parts_summation

divisions = 5.

def Theta_y_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		part1_in = (np.sin(theta_in/2.)/np.cos(theta_in/2.))
		part2_in = (y/constants.y0)**(free_parms['epsilon']-1.)
		theta_y_in = 2.*math.atan(part1_in*part2_in)
		return theta_y_in

def W_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		w_in = w0_in*(y/constants.y0)**free_parms['epsilon']
		return w_in

def W_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		w_in = W_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		w_vprime_in = np.sqrt(w_in**2 - vprime**2)
		return w_vprime_in

def Theta_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		w_in = W_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		w_vprime_in = W_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_y_in = Theta_y_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_in = 2.*math.atan( (w_vprime_in/w_in)*(np.sin(theta_y_in/2.)/np.cos(theta_y_in/2.)) ) 
		return theta_vprime_in

def S1_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		w_vprime_in = W_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_y_in = Theta_y_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		#theta_vprime_in = Theta_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_in = theta_y_in
		if free_parms['epsilon'] ==1:
			s1_vprime_in = w_vprime_in*(np.cos(theta_vprime_in/2.)/np.sin(constants.inc-theta_vprime_in/2.))

		else:
			data = (free_parms['epsilon'],theta_vprime_in,constants.inc)
			theta1_vprime_in = fsolve(angle.th1_solver,theta_vprime_in,args = data)
			dist_from_axis = y*np.sin(theta1_vprime_in/2.)/(np.sin(constants.inc)*np.sin(constants.inc - (theta1_vprime_in/2.)))
			s1_vprime_in = np.sqrt(dist_from_axis**2 - vprime**2)
			#print('dist1',dist_from_axis,vprime)
		return s1_vprime_in


def S2_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		w_vprime_in = W_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_y_in = Theta_y_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		#theta_vprime_in = Theta_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_in = theta_y_in
		if free_parms['epsilon'] ==1:
			s2_vprime_in = w_vprime_in*(np.cos(theta_vprime_in/2.)/np.sin(constants.inc+theta_vprime_in/2.))
		else:
			data = (free_parms['epsilon'],theta_vprime_in,constants.inc)
			theta2_vprime_in = fsolve(angle.th2_solver,theta_vprime_in,args = data)
			dist_from_axis = y*np.sin(theta2_vprime_in/2.)/(np.sin(constants.inc)*np.sin(constants.inc + (theta2_vprime_in/2.)))
			s2_vprime_in = np.sqrt(dist_from_axis**2 - vprime**2)
			#print('dist2',dist_from_axis,vprime)
		return s2_vprime_in


def Theta_y_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		part1_out = (np.sin(theta_out/2.)/np.cos(theta_out/2.))
		part2_out = (y/constants.y0)**(free_parms['epsilon']-1.)
		theta_y_out = 2.*math.atan(part1_out*part2_out)
		return theta_y_out

def W_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		w_out = w0_out*(y/constants.y0)**free_parms['epsilon']
		return w_out

def W_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		w_out = W_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		w_vprime_out = np.sqrt(w_out**2 - vprime**2)
		return w_vprime_out


def Theta_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		w_out = W_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		w_vprime_out = W_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_y_out = Theta_y_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_out = 2.*math.atan( (w_vprime_out/w_out)*(np.sin(theta_y_out/2.)/np.cos(theta_y_out/2.)) ) 
		return theta_vprime_out


def S1_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		w_vprime_out = W_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_y_out = Theta_y_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		#theta_vprime_out = Theta_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_out = theta_y_out
		if free_parms['epsilon'] ==1:
			s1_vprime_out = w_vprime_out*(np.cos(theta_vprime_out/2.)/np.sin(constants.inc- (theta_vprime_out/2.)))
		else:
			data = (free_parms['epsilon'],theta_vprime_out,constants.inc)
			theta1_vprime_out = fsolve(angle.th1_solver,theta_vprime_out,args = data)
			dist_from_axis = y*np.sin(theta1_vprime_out/2.)/(np.sin(constants.inc)*np.sin(constants.inc - (theta1_vprime_out/2.)))
			s1_vprime_out = np.sqrt(dist_from_axis**2 - vprime**2)
		return s1_vprime_out


def S2_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		w_vprime_out = W_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_y_out = Theta_y_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		#theta_vprime_out = Theta_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_out = theta_y_out
		if free_parms['epsilon'] ==1:
			s2_vprime_out = w_vprime_out*(np.cos(theta_vprime_out/2.)/np.sin(constants.inc+ (theta_vprime_out/2.)))
		else:
			data = (free_parms['epsilon'],theta_vprime_out,constants.inc)
			theta2_vprime_out = fsolve(angle.th2_solver,theta_vprime_out,args = data)
			dist_from_axis = y*np.sin(theta2_vprime_out/2.)/(np.sin(constants.inc)*np.sin(constants.inc + (theta2_vprime_out/2.)))
			s2_vprime_out = np.sqrt(dist_from_axis**2 - vprime**2)
		return s2_vprime_out


#constants required : theta0,r0,y0,free_parms['epsilon'],inc,aj,ak,d,free_parms['T0']
#y0 to yamax
def ff1(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):

		s1_vprime_in = 	S1_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s2_vprime_in = S2_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s1_vprime_out = S1_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)		
		s2_vprime_out = S2_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_in = Theta_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_out = Theta_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)

		lower1 = 0.
		upper1 = s1_vprime_out - s1_vprime_in
		lower2 = s2_vprime_in
		upper2 = s2_vprime_out 
		
		tau_ff_inner = jet_parts_summation.free_tau_inner1(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)
		f_ff_inner = jet_parts_summation.free_inner1(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)

		f_behind = jet_parts_summation.flux_b(vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2)
		tau_front = jet_parts_summation.sigma_tau_f(upper1,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1,divisions)
		f_front = jet_parts_summation.flux_f(vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1)

		I_los = f_behind*np.exp(-(tau_ff_inner + tau_front))+ f_ff_inner*np.exp(-tau_front) + f_front
		return I_los

############################################################################################################################################
#yamax to ybmax_in
def ff2(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		
		
		s1_vprime_in = 	S1_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s2_vprime_in = S2_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s1_vprime_out = S1_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)		
		s2_vprime_out = S2_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_in = Theta_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_out = Theta_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
			
		center_dist_in = (ybmax_in-y)/(np.sin(constants.inc)*np.cos(constants.inc))
		center_dist_out = (ybmax_out-y)/(np.sin(constants.inc)*np.cos(constants.inc))
		WE_in =  s1_vprime_in- center_dist_in 
		WE_out =  s1_vprime_out- center_dist_out
		GH_in = (y-ybmax_in)/(np.sin(constants.inc)*np.cos(constants.inc))
		GH_out = (y-ybmax_out)/(np.sin(constants.inc)*np.cos(constants.inc))
		
		lower1 = WE_out
		upper1 = s1_vprime_out - (s1_vprime_in- WE_in)
		lower2 = s2_vprime_in
		upper2 = s2_vprime_out 	
		#print(WE_in,WE_out,y)
		#print('1',upper1,lower1,upper1-lower1)

		tau_ff_inner = jet_parts_summation.free_tau_inner2(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,WE_in,WE_out,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)
		f_ff_inner = jet_parts_summation.free_inner2(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,WE_in,WE_out,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)
		
		f_behind = jet_parts_summation.flux_b(vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2)
		tau_front = jet_parts_summation.sigma_tau_f(upper1,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1,divisions)
		f_front = jet_parts_summation.flux_f(vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1)

		I_los = f_behind*np.exp(-(tau_ff_inner + tau_front))+ f_ff_inner*np.exp(-tau_front) + f_front
		#I_los = f_behind*np.exp(-(tau_ff_inner + tau_front))+ f_ff_inner*np.exp(-tau_front) 

		return I_los


#############################################################################################################################################
#ybmax_in to ybmax_out
def ff3(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		
		s1_vprime_in = 	S1_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s2_vprime_in = S2_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s1_vprime_out = S1_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)		
		s2_vprime_out = S2_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_in = Theta_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_out = Theta_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
			
		center_dist_in = (ybmax_in-y)/(np.sin(constants.inc)*np.cos(constants.inc))
		center_dist_out = (ybmax_out-y)/(np.sin(constants.inc)*np.cos(constants.inc))
		WE_in =  s1_vprime_in- center_dist_in 
		WE_out =  s1_vprime_out- center_dist_out
		GH_in = (y-ybmax_in)/(np.sin(constants.inc)*np.cos(constants.inc))
		GH_out = (y-ybmax_out)/(np.sin(constants.inc)*np.cos(constants.inc))
		
		lower1 = WE_out
		upper1 = s1_vprime_out 
		lower2 = s2_vprime_in
		upper2 = s2_vprime_out		
		lower3 = 0.
		upper3 = GH_in
		#print('3',upper2-lower2)	

		tau_ff_inner = jet_parts_summation.free_tau_inner3(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,WE_in,WE_out,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)
		f_ff_inner = jet_parts_summation.free_inner3(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,WE_in,WE_out,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)

		f_behind = jet_parts_summation.flux_b(vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2)
		tau_front1 = jet_parts_summation.sigma_tau_f(upper1,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1,divisions)
		f_front1 = jet_parts_summation.flux_f(vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1)
		tau_front2 = jet_parts_summation.sigma_tau_b(upper3,vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower3,upper3,divisions)
		f_front2 = jet_parts_summation.flux_b(vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower3,upper3)
		
		I_los = f_behind*np.exp(-(tau_ff_inner + tau_front2 + tau_front1))+ f_ff_inner*np.exp(-(tau_front2 + tau_front1)) + f_front2*np.exp(-tau_front1) + f_front1
		#I_los = f_behind*np.exp(-(tau_ff_inner + tau_front1))+ f_ff_inner*np.exp(- tau_front1) +  f_front1
		#print(I_los)
		return I_los



################################################################################################################################################
#ybmax_out to ycmax_in
def ff4(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		
		s1_vprime_in = 	S1_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s2_vprime_in = S2_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s1_vprime_out = S1_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)		
		s2_vprime_out = S2_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)

		theta_vprime_in = Theta_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_out = Theta_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
			
		GH_in = (y-ybmax_in)/(np.sin(constants.inc)*np.cos(constants.inc))
		GH_out = (y-ybmax_out)/(np.sin(constants.inc)*np.cos(constants.inc))
		
		
		lower1 = GH_out
		upper1 = GH_in
		lower2 = s2_vprime_in
		upper2 = s2_vprime_out	
		#print('4',upper2-lower2)
	
		tau_ff_inner = jet_parts_summation.free_tau_inner4(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)
		f_ff_inner = jet_parts_summation.free_inner4(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms,GH_in,GH_out,s1_vprime_in,s1_vprime_out,s2_vprime_in,s2_vprime_out,theta_vprime_in,theta_vprime_out)

		f_behind = jet_parts_summation.flux_b(vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2)
		tau_front = jet_parts_summation.sigma_tau_f(upper1,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1,divisions)
		f_front = jet_parts_summation.flux_f(vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1)

		I_los = f_behind*np.exp(-(tau_ff_inner + tau_front))+ f_ff_inner*np.exp(-tau_front) + f_front

		return I_los




###############################################################################################################################################
#ycmax_in to ycmax_out
def ff5(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		
		s1_vprime_in = 	S1_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s2_vprime_in = S2_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s1_vprime_out = S1_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)		
		s2_vprime_out = S2_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)

		
		theta_vprime_in = Theta_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_out = Theta_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
	
		GH_out = (y-ybmax_out)/(np.sin(constants.inc)*np.cos(constants.inc))
		
		
		lower2 = GH_out
		upper2 = s2_vprime_out
		#print('5',upper2-lower2)		

		f_behind = jet_parts_summation.flux_b(vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2)

		I_los = f_behind

		return I_los


###############################################################################################################################################
#y0 to yamax
def ff6(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		
		s1_vprime_in = 	S1_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s2_vprime_in = S2_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s1_vprime_out = S1_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)		
		s2_vprime_out = S2_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)

		
		theta_vprime_in = Theta_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_out = Theta_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
	
		lower1 = 0.
		upper1 = s1_vprime_out 
		lower2 = 0.
		upper2 = s2_vprime_out
		
		f_behind = jet_parts_summation.flux_b(vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2)
		tau_front = jet_parts_summation.sigma_tau_f(upper1,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1,divisions)
		f_front = jet_parts_summation.flux_f(vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1)

		I_los = f_behind*np.exp(-tau_front) + f_front

		return I_los



#yamax to ybmax_out
def ff7(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		
		s1_vprime_in = 	S1_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s2_vprime_in = S2_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s1_vprime_out = S1_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)		
		s2_vprime_out = S2_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)

		theta_vprime_in = Theta_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_out = Theta_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		
		center_dist_in = (ybmax_in-y)/(np.sin(constants.inc)*np.cos(constants.inc))
		center_dist_out = (ybmax_out-y)/(np.sin(constants.inc)*np.cos(constants.inc))
		

		WE_out =  s1_vprime_out- center_dist_out 
		lower1 = WE_out
		upper1 = s1_vprime_out 
		lower2 = 0
		upper2 = s2_vprime_out
		#print('7',lower1)		

		f_behind = jet_parts_summation.flux_b(vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2)
		tau_front = jet_parts_summation.sigma_tau_f(upper1,vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1,divisions)
		f_front = jet_parts_summation.flux_f(vprime,y,s1_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower1,upper1)

		I_los = f_behind*np.exp(-tau_front) + f_front

		return I_los

#ybmax_out to ycmax_in
def ff8(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
		
		s1_vprime_in = 	S1_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s2_vprime_in = S2_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		s1_vprime_out = S1_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)		
		s2_vprime_out = S2_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)

		theta_vprime_in = Theta_vprime_in(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
		theta_vprime_out = Theta_vprime_out(vprime,y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms)
			
		WE_out = s1_vprime_out - (ybmax_out-y)/(np.sin(constants.inc)*np.cos(constants.inc))
		GH_out = (y-ybmax_out)/(np.sin(constants.inc)*np.cos(constants.inc))


		lower2 = GH_out
		upper2 = s2_vprime_out
		#print('8',upper2-lower2)
		
		f_behind = jet_parts_summation.flux_b(vprime,y,s2_vprime_out,theta_vprime_out,NU,w0_out,free_parms,lower2,upper2)

		I_los = f_behind

		return I_los


######################################################################################################
######################################################################################################
######################################################################################################

#y0 to yamax
def func1(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
	w_in = w0_in*(y/constants.y0)**free_parms['epsilon']
	ff0 = integrate.quad(ff1,0.,w_in,args=(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
	#print('1',ff0,y)
	return 2.*ff0


#yamax to ybmax_in
def func2(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
	w_in = w0_in*(y/constants.y0)**free_parms['epsilon']
	#print('w_in',w_in,y)
	ff0 = integrate.quad(ff2,0.,w_in,args=(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
	#print(ff0)
	#ff0 = 2.
	return 2.*ff0


#ybmax_in to ybmax_out
def func3(y,theta_in,theta_out,w0_in,w0_out,w_ybmax_in,ybmax_in,ybmax_out,NU,free_parms):
	delta_y = y-ybmax_in
	w_top = np.sqrt(w_ybmax_in**2 - ((delta_y/np.cos(constants.inc))**2)  )
	ff0 = integrate.quad(ff3,0.,w_top,args=(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
	return 2.*ff0

#ybmax_out to ycmax_in
def func4(y,theta_in,theta_out,w0_in,w0_out,w_ybmax_in,ybmax_in,ybmax_out,NU,free_parms):
	delta_y = y-ybmax_in
	w_top = np.sqrt(w_ybmax_in**2 - ((delta_y/np.cos(constants.inc))**2)  )
	ff0 = integrate.quad(ff4,0.,w_top,args=(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
	return 2.*ff0

#ycmax_in to ycmax_out
def func5(y,theta_in,theta_out,w0_in,w0_out,w_ybmax_out,ybmax_in,ybmax_out,NU,free_parms):
	delta_y = y-ybmax_out
	w_top = np.sqrt(w_ybmax_out**2 - ((delta_y/np.cos(constants.inc))**2)  )
	ff0 = integrate.quad(ff5,0.,w_top,args=(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
	return 2.*ff0

#edges--------------------------------------
#y0 to yamax
def func6(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
	w_in = w0_in*(y/constants.y0)**free_parms['epsilon']
	w_out = w0_out*(y/constants.y0)**free_parms['epsilon']
	ff0 = integrate.quad(ff6,w_in,w_out,args=(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
	#print('6',ff0,y)
	return 2.*ff0


#yamax to ybmax_in
def func7(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms):
	w_in = w0_in*(y/constants.y0)**free_parms['epsilon']
	w_out = w0_out*(y/constants.y0)**free_parms['epsilon']
	ff0 = integrate.quad(ff7,w_in,w_out,args=(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
	return 2.*ff0

#ybmax_in to ybmax_out
def func8(y,theta_in,theta_out,w0_in,w0_out,w_ybmax_in,ybmax_in,ybmax_out,NU,free_parms):
	w_in = w0_in*(y/constants.y0)**free_parms['epsilon']
	w_out = w0_out*(y/constants.y0)**free_parms['epsilon']
	delta_y = y-ybmax_in
	w_top_in = np.sqrt(w_ybmax_in**2 - ((delta_y/np.cos(constants.inc))**2)  )
	ff0 = integrate.quad(ff7,w_top_in,w_out,args=(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
	return 2.*ff0

#ybmax_out to ycmax_in
def func9(y,theta_in,theta_out,w0_in,w0_out,w_ybmax_in,w_ybmax_out,ybmax_in,ybmax_out,NU,free_parms):
	w_in = w0_in*(y/constants.y0)**free_parms['epsilon']
	w_out = w0_out*(y/constants.y0)**free_parms['epsilon']
	delta_y = y-ybmax_in
	w_top_in = np.sqrt(w_ybmax_in**2 - ((delta_y/np.cos(constants.inc))**2)  )
	delta_y = y-ybmax_out
	w_top_out = np.sqrt(w_ybmax_out**2 - ((delta_y/np.cos(constants.inc))**2)  )
	ff0 = integrate.quad(ff8,w_top_in,w_top_out,args=(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
	return 2.*ff0

#ybmax_out to ycmax_in
def func10(y,theta_in,theta_out,w0_in,w0_out,w_ybmax_out,ybmax_in,ybmax_out,NU,free_parms):
	w_in = w0_in*(y/constants.y0)**free_parms['epsilon']
	w_out = w0_out*(y/constants.y0)**free_parms['epsilon']
	
	delta_y = y-ybmax_out
	w_top_out = np.sqrt(w_ybmax_out**2 - ((delta_y/np.cos(constants.inc))**2)  )
	ff0 = integrate.quad(ff8,0.,w_top_out,args=(y,theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]
	return 2.*ff0


##############################################################################################################################################
##############################################################################################################################################



#integration over y
def flux_final(yamax,ybmax_in,ybmax_out,ycmax_in,ycmax_out,NU,free_parms):
	theta_in = free_parms['theta_in']	
	theta_out = constants.theta_out
	w0_in = constants.r0*(np.sin(theta_in/2.)/np.cos(theta_in/2.))
	w0_out = constants.r0*(np.sin(theta_out/2.)/np.cos(theta_out/2.))		
	w_ybmax_in = w0_in*(ybmax_in/constants.y0)**free_parms['epsilon']#######
	w_ybmax_out = w0_out*(ybmax_out/constants.y0)**free_parms['epsilon']
	
	flux1 = integrate.quad(func1,constants.y0,yamax,args=(theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]/(constants.d**2)  +   integrate.quad(func6,constants.y0,yamax,args=(theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]/(constants.d**2)

	
	flux2 = integrate.quad(func2,yamax,ybmax_in,args=(theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]/(constants.d**2) + integrate.quad(func7,yamax,ybmax_in,args=(theta_in,theta_out,w0_in,w0_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]/(constants.d**2)
	
	
	flux3 = integrate.quad(func3,ybmax_in,ybmax_out,args=(theta_in,theta_out,w0_in,w0_out,w_ybmax_in,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]/(constants.d**2) + integrate.quad(func8,ybmax_in,ybmax_out,args=(theta_in,theta_out,w0_in,w0_out,w_ybmax_in,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]/(constants.d**2)

	
	flux4 = integrate.quad(func4,ybmax_out,ycmax_in,args=(theta_in,theta_out,w0_in,w0_out,w_ybmax_in,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]/(constants.d**2) + integrate.quad(func9,ybmax_out,ycmax_in,args=(theta_in,theta_out,w0_in,w0_out,w_ybmax_in,w_ybmax_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]/(constants.d**2)

	
	flux5 = integrate.quad(func5,ycmax_in,ycmax_out,args=(theta_in,theta_out,w0_in,w0_out,w_ybmax_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]/(constants.d**2) + integrate.quad(func10,ycmax_in,ycmax_out,args=(theta_in,theta_out,w0_in,w0_out,w_ybmax_out,ybmax_in,ybmax_out,NU,free_parms),epsabs=1.e-1, epsrel=1.e-1)[0]/(constants.d**2)
	


	return flux1+flux2+flux3+flux4+flux5



