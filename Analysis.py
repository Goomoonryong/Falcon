import numpy as np
import pdb
import matplotlib.pyplot as plt
from math import atan

# """ Design Variables """

alpha_w = -8.0
alpha0_t = 0                
cl_alpha = 15.0/1.3                     #""" Airfoils used in wing and Tail """
cl_alpha_t = 2*np.pi**2/180
Cm_0w = -0.2
Cm_af =0 


AR_W = 12.3
S   = 10.375                            # """ Wing and Tail Dimensions """
V_h = 0.497
l_t = 4.0
b = np.sqrt(AR_W*S)
c = S/b
S_t = V_h*S/(l_t*c)
b_t = 2
AR_t = b_t**2/S_t
e = 0.92
eta = 0.9

i_w = 0
i_t = 0
e0 = 0          # epsilon
ea = 0          # d(epsilon)/d(alpha)

x_cg = 0
x_ac = -0.25



# ''' Lift and Moment Coefficients '''
CL_alpha_w = cl_alpha/(1 + 57.3*cl_alpha/(np.pi*AR_W*e))
CL_alpha_t = cl_alpha_t/(1 + 57.3*cl_alpha_t/(np.pi*AR_t*e))

CL_alpha    = CL_alpha_w + eta*(S_t/S)*(1-ea)*CL_alpha_t                        
CL_0        = CL_alpha_w*(i_w-alpha0_w) + eta*(S_t/S)*CL_alpha_t*(i_t)
a0 	    = alpha_w*CL_alpha_w/CL_alpha + eta*(S_t/S)*CL_alpha_t*alpha0_t                   # Alpha(L=0) for wing & tail

print "CL0 ", CL_0

Cm_0        = Cm_0w - eta*V_h*CL_alpha_t*(i_t-e0+(1-ea)*alpha0) + Cm_af*az
x_np        = x_ac - (eta*V_h*CL_alpha_t*(1-ea) + Cm_af)/CL_alpha
Cm_alpha    = (x_np - x_cg)*CL_alpha


# ''' Trim Conditions '''

ae = 2.17/100
CL_delta = eta*(S_t/S)*ae
Cm_delta = -CL_delta*(l_t/c - (x_ac - x_cg))

alpha = np.arange(0,8,0.2)
delta = -Cm_0/Cm_delta - Cm_alpha*alpha/Cm_delta
print 'Delta: ', CL_delta
print 'Cm_delta: ', Cm_delta
print '------------------------------------------------------'
print '------------------------------------------------------'
CL  = CL_alpha*(alpha - a0)


cfe = 0.004
S_fus_vtail = 11.392
S_wet = S_fus_vtail + 2*S + 2*S_t
Cd0 = cfe*S_wet/S
Cdi = CL**2/(np.pi*AR_W*e)
CD = Cdi + Cd0

# ''' Performance Parameters '''
Range = CL/CD
endurance = CL**1.5/CD
W = 160*9.81
rho = 1.225
v_stall = 30*5./18
lift_stall = 0.5*rho*v_stall**2*np.max(CL)*1.75*S
sink = np.sqrt(2*W/S/rho)/endurance
