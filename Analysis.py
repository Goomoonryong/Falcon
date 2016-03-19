import numpy as np
import pdb
import matplotlib.pyplot as plt
from math import atan

rad2deg = 180./np.pi

# AOA at zero lift
alpha0_w = -8.2
alpha0_t = 0  

alpha = np.arange(-5,13,0.2)


# Lift curve slope for airfoils - per deg
cl_alpha_w = 0.1
cl_alpha_t = 2*np.pi**2/180


Cm0_w = -0.2 # got to be changed
Cm0_f = 0

# Wing parameters
b = 10.3
S   = 11.5	
AR_W = b**2/S
c = S/b
Cm_alpha_f  = 2*1.408/S/c/rad2deg


# Tail parameters
l_t = 4.72         
S_t = 1.432      
b_t = 2
AR_t = b_t**2/S_t
V_h = S_t*l_t/(S*c)

# coefficients of efficiency
e = 0.92
eta = 0.9

# wing and tail setting angles
i_w = 0
i_t = -1.98010788483

# locations of CG and AC
x_cg = -0.3928/c
x_ac = -c/(4*c)



# ''' Lift and Moment Coefficients '''
CL_alpha_w = cl_alpha_w/(1 + rad2deg*cl_alpha_w/(np.pi*AR_W*e))
CL_alpha_t = cl_alpha_t/(1 + rad2deg*cl_alpha_t/(np.pi*AR_t*e))

eps0 = CL_alpha_w*(i_w - alpha0_w)*rad2deg/(np.pi*AR_W*e)
eps_alpha = CL_alpha_w*rad2deg/(np.pi*AR_W*e)


CL_alpha = CL_alpha_w + eta*(S_t/S)*CL_alpha_t*(1 - eps_alpha)
CL_0 = CL_alpha_w*(i_w-alpha0_w) + eta*(S_t/S)*CL_alpha_t*(i_t - eps0)

# Alpha(L=0) wing and tail
alpha0 	= -CL_0/CL_alpha     

Cm_alpha = CL_alpha*(x_ac - x_cg) - eta*V_h*(1 - eps_alpha)*CL_alpha_t + Cm_alpha_f
Cm_0     = Cm0_w + CL_0*(x_ac - x_cg) - eta*V_h*(i_t - eps0)*CL_alpha_t + Cm0_f

x_NP = x_ac - (eta*V_h*(1 - eps_alpha)*CL_alpha_t + Cm_alpha_f)/CL_alpha

SM = 0.1*c
Des_x_NP = x_cg*c - SM
Des_Vh = (x_ac - Cm_alpha_f/CL_alpha - Des_x_NP)/(eta*(1 - eps_alpha)*CL_alpha_t/CL_alpha)
Des_lt = Des_Vh*S*c/S_t
print 'x_NP', x_NP
print 'Desired lt', Des_lt

""" Trim Conditons """

CL_delta_t = 0.04
CL_delta = eta*(S_t/S)*CL_delta_t
Cm_delta = -CL_delta*(l_t/c - (x_ac - x_cg))

Cm  = Cm_0 + Cm_alpha*alpha
delta = -Cm/Cm_delta 

CL  = CL_0 + CL_alpha*alpha + CL_delta*delta
CL_tail = eta*(S_t/S)*CL_alpha_t*(alpha*(1-eps_alpha) + i_t - eps0) + CL_delta*delta
CL_wing = CL_alpha_w*(alpha + i_w-alpha0_w)  
cfe = 0.004
S_fus_vtail = 11.392
S_wet = S_fus_vtail + 2*S + 2*S_t
Cd0 = cfe*S_wet/S
Cdi = CL_wing**2/(np.pi*AR_W*e) +(S_t/S)*CL_tail**2/(np.pi*AR_t*e)
CD = Cdi + Cd0


""" Performance Parameters """
Range = CL/CD
endurance = CL**1.5/CD
W = 150*9.81
rho = 1.225
v_stall = 30*5./18
# lift scale factor due to flap
CL_scale = 1.75
lift_stall = 0.5*rho*v_stall**2*np.max(CL)*CL_scale*S
sink = np.sqrt(2*W/S/rho)/endurance

CL_cruise = CL[np.argmax(endurance)]
Cruise_Vel = np.sqrt(2*W/(rho*S*CL_cruise))

print "Sink rate at landing", sink[-1]
print "Stall speed:", np.sqrt(2*W/(rho*S*CL[-1]*CL_scale))*18./5

""" Desired Tail Setting Angle """


print 'Epsilon 0:',eps0;
Des_Cm_0 = -Cm_alpha*alpha[np.argmax(Range)]
Des_it = eps0 - (Des_Cm_0 - Cm0_w - CL_0*(x_ac - x_cg) - Cm0_f)/(eta*V_h*CL_alpha_t)
print 'Desired Tail angle:', Des_it

""" Printing Output """

print 'Gross weight (empty + payload): ', W
print 'Airfoil used: GOE 441 for Re-500,000 to 1,000,000' 
print '-----------------------------------------------'
print 'AR_W: ', AR_W 
print 'S: ', S
print 'Span: ', b
print 'MAC: ', c
print 'tail volume: ', V_h
print 'tail span: ', b_t
print 'tail surface area: ', S_t
print 'total wetted area:', S_wet
print 'eta (downwash ratio):', 0.9
print 'oswald efficiency factor (e):', 0.92
print '----- Performance Parameters--------------------'
print 'Delta', delta
print 'Delta @ max Endurance', delta[np.argmax(endurance)]
print '------------------------------------------------'
print '(CL/CD)_max: ', np.max(Range), 'at AOA: ', alpha[np.argmax(Range)]
print '(CL^1.5/CD)_max: ', np.max(endurance), 'at AOA: ', alpha[np.argmax(endurance)]
print 'CL_max: ', np.max(CL), 'at AOA: ', alpha[np.argmax(CL)]
print 'CL', CL
print '(CL/CD) at max endurance : ', Range[np.argmax(endurance)] 
print 'Sink Rate (m/s): ', np.min(sink), 'at AOA: ', alpha[np.argmax(endurance)]
print 'Glide slope: ', atan(1/np.max(Range))*180/np.pi, 'at AOA: ', alpha[np.argmax(Range)]
print 'lift at stall: ', lift_stall, lift_stall*1.6/1.75
print 'Cruise Velocity', Cruise_Vel

print 'THE END'
