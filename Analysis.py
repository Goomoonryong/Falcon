import numpy as np
import pdb
import matplotlib.pyplot as plt
from math import atan

""" Design Variables """

alpha_w = -8.2
alpha0_t = 0                
cl_alpha = 0.1                 #""" Airfoils used in wing and Tail """
cl_alpha_t = 2*np.pi**2/180
Cm_0w = -0.2
Cm_af = 0
Cm_f  = 0


b = 10.0
S   = 11.0	
l_t = 4.72          # """ Wing and Tail Dimensions """
AR_W = b*b/S
c = S/b
S_t = 1.432      
b_t = 2
AR_t = b_t**2/S_t
V_h = S_t*l_t/(S*c)
e = 0.92
eta = 0.9

i_w = 0
i_t = -1.98010788483


x_cg = -0.3928/c
x_ac = -c/(4*c)



# ''' Lift and Moment Coefficients '''
CL_alpha_w = cl_alpha/(1 + 57.3*cl_alpha/(np.pi*AR_W*e))
CL_alpha_t = cl_alpha_t/(1 + 57.3*cl_alpha_t/(np.pi*AR_t*e))

eps0 = CL_alpha_w*(i_w - alpha_w)*57.3/(np.pi*AR_W*e)         # Epsilon 0
eps_al = CL_alpha_w*57.3/(np.pi*AR_W*e)                       # d(epsilon)/d(alpha)

CL_alpha = CL_alpha_w + eta*(S_t/S)*CL_alpha_t*(1 - eps_al)
CL_0 = CL_alpha_w*(i_w-alpha_w) + eta*(S_t/S)*CL_alpha_t*(i_t - eps0)
alpha0 	= -CL_0/CL_alpha     # Alpha(L=0) wing and tail

Cm_alpha = CL_alpha*(x_ac - x_cg) - eta*V_h*(1 - eps_al)*CL_alpha_t + Cm_af
Cm_0     = Cm_0w + CL_0*(x_ac - x_cg) - eta*V_h*(i_t - eps0)*CL_alpha_t + Cm_f

x_NP = x_ac - (eta*V_h*(1 - eps_al)*CL_alpha_t + Cm_af)/CL_alpha

SM = 0.1*c
Des_x_NP = x_cg*c - SM
Des_Vh = (x_ac - Cm_af/CL_alpha - Des_x_NP)/(eta*(1 - eps_al)*CL_alpha_t/CL_alpha)
Des_lt = Des_Vh*S*c/S_t
print 'x_NP', x_NP
print 'Desired lt', Des_lt
""" Trim Conditons """

CLt_delta = 0.04
CL_delta = eta*(S_t/S)*CLt_delta
Cm_delta = -CL_delta*(l_t/c - (x_ac - x_cg))


alpha = np.arange(-5,13,0.2)

Cm  = Cm_0 + Cm_alpha*alpha
delta = -Cm/Cm_delta 

CL  = CL_0 + CL_alpha*alpha + CL_delta*delta
CL_tail = eta*(S_t/S)*CL_alpha_t*(alpha*(1-eps_al) + i_t - eps0) + CL_delta*delta
CL_wing = CL_alpha_w*(alpha + i_w-alpha_w)  
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
lift_stall = 0.5*rho*v_stall**2*np.max(CL)*1.75*S
sink = np.sqrt(2*W/S/rho)/endurance

CL_cruise = CL[np.argmax(endurance)]
Cruise_Vel = np.sqrt(2*W/(rho*S*CL_cruise))

print "Sink rate at landing", sink[-1]
print "Stall speed:", np.sqrt(2*W/(rho*S*CL[-1]*1.5))*18./5

""" Desired Tail Setting Angle """


print 'Epsilon 0:',epso;
Des_Cm_0 = -Cm_alpha*alpha[np.argmax(Range)]
Des_it = epso - (Des_Cm_0 - Cm_0w - CL_0*(x_ac - x_cg) - Cm_f)/(eta*V_h*CL_alpha_t)
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
