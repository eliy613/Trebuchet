from math import *
from re import A
from tkinter import *
import time
import numpy as np

def reinitialize():
    global using_runge, g, ground, pixels_per_meter, payload, counterweight, major, minor, trebuchet_height, lever_mass, theta_0, theta_f, counterweight_length, mu, fulcrum_x, fulcrum_y, t_step
    global payload_x, payload_y, has_launched, theta, t, theta_prime, phi_prime, phi
    
    g = -9.81 # m/s^2
    
    ground = 572.5
    pixels_per_meter = 50

    payload=0.5 # kg
    counterweight=29 # kg
    major=1.33# m
    minor=0.85 # m
    trebuchet_height=1.45 # m, 4.76 ft
    
    theta_0 = -pi/3 # rad
    theta_f= 1.3# rad
    counterweight_length = 0.5 # m
    mu = 2.1 #mass per unit length of rod
    lever_mass=(major+minor)*mu # kg
    
    fulcrum_x, fulcrum_y = 100, ground - trebuchet_height*pixels_per_meter

    t_step=0.0001 # s

    payload_x, payload_y = fulcrum_x - major*pixels_per_meter*cos(theta_0), fulcrum_y - major*pixels_per_meter*sin(theta_0)

    has_launched = False

    t, theta = 0, theta_0
    theta_prime = 0
    phi_prime = 0
    phi = pi/2 + theta_0
    using_runge = True

reinitialize()

parameters = np.array([major,minor,theta_0,theta_f,counterweight_length])
derivs =      np.array([0,   0,    0,      0,       0,                  ], dtype="float64") #trigger while loop
dv = 0.001
lr=0.0001

def launch_trebuchet():
    # POSITIVE is CLOCKWISE
    global theta, theta_prime, phi, phi_prime, t, has_launched, payload_x, payload_y, trebuchet_height, theta_0, major, minor, payload, counterweight, lever_mass, counterweight_length

    while theta_f - theta > t_step:
        state = np.array([theta, phi, theta_prime, phi_prime])
        derivs = get_derivative([theta, phi, theta_prime, phi_prime]) # t --> theta; p--> phi
        if using_runge:
            k2 = get_derivative(state+derivs*t_step/2)

            k3 = get_derivative(state+k2*t_step/2)

            k4 = get_derivative(state+k3*t_step)

            derivs = 1/6 * (derivs+2*k2+2*k3+k4)

        theta += derivs[0] * t_step #misaligned definitions of theta
        phi += derivs[1] * t_step
        
        theta_prime += derivs[2]*t_step
        phi_prime += derivs[3]*t_step


        t += t_step
    
    has_launched = True
    v_0 = -theta_prime * major * pixels_per_meter
    tru_distance = x_range(v_0/pixels_per_meter,theta,trebuchet_height+sin(theta)*major)
    return tru_distance

def x_range(v_0, theta, y_0):
    t = (-v_0*sin(theta) - sqrt(v_0**2 * sin(theta)**2 + 2*g*(-1)*y_0)) / g
    return v_0*cos(theta)*t

def get_derivative(state):
    '''solved diffeq came from http://www.algobeautytreb.com/trebmath356.pdf'''
    global has_launched
    
    g = 9.81

    theta = state[0]; phi = state[1]; theta_prime= state[2]; phi_prime = state[3]
    
    theta = pi/2 - theta #we define theta differently, rip
    
    m2 = payload if not has_launched else 0.
    m1 = counterweight
    l2 = major
    l1 = minor
    l4 = counterweight_length
    mb = lever_mass
    # theta'' is x, phi'' is y
    a = -l1*l4*m1*2*cos(phi) + l1**2*m1 + l4**2*m1 + l2**2*m2 + 1/3*l1**2*mb - 1/3*l1*l2*mb + 1/3*l2**2*mb
    b = -l1*l4*m1*cos(phi) + l4**2*m1
    c = l1*l4*m1*phi_prime*(phi_prime+2*theta_prime)*sin(phi) + g*l1*m1*sin(theta) - g*l2*m2*sin(theta) + 1/2*g*l1*mb*sin(theta) - 1/2*g*l2*mb*sin(theta) - g*l4*m1*sin(phi+theta)

    d = l4*m1*l4 - l4*m1*l1*cos(phi)
    e = l4*m1*l4
    f = l4*m1*-l1*theta_prime**2*sin(phi) - l4*m1*g*sin(phi+theta)

    left = np.array([[a, b],[d, e]])
    right = np.array([-c, -f])

    theta_double_prime, phi_double_prime = np.linalg.solve(left, right)
    return np.array([-theta_prime, phi_prime, theta_double_prime, phi_double_prime])

def simulate():
    global theta, theta_prime, phi, phi_prime, t, has_launched, payload_x, payload_y, trebuchet_height, theta_0, major, minor, payload, counterweight, lever_mass, counterweight_length, theta_0, theta_f
    global parameters
    major = parameters[0]
    minor = parameters[1]
    theta_0 = parameters[2]
    theta_f = parameters[3]
    counterweight_length = parameters[4]
    return launch_trebuchet()

'''We will keep on optimizing until we get a distance of over 35m. The optimization algorithm
would keep on going, but the results just become nonsensical because the differences between
the simulation and real life build up and our idealizing assumptions (having a point mass for the
hinged counterweight, no friction, etc) become incredibly untrue.'''
#parameters are 
curDist = 0
while curDist<35: 
    print("\n\n\n---------------------------")
    reinitialize()
    lever_mass = mu*(major+minor)

    curDist = simulate()
    print("With these parameters: {} \n We get a distance of: {}".format(parameters,curDist))
    print("\n")

    for i in range(5): #numpy array is stupid
        
        reinitialize()
        parameters[i] += dv
        newDist = simulate()
        derivs[i] = (newDist - curDist)/dv
        print(newDist)
        parameters[i] -= dv
    if np.max(derivs) > 1000:
        print("The derivatives are too high, the optimization at this point is nonsense. Try changing the learning rate or restarting with different parameters")
        break
    print("Derivs: ",derivs)
    parameters += lr*derivs
    parameters[0] = max(parameters[0],0.05); parameters[1] = max(parameters[1],0.05)
    
print("major: {}, minor: {}, theta_0: {}, theta_f: {}, counterweight_length: {}".format(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]))

    
        
