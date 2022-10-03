from math import *
from re import A
#import matplotlib.pyplot as plt
from tkinter import *
import numpy as np

g = -9.81 # m/s^2

ground = 572.5
pixels_per_meter = 50

payload = 0.5 # kg
counterweight = 29 # kg
major = 1.33 # m
minor = 0.86 # m
trebuchet_height = 1.68 # m, 5.5 ft
lever_mass = 2.1*(major+minor) # kg
theta_0 = -1.047 # rad
theta_f = 1.13 # rad
counterweight_length = 0.5 # m

fulcrum_x, fulcrum_y = 100, ground - trebuchet_height*pixels_per_meter

t_step = 0.0005 # s

payload_x, payload_y = fulcrum_x - major*pixels_per_meter*cos(theta_0), fulcrum_y - major*pixels_per_meter*sin(theta_0)

has_launched = False

t, theta = 0, theta_0
theta_prime = 0
phi_prime = 0
phi = pi/2 + theta_0

using_runge = True #if true, the Runge Kutta method is used to solve the diffeqs, and if not Euler's method is used
def x_range(v_0, theta, y_0):
    t = (-v_0*sin(theta) - sqrt(v_0**2 * sin(theta)**2 - 2*g*y_0)) / g
    return v_0*cos(theta)*t

def get_derivative(state):
    '''diffeqs came from http://www.algobeautytreb.com/trebmath356.pdf'''
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

def update():
    global payload_x, payload_y, payload, counterweight, lever_mass, theta_0, theta_f, major, minor, counterweight_length, trebuchet_height

    payload = float(payload_input.get('1.0','end-1c'))
    counterweight = float(counterweight_input.get('1.0','end-1c'))
    major = float(major_input.get('1.0','end-1c'))
    minor = float(minor_input.get('1.0','end-1c'))
    trebuchet_height = float(trebuchet_height_input.get('1.0','end-1c'))
    lever_mass = float(lever_mass_input.get('1.0','end-1c'))
    theta_0 = float(theta_0_input.get('1.0','end-1c'))
    theta_f = float(theta_f_input.get('1.0','end-1c'))
    counterweight_length = float(counterweight_length_input.get('1.0','end-1c'))

    fulcrum_y = ground - trebuchet_height*pixels_per_meter
    canvas.coords(pedestal, fulcrum_x, fulcrum_y, fulcrum_x+pixels_per_meter/5, fulcrum_y+trebuchet_height*pixels_per_meter, fulcrum_x-pixels_per_meter/5, fulcrum_y+trebuchet_height*pixels_per_meter)

    canvas.coords(major_arm_graph, fulcrum_x,fulcrum_y,minor*pixels_per_meter*cos(theta)+fulcrum_x, fulcrum_y+minor*pixels_per_meter*sin(theta))
    canvas.coords(minor_arm_graph, fulcrum_x, fulcrum_y, fulcrum_x - major*pixels_per_meter*cos(theta), fulcrum_y - major*pixels_per_meter*sin(theta))
    canvas.coords(counterweight_graph, minor*pixels_per_meter*cos(theta)+fulcrum_x + counterweight_length*pixels_per_meter*sin(phi - theta - pi/2)-2,
        fulcrum_y+minor*pixels_per_meter*sin(theta) + counterweight_length*pixels_per_meter*cos(phi - theta - pi/2)-2,
        minor*pixels_per_meter*cos(theta)+fulcrum_x + counterweight_length*pixels_per_meter*sin(phi - theta - pi/2)+2,
        fulcrum_y+minor*pixels_per_meter*sin(theta) + counterweight_length*pixels_per_meter*cos(phi - theta - pi/2)+2,
    )
    canvas.coords(counterweight_arm_graph, minor*pixels_per_meter*cos(theta)+fulcrum_x,
    fulcrum_y+minor*pixels_per_meter*sin(theta),
    minor*pixels_per_meter*cos(theta)+fulcrum_x + counterweight_length*pixels_per_meter*sin(phi - theta - pi/2),
    fulcrum_y+minor*pixels_per_meter*sin(theta) + counterweight_length*pixels_per_meter*cos(phi - theta - pi/2))
    
    if not has_launched:
        payload_x, payload_y = fulcrum_x - major*pixels_per_meter*cos(theta), fulcrum_y - major*pixels_per_meter*sin(theta)

    # if we HAVE launched, payload_x and payload_y will be directly set by the launcher function

    canvas.coords(payload_graph, payload_x-2,
        payload_y-2,
        payload_x+2,
        payload_y+2)

    window.update()
    canvas.grid(row=0, column=0)

def reload():
    global t, theta, theta_prime, phi, phi_prime, has_launched, distance_label, efficiency_label
    has_launched = False
    t, theta = 0, theta_0
    theta_prime = 0
    phi_prime = 0
    phi = pi/2 + theta_0

    distance_label.config(text='')
    efficiency_label.config(text='')

    update()

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

        theta += derivs[0] * t_step #misaligned definitions of theta between us and http://www.algobeautytreb.com/trebmath356.pdf
        phi += derivs[1] * t_step
        
        theta_prime += derivs[2]*t_step
        phi_prime += derivs[3]*t_step

        update()

        t += t_step
        
        if t%0.02<0.0001:
            print("Current velocity: "+str(abs(theta_prime*major))+"m/s")
    print(t)
    print("Final velocity: {}".format(abs(theta_prime*major))+"m/s")
    change_payload_kin = 1/2 * payload * (theta_prime * major)**2
    change_cw_pot = counterweight * -g * (minor*sin(-theta_0) - counterweight_length \
        - (-minor*sin(theta) - counterweight_length*cos(phi - theta - pi/2)))

    efficiency = change_payload_kin / change_cw_pot
    efficiency_label.config(text=str(efficiency*100)[:5])

    print("--------------------------------")
    has_launched = True
    v_0 = -theta_prime * major * pixels_per_meter
    v_x, v_y = v_0*sin(theta), -v_0*cos(theta) # switch sin and cos cuz theta is not our lanch angie; it was lever angle which is complement of what we want
    tru_distance = x_range(v_0/pixels_per_meter,theta,trebuchet_height+sin(theta)*major)
    distance_label.config(text = str(tru_distance)[:6])

    while payload_y < 570.5 and has_launched:
        payload_x += v_x*t_step
        v_y += -(g*pixels_per_meter) * t_step # g is negative cuz top-bottom is reversed
        payload_y += v_y * t_step
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

        update()

        t += t_step

    

window = Tk()
window.title('Trebuchet Simulator')
canvas = Canvas(window, width=1000, height=1000, highlightbackground='black')
canvas.grid(row=0, column=0, rowspan=100)

ground_graph = canvas.create_line(0, ground, 1000, ground, fill='green')
minor_arm_graph = canvas.create_line(fulcrum_x,fulcrum_y,minor*pixels_per_meter*cos(theta)+fulcrum_x, fulcrum_y+minor*pixels_per_meter*sin(theta))
major_arm_graph = canvas.create_line(fulcrum_x, fulcrum_y, fulcrum_x - major*pixels_per_meter*cos(theta), fulcrum_y - major*pixels_per_meter*sin(theta))
counterweight_arm_graph = canvas.create_line(minor*pixels_per_meter*cos(theta)+fulcrum_x,
    fulcrum_y+minor*pixels_per_meter*sin(theta),
    minor*pixels_per_meter*cos(theta)+fulcrum_x + counterweight_length*pixels_per_meter*sin(phi - theta - pi/2),
    fulcrum_y+minor*pixels_per_meter*sin(theta) + counterweight_length*pixels_per_meter*cos(phi - theta - pi/2))
payload_graph = canvas.create_oval(payload_x-2,
    payload_y-2,
    payload_x+2,
    payload_y+2,
    fill='red')
counterweight_graph = canvas.create_oval(minor*pixels_per_meter*cos(theta)+fulcrum_x + counterweight_length*pixels_per_meter*sin(phi - theta - pi/2)-2,
    fulcrum_y+minor*pixels_per_meter*sin(theta) + counterweight_length*pixels_per_meter*cos(phi - theta - pi/2)-2,
    minor*pixels_per_meter*cos(theta)+fulcrum_x + counterweight_length*pixels_per_meter*sin(phi - theta - pi/2)+2,
    fulcrum_y+minor*pixels_per_meter*sin(theta) + counterweight_length*pixels_per_meter*cos(phi - theta - pi/2)+2,
    fill='black'
)

pedestal = canvas.create_polygon(fulcrum_x, fulcrum_y, fulcrum_x+pixels_per_meter/5, fulcrum_y+trebuchet_height*pixels_per_meter, fulcrum_x-pixels_per_meter/5, fulcrum_y+trebuchet_height*pixels_per_meter)

person = canvas.create_line(25, ground, 25, ground - pixels_per_meter*1.77)

launch_button = Button(window, command=launch_trebuchet, text='Fire!')
launch_button.grid(row=11, column=1,columnspan=2)

reload_button = Button(window, command=reload, text='Reload')
reload_button.grid(row=12,column=1,columnspan=2)

payload_label = Label(window, text='Payload mass (kg):')
payload_label.grid(row=2, column = 1)
payload_input = Text(window, width=10,height=1)
payload_input.insert(END, payload)
payload_input.grid(row=2, column=2)

counterweight_label = Label(window, text='Counterweight mass (kg):')
counterweight_label.grid(row=3, column = 1)
counterweight_input = Text(window, width=10,height=1)
counterweight_input.insert(END, counterweight)
counterweight_input.grid(row=3, column=2)

major_label = Label(window, text='Left side of lever length (m):')
major_label.grid(row=4, column = 1)
major_input = Text(window, width=10,height=1)
major_input.insert(END, major)
major_input.grid(row=4, column=2)

minor_label = Label(window, text='Right side of lever length (m):')
minor_label.grid(row=5, column = 1)
minor_input = Text(window, width=10,height=1)
minor_input.insert(END, minor)
minor_input.grid(row=5, column=2)

trebuchet_height_label = Label(window, text='Trebuchet height (m):')
trebuchet_height_label.grid(row=6, column = 1)
trebuchet_height_input = Text(window, width=10,height=1)
trebuchet_height_input.insert(END, trebuchet_height)
trebuchet_height_input.grid(row=6, column=2)

lever_mass_label = Label(window, text='Lever mass (kg):')
lever_mass_label.grid(row=7, column = 1)
lever_mass_input = Text(window, width=10,height=1)
lever_mass_input.insert(END, lever_mass)
lever_mass_input.grid(row=7, column=2)

theta_0_label = Label(window, text='Initial angular displacement (rad):')
theta_0_label.grid(row=8, column = 1)
theta_0_input = Text(window, width=10,height=1)
theta_0_input.insert(END, theta_0)
theta_0_input.grid(row=8, column=2)

theta_f_label = Label(window, text='Final angular displacement (rad):')
theta_f_label.grid(row=9, column = 1)
theta_f_input = Text(window, width=10,height=1)
theta_f_input.insert(END, theta_f)
theta_f_input.grid(row=9, column=2)

counterweight_length_label = Label(window, text='Counterweight length (m):')
counterweight_length_label.grid(row=10, column = 1)
counterweight_length_input = Text(window, width=10,height=1)
counterweight_length_input.insert(END, counterweight_length)
counterweight_length_input.grid(row=10, column=2)

Label(window,text='Horizontal Distance (m):').grid(row=13,column=1)
distance_label = Label(window, text='')
distance_label.grid(row=13, column = 2)

Label(window, text='Efficiency (%):').grid(row=14,column=1)
efficiency_label = Label(window, text='')
efficiency_label.grid(row=14, column=2)

window.mainloop()
