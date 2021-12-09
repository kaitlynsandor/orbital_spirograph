import PySimpleGUI as sg
import math
import turtle
import matplotlib.pyplot as plt
from numpy import pi, exp, real, imag, linspace
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from astropy.time import Time
from astroquery.jplhorizons import Horizons

sim_start_date = "2021-01-01"     # simulating a solar system starting from this date
sim_duration = 3 * 365                # (int) simulation duration in days
G = 6.67e-11
G_AU = 1.488e-34

def animate(i):
    return ss.evolve()

class Object:                   # define the objects: the Sun, Earth, Mercury, etc
    def __init__(self, name, rad, color_line, color_spiro, z, r):
        self.name = name
        self.r    = np.array(r, dtype=np.float)
        self.v    = np.array([0,0,0], dtype=np.float)
        self.xs = []
        self.ys = []
        self.line_color = color_line
        self.spiro_xs = []
        self.spiro_ys = []
        self.period = 0
        self.plot = ax.scatter(r[0], r[1], color=color_line, s=rad**2.7, edgecolors=None, zorder=z)
        self.line, = ax.plot([], [], color=color_line, linewidth=3.4)
        self.spiro_line, = ax.plot([], [], color=color_spiro, linewidth=0.7)


class SolarSystem:
    def __init__(self, thesun, mass):
        self.thesun = thesun
        self.planets = []
        self.center_mass = mass
        self.time = None
        self.count = 0

    def planet_report(self):
        length = 70
        beg = "PLANET REPORT"
        print(beg.center(length, '*'))
        print()

        for planet in self.planets:
            name = planet.name
            print(name.center(length, '-'))

            print('Semi Major Axis (AU): ' + str(planet.r[0]))
            print('Initial Velocity (X, Y) (AU/Day): ' + str(planet.v[0]) + ', ' + str(planet.v[1]))
            print('Calculated Period of Planet (Days): ' + str(planet.period))
            print()

    def AU_m(self, distance):
        distance = math.sqrt(distance[0]**2 + distance[1]**2) # get distance to SUN
        return distance*1.5e11 #convert AU to m

    def AU_day(self, acc):
        return acc*5.775e-7

    def calculate_initial_velocity(self):
        for planet in self.planets:
            y = math.sqrt(G_AU*self.center_mass/np.sum(planet.r ** 2)** (1/2))
            planet.v = np.array([0,y], dtype=np.float)

    def calculate_period(self):
        for planet in self.planets:
            result = math.sqrt(4*(3.14**2)/G_AU/self.center_mass*(np.sum(planet.r ** 2)** (3/2)))
            planet.period = result


    def calculate_acceleration(self, planet):
        acc = -G_AU * self.center_mass * planet.r / np.sum(planet.r ** 2)** (3/2)
        return acc

    def add_planet(self, planet):
        self.planets.append(planet)

    def evolve(self):           # evolve the trajectories
        dt = 1
        self.time += dt
        plots = []
        lines = []
        labels = []


        for i in range(len(self.planets)):
            p = self.planets[i]
            p.r[0] += p.v[0] * dt
            p.r[1] += p.v[1] * dt

            acc = self.calculate_acceleration(p)  # in units of AU/day^2

            p.v += acc * dt

            p.xs.append(p.r[0])
            p.ys.append(p.r[1])

            p.plot.set_offsets(p.r[:2])
            p.line.set_xdata(p.xs)
            p.line.set_ydata(p.ys)

            plots.append(p.plot)

            if len(p.xs) > 1 and i > 0 and self.count == 7:
                p.spiro_xs.append(self.planets[i - 1].xs[-1])
                p.spiro_ys.append(self.planets[i-1].ys[-1])
            else:
                p.spiro_xs.append(p.xs[-1])
                p.spiro_ys.append(p.ys[-1])

            p.spiro_line.set_xdata(p.spiro_xs)
            p.spiro_line.set_ydata(p.spiro_ys)

            lines.append(p.spiro_line)
            lines.append(p.line)

        if self.count == 8:
            self.count = 0
        else:
            self.count += 1
        lines = lines[::-1]
        plots.insert(0, self.thesun.plot)
        return labels + plots + lines

# deal with main driver of the program to configure animation
if __name__ == "__main__":
    sg.theme('DefaultNoMoreNagging')  # Add a touch of color
    # All the stuff inside your window.
    layout = [
              [sg.Text('How many planets would you like to simulate?'), sg.InputText()],
        [sg.Text('What is the Mass of the Center Mass/Sun (kg)'), sg.InputText()],
              [sg.Button('Enter'), sg.Button('Quit')]]

    # Create the Window
    window = sg.Window('Spirograph Orbit Simulation Setup', layout)
    new_layout = []
    planets = 4
    planets_to_data = {}
    solar_mass = 0

    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == 'Quit':  # if user closes window or clicks cancel
            break
        elif event == 'Enter':
            if values[0].isnumeric() and 6 > int(values[0]) > 1 and values[1].isnumeric():
                solar_mass = int(values[1])
                planets = int(values[0])
                new_layout.append([sg.Text('SIMULATING '+ str(values[0]) + ' PLANETS...')])
                new_layout.append([sg.Text('DEFAULT SIMULATION EACH COMPARED WITH PLANET BELOW IT')])
                for i in range(int(values[0])):
                    new_layout.append([[sg.Text('Name of planet'), sg.InputText(), sg.Text('Length of Semi-Major Axis (AU)'), sg.InputText(), sg.Text('Color of Spirograph'), sg.InputText()]])
                new_layout.append([sg.Button('Enter'), sg.Button('Quit')])
                break
            elif values[1].isnumeric():
                sg.Popup('You must enter a positive integer <= 5 and >= 2', keep_on_top=True)
            else:
                sg.Popup('Please enter a solar mass', keep_on_top=True)
    window.close()

    window = sg.Window('Spirograph Orbit Simulation Setup', new_layout)

    # Event Loop to process "events" and get the "values" of the inputs
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == 'Quit':  # if user closes window or clicks cancel
            window.close()
            break
        elif all(x != '' for x in values.values()):
            offset = 0
            for i in range(planets):
                planets_to_data[values[offset]] = [float(values[offset+1]), (values[offset+2])] # dict --> planet:{semi, color}
                offset += 3
            break
        else:
            sg.Popup('You must enter a value for each box', keep_on_top=True)

    plt.style.use('dark_background')
    fig = plt.figure(figsize=[6, 6])
    ax = plt.axes([0., 0., 1., 1.], xlim=(-1.8, 1.8), ylim=(-1.8, 1.8))
    ax.set_aspect('equal')
    ax.axis('off')
    ss = SolarSystem(Object("Sun", 12, 'yellow', 'yellow', 15, [0, 0, 0]), solar_mass) # <-- mass of sun hard coded
    ss.time = Time(sim_start_date).jd

    for planet in planets_to_data.keys():
        ss.add_planet(Object(planet, 10, planets_to_data[planet][1], planets_to_data[planet][1], 4,
                         [planets_to_data[planet][0], 0]))

    ss.calculate_initial_velocity()
    ss.calculate_period()
    ss.planet_report()
    ani = animation.FuncAnimation(fig, animate, repeat=True, frames=sim_duration, blit=True, interval=20,)
    plt.show()
