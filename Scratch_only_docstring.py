import matplotlib.pyplot as plt
import matplotlib
import numpy as np
matplotlib.use("TkAgg")

class FreeFlight:
    def __init__(self, origin: list, v_init, angle, time_step=0.01, duration=10, radius=0.1):
        '''
        :param origin: Origin of the free flight as a list of two elements
        :param v_init: Initial velocity of the rocket in m/s
        :param angle: Angle of the rocket at launch in degrees
        :param time_step: Time interval between each step of the rocket for the Euler method in seconds
        :param duration: Duration of the simulation in seconds
        :param radius: Radius of the rocket's body
        '''

    def dvx_dt(self, vx):
        """
        :param vx: Velocity in the x direction
        :return: Accelaration in the x direction
        """

    def dvy_dt(self, vy):
    """
    :param vy: Velocity in the y direction
    :return: Acceleration in the y direction
    """

    def euler(self):
        """
        :return: Aproxiamtion of velocity in the x and y direction using Euler's method
        """

    def velocity(self):
        """
        :return: Velocity of the rocket in m/s
        """

    def plot_velocity(self):
        """
        :return: Plot of the velocity of the rocket in m/s
        """

    def set_path(self):
        """
        :return path: Path of the rocket
        :return d: Distance of the rocket's landing in meters
        """

    def plot_path(self):
        """
        :return: Plot of the path of the rocket
        """

class TakeOff:
    def __init__(self, m_rocket, r_neck, m_water, p_initial, v_bottle, initial_angle, dt=0.0001):
        '''
        :param m_rocket: Mass of the dry rocket in kg
        :param r_neck: Radius of the bottle neck in cm
        :param m_water: Mass of the water in kg
        :param p_initial: Initial Pressure in Bar
        :param v_bottle: Volume of the bottle in liters
        :param initial_angle: Initial angle of the rocket in degrees
        :param dt: Time step of the simulation
        '''


    def force(self):
        '''
        :return: The force of the rocket in N
        '''

    def time_of_launch(self):
        '''
        :return: The time of the launch of the rocket in seconds
        '''

    def mass_flow(self):
        '''
        :return: The mass flow of the rocket in kg/s
        '''

    def dvx_dt(self):
        '''
        :return: The derivative of the velocity in the x direction of the rocket in m/s^2
        '''

    def dvy_dt(self):
        '''
        :return: The derivative of the velocity in the y direction of the rocket in m/s^2
        '''

    def euler_velocity(self):
        '''
        :return: The velocity of the rocket calculated with euler's method
        '''

    def euler_path(self):
        '''
        :return: The path of the rocket using Euler's method
        '''

    def plot(self):
        '''
        :return: Plots the graph of the take-off phase
        '''

# Asking for the necessary input variables for the simulation
p_init = float(input("What's the initial pressure?"))
angle_init = float(input("What's the initial angle?"))*(np.pi/180)

# Simulating take-off
test1 = TakeOff(.2, 1, .5, p_init, 1.5, angle_init)

# Storing the necessary data from take-off (velocities)
test1.euler_path()
v1x = test1.vx
v1y = test1.vy
v_final = np.sqrt(max(v1x) ** 2 + max(v1y) ** 2)

# Storing the necessary data from take-off (path of take-off, final position)
path1 = test1.path
final_dest = [max(path1['x']), max(path1['y'])]
t = test1.time_of_launch()
tt = np.arange(0, t, len(v1x))
angle = np.arctan(max(v1y)/max(v1x))

# Simulating free flight with results from take-off
test2 = FreeFlight(final_dest, v_final, angle)
path2, distance = test2.set_path()
v2x, v2y = test2.vx, test2.vy

def plot_full_path():
    '''
    :return: Plots the full path of the rocket using the variables and paths defined previously
    '''

plot_full_path()
