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
        self.rho = 1.225
        self.g = 9.81
        self.cx = 0.4
        self.m = 0.2
        self.a = np.pi*(radius**2)
        self.dt = time_step
        self.sim_time = duration
        self.iterations = int(self.sim_time/self.dt)
        self.t = np.linspace(0, self.sim_time, self.iterations)
        self.vx = [float(v_init*np.cos(angle))]
        self.vy = [float(v_init*np.sin(angle))]
        self.v_total = []
        if not isinstance(origin, list):
            raise TypeError("origin must be a list of two floats")
        self.origin = origin
        self.path = {'x': [origin[0]],
                     'y': [origin[1]]}
        self.d = 0

    def dvx_dt(self, vx):
        """
        :param vx: Velocity in the x direction
        :return: Accelaration in the x direction
        """
        return (-0.5*self.cx*self.rho*self.a*vx*vx*np.sign(vx))/self.m

    def dvy_dt(self, vy):
        """
        :param vy: Velocity in the y direction
        :return: Acceleration in the y direction
        """
        return ((-0.5*self.cx*self.rho*self.a*vy*vy*np.sign(vy))/self.m) - self.g

    def euler(self):
        """
        :return: Aproxiamtion of velocity in the x and y direction using Euler's method
        """
        for i in range(1, len(self.t)):
            self.vx.append(float(self.vx[i - 1] + self.dt * self.dvx_dt(self.vx[i - 1])))
            self.vy.append(float(self.vy[i - 1] + self.dt * self.dvy_dt(self.vy[i - 1])))

    def velocity(self):
        """
        :return: Velocity of the rocket in m/s
        """
        self.euler()
        for i in range(len(self.vx)):
            self.v_total.append(float(np.sqrt(self.vx[i]*self.vx[i] + self.vy[i]*self.vy[i])))

    def plot_velocity(self):
        """
        :return: Plot of the velocity of the rocket in m/s
        """
        self.velocity()
        plt.plot(self.t, self.v_total, 'b', label="Total Velocity")
        plt.plot(self.t, self.vx, 'r--', label="X Velocity")
        plt.plot(self.t, self.vy, 'g--', label="Y Velocity")
        plt.xlabel("Time")
        plt.ylabel("Velocity")
        plt.legend()
        plt.grid()
        plt.title("Velocity vs Time")
        plt.axis((0, self.sim_time, min(self.vy)*1.2, max(self.v_total) * 1.1))
        plt.show()

    def set_path(self):
        """
        :return path: Path of the rocket
        :return d: Distance of the rocket's landing in meters
        """
        self.velocity()
        for i in range(1, len(self.t)):
            self.path['x'].append(self.path['x'][i-1] + self.dt*self.vx[i-1])
            self.path['y'].append(self.path['y'][i-1] + self.dt*self.vy[i-1])
        for i in range(len(self.path['y'])):
            if self.path['y'][i] > 0:
                self.d = self.path['x'][i]
                i = len(self.path['y'])
        print(f'the rocket will land at d = {self.d}')
        return self.path, self.d

    def plot_path(self):
        """
        :return: Plot of the path of the rocket
        """
        self.set_path()
        plt.plot(self.path['x'], self.path['y'], 'b', label="Flight Path")
        plt.title("Free Flight Path")
        plt.xlim(0, int(self.d*1.1)+1)
        plt.ylim(0, int(max(self.path['y'])*1.1)+1)
        plt.gca().set_aspect('equal')
        plt.ylabel("Height in m")
        plt.xlabel("Distance in m")
        plt.legend()
        plt.tight_layout()
        plt.savefig('Path.png', dpi=1000)
        plt.show()


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
        self.dt = dt
        self.m_rocket = [m_rocket + m_water]
        self.m_water = [m_water]
        self.p_init = p_initial*100000
        self.r_neck = r_neck/100
        self.v_bottle = v_bottle/1000
        self.angle = initial_angle
        self.neck_a = np.pi*self.r_neck**2
        self.ax = []
        self.ay = []
        self.vx = [0]
        self.vy = [0]
        self.rho = 1000
        self.g = 9.81
        self.v_water_init = self.m_water[0] / self.rho
        self.p_final = self.p_init * ((self.v_bottle - self.v_water_init) / self.v_bottle)
        self.v_e = np.sqrt(2*(self.p_init-self.p_final)/self.rho)  # From Bernoulli's equation
        self.m_dot = self.rho * self.neck_a * self.v_e
        self.path = {'x': [0],
                     'y': [0]}

    def force(self):
        '''
        :return: The force of the rocket in N
        '''
        return self.m_dot*self.v_e

    def time_of_launch(self):
        '''
        :return: The time of the launch of the rocket in seconds
        '''
        t = self.m_water[0]/(self.v_e*self.neck_a*self.rho)
        return t

    def mass_flow(self):
        '''
        :return: The mass flow of the rocket in kg/s
        '''
        if len(self.m_rocket) == 1:
            # print('Calculating mass')
            # print(self.m_water)
            while self.m_water[len(self.m_water)-1] > 0:
                self.m_rocket.append(float(self.m_rocket[len(self.m_rocket)-1] - (self.dt * self.m_dot)))
                self.m_water.append(float(self.m_water[len(self.m_water)-1] - (self.dt * self.m_dot)))
            # print(self.m_rocket)
            # print(self.m_water)
        else:
            print('Mass already calculated')

    def dvx_dt(self):
        '''
        :return: The derivative of the velocity in the x direction of the rocket in m/s^2
        '''
        # print('Calculating Acceleration on X')
        self.mass_flow()
        for i in range(len(self.m_rocket)):
            self.ax.append(float((self.force() * np.cos(self.angle)) / self.m_rocket[i]))


    def dvy_dt(self):
        '''
        :return: The derivative of the velocity in the y direction of the rocket in m/s^2
        '''
        # print('Calculating Acceleration on Y')
        self.mass_flow()
        for i in range(len(self.m_rocket)):
            self.ay.append(float((self.force()*np.sin(self.angle))/self.m_rocket[i]) - self.g)

    def euler_velocity(self):
        '''
        :return: The velocity of the rocket calculated with euler's method
        '''
        # print('Calculating Velocities of take-off')
        self.dvx_dt()
        self.dvy_dt()
        for i in range(1, len(self.m_rocket)):
            self.vx.append(self.vx[i - 1] + self.dt*self.ax[i-1])
            self.vy.append(self.vy[i - 1] + self.dt*self.ay[i-1])
        # print('Velocities calculated', '\n')
        return {'vx': self.vx,
                'vy': self.vy}

    def euler_path(self):
        '''
        :return: The path of the rocket using Euler's method
        '''
        # print('Calculating Path of take-off')
        self.euler_velocity()
        for i in range(1, len(self.vx)):
            self.path['x'].append(self.path['x'][i - 1] + self.dt*self.vx[i])
            self.path['y'].append(self.path['y'][i - 1] + self.dt*self.vy[i])
        # print('Path calculated', '\n')
        return self.path

    def plot(self):
        '''
        :return: Plots the graph of the take-off phase
        '''
        self.euler_path()
        plt.plot(self.path['x'], self.path['y'], 'c')
        plt.show()

def calc_distance(p_init, angle_init):

    print(f'calculatig for P={p_init} and Angle={np.degrees(angle_init)}')

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

    # Simulating free flight with results from take-off
    test2 = FreeFlight(final_dest, v_final, angle_init)
    path2, distance = test2.set_path()
    v2x, v2y = test2.vx, test2.vy

    return path1, path2, distance


# Start of full flight calculations
# List of multiple angles and pressures for simulations
list_angles = np.arange(5, 90, 5)
list_pressure = np.arange(1.5, 6, 0.5)

# print(list_angles)
# print(list_pressure)

def calc_all_distances():
    local_res = [[] for _ in range(len(list_angles))]
    local_results = []

    for i, angle in enumerate(list_angles):
        for pressure in list_pressure:
            distance = 0
            _, _, distance = calc_distance(pressure, np.radians(angle))
            local_results.append(f"{pressure} & {angle} & {round(distance, 1)}")
            local_res[i].append(distance)

    return local_res, local_results


def plot_full_path(path1, path2, distance):
    '''
    :return: Plots the full path of the rocket using the variables and paths defined previously
    '''
    plt.plot(path1['x'], path1['y'], 'green', label='Take-off')
    plt.plot(path2['x'], path2['y'], 'orange', label='Free Flight')
    plt.title("Flight Path")
    plt.legend()
    plt.xlim(0, int(distance) + 1)
    plt.ylim(0, int(max(path2['y'])) + 2)
    plt.xlabel("Distance in m")
    plt.ylabel("Height in m")
    plt.gca().set_aspect('equal')
    # plt.savefig('FullPath.png', dpi=200)
    # plt.show()


def plot_all_distances():
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

    x, y = np.meshgrid(list_pressure, list_angles)
    z = np.array(res)
    surf = ax.plot_surface(x, y, z, cmap='plasma')
    ax.set_title('Distance of launch for any angle or pressure')
    ax.set_xlabel('Pressure (Bar)')
    ax.set_ylabel(r'Angle ($^{\circ}$)')
    ax.set_zlabel('Distance (m)')

    fig.colorbar(surf, ax=ax)
    plt.savefig('AllDistances.png')
    plt.show()

def contour():
    a_index = np.where(list_angles == 45)
    print(a_index[0])
    p_index = np.where(list_pressure == 5)
    print(p_index)
    color = 'midnightblue'

    fig1, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(15, 6))
    x, y = np.meshgrid(list_pressure, list_angles)
    print(res)
    z = np.array(res)
    contour = ax1.contourf(x, y, z, cmap='plasma', levels=26)
    ax1.set_title('Contour Plot of Launch Distance')
    ax1.set_xlabel('Pressure (Bar)')
    ax1.set_ylabel(r'Angle ($^{\circ}$)')
    ax1.axvline(x=list_pressure[p_index], color='red', linestyle='--')
    ax1.axhline(y=list_angles[a_index], color='red', linestyle='--')

    cbar = fig1.colorbar(contour, ax=ax1)
    cbar.set_label('Distance (m)', rotation=270, labelpad=15)

    list_ax2 = z[a_index, :]
    list_ax2 = np.reshape(list_ax2, np.shape(list_pressure))
    ax2.plot(list_pressure, list_ax2, color=color)
    ax2.set_title(f'Distance for an angle of {list_angles[a_index]}'+r'$^{\circ}$')
    ax2.set_ylabel('Distance (m)')
    ax2.set_xlabel('Pressure (Bar)')

    list_ax3 = z[:, p_index]
    list_ax3 = np.reshape(list_ax3, np.shape(list_angles))
    ax3.plot(list_angles, list_ax3, color=color)
    ax3.set_title(f'Distance for a pressure of {list_pressure[p_index]} Bar')
    ax3.set_xlabel('Angle '+r'($^{\circ}$)')

    z_max = np.max(z)
    ax2.set_ylim([0, z_max])
    ax3.set_ylim([0, z_max])

    plt.savefig(f'contours.png')
    plt.tight_layout()
    plt.show()


# Asking for the necessary input variables for the simulation
# p_init = float(input("What's the initial pressure?"))
# angle_init = float(input("What's the initial angle?"))*(np.pi/180)
#
# path1, path2, distance = calc_distance(p_init, angle_init)
# plot_full_path(path1=path1, path2=path2, distance=distance)

res, results = calc_all_distances()
plot_all_distances()
contour()

