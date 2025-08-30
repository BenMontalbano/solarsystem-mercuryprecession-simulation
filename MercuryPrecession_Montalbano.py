import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# Author: Ben MontalbanoS
# Purpose: The purpose of this script is to make a computational simulation
# of the Earth falling into the sun.

# Define Gravitational Constant, speed of light, and step size
GC = 6.6743E-11  # SI units
c = 299792458 # meters/second
TIMEINC = 10  # secs

# Create Storage systems for data of bodies
class body:
    def __init__(self, name, mass, Xpos, Ypos, Zpos, Xvel, Yvel, Zvel):
        self.name = name
        self.mass = mass  # kg
        self.Xpos = Xpos  # m
        self.Ypos = Ypos  # m
        self.Zpos = Zpos
        self.Xvel = Xvel  # m/s
        self.Yvel = Yvel  # m/s
        self.Zvel = Zvel  # m/s


#Input Data for Bodies
Sun = body('Sun', 1.989E30, 0, 0, 0, 0, 0, 0)
Mercury = body('Mercury', 3.3E23, 7E10, 0, 0, 0, 39000, 0) 
Venus = body('Venus', 4.87E24, 0, 1.08E11,  0, -35020, 0, 2196.79) 
Earth = body('Earth', 5.97E24, -1.4964E11, 0, 0, 0, -(2 * np.pi * 1.4964E11) / (3.15E7), 0)
Mars = body('Mars', 6.42E23, 2.28E11, 0, 0, 0, 24070, 777.67)  
Jupiter = body('Jupiter', 1.8982E27, 0, -7.78E11, 0, 13070, 0, 559.99) 

# Create List of Bodies
Bodies=[Sun,Mercury]
colors = ['yellow','orange','blue','red', 'pink']
Num_Bodies = len(Bodies)


# Create Arrays for Newobject with Distances to each Body
def Distances(Newobject):
    
    DistanceArray = np.zeros((Num_Bodies, 1), dtype=np.float64)
    i = 0
    for Body in Bodies:
        RadSquare = ((Newobject.Xpos - Body.Xpos) ** 2 +
                     (Newobject.Ypos - Body.Ypos) ** 2 +
                     (Newobject.Zpos - Body.Zpos) ** 2)
        DistanceArray[i] = RadSquare
        i += 1
    return DistanceArray

# Find Gravitational Force by each body
def Gravity(Newobject):
    DistanceArray = Distances(Newobject)
    i = 0
    GForce = np.zeros((3, 1), dtype=np.float64)

    for Body in Bodies:
        if Body == Newobject:
            i += 1
            continue
        
        if Newobject == Mercury and Body == Sun:
            rVector = np.array([Mercury.Xpos - Sun.Xpos,
                                Mercury.Ypos - Sun.Ypos,
                                Mercury.Zpos - Sun.Zpos]).flatten()
            velVector = np.array([Mercury.Xvel, Mercury.Yvel, Mercury.Zvel]).flatten()
            Ang_M_Mercury = np.cross(rVector, velVector)
            L_Squared = np.dot(Ang_M_Mercury, Ang_M_Mercury)
            Radial_G = -(GC * Mercury.mass * Sun.mass / DistanceArray[i]) * \
                       (1 + (3 * L_Squared) / (DistanceArray[i] * c * c))
            Xgrav = Radial_G * (Newobject.Xpos - Body.Xpos) / np.sqrt(DistanceArray[i])
            Ygrav = Radial_G * (Newobject.Ypos - Body.Ypos) / np.sqrt(DistanceArray[i])
            Zgrav = Radial_G * (Newobject.Zpos - Body.Zpos) / np.sqrt(DistanceArray[i])
            GForce[0] += Xgrav
            GForce[1] += Ygrav
            GForce[2] += Zgrav
            i += 1
            continue

        
        Xgrav = -((GC * Body.mass * Newobject.mass) / DistanceArray[i]) * \
                ((Newobject.Xpos - Body.Xpos) / np.sqrt(DistanceArray[i]))
        Ygrav = -((GC * Body.mass * Newobject.mass) / DistanceArray[i]) * \
                ((Newobject.Ypos - Body.Ypos) / np.sqrt(DistanceArray[i]))
        Zgrav = -((GC * Body.mass * Newobject.mass) / DistanceArray[i]) * \
                ((Newobject.Zpos - Body.Zpos) / np.sqrt(DistanceArray[i]))
        GForce[0] += Xgrav
        GForce[1] += Ygrav
        GForce[2] += Zgrav
        i += 1

    return GForce

# Find Motion of Newobject in a Timestep
def Motion(Newobject):
    Gforce = Gravity(Newobject).flatten()
    Xacc = Gforce[0] / Newobject.mass
    Yacc = Gforce[1] / Newobject.mass
    Zacc = Gforce[2] / Newobject.mass

    Newobject.Xvel += Xacc * TIMEINC
    Newobject.Yvel += Yacc * TIMEINC
    Newobject.Zvel += Zacc * TIMEINC

    Newobject.Xpos += Newobject.Xvel * TIMEINC
    Newobject.Ypos += Newobject.Yvel * TIMEINC
    Newobject.Zpos += Newobject.Zvel * TIMEINC

    return [Newobject.Xpos, Newobject.Ypos, Newobject.Zpos,
            Newobject.Xvel, Newobject.Yvel, Newobject.Zvel]

# Setup for Simulation
Duration = 20000000
Timesteps = Duration // TIMEINC
XPaths = np.zeros((Num_Bodies, Timesteps), dtype=np.float64)
YPaths = np.zeros((Num_Bodies, Timesteps), dtype=np.float64)
ZPaths = np.zeros((Num_Bodies, Timesteps), dtype=np.float64)
Mercury_Array = []

# Find X, Y, and Z paths for each body
for j in range(Timesteps):
    Xpoints = np.zeros(Num_Bodies, dtype=np.float64)
    Ypoints = np.zeros(Num_Bodies, dtype=np.float64)
    Zpoints = np.zeros(Num_Bodies, dtype=np.float64)

    for i, ClassObject in enumerate(Bodies):
        if ClassObject == Sun:
            Sun.Xpos, Sun.Ypos, Sun.Zpos = 0.0, 0.0, 0.0
            continue

        new_pos = Motion(ClassObject)
        ClassObject.Xpos, ClassObject.Ypos, ClassObject.Zpos = new_pos[0:3]
        ClassObject.Xvel, ClassObject.Yvel, ClassObject.Zvel = new_pos[3:6]

        Xpoints[i] = new_pos[0]
        Ypoints[i] = new_pos[1]
        Zpoints[i] = new_pos[2]

        if ClassObject == Mercury:
            Mercury_Array.append([
                ClassObject.Xpos - Sun.Xpos,
                ClassObject.Ypos - Sun.Ypos,
                ClassObject.Zpos - Sun.Zpos
            ])

    XPaths[:, j] = Xpoints
    YPaths[:, j] = Ypoints
    ZPaths[:, j] = Zpoints

# Create Empty lists for Precession calculations
Dist_Array = []
Min_Indices = []
AngleExtrema_list=[]
AnnualChange_list=[]
Year_Count=[]
Mercury_Array = np.array(Mercury_Array, dtype=np.float64)

# Calculate angular change of perihelion
for point in Mercury_Array:
    x, y, z = point
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    Dist_Array.append(r)

for i in range(len(Dist_Array) - 2):
    if Dist_Array[i]>Dist_Array[i + 1]<Dist_Array[i+2]:
    	Min_Indices.append(i)

raw_angles = np.array([
    np.arctan2(Mercury_Array[i][1], Mercury_Array[i][0])
    for i in Min_Indices
])
unwrapped_angles = np.unwrap(raw_angles) * 206264.8  # Convert to arcseconds

AngleExtrema_list = unwrapped_angles.tolist()

for Angle in range(len(AngleExtrema_list)-1):
    AnnualChange=AngleExtrema_list[Angle+1]-AngleExtrema_list[Angle]
    AnnualChange_list.append(AnnualChange)
    Year_Count.append(Angle)

# print and plot angular shift of perihelion
print(AnnualChange_list)
print(np.mean(AnnualChange_list))
plt.plot(Year_Count,AnnualChange_list,'o')
plt.xlabel('Orbits')
plt.ylabel('Arcseconds of Change in Perihelion')
#plt.ylim(-20,20)
plt.suptitle('Precession of Mercury')
plt.title('True Value is ~6 arcseconds')
plt.grid(True)
plt.show()


# Create Figure for Solar System animation
fig = plt.figure()
fig.patch.set_facecolor('black')
ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor('black') 
ax.set_xlim(-3E12, 3E12)
ax.set_ylim(-3E12, 3E12)
ax.set_zlim(-3E12, 3E12)
ax.set_axis_off() 
ax.grid(False)


Paths = [ax.plot([], [], [], color=colors[i % len(colors)])[0] for i in range(Num_Bodies)]
Points = [
    ax.scatter([], [], [], 'o', color=colors[i % len(colors)], s = 250 if i == 0 else 60)
    for i in range(Num_Bodies)
]

# Place Bodies
def Place_Planets():
    for Path, Point in zip(Paths, Points):
        empty_array = np.array([])
        Path.set_data(empty_array, empty_array)
        Path.set_3d_properties(empty_array)
        Point._offsets3d = ([], [], [])
    return Paths + Points

# Update Bodies
def Next_Frame(frame):
    ax.view_init(elev=70, azim=45)
    for i in range(Num_Bodies):
        Paths[i].set_data(XPaths[i, :frame], YPaths[i, :frame])
        Paths[i].set_3d_properties(ZPaths[i, :frame])
        Points[i]._offsets3d = ([XPaths[i, frame]], [YPaths[i, frame]], [ZPaths[i, frame]])
    return Paths + Points

# Animate System
ani = animation.FuncAnimation(fig, Next_Frame, frames=range(0, Timesteps, 20),
                              init_func=Place_Planets, blit=False, interval=1)

plt.show()