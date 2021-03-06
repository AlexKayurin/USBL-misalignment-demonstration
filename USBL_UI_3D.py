"""
akayurin@gmail.com
"""
import math
import random # apply to slantrange and rel_az
import numpy as np
from win32api import GetSystemMetrics
# pyplot backend
import matplotlib
matplotlib.use('TkAgg')
#matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Button, RadioButtons, Slider
#print(matplotlib.pyplot.get_backend())

# Get screen resolution and sides ratio
sc_rat = 1.25 / (GetSystemMetrics(0) / GetSystemMetrics(1)) # set original ratio

# Conditions init
offset = 0 # Line offset
dz = 50 # Depth
hdg = 0 # Heading
noise = 1 # Noise
# Corrections init
dS = 1 # Scale factor
dH = 0 # Yaw
dP = 0 # Pitch
dR = 0 # Roll

def ves_shape(hdg): # arrays for vessel shape, water, ground, north arraow  
    Shape_top = np.array([[10.198 * math.sin(hdg + math.pi + 0.197396), 10.198 * math.cos(hdg + math.pi + 0.197396)],
                    [8.246 * math.sin(hdg - 0.244979), 8.246 * math.cos(hdg - 0.244979)],
                    [14 * math.sin(hdg), 14 * math.cos(hdg)],
                    [8.246 * math.sin(hdg + 0.244979), 8.246 * math.cos(hdg + 0.244979)],
                    [10.198 * math.sin(hdg + math.pi - 0.197396), 10.198 * math.cos(hdg + math.pi - 0.197396)],
                    [10.198 * math.sin(hdg + math.pi + 0.197396), 10.198 * math.cos(hdg + math.pi + 0.197396)]])        
    Shape_stern = np.array([[-4, -3], [-4, 2], [-3, 2], [-3, 4], [-2, 4], [-2, 5.4],
                           [-1.2, 5.4], [-1.2, 6.5], [-0.9, 6.5], [-0.9, 5.4], [-0.2, 5.4],
                           [-0.1, 8.5], [0.1, 8.5], [0.2, 5.4],
                           [2, 5.4], [2, 4], [3, 4], [3, 2], [4, 2], [4, -3], [3, -4], [-3, -4], [-4, -3]])    
    Shape_port = np.array([[-14, 2.5], [-12, 2.5], [-12, 5.5], [-11.8, 5.5], [-11.8, 2.5],
                           [-8, 2.5], [-8, 4], [-7, 4], [-7, 5.4], [-4.5, 5.4], [-4.5, 8.5],
                           [-4.3, 8.5], [-4.3, 5.4], [-3, 5.4], [-3, 4],
                           [0, 4], [0, 1.5], [10, 1.5], [10, -1], [9, -2.5], [-9, -2.5], [-14, 2]])
    Water = np.array([[-100, -50], [-100, 0], [100, 0], [100, -50]])
    Ground = np.array([[-100, -90], [-100, -50], [100, -50], [100, -90]])
    North_arrow = np.array([[0, 0], [-3, -5], [0, 15], [3, -5]])
    
    return [Shape_top, Shape_stern, Shape_port, Water, Ground, North_arrow] 
   
def calc_azimuth(dx, dy): # azimuth
    if dy == 0  and dx >= 0:
        azimuth = math.pi / 2
    elif dy == 0 and dx < 0:
        azimuth = 3 * math.pi / 2
    else:
        azimuth = math.atan(dx / dy)
    if (dx >= 0 and dy <0) or (dx < 0 and dy <0):
        azimuth = azimuth + math.pi                   
    return azimuth

def calc_thA(a, b): # vertical angles
    if a == 0:
        thA = 0 
    else:
        thA = math.atan(a / b)
    return thA

def update_slider(val): # update heading slider (in global namespace)
     hdg_slider.ch_val(math.degrees(val))
        
#INIT PLOT
fig = plt.figure(num='akayurin@gmail.com - USBL installation errors demonstration - 2021') 
# Plan
ax_plan = plt.axes([-0.02, 0.19, 0.62 * sc_rat, 0.62], aspect=1, frame_on=True, xlim=(-75, 75), ylim=(-75, 75))
plt.title('Plan view', loc='center', fontsize=8)
plt.tick_params(labelsize=8)
# Back
ax_back = plt.axes([0.55, 0.47, 0.4 * sc_rat, 0.4], aspect=1, xlim=(-100, 100), ylim=(-90, 10))
ax_back.yaxis.tick_right()
plt.tick_params(labelsize=8)
plt.ylabel('Stern view', fontsize=8)
# Port
ax_port = plt.axes([0.55, 0.135, 0.4 * sc_rat, 0.4], aspect=1, xlim=(-100, 100), ylim=(-90, 10))
ax_port.yaxis.tick_right()
plt.tick_params(labelsize=8)
plt.ylabel('Portside view', fontsize=8)

mng = plt.get_current_fig_manager()
mng.window.state('zoomed') # TkAgg
#mng.window.showMaximized() #Qt5Agg

# STATIONARY FIGURES ON PLOTS
# Plan
range_circle1 = ax_plan.add_patch(plt.Circle((0, 0), 25, color='g', alpha=0.1, ls=':', lw=0.5)) # Range 50%
range_circle2 = ax_plan.add_patch(plt.Circle((0, 0), 50, color='g', alpha=0.05, ls=':', lw=0.5)) # Range 100%
beacon_location_top = ax_plan.add_patch(plt.Circle((0, 0), 1, fc='b', alpha=0.7)) # beacon location patch
ax_plan.plot([0, 0], [-75, 750], 'g', lw=0.1) # V Axis 
ax_plan.plot([-75, 75], [0, 0], 'g', lw=0.1) # H Axis
ax_plan.text(25, 0, '50%', style='italic', color='g', fontsize=8)
ax_plan.text(50, 0, '100%', style='italic', color='g', fontsize=8)
north_arrow = ax_plan.add_patch(plt.Polygon(ves_shape(hdg)[5] + (65, 55), fc='k', alpha=0.7))
# Back
beacon_location_back = ax_back.add_patch(plt.Circle((0, -50), 1.5, fc='b', alpha=0.7)) # beacon location patch
ax_back.plot([-100, 100], [-50, -50], 'k:', lw=1) # Seabed 
ax_back.plot([-100, 100], [0, 0], 'b', lw=0.5) # Water
ax_back.plot([0, 0], [-90, 10], 'g', lw=0.1) # V Axis
water_stern = ax_back.add_patch(plt.Polygon(ves_shape(hdg)[3], fc='c', alpha=0.1))
ground_stern = ax_back.add_patch(plt.Polygon(ves_shape(hdg)[4], fc='y', alpha=0.5))
# Port
beacon_location_port = ax_port.add_patch(plt.Circle((0, -50), 1.5, fc='b', alpha=0.7)) # beacon location patch
ax_port.plot([-100, 100], [-50, -50], 'k:', lw=1) # Seabed 
ax_port.plot([-100, 100], [0, 0], 'b', lw=0.5) # Water
ax_port.plot([0, 0], [-90, 10], 'g', lw=0.1) # V Axis
water_port = ax_port.add_patch(plt.Polygon(ves_shape(hdg)[3], fc='c', alpha=0.1))
ground_port = ax_port.add_patch(plt.Polygon(ves_shape(hdg)[4], fc='y', alpha=0.5))

# lists to store x and y axis points init
edata_vesseltrack, ndata_vesseltrack = [], [] # vessel tracklines plan (e, n)
edata_beacontrack, ndata_beacontrack, zdata_beacontrack = [], [], [] # beacon tracklines plan (e, n, z)
xdata_beacontrack, ydata_beacontrack = [], [] # beacon tracklines back/port (x, y, z from above)

def anim_init(): # init moving objects
    global vessel
    # Plan
    beacon = ax_plan.add_patch(plt.Circle((0, 0), 0.5, fc='r')) # beacon track patch init
    vesselshape = ax_plan.add_patch(plt.Polygon(ves_shape(hdg)[0], fc='b', alpha=0.6)) # vessel shape patch init
    line_vesseltrack, = ax_plan.plot([], [], ':b', lw=0.5) # vessel track line init
    beacontrack, = ax_plan.plot([], [], 'r', lw=0.5) # beacon track line init            
    # Back
    vesselshape_back = ax_back.add_patch(plt.Polygon(ves_shape(hdg)[1], fc='b', alpha=0.6)) # vessel shape patch init
    beacon_back = ax_back.add_patch(plt.Circle((0, 0), 0.5, fc='r')) # beacon track patch init
    beacontrack_back, = ax_back.plot([], [], 'r', lw=0.5) # beacon track line init 
    # Port
    vesselshape_port = ax_port.add_patch(plt.Polygon(ves_shape(hdg)[2], fc='b', alpha=0.6))
    beacon_port = ax_port.add_patch(plt.Circle((0, 0), 0.5, fc='r')) # beacon track patch init
    beacontrack_port, = ax_port.plot([], [], 'r', lw=0.5) # beacon track line init         
    # Combine all objects into one list for single animation
    vessel=[line_vesseltrack, vesselshape, beacontrack, beacon, vesselshape_back, beacontrack_back, beacon_back, vesselshape_port, beacontrack_port, beacon_port]   
    return vessel   

def animate_beacon(i): # ANIMATION
    global e_vesseltrack
    global n_vesseltrack   
    global hdg      
# VESSEL CALCULATIONS--------------------------------------------------------------------------------------      
    if spin == True or circulation == True:
        hdg = math.radians(i * 3.6)
        if math.degrees(hdg) >= 360:
            hdg = hdg - 2 * math.pi
        update_slider(hdg) # update heading slider
        
        if spin == False and circulation == True:
            e_vesseltrack = - 50 * math.cos(-hdg)
            n_vesseltrack = - 50 * math.sin(-hdg)
        else:
            e_vesseltrack = 0  
            n_vesseltrack = 0 
    else:
        e_vesseltrack = offset  
        n_vesseltrack = i - 50        
    
    edata_vesseltrack.append(e_vesseltrack) # vessel track arrays 
    ndata_vesseltrack.append(n_vesseltrack) 

# BEACON CALCULATIONS--------------------------------------------------------------------------------------    
    de = -e_vesseltrack
    dn = -n_vesseltrack    
    hordist = math.hypot(de, dn) 
    slantrange = math.hypot(hordist, dz) 
    
    az = calc_azimuth(de, dn)
    rel_az = az - hdg
    
    dx = hordist * math.sin(rel_az)
    dy = hordist * math.cos(rel_az)
    
    thX = calc_thA(dx, dz) # X axis to Z axis   
    thY = calc_thA(dy, dz) # Y axis to Z axis   
    thZ = calc_thA(hordist, dz) # slant range to Z axis
    
    rel_az_corr = rel_az + dH + random.gauss(0, math.radians(noise / 2))
    slantrange_corr = dS * slantrange + random.gauss(0, noise / 2)
    # Correct slant range
    hordist_0 = slantrange_corr * math.sin(thZ)
    dx_0 = hordist_0 * math.sin(rel_az)
    dy_0 = hordist_0 * math.cos(rel_az)
    dz_0 = (slantrange_corr ** 2 - hordist_0 ** 2) ** 0.5  
    # Correct heading
    hordist_1 = hordist_0
    dx_1 = hordist_1 * math.sin(rel_az_corr)
    dy_1 = hordist_1 * math.cos(rel_az_corr)
    dz_1 = dz_0      
    # Correct pitch
    dx_2 = dx_1
    dy_2 = dy_1 - 2 * math.hypot(dy_1, dz_1) * math.sin(dP / 2) * math.cos(thY + dP / 2)
    dz_2 = dz_1 - 2 * math.hypot(dy_1, dz_1) * math.sin(dP / 2) * math.cos(thY + dP / 2) * math.tan(thY + dP / 2)   
    # Correct roll
    dx_3 = dx_2 - 2 * math.hypot(dx_2, dz_2) * math.sin(dR / 2) * math.cos(thX + dR / 2)
    dy_3 = dy_2
    dz_3 = dz_2 - 2 * math.hypot(dx_2, dz_2) * math.sin(dR / 2) * math.cos(thX + dR / 2) * math.tan(thX + dR / 2)
    hordist_3 = math.hypot(dx_3, dy_3) 
    rel_az_corr_3 = calc_azimuth(dx_3, dy_3)   
    # dx, dy into de, dn  
    de_corr = hordist_3 * math.sin(hdg + rel_az_corr_3) 
    dn_corr = hordist_3 * math.cos(hdg + rel_az_corr_3)    
    
    e_beacon = e_vesseltrack + de_corr # beacon track coordinates
    n_beacon = n_vesseltrack + dn_corr
    
    edata_beacontrack.append(e_beacon)  # beacon track arrays
    ndata_beacontrack.append(n_beacon)
    zdata_beacontrack.append(-dz_3)
    xdata_beacontrack.append(dx_3 - dx)
    ydata_beacontrack.append(dy - dy_3)
    
    # BEACON PLOT-------------------------------------------
    # Plan
    vessel[2].set_data(edata_beacontrack, ndata_beacontrack) # beacon track array   
    vessel[3].center = (e_beacon, n_beacon) # beacon shape position 
    # Back
    vessel[5].set_data(xdata_beacontrack, zdata_beacontrack) # beacon track array   
    vessel[6].center = (dx_3 - dx, -dz_3) # beacon shape position     
    # Port
    vessel[8].set_data(ydata_beacontrack, zdata_beacontrack) # beacon track array   
    vessel[9].center = (dy - dy_3, -dz_3) # beacon shape position     
    # VESSEL PLOT----------------------------------------- 
    # Plan
    vessel[0].set_data(edata_vesseltrack, ndata_vesseltrack)  # vessel track array  
    vessel[1].set_xy(ves_shape(hdg)[0] + (e_vesseltrack, n_vesseltrack)) # vessel shape coordinates  
    # Back 
    vessel[4].set_xy(ves_shape(hdg)[1] + (-dx, 0)) # vessel shape coordinates
    # Port
    vessel[7].set_xy(ves_shape(hdg)[2] + (dy, 0)) # vessel shape coordinates
    
#INTERFACE--------------------------------------------------------------------------------------    
def Run(val): # Run button
    global gogogo
    global spin
    global circulation
    # Check mode   
    if mode_button.value_selected == 'Spin':
        spin = True
    if mode_button.value_selected == 'Circle':
        spin = False
        circulation = True
    if mode_button.value_selected == 'Line':
        spin = False
        circulation = False

    if 'vessel' in globals(): # remove previous track if 'vessel' created in anim_init      
        vessel[0].remove()
        vessel[1].remove()
        vessel[2].remove()
        vessel[3].remove()
        vessel[4].remove()
        vessel[5].remove()
        vessel[6].remove()
        vessel[7].remove()
        vessel[8].remove()
        vessel[9].remove()        
        del edata_vesseltrack[:]
        del ndata_vesseltrack[:]
        del edata_beacontrack[:]
        del ndata_beacontrack[:]
        del zdata_beacontrack[:]
        del xdata_beacontrack[:]
        del ydata_beacontrack[:]                        

    gogogo = animation.FuncAnimation(fig, animate_beacon, init_func=anim_init, frames=101, interval=1, repeat=False)

def change_slider(val): # Get changed slider value
    global offset
    global hdg
    global noise
    global dS
    global dH
    global dP
    global dR
           
    offset = off_slider.get_val() / 2
    hdg = math.radians(hdg_slider.get_val())
    noise = noi_slider.get_val()
    dS = scl_slider.get_val()
    dH = math.radians(yaw_slider.get_val())
    dP = math.radians(pitch_slider.get_val())
    dR = math.radians(roll_slider.get_val())
    
class Slider():
    def __init__(self, position, text, v_min, v_max, v_init, v_step):
        self.pl = plt
        self.ax = self.pl.axes(position)        
        self.slid = self.pl.Slider(self.ax, text, valmin=v_min, valmax=v_max, valinit=v_init, dragging=True, valstep=v_step)
        self.slid.label.set_size(8)
        self.slid.valtext.set_size(8)
        self.slid.on_changed(change_slider)
    def ch_val(self, val):
        self.slid.set_val(val)
    def get_val(self):
        return self.slid.val

# Run button
ax_run = plt.axes([0.8, 0.03, 0.1 * sc_rat, 0.1])
b_run = Button(ax_run, 'Run')
b_run.on_clicked(Run)
# Mode selection radiobuttons
ax_mode = plt.axes([0.8, 0.75, 0.1 * sc_rat, 0.3], aspect=1)
mode_button = RadioButtons(ax_mode, ['Line', 'Spin', 'Circle'], active=0, activecolor= 'b')
# Slider
off_slider = Slider([0.11, 0.85, 0.6 * sc_rat, 0.025], 'Offset %WD', -100 , 100, 0, 10)
hdg_slider = Slider([0.11, 0.89, 0.6 * sc_rat, 0.025], 'Heading', 0 , 360, 0, 5)
noi_slider = Slider([0.11, 0.93, 0.6 * sc_rat, 0.025], 'Noise', 0, 1, 0.5, 0.1)
scl_slider = Slider([0.11, 0.12, 0.6 * sc_rat, 0.025], 'Scale error', 0.8, 1.2, 1, 0.001)
yaw_slider = Slider([0.11, 0.09, 0.6 * sc_rat, 0.025], 'Yaw error', -10, 10, 0, 0.1)
pitch_slider = Slider([0.11, 0.06, 0.6 * sc_rat, 0.025], 'Pitch error', -10, 10, 0, 0.1)
roll_slider = Slider([0.11, 0.03, 0.6 * sc_rat, 0.025], 'Roll error', -10, 10, 0, 0.1)

plt.show()
