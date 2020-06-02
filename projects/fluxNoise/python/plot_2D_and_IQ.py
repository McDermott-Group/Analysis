import noiselib
reload(noiselib)
import numpy as np
import matplotlib.pyplot as plt
import tkFileDialog
import os

def _circle(xy, r):
    """Returns points for a circle."""
    n = 100;
    x = xy[0] + r*np.cos(np.linspace(0,2*np.pi,n))
    y = xy[1] + r*np.sin(np.linspace(0,2*np.pi,n))
    return x,y

if 'last_file_opened' in locals():
    last_path = last_file_opened
else:
    last_path = os.getenv('DATA_ROOT').replace('\\','/')
path = tkFileDialog.askopenfilename(initialdir=last_path)
# path = r'Z:\mcdermott-group\data\fluxNoise2\DR1 - 2019-12-17\CorrFar\Q2\General\05-10-20\Qubit_spectroscopy\MATLABData\Qubit_spectroscopy_002.mat'
last_file_opened = path
data = noiselib.loadmat(path)
o = np.array(data['Single_Shot_Occupation'])
i = np.array(data['Is'])
q = np.array(data['Qs'])

fig_o, ax_o = plt.subplots()
ax_o.imshow(o.transpose())
ax_o.set_aspect(1.*o.shape[0]/o.shape[1])
fig_o.suptitle(path)
# fig_o.colorbar(o)
plt.draw()
plt.pause(0.05)

fig_iq, ax_iq = plt.subplots()
iq_plot = ax_iq.scatter([], [], s=0.5)
state0, = ax_iq.plot([],[], 'k', linewidth=2)
state1, = ax_iq.plot([],[], 'k', linewidth=2)
state2, = ax_iq.plot([],[], 'k', linewidth=2)
ax_iq.set_aspect(1.)
plt.draw()
plt.pause(0.05)

coords = []

def onclick(event):
    x,y = event.xdata, event.ydata
    x,y = int(np.round(x)),int(np.round(y))
    # print x,y
    iq_plot.set_offsets( np.dstack((i[x,y],q[x,y]))[0] )
    cx0,cy0 = _circle( (np.array(data['Center_of_0_State_I'])[x,y],
                        np.array(data['Center_of_0_State_Q'])[x,y]),
                        np.array(data['Std_of_0_State'])[x,y] )
    cx1,cy1 = _circle( (np.array(data['Center_of_1_State_I'])[x,y],
                        np.array(data['Center_of_1_State_Q'])[x,y]),
                        np.array(data['Std_of_1_State'])[x,y] )
    cx2,cy2 = _circle( (np.array(data['Center_of_2_State_I'])[x,y],
                        np.array(data['Center_of_2_State_Q'])[x,y]),
                        np.array(data['Std_of_2_State'])[x,y] )
    state0.set_xdata(cx0)
    state0.set_ydata(cy0)
    state1.set_xdata(cx1)
    state1.set_ydata(cy1)
    state2.set_xdata(cx2)
    state2.set_ydata(cy2)
    plt.draw()
    plt.pause(0.05)
cid = fig_o.canvas.mpl_connect('button_press_event', onclick)

