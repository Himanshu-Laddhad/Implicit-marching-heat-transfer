import matplotlib as mt
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from numpy.ma.core import shape, sqrt 


# input values
l = float(10)
rx = int(20)
rb = {
    't' : float(0),
    'name' : "right"
}
lb = {
    't' : float(1),
    'name' : "left"
}
alpha = 1
dt = float(0.001)

#calculation variables
cx = rx
dx = l / rx
print(dx)
a = (dt*alpha)/(dx*dx)
xcells = np.array([(2*i+1)*dx/2 for i in range(cx)])

#weights
ae = -a 
aw = -a
ap = 1 + 2*a

#initialisation
xcells = np.array([(2*i+1)*dx/2 for i in range(rx)])
ot = np.zeros(shape = (rx,1))
nt = np.zeros(shape = (rx,1))
coeff = np.zeros(shape= (rx,cx))
const = np.zeros(shape= (rx,1))
iteration = 0
txtime = []
#matrices formation
while True:
    ot = nt
    for r in range(rx):
        if r == 0: 
            #left
            coeff[r,r] = ap - aw
            coeff[r,r + 1] = ae
            const[r,0] = ot[r,0] - aw*(2*lb['t'])
        elif r == rx - 1: 
            #right
            coeff[r,r] = ap - ae
            coeff[r,r - 1] = aw
            const[r,0] = ot[r,0] - ae*(2*rb['t'])
        else:
            #mid cells
            coeff[r,r] = ap 
            coeff[r,r - 1] = aw
            coeff[r,r + 1] = ae
            const[r,0] = ot[r,0] 

    inv_coeff = np.linalg.inv(coeff)
    nt = np.dot(inv_coeff,const)
    if np.max(nt - ot) < 0.000001:
        break
    iteration += 1
    txtime.append(nt)

net_time = iteration * dt
for i in range(0, int(net_time),int(net_time/10)):
    z = txtime[int(i/dt)]
    plt.plot(xcells,z[:,0])
plt.show()