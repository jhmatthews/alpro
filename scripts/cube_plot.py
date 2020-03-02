import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits import mplot3d

B = np.load("Bfield.npy")

B = np.sqrt(np.sum(B**2, axis=3))

x = np.arange(256)
fig = plt.figure()
ax = plt.axes(projection='3d')

limit = 256
delta = 0.5

cmap = "inferno"

X,Y = np.meshgrid(x,x)
#ax.contourf(X, Y, B[:,:,0,0], 10, zdir='z', offset=256,cmap='afmhot')
ax.contourf(X, Y, B[:,0,:], 10, zdir='z', offset=limit-delta,cmap=cmap)
ax.contourf(X, B[:,0,:], Y, 10, zdir='y', offset=0+delta,  cmap=cmap)
ax.contourf(B[:,0,:], X, Y, 10, zdir='x', offset=limit-delta,cmap=cmap)
#ax.contourf(X, Y, B[0,:,:,0], 10, zdir='x', offset=0,cmap='afmhot')
#ax.contourf(X, Y, B[:,:,0,0], zdir='y', offset=0,cmap='afmhot')
#ax.imshow(B[:,:,0,0])
plt.axis("off")
ax.set_xlim3d(0,limit)
ax.set_ylim3d(0,limit)
ax.set_zlim3d(0,limit)
plt.show()
# ax.pcolormesh()