from scipy.fftpack import fftn, ifftn
import numpy as np 
import matplotlib.pyplot as plt

# def gaussian(A, sigma):
# 	x = np.exp(-(A*A) / 2.0 / sigma / sigma) / np.sqrt(2.0 * np.pi) / sigma
# 	return (x)

# def sigma_k_gauss(k):
# 	return k * k * np.exp(-(k/k0)**2)

# def generate_field(kfunction, kmin, kmax):


# 	A = numpy.random.normal(A, sigma_k)


N = 56
shape = (N,N,N)
kgrid = np.indices(shape, dtype=float).T + 1 - (N/2)

# units are all in kiloparsec
#Lc = 
# domain size 
L = 100
kmin = 2.0 * np.pi / L
dx = L / N
kmax = 2.0 * np.pi / dx 

kgrid *= kmin

# print (np.linalg.norm(kgrid[0,0,0]))
kgrid2 = kgrid * kgrid

kmag = np.sqrt(np.sum(kgrid2, axis=3))

print (kgrid.shape, kgrid2.shape, kmag.shape)

n = (5./3.)+2

#k = np.arange(N, dtype=float)
mod_Ak2 = np.sqrt(np.power(kmag / kmin, -n)).flatten()



#k = np.arange(N, dtype=float)
# first define blank array to keep Atilde in 
Atilde = np.zeros(kgrid.shape, dtype=np.complex_)
amplitudes = np.zeros(kgrid.shape)
phase = np.zeros(kgrid.shape)

# calculate random amplitudes
amplitudes[:,:,:,0] = np.array([np.random.rayleigh(scale=x) for x in mod_Ak2]).reshape(shape)
amplitudes[:,:,:,1] = np.array([np.random.rayleigh(scale=x) for x in mod_Ak2]).reshape(shape)
amplitudes[:,:,:,2] = np.array([np.random.rayleigh(scale=x) for x in mod_Ak2]).reshape(shape)
phase[:,:,:,0] = np.random.random(shape) * 2.0 * np.pi 
phase[:,:,:,1] = np.random.random(shape) * 2.0 * np.pi 
phase[:,:,:,2] = np.random.random(shape) * 2.0 * np.pi 

for i in range(3):
	Atilde[:,:,:,i].real = amplitudes[:,:,:,i] * np.cos(phase[:,:,:,i])
	Atilde[:,:,:,i].imag = amplitudes[:,:,:,i] * np.sin(phase[:,:,:,i])
	#Atilde[:,:,:,i] = np.complex



#Atilde[:,:,:,0].real = np.array([np.random.normal(loc=0.0, scale=x) for x in mod_Ak2]).reshape(shape)
#Atilde[:,:,:,1].real = np.array([np.random.normal(loc=0.0, scale=x) for x in mod_Ak2]).reshape(shape)
#Atilde[:,:,:,2].real = np.array([np.random.normal(loc=0.0, scale=x) for x in mod_Ak2]).reshape(shape)

# calculate random phases
# Atilde[:,:,:,0].imag = np.random.random(shape) * 2.0 * np.pi 
# Atilde[:,:,:,1].imag = np.random.random(shape) * 2.0 * np.pi 
# Atilde[:,:,:,2].imag = np.random.random(shape) * 2.0 * np.pi 




i = np.complex(0,1)
Btilde = np.cross(i * kgrid, Atilde, axisa=3, axisb=3)

# plt.plot(kgrid, Btilde[:,:,:,0])
# plt.show()

# B = ifftn(Btilde).real

# # #plt.imshow(B[:,:,0,0], cmap="Spectral")
# plt.plot(B[:,0,0,0])
# plt.plot(B[:,0,0,1])
# plt.plot(B[:,0,0,2])
# # plt.show()

Brms = np.sqrt(B**2) 
print (np.mean(Brms[:,:,:,0]), np.mean(Brms[:,:,:,1]), np.mean(Brms[:,:,:,2]))

# BB = np.linalg.norm(B, axis=3)

# dBx = np.gradient(B[:,:,:,0], axis=0)
# dBy = np.gradient(B[:,:,:,1], axis=1)
# dBz = np.gradient(B[:,:,:,2], axis=2)
# divB = dBx + dBy + dBz
# plt.imshow(divB[:,:,0]/BB[:,:,0])
# plt.colorbar()
# np.save("Bfield.npy", B)



#print (np.sum(kgrid2, axis=3))
#kmag = np.sqrt(kgrid.dot(kgrid))
#
# from FyeldGenerator import generate_field
# import matplotlib.pyplot as plt
# import numpy as np
# # Helper that generates power-law power spectrum
# def Pkgen(n):
#     def Pk(k):
#         return np.power(k, -n)
#     return Pk
# # Draw samples from a normal distribution
# def distrib(shape):
#     a = np.random.normal(loc=0, scale=1, size=shape)
#     b = np.random.normal(loc=0, scale=1, size=shape)
#     return a + 1j * b
# N = 256
# shape = (N,N,N)
# field = generate_field(distrib, Pkgen(2), shape)
# plt.imshow(field[0], cmap='seismic')
# plt.show()