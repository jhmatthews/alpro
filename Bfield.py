from scipy.fftpack import fftn, ifftn
import numpy as np 
import matplotlib.pyplot as plt

class MagneticField:
	'''
	A container for a magnetic field model
	'''
	def __init__(self, N):
		self.N = N 
		self.shape = (N, N, N)
		self.vec_shape = (N, N, N, 3)
		
	def get_Ak2(self, L):
		self.x = np.linspace(-L,L,self.N)
		self.y = np.linspace(-L,L,self.N)
		self.z = np.linspace(-L,L,self.N)
		kx = 1./self.x
		ky = 1./self.y
		kz = 1./self.z
		kt = np.zeros(3)
		zeta = (11./3.)+2
		kmax = 0.1
		#kmin = L / self.N


		kkx, kky, kkz = np.array(np.meshgrid(kx,ky,kz))
		self.kgrid2 = np.indices(self.shape, dtype=float).T + 1
		kmag = np.sqrt(kkx*kkx + kky*kky + kkz*kkz)
		A_k2 = (kmag)**(-zeta)
		print (np.min(kmag), np.max(kmag))
		A_k2[(kmag < kmin)] = 0.0
		A_k2[(kmag > kmax)] = 0.0
		
		#A_k2 = np.concatenate((A_k2,A_k2,A_k2),axis)

		#print (kkx.shape, kky.shape, kkz.shape)
		print (A_k2)
		self.kgrid = np.zeros(self.vec_shape)
		self.kgrid[:,:,:,0] = kkx
		self.kgrid[:,:,:,1] = kky
		self.kgrid[:,:,:,2] = kkz
		#self.kgrid = np.concatenate((kkx,kky,kkz))
		#print (self.kgrid.shape)

		return (A_k2)

	def get_random_field_tribble(self, spectrum, mod_Ak2):

		# first define blank array to keep Atilde in 
		self.Atilde = np.zeros(self.kgrid.shape, dtype=np.complex_)
		amplitudes = np.zeros(self.kgrid.shape)
		#phase = np.zeros(mod_Ak2.shape)

		#Akx = mod_Ak2[:,:,:,0].flatten()
		#Aky = mod_Ak2[:,:,:,1].flatten()
		#Akz = mod_Ak2[:,:,:,2].flatten()

		# calculate random amplitudes according to a rayleigh distribution 
		amplitudes[:,:,:,0] = (1/(2*np.pi))*np.random.rayleigh(scale=mod_Ak2)
		amplitudes[:,:,:,1] = (1/(2*np.pi))*np.random.rayleigh(scale=mod_Ak2)
		amplitudes[:,:,:,2] = (1/(2*np.pi))*np.random.rayleigh(scale=mod_Ak2)
		#amplitudes[:,:,:,1] = np.array([np.random.rayleigh(scale=x) for x in Aky]).reshape(self.shape)
		#amplitudes[:,:,:,2] = np.array([np.random.rayleigh(scale=x) for x in Akz]).reshape(self.shape)

		# calculate random phases (uniform in 2pi)
		phase = np.random.uniform(size=self.kgrid.shape) * 2.0 * np.pi 
		#phase[:,:,:,0] = np.random.random(self.shape) * 2.0 * np.pi 
		#phase[:,:,:,1] = np.random.random(self.shape) * 2.0 * np.pi 
		#phase[:,:,:,2] = np.random.random(self.shape) * 2.0 * np.pi 

		# turn polar coordinates into rectangular complex form 
		for i in range(3):
			#self.Atilde[:,:,:,i].real = amplitudes[:,:,:,i] * np.cos(phase[:,:,:,i])
			#self.Atilde[:,:,:,i].imag = amplitudes[:,:,:,i] * np.sin(phase[:,:,:,i])
			self.Atilde[:,:,:,i] = amplitudes[:,:,:,i]*np.exp(1j*phase[:,:,:,i])

		i_ = np.complex(0,1)
		#print (self.kgrid.shape, self.Atilde.shape)
		print (self.kgrid)
		self.Btilde = np.cross(i_ * self.kgrid, self.Atilde, axisa=3, axisb=3)
		# self.Btilde *= i
		Ax = self.Atilde[:,:,:,0]
		Ay = self.Atilde[:,:,:,1]
		Az = self.Atilde[:,:,:,2]
		kx = self.kgrid[:,:,:,0]
		ky = self.kgrid[:,:,:,1]
		kz = self.kgrid[:,:,:,2]
		B_x = (Az*1j*ky) - (Ay*1j*kz)
		B_y = -((Az*1j*kx) - (Ax*1j*kz))
		B_z = (Ay*1j*kx) - (Ax*1j*ky)
		B = np.concatenate((B_x,B_y,B_z))

		# finally take Inverse fourier transform using scipy fftpack and get real part 
		B = ifftn(self.Btilde).real
		#print (self.Atilde)
		#print (amplitudes, phase)

		return (B)



N = 256
# units are all in kiloparsec
L = 100  		# size of cluster 
dx = L / N 		# resolution of cluster 
Bfield = MagneticField(N)

# create x,y,z arrays 
# x = np.arange(-L,L,dx) #range in x
# z = np.arange(-L,L,dx) #range in y
# y = np.arange(-L,L,dx) #range in z

kmin = 2.0 * np.pi / L
dx = L / Bfield.N
kmax = 2.0 * np.pi / dx 

# Bfield.kgrid *= kmin

# print (np.linalg.norm(kgrid[0,0,0]))
#kgrid2 = Bfield.kgrid * Bfield.kgrid
#kmag = np.sqrt(np.sum(kgrid2, axis=3))
#
#print (Bfield.kgrid.shape, kgrid2.shape, kmag.shape)

# n = (11./3.)+2

# #k = np.arange(N, dtype=float)
# mod_Ak2 = np.sqrt(np.power(Bfield.kgrid / kmin, -n))
# mod_Ak2[mod_Ak2 == np.inf] = 1e50 
# mod_Ak2[mod_Ak2 <= 0] = 1e-16
# print (mod_Ak2.shape)

mod_Ak2 = Bfield.get_Ak2(L)

print (Bfield.kgrid)

B = Bfield.get_random_field_tribble(None, mod_Ak2)



#Atilde[:,:,:,0].real = np.array([np.random.normal(loc=0.0, scale=x) for x in mod_Ak2]).reshape(shape)
#Atilde[:,:,:,1].real = np.array([np.random.normal(loc=0.0, scale=x) for x in mod_Ak2]).reshape(shape)
#Atilde[:,:,:,2].real = np.array([np.random.normal(loc=0.0, scale=x) for x in mod_Ak2]).reshape(shape)

# calculate random phases
# Atilde[:,:,:,0].imag = np.random.random(shape) * 2.0 * np.pi 
# Atilde[:,:,:,1].imag = np.random.random(shape) * 2.0 * np.pi 
# Atilde[:,:,:,2].imag = np.random.random(shape) * 2.0 * np.pi 




#plt.plot(Bfield.kgrid[:,:,0,0], Bfield.Btilde[:,:,0,0])
#plt.show()





Brms = np.sqrt(B**2) 
print (np.mean(Brms[:,:,:,0]), np.mean(Brms[:,:,:,1]), np.mean(Brms[:,:,:,2]))

BB = np.linalg.norm(B, axis=3)

dBx = np.gradient(B[:,:,:,0], axis=0)
dBy = np.gradient(B[:,:,:,1], axis=1)
dBz = np.gradient(B[:,:,:,2], axis=2)
divB = dBx + dBy + dBz

#plt.imshow(BB[:,:,0])
#plt.colorbar()
#plt.show()
plt.subplot(221)
titles = ["Bx", "By", "Bz"]

for i in range(3):
	plt.subplot(3,2,i+1)
	plt.imshow(B[:,:,0,i])

plt.subplot(3,2,4)
plt.title("B")
plt.imshow(BB[:,:,0])

plt.subplot(3,2,5)
plt.title("div(B)")
plt.imshow(divB[:,:,0])
	#plt.clf


#np.save("Bfield.npy", B)
# plt.clf()
# plt.figure()
plt.subplot(3,2,6)
plt.plot(B[0,:,0,0])
plt.plot(B[0,:,0,1])
plt.plot(B[0,:,0,2])
plt.show()

print (np.mean(divB), np.mean(BB))
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