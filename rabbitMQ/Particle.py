from random import random, seed
from math import *

# Konstante
INIT=-1
AIR=0
WALL=1
TOP=2
BOTTOM=3

# Parametri posla
r = 0.04
h = 0.08
half_path = 6.0e-6
tau = 183/log(2)


# Stampaj status
def print_status(nAir, nWall, nTop, nBottom):
	print ("Bottom=%d Top=%d Wall=%d Air=%d" % (nBottom, nTop, nWall, nAir))


# Distribucija vremena zivota cestice
def distrib(x):
	return (-x*log(random()))

# Gausova distribucija
def Gauss(particle):
	if particle.gauss_kon == 0:
		s1 = random()
		s2 = random()
		temp = sqrt(-2 * log(s1))
		fi = 2 * pi * s2
		value = temp * cos(fi)
		particle.gauss_y = temp * sin(fi)
		particle.gauss_kon = 1
	else:
		value = particle.gauss_y
		particle.gauss_kon = 0

	return value


# Brzina po Maxwell raspodeli
def speed_maxwell(am, te, particle):
	coeff = sqrt(8.314 * te / am)
	return [Gauss(particle)*coeff, Gauss(particle)*coeff, Gauss(particle)*coeff]

# Trenutna pozicija cestice
def in_position(x,y,z):
	result = AIR
	xy = sqrt(x**2 + y**2)

	if xy>=r: # u zidu
		result = WALL
	elif z>=h: # in top
		result = TOP
	elif z<=0: # in bottom
		result = BOTTOM

	return result


def distance(p1, p2):
	return sqrt((p1.x-p2.x)**2+(p1.y-p2.y)**2+(p1.z-p2.z)**2)

#
# Klasa za particle
#
class Particle:
	def __init__(self, id):
		self.id = id
		self.life_time = 0.
		self.half_life = distrib(tau)
		self.status = INIT
		self.gauss_kon = 0
		self.gauss_y = 0

		self.z = 1.e-2 * random()
		rr = r*sqrt(random())
		fir = 2*pi*random()
		self.x = rr*cos(fir)
		self.y = rr*sin(fir)

	def __str__(self):
		return str(self.status)

	def set(self,x,y,z):
		self.x,self.y,self.z = x,y,z