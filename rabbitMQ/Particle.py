from math import *
from random import random, seed, getstate, setstate


# Konstante
INIT=-1
AIR=0
WALL=1
TOP=2
BOTTOM=3

gauss_kon = 0
gauss_y = 0.0

# Parametri posla
r = 0.04
h = 0.08
half_path = 6.0e-6
tau = 183/log(2)

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

# Stampaj status
def print_status(nAir, nWall, nTop, nBottom):
	print ("Bottom=%d Top=%d Wall=%d Air=%d" % (nBottom, nTop, nWall, nAir))

# Distribucija vremena zivota cestice
def distrib(x):
	return (-x*log(random()))

# Gausova distribucija
def Gauss():
	global gauss_kon
	global gauss_y
	if gauss_kon == 0:
		s1 = random()
		s2 = random()
		temp = sqrt(-2 * log(s1))
		fi = 2 * pi * s2
		value = temp * cos(fi)
		gauss_y = temp * sin(fi)
		gauss_kon = 1
	else:
		value = gauss_y
		gauss_kon = 0

	return value

# Brzina po Maxwell raspodeli
def speed_maxwell(am, te):
	coeff = sqrt(8.314 * te / am)
	return [Gauss()*coeff, Gauss()*coeff, Gauss()*coeff]

def distance(p1, p2):
	return sqrt((p1.x-p2.x)**2+(p1.y-p2.y)**2+(p1.z-p2.z)**2)

#
# Klasa za particle
#
class Particle:
	def __init__(self, id, half_life, x, y, z):
		self.id = id
		self.life_time = 0.
		self.half_life = half_life
		self.status = INIT
		self.x = x
		self.y = y
		self.z = z

	def __str__(self):
		return str(self.status)

	def set(self,x,y,z):
		self.x,self.y,self.z = x,y,z
