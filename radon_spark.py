#!/usr/bin/env python

import sys
from random import random, seed
from math import *
from pyspark.sql import SparkSession
from pyspark.context import SparkContext
from Particle import *


def obradiCesticu(p):

    while True:

        impact_distance = distrib(half_path)
        fi0 = 2 * pi * random()
        theta0 = acos(1 - 2 * random())

        # Probne nove koordinate
        x = p.x + impact_distance * sin(theta0) * cos(fi0)
        y = p.y + impact_distance * sin(theta0) * sin(fi0)
        z = p.z + impact_distance * cos(theta0)
        v = speed_maxwell(0.2, 293.0, p)
        v_magnitude = sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

        pos_status = in_position(x, y, z)

        if pos_status == AIR:

            p.life_time += impact_distance / v_magnitude

            if p.life_time >= p.half_life:
                return ("AIR", 1)
            else:
                p.set(x, y, z)
                continue
        # Cestica se nije raspala
        else:
            p.life_time += impact_distance / v_magnitude

            if p.life_time > p.half_life:
                return ("AIR", 1)

            else:
                if pos_status == BOTTOM:
                    return ("BOTTOM", 1)
                elif pos_status == TOP:
                    return ("TOP", 1)
                else:
                    return ("WALL", 1)


def main():

    particles_num = int(sys.argv[1])

    spark = SparkSession.builder.appName("Radon").getOrCreate()
    spark.sparkContext.addPyFile("Particle.py")

    # Generisanje cestica
    particles = []
    particles = [Particle(i) for i in range(particles_num)]

    # Kreiraj RDD od liste cestica
    rddParticles = spark.sparkContext.parallelize(particles)

    # Svaki maper obradjuje jednu cesticu
    rddParticlesReduced = rddParticles.map(lambda x: obradiCesticu(x))

    rezultat = rddParticlesReduced.reduceByKey(lambda a, b : a + b)

    print rezultat.take(4)

    spark.stop()

# Pozovi main() rutinu
if __name__=='__main__':
	main()