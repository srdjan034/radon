import pika
from Particle import *
import json
import sys

particle = None
connection = None
channel = None
particleId = 0

def uzmiPocetneVrednosti(ch, method, properties, pointsJson):
    global particle

    initialValues = json.loads(pointsJson)
    #print(str(initialValues))
    id = int(initialValues['id'])
    half_life = float(initialValues['half_life'])
    x = float(initialValues['x'])
    y = float(initialValues['y'])
    z = float(initialValues['z'])

    particle = Particle(id, half_life, x, y, z)

    # Odjavi se
    channel.stop_consuming()

# funkcija proverava da li bounding box putanje iskace iz cilindra za prosledjene koordinate cestice
def checkBoundingBox(xParticle, yParticle, zParticle, xMin, yMin, zMin, xMax, yMax, zMax, life_time_step_sum):

    global particle

    x = xParticle + xMax
    y = yParticle + yMin
    z = zParticle + zMin

    if in_position(x, y, z) != AIR:
        return True

    x = xParticle + xMax
    y = yParticle + yMin
    z = zParticle + zMax

    if in_position(x, y, z) != AIR:
        return True

    x = xParticle + xMin
    y = yParticle + yMin
    z = zParticle + zMax

    if in_position(x, y, z) != AIR:
        return True

    x = xParticle + xMin
    y = yParticle + yMin
    z = zParticle + zMin

    if in_position(x, y, z) != AIR:
        return True

    x = xParticle + xMax
    y = yParticle + yMax
    z = zParticle + zMin

    if in_position(x, y, z) != AIR:
        return True

    x = xParticle + xMax
    y = yParticle + yMax
    z = zParticle + zMax

    if in_position(x, y, z) != AIR:
        return True

    x = xParticle + xMin
    y = yParticle + yMax
    z = zParticle + zMax

    if in_position(x, y, z) != AIR:
        return True

    x = xParticle + xMin
    y = yParticle + yMax
    z = zParticle + zMax

    if in_position(x, y, z) != AIR:
        return True

    # proveri da li se cestica raspala
    if (particle.life_time + life_time_step_sum) >= particle.half_life:
        return True

def obradiPutanju(ch, method, properties, pointsJson):

    global particle
    global connection
    global channel
    global particleId
    global ii
    message = ""

    pointsData = json.loads(pointsJson)

    xList = pointsData['x']
    yList = pointsData['y']
    zList = pointsData['z']
    life_time_step = pointsData['life_time_step']

    xMin = pointsData['xMin']
    xMax = pointsData['xMax']
    yMin = pointsData['yMin']
    yMax = pointsData['yMax']
    zMin = pointsData['zMin']
    zMax = pointsData['zMax']

    x_last = pointsData['x_last']
    y_last = pointsData['y_last']
    z_last = pointsData['z_last']

    life_time_step_sum = pointsData['life_time_step_sum']

    # proveri da li "bounding box" putanje iskace iz cilindra za trenutne koordinate cestice
    if checkBoundingBox(particle.x, particle.y, particle.z, xMin, yMin, zMin, xMax, yMax, zMax, life_time_step_sum):
        for i in range(len(xList)):
            x = particle.x + xList[i]
            y = particle.y + yList[i]
            z = particle.z + zList[i]

            pos_status = in_position(x, y, z)

            if pos_status == AIR:
                particle.life_time += life_time_step[i]

                if particle.life_time >= particle.half_life:
                    message = "AIR"
                else:
                    particle.set(x, y, z)
                    continue
            # Cestica se nije raspala
            else:
                particle.life_time += life_time_step[i]

                if particle.life_time > particle.half_life:
                    message = "AIR"

                else:
                    if pos_status == BOTTOM:
                        message = "BOTTOM"
                    elif pos_status == TOP:
                        message = "TOP"
                    else:
                        message = "WALL"
            break
    else:
        particle.life_time += life_time_step_sum
        particle.set(particle.x + x_last, particle.y + y_last, particle.z + z_last)

    if message != "":
        channel.queue_declare(queue='particleCounter')

        ch.basic_publish(exchange='',
                         routing_key='particleCounter',
                         body=message,
                         properties=pika.BasicProperties(
                             delivery_mode=2  # make message persistent
                         ))

        ch.basic_ack(delivery_tag=method.delivery_tag)
        ch.basic_cancel(consumer_tag='radonSim')
        print 'Cestica ' + str(particleId) + ' je zavrsila.'
        print "(" + str(particle.x) + "," + str(particle.y) + "," + str(particle.z) + ")"
        sys.exit(0)



if __name__ == "__main__":

    try:
        with open('configuration.json') as data_file:
            conf = json.load(data_file)

        particleId = int(sys.argv[1])

        host = conf['host']
        port = int(conf['port'])
        particles_num = int(conf['brojCestica'])

        credentials = pika.PlainCredentials(conf['user'], conf['pass'])
        connection = pika.BlockingConnection(pika.ConnectionParameters(host=host, port=port, credentials=credentials))
        channel = connection.channel()

        # Uzmi pocetne vrednosti
        channel.queue_declare(queue='pocetne_vrednosti')
        channel.basic_consume(uzmiPocetneVrednosti,
                              queue='pocetne_vrednosti')

        # hack za stop_consuming()
        while channel._consumer_infos:
            channel.connection.process_data_events(time_limit=1)  # 1 second

        # Uzimaj putanje
        args = {'x-max-length': 40}
        channel.queue_declare(queue='path')#, arguments=args)

        channel.basic_consume(obradiPutanju,
                                queue='path')

        channel.start_consuming()
    except KeyboardInterrupt:
        channel.queue_delete(queue='path')
        channel.queue_delete(queue='particleCounter')
        channel.queue_delete(queue='end')
        connection.close()
