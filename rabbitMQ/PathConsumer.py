import pika
from Particle import *
import json
import sys

particle = None
connection = None
channel = None
particleId = 0

def obradiPutanju(ch, method, properties, pointsJson):

    global particle
    global connection
    global channel
    global particleId
    message = ""

    pointsData = json.loads(pointsJson)

    xList = pointsData['x']
    yList = pointsData['y']
    zList = pointsData['z']
    impactDistanceList = pointsData['impactDistance']

    for i in range(len(xList)):
        x = particle.x + xList[i]
        y = particle.y + yList[i]
        z = particle.z + zList[i]

        v = speed_maxwell(0.2, 293.0, particle)
        v_magnitude = sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

        pos_status = in_position(x, y, z)

        if pos_status == AIR:
            particle.life_time += impactDistanceList[i] / v_magnitude

            if particle.life_time >= particle.half_life:
                message = "AIR"
            else:
                particle.set(x, y, z)
                continue
        # Cestica se nije raspala
        else:
            particle.life_time += impactDistanceList[i] / v_magnitude

            if particle.life_time > particle.half_life:
                message = "AIR"

            else:
                if pos_status == BOTTOM:
                    message = "BOTTOM"
                elif pos_status == TOP:
                    message = "TOP"
                else:
                    message = "WALL"

        ch.basic_publish(exchange='',
                              routing_key='particleCounter',
                              body=message,
                              properties=pika.BasicProperties(
                                  delivery_mode=2,  # make message persistent
                              ))

        ch.basic_ack(delivery_tag=method.delivery_tag)
        ch.basic_cancel(consumer_tag='radonSim')
        print 'Cestica ' + str(particleId) + ' je zavrsila.'
        sys.exit(0)

if __name__ == "__main__":

    try:

        with open('configuration.json') as data_file:
            conf = json.load(data_file)

        particleId = int(sys.argv[1])

        host = conf['host']
        port = int(conf['port'])
        particles_num = int(conf['brojCestica'])

        particle = Particle(particleId)

        connection = pika.BlockingConnection(pika.ConnectionParameters(host=host, port=port))
        channel = connection.channel()

        args = {'x-max-length': particles_num * 2}
        channel.queue_declare(queue='path', arguments=args)
        channel.queue_declare(queue='particleCounter')

        channel.basic_consume(obradiPutanju,
                                queue='path')

        channel.start_consuming()
    except KeyboardInterrupt:
        channel.queue_delete(queue='path')
        channel.queue_delete(queue='particleCounter')
        channel.queue_delete(queue='end')
        connection.close()
