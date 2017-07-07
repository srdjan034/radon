import time
import sys
import pika
from Particle import *
from threading import Thread
import json
from pprint import pprint

channel = None
connection = None
running = True

def oznaciKraj(ch, method, properties, poruka):
    global running
    running = False
    sys.exit(0)

def CekajNaPorukuZaKraj():
    global channel
    channel.queue_declare(queue='end')
    channel.basic_consume(oznaciKraj,
                           queue='end',
                           no_ack=True)
    channel.start_consuming()

if __name__ == "__main__":

    with open('configuration.json') as data_file:
        conf = json.load(data_file)

    host = conf['host']
    port = int(conf['port'])
    # broj tacaka u jednoj kreiranoj putanji
    pathPoints_num = int(conf['brojTacakaUJednojPutanji'])
    # broj cestica
    particles_num = int(conf['brojCestica'])

    connection = pika.BlockingConnection(pika.ConnectionParameters(host=host, port=port))
    channel = connection.channel()

    # Maksimalan broj poruka u redu
    args = {'x-max-length' : particles_num * 2}
    q = channel.queue_declare(queue='path', arguments=args)

    # u posebnoj niti cekaj na poruku od RadonCounter.py da je simulacija zavrsena
    thread = Thread(target = CekajNaPorukuZaKraj)
    thread.start()

    try:
        while running:

            xList = []
            yList = []
            zList = []
            impactDistance = []

            # kreiraj putanju od pathPoints_num tacaka
            for i in range(pathPoints_num):
                impact_distance = distrib(half_path)
                fi0 = 2 * pi * random()
                theta0 = acos(1 - 2 * random())

                xList.append(impact_distance * sin(theta0) * cos(fi0))
                yList.append(impact_distance * sin(theta0) * sin(fi0))
                zList.append(impact_distance * cos(theta0))
                impactDistance .append(impact_distance)

            points = '{ "x" : ' + str(xList) + ', "y" : ' + str(yList) + ', "z" : ' + str(zList)  \
                     + ', "impactDistance" : ' + str(impactDistance)  + ' }'

            # posalji putanju u json formatu
            channel.basic_publish(exchange='',
                                routing_key='path',
                                body=points,
                                properties=pika.BasicProperties(
                                    delivery_mode=1
                                ))

        connection.close()
    except KeyboardInterrupt:
        channel.queue_delete(queue='path')
        channel.queue_delete(queue='particleCounter')
        channel.queue_delete(queue='end')
        connection.close()