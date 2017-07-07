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

    while channel._consumer_infos:
        channel.connection.process_data_events(time_limit=1)  # 1 second

if __name__ == "__main__":

    with open('configuration.json') as data_file:
        conf = json.load(data_file)

    host = conf['host']
    port = int(conf['port'])
    # broj tacaka u jednoj kreiranoj putanji
    pathPoints_num = int(conf['brojTacakaUJednojPutanji'])
    # broj cestica
    particles_num = int(conf['brojCestica'])

    credentials = pika.PlainCredentials(conf['user'], conf['pass'])
    connection = pika.BlockingConnection(pika.ConnectionParameters(host=host, port=port, credentials=credentials))
    channel = connection.channel()

    # Posalji pocetne vrednosti cestice
    channel.queue_declare(queue='pocetne_vrednosti')

    for i in range(particles_num):
        half_life = distrib(tau)
        z = 1.e-2 * random()
        rr = r * sqrt(random())
        fir = 2 * pi * random()
        x = rr * cos(fir)
        y = rr * sin(fir)

        pocetneVrednostiCestice = '{"id" : ' + str(i) + ', "half_life" : ' + str(half_life) + ', "z" : ' + str(z) + ', "x" : ' + str(x)  \
                     + ', "y" : ' + str(y)  + ' }'

        channel.basic_publish(exchange='',
                              routing_key='pocetne_vrednosti',
                              body=pocetneVrednostiCestice,
                              properties=pika.BasicProperties(
                                  delivery_mode=1
                              ))
        #print str(x) + ' - ' + str(y) + ' - ' + str(z) + ' - ' + str(half_life)

    args = {'x-max-length' : particles_num * 3} # Maksimalan broj poruka u redu putanja
    channel.queue_declare(queue='path', arguments=args)

    # u posebnoj niti cekaj na poruku od RadonCounter.py da je simulacija zavrsena
    thread = Thread(target = CekajNaPorukuZaKraj)
    thread.start()

    speed_maxwell_List = None
    v = 0.0

    try:
        while running:

            xList = []
            yList = []
            zList = []
            life_time_step_List = []

            # kreiraj putanju od pathPoints_num tacaka
            for i in range(pathPoints_num):
                impact_distance = distrib(half_path)
                fi0 = 2 * pi * random()
                theta0 = acos(1 - 2 * random())

                xList.append(impact_distance * sin(theta0) * cos(fi0))
                yList.append(impact_distance * sin(theta0) * sin(fi0))
                zList.append(impact_distance * cos(theta0))

                v = speed_maxwell(0.2, 293.0)

                v_magnitude = sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

                life_time_step = impact_distance / v_magnitude
                life_time_step_List.append(life_time_step)

            points = '{ "x" : ' + str(xList) + ', "y" : ' + str(yList) + ', "z" : ' + str(zList)  \
                     + ', "life_time_step" : ' + str(life_time_step_List)  + ' }'

            # posalji putanju u json formatu
            channel.basic_publish(exchange='',
                                routing_key='path',
                                body=points,
                                properties=pika.BasicProperties(
                                    delivery_mode=1
                                ))

        channel.queue_delete(queue='pocetne_vrednosti')
        channel.queue_delete(queue='path')
        channel.queue_delete(queue='particleCounter')
        channel.queue_delete(queue='end')
        connection.close()

    except KeyboardInterrupt:
        channel.queue_delete(queue='pocetne_vrednosti')
        channel.queue_delete(queue='path')
        channel.queue_delete(queue='particleCounter')
        channel.queue_delete(queue='end')
        connection.close()
