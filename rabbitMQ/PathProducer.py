import pika
from Particle import *
from threading import Thread
import json
import sys
import time

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

    seed(1)

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

    args = {'x-max-length' : 40} # Maksimalan broj poruka u redu putanja
    channel.queue_declare(queue='path')#, arguments=args)

    # u posebnoj niti cekaj na poruku od RadonCounter.py da je simulacija zavrsena
    thread = Thread(target = CekajNaPorukuZaKraj)
    thread.start()

    speed_maxwell_List = None
    v = 0.0

    try:
        start_time = time.time()
        while running:

            x = 0
            y = 0
            z = 0

            xMin = sys.float_info.max
            xMax = sys.float_info.min
            yMin = sys.float_info.max
            yMax = sys.float_info.min
            zMin = sys.float_info.max
            zMax = sys.float_info.min

            life_time_step_sum = 0

            x_last = 0
            y_last = 0
            z_last = 0

            seed_state_list = getstate()

            # kreiraj putanju od pathPoints_num tacaka
            for j in range(pathPoints_num):

                impact_distance = distrib(half_path)
                fi0 = 2 * pi * random()
                theta0 = acos(1 - 2 * random())

                x = impact_distance * sin(theta0) * cos(fi0)
                y = impact_distance * sin(theta0) * sin(fi0)
                z = impact_distance * cos(theta0)

                x_last += x
                y_last += y
                z_last += z

                if x_last < xMin : xMin = x_last
                if y_last < yMin : yMin = y_last
                if z_last < zMin : zMin = z_last

                if x_last > xMax : xMax = x_last
                if y_last > yMax : yMax = y_last
                if z_last > zMax : zMax = z_last

                v = speed_maxwell(0.2, 293.0)

                v_magnitude = sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

                life_time_step = impact_distance / v_magnitude
                life_time_step_sum += life_time_step

            points = '{ "xMin" : ' + str(xMin) \
                     + ', "yMin" :' + str(yMin) \
                     + ', "zMin" :' + str(zMin) \
                     + ', "xMax" : ' + str(xMax) \
                     + ', "yMax" :' + str(yMax) \
                     + ', "zMax" :' + str(zMax) \
                     + ', "x_last" :' + str(x_last) \
                     + ', "y_last" :' + str(y_last) \
                     + ', "z_last" :' + str(z_last) \
                     + ', "life_time_step_sum" :' + str(life_time_step_sum) \
                     + ', "seed_state_0" :' + str(seed_state_list[0]) \
                     + ', "seed_state_1" : ' + ('[' + ', '.join(str(x) for x in seed_state_list[1]) + ']') \
                     + ', "seed_state_2" :' + (str(-1) if seed_state_list[2] is None else str(seed_state_list[2])) \
                     + ' }'

            # posalji putanju u json formatu
            channel.basic_publish(exchange='',
                               routing_key='path',
                                body=points,
                                properties=pika.BasicProperties(
                                    delivery_mode=2
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
