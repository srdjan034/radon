import pika
import sys
import subprocess
import os
import json

nAir = nWall = nTop = nBottom = 0
particles_num = 0
particles_forced = 0
connection = None
channel = None

def print_status(nAir, nWall, nTop, nBottom):
	print ("Bottom=%d Top=%d Wall=%d Air=%d" % (nBottom, nTop, nWall, nAir))

def counter(ch, method, properties, poruka):

    global nAir
    global nWall
    global nTop
    global nBottom
    global particles_num
    global producers_num

    if(poruka == "AIR"):
        nAir += 1
    elif(poruka == "WALL"):
        nWall += 1
    elif (poruka == "Top"):
        nTop += 1
    else:
        nBottom += 1

    print_status(nAir, nWall, nTop, nBottom)

    # Ako su cestice zavrsile
    if(nAir + nWall + nTop + nBottom == particles_num):
        print "Kraj simulacije"
        channel.queue_declare(queue='end')
        # Posalji svim proizvodjacima da je simulacija zavrsena
        for i in range(producers_num):
            channel.basic_publish(exchange='',
                              routing_key='end',
                              body='kraj',
                              properties=pika.BasicProperties(
                                  delivery_mode=1,  # make message persistent
                              ))

        # Odjavi se sa svog kanala
        ch.basic_ack(delivery_tag=method.delivery_tag)
        ch.basic_cancel(consumer_tag='particleCounter')

        # Obrisi kanale
        ch.queue_delete(queue='particleCounter')
        ch.queue_delete(queue='end')
        sys.exit(0)

if __name__ == "__main__":

    try:

        with open('configuration.json') as data_file:
            conf = json.load(data_file)

        host = conf['host']
        port = int(conf['port'])
        particles_num = int(conf['brojCestica'])
        producers_num = int(conf['brojProizvodjaca'])

        connection = pika.BlockingConnection(pika.ConnectionParameters(host=host, port=port))
        channel = connection.channel()
        channel.queue_declare(queue='particleCounter')

        channel.basic_consume(counter,
                                queue='particleCounter',
                                no_ack=True)

        channel.start_consuming()

    except KeyboardInterrupt:
        channel.queue_delete(queue='particleCounter')
        channel.queue_delete(queue='end')
        connection.close()