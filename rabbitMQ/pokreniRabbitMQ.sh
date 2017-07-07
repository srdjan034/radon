#!/bin/bash

particles_num=`python procitajJsonPolje.py configuration.json brojCestica`       # ukupan broj cestica
producers_num=`python procitajJsonPolje.py configuration.json brojProizvodjaca`  # ukupan broj proizvojaca

# pokreni skripta koja obradjuju po jednu cesticu
for id in `seq 1 $particles_num`
do
    python PathConsumer.py $id &
done

# pokreni proizvodjace putanja
for i in `seq 1 $producers_num`
do
    python PathProducer.py &
done

echo 'Simulacija je u toku...'
# pokreni program koji sabira sve statuse cestica
time python RadonCounter.py


