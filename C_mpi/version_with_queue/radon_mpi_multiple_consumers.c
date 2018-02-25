#include "particle.h"

int main(int argc, char *argv[]) 
{
    Conf * conf = readConf("conf.json");
    
    int rank, size; 

    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if(conf->consumerCount + conf->producerCount + 1 == size)
    {
        initialize();
        
        if(rank == 0) // Queue
        {   
            runQueue(conf, conf->consumerCount, conf->producerCount);
        }
        else if(rank <= conf->consumerCount) // Consumers
        {
            simulateParticles(rank, conf);
        }
        else if(conf->consumerCount < rank) // Producers
        {
            generatePathTrajektories(rank, conf->partialTrajectoryLength);
        }
    }
    else if(rank == 0)
    {
        printf("Simulacija nije izvrsena zbog pogresnih vrednosti parametara u fajlu conf.json. \n");
        printf("Zbir broja potrosaca i broja proizvodjaca mora biti jednak \"broj_procesora - 1\n");
    }

    MPI_Finalize();
    
    return 0;
}

void generatePathTrajektories(int rank, int partialTrajectoryLength)
{    
    MPI_Datatype dt_partial_trajectory;
    createMpiTypeForPartialTrajectory(&dt_partial_trajectory);
    
    struct drand48_data buffer_seed;
    srand48_r(rank, &buffer_seed);
    
    MPI_Status status;
    MPI_Request request;

    int gauss_kon = 1;
    double gauss_y = 0.0;
    
    int pathCount = 0;

    int endOfSimulationFlag = 0;
    
    // Sve dok ne stigne poruka za kraj simulacije
    while(endOfSimulationFlag == 0)
    {
        Partial_trajectory partial_trajectory;
        generatePartialTrajectory(&partial_trajectory, &gauss_kon, &gauss_y, 
                partialTrajectoryLength, &buffer_seed);
        
        MPI_Isend(&partial_trajectory, 1, dt_partial_trajectory, 0, 
                PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &request); 

        pathCount++;
        
        // Proveri da li je stigla poruka za kraj simulacije
        MPI_Iprobe(0, SIMULATION_END, MPI_COMM_WORLD, &endOfSimulationFlag, &status);
    }
    
    writeTrajectoryStatistics(pathCount, 0, rank);
}

void simulateParticles( int rank, Conf * conf)
{
    Partial_trajectory partialTrajectory;
    Particle particle;
    MPI_Request request;
    MPI_Status status;
    
    MPI_Datatype dt_particle;         
    createMpiTypeForParticle(&dt_particle);
    
    MPI_Datatype dt_partial_trajectory;
    createMpiTypeForPartialTrajectory(&dt_partial_trajectory); 
    
    int nAir, nBottom, nWall, nTop;
    nAir = nBottom = nWall = nTop = 0;
    
    int particleCollisionCount = 0;
    
    char fileName[30];
    sprintf(fileName, "rezultat_%d_%d_%d_%d.csv", conf->particleCount, 
            conf->producerCount, conf->consumerCount, rank);
    
    FILE * f = fopen(fileName, "w");
    
    clock_t begin = clock();
    
    while(1)
    {
        //Posalji zahtev za novom cesticom
        MPI_Isend(&rank, 1, MPI_INT, 0, PARTICLE_REQUEST, MPI_COMM_WORLD, &request);

        // Prihvati cesticu
        MPI_Recv(&particle, 1, dt_particle, 0, PARTICLE, MPI_COMM_WORLD, &status);
        
        // Ako je status cestice NONE prekini sa radom
        if(particle.status == NONE) break;
        
        particleCollisionCount = 0;

        while(1)
        {
            // Posalji zahtev za parcijalnom putanjom
            MPI_Isend(&rank, 1, MPI_INT, 0, PARTIAL_TRAJECTORY_REQUEST, 
                    MPI_COMM_WORLD, &request);

            // Prihvati parcijalnu putanju
            MPI_Recv(&partialTrajectory, 1, dt_partial_trajectory, 0, 
                    PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &status);            
             
            if( checkBoundingBox(&particle, &partialTrajectory) != AIR)
            {  
                struct drand48_data seed;
                
                initializeSeed(&seed, &partialTrajectory);
                
                if(checkDepositionPlace(&particle, conf->partialTrajectoryLength, &seed, &particleCollisionCount) == 1)
                {
                    switch(particle.status)
                    {
                        case TOP:
                            printf("TOP, %d, %lf\n", particleCollisionCount, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "TOP, %d, %lf\n", particleCollisionCount, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            nTop++;
                        break;
                        case BOTTOM:
                            printf("BOTTOM, %d, %lf\n", particleCollisionCount, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "BOTTOM, %d,  %lf\n", particleCollisionCount, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            nBottom++;
                        break;
                        case WALL:
                            printf("WALL, %d, %lf\n", particleCollisionCount, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "WALL, %d,  %lf\n", particleCollisionCount, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            nWall++;
                        break;
                        case AIR:
                            printf("AIR, %d, %lf\n", particleCollisionCount, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "AIR, %d,  %lf\n", particleCollisionCount, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            nAir++;
                        break;
                    }
                    
                    fflush(f);
                    
                    break;
                }
                else
                {
                    particle.life_time += partialTrajectory.lifeTimeStepSum;
                    particle.x += partialTrajectory.xLast;
                    particle.y += partialTrajectory.yLast;
                    particle.z += partialTrajectory.zLast;
                    particleCollisionCount += conf->partialTrajectoryLength;
                }
            }
            else
            {
                particle.life_time += partialTrajectory.lifeTimeStepSum;
                particle.x += partialTrajectory.xLast;
                particle.y += partialTrajectory.yLast;
                particle.z += partialTrajectory.zLast;
                particleCollisionCount += conf->partialTrajectoryLength;
            }
        }
    }
    
    fclose(f);
    
    sendDepositionStatistics(nAir, nBottom, nTop, nWall);
    writeTrajectoryStatistics(0, 0, rank);
}

void runQueue(Conf * conf, int consumerCount, int producerCount)
{  
    int running = 1;
    int i, j;
    MPI_Request request;
    int indexProducer = -1;
    int indexConsumer = -1;
    int processedPathCount = 0;

    /*
        Pripremi pocetne pozicije cestica i prijavi se da prihvatas zahteve
        pd simulatora vestica.
     */
    
    MPI_Datatype dt_particle;
    createMpiTypeForParticle(&dt_particle);
    
    int particleCountTmp = conf->particleCount;
    int consumerCountTmp = conf->consumerCount;
    
    int nextParticleIndex = 0;
    
    int bufferParticles[conf->consumerCount];
    MPI_Request particleRequests[conf->consumerCount];
    int outCountParticles;
    int array_of_indices_particles[conf->consumerCount];
    MPI_Status array_of_statuses_particles[conf->consumerCount]; 
    
    struct drand48_data buffer_seed;
    srand48_r(time(NULL), &buffer_seed);
    
    // Kreiraj cestice
    Particle particles[conf->particleCount]; 
    for(i = 0; i < conf->particleCount; i++)
        initialize_particle(&particles[i], i, tau, r, &buffer_seed);
    
    // Prijavi se na zahteve za cesticama
    for(i = 0; i < conf->consumerCount; i++)
    {   
        MPI_Irecv(&bufferParticles[i], 1, MPI_INT, i+1, PARTICLE_REQUEST, 
                MPI_COMM_WORLD, &particleRequests[i]);
    } 
    
    /*
        Prijavi se da prihvatas zahteve za parcijalnim putanjama.
     */
    
    int bufferConsumerTrajectoryRequest[consumerCount];
    MPI_Request consumerRequests[consumerCount];
    MPI_Status  consumerStatuses[consumerCount];
    int outCountConsumers = 0;
    int array_of_indices_consumers[consumerCount];
    
    //
    for(i = 0; i < consumerCount; i++)
    {   
        MPI_Irecv(&bufferConsumerTrajectoryRequest[i], 1, MPI_INT, i + 1, 
                PARTIAL_TRAJECTORY_REQUEST, MPI_COMM_WORLD, &consumerRequests[i]);
    }

    /* 
        Prijavi se da prihvatas parcijalne putanje od proizvodjaca.
     */
    
    MPI_Datatype dt_partial_trajectory;
    createMpiTypeForPartialTrajectory(&dt_partial_trajectory);
    
    Partial_trajectory * partial_trajectory;
    
    Partial_trajectory bufferProducer[producerCount];
    MPI_Request producerRequests[producerCount];
    MPI_Status  producerStatuses[producerCount];
    int outCountProducers = 0;
    int array_of_indices_producers[producerCount];
    
    for(i = consumerCount, j = 0; i < consumerCount + producerCount; i++, j++)
    {   
        MPI_Irecv(&bufferProducer[j], 1, dt_partial_trajectory, i + 1, 
                PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &producerRequests[j]);
    }
    
    /*
        Simulacija reda pocinje ovde.
     */
    
    while(running)
    {
        // Cekaj na putanje od proizvodjaca
        if(outCountProducers == 0 && running)
            MPI_Waitsome(producerCount, producerRequests, &outCountProducers, 
                    array_of_indices_producers, producerStatuses); 
        
        // Sve dok ne posaljes sve putanje koje su pristigle naredbom iznad.
        while(outCountProducers && running)
        {
            // Proveri da li ima zahteva za parcijalnim putanjama
            if(outCountConsumers == 0 && running)
                MPI_Testsome(consumerCount, consumerRequests, &outCountConsumers, 
                        array_of_indices_consumers, consumerStatuses);
            
            // Ako ima zahteva za parcijalnim putanjama, obradi ih
            if(outCountConsumers > 0)
            {
                // Obradi sve zahteve za putanjama
                while(outCountProducers && outCountConsumers && running)
                {   
                    outCountProducers--;
                    outCountConsumers--;

                    indexProducer = array_of_indices_producers[outCountProducers];
                    indexConsumer = array_of_indices_consumers[outCountConsumers];

                    partial_trajectory = &bufferProducer[indexProducer];

                    processedPathCount++;

                    MPI_Isend(partial_trajectory, 1, dt_partial_trajectory, 
                            consumerStatuses[outCountConsumers].MPI_SOURCE, 
                            PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &request);

                    MPI_Irecv(&bufferConsumerTrajectoryRequest[indexConsumer], 
                            1, MPI_INT, consumerStatuses[outCountConsumers].MPI_SOURCE, 
                            PARTIAL_TRAJECTORY_REQUEST, 
                            MPI_COMM_WORLD, &consumerRequests[indexConsumer]);

                    MPI_Irecv(&bufferProducer[indexProducer], 1, dt_partial_trajectory, 
                            producerStatuses[outCountProducers].MPI_SOURCE, 
                            PARTIAL_TRAJECTORY, MPI_COMM_WORLD, 
                            &producerRequests[indexProducer]);
                }
            }
            
            /*
                Kada zavrsis sa obradom zahteva za putanjama, proveri da li ima zahteva za cesticama.
             */
            
            // Proveri da li ima zahteva za cesticama
            MPI_Testsome(conf->consumerCount, particleRequests, &outCountParticles, 
                    array_of_indices_particles, array_of_statuses_particles);
            
            // Obradi sve zahteve za cesticama
            for(i = 0; i < outCountParticles; i++)
            { 
               // Proveri da li jos ima cestica
               if(nextParticleIndex < particleCountTmp)
               {
                    MPI_Isend(&(particles[nextParticleIndex]), 1, dt_particle, 
                            array_of_statuses_particles[i].MPI_SOURCE, PARTICLE, 
                            MPI_COMM_WORLD, &request);

                    MPI_Irecv(&bufferParticles[array_of_indices_particles[i]], 1, MPI_INT, 
                            array_of_statuses_particles[i].MPI_SOURCE, PARTICLE_REQUEST, 
                            MPI_COMM_WORLD, &particleRequests[array_of_indices_particles[i]]);
                    
                    nextParticleIndex++;
               }
               else // Ako su sve cestice poslate, posalji cesticu sa statusom NONE
               {
                   Particle p;
                   p.status = NONE;
                   MPI_Isend(&p, 1, dt_particle, array_of_statuses_particles[i].MPI_SOURCE, 
                           PARTICLE, MPI_COMM_WORLD, &request);
                   
                   consumerCountTmp--;
                   
                   // Ako je svim simulatorima cestica poslat status NONE, zavrsi sa radom
                   if(consumerCountTmp == 0)
                       running = 0;
               }
            }
        }
    }

    // Posalji poruke za kraj
    int n = 1;
    for(i = conf->consumerCount + 1; i <= conf->consumerCount + conf->producerCount; i++)
        MPI_Isend(&n, 1, MPI_INT, i, SIMULATION_END, MPI_COMM_WORLD, &request);
    
    getDepositionStatistics(conf);
    
    writeTrajectoryStatistics(0, processedPathCount, 0);
    
}