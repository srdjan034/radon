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
    int running = 1;
    
    /*
        Wait for message to end in separate thread.
        Thread will set variable "running" to 0.
    */
    
    MPI_Datatype dt_partial_trajectory;
    createMpiTypeForPartialTrajectory(&dt_partial_trajectory);
    
    pthread_t thread;
    pthread_attr_t attr;
    
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        
    pthread_create(&thread, &attr, waitForSignal, (void *) & running);
    
    struct drand48_data buffer_seed;
    srand48_r(rank, &buffer_seed);
    
    MPI_Request request;

    int gauss_kon = 1;
    double gauss_y = 0.0;
    
    int pathCount = 0;

    while(running)
    {
        Partial_trajectory partial_trajectory;
        generatePartialTrajectory(&partial_trajectory, &gauss_kon, &gauss_y, partialTrajectoryLength, &buffer_seed);
        
        MPI_Isend(&partial_trajectory, 1, dt_partial_trajectory, 0, PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &request); 

        pathCount++;
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
    
    int count;
    
    char fileName[30];
    sprintf(fileName, "rezultat_%d_%d_%d_%d.csv", conf->particleCount, conf->producerCount, conf->consumerCount, rank);
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
        
        count = 0;

        while(1)
        {
            // Posalji zahtev za parcijalnom putanjom
            MPI_Isend(&rank, 1, MPI_INT, 0, PARTIAL_TRAJECTORY_REQUEST, MPI_COMM_WORLD, &request);
            
            // Prihvati parcijalnu putanju
            MPI_Recv(&partialTrajectory, 1, dt_partial_trajectory, 0, PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &status);
             
            if( checkBoundingBox(&particle, &partialTrajectory) != AIR)
            {  
                struct drand48_data seed;
                
                initializeSeed(&seed, &partialTrajectory);
                
                if(checkDepositionPlace(&particle, conf->partialTrajectoryLength, &seed, &count) == 1)
                {
                    switch(particle.status)
                    {
                        case TOP:
                            printf("TOP, %d, %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "TOP, %d, %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            nTop++;
                        break;
                        case BOTTOM:
                            printf("BOTTOM, %d, %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "BOTTOM, %d,  %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            nBottom++;
                        break;
                        case WALL:
                            printf("WALL, %d, %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "WALL, %d,  %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            nWall++;
                        break;
                        case AIR:
                            printf("AIR, %d, %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "AIR, %d,  %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
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
                    count += conf->partialTrajectoryLength;
                }
            }
            else
            {
                particle.life_time += partialTrajectory.lifeTimeStepSum;
                particle.x += partialTrajectory.xLast;
                particle.y += partialTrajectory.yLast;
                particle.z += partialTrajectory.zLast;
                count += conf->partialTrajectoryLength;
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
    MPI_Request request;
    
    MPI_Datatype dt_partial_trajectory;
    createMpiTypeForPartialTrajectory(&dt_partial_trajectory);
    
    Partial_trajectory * partial_trajectory;
    
    MPI_Request consumerRequests[consumerCount];
    MPI_Status  consumerStatuses[consumerCount];
    int outCountConsumers = 0;
    int array_of_indices_consumers[consumerCount];
 
    MPI_Request producerRequests[producerCount];
    MPI_Status  producerStatuses[producerCount];
    int outCountProducers = 0;
    int array_of_indices_producers[producerCount];

    int bufferConsumerTrajectoryRequest[consumerCount];
    Partial_trajectory bufferProducer[producerCount];
    
    int i, j;
    
    initialParticlePositionObserver(conf, &running);
    
    // Consumers
    for(i = 0; i < consumerCount; i++)
    {   
        MPI_Irecv(&bufferConsumerTrajectoryRequest[i], 1, MPI_INT, i + 1, PARTIAL_TRAJECTORY_REQUEST, MPI_COMM_WORLD, &consumerRequests[i]);
    }

    // Producers
    for(i = consumerCount, j = 0; i < consumerCount + producerCount; i++, j++)
    {   
        MPI_Irecv(&bufferProducer[j], 1, dt_partial_trajectory, i + 1, PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &producerRequests[j]);
    }
    
    int indexProducer = -1;
    int indexConsumer = -1;
    
    int processedPathCount = 0;
    
    while(running)
    {
       
        // Cekaj na putanje
        while(outCountProducers == 0 && running)
        {   
            MPI_Testsome(producerCount, producerRequests, &outCountProducers, array_of_indices_producers, producerStatuses); 

            if(outCountProducers == 0) 
                usleep(100);
        }
        
        // Sve dok ne posaljes sve putanje
        while(outCountProducers && running)
        {
            // Pokupi sve zahteve za putanjama
            while( outCountConsumers == 0 && running)
            {
                MPI_Testsome(consumerCount, consumerRequests, &outCountConsumers, array_of_indices_consumers, consumerStatuses);
                
                if(outCountConsumers == 0) 
                    usleep(100);
            }
            
            // Posalji sve putanje
            while(outCountProducers && outCountConsumers && running)
            {   
                outCountProducers--;
                outCountConsumers--;
            
                indexProducer = array_of_indices_producers[outCountProducers];
                indexConsumer = array_of_indices_consumers[outCountConsumers];

                partial_trajectory = &bufferProducer[indexProducer];
                
                processedPathCount++;
                
                MPI_Isend(partial_trajectory, 1, dt_partial_trajectory, consumerStatuses[outCountConsumers].MPI_SOURCE, PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &request);
            
                MPI_Irecv(&bufferConsumerTrajectoryRequest[indexConsumer], 1, MPI_INT, consumerStatuses[outCountConsumers].MPI_SOURCE, PARTIAL_TRAJECTORY_REQUEST, MPI_COMM_WORLD, &consumerRequests[indexConsumer]);
            
                MPI_Irecv(&bufferProducer[indexProducer], 1, dt_partial_trajectory, producerStatuses[outCountProducers].MPI_SOURCE, PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &producerRequests[indexProducer]);
            }
        }
    }

    // Posalji poruke za kraj
    int n = 1;
    for(i = conf->consumerCount + 1; i <= conf->consumerCount + conf->producerCount; i++)
        MPI_Send(&n, 1, MPI_INT, i, SIMULATION_END, MPI_COMM_WORLD);
    
    getDepositionStatistics(conf);
    
    writeTrajectoryStatistics(0, processedPathCount, 0);
    
}

void * consumerListener(void * v)
{
    Conf * conf = readConf("conf.json");
    
    int * running = (int *)v;
    
    int nextParticleIndex = 0;
    MPI_Request request;
    
    MPI_Datatype dt_particle;
    createMpiTypeForParticle(&dt_particle);
    
    int buffer[conf->consumerCount];
    MPI_Request requests[conf->consumerCount];
    int outcount;
    int array_of_indices[conf->consumerCount];
    MPI_Status array_of_statuses[conf->consumerCount]; 
     
    int particleCount = conf->particleCount;
    int consumerCount = conf->consumerCount;
    
    Particle particles[conf->particleCount]; 
    
    struct drand48_data buffer_seed;
    srand48_r(time(NULL), &buffer_seed);
    
    int i;
    for(i = 0; i < conf->particleCount; i++)
        initialize_particle(&particles[i], i, tau, r, &buffer_seed);

    for(i = 0; i < conf->consumerCount; i++)
    {   
        MPI_Irecv(&buffer[i], 1, MPI_INT, i+1, PARTICLE_REQUEST, MPI_COMM_WORLD, &requests[i]);
    } 
    
    while(consumerCount > 0)
    {
        MPI_Testsome(conf->consumerCount, requests, &outcount, array_of_indices, array_of_statuses);
        
        // For every completed requests
        for(i = 0; i < outcount; i++)
        {  
           if(nextParticleIndex < particleCount)
           {
                MPI_Isend(&(particles[nextParticleIndex]), 1, dt_particle, array_of_statuses[i].MPI_SOURCE, PARTICLE, MPI_COMM_WORLD, &request);
                MPI_Irecv(&buffer[array_of_indices[i]], 1, MPI_INT, array_of_statuses[i].MPI_SOURCE, PARTICLE_REQUEST, MPI_COMM_WORLD, &requests[array_of_indices[i]]);
                nextParticleIndex++;
           }
           else
           {
               Particle p;
               p.status = NONE;
               MPI_Isend(&p, 1, dt_particle, array_of_statuses[i].MPI_SOURCE, PARTICLE, MPI_COMM_WORLD, &request);
               consumerCount--;;
           }
        }
        
        outcount = 0;
        sleep(1);
    }
    
    *running = 0;
}

/*
    Implementira red koji opsluzuje simulatore pocetnim pozicijama cestica.
 */
void initialParticlePositionObserver(Conf * conf, int * running)
{
    pthread_t thread;
    pthread_attr_t attr;
    
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    pthread_create(&thread, &attr, consumerListener, running);
}