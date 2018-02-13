#include "particle.h"

void * consumerListener(void * p);
void initialParticlePositionObserver(Conf * conf, MPI_Datatype * dt_particle, int * running);
void runQueue(Conf * conf, MPI_Datatype * dt_particle, MPI_Datatype * dt_point, int consumerCount, int producerCount);
void simulateParticles(MPI_Datatype * dt_particle, MPI_Datatype * dt_point, int rank, Conf * conf);
void generatePathTrajektories(int rank, MPI_Datatype  *dt_point,int partialTrajectoryLength);

void main(int argc, char *argv[]) 
{
    Conf * conf = readConf("conf.json");
    // mpi
    int rank, size; 
    MPI_Datatype dt_point;
    MPI_Datatype dt_particle;

    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if(conf->consumerCount + conf->producerCount + 1 == size)
    {
        initialize();

        MPI_Type_contiguous(11, MPI_DOUBLE, &dt_point);
        MPI_Type_commit(&dt_point);
        
        MPI_Type_contiguous(7, MPI_DOUBLE, &dt_particle);
        MPI_Type_commit(&dt_particle);

        Partial_trajectory partial_trajectory;

        if(rank == 0) // Queue
        {   
            runQueue(conf,&dt_particle, &dt_point, conf->consumerCount, conf->producerCount);
        }
        else if(rank <= conf->consumerCount) // Consumers
        {
            simulateParticles(&dt_particle, &dt_point, rank, conf);
        }
        else if(conf->consumerCount < rank) // Producers
        {
            generatePathTrajektories(rank, &dt_point, conf->partialTrajectoryLength);
        }
    }
    else if(rank == 0)
    {
        printf("Simulacija nije izvrsena zbog pogresnih vrednosti parametara u fajlu conf.json. \n");
        printf("Zbir broja potrosaca i broja proizvodjaca mora biti jednak \"broj_procesora - 1\n");
    }

    MPI_Finalize(); 
}

void * consumerListener(void * v)
{
    Conf * conf = readConf("conf.json");
    
    int * running = (int *)v;
    
    int nextParticleIndex = 0;
    MPI_Status status;
    int consumerRank;
    
    MPI_Datatype dt_point;
    MPI_Type_contiguous(11, MPI_DOUBLE, &dt_point);
    MPI_Type_commit(&dt_point);
    
    MPI_Datatype dt_particle;
    MPI_Type_contiguous(7, MPI_DOUBLE, &dt_particle);
    MPI_Type_commit(&dt_particle);
    
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
                MPI_Send(&(particles[nextParticleIndex]), 1, dt_particle, array_of_statuses[i].MPI_SOURCE, PARTICLE, MPI_COMM_WORLD);
                MPI_Irecv(&buffer[array_of_indices[i]], 1, MPI_INT, array_of_statuses[i].MPI_SOURCE, PARTICLE_REQUEST, MPI_COMM_WORLD, &requests[array_of_indices[i]]);
                nextParticleIndex++;
           }
           else
           {
               Particle p;
               p.status = NONE;
               MPI_Send(&p, 1, dt_particle, array_of_statuses[i].MPI_SOURCE, PARTICLE, MPI_COMM_WORLD);
               consumerCount--;;
           }
        }
        
        sleep(1);
    }
    
    *running = 0;
}

/*
    Implementira red koji opsluzuje simulatore pocetnim pozicijama cestica.
 */
void initialParticlePositionObserver(Conf * conf, MPI_Datatype * dt_particle, int * running)
{
    pthread_t thread;
    pthread_attr_t attr;
    
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    pthread_create(&thread, &attr, consumerListener, running);
}

void runQueue(Conf * conf, MPI_Datatype * dt_particle, MPI_Datatype * dt_point, int consumerCount, int producerCount)
{  
    int running = 1;
    
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
    
    initialParticlePositionObserver(conf, dt_particle, &running);
    
    // Consumers
    for(i = 0; i < consumerCount; i++)
    {   
        MPI_Irecv(&bufferConsumerTrajectoryRequest[i], 1, MPI_INT, i + 1, PARTIAL_TRAJECTORY_REQUEST, MPI_COMM_WORLD, &consumerRequests[i]);
    }
    
    // Producers
    for(i = consumerCount, j = 0; i < consumerCount + producerCount; i++, j++)
    {   
        MPI_Irecv(&bufferProducer[j], 1, *dt_point, i + 1, PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &producerRequests[j]);
    }
    
    int indexProducer = -1;
    int indexConsumer = -1;
    
    while(running)
    {
        // Cekaj na putanje
        while(outCountProducers == 0 && running)
            MPI_Testsome(producerCount, producerRequests, &outCountProducers, array_of_indices_producers, producerStatuses); 
        
        // Sve dok ne posaljes sve putanje
        while(outCountProducers && running)
        {
            // Pokupi sve zahteve za putanjama
            while( outCountConsumers == 0 && running)
            {
                MPI_Testsome(consumerCount, consumerRequests, &outCountConsumers, array_of_indices_consumers, consumerStatuses); 
            }
            
            // Posalji sve putanje
            while(outCountProducers && outCountConsumers && running)
            {   
                outCountProducers--;
                outCountConsumers--;
            
                indexProducer = array_of_indices_producers[outCountProducers];
                indexConsumer = array_of_indices_consumers[outCountConsumers];
                        
                partial_trajectory = &bufferProducer[indexProducer];
                
                MPI_Send(partial_trajectory, 1, *dt_point, consumerStatuses[outCountConsumers].MPI_SOURCE, PARTIAL_TRAJECTORY, MPI_COMM_WORLD);

            
                MPI_Irecv(&bufferConsumerTrajectoryRequest[indexConsumer], 1, MPI_INT, consumerStatuses[outCountConsumers].MPI_SOURCE, PARTIAL_TRAJECTORY_REQUEST, MPI_COMM_WORLD, &consumerRequests[indexConsumer]);
            
                MPI_Irecv(&bufferProducer[indexProducer], 1, *dt_point, producerStatuses[outCountProducers].MPI_SOURCE, PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &producerRequests[indexProducer]);
            }
        }
    }

    // Posalji poruke za kraj
    int n = 1;
    for(i = conf->consumerCount + 1; i <= conf->consumerCount + conf->producerCount; i++)
        MPI_Send(&n, 1, MPI_INT, i, SIMULATION_END, MPI_COMM_WORLD);
    
    printf("Deljeni red se odjavljuje. \n");
    
}

void simulateParticles(MPI_Datatype * dt_particle, MPI_Datatype * dt_point, int rank, Conf * conf)
{
    Partial_trajectory partialTrajectory;
    Particle particle;
    MPI_Request request;
    MPI_Status status;
    
    int count;
    
    char fileName[30];
    sprintf(fileName, "rezultat_%d_%d_%d.csv", conf->particleCount, conf->producerCount, conf->consumerCount);
    FILE * f = fopen(fileName, "a");
    
    clock_t begin = clock();
    
    while(1)
    {
        //Posalji zahtev za novom cesticom
        MPI_Send(&rank, 1, MPI_INT, 0, PARTICLE_REQUEST, MPI_COMM_WORLD);

        // Prihvati cesticu
        MPI_Recv(&particle, 1, *dt_particle, 0, PARTICLE, MPI_COMM_WORLD, &status);
        
        // Ako je status cestice NONE prekini sa radom
        if(particle.status == NONE) break;
        
        count = 0;

        while(1)
        {
            // Posalji zahtev za parcijalnom putanjom
            MPI_Send(&rank, 1, MPI_INT, 0, PARTIAL_TRAJECTORY_REQUEST, MPI_COMM_WORLD);
            
            // Prihvati parcijalnu putanju
            MPI_Recv(&partialTrajectory, 1, *dt_point, 0, PARTIAL_TRAJECTORY, MPI_COMM_WORLD, &status);
 
            if( checkBoundingBox(&particle, &partialTrajectory) != AIR)
            {
                if(checkDepositionPlace(&particle, conf->partialTrajectoryLength, &partialTrajectory.seed, &count) == 1)
                {
                    switch(particle.status)
                    {
                        case TOP:
                            printf("TOP, %d, %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "TOP, %d, %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                        break;
                        case BOTTOM:
                            printf("BOTTOM, %d, %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "BOTTOM, %d,  %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                        break;
                        case WALL:
                            printf("WALL, %d, %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "WALL, %d,  %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                        break;
                        case AIR:
                            printf("AIR, %d, %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                            fprintf(f, "AIR, %d,  %lf\n", count, (double)(clock() - begin) / CLOCKS_PER_SEC);
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
    
    printf("Simulator cestice %d. se odjavljuje\n", rank);
    fclose(f);
}


void generatePathTrajektories(int rank, MPI_Datatype  *dt_point,int partialTrajectoryLength)
{
    int running = 1;
    
    /*
        Wait for message to end in separate thread.
        Thread will set variable "running" to 0.
    */
    
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
    
    while(running)
    {
        Partial_trajectory * partial_trajectory = 
                    generatePartialTrajectory(&gauss_kon, &gauss_y, partialTrajectoryLength, &buffer_seed);
            
        MPI_Send(partial_trajectory, 1, *dt_point, 0, PARTIAL_TRAJECTORY, MPI_COMM_WORLD);  
        free(partial_trajectory);
    }
    
    printf("Proizvodjac %d se odjavljuje\n", rank);
}
