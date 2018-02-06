#include "particle.h"

void * consumerListener(void * p);
void initialParticlePositionObserver(Conf * conf, MPI_Datatype * dt_particle, int * running);
void runQueue(Conf * conf, MPI_Datatype * dt_particle, MPI_Datatype * dt_point, int consumerCount, int producerCount);
void simulateParticles(MPI_Datatype * dt_particle, MPI_Datatype * dt_point, int rank, Conf * conf);
void generatePathTrajektories(int rank, MPI_Datatype  *dt_point,int partialTrajectoryLength, struct drand48_data * buffer_seed);

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
        else if(0 < rank && rank <= conf->consumerCount) // Consumers
        {
            simulateParticles(&dt_particle, &dt_point, rank, conf);
        }
        else if(conf->consumerCount < rank) // Producers
        {
            struct drand48_data buffer_seed;
            srand48_r(rank, &buffer_seed);
            generatePathTrajektories(rank, &dt_point, conf->partialTrajectoryLength, &buffer_seed);
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
    
    MPI_Datatype * dt_point;
    MPI_Type_contiguous(11, MPI_DOUBLE, &dt_point);
    MPI_Type_commit(&dt_point);
    
    MPI_Datatype * dt_particle;
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
        MPI_Irecv(&buffer[i], 1, MPI_INT, i+1, 4, MPI_COMM_WORLD, &requests[i]);
    } 
    
    while(consumerCount > 0)
    {
        
        MPI_Testsome(conf->consumerCount, requests, &outcount, array_of_indices, array_of_statuses); 

        // For every completed requests
        for(i = 0; i < outcount; i++)
        {
           consumerRank = buffer[array_of_indices[i]];
           
           if(nextParticleIndex < particleCount)
           {
                MPI_Send(&particles[nextParticleIndex++], 1, dt_particle, consumerRank, 4, MPI_COMM_WORLD);
                MPI_Irecv(&buffer[array_of_indices[i]], 1, MPI_INT, array_of_statuses[i].MPI_SOURCE, 4, MPI_COMM_WORLD, &requests[i]);
           }
           else
           {
               Particle p;
               p.status = NONE;
               MPI_Send(&p, 1, dt_particle, consumerRank, 4, MPI_COMM_WORLD);
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
    initialParticlePositionObserver(conf, dt_particle, &running);
    
    Partial_trajectory * partial_trajectory;
    
    // Consumer
    MPI_Request array_of_consumer_requests[consumerCount];
    MPI_Status array_of_consumer_statuses[consumerCount];
    int outConsumerCount = 0;
    int array_of_consumer_indices[consumerCount];
    int bufferConsumer[consumerCount];
    
    // Producer
    MPI_Request array_of_producer_requests[producerCount];
    MPI_Status array_of_producer_statuses[producerCount];
    int outProducerCount = 0;
    int array_of_producer_indices[producerCount];
    Partial_trajectory bufferProducer[producerCount];
    
    MPI_Request request;
    MPI_Status status;
    int i;
    int flag;
    
    // Consumers 
    for(i = 0; i < consumerCount; i++)
    {   
        MPI_Irecv(&bufferConsumer[i], 1, MPI_INT, i + 1, 1, MPI_COMM_WORLD, &array_of_consumer_requests[i]);
    }
    
    // Producers
    for(i = consumerCount; i < consumerCount + producerCount; i++)
    {   
        MPI_Irecv(&bufferProducer[i-consumerCount], 1, *dt_point, i + 1, 2, MPI_COMM_WORLD, &array_of_producer_requests[i-consumerCount]);
    }
    
    while(running)
    {
        // Cekaj na zahteve za putanjama od potrosaca
        while(outConsumerCount == 0 && running)
        {
            MPI_Testsome(consumerCount, array_of_consumer_requests, &outConsumerCount, array_of_consumer_indices, array_of_consumer_statuses);
        }
        
        // Sve dok ne obradis sve zahteve za putanjama
        while(outConsumerCount > 0 && running)
        {           
            // Cekaj na putanje od proizvodjaca
            while(outProducerCount == 0 && running)
               MPI_Testsome(producerCount, array_of_producer_requests, &outProducerCount, array_of_producer_indices, array_of_producer_statuses); 
            
            // Salji putanje
            while(outConsumerCount > 0 && outProducerCount > 0 && running)
            {
                outConsumerCount--;
                outProducerCount--;
                
                // Preuzmi putanju i posalji je
                partial_trajectory = &bufferProducer[array_of_producer_indices[outProducerCount]];

                int consumerRank = bufferConsumer[array_of_consumer_indices[outConsumerCount]];
                
                MPI_Isend(partial_trajectory, 1, *dt_point, consumerRank, 3, MPI_COMM_WORLD, &request);
                
                // Prijavi se za cekanje potrosaca da se javi da li mu treba jos putanja
                MPI_Irecv(&bufferConsumer[array_of_consumer_indices[outConsumerCount]], 
                        1, MPI_INT, consumerRank, 1, MPI_COMM_WORLD, 
                        &array_of_consumer_requests[outConsumerCount]);
               
                // Prijavi se za prihvatanje jos putanja od proizvodjaca cija je putanja upravo obradjena
                MPI_Irecv(&bufferProducer[array_of_producer_indices[outProducerCount]], 
                            1, *dt_point, array_of_producer_statuses[outProducerCount].MPI_SOURCE, 2, MPI_COMM_WORLD, 
                            &array_of_producer_requests[outProducerCount]);
            }
        }
    }
    
    // Posalji poruke za kraj
    int n = 1;
    for(i = conf->consumerCount + 1; i <= conf->consumerCount + conf->producerCount; i++)
        MPI_Send(&n, 1, MPI_INT, i, 5, MPI_COMM_WORLD);
    
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
    sprintf(fileName, "rezultat.csv", rank);
    FILE * f = fopen(fileName, "a");
    
    clock_t begin = clock();
    
    while(1)
    {
        // Posalji zahtev za novom cesticom
        MPI_Isend(&rank, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, &request);

        // Prihvati cesticu
        MPI_Recv(&particle, 1, *dt_particle, 0, 4, MPI_COMM_WORLD, &status);
        
        // Ako je status cestice NONE prekini saradom
        if(particle.status == NONE) break;
        
        count = 0;

        while(1)
        {
            // Posalji zahtev za parcijalnom putanjom
            MPI_Send(&rank, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            
            // Prihvati parcijalnu putanju
            MPI_Recv(&partialTrajectory, 1, *dt_point, 0, 3, MPI_COMM_WORLD, &status);
 
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


void generatePathTrajektories(int rank, MPI_Datatype  *dt_point,int partialTrajectoryLength, struct drand48_data * buffer_seed)
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
    
    MPI_Request request;
    int gauss_kon = 1;
    double gauss_y = 0.0;
    
    while(running)
    {
        Partial_trajectory * partial_trajectory = 
                    generatePartialTrajectory(&gauss_kon, &gauss_y, partialTrajectoryLength, buffer_seed);
            
        MPI_Isend(partial_trajectory, 1, *dt_point, 0, 2, MPI_COMM_WORLD, &request);
            
        free(partial_trajectory);
    }
    
    printf("Proizvodjac %d se odjavljuje\n", rank);
}
