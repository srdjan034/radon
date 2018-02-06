/* 
 * 
 * File:   radon_mpi.c
 * Author: srdjan
 *
 * Created on January 23, 2018, 1:06 PM
 * 
 */

#include "particle.h"

void main(int argc, char *argv[]) 
{
    Conf * conf = readConf("conf.json");
    
    // Particles
    Particle * particles;
    int particles_num = conf->particleCount; 
    int nPathPoints = conf->partialTrajectoryLength;
    int gauss_kon = 1;
    double gauss_y = 0.0;
    
    // mpi
    int rank, size; 
    double param; 
    MPI_Datatype dt_point;

    int running = 1;

    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    MPI_Request requests[size];
    
    initialize();
    
    struct drand48_data buffer_seed;
    srand48_r(rank, &buffer_seed);

    MPI_Type_contiguous(11, MPI_DOUBLE, &dt_point);
    MPI_Type_commit(&dt_point);

    Partial_trajectory partial_trajectory;

    if(rank == 0)
    {   
        /*
            Particles
         */
        int i, j;
        
        int nAir, nWall, nTop, nBottom;
        nAir = nWall = nTop = nBottom = 0;
        int nDeportedParticles = 0;
        
        Partial_trajectory * partial_trajectory;
        particles = (Particle *)malloc(particles_num * sizeof(Particle));
        for(i = 0; i < particles_num; i++)
            initialize_particle(&particles[i], i, tau, r, &buffer_seed);
        
        Particle * particle = &particles[nDeportedParticles];
        
        enum DEPORTATION_PLACE p;
        int count = 0;
        
        FILE * f = fopen("rezultat.csv", "w");
        clock_t begin = clock();
        
        /*
            MPI
         */
        MPI_Request array_of_requests[size - 1];
        int outcount;
        int array_of_indices[size - 1];
        MPI_Status array_of_statuses[size - 1];
        Partial_trajectory buffer[100];
        
        // Start a standard-mode, nonblocking receive. 
        for(i = 1; i < size; i++)
        {   
            MPI_Irecv(&buffer[i-1], 1, dt_point, i, 0, MPI_COMM_WORLD, &requests[i-1]);
        }
        
        while(nDeportedParticles != particles_num)
        {
            // Tests for some given requests to complete
            MPI_Testsome(size-1, requests, &outcount, array_of_indices, array_of_statuses); 
            
            // For every completed requests
            for(i = 0; i < outcount; i++)
            {
                partial_trajectory = &buffer[array_of_indices[i]];
                
                if( (p = checkBoundingBox(particle, partial_trajectory)) != AIR)
                {
                    
                    if(checkDepositionPlace(particle, nPathPoints, &partial_trajectory->seed, &count) == 1)
                    {
                        fprintf(f, "%d, %d, %lf\n", particle->id, count, (double)(clock() - begin) / CLOCKS_PER_SEC);
                        fflush(f);

                        nDeportedParticles++;

                        switch(particle->status)
                        {
                            case TOP:
                                nTop++;
                                break;
                            case BOTTOM:
                                nBottom++;
                                break;
                            case WALL:
                                nWall++;
                                break;
                            case AIR:
                                nAir++;
                                break;
                        }

                        print_status(nAir, nWall, nTop, nBottom);

                        if(nDeportedParticles == particles_num)
                            break;	
                        else
                        {
                            particle = &particles[nDeportedParticles];
                            count = 0;
                        }
                    }
                    else
                    {
                        particle->life_time += partial_trajectory->lifeTimeStepSum;
                        particle->x += partial_trajectory->xLast;
                        particle->y += partial_trajectory->yLast;
                        particle->z += partial_trajectory->zLast;
                    
                        count += nPathPoints;
                    }
                }
                else
                {		
                    particle->life_time += partial_trajectory->lifeTimeStepSum;
                    particle->x += partial_trajectory->xLast;
                    particle->y += partial_trajectory->yLast;
                    particle->z += partial_trajectory->zLast;
                    
                    count += nPathPoints;
                }
                
                MPI_Irecv(&buffer[array_of_indices[i]], 1, dt_point, array_of_indices[i]+1, 0, MPI_COMM_WORLD, &requests[array_of_indices[i]]);
            }
        }
        
        print_status(nAir, nWall, nTop, nBottom);
        fprintf(f, "nAir = %d, nWall = %d, nTop = %d, nBottom = %d \n", nAir, nWall, nTop, nBottom);
         
        fclose(f);

        // Posalji poruke za kraj
        int n = 1;
        for(i=1; i < size; i++)
            MPI_Send(&n, 1, MPI_INT, i, 5, MPI_COMM_WORLD);
        
    }
    else
    {
        /*
            Wait for message to end in separate thread.
            Thread will set variable "running" to 0.
         */
        pthread_t thread;
        pthread_attr_t attr;
    
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        
        pthread_create(&thread, &attr, waitForSignal, (void *) & running);
        
        /*
            Generate partial trajectories
         */
        while(running)
        {
            Partial_trajectory * partial_trajectory = 
                    generatePartialTrajectory(&gauss_kon, &gauss_y, nPathPoints, &buffer_seed);
            
            MPI_Isend(partial_trajectory, 1, dt_point, 0, 0, MPI_COMM_WORLD, &requests[rank]);
            
            free(partial_trajectory);
        }
    }

    MPI_Finalize(); 
}