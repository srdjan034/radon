
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <float.h>
#include <unistd.h>
#include <mpi.h> 
#include <stddef.h>
#include <string.h>

// Job parameters
double r;
double h;
double half_path;
double tau;

enum DEPORTATION_PLACE
{
    INIT = -1,
    AIR = 0,
    WALL = 1,
    TOP = 2,
    BOTTOM = 3,
    DECOMPOSED = 4,
    NONE = 5
};

enum MESSAGE_TYPE
{
    PARTICLE_REQUEST = 10,
    PARTICLE = 11,
    PARTIAL_TRAJECTORY_REQUEST = 12,
    PARTIAL_TRAJECTORY = 31,
    SIMULATION_END = 14
};

typedef struct
{
    double life_time;
    double half_life;
    double z;
    double x;
    double y;
    
    int id;
    enum DEPORTATION_PLACE status;

} Particle;

typedef struct 
{
    double xMin;
    double xMax;
    double yMin;
    double yMax;
    double zMin;
    double zMax;
    double xLast;
    double yLast;
    double zLast;
    double lifeTimeStepSum;
    
    //drand48_data
    unsigned long long int __a;
    unsigned short int __c;
    unsigned short int __init;
    unsigned short int __old_x1;
    unsigned short int __old_x2;
    unsigned short int __old_x3;
    unsigned short int __x1;
    unsigned short int __x2;
    unsigned short int __x3;
    
} Partial_trajectory;

typedef struct
{
    int nodeCount;
    int processorsPerNode;
    int producerCount;
    int consumerCount;
    int partialTrajectoryLength;
    int particleCount;

} Conf;

void initialize()
{
    // Job parameters
    r = 0.04;
    h = 0.08;
    half_path = 6.0e-6;
    tau = 183/log(2);
}

double randomDoubleWithSeed(struct drand48_data * buffer_seed)
{
    double value;
    drand48_r( buffer_seed, &value);
    return value;
}

double distrib(double x, struct drand48_data * buffer_seed)
{
    return (-x*log(randomDoubleWithSeed(buffer_seed)));
}

double Gauss(int * gauss_kon, double * gauss_y, struct drand48_data * buffer_seed)
{
    double value;
    double s1, s2, temp, fi;

    if(*gauss_kon == 0)
    {
        s1 = randomDoubleWithSeed(buffer_seed);
        s2 = randomDoubleWithSeed(buffer_seed);
        temp = sqrt(-2 * log(s1));
        fi = 2 * M_PI * s2;
        value = temp * cos(fi);
        *gauss_y = temp * sin(fi);
        *gauss_kon = 1;
    }
    else
    {
        value = *gauss_y;
        *gauss_kon = 0;
    }

    return value;
}

double * speed_maxwell(double am, double te, double * values, struct drand48_data * buffer_seed, int * gauss_kon, double * gauss_y)
{	
    double coeff = sqrt(8.314 * te / am);

    values[0] = Gauss(gauss_kon, gauss_y, buffer_seed)*coeff;
    values[1] = Gauss(gauss_kon, gauss_y, buffer_seed)*coeff;
    values[2] = Gauss(gauss_kon, gauss_y, buffer_seed)*coeff;

    return values;
}

enum DEPORTATION_PLACE in_position(double x, double y, double z)
{
    enum DEPORTATION_PLACE result = AIR;

    double xy = sqrt(x*x + y*y);

    if ( xy >= r )
            result = WALL;
    else if ( z >= h )
            result = TOP;
    else if (z <= 0)
            result = BOTTOM;

    return result;
}

double distance(Particle * p1, Particle * p2)
{
    return sqrt ( (p1->x - p2->x) * (p1->x - p2->x) + 
                  (p1->y - p2->y) * (p1->y - p2->y) + 
                  (p1->z - p2->z) * (p1->z - p2->z)
                );
}

void initialize_particle(Particle * particle, int id, double tau, double r, struct drand48_data * buffer_seed)
{
    particle->id = id;
    particle->life_time = 0.0;
    particle->half_life = distrib(tau, buffer_seed);
    particle->status = INIT;

    particle->z = 1.e-2 * randomDoubleWithSeed(buffer_seed);
    double rr = r * sqrt(randomDoubleWithSeed(buffer_seed));
    double fir = 2 * M_PI * randomDoubleWithSeed(buffer_seed);
    particle->x = rr * cos(fir);
    particle->y = rr * sin(fir);
}

void set(Particle * p, double x, double y, double z)
{
    p->x = x;
    p->y = y;
    p->z = z;
}

void print_status(int nAir, int nWall, int nTop, int nBottom)
{
    printf("Bottom=%d Top=%d Wall=%d Air=%d\n", nBottom, nTop, nWall, nAir);
}

int checkDepositionPlace(Particle * particle, int nPathPoints, struct drand48_data * buffer_seed, int * count)
{
    enum DEPORTATION_PLACE pos_status;
    
    double impactDistance = 0.0;
    double fi0, theta0;
    double v_magnitude;
    double * v = (double * )calloc(3, sizeof(double));
    //double lifeTimeStepSum = 0.0;
    double xTemp, yTemp, zTemp;
    int gauss_kon = 1;
    double gauss_y = 0.0;
    
    double xMin = DBL_MAX;
    double xMax = -DBL_MAX;
    double yMin = DBL_MAX;
    double yMax = -DBL_MAX;
    double zMin = DBL_MAX;
    double zMax = -DBL_MAX; 
    
    int i;
    for (i = 0; i < nPathPoints; i++) 
    {
        (*count)++;
        
        impactDistance = distrib(half_path, buffer_seed);
        
        fi0 = 2 * M_PI * randomDoubleWithSeed(buffer_seed);
        theta0 = acos(1 - 2 * randomDoubleWithSeed(buffer_seed));
        
        xTemp = particle->x + impactDistance * sin(theta0) * cos(fi0);
	yTemp = particle->y + impactDistance * sin(theta0) * sin(fi0);
	zTemp = particle->z + impactDistance * cos(theta0);
        
        if(xTemp < xMin) xMin = xTemp;
	if(yTemp < yMin) yMin = yTemp;
	if(zTemp < zMin) zMin = zTemp;
				
	if(xTemp > xMax) xMax = xTemp;
	if(yTemp > yMax) yMax = yTemp;
	if(zTemp > zMax) zMax = zTemp;
        
        speed_maxwell(0.2, 293.0, v, buffer_seed, &gauss_kon, &gauss_y);
        v_magnitude = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        
        pos_status = in_position(xTemp,yTemp,zTemp);
        
        if(pos_status == AIR)
	{
            particle->life_time += (impactDistance / v_magnitude);
            
            if (particle->life_time >= particle->half_life)
            {
                particle->status = AIR;
                break;
            }
            else
            {
                particle->x = xTemp;
                particle->y = yTemp;
                particle->z = zTemp;
                continue;
            }	
	}
	else
	{
            particle->life_time += (impactDistance / v_magnitude);

            if (particle->life_time > particle->half_life)
            {
                particle->status = AIR;
            }
            else
            {
                if(pos_status == BOTTOM)
                {
                    particle->status = BOTTOM;
                }
                else if(pos_status == TOP)
                {
                    particle->status = TOP;
                }
                else
                {
                    particle->status = WALL;
                }
            }
            
            break;
	}
    }
    
    if( i == nPathPoints)
    {
        return 0;
    }
    else
        return 1;
}

void generatePartialTrajectory(Partial_trajectory * partial_trajectory, int * gauss_kon, double * gauss_y, int nPathPoints, struct drand48_data * buffer_seed)
{
    double impactDistance = 0.0;
    double fi0, theta0;
    double v_magnitude;
    double * v = (double * )calloc(3, sizeof(double));
    double lifeTimeStepSum = 0.0;
    
    struct drand48_data seedPathStart = *buffer_seed; 
    
    double xMin = DBL_MAX;
    double xMax = -DBL_MAX;
    double yMin = DBL_MAX;
    double yMax = -DBL_MAX;
    double zMin = DBL_MAX;
    double zMax = -DBL_MAX; 
    
    double xLast = .0; 
    double yLast = .0;
    double zLast = .0;
    
    int i;
    for (i = 0; i < nPathPoints; i++) 
    {
        impactDistance = distrib(half_path, buffer_seed);
        
        fi0 = 2 * M_PI * randomDoubleWithSeed(buffer_seed);
        theta0 = acos(1 - 2 * randomDoubleWithSeed(buffer_seed));
        
        xLast += impactDistance * sin(theta0) * cos(fi0);
	yLast += impactDistance * sin(theta0) * sin(fi0);
	zLast += impactDistance * cos(theta0);
        
	if(xLast < xMin) xMin = xLast;
	if(yLast < yMin) yMin = yLast;
	if(zLast < zMin) zMin = zLast;
				
	if(xLast > xMax) xMax = xLast;
	if(yLast > yMax) yMax = yLast;
	if(zLast > zMax) zMax = zLast;
        
        speed_maxwell(0.2, 293.0, v, buffer_seed, gauss_kon, gauss_y);
        
        v_magnitude = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        
        lifeTimeStepSum += (impactDistance / v_magnitude);
    }
    
    free(v);
    
    partial_trajectory->xMin = xMin;
    partial_trajectory->xMax = xMax;
    partial_trajectory->yMin = yMin;
    partial_trajectory->yMax = yMax;
    partial_trajectory->zMin = zMin;
    partial_trajectory->zMax = zMax;
    partial_trajectory->xLast = xLast;
    partial_trajectory->yLast = yLast;
    partial_trajectory->zLast = zLast;
    partial_trajectory->lifeTimeStepSum = lifeTimeStepSum;
    partial_trajectory->__a = seedPathStart.__a;
    partial_trajectory->__c = seedPathStart.__c;
    partial_trajectory->__init = seedPathStart.__init;
    partial_trajectory->__old_x1 = seedPathStart.__old_x[0];
    partial_trajectory->__old_x2 = seedPathStart.__old_x[1];
    partial_trajectory->__old_x3 = seedPathStart.__old_x[2];
    partial_trajectory->__x1 = seedPathStart.__x[0];
    partial_trajectory->__x2 = seedPathStart.__x[1];
    partial_trajectory->__x3 = seedPathStart.__x[2];
}

enum DEPORTATION_PLACE checkBoundingBox(Particle * particle, Partial_trajectory * partial_trajectory)
{
    enum DEPORTATION_PLACE p = AIR;
    
    double x = particle->x + partial_trajectory->xMax;
    double y = particle->y + partial_trajectory->yMin;
    double z = particle->z + partial_trajectory->zMin;

    if( (p = in_position(x, y, z)) != AIR)
    {
        return p;
    }

    x = particle->x + partial_trajectory->xMax;
    y = particle->y + partial_trajectory->yMin;
    z = particle->z + partial_trajectory->zMax;

    if( (p = in_position(x, y, z)) != AIR)
    {
        return p;
    }

    x = particle->x + partial_trajectory->xMin;
    y = particle->y + partial_trajectory->yMin;
    z = particle->z + partial_trajectory->zMax;

    if( (p = in_position(x, y, z)) != AIR)
    {
        return p;
    }

    x = particle->x + partial_trajectory->xMin;
    y = particle->y + partial_trajectory->yMin;
    z = particle->z + partial_trajectory->zMin;

    if( (p = in_position(x, y, z)) != AIR)
    {
        return p;
    }

   x = particle->x + partial_trajectory->xMax;
   y = particle->y + partial_trajectory->yMax;
   z = particle->z + partial_trajectory->zMin;

   if( (p = in_position(x, y, z)) != AIR)
   {
       return p;
   }

   x = particle->x + partial_trajectory->xMax;
   y = particle->y + partial_trajectory->yMax;
   z = particle->z + partial_trajectory->zMax;

   if( (p = in_position(x, y, z)) != AIR)
   {
       return p;
   }

   x = particle->x + partial_trajectory->xMin;
   y = particle->y + partial_trajectory->yMax;
   z = particle->z + partial_trajectory->zMax;

   if( (p = in_position(x, y, z)) != AIR)
   {
       return p;
   }

   x = particle->x + partial_trajectory->xMin;
   y = particle->y + partial_trajectory->yMax;
   z = particle->z + partial_trajectory->zMin;

   if( (p = in_position(x, y, z)) != AIR)
   {
       return p;
   }

   if((particle->life_time + partial_trajectory->lifeTimeStepSum) >= particle->half_life)
   {
       return DECOMPOSED;
   }
   
   return p;
}

void * waitForSignal(void * t)
{
    MPI_Status status;
    int flag = 0;
    
    while(flag == 0)
    {
        MPI_Iprobe(0, SIMULATION_END, MPI_COMM_WORLD, &flag, &status);
        sleep(1);
    }
    
    int * running = (int * ) t;
    *running = 0;
}

void putPointerAfterColon(char ** s_tmp)
{
    while(**s_tmp != ':')
        (*s_tmp)++;
    (*s_tmp)++;
}

/*
    Lightweight Conf.json parser
 */
int readJsonIntValue(char * s, char * jsonFieldName, int * value)
{
    char * s_tmp = strstr(s, jsonFieldName);
    
    if(s_tmp != NULL)
    {
        int i;
        for(i = 0; i < strlen(jsonFieldName); i++)
            s_tmp++;
        
        putPointerAfterColon(&s_tmp);

        sscanf(s_tmp, "%d", value);
        
        return 1;
    }
    else
        return -1;
}

/*
    Lightweight Conf.json parser
 */
Conf * readConf(char fileName[])
{
    FILE * f = fopen(fileName, "rb");
    
    Conf * conf = (Conf *) malloc(sizeof(Conf)); 
    
    char * buffer = NULL;
    long length;
    
    if (f)
    {
        fseek (f, 0, SEEK_END);
        length = ftell (f);
        fseek (f, 0, SEEK_SET);
        buffer = malloc (length);
        
        if (buffer)
        {
          fread (buffer, 1, length, f);
        }
    }
    
    readJsonIntValue(buffer, "nodeCount", &conf->nodeCount);
    readJsonIntValue(buffer, "processorsPerNode", &conf->processorsPerNode);
    readJsonIntValue(buffer, "producerCount", &conf->producerCount);
    readJsonIntValue(buffer, "consumerCount", &conf->consumerCount);
    readJsonIntValue(buffer, "partialTrajectoryLength", &conf->partialTrajectoryLength);
    readJsonIntValue(buffer, "particleCount", &conf->particleCount);
    
    fclose(f);
    
    return conf;
}

void createMpiTypeForParticle(MPI_Datatype  * dt_particle)
{
    int          blocklengths[3];
    MPI_Aint     offsets[3], intExtent, doubleExtent;
    MPI_Datatype types[3];
    
    offsets[0] = 0;
    types[0] = MPI_DOUBLE;
    blocklengths[0] = 5;
    
    MPI_Type_extent(MPI_INT, &doubleExtent);
    offsets[1] = 5 * doubleExtent;
    types[1] = MPI_INT;
    blocklengths[1] = 1;

    MPI_Type_extent(MPI_DOUBLE, &doubleExtent);
    offsets[2] = intExtent +  5 * doubleExtent;
    types[2] = MPI_SHORT_INT;
    blocklengths[2] = 1;

    MPI_Type_create_struct(3, blocklengths, offsets, types, dt_particle);
    MPI_Type_commit(dt_particle);  
}

void createMpiTypeForParticle(MPI_Datatype  * dt_particle)
{
    int          blocklengths[3];
    MPI_Aint     offsets[3], intExtent, doubleExtent;
    MPI_Datatype types[3];
    
    offsets[0] = 0;
    types[0] = MPI_DOUBLE;
    blocklengths[0] = 5;
    
    MPI_Type_extent(MPI_DOUBLE, &doubleExtent);
    offsets[1] = 5 * doubleExtent;
    types[1] = MPI_INT;
    blocklengths[1] = 1;

    MPI_Type_extent(MPI_INT, &intExtent);
    offsets[2] = intExtent +  5 * doubleExtent;
    types[2] = MPI_SHORT_INT;
    blocklengths[2] = 1;

    MPI_Type_create_struct(3, blocklengths, offsets, types, dt_particle);
    MPI_Type_commit(dt_particle);  
}

void initializeSeed(struct drand48_data * seed, Partial_trajectory * partialTrajectory)
{
    seed->__a = partialTrajectory->__a;
    seed->__c = partialTrajectory->__c;
    seed->__init = partialTrajectory->__init;
    seed->__old_x[0] = partialTrajectory->__old_x1;
    seed->__old_x[1] = partialTrajectory->__old_x2;
    seed->__old_x[2] = partialTrajectory->__old_x3;
    seed->__x[0] = partialTrajectory->__x1;
    seed->__x[1] = partialTrajectory->__x2;
    seed->__x[2] = partialTrajectory->__x3;
}

void printPartialTrajectory(Partial_trajectory * partial_trajectory)
{
    printf("lifeTimeStepSum = %lf \n", partial_trajectory->lifeTimeStepSum);
    printf("seed.__a = %llu \n", partial_trajectory->__a);
    printf("seed.__c = %hu \n", partial_trajectory->__c);
    printf("seed.__init = %hu \n", partial_trajectory->__init);
    
    printf("old_x = %hu \n", partial_trajectory->__old_x1);
    printf("old_x = %hu \n", partial_trajectory->__old_x2);
    printf("old_x = %hu \n", partial_trajectory->__old_x3);
    
    printf("__x = %hu \n", partial_trajectory->__x1);
    printf("__x = %hu \n", partial_trajectory->__x2);
    printf("__x = %hu \n", partial_trajectory->__x3);
}

void sendDepositionStatistics(int nAir, int nBottom, int nTop, int nWall)
{
     MPI_Send(&nAir, 1, MPI_INT, 0, AIR, MPI_COMM_WORLD); 
     MPI_Send(&nBottom, 1, MPI_INT, 0, BOTTOM, MPI_COMM_WORLD); 
     MPI_Send(&nTop, 1, MPI_INT, 0, TOP, MPI_COMM_WORLD); 
     MPI_Send(&nWall, 1, MPI_INT, 0, WALL, MPI_COMM_WORLD); 
}

void getDepositionStatistics(Conf * conf)
{
    int i;
    MPI_Status status;
    int nAir, nBottom, nTop, nWall;
    int nAirTmp, nBottomTmp, nTopTmp, nWallTmp;
    nAir = nBottom = nTop = nWall = 0;
    
    for(i = 1; i <= conf->consumerCount; i++)
    {   
        MPI_Recv(&nAirTmp, 1, MPI_INT, i, AIR, MPI_COMM_WORLD, &status);
        MPI_Recv(&nBottomTmp, 1, MPI_INT, i, BOTTOM, MPI_COMM_WORLD, &status);
        MPI_Recv(&nTopTmp, 1, MPI_INT, i, TOP, MPI_COMM_WORLD, &status);
        MPI_Recv(&nWallTmp, 1, MPI_INT, i, WALL, MPI_COMM_WORLD, &status);
        
        nAir += nAirTmp;
        nBottom += nBottomTmp;
        nTop += nTopTmp;
        nWall += nWallTmp;
    }
    printf("\n****** Statistika deponovanja ****** \n\n");
    printf(" nAir = %d,\n nBottom = %d,\n nTop = %d,\n nWall = %d\n", nAir, nBottom, nTop, nWall);
    printf("\n************************************ \n");
}

void writeTrajectoryStatistics(int producerCount, int processedTrajectoriesCount, int rank)
{
    int globalProducerCount;
    MPI_Reduce( &producerCount, &globalProducerCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    
    if(rank == 0)
    {
        printf("\n****** Statistika putanja  ********* \n\n");
        printf("   Broj proizvedenih putanja :  %d\n", globalProducerCount);
        printf("   Broj obradjenih   putanja :  %d\n", processedTrajectoriesCount);
        printf("\n************************************ \n\n");
    }
}

void * consumerListener(void * p);
void initialParticlePositionObserver(Conf * conf, int * running);
void runQueue(Conf * conf, int consumerCount, int producerCount);
void simulateParticles( int rank, Conf * conf);
void generatePathTrajektories(int rank,int partialTrajectoryLength);
void createMpiTypeForPartialTrajectory(MPI_Datatype  * dt_partial_trajectory);
void createMpiTypeForParticle(MPI_Datatype  * dt_particle);
