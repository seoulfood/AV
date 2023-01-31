#include <mpi.h>
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include <cmath>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include "serialAV.h"
#include <omp.h>
using std::cout;
using std::endl;
using namespace AnisoVoro;
using namespace std::chrono;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int P;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double xLength = 25;
    double yLength = 25;
    double zLength = 0;
    if(argc < 4){
        std::cout << "ERROR CODE 1:" << std::endl;
        std::cout << "\tPlease enter more arguments in ./run <printBool> <voxDegree> <bool useGPU>" << std::endl;
        exit(1);
    }

    double voxDegree = atoi(argv[2]);
    int voxArrSize = (voxDegree*xLength) * (voxDegree*yLength);
    voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize* (voxDegree*zLength));
    std::vector<VoxelBit> voxArr(voxArrSize);
    //pentagon area = 0.25 * sqrt(5*(5 + (2*sqrt(5))))pow(s,2)
    //where s is the length of one of the sides
    //std::vector<Position> unitPentagon(2*pow(voxDegree,2));
    double start;
    if(rank == 0){
        start = MPI_Wtime();
    }
    //SimBox sim = SimBox(xLength, yLength, zLength, voxArr, voxDegree);
    SimBox sim;
    if(atoi(argv[3]) == 0){
        cout << "Setting to use serial" << endl;
        sim = SimBox(xLength, yLength, zLength, voxDegree);
    }
    else if(atoi(argv[3]) == 1){
        cout << "Setting to use MPI" << endl;
        sim = SimBox(xLength, yLength, zLength, voxDegree, 1, P, rank);
    }
    else{
        cout << "Setting to use GPU" << endl;
        sim = SimBox(xLength, yLength, zLength, voxDegree, 2);
    }

    int shapeSize = 0;
    for (double x = 0; x < 1; x = x + (1/voxDegree)){
        for(double y = -2; y < 3; y = y + (1/voxDegree)){
            shapeSize = shapeSize + 1;
        }
    }
    std::vector<Position> unitSquareArr(shapeSize);
    //cout << "Shape size is " << shapeSize << endl;
    shapeSize = 0;
    for (double x = 0; x < 1; x = x + (1/voxDegree)){
        for(double y = -2; y < 3; y = y + (1/voxDegree)){
            unitSquareArr.at(shapeSize) = Position(x, y, 0);
            shapeSize = shapeSize + 1;
        }
    }
    Shape square(unitSquareArr);
    double angle;
    Quaternion q;

    srand(time(NULL));

    angle = 0;

    int shapeNumber = 0;
    double startPlacingShapes;
    double stopPlacingShapes;
    double durationPlacingShapes;
    double stop;
    double duration;
    if(rank == 0){
        startPlacingShapes = MPI_Wtime();
        
    }
    sim.placeShape(square, q, Position(4, 2, 0), 0);
    sim.placeShape(square, q, Position(11, 9, 0), 1);

    if(rank == 0){
        stopPlacingShapes = MPI_Wtime();
        durationPlacingShapes = stopPlacingShapes - startPlacingShapes;
        cout << "Time taken to place shapes: " << durationPlacingShapes << " microseconds" << endl;
        stop = MPI_Wtime();
        //duration = duration_cast<microseconds>(stop - start);
        duration = stop - start;
        cout << "Time taken to initialize voronoi: " << duration 
             << " microseconds" << endl;
    }
    //start = high_resolution_clock::now();
    sim.runVoro();
    //stop = high_resolution_clock::now();
    //duration = duration_cast<microseconds>(stop - start);
    //cout << "Time taken to run voronoi: " << duration.count() 
    //     << " microseconds" << endl;
    int print = atoi(argv[1]);
    if(rank == 0){
        if(print == 1){
            sim.printBox();
            sim.printBoundaries();
            sim.printCells();
        }
    }
    return 0;
}
