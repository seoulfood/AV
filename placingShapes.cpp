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
#include "AV.h"
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

    double xLength = 1000;
    double yLength = 1000;
    double zLength = 0;
    if(argc < 3){
        std::cout << "ERROR CODE 1:" << std::endl;
        std::cout << "\tPlease enter more arguments in ./run <printBool> <voxDegree>" << std::endl;
        exit(1);
    }

    double voxDegree = atoi(argv[2]);
    //int voxArrSize = (voxDegree*xLength) * (voxDegree*yLength);
    //voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize* (voxDegree*zLength));
    //std::vector<VoxelBit> voxArr(voxArrSize);
    //pentagon area = 0.25 * sqrt(5*(5 + (2*sqrt(5))))pow(s,2)
    //where s is the length of one of the sides
    //std::vector<Position> unitPentagon(2*pow(voxDegree,2));
    //SimBox sim = SimBox(xLength, yLength, zLength, voxArr, voxDegree);
    
    SimBox sim;
    sim = SimBox(xLength, yLength, zLength, voxDegree, 1, rank, P);

    int shapeSize = 0;
    for (double x = 0; x < 5; x = x + (1/voxDegree)){
        for(double y = 0; y < 3; y = y + (1/voxDegree)){
            shapeSize = shapeSize + 1;
        }
    }
    std::vector<Position> shapeArr(shapeSize);
    //cout << "Shape size is " << shapeSize << endl;
    shapeSize = 0;
    for (double x = 0; x < 5; x = x + (1/voxDegree)){
        for(double y = 0; y < 3; y = y + (1/voxDegree)){
            shapeArr.at(shapeSize) = Position(x, y, 0);
            shapeSize = shapeSize + 1;
        }
    }
    Shape tallRec(shapeArr);

    shapeSize = 0;
    for (double x = 0; x < 2; x = x + (1/voxDegree)){
        for(double y = 0; y < 10; y = y + (1/voxDegree)){
            shapeSize = shapeSize + 1;
        }
    }
    shapeArr.resize(shapeSize);
    //cout << "Shape size is " << shapeSize << endl;
    shapeSize = 0;
    for (double x = 0; x < 2; x = x + (1/voxDegree)){
        for(double y = 0; y < 10; y = y + (1/voxDegree)){
            shapeArr.at(shapeSize) = Position(x, y, 0);
            shapeSize = shapeSize + 1;
        }
    }
    Shape longRec(shapeArr);
 
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

//    sim.placeShape(tallRec, q, Position(1, 17, 0), 0);
//    sim.placeShape(tallRec, q, Position(5, 28, 0), 1);
//    sim.placeShape(tallRec, q, Position(28, 17, 0), 2);
//    sim.placeShape(tallRec, q, Position(21, 30, 0), 3);


    int particleNum = 0;
    for(int x = 0; x < xLength; x+=60){
        for(int y = 0; y < yLength; y+=50){
            sim.placeShape(tallRec, q, Position(x, y, 0), particleNum);
            particleNum++;
        }
    }

    for(int x = 20; x < xLength; x+=60){
        for(int y = 5; y < yLength; y+=40){
            sim.placeShape(longRec, q, Position(x, y, 0), particleNum);
            particleNum++;
        }
    }
 
 
    /*
    sim.placeShape(square, q, Position(1, 7, 0), 0);
    sim.placeShape(square, q, Position(30, 9, 0), 1);
    sim.placeShape(square, q, Position(14, 30, 0), 2);
    sim.placeShape(square, q, Position(28, 30, 0), 3);
    */


    if(rank == 0){
        stopPlacingShapes = MPI_Wtime();
        durationPlacingShapes = stopPlacingShapes - startPlacingShapes;
        cout << "Time taken to place shapes: " << durationPlacingShapes << " microseconds" << endl;
        cout << "Ranks: " << P << endl;
        cout << "voxDegree: " << voxDegree << endl;
    }
    sim.runVoro();

    //cout << "<<<<<<<<<<<<<<<<<Rank " << rank << " has finished runVoro." << endl;

    //MPI_Barrier(MPI_COMM_WORLD);
    //stop = high_resolution_clock::now();
    //duration = duration_cast<microseconds>(stop - start);
    //cout << "Time taken to run voronoi: " << duration.count() 
    //     << " microseconds" << endl;
    int print = atoi(argv[1]);
    if(print == 1){
        sim.printBox();
        sim.printBoundaries();
        sim.printCells();
        sim.printVoxRank();
    }

    //sim.printVoxRank();
    MPI_Finalize();
    return 0;
}
