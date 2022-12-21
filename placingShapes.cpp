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
using std::cout;
using std::endl;
using namespace AnisoVoro;
using namespace std::chrono;

int main(int argc, char** argv) {
    double xLength = 250;
    double yLength = 250;
    double zLength = 0;
    if(argc < 4){
        std::cout << "ERROR CODE 1:" << std::endl;
        std::cout << "\tPlease enter more arguments in ./run <printBool> <voxDegree>" << std::endl;
        exit(1);
    }

    double voxDegree = atoi(argv[2]);
    int voxArrSize = (voxDegree*xLength) * (voxDegree*yLength);
    voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize* (voxDegree*zLength));
    std::vector<VoxelBit> voxArr(voxArrSize);
    //pentagon area = 0.25 * sqrt(5*(5 + (2*sqrt(5))))pow(s,2)
    //where s is the length of one of the sides
    //std::vector<Position> unitPentagon(2*pow(voxDegree,2));
    auto start = high_resolution_clock::now();
    //SimBox sim = SimBox(xLength, yLength, zLength, voxArr, voxDegree);
    SimBox sim = SimBox(xLength, yLength, zLength, voxDegree);

    int shapeSize = 0;
    for (double x = 0; x < 1; x = x + (1/voxDegree)){
        for(double y = -2; y < 3; y = y + (1/voxDegree)){
            shapeSize = shapeSize + 1;
        }
    }
    std::vector<Position> unitSquareArr(shapeSize);
    cout << "Shape size is " << shapeSize << endl;
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

    //angle = (1.0/(rand() % 100)+15) * M_PI;
    angle = M_PI * atoi(argv[3]);
    
    cout << "Angle: " << angle << endl;
    q = Quaternion(cos(angle/2.0), 0, 0, sin(angle/2.0));
    sim.placeShape(square, q, Position(5, 8, 0), 1);

    angle = M_PI/2;
    cout << "Angle: " << angle << endl;
    q = Quaternion(cos(angle/2.0), 0, 0, sin(angle/2.0));
    sim.placeShape(square, q, Position(13, 17, 0), 2);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken to initialize voronoi: " << duration.count() 
         << " microseconds" << endl;
    //start = high_resolution_clock::now();
    sim.runVoro();
    //stop = high_resolution_clock::now();
    //duration = duration_cast<microseconds>(stop - start);
    //cout << "Time taken to run voronoi: " << duration.count() 
    //     << " microseconds" << endl;
    int print = atoi(argv[1]);
    if(print == 1){
        sim.printBox();
        sim.printBoundaries();
        sim.printCells();
    }
   
    return 0;
}
