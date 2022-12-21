#include <iostream>
#include <memory>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include <cmath>
#include <unistd.h>
#include "serialAV.h"
using std::cout;
using std::endl;
using namespace AnisoVoro;
using namespace std::chrono;

int main(int argc, char** argv) {
    double xLength = 20;
    double yLength = 20;
    double zLength = 0;
    double voxDegree = 10;
    int voxArrSize = (voxDegree*xLength) * (voxDegree*yLength);
    voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize* (voxDegree*zLength));
    std::vector<VoxelBit> voxArr(voxArrSize);
    auto start = high_resolution_clock::now();
    //SimBox sim = SimBox(xLength, yLength, zLength, voxArr, voxDegree);
    SimBox sim = SimBox(xLength, yLength, zLength, voxDegree);

    for (double x = 4; x < 6; x = x + (1/voxDegree)) {
        for (double y = 2; y < 4; y = y + (1/voxDegree)) {
            Position p(x, y, 0);
            sim.setVoxel(p, true, 0);
        }
    }
    for (double x = 4; x < 6; x = x + (1/voxDegree)) {
        for (double y = 7; y < 9; y = y + (1/voxDegree)) {
            Position p(x, y, 0);
            sim.setVoxel(p, true, 1);
        }
    }

    for (double x = 14; x < 19; x = x + (1/voxDegree)) {
        for (double y = 7; y < 10; y = y + (1/voxDegree)) {
            Position p(x, y, 0);
            sim.setVoxel(p, true, 2);
        }
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken to initialize voronoi: " << duration.count() 
         << " microseconds" << endl;
    sim.runVoro();
    int print = atoi(argv[1]);
    if(print == 1){
        sim.printBox();
        sim.printBoundaries();
        sim.printCells();
    }
   
    return 0;
}
