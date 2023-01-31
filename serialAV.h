#include <iostream>
#include <unistd.h>
#include <memory>
#include <stdlib.h>
#include <cmath>
#include <queue>
#include <set>
#include <vector>
using std::cout;
using std::endl;
using namespace std::chrono;

#ifndef AV_H
#define AV_H

namespace AnisoVoro {

    class Queue {
        public:
            Queue();
            Queue(int initialSize);
            int front();
            //bool full();
            bool empty();
            bool pop();
            bool push(int val);
            int at(int i);
            int size();
            void display();
            void emptyQueue();

        private:
            int sz;
            int head;
            int tail;
            int capacity;
            bool full;
            std::unique_ptr<int[]> queueArray;

    };

    class Position {
        public:
            double x, y, z;
            Position();
            Position(double x, double y, double z);
            //friend std::ostream& operator<<(std::ostream &os, const Position &p);
            //Position& operator<<(const Position& p);
            //~Position();
    };

    //Shapes are centered at the origin by default 
    class Shape {
        public:
            std::vector<Position> points;
            //int voxDegree;
            Shape();
            Shape(std::vector<Position>& points);
            Shape(std::vector<int> &shapePoints, int xDim, int yDim, int zDim, int voxDegree);
            //shapeFromVertices(std::vector<Position>& v, int vD);
    };

    class Quaternion {
        public:
            double w;
            double x;
            double y;
            double z;
            Quaternion();
            Quaternion(double w, double x, double y, double z);
    };
   
    Position rotatePoint(Position p, Quaternion q);
    Position rotatePoint(double pd[3], double qd[4]);

    class BoxDim {
        public:
            double x, y, z;
            BoxDim();
            BoxDim(double x, double y, double z);
            void print();
    };
    
    class VoxelIndex{
        public:
            int i, x, y, z;
            VoxelIndex();
            VoxelIndex(Position vbp, BoxDim vbd);
            VoxelIndex(Position vbp, BoxDim vbd, Position refCorner);
            VoxelIndex(int i, BoxDim vbd);
            VoxelIndex(int i, BoxDim vbd, Position refCorner);
            //~VoxelIndex();
    };
    
    class VoxelBit {
        public:
            int layer;
            VoxelIndex index;
            std::set<int> origins;
            bool isParticle, isBoundary;
            //~VoxelBit();
            VoxelBit();
            VoxelBit(int i, bool isParticle, int particleNum, BoxDim vbd); 
            void getNeighborsIndices2D(BoxDim vbd, int (&neighbors)[8]);
            void getNeighborsIndices3D(BoxDim vbd, int (&neighbors)[26]);
    };
    
    class VoxelVector {
        public:
            std::vector<VoxelBit> v;
            VoxelVector();
            VoxelVector(int size);
            VoxelVector(int xLength, int yLength, int zLength, int voxDegree);
    };
    
    class SimBox {
        public:
            std::vector<VoxelBit> pVoxArr; // initial particle voxelization
            std::vector<Shape> shapeArr;
            std::vector<Position> positions;
            int voxDegree; //number of voxels per simulation box unit
            BoxDim simBoxDim;
            BoxDim voxBoxDim; //simUnitBoxDim * voxDegree
            bool is2D;
            int partNum;
            Position center;
            Position voxCenter;
            int mode; 

            SimBox();    
            SimBox(double xLength, double yLength, double zLength, 
                   std::vector<VoxelBit>& pVoxArr, int voxDegree);
            SimBox(double xLength, double yLength, double zLength, int voxDegree);
            SimBox(double xLength, double yLength, double zLength, int voxDegree, int mode);
            void setVoxel(Position p, bool isParticle, int particleNum);
            void placeShape(Shape s, Quaternion q, Position p, int particleNum);
            //void particleTypes(std::vector<int>& type_id, 
            //                   std::vector<std::vector<Position>>& particleTypes,
            //                   std::vector<Position> positions,
            //                   std::vector<Quaternion>& orientations);
            //std::vector<double> localDensity();
            //std::vector<std::vector<int> > densityNeighbors();
            int particleNum();
            void particleNum(int num);
            void printBox();
            void printBoundaries();
            void printCells();
            void setDevice(int mode);
            void setDevice(int mode, int rank);
            void runVoro();
    
        //~SimBox();
    
        private:
            int rank; //MPI Rank
            int mpiWorldSize;//MPI world size
            std::vector<int> mpiNeighbors;
            std::vector<int> voxelTracker;

            std::vector<VoxelBit> neighbors;
            std::queue<int> layerRun;
            std::queue<int> originRun;
            std::queue<int> boundaryIndices;

            void setPVoxelArraySize(double xLength, double yLength, double zLength);
            void initialize();
            void adjustPosition(Position &p);
            int indexFromPosition(Position p); 
            Position positionFromIndex(int i);
            void initializeQueue();
            void runLayerByLayer();
            void divideSimBox();
            void initializeQueueMPI();
            void runLayerByLayerMPI();
            void runLayerByLayerGPU();
            void updateNeighbors(int currentLayer, VoxelBit& v);
            void updateOrigins(int currentLayer);
            void originUpdater(int currentLayer, VoxelBit& v);
    };
    
}

#endif
