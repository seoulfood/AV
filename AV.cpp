#include <omp.h>
#include <cstdio>
#include <iostream>
#include <ostream>
#include <unistd.h>
#include <memory>
#include <stdlib.h>
#include <cmath>
#include <queue>
#include <set>
#include <vector>
#include <chrono>
#include "AV.h"
using std::cout;
using std::endl;
using namespace std::chrono;
using namespace AnisoVoro;

/*
BIG NOTE: ONLY THE SIM BOX EVER DEALS WITH ACTUAL SIM DIMENSIONS
    ALL OTHER CLASSES/STRUCTS AND SUBROUTINE DEAL WITH THE 
    VOXELIZATION BOX AND DIMENSIONS
*/

Queue::Queue(): capacity{100}, queueArray{new int[capacity]}{
    this->head = 0;
    this->tail = 0;
    this->sz = 0;
}
Queue::Queue(int c): capacity{c}, queueArray{new int[capacity]}{
    this->head = 0;
    this->tail = 0;
    this->sz = 0;
}
int Queue::front() {
    return queueArray[head];
}
void Queue::emptyQueue(){
    this->head = 0;
    this->tail = 0;
    this->sz = 0;
}
bool Queue::empty() {
    if(this->sz == 0){
        return true;
    }
    else{
        return false;
    }
}
bool Queue::pop() {
    if(empty()){
        return false;
    }
    if(head == capacity - 1){
        head = 0;
    }
    else{
        head += 1;
    }
    this->sz -= 1;
    return true;
}
bool Queue::push(int val) {
    if(this->sz == capacity) {
        cout << "Queue is full" << endl;
        return false;
    }
    else{
        queueArray[tail] = val;
        if (tail == capacity - 1){
            tail = 0;
        }
        else{
            tail += 1;
        }
        this->sz += 1;
        return true;
    }
}
int Queue::at(int i) {
    int modI = (head + i) % (this->capacity);
    return queueArray[modI];
    /*
    if(tail > head){
        return queueArray[head+i];
    }
    else{
        return queueArray[i - (capacity - 1 - head)];
    }
    */
}
int Queue::size(){
    /*
    if(tail >= head){
        return(tail-head);
    }
    else{
        return(capacity - (head-tail - 1));
    }
    */
    return this->sz;
}
void Queue::display(){
    //if(this->size() == 0){
    if(this->empty()){
        cout << "Queue is empty!";
    }
    else{
        for(int i = 0; i < this->size(); i++){
            cout << this->at(i) << " ";
        }
    }
    cout << "\tSize: " << this->size();
    cout << " head: " << this->head << " tail: " << this->tail;
    cout << endl;

    cout << "\tActual: ";
    for(int i = 0; i < this->capacity; i++){
        cout << this->queueArray[i] << " ";
    }
    cout << endl;
}

Shape::Shape() {
    this->points.reserve(1);
}
Shape::Shape(std::vector<Position>& points) {
    this->points = points; 
}
Shape::Shape(std::vector<int> &shapePoints, int xDim, int yDim, int zDim, int voxDegree){
    int size = xDim * yDim * zDim;
    this->points.reserve(1);
}
//Shape::shapeFromVertices(std::vector<Position>& v, int vD) {
    //this->vertices = v;
    //this->voxDegree = vD;
//}

Quaternion::Quaternion() {
    this->w = 0;
    this->x = 0;
    this->y = 0;
    this->z = 0;

}
Quaternion::Quaternion(double w, double x, double y, double z) {
    double normalize = sqrt(pow(w, 2) + pow(x, 2) + pow(y, 2) + pow(z, 2));
    this->w = w/normalize;
    this->x = x/normalize;
    this->y = y/normalize;
    this->z = z/normalize;
}

Position::Position() {
    x = 0;
    y = 0;
    z = 0;
}
Position::Position(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
}
Position rotatePoint(Position p, Quaternion q) {
    Position rp;


    double r00 = 2 * (q.w * q.w + q.x * q.x) - 1;
    double r01 = 2 * (q.x * q.y + q.w * q.z);
    double r02 = 2 * (q.x * q.z + q.w * q.y);


    double r10 = 2 * (q.x * q.y + q.w * q.z);
    double r11 = 2 * (q.w * q.w + q.y * q.y) - 1;
    double r12 = 2 * (q.y * q.z + q.w * q.x);


    double r20 = 2 * (q.x * q.z + q.w * q.y);
    double r21 = 2 * (q.y * q.z + q.w * q.x);
    double r22 = 2 * (q.w * q.w + q.z * q.z) - 1;

    rp.x = (p.x*r00) + (p.y*r01) + (p.z*r02);
    rp.y = (p.x*r10) + (p.y*r11) + (p.z*r12);
    rp.z = (p.x*r20) + (p.y*r21) + (p.z*r22);

    return rp;
};
Position rotatePoint(double pd[3], double qd[4]) {
    Position rp;

    Position p(pd[0], pd[1], pd[2]);
    Quaternion q(qd[0], qd[1], qd[2], qd[3]);


    double r00 = 2 * (q.w * q.w + q.x * q.x) - 1;
    double r01 = 2 * (q.x * q.y + q.w * q.z);
    double r02 = 2 * (q.x * q.z + q.w * q.y);


    double r10 = 2 * (q.x * q.y + q.w * q.z);
    double r11 = 2 * (q.w * q.w + q.y * q.y) - 1;
    double r12 = 2 * (q.y * q.z + q.w * q.x);


    double r20 = 2 * (q.x * q.z + q.w * q.y);
    double r21 = 2 * (q.y * q.z + q.w * q.x);
    double r22 = 2 * (q.w * q.w + q.z * q.z) - 1;

    rp.x = (p.x*r00) + (p.y*r01) + (p.z*r02);
    rp.y = (p.x*r10) + (p.y*r11) + (p.z*r12);
    rp.z = (p.x*r20) + (p.y*r21) + (p.z*r22);

    return rp;

}

BoxDim::BoxDim() {
    x = 0;
    y = 0;
    z = 0;
}
BoxDim::BoxDim(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
}
void BoxDim::print() {
    cout << "Box Dimensions: [" << x << ", " << y << ", " << z << "]" << endl;
}

VoxelIndex::VoxelIndex() {
    this->i = 0;
    this->x = 0;
    this->y = 0;
    this->z = 0;
}
VoxelIndex::VoxelIndex(Position vbp, BoxDim vbd) {
    this->x = vbp.x;
    this->y = vbp.y;
    this->z = vbp.z;
    this->i = (x) + (y*vbd.x) + (z * vbd.x * vbd.y);
}
VoxelIndex::VoxelIndex(Position vbp, BoxDim vbd, Position refCorner) {
    //TO DO
    //ref corner is the lowest x, y, and z
    this->x = vbp.x;
    this->y = vbp.y;
    this->z = vbp.z;

    double ix = x - refCorner.x;
    double iy = y - refCorner.y;
    double iz = z - refCorner.z;
    
    this->i = (ix) + (iy*vbd.x) + (iz * vbd.x * vbd.y);
}
VoxelIndex::VoxelIndex(int i, BoxDim vbd) {
    this->i = i;
    this->z = floor(i/(vbd.x*vbd.y));
    this->y = floor((i - (z*vbd.x*vbd.y))/vbd.x);
    this->x = i - (z*vbd.x*vbd.y) - (y*vbd.x);
}
VoxelIndex::VoxelIndex(int i, BoxDim vbd, Position refCorner) {
    //TO DO
    this->i = i;

    double iz = floor(i/(vbd.x*vbd.y));
    double iy = floor((i - (z*vbd.x*vbd.y))/vbd.x);
    double ix = i - (z*vbd.x*vbd.y) - (y*vbd.x);

    this->x = ix + refCorner.x;
    this->y = iy + refCorner.y;
    this->z = iz + refCorner.z;

}


VoxelBit::VoxelBit() {
    layer = -1;
    index = VoxelIndex();
    isParticle = false;
    isBoundary = false;
}
VoxelBit::VoxelBit(int i, bool isParticle, int particleNum, BoxDim vbd) {
    index = VoxelIndex(i, vbd);
    this->isParticle = isParticle;
    if (isParticle) {
        layer = 0;
        origins.insert(particleNum);
    }
    else {
        layer = -1;
    }
    isBoundary = false;
}
VoxelBit::VoxelBit(int i, bool isParticle, int particleNum, BoxDim vbd, Position refCorner) {
    //TO DO
    index = VoxelIndex(i, vbd);
    this->isParticle = isParticle;
    if (isParticle) {
        layer = 0;
        origins.insert(particleNum);
    }
    else {
        layer = -1;
    }
    isBoundary = false;
}
void VoxelBit::getNeighborsIndices2D(BoxDim vbd, int (&neighbors)[8]) {
    int vBoxX = vbd.x;
    int vBoxY = vbd.y;
    int i = index.i;

    int xpos = (index.x != vBoxX - 1) ? (1) : (-vBoxX + 1);
    int xneg = (index.x != 0) ? (-1) : (vBoxX - 1);
    int ypos = (index.y != vBoxY - 1) ? (vBoxX) : (-vBoxX * (vBoxY - 1));
    int yneg = (index.y != 0) ? (-vBoxX) : (vBoxX * (vBoxY - 1));

    int n[8] = {i + xpos, i + xneg,
                        i + ypos, i + ypos + xneg, i + ypos + xpos, 
                        i + yneg, i + yneg + xneg, i + yneg + xpos 
                        };

    for(int a = 0; a < 8; a++) {
        neighbors[a] = n[a];
    }

}
void VoxelBit::getNeighborsIndices3D(BoxDim vbd, int (&neighbors)[26]) {
    int vBoxX = vbd.x;
    int vBoxY = vbd.y;
    int vBoxZ = vbd.z;
    int i = index.i;

    int xpos = (index.x != vBoxX - 1) ? (1) : (-vBoxX + 1);
    int xneg = (index.x != 0) ? (-1) : (vBoxX - 1);
    int ypos = (index.y != vBoxY - 1) ? (vBoxX) : (-vBoxX * (vBoxY - 1));
    int yneg = (index.y != 0) ? (-vBoxX) : (vBoxX * (vBoxY - 1));
    int zpos = (index.z != vBoxZ - 1) ? (vBoxX*vBoxY) : (-(vBoxX*vBoxY*(vBoxZ - 1)));
    int zneg = (index.z != 0) ? (-vBoxX*vBoxY) : (vBoxX*vBoxY*(vBoxZ - 1));
    int izpos = i + zpos;
    int izneg = i + zneg;
    
    int n[26] = {i + xpos, i + xneg,
                        i + ypos, i + ypos + xneg, i + ypos + xpos, 
                        i + yneg, i + yneg + xneg, i + yneg + xpos, 
                        izpos, izpos + xpos, izpos + xneg,
                        izpos + ypos, izpos + ypos + xneg, izpos + ypos + xpos, 
                        izpos + yneg, izpos + yneg + xneg, izpos + yneg + xpos, 
                        izneg, izneg + xpos, izneg + xneg,
                        izneg + ypos, izneg + ypos + xneg, izneg + ypos + xpos, 
                        izneg + yneg, izneg + yneg + xneg, izneg + yneg + xpos 
                        };

    for(int a = 0; a < 26; a++) {
        neighbors[a] = n[a];
    }


}

VoxelVector::VoxelVector() {
    this->v.reserve(0);
}
VoxelVector::VoxelVector(int size) {
    this->v.reserve(size);
}
VoxelVector::VoxelVector(int xLength, int yLength, int zLength, int voxDegree) {
    int voxArrSize = (voxDegree*xLength) * (voxDegree*yLength);
    voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize* (voxDegree*zLength));
    this->v.reserve(voxArrSize);
}

SimBox::SimBox() {
    simBoxDim = BoxDim();
    voxBoxDim = BoxDim();
    is2D = false;
    this->pVoxArr = VoxelVector().v;   
    this->voxDegree = 1;
    this->center = Position(0, 0, 0);
}

SimBox::SimBox(double xLength, double yLength, double zLength, std::vector<VoxelBit>& pVoxArr, int voxDegree) {
    simBoxDim = BoxDim(xLength, yLength, zLength);
    voxBoxDim = BoxDim(floor(xLength*voxDegree), 
              floor(yLength*voxDegree),
              floor(zLength*voxDegree));
    is2D = (zLength == 0) ? true : false;
    this->pVoxArr = pVoxArr;
    this->voxDegree = voxDegree;
    this->center = Position(0, 0, 0);
    initialize();
}

SimBox::SimBox(double xLength, double yLength, double zLength, int voxDegree) {
    simBoxDim = BoxDim(xLength, yLength, zLength);
    voxBoxDim = BoxDim(floor(xLength*voxDegree), 
              floor(yLength*voxDegree),
              floor(zLength*voxDegree));
    is2D = (zLength == 0) ? true : false;
    int voxArrSize = (voxDegree*xLength) * (voxDegree*yLength);
    voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize* (voxDegree*zLength));
    std::vector<VoxelBit> tempVec(voxArrSize);
    this->pVoxArr = tempVec;
    this->voxDegree = voxDegree;
    this->center = Position(0, 0, 0);
    this->voxCenter = Position(center.x * voxDegree,
                               center.y * voxDegree, 
                               center.z * voxDegree);
    setDevice(0);
    initialize();
}

SimBox::SimBox(double xLength, double yLength, double zLength, int voxDegree, int mode) {
    simBoxDim = BoxDim(xLength, yLength, zLength);
    voxBoxDim = BoxDim(floor(xLength*voxDegree), 
              floor(yLength*voxDegree),
              floor(zLength*voxDegree));
    is2D = (zLength == 0) ? true : false;
    int voxArrSize = (voxDegree*xLength) * (voxDegree*yLength);
    voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize* (voxDegree*zLength));
    std::vector<VoxelBit> tempVec(voxArrSize);
    this->pVoxArr = tempVec;
    this->voxDegree = voxDegree;
    this->center = Position(0, 0, 0);
    this->voxCenter = Position(center.x * voxDegree,
                               center.y * voxDegree, 
                               center.z * voxDegree);
    setDevice(mode);
    initialize();
}

SimBox::SimBox(double xLength, double yLength, double zLength, int voxDegree, int mode, int rank, int mpiWorldSize) {
    simBoxDim = BoxDim(xLength, yLength, zLength);
    voxBoxDim = BoxDim(floor(xLength*voxDegree), 
              floor(yLength*voxDegree),
              floor(zLength*voxDegree));
    is2D = (zLength == 0) ? true : false;
    int voxArrSize = (voxDegree*xLength) * (voxDegree*yLength);
    voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize* (voxDegree*zLength));
    std::vector<VoxelBit> tempVec(voxArrSize);
    this->pVoxArr = tempVec;
    this->voxDegree = voxDegree;
    this->center = Position(0, 0, 0);
    this->voxCenter = Position(center.x * voxDegree,
                               center.y * voxDegree, 
                               center.z * voxDegree);
    setDevice(mode);
    initialize();
}


void SimBox::setPVoxelArraySize(double xLength, double yLength, double zLength){
    if(mode == 0){
        int voxArrSize = (voxDegree*xLength) * (voxDegree*yLength);
        voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize* (voxDegree*zLength));
        std::vector<VoxelBit> tempVec(voxArrSize);
        this->pVoxArr = tempVec;
    }
    else if(mode == 1){
        int voxArrSize = (voxDegree*xLength) * (voxDegree*yLength);
        voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize* (voxDegree*zLength));
        std::vector<VoxelBit> tempVec(voxArrSize);
        this->pVoxArr = tempVec;
    }
}

void SimBox::setVoxel(Position p, bool isParticle, int particleNum = -1) {
    //adjustPosition(p);
    int i = indexFromPosition(p);
    pVoxArr.at(i) = VoxelBit(i, isParticle, particleNum, voxBoxDim);
}

void SimBox::adjustPosition(Position &p){
    p.x = ((this->simBoxDim.x)/2.0) + (this->voxCenter.x);
    p.y = ((this->simBoxDim.y)/2.0) + (this->voxCenter.y);
    p.z = ((this->simBoxDim.z)/2.0) + (this->voxCenter.z);
}

void SimBox::placeShape(Shape s, Quaternion q, Position p, int particleNum = -1) {
    Position rp;
    Position sp;

    double r00 = pow(q.w, 2) + pow(q.x, 2) - pow(q.y, 2) - pow(q.z, 2); 
    double r01 = (2*q.x*q.y) - (2*q.w*q.z);
    double r02 = (2*q.x*q.z) + (2*q.w*q.y);
    double r10 = (2*q.x*q.y) + (2*q.w*q.z);
    double r11 = pow(q.w, 2) - pow(q.x, 2) + pow(q.y, 2) - pow(q.z, 2);
    double r12 = (2*q.y*q.z) - (2*q.w*q.x);
    double r20 = (2*q.x*q.z) - (2*q.w*q.y);
    double r21 = (2*q.y*q.z) + (2*q.w*q.y);
    double r22 = pow(q.w, 2) - pow(q.x, 2) - pow(q.y, 2) + pow(q.z, 2);

    for (int i = 0; i < s.points.size(); i++) {
        sp = s.points.at(i);
        rp.x = ((sp.x*r00) + (sp.y*r01) + (sp.z*r02)) + p.x;
        rp.y = ((sp.x*r10) + (sp.y*r11) + (sp.z*r12)) + p.y;
        rp.z = ((sp.x*r20) + (sp.y*r21) + (sp.z*r22)) + p.z;
        setVoxel(rp, true, particleNum);
    }
}

void SimBox::printBox(){
    if(mode == 0)
    {
        for (int i = 0; i < static_cast<int>(pVoxArr.size()); i++) {
            Position p = positionFromIndex(i);
            VoxelBit v  = pVoxArr.at(i);
            cout << "[" << p.x << "][" << p.y << "][" << p.z << "]" << " is " << v.layer << endl;  
        }
    }
    else{
        cout << "Cannot return full cell information for non serial methods. Memory limits make this unnecessary.";
    }
}
void SimBox::printBoundaries(){
    if(mode == 0)
    {
        for (int i = 0; i < static_cast<int>(pVoxArr.size()); i++) {
            Position p = positionFromIndex(i);
            VoxelBit v  = pVoxArr.at(i);
            cout << "Boundary[" << p.x << "][" << p.y << "][" << p.z << "]" << " is " << v.isBoundary << endl;  
        }
    }
    else{
        cout << "Cannot return full cell information for non serial methods. Memory limits make this unnecessary.";
    }

}
void SimBox::printCells(){
    if(mode == 0)
    {
        for (int i = 0; i < static_cast<int>(pVoxArr.size()); i++) {
            Position p = positionFromIndex(i);
            VoxelBit v  = pVoxArr.at(i);
            if(v.isParticle) {
                cout << "ParticleCell[" << p.x << "][" << p.y << "][" << p.z << "]" << " is " << *v.origins.begin() << endl;  
            }
            else if(v.isBoundary) {
                cout << "ParticleCell[" << p.x << "][" << p.y << "][" << p.z << "]" << " is " << "B" << endl;  
            }
            else{
                cout << "ParticleCell[" << p.x << "][" << p.y << "][" << p.z << "]" << " is " << "V" << endl;
            }
        }
    }
    else{
        cout << "Cannot return full cell information for non serial methods. Memory limits make this unnecessary.";
    }
}
void SimBox::particleNum(int num) {
    this->partNum = num;
}
int SimBox::particleNum() {
    return partNum;
}

void SimBox::setDevice(int mode){
    this->mode = mode;
    this->rank = 0;
    this->mpiWorldSize = 1;
}
void SimBox::setDevice(int mode, int rank, int mpiWorldSize){
    this->mode = mode;
    this->rank = rank;
    this->mpiWorldSize = mpiWorldSize;
}

void SimBox::runVoro(){
    /*
    double start;
    double stop;
    if(this->rank == 0){
        start = MPI_Wtime();
    }
    */
    if(this->mode == 0)//Using plain Serial
    {
        initializeQueue();
        runLayerByLayer();
    }
    else if(this->mode == 1)//Using Dynamic MPI
    {
        initializeQueue();
        runLayerByLayerMPI();
    }
    else if(this->mode == 2)//OpenMP GPU
    {
        initializeQueue();
        //runLayerByLayerGPU(); still in progress
    }
    /*
    if(this->rank == 0){
        stop = MPI_Wtime();
        double duration = stop - start;
        cout << "Time taken to run voronoi: " << duration 
             << " seconds" << endl;
    }*/
}

void SimBox::divideSimBox(){
    int xLength = this->voxBoxDim.x;
    int yLength = this->voxBoxDim.y;
    int zLength = this->voxBoxDim.z;

    int totalDomains = sqrt(2);

}

void SimBox::initialize() {
    int x, y, z;
    int i;
    int voxX = voxBoxDim.x;
    int voxY = voxBoxDim.y;
    int zMax = (voxBoxDim.z == 0) ? (1) : voxBoxDim.z;

    if(this->mode == 0)
    {
        cout << "Initializing using serial cpu" << endl;
        for(x = 0; x < voxBoxDim.x; x++) {
            for(y = 0; y < voxBoxDim.y; y++) {
                zMax = (voxBoxDim.z == 0) ? (1) : voxBoxDim.z;
                //In the future if it's 2D vbd.z should be 1, not 0
                for(z = 0; z < zMax; z++) {
                    i = ((x) + (y*voxBoxDim.x) + (z * voxBoxDim.x * voxBoxDim.y));
                    pVoxArr.at(i) = VoxelBit(i, false, 0, voxBoxDim);
                }
            }
        }
    }
    else if(this->mode == 1){
        cout << "Initializing using MPI on cpu" << endl;
        divideSimBox();
    }
    else if(this->mode == 2)
    {
        cout << "Using omp to initialize!" << endl;
        int id, np, a;
        #pragma omp parallel for num_threads(32)
        for(x = 0; x < voxX; x++) {
            id = omp_get_thread_num();
            for(y = 0; y < voxY; y++) {
                for(z = 0; z < zMax; z++) {
                    i = ((x) + (y*voxBoxDim.x) + (z * voxBoxDim.x * voxBoxDim.y));
                    pVoxArr.at(i) = VoxelBit(i, false, 0, voxBoxDim);
                }
            }
        }
    }

}

int SimBox::indexFromPosition(Position p) {
    int i = voxDegree * ((p.x) + (p.y*voxBoxDim.x) + (p.z * voxBoxDim.x * voxBoxDim.y));
    return i;
}

Position SimBox::positionFromIndex(int i) {
    double z = floor(i/(voxBoxDim.x*voxBoxDim.y));
    double y = floor((i - (z*voxBoxDim.x*voxBoxDim.y))/voxBoxDim.x);
    double x = i - (z*voxBoxDim.x*voxBoxDim.y) - (y*voxBoxDim.x);
    Position p(x, y, z);
    return (p);
}

void SimBox::initializeQueue() {
    VoxelBit v;
    for (int i = 0; i < static_cast<int>(pVoxArr.size()); i++) {
        v = pVoxArr.at(i);
        if (v.layer == 0) {
            updateNeighbors(0, v);
        }
    }
    updateOrigins(1);
}

void SimBox::runLayerByLayerGPU() {
    int currentLayer = 1;
    VoxelBit v;
    int i;

    for(int i = 0; i < layerRun.size();i++){

    }
    //layerRun.emptyQueue();
    while(!layerRun.empty()) {
        i = layerRun.front();
        layerRun.pop(); 
        v = pVoxArr.at(i);
        if(v.layer != currentLayer) {
            currentLayer += 1;
            updateOrigins(currentLayer);
        }
        updateNeighbors(currentLayer, v);
    }
    updateOrigins(currentLayer);
}

void SimBox::runLayerByLayer() {
    int currentLayer = 1;
    VoxelBit v;
    int i;

    while(!layerRun.empty()) {
        i = layerRun.front();
        layerRun.pop(); 
        v = pVoxArr.at(i);
        if(v.layer != currentLayer) {
            currentLayer += 1;
            updateOrigins(currentLayer);
        }
        updateNeighbors(currentLayer, v);
    }
    updateOrigins(currentLayer);
}

void SimBox::runLayerByLayerMPI() {
    int currentLayer = 1;
    VoxelBit v;
    int i;

    while(!layerRun.empty()) {
        i = layerRun.front();
        layerRun.pop(); 
        v = pVoxArr.at(i);
        if(v.layer != currentLayer) {
            currentLayer += 1;
            updateOrigins(currentLayer);
        }
        updateNeighbors(currentLayer, v);
    }
    updateOrigins(currentLayer);
}

void SimBox::updateNeighbors(int currentLayer, VoxelBit& v) {
    VoxelBit nv;
    if (is2D) {
        int neighbors[8];
        v.getNeighborsIndices2D(voxBoxDim, neighbors);
        for(int n = 0; n < 8; n++) {
            nv = pVoxArr.at(neighbors[n]);
            if(nv.layer == -1) {
                nv.layer = currentLayer + 1;
                if(!v.isBoundary){
                    nv.origins.insert(v.origins.begin(), v.origins.end());
                }
                layerRun.push(nv.index.i);
                originRun.push(nv.index.i);
                pVoxArr.at(neighbors[n]) = nv;
            }
            else if(nv.layer == currentLayer + 1 && !v.isBoundary){
                nv.origins.insert(v.origins.begin(), v.origins.end());
                pVoxArr.at(neighbors[n]) = nv;
            }
        }
    }
    else {//3D Situation
        int neighbors[26];
        v.getNeighborsIndices3D(voxBoxDim, neighbors);

        for(int n = 0; n < 8; n++) {
            nv = pVoxArr.at(neighbors[n]);
            if(nv.layer == -1) {
                nv.layer = currentLayer + 1;
                if(!v.isBoundary){
                    nv.origins.insert(v.origins.begin(), v.origins.end());
                }
                layerRun.push(nv.index.i);
                originRun.push(nv.index.i);
                pVoxArr.at(neighbors[n]) = nv;
            }
            else if(nv.layer == currentLayer + 1 && !v.isBoundary){
                nv.origins.insert(v.origins.begin(), v.origins.end());
                pVoxArr.at(neighbors[n]) = nv;
            }
        }
    }
}

void SimBox::updateOrigins(int currentLayer) {
    VoxelBit w;
    int o;
    while(!originRun.empty()){
        o = originRun.front();
        originRun.pop();
        w = pVoxArr.at(o);
        originUpdater(currentLayer, w);
    }
}

void SimBox::originUpdater(int currentLayer, VoxelBit& v) {
    VoxelBit nv;
    int voidCount = 0;
    bool override = false;
    if (is2D) {
        int num = 8;
        int neighbors[8];
        v.getNeighborsIndices2D(voxBoxDim, neighbors);
        for(int n = 0; n < num; n++) {
            nv = pVoxArr.at(neighbors[n]);
            if(nv.layer == -1) {
            voidCount++;
            }
            if(nv.layer == v.layer) {
                if(nv.origins != v.origins && nv.origins.size() == 1){
                    override = true;
                }
            }
        }
        if ((voidCount == 0 || v.origins.size() > 1) || (override)){
            v.isBoundary = true;
            pVoxArr.at(v.index.i) = v;
        }
    }
    else{
        int num = 26;
        int neighbors[26];
        v.getNeighborsIndices3D(voxBoxDim, neighbors);
        for(int n = 0; n < num; n++) {
            nv = pVoxArr.at(neighbors[n]);
            if(nv.layer == -1) {
            voidCount++;
            }
            if(nv.layer == v.layer) {
                if(nv.origins != v.origins && nv.origins.size() == 1){
                    override = true;
                }
            }
        }
        if ((voidCount == 0 || v.origins.size() > 1) || (override)){
            v.isBoundary = true;
            pVoxArr.at(v.index.i) = v;
        }
    }
}
