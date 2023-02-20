#include <mpi.h>
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
using namespace AnisoVoro;
//using namespace std::chrono;

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
    index = VoxelIndex(i, vbd, refCorner);
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
    this->refCorner = Position(0, 0, 0);
}

SimBox::SimBox(double xLength, double yLength, double zLength, std::vector<VoxelBit>& pVoxArr, int voxDegree) {
    simBoxDim = BoxDim(xLength, yLength, zLength);
    voxBoxDim = BoxDim(floor(xLength*voxDegree), 
              floor(yLength*voxDegree),
              floor(zLength*voxDegree));
    is2D = (zLength == 0) ? true : false;
    this->pVoxArr = pVoxArr;
    this->voxDegree = voxDegree;
    this->refCorner = Position(0, 0, 0);
    initialize();
}

SimBox::SimBox(double xLength, double yLength, double zLength, int voxDegree) {
    setDevice(0);
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
    this->refCorner = Position(0, 0, 0);
    this->voxRefCorner = Position(refCorner.x * voxDegree,
                               refCorner.y * voxDegree, 
                               refCorner.z * voxDegree);
    initialize();
}

SimBox::SimBox(double xLength, double yLength, double zLength, int voxDegree, int mode) {
    setDevice(mode);
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
    this->refCorner = Position(0, 0, 0);
    this->voxRefCorner = Position(refCorner.x * voxDegree,
                               refCorner.y * voxDegree, 
                               refCorner.z * voxDegree);
    initialize();
}

SimBox::SimBox(double xLength, double yLength, double zLength, int voxDegree, int mode, int rank, int mpiWorldSize) {
    this->is2D = (zLength == 0) ? true : false;
    this->voxDegree = voxDegree;
    setDevice(mode, rank, mpiWorldSize);
    cout << (xLength * voxDegree) << ", " << (yLength * voxDegree) << endl;
    if(is2D){
        dcomp.divideSimBox2D(floor(xLength * voxDegree), 
                             floor(yLength * voxDegree));
    }
    else{
        cout << "Not implemented yet!" << endl;
    }
    xLength = (this->dcomp.localXMax) - (this->dcomp.localXMin);
    yLength = (this->dcomp.localYMax) - (this->dcomp.localYMin);
    zLength = (this->dcomp.localZMax) - (this->dcomp.localZMin);
    this->voxBoxDim = BoxDim(xLength, 
                             yLength,
                             zLength);
    voxBoxDim.print();
    this->simBoxDim = BoxDim(xLength/voxDegree, yLength/voxDegree, zLength/voxDegree);
    int voxArrSize = (xLength) * (yLength);
    voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize*zLength);
    std::vector<VoxelBit> tempVec(voxArrSize);
    this->pVoxArr = tempVec;
    this->voxRefCorner = Position(this->dcomp.localXMin,
                               this->dcomp.localYMin, 
                               this->dcomp.localZMin);
    this->refCorner = Position(voxRefCorner.x / voxDegree, 
                               voxRefCorner.y / voxDegree,
                               voxRefCorner.z / voxDegree);
    cout << "Local ref corner is [" << voxRefCorner.x << "][" << voxRefCorner.y << "][" << voxRefCorner.z << "]" << endl;
    initialize();
}


void SimBox::setVoxel(Position p, bool isParticle, int particleNum = -1) {
    adjustInputPosition(p);
    int i = indexFromPosition(p);
    pVoxArr.at(i) = VoxelBit(i, isParticle, particleNum, voxBoxDim, this->voxRefCorner);
}

void SimBox::adjustInputPosition(Position &p){
    p.x = p.x - voxRefCorner.x;
    p.y = p.y - voxRefCorner.y;
    p.z = p.z - voxRefCorner.z;
}

void SimBox::adjustOutputPosition(Position &p){
    p.x = p.x + voxRefCorner.x;
    p.y = p.y + voxRefCorner.y;
    p.z = p.z + voxRefCorner.z;
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
        if(insideVoxBox(rp)){
            setVoxel(rp, true, particleNum);
        }
    }
}

bool SimBox::insideVoxBox(Position p){
    if(p.x >= this->refCorner.x && p.x <= (this->refCorner.x + voxBoxDim.x) &&
       p.y >= this->refCorner.y && p.y <= (this->refCorner.y + voxBoxDim.y) &&
       p.z >= this->refCorner.z && p.z <= (this->refCorner.z + voxBoxDim.z))
    {
        return true;
    }
    else{
        return false;
    }

}

void SimBox::printBox(){
    for (int i = 0; i < static_cast<int>(pVoxArr.size()); i++) {
        Position p = positionFromIndex(i);
        VoxelBit v  = pVoxArr.at(i);
        cout << "[" << p.x << "][" << p.y << "][" << p.z << "]" << " is " << v.layer << endl;  
    }
}
void SimBox::printBoundaries(){
    for (int i = 0; i < static_cast<int>(pVoxArr.size()); i++) {
        Position p = positionFromIndex(i);
        VoxelBit v  = pVoxArr.at(i);
        cout << "Boundary[" << p.x << "][" << p.y << "][" << p.z << "]" << " is " << v.isBoundary << endl;  
    }
}
void SimBox::printCells(){
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
void SimBox::printVoxRank(){
    for (int i = 0; i < static_cast<int>(pVoxArr.size()); i++) {
        Position p = positionFromIndex(i);
        VoxelBit v  = pVoxArr.at(i);
        cout << "VoxCell[" << p.x << "][" << p.y << "][" << p.z << "]" << " is " << dcomp.rank << endl;  
    }
}

void SimBox::setReferenceCorner(Position p){
    this->refCorner.x = this->refCorner.x + this->voxDegree * p.x;
    this->refCorner.y = this->refCorner.y + this->voxDegree * p.y;
    this->refCorner.z = this->refCorner.z + this->voxDegree * p.z;
}

void SimBox::particleNum(int num) {
    this->partNum = num;
}
int SimBox::particleNum() {
    return partNum;
}

void SimBox::setDevice(int mode){
    this->mode = mode;
    this->dcomp = DomainDecomposition(0, 1);
}
void SimBox::setDevice(int mode, int rank, int mpiWorldSize){
    this->mode = mode;
    this->dcomp = DomainDecomposition(rank, mpiWorldSize);
}

void SimBox::runVoro(){
    double start;
    double stop;
    if(this->dcomp.rank == 0){
        start = MPI_Wtime();
    }
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
    if(this->dcomp.rank == 0){
        stop = MPI_Wtime();
        double duration = stop - start;
        cout << "Time taken to run voronoi: " << duration 
             << " seconds" << endl;
    }
}


void SimBox::initialize() {
    int x, y, z;
    int i;
    int voxX = voxBoxDim.x;
    int voxY = voxBoxDim.y;
    int zMax = (voxBoxDim.z == 0) ? (1) : voxBoxDim.z;

    if(this->mode == 0 || this->mode == 1)
    {
        for(x = 0; x < voxBoxDim.x; x++) {
            for(y = 0; y < voxBoxDim.y; y++) {
                zMax = (voxBoxDim.z == 0) ? (1) : voxBoxDim.z;
                //In the future if it's 2D vbd.z should be 1, not 0
                for(z = 0; z < zMax; z++) {
                    i = ((x) + (y*voxBoxDim.x) + (z * voxBoxDim.x * voxBoxDim.y));
                    pVoxArr.at(i) = VoxelBit(i, false, 0, voxBoxDim, this->voxRefCorner);
                }
            }
        }
    }
    else if(this->mode == 2)
    {
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


DomainDecomposition::DomainDecomposition(){
    this->rank = 1;
    this->P = 1;
}
DomainDecomposition::DomainDecomposition(int rank, int P){
    this->rank = rank;
    this->P = P;
}

void DomainDecomposition::divideSimBox2D(int xLength, int yLength){
    this->length = floor(sqrt(this->P));
    this->mainCount = std::pow(length, 2);
    this->mainDomainNumber = this->rank % mainCount;
    this->localNumber = this->rank / mainCount; 


    int xSpace = xLength / length; 
    int ySpace = yLength / length; 
    int zSpace = 0; 

    int xMod = mainDomainNumber % length;
    int yMod = mainDomainNumber / length;

    this->numberInLocal = (mainDomainNumber >= P % mainCount) ? (P/mainCount) : ((P/mainCount) + 1);
    int localxMod = xSpace / numberInLocal;

    //cout << "xSpace is " << xSpace << "\tlocalxMod is " << localxMod << endl;

    cout << "Rank " << this->rank << "\tMain Domain: " << this->mainDomainNumber << "\tlocal rank: " << this->localNumber << "\tlocal number total "<< numberInLocal << endl;

    cout << "Rank " << this->rank << "XMod: " << xMod << endl;

    if(numberInLocal == 1){
        this->localXMin = xSpace * xMod;
        this->localXMax = (xMod == xLength - 1) ? (xLength): (localXMin + xSpace - 1);
        cout << "Rank: " << this->rank << "\tXmin: " << this->localXMin << "\tXMax: " << this->localXMax << endl;
    }
    else{
        cout << "There are subspaces managed by rank " << this->rank << endl;
        this->localXMin = (xSpace * xMod) + (localNumber * localxMod);
        if(localNumber + 1 == numberInLocal){
            this->localXMax = (xMod == xLength - 1) ? (xLength): ((xSpace*xMod) + xSpace - 1);
        }
        else{
            this->localXMax = (localXMin + (localxMod) - 1);
        }
    }

    this->localYMin = ySpace * yMod;
    this->localYMax = (yMod == length - 1) ? (yLength): (localYMin + ySpace - 1);

    this->localZMin = 0;
    this->localZMax = 0;
}

/*
void DomainDecomposition::divideSimBox3D(int xLength, int yLength, int zLength){
    cout << "Not yet implemented!" << endl;
}

int DomainDecomposition::PlusXRank(){
    if(this->localNumber + 1 == this->numberInLocal){
        return mainDomainNumber+1;
    }
    else{
        return mainCount * (localNumber + 1) + mainDomainNumber
    }
}
int DomainDecomposition::PlusY(){

}
int DomainDecomposition::PlusZ(){

}
int DomainDecomposition::MinusX(){
    if(this->localNumber == 1){
        if{}
        return (mainDomainNumber-1);
    }
    else{
        return mainCount * (localNumber - 1) + mainDomainNumber
    }

}
int DomainDecomposition::MinusY(){

}
int DomainDecomposition::MinusZ(){

}

void DomainDecomposition::sendData(&int data, int recvRank){

}

void DomainDecomposition::ReceiveData(&int data, int sendRank){

}
*/