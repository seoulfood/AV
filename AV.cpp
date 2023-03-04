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
VoxelIndex::VoxelIndex(Position vbp, BoxDim vbd, Position voxRefCorner) {
    this->x = vbp.x;
    this->y = vbp.y;
    this->z = vbp.z;

    double ix = x - voxRefCorner.x;
    double iy = y - voxRefCorner.y;
    double iz = z - voxRefCorner.z;
    
    this->i = (ix) + (iy*vbd.x) + (iz * vbd.x * vbd.y);
}
VoxelIndex::VoxelIndex(int i, BoxDim vbd) {
    this->i = i;
    this->z = floor(i/(vbd.x*vbd.y));
    this->y = floor((i - (z*vbd.x*vbd.y))/vbd.x);
    this->x = i - (z*vbd.x*vbd.y) - (y*vbd.x);
}
VoxelIndex::VoxelIndex(int i, BoxDim vbd, Position voxRefCorner) {
    this->i = i;
 
    double iz = floor(i/(vbd.x*vbd.y));
    double iy = floor((i - (iz*vbd.x*vbd.y))/vbd.x);
    double ix = i - (iz*vbd.x*vbd.y) - (iy*vbd.x);

    this->x = ix + voxRefCorner.x;
    this->y = iy + voxRefCorner.y;
    this->z = iz + voxRefCorner.z;

    /*
    if(ix == 0 && iy == 0 && iz == 0){
        cout << "\t--------->p started as [" << ix << "][" << iy << "] and ended as [" << this->x << "][" << this->y << "]" << endl;
    }
    cout << "VoxCell[" << this->x << "][" << this->y << "][" << '0' << "] is 0" << endl;
    */
}

VoxelBit::VoxelBit() {
    this->index = VoxelIndex();
}
VoxelBit::VoxelBit(int i, bool isParticle, int particleNum, BoxDim vbd) {
    this->index = VoxelIndex(i, vbd);
    this->isParticle = isParticle;
    if (isParticle) {
        this->layer = 0;
        origins.insert(particleNum);
    }
    else {
        this->layer = -1;
    }
}
VoxelBit::VoxelBit(int i, bool isParticle, int particleNum, BoxDim vbd, Position voxRefCorner) {
    this->index = VoxelIndex(i, vbd, voxRefCorner);
    this->isParticle = isParticle;
    if (isParticle) {
        this->layer = 0;
        origins.insert(particleNum);
    }
    else {
        this->layer = -1;
    }
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
    totalBoxDim = BoxDim();
    voxBoxDim = BoxDim();
    is2D = false;
    this->pVoxArr = VoxelVector().v;   
    this->voxDegree = 1;
    this->mainRefCorner = Position(0, 0, 0);
}

/*
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
*/

SimBox::SimBox(double xLength, double yLength, double zLength, int voxDegree, int mode, int rank, int mpiWorldSize) {

    //converting everything to vox units
    xLength = floor(xLength * voxDegree);
    yLength = floor(yLength * voxDegree);
    zLength = floor(zLength * voxDegree);
    
    this->totalBoxDim = BoxDim(xLength, yLength, zLength);

    //unsure why we need this?
    //Eventually can make the box moveable but add that later
    this->mainRefCorner = Position(0, 0, 0);

    this->is2D = (zLength == 0) ? true : false;
    this->voxDegree = voxDegree;
    setDevice(mode, rank, mpiWorldSize);
    if(is2D){
        dcomp.divideSimBox2D(xLength, 
                             yLength);

        for(int x = dcomp.localXMin; x <= dcomp.localXMax; x++)
        {
            for(int y = dcomp.localYMin; y <= dcomp.localYMax; y++)
            {
                //cout << "VoxCell[" << x << "][" << y << "][" << "0" << "]" << " is " << dcomp.rank << endl;  
                int a = 5;
            }
        }
    }
    else{
        cout << "Not implemented yet!" << endl;
    }

    xLength = (dcomp.localXMax != xLength) ? 1 : 0;
    yLength = (dcomp.localYMax != yLength) ? 1 : 0;
    zLength = (dcomp.localZMax != zLength) ? 1 : 0;

    xLength += (this->dcomp.localXMax) - (this->dcomp.localXMin);
    yLength += (this->dcomp.localYMax) - (this->dcomp.localYMin);
    zLength += (is2D) ? 0 : (this->dcomp.localZMax) - (this->dcomp.localZMin);

    if(mpiWorldSize > 1){
        adjustForGhostCells(&xLength, &yLength, &zLength);
    }

    this->voxBoxDim = BoxDim(xLength, 
                             yLength,
                             zLength);
    voxBoxDim.print();
    int voxArrSize = (xLength) * (yLength);
    voxArrSize = (zLength==0) ? (voxArrSize) : (voxArrSize*zLength);
    std::vector<VoxelBit> tempVec(voxArrSize);
    this->pVoxArr = tempVec;
    this->voxRefCorner = Position(this->dcomp.localXMin,
                               this->dcomp.localYMin, 
                               this->dcomp.localZMin);
    if(mpiWorldSize > 1){
        adjustVoxRefCorner();
    }

    //cout << "Local vox ref corner is [" << voxRefCorner.x << "][" << voxRefCorner.y << "][" << voxRefCorner.z << "]" << endl;

    initialize();
}


void SimBox::setVoxel(Position p, bool isParticle, int particleNum = -1) {
    Position a = p;
    adjustInputPosition(p);
    int i = indexFromPosition(p);
    pVoxArr.at(i) = VoxelBit(i, isParticle, particleNum, voxBoxDim, this->voxRefCorner);
}


void SimBox::adjustGhostCells(){
    VoxelBit v;
    VoxelIndex a;
    int pVoxSize = static_cast<int>(pVoxArr.size());
    for(int i = 0; i < pVoxSize; i += voxBoxDim.x){
        (pVoxArr.at(i)).isGhost = true;
        if(this->dcomp.localXMin == mainRefCorner.x){
            ((pVoxArr.at(i)).index).x = totalBoxDim.x;
        }
    }
    for(int i = voxBoxDim.x-1; i < pVoxSize; i += voxBoxDim.x){
        (pVoxArr.at(i)).isGhost = true;
        if(this->dcomp.localXMax + 1 == totalBoxDim.x){
            ((pVoxArr.at(i)).index).x = mainRefCorner.x;
        }
    }
    for(int i = 0; i < voxBoxDim.x; i++){
        (pVoxArr.at(i)).isGhost = true;
        if(this->dcomp.localYMin == mainRefCorner.y){
            ((pVoxArr.at(i)).index).y = totalBoxDim.y;
        }
    }
    for(int i = pVoxSize - voxBoxDim.x; i < pVoxSize ; i++){
        (pVoxArr.at(i)).isGhost = true;
        if(this->dcomp.localYMax + 1 == totalBoxDim.y){
            ((pVoxArr.at(i)).index).y = mainRefCorner.y;
        }
    }
    this->ghostsToSend[0].reserve(voxBoxDim.x);
    this->ghostsToSend[1].reserve(voxBoxDim.x);
    this->ghostsToSend[2].reserve(voxBoxDim.y);
    this->ghostsToSend[3].reserve(voxBoxDim.y);
    this->ghostsToSend[4].reserve(voxBoxDim.z);
    this->ghostsToSend[5].reserve(voxBoxDim.z);

    this->ghostsToRecv[0].reserve(voxBoxDim.x);
    this->ghostsToRecv[1].reserve(voxBoxDim.x);
    this->ghostsToRecv[2].reserve(voxBoxDim.y);
    this->ghostsToRecv[3].reserve(voxBoxDim.y);
    this->ghostsToRecv[4].reserve(voxBoxDim.z);
    this->ghostsToRecv[5].reserve(voxBoxDim.z);
}

void SimBox::adjustForGhostCells(double *xLength, double *yLength, double *zLength){
    *xLength += 2;
    *yLength += 2;
    if(!this->is2D){*zLength += 2;}
}

void SimBox::adjustVoxRefCorner(){
    this->voxRefCorner.x -= 1;
    this->voxRefCorner.y -= 1;
    if(!this->is2D){this->voxRefCorner.z -= 1;}
}

void SimBox::adjustInputPosition(Position &p){
    p.x -= voxRefCorner.x;
    p.y -= voxRefCorner.y;
    p.z -= voxRefCorner.z;
}

void SimBox::adjustOutputPosition(Position &p){
    p.x += voxRefCorner.x;
    p.y += voxRefCorner.y;
    p.z += voxRefCorner.z;
}

void SimBox::placeShape(Shape s, Quaternion q, Position p, int particleNum = -1) {
    Position rp;
    Position sp;

    /*
    double r00 = pow(q.w, 2) + pow(q.x, 2) - pow(q.y, 2) - pow(q.z, 2); 
    double r01 = (2*q.x*q.y) - (2*q.w*q.z);
    double r02 = (2*q.x*q.z) + (2*q.w*q.y);
    double r10 = (2*q.x*q.y) + (2*q.w*q.z);
    double r11 = pow(q.w, 2) - pow(q.x, 2) + pow(q.y, 2) - pow(q.z, 2);
    double r12 = (2*q.y*q.z) - (2*q.w*q.x);
    double r20 = (2*q.x*q.z) - (2*q.w*q.y);
    double r21 = (2*q.y*q.z) + (2*q.w*q.y);
    double r22 = pow(q.w, 2) - pow(q.x, 2) - pow(q.y, 2) + pow(q.z, 2);
    */

    for (int i = 0; i < s.points.size(); i++) {
        sp = s.points.at(i);
        /*
        rp.x = ((sp.x*r00) + (sp.y*r01) + (sp.z*r02)) + p.x;
        rp.y = ((sp.x*r10) + (sp.y*r11) + (sp.z*r12)) + p.y;
        rp.z = ((sp.x*r20) + (sp.y*r21) + (sp.z*r22)) + p.z;
        */
        rp.x = sp.x + p.x;
        rp.y = sp.y + p.y;
        rp.z = sp.z + p.z;
        if(insideVoxBox(rp)){
            setVoxel(rp, true, particleNum);
        }
    }
}

bool SimBox::insideVoxBox(Position p){
    //NOTE: FixThis!!!!!
    /*
    cout << "Is [" << p.x << "][" << p.y << "] inside box [" 
    << (voxRefCorner.x + 1) << "][" << (voxRefCorner.y + 1) << "], [" 
    << (voxRefCorner.x + voxBoxDim.x - 1) << "][" << (voxRefCorner.y + 1) << "], [" 
    << (voxRefCorner.x + 1) << "][" << (voxRefCorner.y + voxBoxDim.y - 1) << "], ["
    << (voxRefCorner.x + voxBoxDim.x - 1) << "][" << (voxRefCorner.y + voxBoxDim.y - 1) << "]?";
    */
    if(p.x >= this->voxRefCorner.x + 1 && p.x <= (this->voxRefCorner.x + voxBoxDim.x - 1) &&
       p.y >= this->voxRefCorner.y + 1 && p.y <= (this->voxRefCorner.y + voxBoxDim.y - 1))
       // && p.z >= this->voxRefCorner.z && p.z <= (this->voxRefCorner.z + voxBoxDim.z))
    {
        //cout << "\t ---> Yes! " << "\tRank: " << this->dcomp.rank << endl;
        return true;
    }
    else{
        //cout << "\t ---> No! " << "\tRank: " << this->dcomp.rank << endl;
        return false;
    }

}

void SimBox::printBox(){
    int pVoxSize = static_cast<int>(pVoxArr.size());
    if(this->dcomp.P == 1){
        for (int i = 0; i < pVoxSize; i++) {
            VoxelBit v  = pVoxArr.at(i);
            VoxelIndex a = v.index;
            cout << "[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << v.layer << endl;  
        }
    }
    else{
        bool g1;
        bool g2;
        bool g3;
        bool g4;

        for (int i = 0; i < pVoxSize; i++) {
            VoxelBit v  = pVoxArr.at(i);
            VoxelIndex a = v.index;
            g1 = (i % int(voxBoxDim.x) == 0);
            g2 = ((i+1) % int(voxBoxDim.x) == 0);
            g3 = (i <= int(voxBoxDim.x));
            g4 = ((pVoxSize - int(voxBoxDim.x)) <= i);
            if(g1 || g2 || g3 || g4)
            {
                continue;
            }
            cout << "[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << v.layer << endl;  
        }
    }
}
void SimBox::printBoundaries(){
    int pVoxSize = static_cast<int>(pVoxArr.size());
    if(this->dcomp.P == 1){
        for (int i = 0; i < pVoxSize; i++) {
            VoxelBit v  = pVoxArr.at(i);
            VoxelIndex a = v.index;
            cout << "Boundary[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << v.isBoundary << endl;  
        }
    }
    else{
        bool g1;
        bool g2;
        bool g3;
        bool g4;

        for (int i = 0; i < pVoxSize; i++) {
            VoxelBit v  = pVoxArr.at(i);
            VoxelIndex a = v.index;
            g1 = (i % int(voxBoxDim.x) == 0);
            g2 = ((i+1) % int(voxBoxDim.x) == 0);
            g3 = (i <= int(voxBoxDim.x));
            g4 = ((pVoxSize - int(voxBoxDim.x)) <= i);
            if(g1 || g2 || g3 || g4)
            {
                continue;
            }
            cout << "Boundary[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << v.isBoundary << endl;  
        }
    }

}
void SimBox::printCells(){
    int pVoxSize = static_cast<int>(pVoxArr.size());
    if(this->dcomp.P == 1){
        for (int i = 0; i < pVoxSize; i++) {
            VoxelBit v  = pVoxArr.at(i);
            VoxelIndex a = v.index;
    
            if(v.isParticle) {
                cout << "ParticleCell[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << *v.origins.begin() << endl;  
            }
            else if(v.isBoundary) {
                cout << "ParticleCell[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << "B" << endl;  
            }
            else{
                cout << "ParticleCell[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << "V" << endl;
            }
        }
    }
    else{
        bool g1;
        bool g2;
        bool g3;
        bool g4;

        for (int i = 0; i < pVoxSize; i++) {
            VoxelBit v  = pVoxArr.at(i);
            VoxelIndex a = v.index;
            g1 = (i % int(voxBoxDim.x) == 0);
            g2 = ((i+1) % int(voxBoxDim.x) == 0);
            g3 = (i <= int(voxBoxDim.x));
            g4 = ((pVoxSize - int(voxBoxDim.x)) <= i);
            if(g1 || g2 || g3 || g4)
            {
                continue;
            }

            if(v.isParticle) {
                cout << "ParticleCell[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << *v.origins.begin() << endl;  
            }
            else if(v.isBoundary) {
                cout << "ParticleCell[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << "B" << endl;  
            }
            else{
                cout << "ParticleCell[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << "V" << endl;
            }
        }
    }
}
void SimBox::printVoxRank(){
    int pVoxSize = static_cast<int>(pVoxArr.size());
    if(this->dcomp.P == 1){
        for (int i = 0; i < pVoxSize; i++) {
            VoxelBit v  = pVoxArr.at(i);
            VoxelIndex a = v.index;
            cout << "VoxCell[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << dcomp.rank << endl;  
        }
    }
    else{
        bool g1;
        bool g2;
        bool g3;
        bool g4;

        for (int i = 0; i < pVoxSize; i++) {
            VoxelBit v  = pVoxArr.at(i);
            VoxelIndex a = v.index;

            g1 = (i % int(voxBoxDim.x) == 0);
            g2 = ((i+1) % int(voxBoxDim.x) == 0);
            g3 = (i <= int(voxBoxDim.x));
            g4 = ((pVoxSize - int(voxBoxDim.x)) <= i);

            if(g1 || g2 || g3 || g4)
            {
                //cout << "VoxCell[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << "G" << endl;  
                continue;
            }
            if(!v.isParticle){
                cout << "VoxCell[" << a.x << "][" << a.y << "][" << a.z << "]" << " is " << dcomp.rank << endl;  
            }
            else{
                cout << "VoxCell[" << a.x << "][" << a.y << "][" << a.z << "]" << " is P" << endl;  
            }
        }
    }
}

void SimBox::setReferenceCorner(Position p){
    this->mainRefCorner.x = this->mainRefCorner.x + this->voxDegree * p.x;
    this->mainRefCorner.y = this->mainRefCorner.y + this->voxDegree * p.y;
    this->mainRefCorner.z = this->mainRefCorner.z + this->voxDegree * p.z;
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
    if(this->mode == 0 || this->dcomp.P == 1)//Using plain Serial
    {
        initializeQueue();
        runLayerByLayer();
    }
    else if(this->dcomp.P > 1)//Using Dynamic MPI
    {
        initializeQueue();
        //runLayerByLayerMPI();
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

    if(this->dcomp.P > 1){
        adjustGhostCells();
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
    adjustOutputPosition(p);
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

    if(this->dcomp.P > 1){
        updateGhosts(1);
        int a = 1;
    }
}

void SimBox::runLayerByLayerGPU() {
    cout << "NOT IMPLEMENTED!!" << endl;
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
            updateGhosts(currentLayer);
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

        for(int n = 0; n < 26; n++) {
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

void SimBox::updateGhosts(int layer){

    //In here we need to send, receive all ghosts

    sendRecvGhosts2D(layer);
    //int i = 1;

    //And then we need to update the origins if needed
    //Then add the ghosts into the layerRun for the next step with the same rules as the regular cells

}

void SimBox::integrateGhosts(std::vector<int> *ghosts){
    cout << "Suffer!!!" << endl;
}

//HORRENDOUS LONG SUBROUTINE
void SimBox::sendRecvGhosts2D(int layer){
    int pVoxSize = static_cast<int>(pVoxArr.size());
    VoxelBit v;
    int m;
    int n;
    int s;
    int e; 
    int w; 

    int northSendSize;
    int southSendSize;
    int eastSendSize;
    int westSendSize;

    int northRecvSize;
    int southRecvSize;
    int eastRecvSize;
    int westRecvSize;

    //North vector creation
    n = 0;
    for(int i = pVoxSize - ((voxBoxDim.x*2)-2); i < pVoxSize - (voxBoxDim.x+1); i++){
        v = pVoxArr.at(i);
        if(v.layer == layer){
            ghostsToSend[0].push_back(v.index.x);
            ghostsToSend[0].push_back(v.index.y);
            ghostsToSend[0].push_back(v.index.z);

            if(v.isParticle || v.isBoundary){
                m = v.isParticle ? 1 : 2;
                ghostsToSend[0].push_back(m);
            }
            else{
                ghostsToSend[0].push_back(0);
            }

            for(auto& ori : v.origins){
                ghostsToSend[0].push_back(ori);
            }
            ghostsToSend[0].push_back(-1);
            n++;
        }
    }
    northSendSize = static_cast<int>(ghostsToSend[0].size());
    //South vector creation 
    s = 0;
    for(int i = voxBoxDim.x + 1; i < (voxBoxDim.x*2) - 2; i++){
        v = pVoxArr.at(i);
        if(v.layer == layer){
            ghostsToSend[1].push_back(v.index.x);
            ghostsToSend[1].push_back(v.index.y);
            ghostsToSend[1].push_back(v.index.z);

            if(v.isParticle || v.isBoundary){
                m = v.isParticle ? 1 : 2;
                ghostsToSend[1].push_back(m);
            }
            else{
                ghostsToSend[1].push_back(0);
            }

            for(auto& ori : v.origins){
                ghostsToSend[1].push_back(ori);
            }
            ghostsToSend[1].push_back(-1);
            s++;
        }
    }
    southSendSize = static_cast<int>(ghostsToSend[1].size()); 
    //North/South sending along main bits
    if(this->dcomp.mainPlusY != this->dcomp.mainDomainNumber){
        cout << "We want to send " << n << " voxels to the north from rank " << this->dcomp.rank << " to rank " << this->dcomp.mainPlusY << endl;
        cout << "We want to send " << s << " voxels to the north from rank " << this->dcomp.rank << " to rank " << this->dcomp.mainMinusY << endl;
        MPI_Send(&northSendSize, 1, MPI_INT, this->dcomp.mainPlusY, 0, MPI_COMM_WORLD);
        MPI_Recv(&southRecvSize, 1, MPI_INT, this->dcomp.mainMinusY, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(northSendSize != 0){
            MPI_Send(&(ghostsToSend[0])[0], northSendSize, MPI_INT, this->dcomp.mainPlusY, 0, MPI_COMM_WORLD);
        }
        if(southRecvSize != 0){
            ghostsToRecv[1].reserve(southRecvSize);
            MPI_Recv(&(ghostsToRecv[1])[0], southRecvSize, MPI_INT, this->dcomp.mainMinusY, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if(northSendSize != 0){
            MPI_Send(&(ghostsToSend[0])[0], northSendSize, MPI_INT, this->dcomp.mainPlusY, 0, MPI_COMM_WORLD);
        }
        if(southRecvSize != 0){
            ghostsToRecv[1].reserve(southRecvSize);
            MPI_Recv(&(ghostsToRecv[1])[0], southRecvSize, MPI_INT, this->dcomp.mainMinusY, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Send(&southSendSize, 1, MPI_INT, this->dcomp.mainMinusY, 0, MPI_COMM_WORLD);
        MPI_Recv(&northRecvSize, 1, MPI_INT, this->dcomp.mainPlusY, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(southSendSize != 0){
            MPI_Send(&(ghostsToSend[1])[0], southSendSize, MPI_INT, this->dcomp.mainMinusY, 0, MPI_COMM_WORLD);
        }
        if(northRecvSize != 0){
            ghostsToRecv[0].reserve(northRecvSize);
            MPI_Recv(&(ghostsToRecv[0])[0], northRecvSize, MPI_INT, this->dcomp.mainPlusY, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    //integrateGhosts(&ghostsToRecv[0]);
    //integrateGhosts(&ghostsToRecv[1]);

    //East vector creation
    e = 0;
    for(int i = voxBoxDim.x - 2; i < pVoxSize; i += voxBoxDim.x){
        v = pVoxArr.at(i);
        if(v.layer == layer){
            ghostsToSend[2].push_back(v.index.x);
            ghostsToSend[2].push_back(v.index.y);
            ghostsToSend[2].push_back(v.index.z);

            if(v.isParticle || v.isBoundary){
                m = v.isParticle ? 1 : 2;
                ghostsToSend[2].push_back(m);
            }
            else{
                ghostsToSend[2].push_back(0);
            }

            for(auto& ori : v.origins){
                ghostsToSend[2].push_back(ori);
            }
            ghostsToSend[2].push_back(-1);
            e++;
        }
    }
    eastSendSize = static_cast<int>(ghostsToSend[2].size());
    //West vector creation
    w = 0;
    for(int i = 1; i < pVoxSize; i += voxBoxDim.x){
        v = pVoxArr.at(i);
        if(v.layer == layer){
            ghostsToSend[3].push_back(v.index.x);
            ghostsToSend[3].push_back(v.index.y);
            ghostsToSend[3].push_back(v.index.z);

            if(v.isParticle || v.isBoundary){
                m = v.isParticle ? 1 : 2;
                ghostsToSend[3].push_back(m);
            }
            else{
                ghostsToSend[3].push_back(0);
            }

            for(auto& ori : v.origins){
                ghostsToSend[3].push_back(ori);
            }
            ghostsToSend[3].push_back(-1);
            w++;
        }
    }
    westSendSize = static_cast<int>(ghostsToSend[3].size());
    //NOTE: in future make it so before this we split along subsections and also send only to the local x neighbors so it works with non perfect squares
    if(this->dcomp.mainPlusX != this->dcomp.mainDomainNumber){
        cout << "We want to send " << e << " voxels to the east from rank " << this->dcomp.rank << " to rank " << this->dcomp.mainPlusX << endl;
        cout << "We want to send " << w << " voxels to the west from rank " << this->dcomp.rank << " to rank " << this->dcomp.mainMinusX << endl;
        MPI_Send(&eastSendSize, 1, MPI_INT, this->dcomp.mainPlusX, 0, MPI_COMM_WORLD);
        MPI_Recv(&westRecvSize, 1, MPI_INT, this->dcomp.mainMinusX, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(eastSendSize != 0){
            MPI_Send(&(ghostsToSend[2])[0], eastSendSize, MPI_INT, this->dcomp.mainPlusX, 0, MPI_COMM_WORLD);
        }
        if(westRecvSize != 0){
            ghostsToRecv[3].reserve(westRecvSize);
            MPI_Recv(&(ghostsToRecv[3])[0], westRecvSize, MPI_INT, this->dcomp.mainMinusX, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Send(&westSendSize, 1, MPI_INT, this->dcomp.mainMinusX, 0, MPI_COMM_WORLD);
        MPI_Recv(&eastRecvSize, 1, MPI_INT, this->dcomp.mainPlusX, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(westSendSize != 0){
            MPI_Send(&(ghostsToSend[3])[0], westSendSize, MPI_INT, this->dcomp.mainMinusX, 0, MPI_COMM_WORLD);
        }
        if(eastRecvSize != 0){
            ghostsToRecv[2].reserve(eastRecvSize);
            MPI_Recv(&(ghostsToRecv[2])[0], eastRecvSize, MPI_INT, this->dcomp.mainPlusX, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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

    if(numberInLocal == 1){
        this->localXMin = xSpace * xMod;
        this->localXMax = (xMod == xLength - 1) ? (xLength): (localXMin + xSpace - 1);
    }
    else{
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

    mainNeighbors2D();
}


void DomainDecomposition::mainNeighbors2D(){
    int dim = floor(sqrt(mainCount));
    int z = floor(mainDomainNumber/(mainCount));
    int y = floor((mainDomainNumber- (z*mainCount))/dim);
    int x = mainDomainNumber - (z*mainCount) - (y*dim);

    int minusXtemp, plusXtemp;
    int minusYtemp, plusYtemp;
    int minusZtemp, plusZtemp;
 
    plusXtemp = (x != dim - 1) ? (1) : (-dim + 1);
    minusXtemp = (x != 0) ? (-1) : (dim - 1);
    minusYtemp = (y != dim - 1) ? (dim) : (-dim * (dim - 1));
    plusYtemp = (y != 0) ? (-dim) : (dim * (dim - 1));
    minusZtemp = -1;
    plusZtemp = -1;

    this->mainMinusX = mainDomainNumber + plusXtemp;
    this->mainPlusX = mainDomainNumber + minusXtemp;
    this->mainMinusY = mainDomainNumber + plusYtemp;
    this->mainPlusY = mainDomainNumber + minusYtemp;

    this->mainPlusYPlusX = mainDomainNumber + plusYtemp + plusXtemp;
    this->mainPlusYMinusX = mainDomainNumber + plusYtemp + minusXtemp;
    this->mainMinusYPlusX = mainDomainNumber + minusYtemp + plusXtemp;
    this->mainMinusYMinusX = mainDomainNumber + minusYtemp + minusXtemp;

}

/*
void DomainDecomposition::sendData(&int data, int recvRank){

}

void DomainDecomposition::ReceiveData(&int data, int sendRank){

}
*/