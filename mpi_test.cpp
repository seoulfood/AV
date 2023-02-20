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
//using std::chrono;

class Functions{
    public:
        int rank;
        int P;
        DomainDecomposition dcomp;
    
        Functions(int rank, int P);
    
        void printRank();
        void receiveLeftNeighbor();
        void sendRightNeighbor();
        void splitUpDomain(int xLength, int yLength, int zLength, int print);

};


void Functions::splitUpDomain(int xLength, int yLength, int zLength, int print){
    this->dcomp = DomainDecomposition(this->rank, this->P);

    dcomp.divideSimBox2D(xLength, yLength);

    for(int x = dcomp.localXMin; x <= dcomp.localXMax; x++)
    {
        for(int y = dcomp.localYMin; y <= dcomp.localYMax; y++)
        {
            if(print == 1){
                cout << "[" << x << "][" << y << "][" << "0" << "]" << " is " << this->rank << endl;  
            }
            int a = 5;
        }
    }
}

Functions::Functions(int rank, int P){
    this->rank = rank;
    this->P = P;
}

void Functions::printRank()
{
    if(this->rank==0){
        cout << "Start: " << MPI_Wtime() << endl;
    }

    cout << "Hello from rank " << this->rank << endl;

}

void Functions::receiveLeftNeighbor()
{
    int left;
    int recv_from = (this->rank != 0) ? ((this->rank - 1)) : this->P - 1;
    MPI_Recv(&left, 1, MPI_INT, recv_from, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    cout << "Rank " << this->rank << "'s left neighbor is " << left << endl;

}

void Functions::sendRightNeighbor()
{
    int send_to = (this->rank != this->P - 1) ? ((this->rank + 1)) : 0;
    MPI_Send(&this->rank, 1, MPI_INT, send_to, 1, MPI_COMM_WORLD);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int P;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Functions func = Functions(rank, P);
    //func.printRank();
    //func.sendRightNeighbor();
    //func.receiveLeftNeighbor();
    func.splitUpDomain(18, 18, 0, atoi(argv[1]));

    MPI_Finalize();

    return 0;
}

