#include <vector>
#include <iostream>
#include <iomanip>

#include <chrono>
#include "mpi.h"
#include <math.h>


#include <list>
#include <cstring>    /* memset & co. */
#include <ctime>
#include <cassert>
#include <cuda.h>
#include <cuda_runtime.h>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <iostream>
#include <numeric>

#include<unistd.h>
#include <iostream>

//#include "input.hpp"
#include <iostream>
#include <fstream>

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <chrono>
#include <cuda.h>
#include <cuda_runtime.h>
#include <list>
MPI_Datatype type;
int rank, num_procs;
bool  server;
void pingpong(void *buffer , int i) {

    int size_s =0;

    std::list<double> times0;
    for (int k = 0; k < 10; ++k) {

        double t0 = MPI_Wtime();
        for (int j = 0; j < 10; ++j) {
            if (rank == 0) {
                MPI_Send(buffer, i, type, 1, 6, MPI_COMM_WORLD);
//                MPI_Recv(buffer, i, type, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank == 1) {
                MPI_Recv(buffer, i, type, 0, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//                MPI_Send(buffer, i, type, 0, 0, MPI_COMM_WORLD);
            }
        }
        double tfinal = (MPI_Wtime() - t0) / (10);

        times0.push_back(tfinal);


    }

    times0.sort();
    int count = 0;
    double low = 0;
    double med = 0;
    double high = 0;
    double total = 0;
    for (const auto &item: times0) {
        total+=item;
        if (count == 0) {
            low = item;
        }
        if (count == (10)/2) {
            med = item;
        }
        if (count == 10-1) {
            high = item;
        }

        count++;

    }

    double mean0 =total/times0.size();

    if (rank==1){
        printf("pingpong--%i,%15.9f,%15.9f,%15.9f,%15.9f\n", i * size_s, mean0, low, med, high);
        fflush(stdout);

    }



}
static void *mpi_cuda_malloc(size_t count) {
    void *d_data1;
    cudaMalloc((void **) &d_data1, count * sizeof(float));
//    cudaMalloc((void **) &d_data0, data_size * sizeof(float));

    float *h_data = (float *) malloc(count * sizeof(float));

    for (int i = 0; i < count; ++i) {
        h_data[i] = i * 1.0f;
    }


    cudaMemcpy(d_data1, h_data, count * sizeof(float), cudaMemcpyHostToDevice);
//    free(h_data);
    return d_data1;
}




int main(int argc, char *argv[]) {




    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


    void* buffer = mpi_cuda_malloc(10000000);

    server = rank == 0;

    int k = 2;
    int j = 3;
    int i = 5;
    int amount = 60960;

    MPI_Datatype oldType = MPI_FLOAT;

   int ret = MPI_Type_vector(k, j, i, oldType, &type);
    assert(ret == MPI_SUCCESS);
    ret = MPI_Type_commit(&type);

    pingpong(buffer, amount);

    return 0;
}

