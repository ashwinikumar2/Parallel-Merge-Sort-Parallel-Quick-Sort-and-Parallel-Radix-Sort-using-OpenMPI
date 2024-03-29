#include <iostream>
#include <ctime> 
#include <time.h> 
#define LOADSIZE	4
#include <mpi.h>
#include "sort.h"

const int MINCOUNT = 1000000;
const int MAXCOUNT = 1000000;

inline int randomCount() { return MINCOUNT + rand()%(MAXCOUNT-MINCOUNT+1); }
char randomChar() {
    
    std::string str = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    int pos = rand()%str.size();
    return str[pos];
}

pSort::dataType* generate(long num_of_records)
{
    //Creating an array of structure of type "dataType" declared in sort.h

    pSort::dataType *data = new pSort::dataType[num_of_records];

    srand(time(0));
    for (int i = 0; i < num_of_records; i++) {
        data[i].key = rand();
        for(int j=0; j<LOADSIZE; j++) data[i].payload[j] = randomChar();
        // printf("(%d: %d %c%c%c%c) ",  i,  data[i].key,  data[i].payload[0],  data[i].payload[1], data[i].payload[2],  data[i].payload[3]);
    }

    return data;
}

bool check_sorted(pSort::dataType *test_data,int num_of_records)
{
   int rank, size;

   MPI_Comm_rank(MPI_COMM_WORLD , &rank);
   MPI_Comm_size(MPI_COMM_WORLD , &size);
    // I expect the same number num_of_records here as the input set at this rank
   if(rank != size-1) {
      if(MPI_SUCCESS != MPI_Send(&test_data[num_of_records-1], 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD)) {
         std::cerr << "End Data sending fails" << rank << std::endl;
         return false;
      }
   }
   if(rank != 0) {
      int previouskey;
      MPI_Status status;
      if(MPI_SUCCESS != MPI_Recv(&previouskey, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, &status)) {
         std::cerr << "Previous receiving fails " << rank << std::endl;
         return false;
      }
      if(test_data[0].key < previouskey) return false;
   }
   for(int i=1; i<num_of_records; i++)
      if(test_data[i].key < test_data[i-1].key) return false;

   // Also check that payload intact for each key -- TO BE IMPLEMENTED

   return true;
}

void runExperiment(pSort sorter, int num_of_records=0, pSort::SortType type = pSort::BEST, bool term=true)
{
    /*Processing command line arguments supplied with mpirun*/
    if(num_of_records == 0) num_of_records = randomCount();
    
    //Creating an array of structure of type "dataType" declared in sort.h
    pSort::dataType *test_data  = generate(num_of_records);

    /*Calling functions defined in pSort library to sort records stored in test_data[]*/
    time_t begin,end;
    time(&begin);
    sorter.sort(test_data, num_of_records, type);
    time(&end);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    double timetaken = difftime(end,begin);
    if(check_sorted(test_data,num_of_records)) {
       std::cout << type << " Successful in " << timetaken <<" processor= "<< my_rank<<std::endl;
    }
    if(term) delete test_data;
}

int main(int argc, char *argv[]){
    

    pSort sorter;

    //Calling your init() to set up MPI	
    sorter.init();

    /*=================================================================*/
    runExperiment(sorter, 0, pSort::RADIX); // For example

    //Calling your close() to finalize MPI 
    sorter.close();

    return 0;
}
