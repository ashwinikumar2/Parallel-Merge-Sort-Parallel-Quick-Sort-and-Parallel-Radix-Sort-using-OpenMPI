#include <mpi.h>
#include "sort.h"
#include<iostream>
#include <vector>
#include <cmath>
#define LOADSIZE 4
using namespace std;

MPI_Datatype newtype;           //A Custom data type representing psort::datatype

void pSort::init()              //Used to initialize the program
{
	MPI_Init(NULL, NULL);
	int size_of_type;
    MPI_Type_size(MPI_CHAR, &size_of_type);
    MPI_Aint displacements[2]  = {0,size_of_type};
    int block_lengths[2]  = {1, LOADSIZE};
    MPI_Datatype type_of_blocks[2] = {MPI_INT,MPI_CHAR};
   
    MPI_Type_create_struct(2,block_lengths,displacements,type_of_blocks,&newtype);      //creating a new struct for psort::datatype consisting of two blocks of lengths 1 and LOADSIZE
    MPI_Type_commit(&newtype);                      //commiting the new custom data type named as 'newtype'
}

void pSort::close()             //Used to close the MPI program
{
	MPI_Finalize();
}

/////////////////////////////////////////////////////////////////////////////////////MERGE SORT/////////////////////////////////////////////////////////////////////////
pSort::dataType* merge(pSort::dataType *arr, pSort::dataType *arr1 , int l, int m, int r)                       //used to merge two arrays by comparing the first elements(in iteration) of left portion as compared to right portion
{                                                                                                               //or we can say the usual merge used in serial merge sort algorithm    

    for(int i = 0; i < m-l+1; i++)
    {    
        arr1[i] = arr[l + i];
    }
    for(int j = 0; j < r-m; j++)
    {
        arr1[j+m-l+1] = arr[m + 1 + j];
    }
    int i = 0,j=0,k=l; 
    
    while (i < m - l + 1 && j < r - m)                  //while we are iterating in the left portion of arr i.e. from 0 to m-l and iterating in the right portion of arr i.e. m-l+1 to r-m
    {
        if (arr1[i].key <= arr1[j+m-l+1].key)                       
        {
            arr[k] = arr1[i];
            i=i+1;
        }
        else
        {
            arr[k] = arr1[j+m-l+1];
            j=j+1;
        }
        k=k+1;
    }

    while (i < m - l + 1) 
    {
        arr[k] = arr1[i];
        i=i+1;
        k=k+1;
    }
    while (j < r - m)
    {
        arr[k] = arr1[j+m-l+1];
        j=j+1;
        k=k+1;
    }

    return arr;                                         
}
void mergeSort(pSort::dataType *arr, int l, int r)                                                  //serial merge sort used for doing merge sort serially on an array by the usual merge sort algorithm
{
    if (l < r)
    {   pSort:: dataType *temp=arr;

        int m = (l+r) / 2;
        pSort::dataType *arr1=new pSort::dataType[r-l+1];
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);

        pSort::dataType *output=merge(arr,arr1, l, m, r);
        for (int i = l; i < r; i++)                                                                 //used for setting the pointer values of what was obtained as data at initialization of this function
        {
            /* code */
            temp[i]=output[i];
        }
    }
}

// void parallel_merge(int *data, int ndata)
void parallel_merge(pSort::dataType *data, int ndata)                                               //parallel version of merge sort algorithm 
{
    pSort::dataType *temp=data;
    mergeSort(data, 0 , ndata-1);                                                                   //locally use serial merge sort algorithm         

    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    
    int world_size;

    MPI_Comm_size( MPI_COMM_WORLD, &world_size ); 
    if (my_rank == 0)                                                                               //if we are at the first processor(i.e. rank 0) then gather all the counts of how much data every other processor
    {                                                                                               //is going to send to it
        int *rcounts =new int[world_size];
        MPI_Gather( &ndata, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);                    //gather command for knowing the size that is going to be recieved by 0th processor

        int count_of_collect=0;
        for(int i=0;i<world_size;i++)
        {
            count_of_collect+=rcounts[i];
        }
        int *displs=new int[world_size];                    //displacement array
        // int x=0;
        for(int i=0;i<world_size;i++)
        {
            
            displs[i]=0;
        }
        for(int i=0;i<world_size;i++)
        {
            if(i==0)
                displs[i]=0;
            else{
            displs[i]=displs[i-1]+rcounts[i-1];
            }
        }
        pSort::dataType *collect=(pSort::dataType *)malloc(count_of_collect*sizeof(pSort::dataType));           //array used for collecting the elements sent by every processor to 0th processor

        MPI_Gatherv(data, ndata, newtype, collect, rcounts, displs, newtype, 0, MPI_COMM_WORLD);                //using gatherv to obtain elements from every other processors to parallely merge them  
                
        if(world_size>1)                                          //if number of processors are greater than 1 
        {
            int low=0,high=0,mid=rcounts[0]-1;
            pSort::dataType *arr1=new pSort::dataType[count_of_collect];
            for(int i=0;i<world_size;i++)
            {

                if(i+1<world_size)
                {

                    high=rcounts[i+1]+mid;
                    collect=merge(collect, arr1,low, mid, high);
                    mid=high;
                }
            }

            //scattering the data
            
            int *displs=new int[world_size];                        //displacement array
            int x=0;
            for(int i=0;i<world_size;i++)
            {
                
                displs[i]=0;
            }
            for(int i=0;i<world_size;i++)
            {
                if(i==0)
                    displs[i]=0;
                else{
                displs[i]=displs[i-1]+rcounts[i-1];
                }
            } 

            MPI_Scatterv(collect, rcounts, displs, newtype , data, ndata, newtype , 0, MPI_COMM_WORLD) ;            //scatter the sorted data to every other processor according to what was supposed to be their original size
                                                                                                                    //i.e. every processor contains the same number of elements it contained originally
            for(int i=0;i<ndata;i++)                                                                                //used for setting the pointer values of what was obtained as data at initialization of this function
            {
                temp[i]=data[i];
            }
        }
        
    }

    else
    {
        MPI_Gather( &ndata, 1, MPI_INT, NULL, 10, MPI_INT, 0, MPI_COMM_WORLD);                  //gather command for sending the size of elements that current processor own

        MPI_Gatherv(data, ndata, newtype, NULL, NULL, 0, newtype, 0, MPI_COMM_WORLD);           //using gatherv to send elements from current processor to 0 for parallel merge  

        MPI_Scatterv(NULL, NULL, NULL, newtype , data, ndata, newtype, 0, MPI_COMM_WORLD) ;

        for(int i=0;i<ndata;i++)                                                                //used for setting the pointer values of what was obtained as data at initialization of this function
        {
            temp[i]=data[i];
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////MERGE SORT//////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////QUICK SORT/////////////////////////////////////////////////////////////////////////
pSort::dataType* merge_arrays(pSort::dataType *left, pSort::dataType *right, int leftsize, int rightsize)           //used to merge two arrays as recieved. It does not change the order of elements just appends the right arrau ahead of left array
{
    pSort::dataType *output=new pSort::dataType[leftsize+rightsize];
    for(int i=0;i<leftsize;i++)
    {
        output[i]= left[i];

    }

    for(int i=0;i<rightsize;i++)
    {
        output[i+leftsize]= right[i];
        
    }

    return output;
}
pSort::dataType* sequential_quick_sort(pSort::dataType *data, int ndata, int low, int high)                         //serial quick sort implementing the usual quick sort algorithm by finding the partition index and 
{                                                                                                       //then recursively use this on left and right portions seperated by partition index
    if(low<high)
    {
        // pSort::dataType *temp=data;

        pSort::dataType pivot = data[high];                         //taking pivot to be the value of key of last element in the recieved array  
        int i = (low - 1); 
      
        for (int j=low;j<high;j++)  
        {  
            if (data[j].key < pivot.key)  
            {  
                i++; 
                pSort::dataType x=data[i];
                data[i]=data[j];                                    //swapping values of array elements at indexes i and j
                data[j]=x; 
            }  
        }  
        pSort::dataType x=data[i+1];
        data[i+1]=data[high];                                       //swapping values of array elements at indexes i+1 and high
        data[high]=x;

        int pi=i+1;
        data= sequential_quick_sort(data, ndata, low, pi-1);
        data=sequential_quick_sort(data, ndata, pi+1, high);
        return data;
        
    }
    return data;
}

void quick_sort(pSort::dataType *data, int &ndata, MPI_Comm comm)                                       //parallel version of quick sort
{
    pSort::dataType *temp=data;
    int my_rank,ndata_new;
    int recv_count,  i;
    MPI_Status status;
    int numProcess;    
    MPI_Comm_size(comm, &numProcess);
    MPI_Comm_rank(comm, &my_rank);
    pSort::dataType pivot1;             
    if(numProcess==1){                                                                                  //if number of processors are 1 in this communicator then apply serial quick sort                                      
    	if(ndata!=0)
        data=sequential_quick_sort(data, ndata, 0, ndata-1); 

        for(int i=0;i<ndata;i++)                                                    //used for setting the pointer values of what was obtained as data at initialization of this function
        {
            temp[i]=data[i];
        }
        return;
    }

    int low=0,high=numProcess-1;
    if (my_rank==(numProcess-1))                                                    //if this is the last processor 
    {
        pSort::dataType pivot=data[ndata/2];
        MPI_Bcast(&pivot, 1,newtype,my_rank, comm);                                 //then broadcast to every other processor the pivot value i.e. at ndata/2 index 
        int x1=0;
        vector<pSort::dataType> left1;
        vector<pSort::dataType> right1;
        

        for(int i=0;i<ndata;i++)
        {
            if(data[i].key<=pivot.key)
              left1.push_back(data[i]);
            else
              right1.push_back(data[i]);
        }

        pSort::dataType left[left1.size()];
        pSort::dataType *right= new pSort::dataType[right1.size()];
        for(int i=0;i<left1.size();i++)
        {            
            left[i]=left1[i];
        }
        for(int i=0;i<right1.size();i++)
        {
            right[i]=right1[i];
        }

        MPI_Send( left, left1.size(), newtype, x1, 10, comm);                           //send the values less than pivot to the partner process which belongs to lower half list of processors

        MPI_Status statusx;
        MPI_Probe(x1,10,comm,&statusx);
        int recv_count;
        MPI_Get_count(&statusx, newtype, &recv_count);
        pSort::dataType *recieving_high_list_from_partner=new pSort::dataType[recv_count];
        
        MPI_Recv( recieving_high_list_from_partner, recv_count, newtype, x1,  10, comm, &status);     //recieve the values greater than pivot from the partner process which belongs to lower half list of processors         
        data= merge_arrays(right, recieving_high_list_from_partner, right1.size(), recv_count);

        ndata_new=right1.size()+recv_count;
        
        pivot1=pivot;

        if(my_rank==(numProcess/2)+1 && numProcess%2!=0)
        {

            MPI_Status statusx1;
            MPI_Probe(numProcess/2,10,comm,&statusx1);
            int recv_count1;
            MPI_Get_count(&statusx1, newtype, &recv_count);
            pSort::dataType *recieving_high_list_from_partner1=new pSort::dataType[recv_count];
            MPI_Recv( recieving_high_list_from_partner1, recv_count, newtype, numProcess/2,  MPI_ANY_TAG, comm, &status); //recieve the values greater than pivot from the process at rank= number_of_processors/2 which does not have any partner process

            data= merge_arrays(recieving_high_list_from_partner1,data, recv_count, ndata_new);
            ndata_new+=recv_count;
            delete[] recieving_high_list_from_partner1;
        }

        left1.clear();
        right1.clear();
        delete[] right;
        delete[] recieving_high_list_from_partner;
        
    }
    //reorganising data from every processor such that upper processors contain data greater than pivot and lower processors contains data less than or equal to pivot
    else 
    {
        if(my_rank>((numProcess-1)/2) )     //upper half listed processor send data to partner processor in lower half and recieve high list data(data greater than pivot) from partner processor
        {
            pSort::dataType pivot;
            MPI_Bcast(&pivot, 1,newtype, numProcess-1, comm);                   //calling broadcast so that this processor can recieve the pivot value from last processor
            vector<pSort::dataType> left1;
            vector<pSort::dataType> right1;
            for(int i=0;i<ndata;i++)
            {

                if(data[i].key<=pivot.key)
                  left1.push_back(data[i]);
                else
                  right1.push_back(data[i]);
            }

            pSort::dataType left[left1.size()];
            pSort::dataType *right=new pSort::dataType[right1.size()];
            for(int i=0;i<left1.size();i++)
            {
                left[i]=left1[i];
            }
            for(int i=0;i<right1.size();i++)
            {
                right[i]=right1[i];
            }
            data=right;
            int x1=numProcess-1-my_rank;

            MPI_Send( left, left1.size(), newtype, x1, 10, comm);                   //send the values less than pivot to the partner process which belongs to lower half list of processors

            MPI_Status statusx;
            MPI_Probe(x1,10,comm,&statusx);
            int recv_count;
            MPI_Get_count(&statusx, newtype, &recv_count);
            pSort::dataType *recieving_high_list_from_partner=new pSort::dataType[recv_count];
            MPI_Recv( recieving_high_list_from_partner, recv_count, newtype, x1,  MPI_ANY_TAG, comm, &status);    //recieve the values greater than pivot from the partner process which belongs to lower half list of processors         

            data= merge_arrays(data, recieving_high_list_from_partner, right1.size(),recv_count);
            ndata_new=right1.size()+recv_count;
            pivot1=pivot;

            if(my_rank==(numProcess/2)+1 && numProcess%2!=0)                        
            {
                MPI_Status statusx1;
                MPI_Probe(numProcess/2,10,comm,&statusx1);
                int recv_count1;
                MPI_Get_count(&statusx1, newtype, &recv_count);
                pSort::dataType *recieving_high_list_from_partner1=new pSort::dataType[recv_count];
                MPI_Recv( recieving_high_list_from_partner1, recv_count, newtype, numProcess/2,  MPI_ANY_TAG, comm, &status);           //recieve the values greater than pivot from the process at rank= number_of_processors/2 which does not have any partner process

                data= merge_arrays(recieving_high_list_from_partner1,data, recv_count, ndata_new);
                ndata_new+=recv_count;
                delete[] recieving_high_list_from_partner1;
            }
            left1.clear();
            right1.clear();
            delete[] right;
            delete[] recieving_high_list_from_partner;
          

        }
        else                        //lower half processor receive low list data from partner processor (in upper half of list of processors) and sends a high list of data to partner processor
        {
            pSort::dataType pivot;
            MPI_Bcast(&pivot, 1,newtype, numProcess-1, comm);                       //calling broadcast so that this processor can recieve the pivot value from last processor
            
            int recieving_rank=(numProcess-1-my_rank);
            
            if(my_rank==numProcess/2 && numProcess%2!=0)                            //if the number of processors are odd which leaves the processor at rank=number_of_processors/2 with no partner processor
            {
                vector<pSort::dataType> left1;
                vector<pSort::dataType> right1;
                
                for(int i=0;i<ndata;i++)
                {
                    
                    if(data[i].key<=pivot.key)
                      left1.push_back(data[i]);
                    else
                      right1.push_back(data[i]);
                }

                pSort::dataType *left=new pSort::dataType[left1.size()];
                pSort::dataType right[right1.size()];
                for(int i=0;i<left1.size();i++)
                {
                    left[i]=left1[i];
                }
                for(int i=0;i<right1.size();i++)
                {
                    right[i]=right1[i];
                }
                data=left;
                ndata_new=left1.size();
                MPI_Send( right, right1.size(), newtype, (my_rank+1)%numProcess, 10, comm);         //send values greaer than pivot to the next processor in this special case
                
                left1.clear();
                right1.clear();
                
            }
            else
            {
                MPI_Status statusx;
                MPI_Probe(recieving_rank,10,comm,&statusx);
                int recv_count;
                MPI_Get_count(&statusx, newtype, &recv_count);
                
                pSort::dataType *numbertoreceive=new pSort::dataType[recv_count];
                MPI_Recv( numbertoreceive, recv_count, newtype, recieving_rank,  MPI_ANY_TAG, comm, &status);               //recieve the values less than pivot from the partner process in the upper half 
                    
                vector<pSort::dataType> left1;
                vector<pSort::dataType> right1;
                
                for(int i=0;i<ndata;i++)
                {
                    
                    if(data[i].key<=pivot.key)
                      left1.push_back(data[i]);
                    else
                      right1.push_back(data[i]);
                }

                pSort::dataType *left=new pSort::dataType[left1.size()];
                pSort::dataType right[right1.size()];
                for(int i=0;i<left1.size();i++)
                {
                    left[i]=left1[i];
                    // if(i>left1.size()-1)
        		// cout<<"Left IS WRONG at line 438"<<endl;
                }
                for(int i=0;i<right1.size();i++)
                {
                    right[i]=right1[i];
                }
                data=merge_arrays(left, numbertoreceive, left1.size(),recv_count);
                ndata_new=left1.size()+recv_count;
                
                MPI_Send( right, right1.size(), newtype, recieving_rank, 10, comm);                                 //send values greaer than pivot to the partner processor in the upper half
                left1.clear();
                right1.clear();
                // free(right);
                delete[] left;
                delete[] numbertoreceive;

            }

            pivot1=pivot;
            //free(array);
        }
    }

    MPI_Barrier(comm);
    int color = my_rank/((numProcess+1)>>1);
    MPI_Comm row_comm;
    
    MPI_Comm_split(comm, color,my_rank, &row_comm);                                     //split the current communicator to 2 communicators
    // cout<<"called for recursion in my_rank= "<<my_rank<<endl;
    quick_sort(data,ndata_new,row_comm);                                                //call recursion on the new communicators formed
    

    //every processor need not contain same number of elements it originally used to so moving the elements 
    if(my_rank==0)                              
    {
        int extra_size=ndata_new-ndata;
        if(extra_size<0)
            extra_size=0;
        pSort::dataType *send_extra_data=new pSort::dataType[extra_size];

        for(int i=ndata;i<ndata+extra_size;i++)
        {
            send_extra_data[i-(ndata)]=data[i];
        }

        MPI_Send( send_extra_data, extra_size, newtype, (my_rank+1)%numProcess, 10, comm);              //send extra elements(in case there are extra) to the next processor
        if(extra_size!=0)
        {
            pSort::dataType *restore_data=new pSort::dataType[ndata];
            for(int i=0;i<ndata;i++)
            {
                restore_data[i]=data[i];
            }

            data= restore_data;
            ndata_new=ndata;
        }

        MPI_Status statusx;
        MPI_Probe((my_rank+1)%numProcess,10,comm,&statusx);
        int recv_count;
        MPI_Get_count(&statusx, newtype, &recv_count);
        pSort::dataType *numbertoreceive=new pSort::dataType[recv_count];
        MPI_Recv( numbertoreceive, recv_count, newtype, (my_rank+1)%numProcess,  MPI_ANY_TAG, comm, &status);           //recieve the extra elements from next processor
        
        data=merge_arrays(data,numbertoreceive, ndata_new,recv_count);
        ndata_new+=recv_count;
        for(int i=0;i<ndata_new;i++)
        {
            temp[i]=data[i];
        }


        delete[] send_extra_data;
        delete[] numbertoreceive;

    }
    else
    {

        int recieving_rank=(my_rank-1);
        MPI_Status statusx;
        MPI_Probe(recieving_rank,10,comm,&statusx);
        int recv_count;
        MPI_Get_count(&statusx, newtype, &recv_count);
        pSort::dataType *numbertoreceive=new pSort::dataType[recv_count];
        MPI_Recv( numbertoreceive, recv_count, newtype, recieving_rank,  10, comm, &status);                    //recieve the extra elements from previous processor

        data=merge_arrays(numbertoreceive,data, recv_count, ndata_new);
        ndata_new+=recv_count;

        if(my_rank!=numProcess-1)
        {
            int extra_size=ndata_new-ndata;
            if(extra_size<0)
                extra_size=0;
            pSort::dataType *send_extra_data=new pSort::dataType[extra_size];

            for(int i=0;i<extra_size;i++)
            {
                send_extra_data[i]=data[i+ndata];
            }
            
            MPI_Send( send_extra_data, extra_size, newtype, (my_rank+1)%numProcess, 10, comm);          //send extra elements(in case there are extra) to the next processor
            if(extra_size!=0)
            {
                pSort::dataType *restore_data=new pSort::dataType[ndata];
                for(int i=0;i<ndata;i++)
                {
                    restore_data[i]=data[i];
                }

                data= restore_data;
                ndata_new=ndata;
            }
            delete[] send_extra_data;
            delete[] numbertoreceive;
        }
        else                                                                                //i.e. it is the last processor
        {
            int extra_size=ndata_new-ndata;
            if(extra_size<0)
                extra_size=0;
            pSort::dataType *send_extra_data=new pSort::dataType[extra_size];

            for(int i=0;i<extra_size;i++)
            {
                send_extra_data[i]=data[i];
            }

            MPI_Send( send_extra_data, extra_size, newtype, (my_rank-1), 10, comm);             //send extra elements to previous processor 
            
            if(extra_size!=0)
            {
                pSort::dataType *restore_data=new pSort::dataType[ndata];
                for(int i=0;i<ndata;i++)
                {
                    restore_data[i]=data[i+extra_size];
                }

                data= restore_data;
                ndata_new=ndata;
            }
            delete[] send_extra_data;
            delete[] numbertoreceive;

        }

        if(my_rank!=numProcess-1)                                           //if it is not the last processor 
        {
        	int recieving_rank=(my_rank+1);
	        MPI_Status statusx;
	        MPI_Probe(recieving_rank,10,comm,&statusx);
	        int recv_count;
	        MPI_Get_count(&statusx, newtype, &recv_count);
	        pSort::dataType *numbertoreceive=new pSort::dataType[recv_count];
	        MPI_Recv( numbertoreceive, recv_count, newtype, recieving_rank,  10, comm, &status);           //recieve extra elements from next processor

	        data=merge_arrays(data,numbertoreceive, ndata_new,recv_count);
	        ndata_new+=recv_count;

	        // sequential_quick_sort(data, ndata_new, 0, ndata_new-1);
	        int extra_size=ndata_new-ndata;
	        if(extra_size<0)
	            extra_size=0;
	        pSort::dataType *send_extra_data=new pSort::dataType[extra_size];

	        for(int i=0;i<extra_size;i++)
	        {
	            send_extra_data[i]=data[i];
	        }

	        // if(my_rank!=numProcess-1)
	        MPI_Send( send_extra_data, extra_size, newtype, (my_rank-1), 10, comm);                    //send extra elements to previous processor
	        
	        if(extra_size!=0)
	        {
	            pSort::dataType *restore_data=new pSort::dataType[ndata];
	            for(int i=0;i<ndata;i++)
	            {
	                restore_data[i]=data[i+extra_size];
	            }

	            data= restore_data;
	            ndata_new=ndata;
	            // delete[] restore_data;
	        }

            for(int i1=0;i1<extra_size;i1++)
            {
                send_extra_data[i1]=numbertoreceive[i1+recv_count-extra_size];
            }
            pSort::dataType *collect_to_keep=new pSort::dataType[recv_count- extra_size];
            for(int i1=0;i1<recv_count- extra_size;i1++)
            {
                collect_to_keep[i1]=numbertoreceive[i1];
                // v.push_back(collect_to_keep[i1]);
            
            }

            delete[] numbertoreceive;
            delete[] send_extra_data;
            delete[] collect_to_keep;

           
        }

         for(int i=0;i<ndata_new;i++)                                                       //used for setting the pointer values of what was obtained as data at initialization of this function
        
        
        {
            temp[i]=data[i];
            
        }
        
    }
    
        

    for(int i=0;i<ndata;i++)                                                                    //used for setting the pointer values of what was obtained as data at initialization of this function
    {
        temp[i]=data[i];
    }
}

/////////////////////////////////////////////////////////////////////////////////////////QUICK SORT///////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////RADIX SORT///////////////////////////////////////////////////////////////////////
int number_of_digits(int n)                                     //returns the number of digits in an integer(without sign)
{
    if(n>0)
    return floor(log10(n)+1);
	else if (n==0)
		return 1;
    else{
        int x= floor(log10(-n)+1);
        if(x<0)
            return (-x);
        else
            return x;
    }   
}

void radix_sort(pSort::dataType *data, int ndata)                               //parallel version of radix sort
{

	
    pSort::dataType *temp=data;
    pSort::dataType *data_duplicate=data;
    int my_rank;
    int recv_count,  i;
    MPI_Status status;
    int numProcess;    
    MPI_Comm_size(MPI_COMM_WORLD, &numProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int ndata_new=ndata;

    for(int i=1;i<=15;i++)                                                      //iterate till 15th digit
    {

        map<int, vector<pSort::dataType>> m;                                    //a map containing 10 buckets for each digit 0-9
        vector<pSort::dataType> numbers_not_in_map;
        for(int j=0;j<ndata_new;j++)
        {
            int y=data_duplicate[j].key;
            if(number_of_digits(y)>=i)
            {
                y=y/((int)(pow(10, i-1)));
                if(y<0)
                    y=-y;
                int x=(y)%(10);
                (m[x]).push_back(data_duplicate[j]);
            }
            else
            {
            	int x=0;
                (m[x]).push_back(data_duplicate[j]);
            }
        }
       
        
        for(int i1=0;i1<10;i1++)       //find out which digit is not there
        {
            if(!(m.find(i1)!=m.end()))       //table doesn't contain the char
            {
                vector<pSort::dataType> v;
                m[i1]=v;
            }
        }

        ndata_new=numbers_not_in_map.size();
        vector<pSort::dataType> v;
        v=numbers_not_in_map;
        for(auto itr=m.begin();itr!=m.end();++itr)                          //for each bucket send elements to 0 and then 0 will decide whether to move the extra elements to next processor or not
        {

            const int size_of_itr_th_digit_numbers=((*itr).second).size();
           
            int x=v.size();
            pSort::dataType *numbers_not_in_map_array=new pSort::dataType[x];
            for(int j=0;j<x;j++)
            {
                numbers_not_in_map_array[j]=v[j];
                
            }

            if(0 == my_rank)         //collect on *itr.first ranked process
            {
                int *rcounts=new int[numProcess];
                MPI_Allgather(&size_of_itr_th_digit_numbers, 1, MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);                 //every processor knows the elements on each processor in the current iteration of the buckets
                
                int count_of_collect=0;
                for(int k=0;k<numProcess;k++)
                {
                    count_of_collect+=rcounts[k];
                }
                
                pSort::dataType *numbers_to_be_sent=new pSort::dataType[size_of_itr_th_digit_numbers];
                for(int p=0;p<size_of_itr_th_digit_numbers;p++)
                {
                    numbers_to_be_sent[p]=((*itr).second)[p];
                }
                int *displs=new int[numProcess];
                // int x=0;
                for(int i=0;i<numProcess;i++)
                {
                    
                    displs[i]=0;
                }
                for(int i=0;i<numProcess;i++)
                {
                    if(i==0)
                        displs[i]=0;
                    else{
                    displs[i]=displs[i-1]+rcounts[i-1];
                    }
                }
                pSort::dataType *collect=new pSort::dataType[count_of_collect];
                MPI_Gatherv(numbers_to_be_sent, size_of_itr_th_digit_numbers, newtype, collect, rcounts, displs, newtype, my_rank, MPI_COMM_WORLD);             //collect this iteration's bucket on 0

                int extra_size=v.size()+count_of_collect-ndata;
                if(extra_size<0)
                    extra_size=0;
                pSort::dataType *send_extra_data=new pSort::dataType[extra_size];
                for(int i1=0;i1<extra_size;i1++)
                {
                    send_extra_data[i1]=collect[i1+count_of_collect- extra_size];
                    
                }
                pSort::dataType *collect_to_keep=new pSort::dataType[count_of_collect- extra_size];
                for(int i1=0;i1<count_of_collect- extra_size;i1++)
                {
                    collect_to_keep[i1]=collect[i1];

                }

				MPI_Send( send_extra_data, extra_size, newtype, (my_rank+1)%numProcess, 10, MPI_COMM_WORLD);                                    //send extra elements to next process
                data_duplicate=merge_arrays(numbers_not_in_map_array, collect_to_keep, v.size(), count_of_collect- extra_size);
                ndata_new=v.size()+ count_of_collect- extra_size;
                
                
                for(int i1=0;i1<count_of_collect- extra_size;i1++)
                {
                    v.push_back(collect_to_keep[i1]);
                }

                delete[] rcounts;
                delete[] numbers_to_be_sent;
                delete[] displs;
                delete[] collect;
                delete[] send_extra_data;
                delete[] collect_to_keep;
            }
            else
            {
                int *rcounts=new int[numProcess];
                MPI_Allgather(&size_of_itr_th_digit_numbers, 1, MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);             //every processor knows the elements on each processor in the current iteration of the buckets
                int count_of_collect=0;
                for(int k=0;k<numProcess;k++)
                {
                    count_of_collect+=rcounts[k];
                }
                
                pSort::dataType *numbers_to_be_sent=new pSort::dataType[size_of_itr_th_digit_numbers];
                for(int p=0;p<size_of_itr_th_digit_numbers;p++)
                {
                    numbers_to_be_sent[p]=((*itr).second)[p];
                }
                MPI_Gatherv(numbers_to_be_sent, size_of_itr_th_digit_numbers, newtype, NULL, NULL, NULL, newtype, 0, MPI_COMM_WORLD); //collect this iteration's bucket on 0

                int recieving_rank=(my_rank-1);
                

                MPI_Status statusx;
                MPI_Probe(recieving_rank,10,MPI_COMM_WORLD,&statusx);
                int recv_count;
                MPI_Get_count(&statusx, newtype, &recv_count);
                pSort::dataType *numbertoreceive=new pSort::dataType[recv_count];
                MPI_Recv( numbertoreceive, recv_count, newtype, recieving_rank,  MPI_ANY_TAG, MPI_COMM_WORLD, &status);         //recieve extra elements from previous processor

                int extra_size=v.size()+recv_count-ndata;
                if(extra_size<0)
                    extra_size=0;
                pSort::dataType *send_extra_data=new pSort::dataType[extra_size];
                
                for(int i1=0;i1<extra_size;i1++)
                {
                    send_extra_data[i1]=numbertoreceive[i1+recv_count-extra_size];
                }
                pSort::dataType *collect_to_keep=new pSort::dataType[recv_count- extra_size];
                for(int i1=0;i1<recv_count- extra_size;i1++)
                {
                    collect_to_keep[i1]=numbertoreceive[i1];
                
                }
                if(my_rank!=numProcess-1){
	                MPI_Send( send_extra_data, extra_size, newtype, (my_rank+1)%numProcess, 10, MPI_COMM_WORLD);           //send extra elements if it is not the last processor
	                
            	}

            	data_duplicate=merge_arrays(numbers_not_in_map_array, collect_to_keep, v.size(), recv_count- extra_size);
	                ndata_new=recv_count- extra_size + v.size();
                for(int i1=0;i1<recv_count-extra_size;i1++)
                {
                    v.push_back(collect_to_keep[i1]);

                }


                delete[] rcounts;
                delete[] numbers_to_be_sent;
                delete[] numbertoreceive;
                delete[] send_extra_data;
                delete[] collect_to_keep;
            }

            delete[] numbers_not_in_map_array;
        }  
        
        data_duplicate=new pSort::dataType[v.size()];
        for(int i1=0;i1<v.size();i1++)
        {
            data_duplicate[i1]=v[i1];
        }
        ndata_new=v.size();


        v.clear(); 
        m.clear(); 
        numbers_not_in_map.clear();
    }

    
    MPI_Barrier(MPI_COMM_WORLD);
    //now absolute values of every element's key has been sorted so now sorting based on signs
    map<int, vector<pSort::dataType>> m;                //contains two buckets- one for negative numbers and other for positive numbers
    vector<pSort::dataType> numbers_not_in_map;
    for(int j=0;j<ndata_new;j++)
    {
        pSort::dataType y=data_duplicate[j];
        if(y.key<0){
            (m[0]).push_back(y);
        }
        else
        {
            m[1].push_back(y);
        }
    }

    for(int i1=0;i1<2;i1++)       //find out which digit is not there
    {
        if(!(m.find(i1)!=m.end()))       //table doesn't contain the char
        {
            vector<pSort::dataType> v;
            m[i1]=v;
        }
    }
    ndata_new=0;
    vector<pSort::dataType> v;
    for(auto itr=m.begin();itr!=m.end();++itr)
    {
        const int size_of_itr_th_digit_numbers=((*itr).second).size();
           
        int x=v.size();
        pSort::dataType *numbers_not_in_map_array=new pSort::dataType[x];
        for(int j=0;j<x;j++)
        {
            numbers_not_in_map_array[j]=v[j];
            
        }
        if(0 == my_rank)         //collect on *itr.first ranked process
        {
            int *rcounts=new int[numProcess];
            MPI_Allgather(&size_of_itr_th_digit_numbers, 1, MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);//every processor knows the elements on each processor in the current iteration of the buckets
            
            int count_of_collect=0;
            for(int k=0;k<numProcess;k++)
            {
                count_of_collect+=rcounts[k];
            }
            
            pSort::dataType *numbers_to_be_sent=new pSort::dataType[size_of_itr_th_digit_numbers];
            for(int p=0;p<size_of_itr_th_digit_numbers;p++)
            {
                numbers_to_be_sent[p]=((*itr).second)[p];
            }
            int *displs=new int[numProcess];
            // int x=0;
            for(int i=0;i<numProcess;i++)
            {
                
                displs[i]=0;
            }
            for(int i=0;i<numProcess;i++)
            {
                if(i==0)
                    displs[i]=0;
                else{
                displs[i]=displs[i-1]+rcounts[i-1];
                }
            }
            pSort::dataType *collect=new pSort::dataType[count_of_collect];
            if((*itr).first==0){
                pSort::dataType *collect_opposite=new pSort::dataType[count_of_collect];
                MPI_Gatherv(numbers_to_be_sent, size_of_itr_th_digit_numbers, newtype, collect_opposite, rcounts, displs, newtype, my_rank, MPI_COMM_WORLD);//collect this iteration's bucket on 0
                
                for(int i1=0;i1< count_of_collect;i1++)                         //take opposite ordering of numbers if they belong to the negative numbers' bucket
                {
                    collect[i1]=collect_opposite[count_of_collect-1-i1];
                }
                delete[] collect_opposite;
            }
            else
            {
                MPI_Gatherv(numbers_to_be_sent, size_of_itr_th_digit_numbers, newtype, collect, rcounts, displs, newtype, my_rank, MPI_COMM_WORLD);//collect this iteration's bucket on 0
                
            }

            int extra_size=v.size()+count_of_collect-ndata;
            if(extra_size<0)
                extra_size=0;
            pSort::dataType *send_extra_data=new pSort::dataType[extra_size];
            for(int i1=0;i1<extra_size;i1++)
            {
                send_extra_data[i1]=collect[i1+count_of_collect- extra_size];
                
            }
            pSort::dataType *collect_to_keep=new pSort::dataType[count_of_collect- extra_size];
            for(int i1=0;i1<count_of_collect- extra_size;i1++)
            {
                collect_to_keep[i1]=collect[i1];

            }
            MPI_Send( send_extra_data, extra_size, newtype, (my_rank+1)%numProcess, 10, MPI_COMM_WORLD);            //send extra data to next processor

            data_duplicate=merge_arrays(numbers_not_in_map_array, collect_to_keep, v.size(), count_of_collect- extra_size);
            ndata_new=v.size()+ count_of_collect- extra_size;
            
            
            for(int i1=0;i1<count_of_collect- extra_size;i1++)
            {
                v.push_back(collect_to_keep[i1]);
            }

            delete[] rcounts;
            delete[] numbers_to_be_sent;
            delete[] displs;
            delete[] collect;
            delete[] send_extra_data;
            delete[] collect_to_keep;
                        
        }
        else
        {

            int *rcounts=new int[numProcess];
            MPI_Allgather(&size_of_itr_th_digit_numbers, 1, MPI_INT, rcounts, 1 , MPI_INT, MPI_COMM_WORLD);//every processor knows the elements on each processor in the current iteration of the buckets
            int count_of_collect=0;
            for(int k=0;k<numProcess;k++)
            {
                count_of_collect+=rcounts[k];
            }
            
            pSort::dataType *numbers_to_be_sent=new pSort::dataType[size_of_itr_th_digit_numbers];
            for(int p=0;p<size_of_itr_th_digit_numbers;p++)
            {
                numbers_to_be_sent[p]=((*itr).second)[p];
            }
            MPI_Gatherv(numbers_to_be_sent, size_of_itr_th_digit_numbers, newtype, NULL, NULL, NULL, newtype, 0, MPI_COMM_WORLD);//collect this iteration's bucket on 0
            int recieving_rank=(my_rank-1);
                

            MPI_Status statusx;
            MPI_Probe(recieving_rank,10,MPI_COMM_WORLD,&statusx);
            int recv_count;
            MPI_Get_count(&statusx, newtype, &recv_count);
            pSort::dataType *numbertoreceive=new pSort::dataType[recv_count];
         
            MPI_Recv( numbertoreceive, recv_count, newtype, recieving_rank,  MPI_ANY_TAG, MPI_COMM_WORLD, &status);     //recieve extra elements from next processor
         
            int extra_size=v.size()+recv_count-ndata;
            if(extra_size<0)
                extra_size=0;
            pSort::dataType *send_extra_data=new pSort::dataType[extra_size];       

            for(int i1=0;i1<extra_size;i1++)
            {
                send_extra_data[i1]=numbertoreceive[i1+recv_count-extra_size];
            }
            pSort::dataType *collect_to_keep=new pSort::dataType[recv_count- extra_size];
            for(int i1=0;i1<recv_count- extra_size;i1++)
            {
                collect_to_keep[i1]=numbertoreceive[i1];
                // v.push_back(collect_to_keep[i1]);
            
            }
            if(my_rank!=numProcess-1){
                MPI_Send( send_extra_data, extra_size, newtype, (my_rank+1)%numProcess, 10, MPI_COMM_WORLD);            //send extra elements to next processor if it is not the last processor
                
        	}

        	data_duplicate=merge_arrays(numbers_not_in_map_array, collect_to_keep, v.size(), recv_count- extra_size);
                ndata_new=recv_count- extra_size + v.size();
            for(int i1=0;i1<recv_count-extra_size;i1++)
            {
                v.push_back(collect_to_keep[i1]);

            }
			
            delete[] rcounts;
            delete[] numbers_to_be_sent;
            delete[] numbertoreceive;
            delete[] send_extra_data;
            delete[] collect_to_keep;            
        }

    }
    v.clear(); 
    m.clear(); 
    numbers_not_in_map.clear();
    for(int i1=0;i1<ndata;i1++)
    {
        temp[i1]=data_duplicate[i1];
    }

    delete[] data_duplicate;
    
}
////////////////////////////////////////////////////////////////////////////////////////RADIX SORT///////////////////////////////////////////////


void pSort::sort(dataType *data, int ndata, SortType sorter)
{
	if(sorter== BEST)
		parallel_merge(data, ndata);  
	else if (sorter== QUICK)
		quick_sort(data, ndata, MPI_COMM_WORLD);
	else if(sorter == MERGE)
		parallel_merge(data, ndata);
	else 
    radix_sort(data, ndata);  
}
    