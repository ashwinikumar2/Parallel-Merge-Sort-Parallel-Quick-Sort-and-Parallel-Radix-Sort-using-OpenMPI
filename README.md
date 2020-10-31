# Parallel-Merge-Sort-Parallel-Quick-Sort-and-Parallel-Radix-Sort-using-OpenMPI
**Implementation  of parallel versions of merge sort, quick sort and radix sort**.
Among these, merge sort show the best performance on 64 processors.

In this, we have implemented parallel versions of quick sort, merge sort,
and radix sort. Input is given on each processor and the
output is the soted version of data on all processors such that elements at
processor i are less than elements at processor j for i<j. Let us
go through every type of sort one by one:

**RADIX SORT:**

We start from the right most digit of every data's key in all the processors.
For each position of digits, 10 buckets ranging from (0-9) are made on every
processor. All elements of data on a processor get in the bucket
of their digit on which iteration is going on. The important thing to note here is
that negative numbers are also treated as positive number
while allocating them into a bucket but there sign has not been changed in the
whole program.

The buckets are then transferred to processor with rank 0. Processor 0 decides
whether to keep the elements of buckets within itself or transfer
it to another processor. This decision is based upon the size that gets increased
as compared to the original size of processor 0. If the size
does not exceed the original size of processor 0 then it keeps all the collected
elements wihtin itself. Else send to it to the next processor
available.

After successfully allocating all the collected elements of buckets have been
allocated to either 0th processor or many processors, we search
for the numbers having number of digits greater than the iteration's digit number.
These elemnents are all sent to the last processor and deleted
from the processors. The last processors appends them at the last and moves the
elements contained already in it to previous processor. This way
all elements with less digit than the iteration's digit are filled from 0th
processor to some p'th processor and then numbers with digits more
than the iteration's digit are filled from p'th processor to the last processor.
Once all the elements are sorted based on this classification, all elements
magnitude wise are sorted on each processor which treats negative and
positive numbers to be the same.

So, again two buckets are made on each processor. One of positive numbers and other
of negative numbers. Likewise before, all buckets of all processors
are collected on 0th processor and then redistributed to other processors based on
the same criteria as stated above.
Finally, we get that all elements are sorted globally and that every element of
i'th processor is less than j'th processor if i<j.

**QUICK SORT:**

In this we first chose a pivot to be the mid element of last processor and then
broadcast it so that every other processors know the value of pivot.
Initially the communicator is MPI_COMM_WORLD.

Then every processor divides its elements into two categories: the ones lower than
pivot and the other greater than the pivot.

Every processor has a partner processor. Simply we can say that processor 0 is a
partner processor of the last processor and so on.

If the processor belongs to upper half of the processors in the current
communicator then it sends its elements having key less than the pivot's key
value to the partner processor in the lower half of processors in this
communicator.

Else it sends the elements greater than pivot to partner processor in the upper
half of the processor's list.

This way all upper half processors contain the values greater than pivot and lower
half processors contain the elements less than pivot.

Now, this communicator is split into two and a recursive call is made to sort in
the same way. Note that the recursion is on the number of processors.

After the end of recursion we get all sorted elements but it is not necessary that
every processor now contains the same number of elements as it
originally contained.

For this, we send extra elements from processor 0 to the last processor. The last
processor sends extra elements to the second last element and so on.

Now, every processor contains the same number of elements it originally had and
every element of i'th processor is less than j'th processor if i<j.

**MERGE SORT:**

For this we start by a sequential sort locally in each processor by the usual
sequential merge sort algorithm.

The algorithm is that we divide an array into two halves adn then recursively call
merge sort function on both half of the array.

Then we perform the merge operation which only on both the arrays.
This way all elements are sorted locally.

Now we merge all the elements by gathering them all on 0th processor and then
distributing it to every other processor by using MPI_Scatterv and based
on the original size of every processor.

So, now all elemenets are sorted such that every element of i'th processor is less
than j'th processor if i<j.

The major decisions made were that:

(1) Sequential algorithms are much more time taking in comparison to parallel
algorithms. Still all parallel algorithms use some kind of serialization because it
was not possible to sort parallely when there is only one processor available such
as in the case of quick sort. Parallel algorithms proved to be much faster as shown
in the tables below.

(2) In Radix sort, at first I used many send-receive matching pairs after each
bucket sort so that the total size on each processor must remain the same as
originally it had. This led to a very high time taken for program to execute. So, I
did not make any transfer among all the processors after each bucket. This reduced
the time taken by many folds. In the later case, the program was very fast.

(3) In quick sort, the data was sorted but the sizes were not equal to the original
size every processor had. So, to do this I had to transfer data in a chain like
fashion from 0th processor to last processor such that every processor now has the
necessary size. But, there was still a case where some processor with less data
than it should have. So, a there was another transfer of data in a backward chain.
Otherwise, the data could be left out at some processor.

(4) In merge function used for merge sort, a duplicate array was passed. Since, it
was giving a segmentation fault when we do not pass a duplicate array. It was
previously creating two arrays of large size which was taking up extra memory. When
a duplicate array was passed instead of creating two new arrays, it was working
properly and with good efficiency.

(5) In Quick sort, when number of processors were odd then there was a processor in
mid which was left out because there was no partner process for it. So, it should
perform an extra send recieve to the next processor exceptionally.
SCALABILITY: All sort versions show scalability(almost linear in case of merge &
quick sort while it is not in radix).

Experiments:

For number of elements(N)= 2^20(approx. 10,00,000) on each processor

**MERGE SORT:**\
|NUMBER OF PROCESSORS| Time Taken|
|-----------------------|-----------------|
|1| 1|
|4 |2|
|16| 8|
|32| 17|
|64 |38|

**QUICK SORT:**\

|NUMBER OF PROCESSORS| Time Taken|
|----|-----|
|1| 1|
|4 |2|
|16| 6|
|32| 13|
|64| 110|

**RADIX SORT:**\

|NUMBER OF PROCESSORS |Time Taken|
|------|------|
|1 |4|
|4| 7|
|16| 67|
|32 |180|
|64| 248|

**BEST:** Merge Sort proved to be fast as compared to other two sorts. So , merge sort
is kept in BEST.

**HOW TO RUN:**\
make\
make run

**Note:** make run will run on 6 processors else we can use\
mpirun -np x  ./output\
where x is number of processors
