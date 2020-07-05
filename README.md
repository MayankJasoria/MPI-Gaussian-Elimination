# Gaussian Elimination using MPI
This program solves a system of linear equations using the Gaussian Elimination technique. The normal sequential method, as applied on an augmented matrix, consists of:
 * Formation of an Upper Triangular Matrix
 * Back-substitution from bottom to top
 
In the parallel implementation, since the data is split into multiple processes, the technique becomes slightly more involved, due to the communication design and ordering of the data. The technique now additionally involves data partitioning before the sequential steps, and gathering the results at the end, and is coupled with ultiple communication stpes in each phase to achieve the desired result.

## Data Parititioning
The data is divided in such a manner that each processes gets approximately the same number of rows of the matrix. For simplification of the communication, the matrix is maintained as a flattened 1-D array. The partition is done in a cyclic manner, such that each process receives rows such that the indices of the rows satisfy the equation `row_index % number_of_processes = process_rank`. That is to say, if 10 rows are to be divided among 3 processes (with both indexing and rank starting from 0), then the division would be as follows:

| Process Rank 	| Row indices 	|
|--------------	|-------------	|
| 0            	| 0,3,6,9     	|
| 1            	| 1,4,7       	|
| 2            	| 2,5,8       	|

A direct partitioning scheme could also have been used, by allotting rows to processes in sequentially ordered chunks from the top, leading to a partition as follows:

| Process Rank 	| Row indices 	|
|--------------	|-------------	|
| 0            	| 0,1,2,3     	|
| 1            	| 4,5,6       	|
| 2            	| 7,8,9       	|

However, such direct partitioning increases the overall waiting time, since processes 1 and 2 would be idle till the task of process 1 is completed. Cyclic partitioning overcomes this issue to some extent by assigning almost equal workloads to all processes.

## Formation of an Upper Triangular Matrix
In this step, the input matrix is converted into an upper triangular matrix, where the diagonal of the matrix is 1 and all elements below the diagonal are 0. To acheive this, the following steps are performed in a pipelined manner:
1. **Operations with rows having lower index**
 * **Rearrange the order of elements of current row according to the indices of the pivots of previous rows**: Simply swap the element at the position below the diagonal of the previous row with the element that appears at the original index of the pivot of previous row.
 * **Make element below diagonal of previous rows 0:** This is achieved by subtracting the previous row multiplied by the element just below diagonal of previous row, with the current row.
2. **Finding the pivot element:** The pivot element of any row is the largest number (by magnitude) in that row. This is placed as the diagonal element of that row.
3. **Divide row by pivot:** This makes the value of the diagonal 1.
4. **Pass on row to repeat step 1 in rows having higher index**

Assuming that the processes are arranged curcularly in sorted orer of their rank, then each processes receives rows for operations only from the process that appears immediately before the current process, and sends data only to the process immediately after the current process. Thus, each row is communicated only once to every process.

## Back-substitution
At this stage, the reverse idea is applied, where the solution of the rows below is propagated upward, multiplied with the value at that index in the current row, then subtracted from current row. This process is continued till the solution of all variables is found.

## Gathering the results
Finally, the results of all processes are gathered into a single process. This process then reshuffles the results to rearrange them from the cyclically-ordered chunks to the sequential manner, then writes the final result as output.
