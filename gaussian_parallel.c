#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Subtracts an earlier row of the matrix with all later rows that are present in the current process
 * @param recv_id       The ID of the process from which the row was received
 * @param id            ID of current process
 * @param num_eq        The number of equations in the system
 * @param proc_rows     The chunk of the coefficient matrix contained within current process
 * @param proc_vals     The chink of the values vector contained within the current process
 * @param curr          The row to be processed next by the current process
 * @param recvd_row     The complete row that was received
 * @param rows_per_proc The total number of rows present within current process
 * @param num_proc      The total number of processes
 */
void perform_elimination(int recv_id, int id, int num_eq, double* proc_rows, double* proc_vals, int curr, double* recvd_row, int rows_per_proc, int num_proc) {
    
    for(int i = curr / num_proc; i < rows_per_proc; i++) {
        int piv = recvd_row[num_eq];
        int tmp = proc_rows[(i * num_eq) + recv_id];
        proc_rows[(i*num_eq) + recv_id] = proc_rows[(i*num_eq) + piv];
        proc_rows[(i*num_eq) + piv] = tmp;
        for(int j = 0; j < num_eq; j++) {
            proc_rows[(i * num_eq) + j] -= (recvd_row[j] * proc_rows[(i*num_eq) + recv_id]);
        }
        proc_vals[i] -= recvd_row[num_eq + 1] * proc_rows[(i*num_eq) + recv_id];
    }
}

/**
 * Computes the pivot element of a given row
 * @param curr      The row whose pivot element is to be computed
 * @param proc_rows The number of rows assigned to current process
 * @param num_proc  The total number of processes
 * @param num_eq    The number of equations in the given system
 *
 * @return The index of the pivot element
 */
int compute_pivot(int curr, int num_proc, int num_eq, double* proc_rows) {
    int mx = -1;
    int pivot;
    int row_id = curr / num_proc;
    for(int i = curr; i < num_eq; i++) {
        int val = proc_rows[(row_id * num_eq) + i];
        if(abs(val) > mx) {
            mx = abs(val);
            pivot = i;
        }
    }
    return pivot;
}

/**
 * Divides a complete row by the pivot element
 * @param id        The ID of the current process
 * @param curr      The row in which division is to be performed
 * @param proc_rows The number of rows assigned to current process
 * @param pivot     Index of the pivot element
 * @param num_proc  The total number of processes
 * @param num_eq    The number of equations in the given systemc_vals The chink of the values vector contained within the current process
 */
void perform_division(int id, int curr, double* proc_rows, int pivot, int num_proc, int num_eq, int rows_per_proc, double * proc_vals) {
    
    int row_id = curr / num_proc;
    int tmp = proc_rows[(row_id * num_eq) + pivot]; 
    proc_rows[(row_id*num_eq)+ pivot] = proc_rows[(row_id*num_eq) + curr];
    proc_rows[(row_id*num_eq) + curr] = tmp;
    
    for(int i = curr; i < num_eq; i++) {
        proc_rows[(row_id*num_eq) + i] /= proc_rows[(row_id * num_eq) + curr];
    }
    proc_vals[row_id] /= proc_rows[(row_id*num_eq) + curr];
}


int main(int argc, char** argv) {
    int num_proc;
    int id;
    int num_eq = 0;
    int rows_per_proc;
    double time_taken;
    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    /* Get number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    /* Get current process id */
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    FILE* inputfp = NULL;
    if(id == 0) {
        /* read input from file specified in command line arguments */
        inputfp = fopen(argv[1], "r");
        
        /* read number of equations */
        fscanf(inputfp, "%d", &num_eq);
        
    }
    /*
        Declaring static arrays for equation coefficients.
        NOTE: num_eq = 0 for processes with rank other than 0, so no memory assigned to them,
        The memory is assigned only to the root process.
    */
    //double eq_mat[num_eq][num_eq];
    double eq_mat[num_eq * num_eq];
    double val_mat[num_eq];
    int divs[num_proc];
    int displs[num_proc];
    MPI_Bcast(&num_eq, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int rpp = num_eq / num_proc;                             
        
    /* Number of processes with rows = rpp (np0), rows = rpp + 1(np1) */
    int np0 = num_proc - (num_eq % num_proc);
    int np1 = num_proc - np0;
    if(id >= np1)
        rows_per_proc = rpp;
    else
        rows_per_proc = rpp + 1;
    
    if(id == 0) {    
        
        for(int i = 0; i < num_proc; i++) {
            divs[i] = 0;
        }
        
        /* reading coefficient matrix */
        for(int i = 0; i < num_eq; i++) {
            
            int assigned_proc = i % num_proc;
            int rows_before = 0;
            if(assigned_proc <= np1) {
                rows_before = (rpp + 1) * (assigned_proc);
            } else {
                rows_before = (rpp + 1) * (np1);
                rows_before += (assigned_proc - np1)*(rpp);
            }
            int eff_row = rows_before + (divs[assigned_proc]++);
            if(i >= num_eq)
                continue;
            for(int j = 0; j < num_eq; j++) {
                fscanf(inputfp, "%lf", &eq_mat[eff_row*num_eq + j]);
            }
            fscanf(inputfp, "%lf", &val_mat[eff_row]);
        }
        /* calculating the displacement vector */
        displs[0] = 0;
        for(int i = 1; i < num_proc; i++) {
            divs[i-1] *= num_eq;
            displs[i] = displs[i - 1] + divs[i - 1];
            /// divs[i-1] *= num_eq;
        } 
        divs[num_proc-1] *= num_eq;
        
        fclose(inputfp);
    }
    /* input read, synchronize at this point */
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* chunk of coefficient matrix per process */
    double proc_rows[num_eq * rows_per_proc];
    MPI_Scatterv(eq_mat, divs, displs, MPI_DOUBLE, proc_rows, rows_per_proc * num_eq, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    
    /* chunk of values vector per process */
    double proc_vals[rows_per_proc];
    for(int i = 0; i < num_proc; i++) {
        divs[i] /= num_eq;
        displs[i] /= num_eq;
    }
    MPI_Scatterv(val_mat, divs, displs, MPI_DOUBLE, proc_vals, rows_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* DEBUG */
    char fname[11];
    strcpy(fname, "debug");
    fname[5] = ('0' + (char)id);
    fname[6] = '.';
    fname[7] = 't';
    fname[8] = 'x';
    fname[9] = 't';
    fname[10] = '\0';
    FILE* debugfp = fopen(fname, "w");
    for(int i = 0; i < rows_per_proc; i++) {
        for(int j = 0; j < num_eq; j++) {
            fprintf(debugfp, "%lf ", proc_rows[(i * num_eq) + j]);
        }
      fprintf(debugfp, "%lf\n", proc_vals[i]);
    //    for(int j = 0; j < num_eq; j++) {
    //        fprintf(debugfp, "%lf ", eq_mat[i][j]);
    //    }
        // fprintf(debugfp, "\n");
    }
    fclose(debugfp);
    MPI_Barrier(MPI_COMM_WORLD);
    /* DEBUG ENDS */
    
    /* The process which receives data from the current proc */
    int next_proc = (id + 1) % num_proc;

    /* The process which sends the data to the current process */
    int prev_proc = (id - 1 + num_proc) % num_proc;

    /* The row which is being processes currently */
    int curr = id;
    int prev_curr = -1;
    int piv; 
    double recvd_row[num_eq + 2]; // recvd_row[num_eq] = pivot, recvd_row[num_eq + 1] = value
    MPI_Status st;  

    for(int cnt = 0; cnt < rows_per_proc; cnt++) {
        for(int i = prev_curr + 1; i < curr; i++) {
            // pipelined send receive of all rows between previoud row and current row
            MPI_Recv(recvd_row, num_eq + 2, MPI_DOUBLE, prev_proc, i, MPI_COMM_WORLD, &st);
            MPI_Send(recvd_row, num_eq + 2, MPI_DOUBLE, next_proc, i, MPI_COMM_WORLD);
            perform_elimination(i, id, num_eq, proc_rows, proc_vals, curr, recvd_row, rows_per_proc, num_proc);  // eliminates from current row till last row in process
        }
        piv = compute_pivot(curr, num_proc, num_eq, proc_rows); 
        perform_division(id, curr, proc_rows, piv, num_proc, num_eq, rows_per_proc, proc_vals);
        double send_buf[num_eq + 2];
        for(int j = 0; j < num_eq; j++) {
            send_buf[j] = proc_rows[(cnt * num_eq) + j];
        }
        send_buf[num_eq] = piv;
        send_buf[num_eq + 1] = proc_vals[cnt];
        
        if(curr < num_eq) {
            
            /* We do not want to send last row */
            MPI_Send(send_buf, num_eq + 2, MPI_DOUBLE, next_proc, curr, MPI_COMM_WORLD);
        }
        
        prev_curr = curr;
        curr += num_proc;
        
        if(curr < num_eq) {
            perform_elimination(prev_curr, id, num_eq, proc_rows, proc_vals, curr, send_buf, rows_per_proc, num_proc);
        }
    }
        
    MPI_Finalize();
    return 0;
}