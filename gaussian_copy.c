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
 * @param var_perm  
 */
void perform_elimination(int recv_id, int id, int num_eq, double* proc_rows, double* proc_vals, int curr, double* recvd_row, int rows_per_proc, int num_proc, int * var_perm) {
	/* update indices of each variable */
	int piv_idx = (int) recvd_row[num_eq];
	int tmp = var_perm[piv_idx];
	var_perm[piv_idx] = var_perm[recv_id];
	var_perm[recv_id] = tmp;
	
	for(int i = 0; i < rows_per_proc; i++) {
		if(i*num_proc + id == recv_id) {
			continue;
		}
		int piv = (int) recvd_row[num_eq];
		double tmp = proc_rows[(i * num_eq) + recv_id];
		proc_rows[(i*num_eq) + recv_id] = proc_rows[(i*num_eq) + piv];
		proc_rows[(i*num_eq) + piv] = tmp;
	  
		if(i >= curr / num_proc) {
			double piv_val = proc_rows[(i*num_eq) + recv_id];
			for(int j = 0; j < num_eq; j++) {
				proc_rows[(i * num_eq) + j] -= (recvd_row[j] * piv_val);
			}
			proc_vals[i] -= recvd_row[num_eq + 1] * piv_val;
		}
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
	double mx = -1;
	int pivot;
	int row_id = curr / num_proc;

	for(int i = curr; i < num_eq; i++) {
		double val = proc_rows[(row_id * num_eq) + i];
		if(val < 0)
			val *= -1;
		if(val > mx) {
			mx = val;
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
	/* swapping pivot element with element at pivot index */
	int row_id = curr / num_proc;
	double tmp = proc_rows[(row_id * num_eq) + pivot]; 
	proc_rows[(row_id*num_eq)+ pivot] = proc_rows[(row_id*num_eq) + curr];
	proc_rows[(row_id*num_eq) + curr] = tmp;
	
	double piv_val = proc_rows[(row_id * num_eq) + curr];
	for(int i = curr; i < num_eq; i++) {
		proc_rows[(row_id*num_eq) + i] /= piv_val;
	}
	proc_vals[row_id] /= piv_val;
}

void print_arr(double* arr, int size, FILE* fptr) {
	for(int i = 0; i < size; i++) {
		fprintf(fptr, "%lf ", arr[i]);
	}
	fprintf(fptr, "\n");
}


int main(int argc, char** argv) {
	int num_proc; // number of processors
	int id; // process rank
	int num_eq = 0; // number of equations
	int rows_per_proc; // rows allocated to the current process
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
	
	double eq_mat[num_eq * num_eq]; //double eq_mat[num_eq][num_eq]; flattened to 1d
	double val_mat[num_eq]; // Rhs values (y)
	
	/* for flattening */
	int divs[num_proc]; // Number of coefficients per proc
	int displs[num_proc]; // Begining index of coefficients of each process

	time_taken = MPI_Wtime(); // recording start time

	MPI_Bcast(&num_eq, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcasted number of equations to all processes
	int rpp = num_eq / num_proc; // temporary variables

	 
	int np0 = num_proc - (num_eq % num_proc); // Number of processes with rows = rpp (np0),
	int np1 = num_proc - np0; // number of processes with rows = rpp + 1(np1)

	/* assigning rows_per_proc based on the rank */
	if(id >= np1)
		rows_per_proc = rpp;
	else
		rows_per_proc = rpp + 1;
	
	/* The permutation of the columns of the current row */
	int var_perm[num_eq];
	for(int i = 0; i < num_eq; i++) {
		var_perm[i] = i;
	}
	
	/* FILE input */
	if(id == 0) {    
		
		for(int i = 0; i < num_proc; i++) {
			divs[i] = 0;
		}
		
		/* iterating over all the rows */
		for(int i = 0; i < num_eq; i++) {
			
			int assigned_proc = i % num_proc; // process assigned for the current row to be read
			
			/* calculating the actual position (eff_row) of the row in the partitioned view of rows */
			int rows_before = 0; 
			if(assigned_proc <= np1) {
				rows_before = (rpp + 1) * (assigned_proc);
			} else {
				rows_before = (rpp + 1) * (np1);
				rows_before += (assigned_proc - np1)*(rpp);
			}
			/* adding the position of the row in the assigned processor */
			int eff_row = rows_before + (divs[assigned_proc]++);
			
			/* reading the corresponding coefficients (num_eq + 1) */
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
		}
		divs[num_proc-1] *= num_eq;
		
		fclose(inputfp);
	}

	/* input read, synchronize at this point */
	MPI_Barrier(MPI_COMM_WORLD);
	
	/* chunk of coefficient matrix per process */
	double proc_rows[num_eq * rows_per_proc]; // data for the current process
	MPI_Scatterv(eq_mat, divs, displs, MPI_DOUBLE, proc_rows, rows_per_proc * num_eq, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// MPI_Barrier(MPI_COMM_WORLD); /* TODO: remove */
	
	/* chunk of values vector per process */
	double proc_vals[rows_per_proc]; // y values for the current process
	for(int i = 0; i < num_proc; i++) {
		divs[i] /= num_eq;
		displs[i] /= num_eq;
	}
	MPI_Scatterv(val_mat, divs, displs, MPI_DOUBLE, proc_vals, rows_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// MPI_Barrier(MPI_COMM_WORLD); /* TODO: remove */
	
	/* The process which receives data from the current proc */
	int next_proc = (id + 1) % num_proc;

	/* The process which sends the data to the current process */
	int prev_proc = (id - 1 + num_proc) % num_proc;

	
	int curr = id; // The row which is being processes currently
	int prev_curr = -1; // Id of the previous row
	int piv; // Id of the pivot element in the current row
	double recvd_row[num_eq + 2]; // recvd_row[num_eq] = pivot, recvd_row[num_eq + 1] = y-value
	MPI_Status st; // temoporary status variable

	/*************** Upper Triangualar Matrix Generation Phase ***************/

	char fname[15];
	sprintf(fname, "out%d.txt", id);
	FILE* fptr = fopen(fname, "w");

	/*   
		Every process herein will iterate on the rows assigned to it in the sorted order,
		performing relevant communications in a pipelined manner.
	*/
	fprintf(fptr, "Upper Triangular Matrix:\n");
	for(int cnt = 0; cnt < rows_per_proc; cnt++) {
		fprintf(fptr, "cnt = %d\n", cnt);
		/* Iterating over all the rows before the current row and previously processed row, to perform elimination corresponding to that row */
		for(int i = prev_curr + 1; i < curr; i++) {
			// pipelined send receive of all rows between previous row and current row
			MPI_Recv(recvd_row, num_eq + 2, MPI_DOUBLE, prev_proc, i, MPI_COMM_WORLD, &st);
			fprintf(fptr, "i = %d, curr = %d, prev_curr = %d, Array received from %d: ", i, curr, prev_curr, prev_proc);
			print_arr(recvd_row, num_eq + 2, fptr);
			/* 
				The recieved row has traversed (curr - i) number of times, so it it has traversed num_proc or more times, 
				it has been delivered to all processes 
				TODO: check this
			*/
			if(curr - i < num_proc - 1) {
				MPI_Send(recvd_row, num_eq + 2, MPI_DOUBLE, next_proc, i, MPI_COMM_WORLD);
				fprintf(fptr, "i = %d, curr = %d, prev_curr = %d, Array sent to %d: ", i, curr, prev_curr, next_proc);
				print_arr(recvd_row, num_eq + 2, fptr);
			}

			/* preforming elimination corresponding to the recieved row */
			perform_elimination(i, id, num_eq, proc_rows, proc_vals, curr, recvd_row, rows_per_proc, num_proc, var_perm);

			fprintf(fptr, "i = %d, Post Elimination:\n", i);
			for(int m = 0; m < rows_per_proc; m++) {
				print_arr(proc_rows + (m * num_eq), num_eq, fptr);
			}
			fprintf(fptr, "Proc Vals:\n");
			for(int m = 0; m < rows_per_proc; m++) {
				fprintf(fptr, "%lf ", proc_vals[m]);
			}
			fprintf(fptr, "\n");
		}
		piv = compute_pivot(curr, num_proc, num_eq, proc_rows); 
		int tmp = var_perm[piv];
		var_perm[piv] = var_perm[curr];
		var_perm[curr] = tmp;
		fprintf(fptr, "\nComputed pivot: %d\n", piv);

		perform_division(id, curr, proc_rows, piv, num_proc, num_eq, rows_per_proc, proc_vals);
		fprintf(fptr, "Post division:\n");
		for(int i = 0; i < rows_per_proc; i++) {
			print_arr(proc_rows + (i * num_eq), num_eq, fptr);
		}
		fprintf(fptr, "Proc Vals:\n");
		for(int i = 0; i < rows_per_proc; i++) {
			fprintf(fptr, "%lf ", proc_vals[i]);
		}
		fprintf(fptr, "\n");
		
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
		
		perform_elimination(prev_curr, id, num_eq, proc_rows, proc_vals, curr, send_buf, rows_per_proc, num_proc, var_perm);
		fprintf(fptr, "post elimination with current (%d):\nProc Rows:", curr);
		for(int i = 0; i < rows_per_proc; i++) {
			print_arr(proc_rows + (i * num_eq), num_eq, fptr);
		}
		fprintf(fptr, "Proc Vals:\n");
		for(int i = 0; i < rows_per_proc; i++) {
			fprintf(fptr, "%lf ", proc_vals[i]);
		}
		
		fprintf(fptr, "\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	fprintf(fptr, "\n\nPost UTM\n");
	for(int i = 0; i < rows_per_proc; i++) {
		for (int j = 0; j < num_eq; j++) {
			fprintf(fptr, "%lf ", proc_rows[i*num_eq + j]);
		}
		fprintf(fptr, "%lf\n", proc_vals[i]);
	}
	
	fprintf(fptr, "Var perm\n");
	for (int j = 0; j < num_eq; j++) {
		fprintf(fptr, "%d ", var_perm[j]);
	}

	/*************** Back Substitution Phase ***************/
	
	int lst = id + (rows_per_proc - 1) * num_proc; // actual id of the last row of the current process
	int prev_lst = num_eq;
	double x_val;
	int count = 1;

	fprintf(fptr, "\nBack Substitution Steps\n");
	for(int cnt = rows_per_proc - 1; cnt >= 0; cnt--) {
		fprintf(fptr, "cnt: %d\n", cnt);
		for(int i = prev_lst - 1; i > lst; i--) {
			MPI_Recv(&x_val, count, MPI_DOUBLE, next_proc, i + num_eq, MPI_COMM_WORLD, &st);
			fprintf(fptr, "i = %d, prev_lst = %d, lst = %d, Received %lf from %d\n", i, prev_lst, lst, x_val, next_proc);
			if(i - lst < num_proc - 1) {
				MPI_Send(&x_val, count, MPI_DOUBLE, prev_proc, i + num_eq, MPI_COMM_WORLD);
				fprintf(fptr, "i = %d, prev_lst = %d, lst = %d, Sent %lf to %d\n", i, prev_lst, lst, x_val, prev_proc);
			}
			for(int k = cnt; k >= 0; k--) {
				proc_vals[k] -= proc_rows[(k * num_eq) + i] * x_val;
			}
			fprintf(fptr, "i = %d, proc_vals = ", i);
			print_arr(proc_vals, rows_per_proc, fptr);
		}
		double ans = proc_vals[cnt];
		for(int k = cnt - 1; k >= 0; k--) {
			proc_vals[k] -= proc_rows[(k * num_eq) + lst] * ans;
		}
		fprintf(fptr, "cnt = %d, proc_vals = ", cnt);
		print_arr(proc_vals, rows_per_proc, fptr);

		if(lst > 0) {
			/* We don't want to send first pivot */
			MPI_Send(&ans, count, MPI_DOUBLE, prev_proc, lst + num_eq, MPI_COMM_WORLD);
		}
		prev_lst = lst;
		lst -= num_proc;
	}

	double res[num_eq];

	MPI_Gatherv(proc_vals, rows_per_proc, MPI_DOUBLE, res, divs, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	int proc_with_last_row = (num_eq - 1) % num_proc;

	if(proc_with_last_row != 0 && id == proc_with_last_row) {
		MPI_Send(var_perm, num_eq, MPI_INT, 0, num_eq * 4 + 4, MPI_COMM_WORLD);
	}

	if(id == 0) {
		if(proc_with_last_row != 0)
			MPI_Recv(var_perm, num_eq, MPI_INT, proc_with_last_row, num_eq * 4 + 4, MPI_COMM_WORLD, &st);
		
		int k = 0;
		double solution[num_eq];
		
		/* arrange data in the sequential format */
		for(int i = 0; i < num_proc; i++) {
			for(int j = i; j < num_eq; j += num_proc) {
				solution[j] = res[k++];
			}
		}
		
		/* swap data according to the variable permutations */
		// for(int i = 0; i < num_eq; i++) {
		// 	double temp = solution[i];
		// 	solution[i] = solution[var_perm[i]];
		// 	solution[var_perm[i]] = temp;
		// }

		double sol[num_eq];
		for(int i = 0; i < num_eq; i++) 
			sol[var_perm[i]] = solution[i];
		
		/* print output */
		FILE* outfile = fopen("output.txt", "w");
		fprintf(outfile, "The solution of the system of equations is \n{");
		for(int i = 0; i < num_eq; i++) {
			fprintf(outfile, "%lf ", sol[i]);
		}
		fprintf(outfile, "}\n");
		fclose(outfile);
	}
	  
	/* Ensuring that MPI_Finalize() is not called prematurely by any process */
	MPI_Barrier(MPI_COMM_WORLD);  

	time_taken = MPI_Wtime() - time_taken;

	if(id == 0) {
		printf("Time taken for execution is %lf\n", time_taken);
	}
	 
	MPI_Finalize();
	return 0;
}
