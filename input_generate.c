#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv) {
    if(argc < 3) {
        fprintf(stderr, "Expected 2 arguments, found none\nUsage: %s <no_equations> <output_file>\n", argv[0]);
        return 0;
    }

    int num_eq = atoi(argv[1]);

    FILE* fptr = fopen(argv[2], "w");

    fprintf(fptr, "%d\n", num_eq);

    srand(time(0));

    for(int j = 0; j < num_eq; j++) {
        for(int i = 0; i <= num_eq; i++) {
            int val = rand() % num_eq + 1;
            fprintf(fptr, "%d", val);
            if(i < num_eq) {
                fprintf(fptr, " ");
            } else {
                fprintf(fptr, "\n");
            }
        }
    }

    fclose(fptr);

    return 0;
}