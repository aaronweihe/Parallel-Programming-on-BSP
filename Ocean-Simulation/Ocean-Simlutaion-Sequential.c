//
//  main.c
//  Ocean-Prgramming
//
//  Created by Wei He on 3/20/12.
//  Copyright (c) 2012 Suffolk University. All rights reserved.
//

#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>

#define N 9999
#define TOL 0.001

int main(int argc, const char * argv[])
{
    int i, j, done, counter;
    float diff =0.0, temp = 0.0;
    
    struct timeval start, end;
    float elapsed_time = 0.0;
    gettimeofday(&start, NULL);
    srand(1452764);
    //Matrix initilization
    float **matrix = (float**)calloc(N+2, sizeof(float*));
    for (i=0; i<N+2; i++) {
        matrix[i] = (float*)calloc(N+2, sizeof(float));
    }
    for (i=0; i<N+2; i++) {
        for (j=0; j<N+2; j++) {
            matrix[i][j] = (float)rand()/(float)RAND_MAX;
            //printf("row %d, coloum %d, element: %f", i, j, matrix[i][j]);
        }
    }
    
    //Ocean Simulation
    done = counter = 0;
    //Computation
    while (!done) {
        //printf("Now, %d th roud:\n", ++counter);
        diff=0.0;
        
        for (i=1; i<N+1; i++) {
            for (j=1; j<N+1; j++) {
                temp = matrix[i][j];
                matrix[i][j] = 0.2 * (matrix[i][j] 
                                      + matrix[i][j-1] + matrix[i-1][j] 
                                      + matrix[i][j+1] + matrix[i+1][j]);
                diff += fabs(matrix[i][j] - temp);
            }
        }
        temp = diff/(float)(N*N);
        //printf("temp: %f\n", temp);
        if (temp<TOL) {
            done=1;
        }
        //printf("Current difference is %lf .\nAnd convergence is %lf \n", diff, temp);
        //printf("*****************************\n");
    }
    gettimeofday(&end, NULL);
    elapsed_time = (end.tv_sec - start.tv_sec)
                        + (end.tv_usec - start.tv_usec) / 1000000000.0;
    printf("Total run time is %.3f seconds.\n", elapsed_time);
    
    //Deallocation
    for (i=0; i<N+2; i++) {
        free(matrix[i]);
    }
    free(matrix);
    
    return 0;
}

