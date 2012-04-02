//
//  Ocean-Programming.c
//  Ocean-Programming
//
//  Created by Wei He on 3/21/12.
//  Copyright (c) 2012 Suffolk University. All rights reserved.
//

#include "bsp.h"

#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define N 10000
#define TOL 0.001

void parallel_part()
{
    int i, j;    
    srand(1452764);
    
    //Matrix initilization
    float **matrix = (float**)calloc(N+2, sizeof(float*));
    for (i=0; i<N+2; i++) {
        matrix[i] = (float*)calloc(N+2, sizeof(float));
    }
    for (i=0; i<N+2; i++) {
        for (j=0; j<N+2; j++) {
            matrix[i][j] = (float)rand()/(float)RAND_MAX;
            //printf("row %d, coloum %d, element: %f\n", i, j, matrix[i][j]);
        }
    }
    
    //Parallel part
    bsp_begin(bsp_nprocs());
    int pid, x, y, done;
    pid=x=y=done=0;
    int sqroot = (int)(sqrt(bsp_nprocs()));
    int size = (int)(N/sqroot);    //side
    float Ai_jm1, Aim1_j, Ai_jp1, Aip1_j;
    Ai_jm1 = Aim1_j = Ai_jp1 = Aip1_j = 0.0;
    float temp, diff, convergence, total_diff;
    temp = convergence = 0.0;
    float *diffs = (float*)calloc(bsp_nprocs(), sizeof(float));
    int counter= 0;
    
    //(N/sqrt(p)) is an integer assurance
    if ( N%sqroot!=0) {
        bsp_abort("N/sqrt(p) is not an integer.\nProgram Aborted.\n");
    }
    
    //Initiliaze a piece of martix in decomposition
    float **sub_martix = (float**)calloc(size, sizeof(float*));
    for (i=0; i<size; i++) {
        sub_martix[i] = (float*) calloc(size, sizeof(float));
    }
    //Initiliaze borders
    float *upper = (float*)calloc(size, sizeof(float));
    float *lower = (float*)calloc(size, sizeof(float));
    float *left = (float*)calloc(size, sizeof(float));
    float *right = (float*)calloc(size, sizeof(float));
    float *overlap = (float*)calloc(size, sizeof(float));
    
    bsp_push_reg(&diff, sizeof(float));
    bsp_push_reg(upper, size*sizeof(float));
    bsp_push_reg(lower, size*sizeof(float));
    bsp_push_reg(left, size*sizeof(float));
    bsp_push_reg(right, size*sizeof(float));
    
    //Make each matrix and border available globally
    for (i=0; i<size; i++) {
        bsp_push_reg(sub_martix[i], size*sizeof(float));
    }
    bsp_sync();
    /*Processor 0 distributes the data*/
    if (bsp_pid()==0) {
        for (pid = 0; pid<bsp_nprocs(); pid++) {
            //Determine which part of the original matrix
            x = pid/sqroot;
            y = pid%sqroot;
            //Then the processor 0 copy the data to each processor
            for (i=0; i<size; i++) {
                for (j=0; j<size; j++) {
                    sub_martix[i][j] = matrix[x*size+i+1][y*size+j+1];
                }
            }
            if (pid!=0) {
                for (i=0; i<size; i++) {
                    bsp_put(pid, sub_martix[i], sub_martix[i], 0, size*sizeof(float));
                }
            }
        }
    }
    bsp_sync();
    
    if (bsp_pid()==0) {
        for (pid=0; pid<bsp_nprocs(); pid++) {
            x=pid/sqroot;
            x=pid%sqroot;
            
            //if the part is in 1st row
            if (x==0) {
                for (i=0; i<size; i++) {
                    upper[i] = matrix[0][y*size+1+i];
                }
            }
            //if the part is in leftmost column
            if (y==0) {
                for (i=0; i<size; i++) {
                    left[i] = matrix[x*size+1+i][0];
                }
            }
            //if the part is in last row
            if (x==sqroot-1) {
                for (i=0; i<size; i++) {
                    lower[i] = matrix[N+1][y*size+1+i];
                }
            }
            //if the part is in rightmost column
            if (y==1) {
                for (i=0; i<size; i++) {
                    right[i] = matrix[x*size+1+i][N+1];
                }
            }
            
            if (pid!=0) {
                bsp_put(pid, upper, upper, 0, size*sizeof(float));
                bsp_put(pid, lower, lower, 0, size*sizeof(float));
                bsp_put(pid, left, left, 0, size*sizeof(float));
                bsp_put(pid, right, right, 0, size*sizeof(float));
            }
        }
    }
    bsp_sync();
    
    /* Computation */
    while (!done) {
        pid = bsp_pid();
        diff=0.0;
        total_diff=0.0;
        x = pid/sqroot;
        y = pid%sqroot;
        //printf("Now %d th round:", ++counter);
        
        if (x<sqroot-1) {
            for (i=0; i<size; i++) {
                overlap[i] = sub_martix[size-1][i];
            }
            bsp_put(bsp_pid()+sqroot, overlap, upper, 0, size*sizeof(float));
        }
        if (y<sqroot-1) {
            for (i=0;i<size;i++) {
                overlap[i]=sub_martix[i][size-1];
            }
            bsp_put(bsp_pid()+1, overlap, left, 0, size*sizeof(float));
        }
        if (x>0) {
            for (i=0; i<size; i++) {
                overlap[i]=sub_martix[0][i];
            }
            bsp_put(bsp_pid()-sqroot, overlap, lower, 0, size*sizeof(float));
        }
        if (y>0) {
            for (i=0;i<size; i++) {
                overlap[i]=sub_martix[i][0];
            }
            bsp_put(bsp_pid()-1, overlap, right, 0, size*sizeof(float));
        }
        bsp_sync();
        
        for (i=0; i<size; i++) {
            for (j=0; j<size; j++) {
                temp = sub_martix[i][j];
                if (i-1<0) {
                    Aim1_j=upper[j];
                }
                else {
                    Aim1_j=sub_martix[i-1][j];
                }
                if (i+1>size-1) {
                    Aip1_j=lower[j];
                }
                else {
                    Aip1_j=sub_martix[i+1][j];
                }
                if (j-1<0) {
                    if (y!=0) {
                        Ai_jm1 = left[size-1];
                    }
                    else {
                        Ai_jm1 = left[i];
                    }
                }
                else {
                    Ai_jm1 = sub_martix[i][j-1];
                }
                if (j+1>size-1) {
                    if (y!=sqroot-1) {
                        Ai_jp1 = right[0];
                    }
                    else {
                        Ai_jp1 = right[i];
                    }
                }
                else {
                    Ai_jp1 = sub_martix[i][j+1];
                }
                sub_martix[i][j] = 0.2*(sub_martix[i][j]
                                        + Ai_jm1
                                        + Aim1_j
                                        + Ai_jp1
                                        + Aip1_j);
                //printf("data is %f\n", sub_martix[i][j]);
                diff += fabs(sub_martix[i][j]-temp);
            }
        }
        
        //printf("Result from pid: %d: difference= %f \n", bsp_pid(), diff);
        bsp_sync();
        
        for (i=0; i<bsp_nprocs(); i++) {
            bsp_get(i, &diff, 0, &diffs[i], sizeof(float));
        }
        bsp_sync();
        
        for (i=0; i<bsp_nprocs(); i++) {
            total_diff += diffs[i];
        }
        bsp_sync();
        convergence = (total_diff)/(float)(N*N);
        //printf("Current Convergence is %f\n", convergence);
        if (convergence<TOL) {
            done = 1;
        }
        bsp_sync();
    }
    
    for (i=0;i<size;i++) {
        bsp_pop_reg(sub_martix[i]);
    }
    bsp_pop_reg(&diff);
    bsp_pop_reg(lower);
    bsp_pop_reg(upper);
    bsp_pop_reg(left);
    bsp_pop_reg(right);
    bsp_sync();
    
    for (i=0; i<size; i++) {
        free(sub_martix[i]);
    }
    free(sub_martix);
    free(diffs);
    free(lower);
    free(upper);
    free(left);
    free(right);
    free(overlap);
    bsp_sync();
    bsp_end();
    for (i=0; i<N+2; i++) {
        free(matrix[i]);
    }
    free(matrix);
    
}

int main(int argc, char **argv)
{
    struct timeval start, end;
    float elapsed_time = 0.0;
    gettimeofday(&start, NULL);
    bsp_init(parallel_part, argc, argv);
    parallel_part();
    gettimeofday(&end, NULL);
    elapsed_time = (end.tv_sec - start.tv_sec)
    + (end.tv_usec - start.tv_usec) / 1000000000.0;
    printf("Total run time is %f seconds.\n", elapsed_time);
    return 0;
}