#include <stdio.h>
#include <stdlib.h>
#include "mmio.c"
#include <time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>


struct timespec diff(struct timespec start, struct timespec end){
        struct timespec temp;
        if ((end.tv_nsec - start.tv_nsec) < 0) 
        {
                temp.tv_sec = end.tv_sec - start.tv_sec - 1;
                temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
        } 
        else 
        {
                temp.tv_sec = end.tv_sec - start.tv_sec;
                temp.tv_nsec = end.tv_nsec - start.tv_nsec;
        }
        return temp;
}

/*matrix in coo format is stores as array of coo structs
Useful for qsort() */
typedef struct coo{
    int i;
    int j;
    int value;
}coo;


/*compare to be used with qsort*/
int compareCoo(const void *a, const void *b);

/*Given the .mtx file this function extracts coordinates of non zeroes*/
int* getMatrixData(char* mat, int** I, int** J);

void cooToCsr(int dim,int nz,coo ** pairs,int ** csr_row,int ** csr_col );

/*csr_col length is 2*nz and csr_row dim+1*/
int rowMul(int ri, int rj, int **csr_col, int **csr_row);


int main(int argc, char *argv[]){
    int i, *data, *I, *J;

    data = getMatrixData(argv[1], &I, &J);
    int dim = data[0];
    int nz = data[1];

    /* non zero elements are only given once in .mtx file */
    printf("\nNumber of rows is %d and number of non zeroes is: %d \n", dim, 2*nz);

    /*After getting data create coo form of array. Value is 1 and will be ignored for now.*/
    coo* pairs ;
    pairs = malloc(2*nz*sizeof(coo));

    /*Since matrix is symmetrical we need to account for twice as many elements*/
    for (int k=0; k<nz; k++){
        pairs[k].i = J[k];
        pairs[k].j = I[k];
        pairs[2*nz-k-1].i = I[k];
        pairs[2*nz-k-1].j = J[k];
    }

    /*Always free malloced memory*/
    free(I);
    free(J);

    /*using a custom comparator function coo is (i, j) sorted*/
    qsort((void*)pairs, 2*nz, sizeof(pairs[0]), compareCoo);

    printf("\n-------------------------------\n");
    
    /*COO to CSR conversion*/
    int *A_row, *A_col;
    A_col = malloc(2*nz*sizeof(int));
    A_row = malloc((dim+1)*sizeof(int));

    for(i=0; i<dim+1; i++){
        A_row[i] = 0;
    }

    cooToCsr(dim,nz,&pairs,&A_row,&A_col);

    free(pairs);

    /*-------------------------------------------*/
    char answer;
    printf("Do you want to print CSR? [Y/N]: ");
    scanf("%c", &answer);

    if(answer == 'Y'){
        printf("Col_index is: ");
        for(i=0;i<2*nz;i++){
            printf(" %d ", A_col[i]);
        }
        printf("\nRow_index is: ");
        for(i=0;i<dim+1;i++)
            printf(" %d ", A_row[i]);
    }
    printf("\n-------------------------------\n");

    /*----------------------------------------------*/


    struct timespec tStart,tEnd;
    double minTime = 100;

    int* C_val;
    C_val = malloc(2*nz*sizeof(int));
    for (i = 0; i < 2*nz; ++i)  C_val[i] = 0 ;


    int dif,starti, startj, endi, endj ; 
    clock_gettime(CLOCK_MONOTONIC, &tStart);
    
    /*only calculate square for non zeros of matrix(skips Hadamard)*/
    cilk_for(int row = 0; row < dim; row++){
        cilk_for(int col = A_row[row] ; col < A_row[row+1]; col++){
            /*current index to compare*/
            starti = A_row[row];
            startj = A_row[ A_col[col] ];
            
            /*end of rows*/
            endi = A_row[row+1];
            endj = A_row[A_col[col]+1];

            /*find common elements*/
            while(starti < endi && startj < endj){
                dif = A_col[starti] - A_col[startj]; 

                /*start with comparisons most likely to occur*/
                if(dif<0){   
                     starti++;  
                }
                else if(dif>0){
                    startj++;
                }
                else{
                    C_val[col]++;
                    starti++;
                    startj++;
                }     
            }
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &tEnd);

    FILE *fp;
    fp = fopen("cilkS.txt", "w");
    for (int i = 0; i < 2*nz; ++i)
    {
        fprintf(fp, "%d", C_val[i]);
    }
    fclose(fp);

    printf("\nEnd of calculations\n");
        
 /*------------------------------------------------------------------------*/
    struct timespec tResult = diff(tStart,tEnd);
    double timeR = (double)(tResult.tv_sec+(double)tResult.tv_nsec/1000000000);
    printf( "Benchmark = %f s\n", timeR);

    return 0;
}

int* getMatrixData(char* mat, int** I, int** J){
    FILE *f;
    int ret_code, M, N, nz, i;   
    int* ret = malloc(2*sizeof(int));


    if ((f = fopen(mat, "r")) == NULL) {
            printf("Unexpected argument\n");
            exit(1);
    }
    
    /* find out size of sparse matrix */
    ret_code = mm_read_mtx_crd_size(f, &ret[0], &N, &ret[1]);

    /* reseve memory for matrices */
    *I = malloc(ret[1] * sizeof(int));
    *J = malloc(ret[1] * sizeof(int));

    for (i=0; i<ret[1]; i++){
        fscanf(f, "%d %d \n", &(*I)[i], &(*J)[i]);
        (*I)[i]--;  /* adjust from 1-based to 0-based */
        (*J)[i]--;    
    }
    return ret;
}

int compareCoo(const void *c1, const void *c2){
    const struct coo *a = c1;
    const struct coo *b = c2;

    if(a->i < b->i)
        return -1;
    else if(a->i > b->i)
        return 1;
    else if(a->i == b->i){
        if(a->j < b->j)
            return -1;
        else if(a->j > b->j)
            return 1;
        else 
            return 0;
    } 
    return 2;
}

void cooToCsr(int dim,int nz,coo ** pairs,int ** csr_row,int ** csr_col ){
    int temp = 0;
    (*csr_row)[dim] = 2*nz;
    for(int k=0; k<2*nz; k++){
        (*csr_col)[k] = (*pairs)[k].j;

        if((*pairs)[k].i != temp){
            (*csr_row)[temp+1] = k;
            temp++;
        }
    }
}

/*Calculates vector*vector (product) of row i and row j of csr matrix
 using the fact that vectors are ordered*/
int rowMul(int ri, int rj, int **csr_col, int **csr_row){

    int k = 0, l = 0, result = 0;
    int col_starti = (*csr_row)[ri]; 
    int col_endi = (*csr_row)[ri+1] ; 
    int diffi = (*csr_row)[ri+1] - (*csr_row)[ri];

    int col_startj = (*csr_row)[rj]; 
    int col_endj = (*csr_row)[rj+1] ; 
    int diffj = (*csr_row)[rj+1] - (*csr_row)[rj];
    int dif =0;

    /*find common elements*/
    while(k < diffi && l < diffj){
        dif = (*csr_col)[col_starti+k] - (*csr_col)[col_startj+l]; 
        if(dif<0){   
             k++;  
        }
        else if(dif>0){
            l++;
        }
        else{
            result++;
            k++;
            l++;
        }     
    }
    return result;
}

