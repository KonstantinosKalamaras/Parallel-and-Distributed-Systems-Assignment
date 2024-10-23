#include <stdio.h>
#include <stdlib.h>
#include "mmio.c"
#include <time.h>


struct timespec diff(struct timespec start, struct timespec end)
{
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
    int *csr_row, *csr_col;
    csr_col = malloc(2*nz*sizeof(int));
    csr_row = malloc((dim+1)*sizeof(int));

    for(i=0; i<dim+1; i++){
        csr_row[i] = 0;
    }

    cooToCsr(dim,nz,&pairs,&csr_row,&csr_col);

    free(pairs);

    /*-------------------------------------------*/
    char answer;
    printf("Do you want to print CSR? [Y/N]: ");
    scanf("%c", &answer);

    if(answer == 'Y'){
        printf("Col_index is: ");
        for(i=0;i<2*nz;i++){
            printf(" %d ", csr_col[i]);
        }
        printf("\nRow_index is: ");
        for(i=0;i<dim+1;i++)
            printf(" %d ", csr_row[i]);
    }
    printf("\n-------------------------------\n");

    /*----------------------------------------------*/

    int val;
    struct timespec tStart,tEnd;
    double minTime = 100;

    int *triangles = malloc(dim*sizeof(int));
    for(i=0; i<dim ; i++) triangles[i]=0;

    ///////custom matrix with one triangle/////////
/*    nz = 6;
    int* csr_col2 = malloc(12*sizeof(int)); 
    int* csr_row2 = malloc(6*sizeof(int));
    int gakldlj[12] = {1,3,4,0,2,1,4,0,4,0,2,3}; 
    int paapapa[6] =     {0,3,5,7,9,12};
    for(i=0;i<12;i++)  csr_col2[i] = gakldlj[i];
    for(i=0;i<6;i++)  csr_row2[i] = paapapa[i];
*/
    /////////////////////////////

    clock_gettime(CLOCK_MONOTONIC, &tStart);
    
    /*only calculate square for non zeros of matrix(skips Hadamard)*/
    int currentRow = 0;

    for(int k = 0; k < 2*nz ; k++){
        if(k>=csr_row[currentRow+1]){
            triangles[currentRow]/=2;
            currentRow++;
        }
        val = rowMul(currentRow, csr_col[k], &csr_col, &csr_row);
        triangles[currentRow] += val;
    }
    triangles[dim-1]/=2;
    
    int all_triangles = 0 ;
    printf("\nPrinting triangle vector\n");
    for(i=0; i<dim ; i++){
        all_triangles += triangles[i];
        //printf("%d ",triangles[i]);
    }
    printf("\nTotal triangles are:%d \n", all_triangles/3);


 /*------------------------------------------------------------------------*/
    clock_gettime(CLOCK_MONOTONIC, &tEnd);
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
//NOT USED RIGHT NOW 
/*Calculates vector*vector (product) of row i and row j of csr matrix
 using the fact that vectors are ordered*/
int rowMul(int ri, int rj, int **csr_col, int **csr_row){

    int k = 0, l = 0, result = 0;
    int col_starti = (*csr_row)[ri]; 
    int col_endi = (*csr_row)[ri+1] ; 
    int diffi = col_endi - col_starti;

    int col_startj = (*csr_row)[rj]; 
    int col_endj = (*csr_row)[rj+1] ; 
    int diffj = col_endj - col_startj;
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

