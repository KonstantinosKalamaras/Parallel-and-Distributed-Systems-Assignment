//used to transform coo -> csr
void cooToScr(int dim,int nz,coo ** pairs,int ** csr_row,int ** csr_col );

/*COO to CSR conversion*/
    cooToScr(dim,nz,&pairs,&csr_row,&csr_col);

void cooToScr(int dim,int nz,coo ** pairs,int ** csr_row,int ** csr_col ){
    int temp = 0;
    (*csr_row)[dim] = 2*nz-1;
    for(int k=0; k<2*nz; k++){
        (*csr_col)[k] = (*pairs)[k].j;

        if((*pairs)[k].i != temp){
            (*csr_row)[temp+1] = k;
            temp++;
        }
    }
}
