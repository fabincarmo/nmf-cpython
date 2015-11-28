#include <cblas.h>

#define SUMESC(A,b,C,n,m) for(i=0; i<n;i++){for(j=0; j<m;j++){ C[j+i*m] = A[j+i*m] + b;}}
#define SUMM(A,B,C,n,m) for(i=0; i<n;i++){for(j=0; j<m;j++){ C[j+i*m] = A[j+i*m] + B[j+i*m];}}
#define MULT(A,B,C,n,m) for(i=0;i<n;i++){for(j=0;j<m;j++){C[i*m+j] = A[i*m+j] * B[i*m+j];}}
#define DIV(A,B,C,n,m) for(i=0;i<n;i++){for(j=0;j<m;j++){C[i*m+j] = A[i*m+j] / B[i*m+j];}}
#define MULTM(A,B,C,n,r,m) cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m, r, 1.0, A, r, B, m, 0.0, C, m);
#define MULTMtnt(A,B,C,n,r,m) cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, r, 1.0, A, n, B, m, 0.0, C, m);
#define MULTMntt(A,B,C,n,r,m) cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, r, 1.0, A, r, B, r, 0.0, C, m);
#define PRINTM(M,n,m) for(i=0;i<n;i++){for(j=0;j<m;j++){printf("%.8f ",M[j+i*m]);}printf("\n");}
#define EPS 1E-20
#define VEPS(M,n,m) for(i=0;i<n;i++){for(j=0;j<m;j++){if(M[j+i*m]<EPS){M[j+i*m] = EPS;}}}

static double * NormaColuna(double * M, const int n, const int m){
	int i, j;
	double norma;
	double *Mn = malloc(sizeof(double) * n*m);
	for(i=0;i<m;i++){
		norma = cblas_dnrm2(n,&M[i],m);
		for(j=0;j<n;j++){
			Mn[i+j*m] = M[i+j*m] / norma;
		}
	}
	return Mn;
}

static void nmf_kl(double *V, double *W, double *H, int n, int r, int m, int itmax){

	int it, i ,j;
	double *R, *aux1, *aux2, *Wn, *ones;

	R = malloc(sizeof(double) * n*m);
	Wn = malloc(sizeof(double) * n*r);
	ones = malloc(sizeof(double) * n*m);
	for(i=0; i<n;i++){for(j=0; j<m;j++){ ones[j+i*m] = 1.0;}}
    Wn = NormaColuna(W,n,r);
	for(it=0; it<itmax; it++){
		MULTM(Wn, H, R, n, r, m);

		aux1 = malloc(sizeof(double) * n*m);
		aux2 = malloc(sizeof(double) * r*m);
		VEPS(R,n,m);
		DIV(V,R,aux1,n,m);

		MULTMtnt(Wn,aux1,aux2,r,n,m);

		free(aux1);
		aux1 = malloc(sizeof(double) * r*m);
		MULTMtnt(Wn,ones,aux1,r,n,m);

		VEPS(aux1,r,m);
		DIV(aux2,aux1,aux2,r,m);
		MULT(H,aux2,H,r,m);
		free(aux1);
		free(aux2);
		MULTM(Wn, H, R, n, r, m);
		aux1 = malloc(sizeof(double) * n*m);
		aux2 = malloc(sizeof(double) * n*r);
		VEPS(R,n,m);
		DIV(V,R,aux1,n,m);
		MULTMntt(aux1,H,aux2,n,m,r);
		free(aux1);
		aux1 = malloc(sizeof(double) * n*r);
		MULTMntt(ones,H,aux1,n,m,r);
		VEPS(aux1,n,r);
		DIV(aux2,aux1,aux2,n,r);
		MULT(Wn,aux2,W,n,r);
		free(aux1);
		free(aux2);
		Wn = NormaColuna(W,n,r);
	}
	W = Wn;
	return;
}

static void nmf_kl_h(double *V, double *W, double *H, int n, int r, int m, int itmax){

	int it, i ,j;
	double *R, *aux1, *aux2, *ones;

	R = malloc(sizeof(double) * n*m);
	ones = malloc(sizeof(double) * n*m);
	for(i=0; i<n;i++){for(j=0; j<m;j++){ ones[j+i*m] = 1.0;}}
	for(it=0; it<itmax; it++){
		MULTM(W, H, R, n, r, m);

		aux1 = malloc(sizeof(double) * n*m);
		aux2 = malloc(sizeof(double) * r*m);
		VEPS(R,n,m);
		DIV(V,R,aux1,n,m);

		MULTMtnt(W,aux1,aux2,r,n,m);

		free(aux1);
		aux1 = malloc(sizeof(double) * r*m);
		MULTMtnt(W,ones,aux1,r,n,m);

		VEPS(aux1,r,m);
		DIV(aux2,aux1,aux2,r,m);
		MULT(H,aux2,H,r,m);
		free(aux1);
		free(aux2); 
	}
	return;
}
