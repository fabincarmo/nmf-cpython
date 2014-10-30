#include <cblas.h>

#define SUMESC(A,b,C,n,m) for(i=0; i<n;i++){for(j=0; j<m;j++){ C[j+i*m] = A[j+i*m] + b;}}
#define SUMM(A,B,C,n,m) for(i=0; i<n;i++){for(j=0; j<m;j++){ C[j+i*m] = A[j+i*m] + B[j+i*m];}}
#define MULT(A,B,C,n,m) for(i=0;i<n;i++){for(j=0;j<m;j++){C[i*m+j] = A[i*m+j] * B[i*m+j];}}
#define DIV(A,B,C,n,m) for(i=0;i<n;i++){for(j=0;j<m;j++){C[i*m+j] = A[i*m+j] / B[i*m+j];}}
#define MULTM(A,B,C,n,r,m) cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, m, r, 1.0, A, r, B, m, 0.0, C, m);
#define MULTMtnt(A,B,C,n,r,m) cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, m, r, 1.0, A, n, B, m, 0.0, C, m);
#define MULTMntt(A,B,C,n,r,m) cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, r, 1.0, A, r, B, r, 0.0, C, m);
#define PRINTM(M,n,m) for(i=0;i<n;i++){for(j=0;j<m;j++){printf("%.5f ",M[j+i*m]);}printf("\n");}
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

static double err(double *A, double *B, double *C, int n, int r, int m){
	int i, j;
	double e = 0.;
	double *R = malloc(sizeof(double) * n*m);
	MULTM(A,B,R,n,r,m);
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			e = e + (C[j+i*m] - R[j+i*m]) * (C[j+i*m] - R[j+i*m]);
		}
	}
	return e;
}

static void nmf_euc(double *V, double *W, double *H, int n, int r, int m, int itmax){

	int it, i ,j;
	double *R, *aux1, *aux2, *Wn;

	R = malloc(sizeof(double) * n*m);
	Wn = malloc(sizeof(double) * n*r);

	for(it=0; it<itmax; it++){
		Wn = NormaColuna(W,n,r);
		MULTM(Wn, H, R, n, r, m);
		aux1 = malloc(sizeof(double) * r*m);
		aux2 = malloc(sizeof(double) * r*m);
		MULTMtnt(Wn, V, aux1, r,n,m);
		MULTMtnt(Wn, R, aux2, r,n,m);
		VEPS(aux2,r,m);
		DIV(aux1, aux2, aux1, r,m);
		MULT(H, aux1, H, r,m);
		free(aux1);
		free(aux2);

		aux1 = malloc(sizeof(double) * n*r);
		aux2 = malloc(sizeof(double) * n*r);
		MULTM(Wn, H, R, n,r,m);
		MULTMntt(V, H, aux1, n,m,r);
		MULTMntt(R, H, aux2, n,m,r);
		VEPS(aux2,n,r);
		DIV(aux1, aux2, aux1, n,r);
		MULT(W, aux1, W, n,r);
		free(aux1);
		free(aux2);
	}

	W = NormaColuna(W,n,r);
	printf("Erro: %g\n ",err(W,H,V,n,r,m));
	return;
}

static void nmf_euc_sparse(double *V, double *W, double *H, int n, int r, int m, int itmax, double lambda){

	int i, j, it;
	double *R, *Wn, *Ones, *aux1, *aux2, *aux3;

	R = malloc(sizeof(double) * n*m);
	Wn = malloc(sizeof(double) * n*r);
	Ones = malloc(sizeof(double) * n*n);
	for(i=0; i<n;i++){for(j=0; j<n;j++){ Ones[j+i*n] = 1.0;}}

	for(it=0; it<itmax; it++){
		Wn = NormaColuna(W,n,r);
		MULTM(Wn, H, R, n, r, m);
		aux1 = malloc(sizeof(double) * r*m);
		aux2 = malloc(sizeof(double) * r*m);
		MULTMtnt(Wn, V, aux1, r,n,m);
		MULTMtnt(Wn, R, aux2, r,n,m);
		SUMESC(aux2,lambda,aux2,r,m);
		VEPS(aux2,r,m);
		DIV(aux1, aux2, aux1, r,m);
		MULT(H, aux1, H, r,m);
		free(aux1);
		free(aux2);

		aux1 = malloc(sizeof(double) * n*r);
		aux2 = malloc(sizeof(double) * n*r);
		aux3 = malloc(sizeof(double) * n*r);
		MULTM(Wn, H, R, n,r,m);
		MULTMntt(R, H, aux1, n,m,r);
		MULT(aux1,Wn,aux1,n,r);
		MULTM(Ones,aux1,aux1,n,n,r);
		MULTMntt(V,H,aux2,n,m,r);
		SUMM(aux1,aux2,aux1,n,r);
		MULTM(Ones,aux2,aux2,n,n,r);
		MULT(aux2,Wn,aux2,n,r);
		MULT(Wn,aux2,aux2,n,r);
		MULTMntt(R, H, aux3, n,m,r);
		SUMM(aux2,aux3,aux2,n,r);
		VEPS(aux2,n,r);
		DIV(aux1, aux2, aux1, n,r);
		MULT(Wn, aux1, W, n,r);
		free(aux1);
		free(aux2);
		free(aux3);
	}
	W = NormaColuna(W,n,r);
	printf("Erro: %g\n ",err(W,H,V,n,r,m));

	return;
}

