#include <cstdio>
#include <cstring>
#include <mkl.h>
using namespace std;
void imprimeMatriz(double *A, int fil, int col);
int main(int argc, char *argv[]) {
	double A[16] = {
		5, 5, 7, 9,
		4, 8, 8, 10,
		5, 7, 11, 11,
		6, 8, 10, 14 };
	int pivot[4];
	double A2[16];
	memcpy(A2, A, 16 * sizeof(double));
	//Obtención de L y U
	int result = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 4, 4, A2, 4, pivot);
	if (result != 0) {
		fprintf(stderr, "Fallo en dgetrf: %d\n", result);
		return 1;
	}
	//Obtención del determinante
	double determinante = 1;
	for (int i = 0; i < 4; i++) {
		determinante *= pivot[i] == i + 1 ? A2[i * 4 + i] : -A2[i * 4 + i];
	}
	//Obtención de la inversa 1
	double B1[16];
	memset(B1, 0, 16 * sizeof(double));
	for (int i = 0; i < 4; i++) B1[i * 4 + i] = 1;
	result = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', 4, 4, A2, 4, pivot, B1, 4);
	if (result != 0) {
		fprintf(stderr, "Fallo en dgetrs: %d\n", result);
		return 2;
	}
	//Obtención de la inversa 2
	double B2[16];
	memcpy(B2, A2, 16 * sizeof(double));
	result = LAPACKE_dgetri(LAPACK_ROW_MAJOR, 4, B2, 4, pivot);
	if (result != 0) {
		fprintf(stderr, "Fallo en dgetri: %d\n", result);
		return 3;
	}
	printf("Matriz con L y U mezcladas:\n");
	imprimeMatriz(A2, 4, 4);
	printf("\nVector de pivotamiento:\n");
	for (int i = 0; i < 4; i++) {
		printf("%d ", pivot[i]);
	}
	printf("\n\nDeterminate: %lf\n", determinante);
	printf("\nMatriz inversa metodo 1:\n");
	imprimeMatriz(B1, 4, 4);
	printf("\nMatriz inversa metodo 2:\n");
	imprimeMatriz(B2, 4, 4);
	getchar();
	return 0;
}
void imprimeMatriz(double *A, int fil, int col) {
	for (int i = 0; i < fil; i++) {
		for (int j = 0; j < col; j++) printf("%lf ", A[fil*i + j]);
		printf("\n");
	}
}