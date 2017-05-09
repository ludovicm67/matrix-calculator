#ifndef __MATRIX_H__
#define __MATRIX_H__

typedef float E;
typedef struct matrix {
    E * mat;
    unsigned int nb_rows;
    unsigned int nb_columns;
} * Matrix;

typedef struct {
    Matrix P;
    Matrix L;
    Matrix U;
} PLU;

Matrix new_matrix_copy(Matrix m);
Matrix newMatrix(unsigned int nb_rows, unsigned int nb_columns);
Matrix newMatrix_tab(unsigned int nb_rows, unsigned int nb_columns, E * tab);
E getElt(Matrix m, unsigned int row, unsigned int column);
void setElt(Matrix m, unsigned int row, unsigned int column, E val);
void deleteMatrix(Matrix m);
int isSquare(Matrix m);
int sameSize(Matrix a, Matrix b);
int isTriangulaire(Matrix m);
int isSymetric(Matrix m);
Matrix transpose(Matrix m);
void printMatrix(Matrix m);
Matrix matrix_identite(unsigned int n);
Matrix addition(Matrix m1, Matrix m2);
Matrix mult_scalar(E s, Matrix m);
Matrix multiplication(Matrix a, Matrix b);
Matrix extraction(Matrix m, unsigned int row, unsigned int column);
E det(Matrix m);
Matrix inversion(Matrix m);
void multiplier_ligne(Matrix m, unsigned int i, E k);
void permuter_ligne(Matrix m, unsigned int i, unsigned int j);
E m_determinant(Matrix m);
void copy_matrix(Matrix source, Matrix dest);
void addition_multiplication(Matrix m, unsigned int i, unsigned int j, E k);
Matrix triangulariser(Matrix m);
Matrix inversion_gauss(Matrix m);
PLU decomposition_PLU(Matrix m);
void m_PLU(Matrix m);
Matrix m_PLU_p(Matrix m);
Matrix m_PLU_l(Matrix m);
Matrix m_PLU_u(Matrix m);
E * valeurs_propres(Matrix m);

#endif
