#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

typedef float E;
typedef struct matrix {
    E * mat;
    unsigned int nb_rows;
    unsigned int nb_columns;
} * Matrix;

typedef struct {
    Matrix L;
    Matrix U;
} PLU;

// Permet de générer une nouvelle matrice
Matrix newMatrix(unsigned int nb_rows, unsigned int nb_columns) {
    unsigned int i;
    Matrix m = malloc(sizeof(Matrix));
    if(!m) return m;
    m->nb_rows = nb_rows;
    m->nb_columns = nb_columns;
    m->mat = malloc(nb_rows * nb_columns * sizeof(E));
    for(i = 0; i < nb_rows * nb_columns; i++) m->mat[i] = 0;
    return m;
}

// Permet de récupérer un élément d'une matrice
E getElt(Matrix m, unsigned int row, unsigned int column) {
    return m->mat[m->nb_columns*row + column];
}

// Permet de définir un élément d'une matrice
void setElt(Matrix m, unsigned int row, unsigned int column, E val) {
    m->mat[m->nb_columns*row + column] = val;
}

// Permet de supprimer une matrice
void deleteMatrix(Matrix m) {
    free(m->mat);
    free(m);
}

// Teste si une matrice est carré
int isSquare(Matrix m) {
    return m->nb_rows == m->nb_columns;
}

// Teste si une matrice est symétrique
int isSymetric(Matrix m) {
    unsigned int i, j, n = 0;
    if(!isSquare(m)) return 0;
    for(i = 0; i < m->nb_rows; i++) {
        for(j = 0; j < m->nb_columns; j++) {
            if(getElt(m, i, j) != getElt(m, j, i)) {
                return 0;
            }
            n++;
        }
    }
    return 1;
}

// Retourne la transposée de la matrice
Matrix transpose(Matrix m) {
    unsigned int i, j;

    if(!isSquare(m)) {
        fprintf(stderr, "La matrice n'est pas carré\n");
        exit(EXIT_FAILURE);
    }

    Matrix r = newMatrix(m->nb_columns, m->nb_rows);

    for(i = 0; i < m->nb_rows; i++) {
        for(j = 0; j < m->nb_columns; j++) {
            setElt(r, i, j, getElt(m, j, i));
        }
    }

    return r;
}

// Permet d'afficher une matrice
void printMatrix(Matrix m) {

    if(!m) {
        printf("Y'a comme un bug dans la matrice (:\n");
        return;
    }

    unsigned int i, j;

    printf("\nAFFICHAGE DE LA MATRICE %dx%d :\n", m->nb_rows, m->nb_columns);
    for(i = 0; i < m->nb_rows; i++) {
        for(j = 0; j < m->nb_columns; j++) {
            printf("\t%.0f", getElt(m, i, j));
        }
        printf("\n");
    }
    printf("\n");

}

Matrix matrix_identite(unsigned int n) {
    unsigned int i;
    Matrix m = newMatrix(n, n);
    for (i = 0; i < n; i++) setElt(m, i, i, 1);
    return m;
}

// Additionne deux matrices
// == pré-condition : m1 et m2 doivent être de même dimension
Matrix addition(Matrix m1, Matrix m2) {
    Matrix m = newMatrix(m1->nb_rows, m1->nb_columns);
    unsigned int i, j;

    for(i = 0; i < m->nb_rows; i++) {
        for(j = 0; j < m->nb_columns; j++) {
            setElt(m, i, j, getElt(m1, i, j) + getElt(m2, i, j));
        }
    }

    return m;
}

// Multiplie une matrice par un scalaire
Matrix mult_scalar(E s, Matrix m) {
    unsigned int i;
    for(i = 0; i < m->nb_rows * m->nb_columns; i++) {
        m->mat[i] *= s;
    }
    return m;
}

// Permet d'effectuer la multiplicaton de deux matrices
Matrix multiplication(Matrix a, Matrix b) {
    if(a->nb_columns != b->nb_rows) return NULL;

    Matrix r = newMatrix(a->nb_rows, b->nb_columns);
    unsigned int i, j, k;
    E val;

    for(i = 0; i < r->nb_rows; i++) {
        for(j = 0; j < r->nb_columns; j++) {
            val = 0;
            for(k = 0; k < a->nb_columns; k++) {
                val += getElt(b, k, j) * getElt(a, i, k);
            }
            setElt(r, i, j, val);
        }
        printf("\n");
    }

    return r;
}

Matrix extraction(Matrix m, unsigned int row, unsigned int column) {
    if(m->nb_rows <= row || m->nb_columns <= column) return m;

    Matrix r = newMatrix(m->nb_rows-1, m->nb_columns-1);
    unsigned int i, j, x, y;

    for(i = 0; i < r->nb_rows; i++) {
        for(j = 0; j < r->nb_columns; j++) {
            x = i;
            y = j;
            if(i >= row) x++;
            if(j >= column) y++;
            setElt(r, i, j, getElt(m, x, y));
        }
    }

    return r;
}

double det(Matrix m) {
	if(isSquare(m)) {
		if(m->nb_rows == 0) return 0;
		if(m->nb_rows == 1) return getElt(m, 0, 0);
		int i;
		double somme = 0;
		for(i = 0; i < m->nb_rows; i++) {
			somme += pow((-1), i)* getElt(m, 0, i)* (det(extraction(m, 0, i)));
		}
		return somme;
	}

}

// multiplier la ligne i par un facteur k
void multiplier_ligne(Matrix m, unsigned int i, E k) {

    unsigned int j;

    if(i < 0 || i >= m->nb_rows) {
        fprintf(stderr, "La ligne %d n'existe pas dans cette matrice\n", i);
        exit(EXIT_FAILURE);
    }

    for(j = 0; j < m->nb_columns; j++)
        setElt(m, i, j, k * getElt(m, i, j));
}

// permute les lignes i et j de la matrice m
Matrix permuter_ligne(Matrix m, unsigned int i, unsigned int j) {
    int k;
    // check : si a et b < m->nb_rows
    if((i < m->nb_rows)&&(j < m->nb_rows)) {
        for(k = 0; k < m->nb_columns; k++) {
            E tmp_i = getElt(m, i, k);
            E tmp_j = getElt(m, j, k);
            setElt(m, i, k, tmp_j);
            setElt(m, j, k, tmp_i);
        }
    }
    else{
        return NULL;
    }

    return m;
}

// prec : m doit être triangulariser
E m_determinant(Matrix m) {
    unsigned int i;
    E determinant = 1;

    if(!isSquare(m)) {
        fprintf(stderr, "La matrice n'est pas carré\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m->nb_rows; i++) {
        determinant *= getElt(m, i, i);
    }

    return determinant;
}

// pre-cond : dest et in de même dimensions
void copy_matrice(Matrix in, Matrix dest) {
    unsigned int i, j;
    for(i = 0; i < in->nb_rows; i++) {
        for(j = 0; j < in->nb_columns; j++) {
            setElt(dest, i, j, getElt(in, i, j));
        }
    }
}

void addition_multiplication(Matrix m, unsigned int i, unsigned int j, E k){
    int n;
    Matrix m1, m2 = newMatrix(m->nb_rows, m->nb_columns);
    for(n = 0; n < m->nb_columns; n++){
        E tmp = getElt(m, j, n);
        setElt(m2, i, n, tmp);
    }
    multiplier_ligne(m2, i, k);
    m1 = addition(m, m2);
    copy_matrice(m1, m);

    deleteMatrix(m1);
    deleteMatrix(m2);
}

void triangulariser(Matrix m) {

    unsigned int i, j;
    E k;

    if(!isSquare(m)) {
        fprintf(stderr, "La matrice n'est pas carré\n");
        exit(EXIT_FAILURE);
    }

    for(i = 0; i < m->nb_rows-1; i++) {
        for(j = i+1; j < m->nb_rows; j++) {
            k = -getElt(m, j, i) / getElt(m, i, i);
            addition_multiplication(m, j, i, k);
        }
    }

}


int main() {

    // Matrix m = newMatrix(3, 4);
    // m->mat[0] = 0;
    // m->mat[1] = 1;
    // m->mat[2] = 2;
    // m->mat[3] = 3;
    // m->mat[4] = 4;
    // m->mat[5] = 5;
    // m->mat[6] = 6;

    // Matrix w = newMatrix(3, 4);
    // w->mat[0] = 40;
    // w->mat[1] = 41;
    // w->mat[2] = 42;
    // w->mat[3] = 43;
    // w->mat[4] = 44;
    // w->mat[5] = 45;
    // w->mat[6] = 46;

    // setElt(m, 1, 2, 42);
    // assert(getElt(m, 1, 2) == 42);
    // printMatrix(m);
    // printMatrix(w);
    // Matrix t = addition(m, w);
    // printMatrix(t);

    // Matrix m2 = newMatrix(3, 3);
    // for (int i = 0; i < 3 * 3; i++) {
    //  m2->mat[i] = i;
    // }

    // printMatrix(transpose(m2));


    // Matrix m1 = newMatrix (2, 2);
 //    setElt (m1, 0, 0, 2);
 //    setElt (m1, 0, 1, -3);
 //    setElt (m1, 1, 0, -1);
 //    setElt (m1, 1, 1, 2);
    // printMatrix(m1);

 //    m2 = newMatrix (2, 3);
 //    setElt (m2, 0, 0, 3);
 //    setElt (m2, 0, 1, 1);
 //    setElt (m2, 0, 2, 2);
 //    setElt (m2, 1, 0, 1);
 //    setElt (m2, 1, 1, 0);
 //    setElt (m2, 1, 2, 2);
    // printMatrix(m2);
    // printf("==");
    // printMatrix(multiplication(m1, m2));
    // printf("===WOLOLO==\n");
    // printMatrix(extraction(w, 0, 0));


    Matrix m = newMatrix(4, 4);
    int n = 1;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            setElt(m, i, j, n++);
        }
    }

    multiplier_ligne(m, 0, 2);
    permuter_ligne(m, 0, 1);
    printMatrix(m);
    // m = addition_multiplication(m, 1, 0, 2);
    // printMatrix(m);

    triangulariser(m);
    printMatrix(m);

    printMatrix(matrix_identite(4));

    printf("DET = %f\n", det(m));
    printf("DET = %f\n", m_determinant(m));

    return EXIT_SUCCESS;
}
