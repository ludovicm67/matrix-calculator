#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "system.h"

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

// Permet de générer une nouvelle matrice
Matrix newMatrix(unsigned int nb_rows, unsigned int nb_columns) {
    Matrix m = malloc(sizeof(struct matrix));
    if (!m) {
        print_error("Cannot allocate memory for the matrix");
        exit(EXIT_FAILURE);
    }
    m->nb_rows = nb_rows;
    m->nb_columns = nb_columns;
    m->mat = calloc(nb_rows * nb_columns, sizeof(E));
    return m;
}

// Permet de générer une nouvelle matrice depuis un tableau
Matrix newMatrix_tab(unsigned int nb_rows, unsigned int nb_columns, E * tab) {
    Matrix m = newMatrix(nb_rows, nb_columns);
    memcpy(m->mat, tab, nb_rows * nb_columns * sizeof(E));
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
    if (!m) return;
    free(m->mat);
    free(m);
}

// pre-cond : dest et source de même dimensions
void copy_matrix(Matrix source, Matrix dest) {
    unsigned int i, j;
    for (i = 0; i < source->nb_rows; i++) {
        for (j = 0; j < source->nb_columns; j++) {
            setElt(dest, i, j, getElt(source, i, j));
        }
    }
}

// Teste si une matrice est carré
int isSquare(Matrix m) {
    return m->nb_rows == m->nb_columns;
}

// Teste si une matrice est symétrique
int isSymetric(Matrix m) {
    unsigned int i, j, n = 0;
    if (!isSquare(m)) return 0;
    for (i = 0; i < m->nb_rows; i++) {
        for (j = 0; j < m->nb_columns; j++) {
            if (getElt(m, i, j) != getElt(m, j, i)) {
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

    Matrix r = newMatrix(m->nb_columns, m->nb_rows);

    for (i = 0; i < m->nb_rows; i++) {
        for (j = 0; j < m->nb_columns; j++) {
            setElt(r, j, i, getElt(m, i, j));
        }
    }

    return r;
}

// Permet d'afficher une matrice
void printMatrix(Matrix m) {

    if (!m) {
        print_error("No matrix to print");
        return;
    }

    unsigned int i, j;

    printf("╭");
    for (i = 0; i < m->nb_columns+1; i++) printf("        ");
    printf("       ╮\n");
    for (i = 0; i < m->nb_rows; i++) {
        printf("│   ");
        for (j = 0; j < m->nb_columns; j++) {
            printf("\t%.2f", getElt(m, i, j));
        }
        printf("\t\t│\n");
    }
    printf("╰");
    for (i = 0; i < m->nb_columns+1; i++) printf("        ");
    printf("       ╯\n");
}

// Matrice identité
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

    for (i = 0; i < m->nb_rows; i++) {
        for (j = 0; j < m->nb_columns; j++) {
            setElt(m, i, j, getElt(m1, i, j) + getElt(m2, i, j));
        }
    }

    return m;
}

// Multiplie une matrice par un scalaire
Matrix mult_scalar(E s, Matrix m) {
    Matrix r = newMatrix(m->nb_rows, m->nb_columns);
    unsigned int i;

    copy_matrix(m, r);
    for (i = 0; i < m->nb_rows * m->nb_columns; i++) {
        r->mat[i] *= s;
    }

    return r;
}

// Permet d'effectuer la multiplicaton de deux matrices
Matrix multiplication(Matrix a, Matrix b) {
    if (a->nb_columns != b->nb_rows) return NULL;

    Matrix r = newMatrix(a->nb_rows, b->nb_columns);
    unsigned int i, j, k;
    E val;

    for (i = 0; i < r->nb_rows; i++) {
        for (j = 0; j < r->nb_columns; j++) {
            val = 0;
            for (k = 0; k < a->nb_columns; k++) {
                val += getElt(b, k, j) * getElt(a, i, k);
            }
            setElt(r, i, j, val);
        }
    }

    return r;
}

Matrix extraction(Matrix m, unsigned int row, unsigned int column) {
    if (m->nb_rows <= row || m->nb_columns <= column) return m;

    Matrix r = newMatrix(m->nb_rows-1, m->nb_columns-1);
    unsigned int i, j, x, y;

    for (i = 0; i < r->nb_rows; i++) {
        for (j = 0; j < r->nb_columns; j++) {
            x = i;
            y = j;
            if(i >= row) x++;
            if(j >= column) y++;
            setElt(r, i, j, getElt(m, x, y));
        }
    }

    return r;
}

E det(Matrix m) {
    unsigned int i;
    E somme = 0;

    if (!isSquare(m)) {
        fprintf(stderr, "La matrice n'est pas carrée\n");
        exit(EXIT_FAILURE);
    }

    if (m->nb_rows == 0) return 0;
    else if (m->nb_rows == 1) return getElt(m, 0, 0);

    // Cas de la matrice 2x2
    if ((m->nb_rows == 2) && (m->nb_columns == 2)) {
        return getElt(m, 0, 0) * getElt(m, 1, 1) - (getElt(m, 0, 1) * getElt(m, 1, 0));
    }

    for (i = 0; i < m->nb_rows; i++) {
        somme += (i%2 == 0 ? 1 : -1) * getElt(m, 0, i) * (det(extraction(m, 0, i)));
    }

    return somme;
}

Matrix inversion(Matrix m) {
    int i, j;

    if (!det(m)) {
        return NULL;
    }

    Matrix inverse = newMatrix(m->nb_rows, m->nb_columns);
    for (i = 0; i < (int) m->nb_rows; i++) {
        for (j = 0; j < (int) m->nb_columns; j++) {
            if (((-1 * (j % 2)) -1 * (i % 2)) == -1) {
                setElt(inverse, i, j, -det(extraction(m, i, j)));
            } else {
                setElt(inverse, i, j, det(extraction(m, i, j)));
            }
        }
    }
    mult_scalar((1/det(m)), inverse);
    return inverse;
}

// multiplier la ligne i par un facteur k
void multiplier_ligne(Matrix m, unsigned int i, E k) {

    unsigned int j;

    if (i >= m->nb_rows) {
        fprintf(stderr, "La ligne %d n'existe pas dans cette matrice\n", i);
        exit(EXIT_FAILURE);
    }

    for (j = 0; j < m->nb_columns; j++)
        setElt(m, i, j, k * getElt(m, i, j));
}

// permute les lignes i et j de la matrice m
void permuter_ligne(Matrix m, unsigned int i, unsigned int j) {
    unsigned int k;
    // check : si a et b < m->nb_rows
    if ((i < m->nb_rows) && (j < m->nb_rows)) {
        for(k = 0; k < m->nb_columns; k++) {
            E tmp_i = getElt(m, i, k);
            E tmp_j = getElt(m, j, k);
            setElt(m, i, k, tmp_j);
            setElt(m, j, k, tmp_i);
        }
    } else return;
}

// prec : m doit être triangulariser
E m_determinant(Matrix m) {
    unsigned int i;
    E determinant = 1;

    if (!isSquare(m)) {
        fprintf(stderr, "La matrice n'est pas carrée\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m->nb_rows; i++) {
        determinant *= getElt(m, i, i);
    }

    return determinant;
}

void addition_multiplication(Matrix m, unsigned int i, unsigned int j, E k) {
    unsigned int n;
    Matrix m1, m2 = newMatrix(m->nb_rows, m->nb_columns);
    for (n = 0; n < m->nb_columns; n++){
        E tmp = getElt(m, j, n);
        setElt(m2, i, n, tmp);
    }
    multiplier_ligne(m2, i, k);
    m1 = addition(m, m2);
    copy_matrix(m1, m);

    deleteMatrix(m1);
    deleteMatrix(m2);
}

void triangulariser(Matrix m) {

    unsigned int i, j;
    E k;

    if (!isSquare(m)) {
        fprintf(stderr, "La matrice n'est pas carrée\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < m->nb_rows-1; i++) {
        for (j = i+1; j < m->nb_rows; j++) {
            k = -getElt(m, j, i) / getElt(m, i, i);
            addition_multiplication(m, j, i, k);
        }
    }

}

void inversion_gauss(Matrix m) {

    unsigned int i, j;
    E k;

    if (!isSquare(m)) {
        fprintf(stderr, "La matrice n'est pas carrée\n");
        exit(EXIT_FAILURE);
    }

    Matrix id = matrix_identite(m->nb_rows);

    for (i = 0; i < m->nb_rows-1; i++) {
        for (j = i+1; j < m->nb_rows; j++) {
            k = -getElt(m, j, i) / getElt(m, i, i);
            addition_multiplication(m, j, i, k);
            addition_multiplication(id, j, i, k);
        }
    }



    printf("\n\n\nINVERSION DE LA MATRICE AVEC ALGO DE GAUSS\n");
    printMatrix(id);

    for (i = 0; i < m->nb_rows-1; i++) {
        multiplier_ligne(m, i, 1/getElt(m, i, i));
        multiplier_ligne(id, i, 1/getElt(id, i, i));
    }


    printf("\n\n\nINVERSION DE LA MATRICE AVEC ALGO DE GAUSS\n");
    printMatrix(id);

}

PLU decomposition_PLU(Matrix m) {
    Matrix P = matrix_identite(m->nb_rows);
    Matrix L = matrix_identite(m->nb_rows);
    Matrix U = newMatrix(m->nb_rows, m->nb_rows);
    int i, j, k, n = m->nb_rows;
    E l, p, temp, q;
    PLU m2;
    m2.P = P;
    m2.U = U;
    m2.L = L;
    copy_matrix(m, m2.U);

    for (k = 0; k < n; k++) {
        p = getElt(m2.U, k, k);
        l = k;
        for (i = k; i < n; i++) {
            if (abs(getElt(m2.U, i, k)) > p) {
                p = getElt(m2.U, i, k);
                l = i;
            }
        }
        if (l != k) {
            for (j = 0; j < n; j++) {
                temp = getElt(m2.U, k, j);
                setElt(m2.U, k, j, getElt(m2.U, l, j));
                setElt(m2.U, l, j, temp);
                if (j < k) {
                    temp = getElt(m2.L, k, j);
                    setElt(m2.L, k, j, getElt(m2.L, l, j));
                    setElt(m2.L, l, j, temp);
                    temp = getElt(m2.P, k, j);
                    setElt(m2.P, k, j, getElt(m2.P, l, j));
                    setElt(m2.P, l, j, temp);
                }
            }
            for (i = k + 2; i < n; i++){
                q = getElt(m2.U, i, k);
                setElt(m2.U, i, k, 0);
                setElt(m2.L, i, k, (q/p));
                for (j = k + 2; j < n; j++) {
                    setElt(m2.U, i, j, getElt(m2.U, i, j) - (getElt(m2.U, k, j) * (q/p)));
                }
            }
        }
    }

    return m2;
}

E * valeurs_propres(Matrix m) {
    E* val = malloc(2*sizeof(E));
    E b, c, delta;

    if (!isSquare(m)) {
        fprintf(stderr, "La matrice n'est pas carrée\n");
        exit(EXIT_FAILURE);
    }

    if (m->nb_rows != 2) {
        fprintf(stderr, "La matrice n'est pas de taille 2\n");
    }

    if (m->nb_rows == 2 && m->nb_columns == 2) {
        b = -(getElt(m,0,0) + getElt(m,1,1));
        c = m_determinant(m);
        delta = b*b - 4*c;
        if (delta < 0) {
            fprintf(stderr, "Les valeurs propres ne sont pas réelles\n");
            return NULL;
        } else if (delta == 0) {
            val[0] = -(b/2);
            val[1] = val[0];
        } else {
            val[0] = (-b- sqrtf(delta))/2;
            val[1] = (-b+ sqrtf(delta))/2;
        }
    }
    return val;
}
