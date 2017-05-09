#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "system.h"
#include "matrix.h"


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

Matrix new_matrix_copy(Matrix m) {
    Matrix r = newMatrix(m->nb_rows, m->nb_columns);
    copy_matrix(m, r);
    return r;
}

// Teste si une matrice est carré
int isSquare(Matrix m) {
    return m->nb_rows == m->nb_columns;
}

// Teste si deux matrices sont de même dimension
int sameSize(Matrix a, Matrix b) {
    return a->nb_rows == b->nb_rows && a->nb_columns == b->nb_columns;
}

// Teste si une matrice est triangulaire
int isTriangulaire(Matrix m) {
    unsigned int ligne, col;
    for (ligne = 1; ligne < m->nb_rows; ligne++) {
        for (col = 0; col < ligne; col++) {
            if (getElt(m, ligne, col) != 0) return 0;
        }
    }
    return 1;
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
            if (i >= row) x++;
            if (j >= column) y++;
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

    if (!det(m)) return NULL;

    Matrix inverse = newMatrix(m->nb_rows, m->nb_columns);
    for (i = 0; i < (int) m->nb_rows; i++) {
        for (j = 0; j < (int) m->nb_columns; j++) {
            if ((i + j) % 2 == 0) {
                setElt(inverse, i, j, det(extraction(m, i, j)));
            } else {
                setElt(inverse, i, j, -det(extraction(m, i, j)));
            }
        }
    }
    inverse = transpose(inverse);
    inverse = mult_scalar((1/det(m)), inverse);
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

// Une fonction permettant d’additionner à la ligne i le résultat de la multiplication de la ligne j par un facteur k
void addition_multiplication(Matrix m, unsigned int i, unsigned int j, E k) {
    unsigned int n;
    Matrix m1, m2 = newMatrix(m->nb_rows, m->nb_columns);
    for (n = 0; n < m->nb_columns; n++){
        setElt(m2, i, n, getElt(m, j, n));
    }
    multiplier_ligne(m2, i, k);
    m1 = addition(m, m2);
    copy_matrix(m1, m);

    deleteMatrix(m1);
    deleteMatrix(m2);
}

Matrix triangulariser(Matrix in) {
    Matrix m = new_matrix_copy(in);

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

    return m;
}

// calcul du déterminant avec matrice triangulaire sup.
// prec : m doit être carrée
E m_determinant(Matrix m) {
    unsigned int i;
    E determinant = 1;

    if (!isTriangulaire(m)) {
        triangulariser(m);
    }

    for (i = 0; i < m->nb_rows; i++) {
        determinant *= getElt(m, i, i);
    }

    return determinant;
}

// prec : La matrice doit être carrée
Matrix inversion_gauss(Matrix m) {

    unsigned int i, j;
    E k;

    if (!m_determinant(m)) return NULL;

printf("Matrice entée :\n");
printMatrix(m);

    Matrix id = matrix_identite(m->nb_rows);
    Matrix r = new_matrix_copy(m);

    // On commence par triangulariser la matice
    for (i = 0; i < r->nb_rows-1; i++) {
        for (j = i+1; j < r->nb_rows; j++) {
            k = -getElt(r, j, i) / getElt(r, i, i);
            addition_multiplication(r, j, i, k);
            addition_multiplication(id, j, i, k);
        }
    }


    k = 1/getElt(r, r->nb_rows-1, r->nb_rows-1);
    multiplier_ligne(r, r->nb_rows-1, k);
    multiplier_ligne(id, r->nb_rows-1, k);

printf("Matrice après trig :\n");
printMatrix(m);
printf("Matrice id après trig :\n");
printMatrix(id);
    // addition_multiplication(r, r->nb_rows-1, r->nb_rows-1, k);
    // addition_multiplication(id, r->nb_rows-1, r->nb_rows-1, k);

    // if (r->nb_rows >= 2) {
    //     for (j = r->nb_rows-1; (int) j >= 0; j--) {
    //         for (i = r->nb_rows-2; (int) i >= 0; i--) {
    //             k = -getElt(r, i, j);
    //             addition_multiplication(r, i, j, k);
    //             addition_multiplication(id, i, j, k);
    //             printf("%d %d %f\n", i, j, k);
    //         }
    //         k = 1/getElt(r, j, j);
    //         multiplier_ligne(r, j, k);
    //         multiplier_ligne(id, j, k);
    //     }
    // }

    // j = r->nb_rows-1;
    // for (i = r->nb_rows-2; (int) i >= 0; i--) {
    //     k = -getElt(r, i, j);
    //     addition_multiplication(r, i, j, k);
    //     addition_multiplication(id, i, j, k);
    //     printf("%d %d %f\n", i, j, k);
    // }
    // k = 1/getElt(r, j, j);
    // multiplier_ligne(r, j, k);
    // multiplier_ligne(id, j, k);
    // j = r->nb_rows-2;
    // for (i = r->nb_rows-2; (int) i >= 0; i--) {
    //     k = -getElt(r, i, j);
    //     addition_multiplication(r, j, i, k);
    //     addition_multiplication(id, i, j, k);
    //     printf("%d %d %f\n", i, j, k);
    // }
    // k = 1/getElt(r, j, j);
    // multiplier_ligne(r, j, k);
    // multiplier_ligne(id, j, k);


    printf("\n\n\nINVERSION DE LA MATRICE AVEC ALGO DE GAUSS\n");
    printMatrix(r);

    return id;

}

PLU decomposition_PLU(Matrix m) {
    unsigned int i, j, k, l;
    E p = 0, q = 0, tmp = 0;

    PLU M;
    Matrix P = matrix_identite(m->nb_rows);
    Matrix L = matrix_identite(m->nb_rows);
    Matrix U = m;

    for(k = 0; k < U->nb_rows; k++) {
        p = getElt(U, k, k);
        l = k;
        for(i = k; i < U->nb_rows; i++) {
            if(abs(getElt(U, i, k)) > p) {
                p = getElt(U, i, k);
                l = i;
            }
        }
        if(l != k) {
            for(j = 0; j < U->nb_rows; j++) {
                tmp = getElt(U, k, j);
                setElt(U, k, j, getElt(U, l, j));
                setElt(U, l, j, tmp);
                if(j < k) {
                    tmp = getElt(L, k, j);
                    setElt(L, k, j, getElt(L, l, j));
                    setElt(L, l, j, tmp);
                }
                tmp = getElt(P, k, j);
                setElt(P, k, j, getElt(P, l, j));
                setElt(P, l, j, tmp);
            }
        }
        for(i = k + 1; i < U->nb_rows; i++) {
            q = getElt(U, i, k);
            setElt(U, i, k, 0);
            setElt(L, i, k, (q / p));
            for(j = k + 1; j < U->nb_rows; j++) {
                setElt(U, i, j, getElt(U, i, j) - getElt(U, k, j) * (q / p));
            }
        }
    }

    M.P = P;
    M.L = L;
    M.U = U;

    return M;
}

void m_PLU(Matrix m) {
    PLU p = decomposition_PLU(m);
    printf("Matrice P :\n");
    printMatrix(p.P);
    printf("Matrice L :\n");
    printMatrix(p.L);
    printf("Matrice U :\n");
    printMatrix(p.U);
}

Matrix m_PLU_p(Matrix m) {
    PLU p = decomposition_PLU(m);
    return p.P;
}

Matrix m_PLU_l(Matrix m) {
    PLU p = decomposition_PLU(m);
    return p.L;
}

Matrix m_PLU_u(Matrix m) {
    PLU p = decomposition_PLU(m);
    return p.U;
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
