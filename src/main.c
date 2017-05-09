#include <stdlib.h>
#include "parser.h"

int main() {
    Matrix A = newMatrix(3, 3);
    setElt(A, 0, 0, 0);
    setElt(A, 0, 1, 3);
    setElt(A, 0, 2, 1);
    setElt(A, 1, 0, 4);
    setElt(A, 1, 1, 1);
    setElt(A, 1, 2, 1);
    setElt(A, 2, 0, 2);
    setElt(A, 2, 1, 2);
    setElt(A, 2, 2, 4);
    PLU decompo = decomposition_PLU(A);
    printMatrix(decompo.P);
    printMatrix(multiplication(decompo.P, A));
    printMatrix(multiplication(decompo.L, decompo.U));
    run_parser();
    return EXIT_SUCCESS;
}
