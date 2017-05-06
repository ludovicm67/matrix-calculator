#include <stdio.h>

// Affiche un message d'erreur
void print_error(char * msg) {
    fprintf(stderr, "\033[1;31mERROR : %s\033[0m\n", msg);
}
