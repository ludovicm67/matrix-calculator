#ifndef __PARSER_H__
#define __PARSER_H__

#include "mpc.h"
#include "matrix.h"

typedef struct s_assign {
    char * symbol;
    struct s_expression * e;
    struct s_assign * next;
} * assign;

// ligne d'une matrice
typedef struct s_matrix_row {
    unsigned int size;
    float * row;
} matrix_row;

// structure de base permettant de générer ensuite une vraie matrice
typedef struct s_matrix_raw {
    unsigned int size;
    struct s_matrix_row * raw;
} matrix_raw;

typedef struct s_expression {
    enum {
        UNKNOWN = 0,
        ERROR,
        MATRIX,
        MATRIX_ROW,
        MATRIX_RAW,
        SCALAR,
        ASSIGN,
        IDENT,
        CALL
    } type;
    union {
        Matrix m;
        matrix_row mr;
        matrix_raw mra;
        float s;
        assign a;
        char * str;
    } c;
} * Expression;

void print_expression(Expression e);
Expression new_expression();
Expression new_expression_error(char * msg);
mpc_val_t* val_to_expr(mpc_val_t* val);
mpc_val_t* ident_to_expr(int n, mpc_val_t ** xs);
mpc_val_t *fold_sum(int n, mpc_val_t ** xs);
mpc_val_t *fold_prod(int n, mpc_val_t ** xs);
mpc_val_t *fold_assign(int n, mpc_val_t ** xs);
mpc_val_t *fold_mat_row_first(int n, mpc_val_t ** xs);
mpc_val_t *fold_mat_row(int n, mpc_val_t ** xs);
mpc_val_t *fold_mat_first(int n, mpc_val_t ** xs);
mpc_val_t *fold_mat(int n, mpc_val_t ** xs);
mpc_val_t *fold_value(int n, mpc_val_t ** xs);
void catch_segfault(int signum);
void run_parser();

#endif
