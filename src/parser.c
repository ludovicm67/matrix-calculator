#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "mpc.h"
#include "system.h"
#include "matrix.h"
#include "parser.h"

void print_expression(Expression e) {
    if (!e) print_error("Empty expression");

    switch (e->type) {
        case UNKNOWN:
            print_error("Unknown expression");
            break;
        case ERROR:
            print_error(e->c.str);
            break;
        case MATRIX:
            printMatrix(e->c.m);
            deleteMatrix(e->c.m);
            break;
        case SCALAR:
            printf("%f\n", e->c.s);
            break;
        case ASSIGN:
            printf("La variable '%s' vaut désormais :\n", e->c.a->symbol);
            print_expression(e->c.a->e);
            break;
        case IDENT:
            print_error("Undefined variable");
            break;
        case CALL:
            print_error("Undefined function");
            break;
        case NOTHING:
            break;
        default:
            print_error("Unknown expression type");
    }
}

Expression new_expression() {
    Expression e = malloc(sizeof(struct s_expression));
    if (e == NULL) {
        print_error("Impossible d'allouer de la mémoire !");
        exit(EXIT_FAILURE);
    }
    e->type = UNKNOWN;
    return e;
}

Expression new_expression_error(char * msg) {
    Expression e = new_expression();
    e->type = ERROR;
    e->c.str = msg;
    return e;
}

mpc_val_t* val_to_expr(mpc_val_t* val) {
    Expression e = new_expression();
    e->type = SCALAR;
    e->c.s = *(float *) val;
    free(val);
    return e;
}

mpc_val_t* ident_to_expr(int n, mpc_val_t ** xs) {
    (void) n;
    char * name = (char *) xs[0];
    assign env = (assign) xs[1];

    Expression e = new_expression();

    while (env) {
        if (!strcmp(env->symbol, name)) {
            if (e->type == MATRIX) deleteMatrix(e->c.m);
            memcpy(e, env->e, sizeof(struct s_expression));
            if (env->e->type == MATRIX) {
                e->c.m = new_matrix_copy(env->e->c.m);
            }
        }
        env = env->next;
    }

    if (e->type == UNKNOWN) e = new_expression_error("variable inconnue");


    return e;
}

mpc_val_t* call_to_expr(int n, mpc_val_t ** xs) {
    Expression e = new_expression();
    char * name = (char *) xs[0];
    Expression param = (Expression) xs[1];
    unsigned int has_param = param != NULL;

    if (n == 2) {
        if (has_param) {
            if (!strcmp(name, "id")) {
                if (param->type == SCALAR) {
                    e->type = MATRIX;
                    e->c.m = matrix_identite(param->c.s);
                }
            }

            else if (!strcmp(name, "det")) {
                if (param->type == MATRIX) {
                    if (!isSquare(param->c.m)) {
                        e->type = ERROR;
                        e->c.str = "La matrice doit être carrée !";
                    } else {
                        e->type = SCALAR;
                        e->c.s = det(param->c.m);
                    }
                }
            }

            else if (!strcmp(name, "det_tri")) {
                if (param->type == MATRIX) {
                    if (!isSquare(param->c.m)) {
                        e->type = ERROR;
                        e->c.str = "La matrice doit être carrée !";
                    } else {
                        e->type = SCALAR;
                        e->c.s = m_determinant(param->c.m);
                    }
                }
            }

            else if (!strcmp(name, "tr")) {
                if (param->type == MATRIX) {
                    e->type = MATRIX;
                    e->c.m = transpose(param->c.m);
                    deleteMatrix(param->c.m);
                }
            }

            else if (!strcmp(name, "inv")) {
                if (param->type == MATRIX) {
                    if (!isSquare(param->c.m)) {
                        e->type = ERROR;
                        e->c.str = "La matrice doit être carrée !";
                    } else {
                        e->c.m = inversion(param->c.m);
                        if (!e->c.m) {
                            e->type = ERROR;
                            e->c.str = "La matrice n'est pas inversible.";
                        } else e->type = MATRIX;
                    }
                    deleteMatrix(param->c.m);
                }
            }

            else if (!strcmp(name, "invg")) {
                if (param->type == MATRIX) {
                    if (!isSquare(param->c.m)) {
                        e->type = ERROR;
                        e->c.str = "La matrice doit être carrée !";
                    } else {
                        e->c.m = inversion_gauss(param->c.m);
                        if (!e->c.m) {
                            e->type = ERROR;
                            e->c.str = "La matrice n'est pas inversible.";
                        } else e->type = MATRIX;
                    }
                    deleteMatrix(param->c.m);
                }
            }

            else if (!strcmp(name, "plu")) {
                if (param->type == MATRIX) {
                    m_PLU(param->c.m);
                    e->type = NOTHING;
                    deleteMatrix(param->c.m);
                }
            }

            else if (!strcmp(name, "plu_p")) {
                if (param->type == MATRIX) {
                    e->c.m = m_PLU_p(param->c.m);
                    if (!e->c.m) {
                        e->type = ERROR;
                        e->c.str = "Aucune matrice P trouvée.";
                    } else e->type = MATRIX;
                    deleteMatrix(param->c.m);
                }
            }

            else if (!strcmp(name, "plu_l")) {
                if (param->type == MATRIX) {
                    e->c.m = m_PLU_l(param->c.m);
                    if (!e->c.m) {
                        e->type = ERROR;
                        e->c.str = "Aucune matrice P trouvée.";
                    } else e->type = MATRIX;
                    deleteMatrix(param->c.m);
                }
            }

            else if (!strcmp(name, "plu_u")) {
                if (param->type == MATRIX) {
                    e->c.m = m_PLU_u(param->c.m);
                    if (!e->c.m) {
                        e->type = ERROR;
                        e->c.str = "Aucune matrice P trouvée.";
                    } else e->type = MATRIX;
                    deleteMatrix(param->c.m);
                }
            }

        }

    }

    return e;
}


mpc_val_t *fold_sum(int n, mpc_val_t ** xs) {
    int i;
    Expression * e = (Expression *) xs;

    // valeur par défaut
    if (!n) {
        Expression e = new_expression();
        e->c.s = 0;
        return e;
    }

    for (i = 1; i < n; i++) {
        if (e[0]->type == MATRIX && e[i]->type == MATRIX) {
            if(sameSize(e[0]->c.m, e[i]->c.m)) {
                e[0]->c.m = addition(e[0]->c.m, e[i]->c.m);
            } else {
                e[0]->type = ERROR;
                e[0]->c.str = "Les matrices doivent être de même dimensions.";
            }
        } else if ((e[0]->type == MATRIX && e[i]->type == SCALAR) || (e[0]->type == SCALAR && e[i]->type == MATRIX)) {
            e[0]->type = ERROR;
            e[0]->c.str = "Impossible d'additioner un scalaire avec une matrice.";
        } else {
            e[0]->c.s += e[i]->c.s;
        }
        free(xs[i]);
    }

    return xs[0];
}

mpc_val_t *fold_prod(int n, mpc_val_t ** xs) {
    int i;
    Expression * e = (Expression *) xs;
    Matrix m;

    // valeur par défaut
    if (!n) {
        Expression e = new_expression();
        e->c.s = 1;
        return e;
    }

    for (i = 1; i < n; i++) {
        if (e[0]->type == MATRIX && e[i]->type == MATRIX) {
            e[0]->c.m = multiplication(e[0]->c.m, e[i]->c.m);
        } else if (e[0]->type == SCALAR && e[i]->type == MATRIX) {
            e[0]->type = MATRIX;
            e[0]->c.m = mult_scalar(e[0]->c.s, e[i]->c.m);
            deleteMatrix(e[i]->c.m);
        } else if (e[0]->type == MATRIX && e[i]->type == SCALAR) {
            m = mult_scalar(e[i]->c.s, e[0]->c.m);
            deleteMatrix(e[0]->c.m);
            e[0]->c.m = m;
        } else if (e[0]->type == SCALAR && e[i]->type == SCALAR) {
            e[0]->c.s *= e[i]->c.s;
        }
        free(xs[i]);
    }

    return xs[0];
}

mpc_val_t *fold_assign(int n, mpc_val_t ** xs) {
    char * name = (char *) xs[0];
    Expression e = (Expression) xs[2];

    Expression assign = new_expression();
    assign->type = ASSIGN;
    assign->c.a = malloc(sizeof(struct s_assign));
    assign->c.a->e = e;
    assign->c.a->symbol = name;

    (void) n;

    return assign;
}

mpc_val_t *fold_solve(int n, mpc_val_t ** xs) {
    Matrix a = ((Expression) xs[0])->c.m;
    Matrix b = ((Expression) xs[3])->c.m;
    Expression e = new_expression();
    e->c.m = multiplication(inversion(a), b);
    e->type = MATRIX;

    (void) n;

    return e;
}

mpc_val_t *fold_mat_row_first(int n, mpc_val_t ** xs) {
    (void) n;
    Expression head = (Expression) xs[0];
    Expression rest = (Expression) xs[1];

    if (head->type == SCALAR) {
        rest->c.mr.row[0] = head->c.s;
    }

    free(head);

    return rest;
}

mpc_val_t *fold_mat_row(int n, mpc_val_t ** xs) {
    int i;
    Expression row = new_expression();
    Expression * e = (Expression *) xs;
    row->type = MATRIX_ROW;
    if (!n) {
        row->c.mr.size = 1;
        row->c.mr.row = malloc(sizeof(float));
        return row;
    }
    row->c.mr.size = n+1;
    row->c.mr.row = malloc((n+1) * sizeof(float));
    for (i = 0; i < n; i++) {
        if (e[i]->type == SCALAR) {
            row->c.mr.row[i+1] = e[i]->c.s;
            free(xs[i]);
        }
    }
    return row;
}

mpc_val_t *fold_mat_first(int n, mpc_val_t ** xs) {
    (void) n;
    unsigned int i, j, max_cols = 0;
    Expression head = (Expression) xs[0];
    Expression rest = (Expression) xs[1];

    if (head->type == MATRIX_ROW) {
        rest->c.mra.raw[0] = head->c.mr;
    }

    for (i = 0; i < rest->c.mra.size; i++) {
        if (rest->c.mra.raw[i].size > max_cols) {
            max_cols = rest->c.mra.raw[i].size;
        }
    }

    Matrix m = newMatrix(rest->c.mra.size, max_cols);
    for (i = 0; i < rest->c.mra.size; i++) {
        for (j = 0; j < rest->c.mra.raw[i].size; j++) {
            setElt(m, i, j, rest->c.mra.raw[i].row[j]);
        }
        free(rest->c.mra.raw[i].row);
    }
    free(rest->c.mra.raw);

    free(head);

    rest->type = MATRIX;
    rest->c.m = m;
    return rest;
}

mpc_val_t *fold_mat(int n, mpc_val_t ** xs) {
    int i;
    Expression raw = new_expression();
    Expression * e = (Expression *) xs;
    raw->type = MATRIX_RAW;
    if (!n) {
        raw->c.mra.size = 1;
        raw->c.mra.raw = malloc(sizeof(struct s_matrix_row));
        return raw;
    }
    raw->c.mra.size = n+1;
    raw->c.mra.raw = malloc((n+1) * sizeof(struct s_matrix_row));
    for (i = 0; i < n; i++) {
        if (e[i]->type == MATRIX_ROW) {
            raw->c.mra.raw[i+1] = e[i]->c.mr;
            free(xs[i]);
        }
    }
    return raw;
}


mpc_val_t *fold_value(int n, mpc_val_t ** xs) {

    char * op = (char *) xs[0];
    Expression e = (Expression) xs[1];

    switch (*op) {
        case '-':
            // Si l'opérateur était un -, on multiplie par -1
            if (e->type == SCALAR) e->c.s *= -1;
            else if (e->type == MATRIX) {
                mult_scalar(-1, e->c.m);
            }
            break;
        case '/':
            // Si l'opérateur était un /, on  prend l'inverse
            if (e->type == SCALAR) e->c.s = 1 / e->c.s;
            else if (e->type == MATRIX) {
                e->c.m = inversion(e->c.m);
            }
            break;
        default:
            break;
    }

    free(xs[0]);

    (void) n;

    return xs[1];
}


void catch_segfault(int signum) {
    (void) signum;
    printf("\n");
    print_error("Problème de mémoire ! Arrêt du programme.");
    exit(EXIT_FAILURE);
}

void free_env(assign env) {
    if (!env) return;
    if (env->next) free_env(env->next);
    if (env->e) {
        if (env->e->type == MATRIX) deleteMatrix(env->e->c.m);
        free(env->e);
    }
    free(env);
}


void run_parser() {
    signal(SIGSEGV, catch_segfault);

    assign environnement = malloc(sizeof(struct s_assign));
    environnement->symbol = "pi";
    environnement->e = new_expression();
    environnement->e->type = SCALAR;
    environnement->e->c.s = 3.141593;
    environnement->next = NULL;

    Expression e;

    int is_tty = isatty(0);
    if (is_tty) printf("\033[1mBonjour !\033[0m\n"); // convivialité !

    mpc_parser_t *Expr     = mpc_new("expression");
    mpc_parser_t *Prod     = mpc_new("product");
    mpc_parser_t *Constant = mpc_new("constant");
    mpc_parser_t *Value    = mpc_new("value");
    mpc_parser_t *Line     = mpc_new("line");
    mpc_parser_t *Input    = mpc_new("input");
    mpc_parser_t *Ident    = mpc_new("ident");
    mpc_parser_t *Assign   = mpc_new("assign");
    mpc_parser_t *Mat      = mpc_new("mat");
    mpc_parser_t *MatRow   = mpc_new("mat-row");
    mpc_parser_t *Row      = mpc_new("row");
    mpc_parser_t *Call     = mpc_new("call");
    mpc_parser_t *Solve    = mpc_new("solve");

    mpc_define(Ident, mpc_ident());

    mpc_define(Constant, mpc_apply(mpc_re("-?([0-9]*\\.[0-9]+|[0-9]+\\.?)"), mpcf_float));

    mpc_define(Expr, mpc_and(2, fold_sum,
        Prod, mpc_many(fold_sum, mpc_and(2, fold_value,
            mpc_oneof("+-"), Prod,
            free
        )),
        free
    ));

    mpc_define(Prod, mpc_and(2, fold_prod,
        Value, mpc_many(fold_prod, mpc_and(2, fold_value,
            mpc_oneof("*/"), Value,
            free
        )),
        free
    ));

    mpc_define(Assign, mpc_and(3, fold_assign,
        Ident, mpc_strip(mpc_char('=')), Expr,
        free, free
    ));

    mpc_define(MatRow, mpc_and(2, fold_mat_row_first,
        Expr, mpc_many(fold_mat_row, mpc_and(2, mpcf_snd,
            mpc_char(','), Expr,
            free
        )),
        free
    ));

    mpc_define(Mat, mpc_tok_squares(mpc_and(2, fold_mat_first,
        MatRow, mpc_many(fold_mat, mpc_and(2, mpcf_snd,
            mpc_char(';'), MatRow,
            free
        )),
        free
    ), free));

    mpc_define(Call, mpc_and(2, call_to_expr,
        Ident,
        mpc_parens(mpc_maybe(Expr), free),
        free
    ));

    mpc_define(Value, mpc_strip(mpc_or(5,
        Call,
        mpc_and(2, ident_to_expr, Ident, mpc_lift_val(environnement), free),
        mpc_apply(Constant, val_to_expr),
        Mat,
        mpc_parens(Expr, free)
    )));

    mpc_define(Solve, mpc_and(4, fold_solve,
        Mat, mpc_char('X'), mpc_strip(mpc_char('=')), Mat,
        free, free, free
    ));

    mpc_define(Line, mpc_strip(mpc_or(3,
        Solve, Assign, Expr
    )));

    mpc_define(Input, mpc_whole(Line, free));

    mpc_optimise(Expr);
    mpc_optimise(Prod);
    mpc_optimise(Value);
    mpc_optimise(Call);
    mpc_optimise(Ident);
    mpc_optimise(Constant);
    mpc_optimise(Assign);
    mpc_optimise(Line);
    mpc_optimise(Input);
    mpc_optimise(Row);
    mpc_optimise(Mat);
    mpc_optimise(MatRow);
    mpc_optimise(Solve);


    mpc_result_t r;

    // pour le getline
    char * line = NULL;
    size_t len = 0;

    if (is_tty) printf("\033[1;34m>>> \033[0m");
    while ((getline(&line, &len, stdin)) != -1) {
        line[strcspn(line, "\r\n#")] = 0;

        if (strlen(line) > 0) {
            if (mpc_parse("input", line, Input, &r)) {
                e = (Expression) r.output;
                if (e->type == ASSIGN) {
                    assign new_assign = malloc(sizeof(struct s_assign));
                    new_assign->symbol = e->c.a->symbol;
                    new_assign->e = new_expression();
                    memcpy(new_assign->e, e->c.a->e, sizeof(struct s_expression));
                    if (e->c.a->e->type == MATRIX) {
                        new_assign->e->c.m = new_matrix_copy(e->c.a->e->c.m);
                    }
                    new_assign->next = NULL;
                    environnement->next = new_assign;
                }
                print_expression(r.output);
                free(r.output);
            } else {
                if (!is_tty) fprintf(stderr, "%s\n", line);
                printf("%*s", (int) (is_tty ? r.error->state.col+4 : r.error->state.col), "");
                printf("\033[1;31m^\033[0m\n");
                print_error(err_msg_only(r.error));
                mpc_err_delete(r.error);
            }
        }

        if (is_tty) printf("\033[1;34m>>> \033[0m");
    }

    if (line) free(line);

    mpc_cleanup(13, Assign, Call, Constant, Ident, Expr, Prod, Value, Line, Input, Row, Mat, MatRow, Solve);

    free_env(environnement);

    if (is_tty) printf("\n\033[1mAu revoir ! :)\033[0m\n"); // convivialité !
}
