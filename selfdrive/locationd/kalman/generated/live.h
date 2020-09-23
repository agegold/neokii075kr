/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_8733412135913564155);
void inv_err_fun(double *nom_x, double *true_x, double *out_8163023127837648523);
void H_mod_fun(double *state, double *out_6055917334837167653);
void f_fun(double *state, double dt, double *out_8020493960037440266);
void F_fun(double *state, double dt, double *out_7540607979955595376);
void h_3(double *state, double *unused, double *out_6724944771315514963);
void H_3(double *state, double *unused, double *out_3669764119120110885);
void h_4(double *state, double *unused, double *out_4154834335693279250);
void H_4(double *state, double *unused, double *out_4947678365094884967);
void h_9(double *state, double *unused, double *out_8434469681780947925);
void H_9(double *state, double *unused, double *out_447942525258091791);
void h_10(double *state, double *unused, double *out_4250368256891693028);
void H_10(double *state, double *unused, double *out_1691537069964150712);
void h_12(double *state, double *unused, double *out_5379999212265071567);
void H_12(double *state, double *unused, double *out_1484137471320515425);
void h_13(double *state, double *unused, double *out_9126168919793352425);
void H_13(double *state, double *unused, double *out_8674164924223822620);
void h_14(double *state, double *unused, double *out_8434469681780947925);
void H_14(double *state, double *unused, double *out_447942525258091791);
void h_19(double *state, double *unused, double *out_8436829164367959717);
void H_19(double *state, double *unused, double *out_9190153002668699627);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);