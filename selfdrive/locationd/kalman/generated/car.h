/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6983171995119894410);
void inv_err_fun(double *nom_x, double *true_x, double *out_3289406534228174586);
void H_mod_fun(double *state, double *out_5597874357973578331);
void f_fun(double *state, double dt, double *out_7992265079146681318);
void F_fun(double *state, double dt, double *out_4136008500583111326);
void h_25(double *state, double *unused, double *out_7858839603282072417);
void H_25(double *state, double *unused, double *out_3243248075051422251);
void h_24(double *state, double *unused, double *out_434664235133717623);
void H_24(double *state, double *unused, double *out_6950407459951144513);
void h_30(double *state, double *unused, double *out_8558938070953574819);
void H_30(double *state, double *unused, double *out_542288563727397583);
void h_26(double *state, double *unused, double *out_6698867777079951126);
void H_26(double *state, double *unused, double *out_1744070034470562349);
void h_27(double *state, double *unused, double *out_5004108522969749662);
void H_27(double *state, double *unused, double *out_3356740395752664785);
void h_29(double *state, double *unused, double *out_6068435327800607520);
void H_29(double *state, double *unused, double *out_2500405786000066171);
void h_28(double *state, double *unused, double *out_6603633576743408955);
void H_28(double *state, double *unused, double *out_7882107464910371989);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
