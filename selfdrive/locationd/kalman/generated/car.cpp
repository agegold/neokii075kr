
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6983171995119894410) {
   out_6983171995119894410[0] = delta_x[0] + nom_x[0];
   out_6983171995119894410[1] = delta_x[1] + nom_x[1];
   out_6983171995119894410[2] = delta_x[2] + nom_x[2];
   out_6983171995119894410[3] = delta_x[3] + nom_x[3];
   out_6983171995119894410[4] = delta_x[4] + nom_x[4];
   out_6983171995119894410[5] = delta_x[5] + nom_x[5];
   out_6983171995119894410[6] = delta_x[6] + nom_x[6];
   out_6983171995119894410[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3289406534228174586) {
   out_3289406534228174586[0] = -nom_x[0] + true_x[0];
   out_3289406534228174586[1] = -nom_x[1] + true_x[1];
   out_3289406534228174586[2] = -nom_x[2] + true_x[2];
   out_3289406534228174586[3] = -nom_x[3] + true_x[3];
   out_3289406534228174586[4] = -nom_x[4] + true_x[4];
   out_3289406534228174586[5] = -nom_x[5] + true_x[5];
   out_3289406534228174586[6] = -nom_x[6] + true_x[6];
   out_3289406534228174586[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_5597874357973578331) {
   out_5597874357973578331[0] = 1.0;
   out_5597874357973578331[1] = 0.0;
   out_5597874357973578331[2] = 0.0;
   out_5597874357973578331[3] = 0.0;
   out_5597874357973578331[4] = 0.0;
   out_5597874357973578331[5] = 0.0;
   out_5597874357973578331[6] = 0.0;
   out_5597874357973578331[7] = 0.0;
   out_5597874357973578331[8] = 0.0;
   out_5597874357973578331[9] = 1.0;
   out_5597874357973578331[10] = 0.0;
   out_5597874357973578331[11] = 0.0;
   out_5597874357973578331[12] = 0.0;
   out_5597874357973578331[13] = 0.0;
   out_5597874357973578331[14] = 0.0;
   out_5597874357973578331[15] = 0.0;
   out_5597874357973578331[16] = 0.0;
   out_5597874357973578331[17] = 0.0;
   out_5597874357973578331[18] = 1.0;
   out_5597874357973578331[19] = 0.0;
   out_5597874357973578331[20] = 0.0;
   out_5597874357973578331[21] = 0.0;
   out_5597874357973578331[22] = 0.0;
   out_5597874357973578331[23] = 0.0;
   out_5597874357973578331[24] = 0.0;
   out_5597874357973578331[25] = 0.0;
   out_5597874357973578331[26] = 0.0;
   out_5597874357973578331[27] = 1.0;
   out_5597874357973578331[28] = 0.0;
   out_5597874357973578331[29] = 0.0;
   out_5597874357973578331[30] = 0.0;
   out_5597874357973578331[31] = 0.0;
   out_5597874357973578331[32] = 0.0;
   out_5597874357973578331[33] = 0.0;
   out_5597874357973578331[34] = 0.0;
   out_5597874357973578331[35] = 0.0;
   out_5597874357973578331[36] = 1.0;
   out_5597874357973578331[37] = 0.0;
   out_5597874357973578331[38] = 0.0;
   out_5597874357973578331[39] = 0.0;
   out_5597874357973578331[40] = 0.0;
   out_5597874357973578331[41] = 0.0;
   out_5597874357973578331[42] = 0.0;
   out_5597874357973578331[43] = 0.0;
   out_5597874357973578331[44] = 0.0;
   out_5597874357973578331[45] = 1.0;
   out_5597874357973578331[46] = 0.0;
   out_5597874357973578331[47] = 0.0;
   out_5597874357973578331[48] = 0.0;
   out_5597874357973578331[49] = 0.0;
   out_5597874357973578331[50] = 0.0;
   out_5597874357973578331[51] = 0.0;
   out_5597874357973578331[52] = 0.0;
   out_5597874357973578331[53] = 0.0;
   out_5597874357973578331[54] = 1.0;
   out_5597874357973578331[55] = 0.0;
   out_5597874357973578331[56] = 0.0;
   out_5597874357973578331[57] = 0.0;
   out_5597874357973578331[58] = 0.0;
   out_5597874357973578331[59] = 0.0;
   out_5597874357973578331[60] = 0.0;
   out_5597874357973578331[61] = 0.0;
   out_5597874357973578331[62] = 0.0;
   out_5597874357973578331[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_7992265079146681318) {
   out_7992265079146681318[0] = state[0];
   out_7992265079146681318[1] = state[1];
   out_7992265079146681318[2] = state[2];
   out_7992265079146681318[3] = state[3];
   out_7992265079146681318[4] = state[4];
   out_7992265079146681318[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_7992265079146681318[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_7992265079146681318[7] = state[7];
}
void F_fun(double *state, double dt, double *out_4136008500583111326) {
   out_4136008500583111326[0] = 1;
   out_4136008500583111326[1] = 0;
   out_4136008500583111326[2] = 0;
   out_4136008500583111326[3] = 0;
   out_4136008500583111326[4] = 0;
   out_4136008500583111326[5] = 0;
   out_4136008500583111326[6] = 0;
   out_4136008500583111326[7] = 0;
   out_4136008500583111326[8] = 0;
   out_4136008500583111326[9] = 1;
   out_4136008500583111326[10] = 0;
   out_4136008500583111326[11] = 0;
   out_4136008500583111326[12] = 0;
   out_4136008500583111326[13] = 0;
   out_4136008500583111326[14] = 0;
   out_4136008500583111326[15] = 0;
   out_4136008500583111326[16] = 0;
   out_4136008500583111326[17] = 0;
   out_4136008500583111326[18] = 1;
   out_4136008500583111326[19] = 0;
   out_4136008500583111326[20] = 0;
   out_4136008500583111326[21] = 0;
   out_4136008500583111326[22] = 0;
   out_4136008500583111326[23] = 0;
   out_4136008500583111326[24] = 0;
   out_4136008500583111326[25] = 0;
   out_4136008500583111326[26] = 0;
   out_4136008500583111326[27] = 1;
   out_4136008500583111326[28] = 0;
   out_4136008500583111326[29] = 0;
   out_4136008500583111326[30] = 0;
   out_4136008500583111326[31] = 0;
   out_4136008500583111326[32] = 0;
   out_4136008500583111326[33] = 0;
   out_4136008500583111326[34] = 0;
   out_4136008500583111326[35] = 0;
   out_4136008500583111326[36] = 1;
   out_4136008500583111326[37] = 0;
   out_4136008500583111326[38] = 0;
   out_4136008500583111326[39] = 0;
   out_4136008500583111326[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_4136008500583111326[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_4136008500583111326[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4136008500583111326[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4136008500583111326[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_4136008500583111326[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_4136008500583111326[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_4136008500583111326[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_4136008500583111326[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_4136008500583111326[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_4136008500583111326[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4136008500583111326[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4136008500583111326[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_4136008500583111326[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_4136008500583111326[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_4136008500583111326[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4136008500583111326[56] = 0;
   out_4136008500583111326[57] = 0;
   out_4136008500583111326[58] = 0;
   out_4136008500583111326[59] = 0;
   out_4136008500583111326[60] = 0;
   out_4136008500583111326[61] = 0;
   out_4136008500583111326[62] = 0;
   out_4136008500583111326[63] = 1;
}
void h_25(double *state, double *unused, double *out_7858839603282072417) {
   out_7858839603282072417[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3243248075051422251) {
   out_3243248075051422251[0] = 0;
   out_3243248075051422251[1] = 0;
   out_3243248075051422251[2] = 0;
   out_3243248075051422251[3] = 0;
   out_3243248075051422251[4] = 0;
   out_3243248075051422251[5] = 0;
   out_3243248075051422251[6] = 1;
   out_3243248075051422251[7] = 0;
}
void h_24(double *state, double *unused, double *out_434664235133717623) {
   out_434664235133717623[0] = state[4];
   out_434664235133717623[1] = state[5];
}
void H_24(double *state, double *unused, double *out_6950407459951144513) {
   out_6950407459951144513[0] = 0;
   out_6950407459951144513[1] = 0;
   out_6950407459951144513[2] = 0;
   out_6950407459951144513[3] = 0;
   out_6950407459951144513[4] = 1;
   out_6950407459951144513[5] = 0;
   out_6950407459951144513[6] = 0;
   out_6950407459951144513[7] = 0;
   out_6950407459951144513[8] = 0;
   out_6950407459951144513[9] = 0;
   out_6950407459951144513[10] = 0;
   out_6950407459951144513[11] = 0;
   out_6950407459951144513[12] = 0;
   out_6950407459951144513[13] = 1;
   out_6950407459951144513[14] = 0;
   out_6950407459951144513[15] = 0;
}
void h_30(double *state, double *unused, double *out_8558938070953574819) {
   out_8558938070953574819[0] = state[4];
}
void H_30(double *state, double *unused, double *out_542288563727397583) {
   out_542288563727397583[0] = 0;
   out_542288563727397583[1] = 0;
   out_542288563727397583[2] = 0;
   out_542288563727397583[3] = 0;
   out_542288563727397583[4] = 1;
   out_542288563727397583[5] = 0;
   out_542288563727397583[6] = 0;
   out_542288563727397583[7] = 0;
}
void h_26(double *state, double *unused, double *out_6698867777079951126) {
   out_6698867777079951126[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1744070034470562349) {
   out_1744070034470562349[0] = 0;
   out_1744070034470562349[1] = 0;
   out_1744070034470562349[2] = 0;
   out_1744070034470562349[3] = 0;
   out_1744070034470562349[4] = 0;
   out_1744070034470562349[5] = 0;
   out_1744070034470562349[6] = 0;
   out_1744070034470562349[7] = 1;
}
void h_27(double *state, double *unused, double *out_5004108522969749662) {
   out_5004108522969749662[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3356740395752664785) {
   out_3356740395752664785[0] = 0;
   out_3356740395752664785[1] = 0;
   out_3356740395752664785[2] = 0;
   out_3356740395752664785[3] = 1;
   out_3356740395752664785[4] = 0;
   out_3356740395752664785[5] = 0;
   out_3356740395752664785[6] = 0;
   out_3356740395752664785[7] = 0;
}
void h_29(double *state, double *unused, double *out_6068435327800607520) {
   out_6068435327800607520[0] = state[1];
}
void H_29(double *state, double *unused, double *out_2500405786000066171) {
   out_2500405786000066171[0] = 0;
   out_2500405786000066171[1] = 1;
   out_2500405786000066171[2] = 0;
   out_2500405786000066171[3] = 0;
   out_2500405786000066171[4] = 0;
   out_2500405786000066171[5] = 0;
   out_2500405786000066171[6] = 0;
   out_2500405786000066171[7] = 0;
}
void h_28(double *state, double *unused, double *out_6603633576743408955) {
   out_6603633576743408955[0] = state[5];
   out_6603633576743408955[1] = state[6];
}
void H_28(double *state, double *unused, double *out_7882107464910371989) {
   out_7882107464910371989[0] = 0;
   out_7882107464910371989[1] = 0;
   out_7882107464910371989[2] = 0;
   out_7882107464910371989[3] = 0;
   out_7882107464910371989[4] = 0;
   out_7882107464910371989[5] = 1;
   out_7882107464910371989[6] = 0;
   out_7882107464910371989[7] = 0;
   out_7882107464910371989[8] = 0;
   out_7882107464910371989[9] = 0;
   out_7882107464910371989[10] = 0;
   out_7882107464910371989[11] = 0;
   out_7882107464910371989[12] = 0;
   out_7882107464910371989[13] = 0;
   out_7882107464910371989[14] = 1;
   out_7882107464910371989[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
