// C++ MEX code for SOAV minimization with Linear and Box Constraints
// Mitsuru Toyoda (Tokyo Metropolitan University)
// 2020/**/**
//
// #define PDOptimalityCriteria_ON

#include <Eigen/Cholesky>
// #include <Eigen/QR>
// #include <Eigen/LU>
#include <Eigen/Core>
#include <chrono>  // computation time
#include <iostream>

#include "mex.h"

template <typename Derived>
double interp1(                           //
    const Eigen::MatrixBase<Derived>& x,  //
    const Eigen::MatrixBase<Derived>& v,  //
    double xq,                            //
    int n_point) {
  // extrapolation
  if (xq <= x(0)) {
    return v(0);
  } else if (xq >= x(n_point - 1)) {
    return v(n_point - 1);
  }
  //  start bisection search
  int iL = 0;
  int iR = n_point - 1;
  int iM;
  while (iR - iL > 1) {
    iM = int(std::floor(double(iL + iR) / 2.0));
    if (x(iM) > xq) {
      iR = iM;
    } else if (x(iM) < xq) {
      iL = iM;
    } else {
      return v(iM);
    }
  }
  return (v(iR) * (xq - x(iL)) + v(iL) * (x(iR) - xq)) / (x(iR) - x(iL));
}

template <typename Derived>
int find_i_ubar(const Eigen::MatrixBase<Derived>& x, double xq) {
  // find an index "iubar" s.t. x(iubar-1) < xq <= x(iubar) by bisection search
  int iL = 0;
  int iR = int(x.size()) - 1;
  int iM;
  while (iR - iL > 1) {
    iM = int(std::floor(double(iL + iR) / 2.0));
    if (x(iM) >= xq) {
      iR = iM;
    } else {  // x(iM) < xq
      iL = iM;
    }
  }
  return iL + 1;  // iL = iubar -1
}

template <typename Derived>
int find_i_bar(const Eigen::MatrixBase<Derived>& x, double xq) {
  // find an index "ibar" s.t. x(ibar) <= xq < x(ibar+1) by bisection search
  int iL = 0;
  int iR = int(x.size()) - 1;
  int iM;
  while (iR - iL > 1) {
    iM = int(std::floor(double(iL + iR) / 2.0));
    if (x(iM) > xq) {
      iR = iM;
    } else {  // x(iM) <= xq
      iL = iM;
    }
  }
  return iL;  // iL = ibar
}

// soav_bisec(u,w,x_ubar,x_bar,Aeq,beq,Q,c,gamma,y0,z0);
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  // start timer for computation time
  std::chrono::system_clock::time_point t_start, t_end;
  t_start = std::chrono::system_clock::now();
  // check the number of arguments
  if (nrhs != 12) {
    mexErrMsgIdAndTxt("MATLAB:soav_interp:nrhs", "the number of arguments are 12. ");
  }
  // exit flag
  int exitflag = 0;
  // obtain size
  size_t n_w = mxGetM(prhs[0]);
  size_t n_x = mxGetN(prhs[4]);  // number of columns, i.e., width
  size_t n_b = mxGetM(prhs[4]);  // number of rows, i.e., height
  // read given arguments
  Eigen::Map<Eigen::MatrixXd> umat(mxGetPr(prhs[0]), n_w, n_x);
  Eigen::Map<Eigen::VectorXd> w(mxGetPr(prhs[1]), n_w);
  Eigen::Map<Eigen::VectorXd> x_ubar_table(mxGetPr(prhs[2]), n_x);
  Eigen::Map<Eigen::VectorXd> x_bar_table(mxGetPr(prhs[3]), n_x);
  Eigen::Map<Eigen::MatrixXd> A(mxGetPr(prhs[4]), n_b, n_x);
  Eigen::Map<Eigen::VectorXd> b(mxGetPr(prhs[5]), n_b);
  Eigen::Map<Eigen::MatrixXd> Q(mxGetPr(prhs[6]), n_x, n_x);
  Eigen::Map<Eigen::VectorXd> c(mxGetPr(prhs[7]), n_x);
  double gamma = mxGetScalar(prhs[8]);
  Eigen::Map<Eigen::VectorXd> y0(mxGetPr(prhs[9]), n_x);
  Eigen::Map<Eigen::VectorXd> z0(mxGetPr(prhs[10]), n_x);
#ifdef PDOptimalityCriteria_ON  // [eps_rel; eps_abs]
  Eigen::Map<Eigen::VectorXd> StoppingCriteriaPara(mxGetPr(prhs[11]), 2);
#else  // [ObjValTrue; ObjValTol; ConstraintTol]
  Eigen::Map<Eigen::VectorXd> StoppingCriteriaPara(mxGetPr(prhs[11]), 3);
#endif

  // preparation for matrix multiplication
  Eigen::MatrixXd Inx = Eigen::MatrixXd::Identity(n_x, n_x);
  Eigen::MatrixXd Q_gamma = (1. / gamma) * Inx + Q;
  Eigen::MatrixXd bmat = Q_gamma.llt().solve(((A * Q_gamma.llt().solve(A.transpose())).llt().solve(A)).transpose());
  Eigen::MatrixXd S_gamma = Q_gamma.llt().solve(Inx - A.transpose() * bmat.transpose()) * (1. / gamma);
  Eigen::VectorXd b_gamma = bmat * b;
  Eigen::VectorXd c_gamma = S_gamma * gamma * c;
  // variables for sort
  Eigen::VectorXi i_ubar_table(n_x);
  Eigen::VectorXi i_bar_table(n_x);
  Eigen::VectorXd sigma(n_w + 2);
  // state variables
  Eigen::VectorXd xk(n_x);
  Eigen::VectorXd yk(n_x);
  Eigen::VectorXd zk(n_x);
  Eigen::VectorXd ykp1(n_x);
  Eigen::VectorXd zkp1(n_x);
  yk = y0;
  zk = z0;

  // variables for output of function
  int MaxIter = int(1e8);
  int n_iter = MaxIter;
  double t_iter = INFINITY;
  double myInf = 1e20;
  double fval = myInf, ConstraintViolation = myInf;
  Eigen::VectorXd ConstraintViolationVecX(n_b);

  // variables for stopping criteria
#ifdef PDOptimalityCriteria_ON
  double rkp1_norm, skp1_norm;
  double eps_pri, eps_dual;
  double eps_rel = StoppingCriteriaPara(0);
  double eps_abs = StoppingCriteriaPara(1);
#else
  double ObjValTrue = StoppingCriteriaPara(0);
  double ObjValTol = StoppingCriteriaPara(1);
  double ConstraintTol = StoppingCriteriaPara(2);
  double fkX = myInf, fkp1Y = myInf;
  Eigen::VectorXd ConstraintViolationVecY(2 * n_x);
  // std::cout << "ObjValTrue = " <<  ObjValTrue << " ConstraintTol = " << ConstraintTol  << std::endl;
#endif

  // define permutation matrix
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> sorted_ind_mat(n_w + 2);
  size_t n_v_max = 2 * (n_w - 0 + 1) + 2;  // maximum width of table
  Eigen::VectorXd u_j(n_w + 2), w_ex(n_w + 2);
  Eigen::MatrixXd v_mat(n_v_max, n_x);
  v_mat.fill(0);
  Eigen::VectorXd v_mat_j(n_v_max);
  Eigen::VectorXi n_point(n_x);  // number of table points for each "x_j"
  Eigen::MatrixXd x_mat(n_v_max, n_x);
  x_mat.fill(0);
  w_ex << 0, w, 0;
  for (int j = 0; j < n_x; j++) {
    sorted_ind_mat.setIdentity();  // make diagonal elements ones
    u_j << -myInf, umat.col(j), myInf;
    std::sort(sorted_ind_mat.indices().data(), sorted_ind_mat.indices().data() + (n_w + 2), [&u_j](size_t i1, size_t i2) { return u_j(i1) < u_j(i2); });
    w_ex = sorted_ind_mat.transpose() * w_ex.eval();  // transpose and multiply from right side
    u_j = sorted_ind_mat.transpose() * u_j.eval();
    i_ubar_table(j) = find_i_ubar(u_j, x_ubar_table(j));
    i_bar_table(j) = find_i_ubar(u_j, x_bar_table(j));
    n_point(j) = 2 * (i_bar_table(j) - i_ubar_table(j) + 1) + 2;  // table width for index "j"
    for (int i = i_ubar_table(j) - 1; i < i_bar_table(j) + 1; i++) {
      sigma(i) = w_ex.head(i + 1).sum() - w_ex.tail(n_w + 1 - i).sum();
    }
    //
    Eigen::MatrixXd mat_xtmp(2, i_bar_table(j) - i_ubar_table(j) + 1);      // temporary matrix to build table X
    Eigen::MatrixXd mat_vtmp(2, i_bar_table(j) - i_ubar_table(j) + 1 + 1);  // temporary matrix to build table V
    //
    mat_xtmp = u_j.segment(i_ubar_table(j), i_bar_table(j) - i_ubar_table(j) + 1).transpose().replicate(2, 1);
    Eigen::Map<Eigen::VectorXd> mat_xtmp_v(mat_xtmp.data(), mat_xtmp.size());
    // calculate x(:,j)
    x_mat(0, j) = x_ubar_table(j);
    for (int i = 0; i < n_point(j) - 2; i++) {
      x_mat(i + 1, j) = mat_xtmp_v(i);
    }
    x_mat(n_point(j) - 1, j) = x_bar_table(j);
    // calculate v(:,j)
    mat_vtmp = sigma.segment(i_ubar_table(j) - 1, i_bar_table(j) - i_ubar_table(j) + 1 + 1).transpose().replicate(2, 1);
    Eigen::Map<Eigen::VectorXd> mat_vtmp_v(mat_vtmp.data(), mat_vtmp.size());
    for (int i = 0; i < n_point(j); i++) {
      v_mat(i, j) = x_mat(i, j) + gamma * mat_vtmp_v(i);
    }
  }
  //     std::cout << x_mat << std::endl; std::cout << v_mat << std::endl;
  for (int k = 0; k < MaxIter; k++) {
    // x[k]: proximal operation for $f$
    for (int j = 0; j < n_x; j++) {
      xk(j) = interp1(v_mat.col(j).head(n_point(j)), x_mat.col(j).head(n_point(j)), yk(j) - zk(j), n_point(j));
    }

#ifndef PDOptimalityCriteria_ON
    // objective function value and constraint violation
    fkX = ((-umat).transpose().colwise() + xk).cwiseAbs().colwise().sum() * w;
    fkX = fkX + (1. / 2.) * xk.transpose() * Q * xk + c.transpose() * xk;
    ConstraintViolationVecX = (b - A * xk).cwiseAbs();
    // std::cout << "X -> k: " << k << " ObjValX = " <<  fkX << " ConstraintViolationX = " << ConstraintViolationVecX.maxCoeff() << std::endl;
    if (ConstraintViolationVecX.maxCoeff() <= ConstraintTol && std::abs(fkX - ObjValTrue) <= ObjValTol) {
      // std::cout << "X -> k: " << k << " ObjValX = " << fkX << " ConstraintViolationX = " << ConstraintViolationVecX.maxCoeff() << std::endl;
      n_iter = k;
      exitflag = 1;
      fval = fkX;
      ConstraintViolation = ConstraintViolationVecX.maxCoeff();
      break;
    }
#endif

    // y[k]; proximal operation for $g$
    ykp1 = S_gamma * (xk + zk) + b_gamma - c_gamma;  // y(k+1)
    zkp1 = zk + xk - ykp1;                           // z(k+1)

#ifdef PDOptimalityCriteria_ON
    rkp1_norm = (xk - ykp1).norm();                  // x(k+1) - y(k+1)
    skp1_norm = (-1. / gamma * (ykp1 - yk)).norm();  // y(k+1) - y(k)
    eps_pri = sqrt(n_x) * eps_abs + eps_rel * std::max(xk.norm(), ykp1.norm());
    eps_dual = sqrt(n_x) * eps_abs + eps_rel * (zkp1 / gamma).norm();
    // std::cout << "k+1: " << k+1 << " eps_pri = " <<  eps_pri << " eps_dual = " <<  eps_dual << std::endl;
    if (rkp1_norm <= eps_pri && skp1_norm <= eps_dual) {
      n_iter = k + 1;
      exitflag = 1;
      fval = ((-umat).transpose().colwise() + xk).cwiseAbs().colwise().sum() * w;
      fval = fval + (1. / 2.) * xk.transpose() * Q * xk + c.transpose() * xk;
      ConstraintViolationVecX = (b - A * xk).cwiseAbs();
      ConstraintViolation = ConstraintViolationVecX.maxCoeff();
      break;
    }
#else
    // objective function value and constraint violation
    fkp1Y = ((-umat).transpose().colwise() + ykp1).cwiseAbs().colwise().sum() * w;
    fkp1Y = fkp1Y + (1. / 2.) * ykp1.transpose() * Q * ykp1 + c.transpose() * ykp1;
    ConstraintViolationVecY << (ykp1 - x_bar_table), (x_ubar_table - ykp1);
    // std::cout << "Y -> k+1: " << k+1 << " ObjValY = " << fkp1Y << " ConstraintViolationY = " << ConstraintViolationVecY.maxCoeff() << std::endl;
    if (ConstraintViolationVecY.maxCoeff() <= ConstraintTol && std::abs(fkp1Y - ObjValTrue) <= ObjValTol) {
      // std::cout << "Y -> k+1: " << k+1 << " ObjValY = " << fkp1Y << " ConstraintViolationY = " << ConstraintViolationVecY.maxCoeff() << std::endl;
      n_iter = (k + 1);  // the number of iterations
      exitflag = 1;
      xk = ykp1;
      fval = fkp1Y;
      ConstraintViolation = ConstraintViolationVecY.maxCoeff();
      break;
    }
#endif
    // move to next step
    yk = ykp1;
    zk = zkp1;
  }
  // if iteration terminates by MaxIter setting
  if (exitflag != 1) {
    fval = ((-umat).transpose().colwise() + xk).cwiseAbs().colwise().sum() * w;
    ConstraintViolationVecX = (b - A * xk).cwiseAbs();
    ConstraintViolation = ConstraintViolationVecX.maxCoeff();
  }
  // solution
  plhs[0] = mxCreateDoubleMatrix(n_x, 1, mxREAL);
  Eigen::Map<Eigen::VectorXd> Out0(mxGetPr(plhs[0]), n_x);
  Out0 << xk;
  // objective function value
  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
  Eigen::Map<Eigen::VectorXd> Out1(mxGetPr(plhs[1]), 1);
  Out1(0) = fval;
  // exitflag
  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
  Eigen::Map<Eigen::VectorXd> Out2(mxGetPr(plhs[2]), 1);
  Out2(0) = exitflag;
  // options
  plhs[3] = mxCreateDoubleMatrix(3, 1, mxREAL);
  Eigen::Map<Eigen::VectorXd> Out3(mxGetPr(plhs[3]), 3);
  t_end = std::chrono::system_clock::now();  // stop timer
  if (exitflag == 1) {
    t_iter = double(std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count()) * 1e-9;  // convert micro-seconds to seconds
  }
  Out3(0) = t_iter;
  Out3(1) = n_iter;
  Out3(2) = ConstraintViolation;
}
