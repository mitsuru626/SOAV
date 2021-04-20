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
#include <unsupported/Eigen/KroneckerProduct>

#include "mex.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  // start timer for computation time
  std::chrono::system_clock::time_point t_start, t_end;
  t_start = std::chrono::system_clock::now();
  // check the number of arguments
  if (nrhs != 10) {
    mexErrMsgIdAndTxt("MATLAB:soav_conventional:nrhs", "the number of arguments are 10. ");
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
  double gamma = mxGetScalar(prhs[6]);
  Eigen::Map<Eigen::VectorXd> y0(mxGetPr(prhs[7]), n_x * (n_w + 1) + n_b);
  Eigen::Map<Eigen::VectorXd> z0(mxGetPr(prhs[8]), n_x * (n_w + 1) + n_b);
#ifdef PDOptimalityCriteria_ON  // [eps_rel; eps_abs]
  Eigen::Map<Eigen::VectorXd> StoppingCriteriaPara(mxGetPr(prhs[9]), 2);
#else  // [ObjValTrue; ObjValTol; ConstraintTol]
  Eigen::Map<Eigen::VectorXd> StoppingCriteriaPara(mxGetPr(prhs[9]), 3);
#endif

  Eigen::MatrixXd L(n_x * (n_w + 1) + n_b, n_x);
  L << Eigen::kroneckerProduct(Eigen::MatrixXd::Ones(n_w + 1, 1), Eigen::MatrixXd::Identity(n_x, n_x)), A;
  // std::cout << "L = " <<  L << std::endl;
  // std::cout << "#(rows of L) = " <<  L.rows() << " #(cols of L) = " << L.cols() << std::endl;
  Eigen::MatrixXd Lplus = ((n_w + 1) * Eigen::MatrixXd::Identity(n_x, n_x) + A.transpose() * A).llt().solve(L.transpose());
  // std::cout << "#(rows of Lplus) = " <<  Lplus.rows() << " #(cols of Lplus) = " << Lplus.cols() << std::endl;
  Eigen::VectorXd xk(n_x);
  Eigen::VectorXd yk(n_x * (n_w + 1) + n_b);
  Eigen::VectorXd zk(n_x * (n_w + 1) + n_b);
  Eigen::VectorXd ykp1(n_x * (n_w + 1) + n_b);
  Eigen::VectorXd zkp1(n_x * (n_w + 1) + n_b);
  Eigen::VectorXd Lxpz(n_x * (n_w + 1) + n_b);
  Eigen::VectorXd Lxpz_seg(n_x);
  yk = y0;
  zk = z0;

  // variables for output of function
  int MaxIter = int(1e8);
  int n_iter = MaxIter;
  double t_iter = INFINITY;
  double TimeLimit = 30 * 60;
  double myInf = 1e20;
  double fval = myInf, ConstraintViolation = myInf;
  Eigen::VectorXd ConstraintViolationVecX(2 * n_x + n_b);

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
  Eigen::VectorXd yjkp1(n_x);
  Eigen::VectorXd ConstraintViolationVecY(2 * n_x + n_b);
// std::cout << "ObjVal = " <<  ObjValTrue << " ConstraintTol = " << ConstraintTol << std::endl;
#endif

  for (int k = 0; k < MaxIter; k++) {
    // TimeLimit
    if (k % 10000 == 0) {
      t_end = std::chrono::system_clock::now();
      t_iter = double(std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count()) * 1e-9;
      if (t_iter > TimeLimit) {
        exitflag = -1;
        break;
      }
    }
    // x[k]: proximal operation for $f$
    xk = Lplus * (yk - zk);

#ifndef PDOptimalityCriteria_ON
    // objective function value and constraint violation
    fkX = ((-umat).transpose().colwise() + xk).cwiseAbs().colwise().sum() * w;
    ConstraintViolationVecX.segment(0, n_x) = xk - x_bar_table;
    ConstraintViolationVecX.segment(n_x, n_x) = x_ubar_table - xk;
    ConstraintViolationVecX.segment(2 * n_x, n_b) = (b - A * xk).cwiseAbs();
    // std::cout << "X -> k: " << k << " ObjValX = " << fkX << " ConstraintViolationX = " << ConstraintViolationVecX.maxCoeff() << std::endl;
    // stopping criteria
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
    Lxpz = L * xk + zk;
    for (int j = 0; j < n_w; j++) {
      Lxpz_seg = Lxpz.segment(j * n_x, n_x) - umat.row(j).transpose();
      ykp1.segment(j * n_x, n_x) = umat.row(j).transpose().array() + Lxpz_seg.array().sign() * (Lxpz_seg.array().abs() - gamma * w(j)).cwiseMax(0);
    }
    Lxpz_seg = Lxpz.segment(n_w * n_x, n_x);
    ykp1.segment(n_w * n_x, n_x) = (Lxpz_seg.cwiseMax(x_ubar_table)).cwiseMin(x_bar_table);
    ykp1.segment((n_w + 1) * n_x, n_b) = b;
    zkp1 = Lxpz - ykp1;  // z(k+1)

#ifdef PDOptimalityCriteria_ON
    rkp1_norm = (L * xk - ykp1).norm();                               // x(k+1) - y(k+1)
    skp1_norm = (-1 / gamma * L.transpose() * (ykp1 - yk)).norm();  // y(k+1) - y(k)
    eps_pri = sqrt(n_x * (n_w + 1) + n_b) * eps_abs + eps_rel * std::max((L * xk).norm(), ykp1.norm());
    eps_dual = sqrt(n_x) * eps_abs + eps_rel * (zkp1 / gamma).norm();
    // std::cout << "k+1: " << k+1 << " eps_pri = " <<  eps_pri << " eps_dual = " <<  eps_dual << std::endl;
    if (rkp1_norm <= eps_pri && skp1_norm <= eps_dual) {
      n_iter = k + 1;
      exitflag = 1;
      fval = ((-umat).transpose().colwise() + xk).cwiseAbs().colwise().sum() * w;
      ConstraintViolationVecX.segment(0, n_x) = xk - x_bar_table;
      ConstraintViolationVecX.segment(n_x, n_x) = x_ubar_table - xk;
      ConstraintViolationVecX.segment(2 * n_x, n_b) = (b - A * xk).cwiseAbs();
      ConstraintViolation = ConstraintViolationVecX.maxCoeff();
      break;
    }
#else
    for (int j = 0; j < n_w + 1; j++) {
      yjkp1 = ykp1.segment(j * n_x, n_x);
      fkp1Y = ((-umat).transpose().colwise() + yjkp1).cwiseAbs().colwise().sum() * w;
      ConstraintViolationVecY.segment(0, n_x) = yjkp1 - x_bar_table;
      ConstraintViolationVecY.segment(n_x, n_x) = x_ubar_table - yjkp1;
      ConstraintViolationVecY.segment(2 * n_x, n_b) = (b - A * yjkp1).cwiseAbs();
      // std::cout << "Y" << j << "-> k+1: " << k+1 << " ObjValY = " << fkp1Y << " ConstraintViolationY = " << ConstraintViolationVecY.maxCoeff() << std::endl;
      if (ConstraintViolationVecY.maxCoeff() <= ConstraintTol && std::abs(fkp1Y - ObjValTrue) <= ObjValTol) {
        // std::cout << "Y" << j << "-> k+1: " << k+1 << " ObjValY = " << fkp1Y << " ConstraintViolationY = " << ConstraintViolationVecY.maxCoeff() << std::endl;
        n_iter = (k + 1);  // the number of iterations
        exitflag = 1;
        xk = yjkp1;
        fval = fkp1Y;
        ConstraintViolation = ConstraintViolationVecY.maxCoeff();
        break;
      }
    }
    if (exitflag == 1) {
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
    ConstraintViolationVecX.segment(0, n_x) = xk - x_bar_table;
    ConstraintViolationVecX.segment(n_x, n_x) = x_ubar_table - xk;
    ConstraintViolationVecX.segment(2 * n_x, n_b) = (b - A * xk).cwiseAbs();
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
