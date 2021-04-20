function [x_opt, fval, output] = soav_LP_QP(u, w, lb, ub, Aeq, beq, Q, c, options)
% soav_LP_QP - SOAV minimization by LP/QP solver provided by Gurobi
%
% [x_opt,fval,output] = soav_LP_QP(u,w,lb,ub,Aeq,beq,Q,c,options)
%
% Arguments:
%   u,w,lb,ub,Aeq,beq,Q,c: data of considered SOAV minimization problem
%       note u = [u_1'; ... u_{n_w}'] \in \R^{n_w \times n_x}
%   options: options.Method indicates algorithm
%
% Mitsuru Toyoda (Tokyo Metropolitan University)
% 

n_w = length(w);
n_x = size(Aeq, 2);
n_b = length(beq);
Q_LP = sparse((2*n_w+1)*n_x, (2 * n_w + 1)*n_x);
Q_LP(1:n_x, 1:n_x) = Q;
f_LP = [c; kron(sparse(w), spones(ones(2 * n_x, 1)))];
Aeq_LP = ...
    [Aeq, sparse(n_b, 2 * n_w * n_x); repmat(-speye(n_x, n_x), [n_w, 1]), kron(speye(n_w, n_w), [speye(n_x, n_x), -speye(n_x, n_x)])];
beq_LP = [beq; -reshape(u', [n_w * n_x, 1])];
A_LP = blkdiag([-speye(n_x, n_x); speye(n_x, n_x)], -speye(2 * n_w * n_x, 2 * n_w * n_x));
b_LP = [-lb; ub; sparse(2 * n_w * n_x, 1)];

% Gurobi model
model = struct();
if norm(Q) ~= 0 % QP problem
    model.Q = Q_LP / 2; % gurobi requires Q as a form x'*Q*x
end
model.obj = full(f_LP); % f must be dense
model.A = [sparse(A_LP); sparse(Aeq_LP)]; % A must be sparse
model.sense = [repmat('<', size(A_LP, 1), 1); repmat('=', size(Aeq_LP, 1), 1)];
model.rhs = full([b_LP; beq_LP]); % rhs must be dense
model.lb = -inf(size(model.A, 2), 1);
model.ub = inf(size(model.A, 2), 1);

% Extract relevant Gurobi parameters from (subset of) options
params = struct();
params.OutputFlag = 0; % display of iteration is inactivated
params.Crossover = 0;
params.FeasibilityTol = 1e-8;
params.OptimalityTol = 1e-8;
params.Method = options.Method; % optimization method
% -1=automatic,
%  0=primal simplex,
%  1=dual simplex,
%  2=barrier,

% Solve model with Gurobi
tic;
result = gurobi(model, params);
MfileTime = toc; % Measured by Mfile
GurobiTime = result.runtime; % Measured by Gurobi MEX interface

output = struct();
output.MfileTime = MfileTime;
output.GurobiTime = GurobiTime;
output.IterCount = result.itercount;
output.BarIterCount = result.baritercount;
% first part (original x)
x_opt = result.x(1:n_x);
fval = w' * sum(abs(repmat(x_opt', [n_w, 1]) - u), 2) + 1 / 2 * x_opt' * Q * x_opt + c' * x_opt;
ConstraintViolation = max([0; abs(beq - Aeq * x_opt); x_opt - ub; lb - x_opt]);
% calculate x_i = v_i^+ - v_i^ + u_i
fvalList = zeros(n_w+1, 1);
fvalList(1) = fval;
ConstraintViolationList = zeros(n_w+1, 1);
ConstraintViolationList(1) = ConstraintViolation;
for i = 1:n_w
    LastInd = n_x + 2 * n_x * (i - 1);
    x_i = result.x(LastInd+(1:n_x)) - result.x(LastInd+((n_x + 1):2 * n_x)) + u(i, :)';
    fvalList(i+1, 1) = w' * sum(abs(repmat(x_i', [n_w, 1]) - u), 2) + 1 / 2 * x_i' * Q * x_i + c' * x_i;
    ConstraintViolationList(i+1, 1) = max([0; abs(beq - Aeq * x_opt); x_opt - ub; lb - x_opt]);
end
output.fvalList = fvalList;
output.ConstraintViolationList = ConstraintViolationList;
end