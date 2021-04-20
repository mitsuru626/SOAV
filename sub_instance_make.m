function [] = sub_instance_make(nwTable, NTable, ConsiderQuadraticCost)
% sub_instance_make - Save problem instance in MAT-file
%
% [] = sub_instance_make(nwTable,NTable)
%
% Arguments:
%   nwTable: [n_w(1); ,,,; n_w(n_instance)]
%   NTable:  [N(1)  ; ,,,; N(n_instance)  ]
%   ConsiderQuadraticCost: if true, quadratic cost is included
%
% Mitsuru Toyoda (Tokyo Metropolitan University)
% 

% parameters
Kt = 0.0185;
Ke = 0.0285;
R = 1.15;
J = 0.0021;
ks = 0.087;
B = 0.0088;
theta0 = deg2rad(12);

A_c = [0, 1; -ks / J, -B / J - Kt * Ke / (R * J)];
B_c = [0; Kt / (R * J)];
d_c = [0; ks * theta0 / J];
x0 = [deg2rad(20); 0];
x_target = [deg2rad(80); 0];
lb = -6;
ub = 6;
n_x = size(B_c, 1);
n_u = size(B_c, 2);
T = 0.5;

AeqTable = cell(length(nwTable), 1);
beqTable = cell(length(nwTable), 1);
lbTable = cell(length(nwTable), 1);
ubTable = cell(length(nwTable), 1);
uTable = cell(length(nwTable), 1);
wTable = cell(length(nwTable), 1);
QTable = cell(length(nwTable), 1);
cTable = cell(length(nwTable), 1);

for i_point = 1:length(nwTable)
    n_w = nwTable(i_point);
    N = NTable(i_point);
    fprintf('n_w = %g\n', n_w);
    h = T / N;

    if ConsiderQuadraticCost
        QTable{i_point, 1} = (h / R) * eye(N, N);
        cTable{i_point, 1} = zeros(N, 1);
    else
        QTable{i_point, 1} = zeros(N, N);
        cTable{i_point, 1} = zeros(N, 1);
    end

    if n_w == 1
        w = h;
        u_div = (lb + ub) / 2;
    else
        w = ones(n_w, 1) * h;
        u_div = linspace(lb, ub, n_w)';
    end

    A = expm(A_c*h);
    IeAh = A_c \ (expm(A_c * h) - eye(n_x, n_x));
    B = IeAh * B_c;
    d = IeAh * d_c;
    Delta = zeros(n_x, 1);
    for j = 1:N
        Phi(:, (j - 1)*n_u+1:j*n_u) = A^(N - j) * B;
        Delta = Delta + A^(N - j) * d;
    end
    beq = x_target - A^(N) * x0 - Delta;
    Aeq = Phi;
    AeqTable{i_point, 1} = Aeq;
    beqTable{i_point, 1} = beq;
    lbTable{i_point, 1} = repmat(lb, [N, 1]);
    ubTable{i_point, 1} = repmat(ub, [N, 1]);
    uTable{i_point, 1} = repmat(u_div, [1, N]);
    wTable{i_point, 1} = w;
end

% save as mat-file
SaveFilenameSuffix = ['N_', num2str(min(NTable)), '_', num2str(max(NTable)), '_nw_', num2str(min(nwTable)), '_', num2str(max(nwTable))];
save(['Probrem_', SaveFilenameSuffix], ...
    'nwTable', 'NTable', 'AeqTable', 'beqTable', 'lbTable', 'ubTable', 'uTable', 'wTable', 'QTable', 'cTable');
end