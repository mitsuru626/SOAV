% Numerical Experiment of SOAV Minimization
% Mitsuru Toyoda (Tokyo Metropolitan University)
%

%% 
clear;
clc;
ConsiderQuadraticCost = true; % if true, quadratic cost is included
CompileMex = false; % if true, compile MEX files
n_trial = 10; % number of trials of ADMM and LP calculation
n_trial_conventional = 0; %   number of trials of conventional ADMM
gamma = 1; % parameter $gamma$ in ADMM

%% Complie Mex functions
if CompileMex
    clear soav_bisec_fixed.mexw64;
    delete soav_bisec_fixed.mexw64;
    clear soav_bisec_adaptive.mexw64;
    delete soav_bisec_adaptive.mexw64;
    clear soav_interp.mexw64;
    delete soav_interp.mexw64;

    % proposed method (bisection search, fixed gamma)
    mex -I"C:\Library\cpp\eigen-3.3.8" soav_bisec.cpp
    movefile soav_bisec.mexw64 soav_bisec_fixed.mexw64 f
    % proposed method (bisection search, adaptive gamma)
    mex -I"C:\Library\cpp\eigen-3.3.8" -DADAPTIVE_ON soav_bisec.cpp
    movefile soav_bisec.mexw64 soav_bisec_adaptive.mexw64 f
    % proposed method (table)
    mex -I"C:\Library\cpp\eigen-3.3.8" soav_interp.cpp
    if not(ConsiderQuadraticCost)
        clear soav_conventional.mexw64;
        delete soav_conventional.mexw64;
        % conventional method
        mex -I"C:\Library\cpp\eigen-3.3.8" soav_conventional.cpp
    end
end

%% make problems data as MAT file
nwTable = [1:20:201];
% nwTable = [1:20:201];
NTable = 10 * ones(size(nwTable));sub_instance_make(nwTable,NTable,ConsiderQuadraticCost);
SuffixList{1} = ['N_', num2str(min(NTable)), '_', num2str(max(NTable)), '_nw_', num2str(min(nwTable)), '_', num2str(max(nwTable))];
NTable = 20 * ones(size(nwTable));sub_instance_make(nwTable,NTable,ConsiderQuadraticCost);
SuffixList{2} = ['N_', num2str(min(NTable)), '_', num2str(max(NTable)), '_nw_', num2str(min(nwTable)), '_', num2str(max(nwTable))];
NTable = 50 * ones(size(nwTable));sub_instance_make(nwTable,NTable,ConsiderQuadraticCost);
SuffixList{3} = ['N_', num2str(min(NTable)), '_', num2str(max(NTable)), '_nw_', num2str(min(nwTable)), '_', num2str(max(nwTable))];
NTable = 100 * ones(size(nwTable));sub_instance_make(nwTable,NTable,ConsiderQuadraticCost);
SuffixList{4} = ['N_', num2str(min(NTable)), '_', num2str(max(NTable)), '_nw_', num2str(min(nwTable)), '_', num2str(max(nwTable))];

%% start simulation
SaveFileNameList = cell(4, 1);
LoadFileName = ['Probrem_', SuffixList{1}, '.mat'];
SaveFileName = ['Result_', SuffixList{1}, '.mat'];
SaveFileNameList{1} = SaveFileName;
run('sub_simulate_LP_ADMM.m');

LoadFileName = ['Probrem_', SuffixList{2}, '.mat'];
SaveFileName = ['Result_', SuffixList{2}, '.mat'];
SaveFileNameList{2} = SaveFileName;
run('sub_simulate_LP_ADMM.m');

LoadFileName = ['Probrem_', SuffixList{3}, '.mat'];
SaveFileName = ['Result_', SuffixList{3}, '.mat'];
SaveFileNameList{3} = SaveFileName;
run('sub_simulate_LP_ADMM.m');

LoadFileName = ['Probrem_', SuffixList{4}, '.mat'];
SaveFileName = ['Result_', SuffixList{4}, '.mat'];
SaveFileNameList{4} = SaveFileName;
run('sub_simulate_LP_ADMM.m');

%% load data
clearvars -except SaveFileNameList n_trial_conventional ConsiderQuadraticCost;
load(SaveFileNameList{1});
load(SaveFileNameList{2});
load(SaveFileNameList{3});
load(SaveFileNameList{4});

%% plot results
figure();

subplot(411);
hold('on');
box('on');
plot(nwTable, solver_time_LP_interior_N_10_10_nw_1_201, 'ko-', 'MarkerSize', 3);
plot(nwTable, solver_time_LP_psimplex_N_10_10_nw_1_201, 'k^-', 'MarkerSize', 3);
plot(nwTable, solver_time_LP_dsimplex_N_10_10_nw_1_201, 'kv-', 'MarkerSize', 3);
plot(nwTable, solver_time_bisec_fixed_N_10_10_nw_1_201, 'bo-', 'MarkerSize', 3);
plot(nwTable, solver_time_bisec_adaptive_N_10_10_nw_1_201, 'go-', 'MarkerSize', 3);
plot(nwTable, solver_time_interp_N_10_10_nw_1_201, 'ro-', 'MarkerSize', 3);
xlim([min(nwTable), max(nwTable)]);
ax = gca();
ax.XAxis.TickValues = [];
ylabel(['Time [s]', char(10), '{\itN}=10']);
ylim([0, 0.1]); %ylim([0,0.05]);

subplot(412);
hold('on');
box('on');
plot(nwTable, solver_time_LP_interior_N_20_20_nw_1_201, 'ko-', 'MarkerSize', 3);
plot(nwTable, solver_time_LP_psimplex_N_20_20_nw_1_201, 'k^-', 'MarkerSize', 3);
plot(nwTable, solver_time_LP_dsimplex_N_20_20_nw_1_201, 'kv-', 'MarkerSize', 3);
plot(nwTable, solver_time_bisec_fixed_N_20_20_nw_1_201, 'bo-', 'MarkerSize', 3);
plot(nwTable, solver_time_bisec_adaptive_N_20_20_nw_1_201, 'go-', 'MarkerSize', 3);
plot(nwTable, solver_time_interp_N_20_20_nw_1_201, 'ro-', 'MarkerSize', 3);
xlim([min(nwTable), max(nwTable)]);
ax = gca();
ax.XAxis.TickValues = [];
ylabel(['Time [s]', char(10), '{\itN}=20']);
ylim([0, 0.1]);

subplot(413);
hold('on');
box('on');
plot(nwTable, solver_time_LP_interior_N_50_50_nw_1_201, 'ko-', 'MarkerSize', 3);
plot(nwTable, solver_time_LP_psimplex_N_50_50_nw_1_201, 'k^-', 'MarkerSize', 3);
plot(nwTable, solver_time_LP_dsimplex_N_50_50_nw_1_201, 'kv-', 'MarkerSize', 3);
plot(nwTable, solver_time_bisec_fixed_N_50_50_nw_1_201, 'bo-', 'MarkerSize', 3);
plot(nwTable, solver_time_bisec_adaptive_N_50_50_nw_1_201, 'go-', 'MarkerSize', 3);
plot(nwTable, solver_time_interp_N_50_50_nw_1_201, 'ro-', 'MarkerSize', 3);
xlim([min(nwTable), max(nwTable)]);
ax = gca();
ax.XAxis.TickValues = [];
ylabel(['Time [s]', char(10), '{\itN}=50']);
ylim([0, 0.5]);

subplot(414);
hold('on');
box('on');
plot(nwTable, solver_time_LP_interior_N_100_100_nw_1_201, 'ko-', 'MarkerSize', 3);
plot(nwTable, solver_time_LP_psimplex_N_100_100_nw_1_201, 'k^-', 'MarkerSize', 3);
plot(nwTable, solver_time_LP_dsimplex_N_100_100_nw_1_201, 'kv-', 'MarkerSize', 3);
plot(nwTable, solver_time_bisec_fixed_N_100_100_nw_1_201, 'bo-', 'MarkerSize', 3);
plot(nwTable, solver_time_bisec_adaptive_N_100_100_nw_1_201, 'go-', 'MarkerSize', 3);
plot(nwTable, solver_time_interp_N_100_100_nw_1_201, 'ro-', 'MarkerSize', 3);
xlim([min(nwTable), max(nwTable)]);
ax = gca();
ax.XAxis.TickValues = nwTable(1:end);
ylabel(['Time [s]', char(10), '{\itN}=100']);
ylim([0, 0.5]);
xlabel('\itn_{\itw}');
if ConsiderQuadraticCost
    leg = legend({'QP (interior point)', 'QP (primal simplex)', 'QP (dual simplex)', ...
        'Proposed (fixed \gamma)', 'Proposed (adaptive \gamma)', 'Proposed (table)'}, 'Orientation', 'vertical');
else
    leg = legend({'LP (interior point)', 'LP (primal simplex)', 'LP (dual simplex)', ...
        'Proposed (fixed \gamma)', 'Proposed (adaptive \gamma)', 'Proposed (table)'}, 'Orientation', 'vertical');
end
leg.NumColumns = 2;

%% plot results of conventional ADMM
if n_trial_conventional ~= 0
    figure();
    subplot(212);
    hold('on');
    box('on');
    % 1
    data = solver_time_conventional_N_10_10_nw_1_201;
    ind_terminated1 = find(data == -1);
    data(ind_terminated1) = 1800;
    plot(nwTable, data, 'ko-', 'MarkerSize', 3);
    % 2
    data = solver_time_conventional_N_20_20_nw_1_201;
    ind_terminated2 = find(data == -1);
    data(ind_terminated2) = 1800;
    plot(nwTable, data, 'bo-', 'MarkerSize', 3);
    % 3
    data = solver_time_conventional_N_50_50_nw_1_201;
    ind_terminated3 = find(data == -1);
    data(ind_terminated3) = 1800;
    plot(nwTable, data, 'co-', 'MarkerSize', 3);
    % 4
    data = solver_time_conventional_N_100_100_nw_1_201;
    ind_terminated4 = find(data == -1);
    data(ind_terminated4) = 1800;
    plot(nwTable, data, 'go-', 'MarkerSize', 3);
    ylim([0, 1800]);
    ylabel(['Time [s]']);
    xlim([min(nwTable), max(nwTable)]);
    xlabel('\itn_{\itw}');
    ax = gca();
    ax.XAxis.TickValues = nwTable(1:end);
    ax.YAxis.TickValues = [0, 1000, 1800];
    leg = legend({'{\itN}=10', '{\itN}=20', '{\itN}=50', '{\itN}=100'}, 'Orientation', 'horizontal');
    leg.NumColumns = 4;
end