% matlab settings

%% clear workspace, close open figures
clearvars;
close all;
clc;

%% addpath...
fprintf('addpath casadi\n');
addpath('C:\Program Files\MATLAB\R2018b\toolbox\casadi-windows-matlabR2016a-v3.5.3');

fprintf('addpath mpctools\n');
addpath('C:\Program Files\MATLAB\R2018b\toolbox\mpctools');

% import mpctools
fprintf('import mpctools\n');
mpc = import_mpctools();