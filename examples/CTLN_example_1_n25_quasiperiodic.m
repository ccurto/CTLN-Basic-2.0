% CTLN_example_1_n25_quasiperiodic.m
%
% comments: n25 graph with quasiperiodic attractor.
% This corresponds to Figure 2C of the paper https://arxiv.org/abs/1605.04463
%
% Option at end of file to save data as: examples/CTLN_example_1_n25_quasiperiodic.mat

% adjacency matrix 
sA = [ 0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1; ...
       1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1; ...
       1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1; ...
       1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1; ...
       0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1; ...
       1  1  1  1  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1; ...
       1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1; ...
       1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1; ...
       1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1; ...
       1  1  1  1  1  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1; ...
       1  1  1  1  1  1  1  1  1  1  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0; ...
       1  1  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0; ...
       1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0; ...
       1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0; ...
       1  1  1  1  1  1  1  1  1  1  0  1  1  1  0  0  0  0  0  0  0  0  0  0  0; ...
       0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  0  0  0  0  1  0  0  0  0  0; ...
       0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0  0; ...
       0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0  0; ...
       0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0  0  0  0; ...
       0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  0  1  1  1  0  0  0  0  0  0; ...
       0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  0  0  0  0  1; ...
       0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  0  0  0  0; ...
       0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  0  0  0; ...
       0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0; ...
       0  0  0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  0  1  1  1  0 ]; 

% simulation parameters 
n = 25;    % number of neurons 
e = 0.25;  % epsilon value 
d = 0.5;   % delta value 
theta = 1;   % external drive
T = 600;   % time in units of membrane time constant

% initial conditions 
X0cell=[]; % initialize cell array of interesting initial conditions 
X0cell{1} = [ 0.071   0.003   0.028   0.005   0.010   0.082   0.069   0.032   0.095   ...
              0.003   0.044   0.038   0.077   0.080   0.019   0.049   0.045   0.065   ...
              0.071   0.075   0.028   0.068   0.066   0.016   0.012];

X0 = X0cell{1};   % for other initial conditions, modify command to access other entries in cell array 

% projection - load some nice ones, and pick one
load('example_1_projections.mat','proj1','proj2','proj3');
proj = proj3;

% proj = rand(n,2); % to generate a random one

% colors 
colors = [ ]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% To save .mat file, uncomment commands below 

% ex_num = 1;
% ex_name = 'n25_quasiperiodic';
% graph_comments = 'n25 graph with quasiperiodic attractor from Figure 2 of CTLN preprint';
% 
% save('examples/CTLN_example_1_n25_quasiperiodic.mat', 'sA', 'n', 'e', 'd','theta','T','X0','X0cell','proj','colors','ex_num','ex_name','graph_comments')
