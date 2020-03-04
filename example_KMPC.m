% example_KMPC.m
%   This script identifies a linear Koopman model from data, constructs an
%   MPC controller, and simulates the performance of the controller on a
%   trajectory following task.
%
%   The file "/datafiles/softrobot_train-13_val-4.mat" contains the data
%   used to train the soft robot model in the paper entitled "Data-driven
%   control of soft robots using Koopman operatory theory" by Daniel
%   Bruder, Xun Fu, R. Brent Gillespie, C. David Remy, and Ram Vasudevan.
%
%   The files found in "/ref-trajectories" are the trajectories used for
%   the experiments in the paper as well.
%
%   Daniel Bruder
%   bruderd@umich.edu
%   2019-11-23

%%
%-------------------------------------------------------------------------%
%----------------- Identify Koopman Model from Data ----------------------%
%-------------------------------------------------------------------------%

% load in data from file
[ datafile_name , datafile_path ] = uigetfile( 'datafiles/*.mat' , 'Choose data file for sysid...' );
data4sysid = load( [datafile_path , datafile_name] );

% construct sysid class
ksysid = ksysid( data4sysid, ...
        'model_type' , 'linear' ,...    % model type (linear or nonlinear)
        'obs_type' , { 'poly' } ,...    % type of basis functions
        'obs_degree' , [ 2 ] ,...       % "degree" of basis functions
        'snapshots' , Inf ,...          % Number of snapshot pairs
        'lasso' , [ Inf ] ,...           % L1 regularization term
        'delays' , 1 );                 % Numer of state/input delays
    
% train linear Koopman model(s)
ksysid = ksysid.train_models;

% validate model(s)
results = cell( size(ksysid.candidates) );    % store results in a cell array
err = cell( size(ksysid.candidates) );    % store error in a cell array 

if iscell(ksysid.candidates)
    for i = 1 : length(ksysid.candidates)
        [ results{i} , err{i} ] = ksysid.valNplot_model( i );
    end
else
    [ results{1} , err{1} ] = ksysid.valNplot_model;
end


%%
%-------------------------------------------------------------------------%
%----------------- Construct Koopman MPC Controller ----------------------%
%-------------------------------------------------------------------------%

kmpc = kmpc( ksysid ,...
        'horizon' , 25 ,...
        'input_bounds' , [ 0 , 10 ],... 
        'input_slopeConst' , [0.5e-2],... 
        'input_smoothConst' , [1e-1],... 
        'state_bounds' , [] ,...
        'cost_running' , 0.1 ,...
        'cost_terminal' , 100 ,...
        'cost_input' , 0 ); 
    
%%
%-------------------------------------------------------------------------%
%-------------------------- Run Simulation -------------------------------%
%-------------------------------------------------------------------------%

% load in reference trajectory from file
[ reffile_name , reffile_path ] = uigetfile( 'ref-trajectories/*.mat' , 'Choose reference trajectory file...' );
temp = load( [reffile_path , reffile_name] );
ref = temp.ref;

% run simulation
y0 = [ 1 ,1 ];      % initial laser dot position
u0 = [ 0 , 0 , 0];  % initial input
sim = kmpc.run_simulation( ref , y0 , u0 );


%%
%-------------------------------------------------------------------------%
%---------------------------- Plot Results -------------------------------%
%-------------------------------------------------------------------------%

% Define color(s)
cb_blue = [55,126,184] ./ 255;  % blue

figure;
hold on;
plot( sim.R(:,1) , sim.R(:,2) , 'color' , [0.25 0.25 0.25 1] , 'LineWidth' , 2);
plot( sim.Y(:,1) , sim.Y(:,2) , 'color' , cb_blue , 'LineWidth' , 2);
xlabel('$x_1$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 20);
ylabel('$x_2$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 20);
hold off;
grid on; box on;





