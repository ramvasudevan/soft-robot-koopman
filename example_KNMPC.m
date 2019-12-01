% example_KNMPC.m
%   This script identifies a nonlinear Koopman model from data, constructs
%   an open-loop NMPC controller, and simulates the performance of the 
%   controller on a trajectory following task.
%
%   The file "/datafiles/softrobot_train-13_val-4.mat" contains the data
%   used to train the soft robot model in the paper entitled "Data-driven
%   control of soft robots using Koopman operatory theory" by Daniel
%   Bruder, Xun Fu, R. Brent Gillespie, C. David Remy, and Ram Vasudevan.
%
%   The files found in "/ref-trajectories" are the trajectories used for
%   the experiments in the paper as well.
%
%   NOTE: Constructing the NMPC controller requires GPOPS-II to be
%   installed. See http://www.gpops2.com/ for details.
%
%   Xun Fu (xunfu@umich.edu) 
%   Daniel Bruder (bruderd@umich.edu)
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
        'model_type' , 'nonlinear' ,... % model type (linear or nonlinear)
        'obs_type' , { 'poly' } ,...    % type of basis functions
        'obs_degree' , [ 3 ] ,...       % "degree" of basis functions
        'snapshots' , Inf ,...          % Number of snapshot pairs
        'lasso' , [ 10 ] ,...           % L1 regularization term
        'delays' , 0 );                 % Numer of state/input delays

    
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

% Save the nonlinear system vector field F as a function in local directory
matlabFunction( ksysid.model.F_sym , 'File' , 'vf_koopman' , 'Vars' , ...
    { ksysid.params.x , ksysid.params.u } );


%%
%-------------------------------------------------------------------------%
%----------- Construct Koopman NMPC controller (open-loop) ---------------%
%-------------------------------------------------------------------------%

% Choose reference trajectory for the controller to follow
[t_sample,t_ref, traj_ref, traj_ref_sc] = load_files( ksysid );

steps =(t_ref(end)-t_ref(1))/t_sample; 
dt = t_sample;

x0 = traj_ref_sc(1,:);
xf = traj_ref_sc(end,:);

xmin = -1;
xmax = 1;

slp_min = -0.1;
slp_max = 0.1;

%-------------------------------------------------------------------------%
%----------------- Provide All Bounds for Problem ------------------------%
%-------------------------------------------------------------------------%
t0 = t_ref(1);
tf = t_ref(end);

whole_time = t0:dt:tf;
phase_num = (tf-t0)/dt;

xMin = xmin*ones(1,2);
xMax = xmax*ones(1,2);

pMin = 0*ones(1,3*phase_num);
pMax = 1*ones(1,3*phase_num);
pGuess = zeros(1,3*phase_num);

auxdata.phase_num = phase_num;
auxdata.order = length(x0);
auxdata.xf = xf;
auxdata.dt = dt;

%-------------------------------------------------------------------------%
%------------------ Setup for Problem Bounds and Guess -------------------%
%-------------------------------------------------------------------------%
for i = 1:phase_num
    
    if i == 1
        bounds.phase(i).initialtime.lower = whole_time(i);
        bounds.phase(i).initialtime.upper = whole_time(i);
        bounds.phase(i).finaltime.lower = whole_time(i+1);
        bounds.phase(i).finaltime.upper = whole_time(i+1);
        bounds.phase(i).initialstate.lower = x0; 
        bounds.phase(i).initialstate.upper = x0;
        bounds.phase(i).state.lower = xMin;
        bounds.phase(i).state.upper = xMax;
        bounds.phase(i).finalstate.lower = traj_ref_sc(i+1,:); 
        bounds.phase(i).finalstate.upper = traj_ref_sc(i+1,:);
        bounds.phase(i).integral.lower = 0;
        bounds.phase(i).integral.upper = 10000;
        bounds.phase(i).path.lower = slp_min*ones(1,3);
        bounds.phase(i).path.upper = slp_max*ones(1,3);
        bounds.parameter.lower = pMin;
        bounds.parameter.upper = pMax;
        
        guess.phase(i).time     = [whole_time(i); whole_time(i+1)]; 
        guess.phase(i).state    = [x0; zeros(1,2)];
        guess.phase(i).integral = 0;
        guess.parameter = pGuess;
    
    elseif i == phase_num 
        bounds.phase(i).initialtime.lower = whole_time(i);
        bounds.phase(i).initialtime.upper = whole_time(i);
        bounds.phase(i).finaltime.lower = whole_time(i+1);
        bounds.phase(i).finaltime.upper = whole_time(i+1);
        bounds.phase(i).initialstate.lower = traj_ref_sc(i,:); 
        bounds.phase(i).initialstate.upper = traj_ref_sc(i,:);
        bounds.phase(i).state.lower = xMin;
        bounds.phase(i).state.upper = xMax;
        bounds.phase(i).finalstate.lower = xf; 
        bounds.phase(i).finalstate.upper = xf;
        bounds.phase(i).integral.lower = 0;
        bounds.phase(i).integral.upper = 10000;
        bounds.phase(i).path.lower = slp_min*ones(1,3);
        bounds.phase(i).path.upper = slp_max*ones(1,3);
        bounds.parameter.lower = pMin;
        bounds.parameter.upper = pMax;

        guess.phase(i).time     = [whole_time(i); whole_time(i+1)];
        guess.phase(i).state    = [zeros(1,2); xf];
        guess.phase(i).integral = 0;
        guess.parameter = pGuess;
        
    else
        bounds.phase(i).initialtime.lower = whole_time(i);
        bounds.phase(i).initialtime.upper = whole_time(i);
        bounds.phase(i).finaltime.lower = whole_time(i+1);
        bounds.phase(i).finaltime.upper = whole_time(i+1);
        bounds.phase(i).initialstate.lower = traj_ref_sc(i,:); 
        bounds.phase(i).initialstate.upper = traj_ref_sc(i,:);
        bounds.phase(i).state.lower = xMin;
        bounds.phase(i).state.upper = xMax;
        bounds.phase(i).finalstate.lower = traj_ref_sc(i+1,:); 
        bounds.phase(i).finalstate.upper = traj_ref_sc(i+1,:);
        bounds.phase(i).integral.lower = 0;
        bounds.phase(i).integral.upper = 10000;
        bounds.phase(i).path.lower = slp_min*ones(1,3);
        bounds.phase(i).path.upper = slp_max*ones(1,3);
        bounds.parameter.lower = pMin;
        bounds.parameter.upper = pMax;

        guess.phase(i).time     = [whole_time(i); whole_time(i+1)]; 
        guess.phase(i).state    = [zeros(1,2); zeros(1,2)];
        guess.phase(i).integral = 0;
        guess.parameter = pGuess;
        
    end
   

    
end 

%-------------------------------------------------------------------------%
%------------- Set up Event Constraints That Link Phases -----------------%
%-------------------------------------------------------------------------%
for i = 1:phase_num - 1 
    
    bounds.eventgroup(i).lower = zeros(1,3);
    bounds.eventgroup(i).upper = zeros(1,3);
   
end
    
%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method          = 'hp-LiuRao-Legendre';
mesh.tolerance       = 1e-3;
mesh.maxiterations   = 20;
mesh.colpointsmin    = 3;
mesh.colpointsmax    = 3;

for i = 1:phase_num
    mesh.phase(i).colpoints = 4*ones(1,10);
    mesh.phase(i).fraction  = 0.1*ones(1,10);
end

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name                           = 'GPOPS-EVENT-Problem';
setup.functions.continuous           = @gpops_Continuous;
setup.functions.endpoint             = @gpops_Endpoint;
setup.displaylevel                   = 0;
setup.bounds                         = bounds;
setup.guess                          = guess;
setup.mesh                           = mesh;
setup.nlp.solver                     = 'ipopt';
setup.nlp.snoptoptions.tolerance     = 1e-10;
setup.nlp.snoptoptions.maxiterations = 200;
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.nlp.ipoptoptions.tolerance     = 1e-5;
setup.derivatives.supplier           = 'sparseCD';
setup.derivatives.derivativelevel    = 'second';
setup.method                         = 'RPM-Differentiation';
setup.auxdata                        = auxdata;
%-------------------------------------------------------------------------%
%---------------------- Solve Problem Using GPOPS2 -----------------------%
%-------------------------------------------------------------------------%

tic;
output = gpops2(setup);
elapsed_time = toc;

%-------------------------------------------------------------------------%
%--------------------------------- Plots ---------------------------------%
%-------------------------------------------------------------------------%

% Define color(s)
cb_orange = [255,127,0] ./ 255;   % orange
cb_grey = [0.25 0.25 0.25 1]; % grey

solution = output.result.solution;
time = [];
state =[];
for i = 1:phase_num
    
    time = [time; solution.phase(i).time];
    state = [state; solution.phase(i).state];
    
end

control = solution.parameter;
control = reshape(control,[3,steps]);

% scale trajectories and inputs back up to real-world values
state = ksysid.scaleup.y( state );
traj_ref = ksysid.scaleup.y( traj_ref_sc );
x0 = ksysid.scaleup.y( x0 );
xf = ksysid.scaleup.y( xf );
control = ksysid.scaleup.u( control' )';

% plot control vs. reference trajectory
figure;
plot(state(:,1),state(:,2),'color' , cb_orange , 'LineWidth' , 2);
hold on;
plot(traj_ref(:,1),traj_ref(:,2),'color' , cb_grey , 'LineWidth' , 2);
hold on;
plot(x0(1),x0(2),'o','MarkerSize',8);
hold on;
plot(xf(1),xf(2),'*','MarkerSize',8);
xlabel('$x_1$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 20);
ylabel('$x_2$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 20);
legend('trajectory','reference_traj','starting point','ending point');
title('trajectory comparison');
grid on; box on;

% plot open-loop inputs
figure;
plot(control(1,:),'Linewidth',1.25);
hold on;
plot(control(2,:),'Linewidth',1.25);
hold on;
plot(control(3,:),'Linewidth',1.25);
xlabel('time');
ylabel('control input');
legend('control input 1','control input 2','control input 3');
title('control input');
