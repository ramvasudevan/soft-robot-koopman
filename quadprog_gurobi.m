function [ x , results ] = quadprog_gurobi( H, f, A, b )
%quadprog: Solves a quadratic program using Gurobi
%   Detailed explanation goes here

%% solve QP using gurobi

% %YOU MUST REMOVE THIS !!!!!!!!!!!!!! DONT FORGET!!!!!!!!!!!!!!!!!!!!!!!!!
% H = H + eye( size(H,1) ) * 1e-6;

model.Q = sparse(H);  % removed 0.5 on 2020-02-25 (was added on 2019-10-10)
model.A = sparse(A);
model.obj = f';
model.rhs = b';
model.sense = '<';
model.lb = -Inf * ones( size(H,1) , 1 );    % change from default lower bound which is zero

% solve using gurobi
params.outputflag = 1;  % if 0 will turn off the display while running
results = gurobi(model , params);

% extract the solution from the results struct
if strcmp( results.status , 'OPTIMAL' )
    x = results.x;
else
    disp( [ 'NOT SOLVED BECAUSE MODEL IS ' , results.status ] );
    x = nan( size(A,2) , 1 );
end

end

