classdef kmpc
    %kmpc: Model predictive controller class
    %   Detailed explanation goes here
    
    properties
        params; % paramaters of the system
        model;  % linear model of the system
        lift;   % lifting functions for system
        basis;  % symbolic basis set of observables
        horizon;
        input_bounds;
        input_slopeConst;
        input_smoothConst;
        state_bounds;
        cost_running;
        cost_terminal;
        cost_input;
        projmtx; % projection matrix from liftes state (z) to reference state   
        cost;   % stores cost matrices
        constraints;    % stores constraint matrices
        set_constRHS;  % function that sets the value of the RHS of constraints
        get_zeta;   % function that constructs zeta from state and input data in time
        
        scaledown;  % functions for scaling to [-1,1]
        scaleup;    % functions for scaling from [-1,1]
    end
    
    methods
        % CLASS CONSTRUCTOR
        function obj = kmpc( sysid_class , varargin )
            %kmpc: Construct an instance of this class
            %   sysid_class - sysid class object with a model and params
            %    properties
            %   varargin - Name, Value pairs for class properties
            
            % take some properties/methods from the sysid class
            obj.params = sysid_class.params;
            obj.model = sysid_class.model;
            obj.lift = sysid_class.lift;
            obj.basis = sysid_class.basis;
            obj.get_zeta = @sysid_class.get_zeta;   % copies this method for convenience
            obj.scaledown = sysid_class.scaledown;
            obj.scaleup = sysid_class.scaleup;
            
            % define default values of properties
            obj.horizon = floor( 1 / obj.params.Ts );
            obj.input_bounds = [];  % [min,max] can be 1x2 or mx2
            obj.input_slopeConst = [];
            obj.input_smoothConst = [];
            obj.state_bounds = []; % [min,max] can be 1x2 or nx2 
            obj.cost_running = 0.1;
            obj.cost_terminal = 100;
            obj.cost_input = 0;
            obj.projmtx = obj.model.C;   % recovers measured state (could also use Cshape)
            obj.cost = [];
            obj.constraints = [];
            
            % replace default values with user input values
            obj = obj.parse_args( varargin{:} );
            
            % resize some properties if they aren't full size already
            obj = obj.expand_props;
            
            % get cost and constraint matrices
            obj = obj.get_costMatrices;
            obj = obj.get_constraintMatrices;   
        end
        
        % parse_args: Parses the Name, Value pairs in varargin
        function obj = parse_args( obj , varargin )
            %parse_args: Parses the Name, Value pairs in varargin of the
            % constructor, and assigns property values
            for idx = 1:2:length(varargin)
                obj.(varargin{idx}) = varargin{idx+1} ;
            end
        end
        
        % expand_props: Converts props from shorthand to fully defined
        function obj = expand_props( obj )
            %expand_props: Converts props from shorthand to fully defined
            %   e.g. input_bounds = [ -Inf , Inf ] but params.m = 3,
            %   ==> [ -Inf , Inf ; -Inf , Inf ; -Inf , Inf ]
            
            % input_bounds
            if ~isempty( obj.input_bounds ) && size( obj.input_bounds , 1 ) ~= obj.params.m
                obj.input_bounds = kron( ones( obj.params.m , 1 ) , obj.input_bounds );
            end
            
            % state_bounds
            if ~isempty( obj.state_bounds ) && size( obj.state_bounds , 1 ) ~= obj.params.n
                obj.state_bounds = kron( ones( obj.params.n , 1 ) , obj.state_bounds );
            end     
        end
        
        % get_costMatrices: Contructs the matrices for the mpc optim. problem
        function obj = get_costMatrices( obj )
            %get_costMatrices: Constructs cost the matrices for the mpc 
            % optimization problem.
            %   obj.cost has fields H, G, D, A, B, C, Q, R
            
            % define cost function matrices
            % Cost function is defined: U'HU + ( z0'G + Yr'D )U
            
            model = obj.model;
            
            % A
            N = size(model.A,1);
            A = sparse( N*(obj.horizon+1) , N );
            for i = 0 : obj.horizon
                A( (N*i + 1) : N*(i+1) , : ) = model.A^i ;
            end
            
            % B
            Bheight = N*(obj.horizon+1);
            Bcolwidth = size(model.B,2);
            Bcol = sparse( Bheight , Bcolwidth );    % first column of B matrix
            for i = 1 : obj.horizon
                Bcol( (N*i + 1) : N*(i+1) , : ) = model.A^(i-1) * model.B ;
            end
            
            Lshift = spdiags( ones( N*obj.horizon , 1 ) , -N , N*(obj.horizon+1) , N*(obj.horizon+1) );    % lower shift operator
            
            Bwidth = size(model.B,2)*(obj.horizon);    % total columns in B matrix
            Bblkwidth = obj.horizon;   % how many Bcol blocks wide B matrix is
            B = spalloc( Bheight , Bwidth , floor(Bheight * Bwidth / 2) ); % initialze sparse B matrix
            B(: , 1:Bcolwidth) = Bcol;
            for i = 2 : Bblkwidth
                B(: , (i-1)*Bcolwidth+1 : i*Bcolwidth) = Lshift * B(: , (i-2)*Bcolwidth+1 : (i-1)*Bcolwidth);
            end
            
            % C: matrix that projects lifted state into reference trajectory space
            C = kron( speye(obj.horizon+1) , obj.projmtx);
            nproj = size( obj.projmtx , 1 );
            
            % Q: Error magnitude penalty
            Q = kron( speye(obj.horizon+1) , eye(nproj) * obj.cost_running); % error magnitude penalty (running cost) (default 0.1)
            Q(end-nproj+1 : end , end-nproj+1 : end) = eye(nproj) * obj.cost_terminal;    % (terminal cost) (default 100)
            
            % R: Input magnitude penalty
            R = kron( speye(obj.horizon) , eye(model.params.m) * obj.cost_input );  % input magnitude penalty (for flaccy use 0.5e-2) (new videos used 0.5e-3)

            % H, G, D
            H = B' * C' * Q * C * B + R;
            G = 2 * A' * C' * Q * C * B;
            D = -2 * Q * C * B;
            
            % set outputs
            obj.cost.H = H; obj.cost.G = G; obj.cost.D = D; % constructed matrices
            obj.cost.A = A; obj.cost.B = B; obj.cost.C = C; obj.cost.Q = Q; obj.cost.R = R; % component matrices
        end
        
        % get_constraintMatrices: Constructs the constraint matrices
        function obj = get_constraintMatrices( obj )
            %get_constraintMatrices: Constructs the constraint matrices for
            % the mpc optimization problem.
            %   obj.constraints has fields L, M, F, E, (c?)
            %   F is for input constraints
            %   E is for state constraints
            
            % shorten some variable names
            Np = obj.horizon;     % steps in horizon
            params = obj.params;    % linear system model parameters
            cost = obj.cost;     % cost matrices
            
            F = []; E = [];     % initialize empty matrices
            c = [];
            
            % input_bounds
            if ~isempty( obj.input_bounds )
                num = 2*params.m;     % number of input bound constraints
                
                % F: input_bounds
                Fbounds_i = [ -speye(params.m) ; speye(params.m) ];    % diagonal element of F, for bounded inputs
                Fbounds = sparse( num * (Np+1) , size(cost.B,2) );  % no constraints, all zeros
                Fbounds( 1:num*Np , 1:Np*params.m ) = kron( speye(Np) , Fbounds_i );     % fill in nonzeros
                F = [ F ; Fbounds ];    % append matrix
                
                % E: input_bounds (just zeros)
                Ebounds = sparse( num * (Np+1) , size(cost.B,1) );  % no constraints, all zeros
                E = [ E ; Ebounds ];    % append matrix
                
                % c: input_bounds
                if isfield( obj.params , 'NLinput' )    % don't scale input if it's nonlinear
                    input_bounds_sc = obj.input_bounds;
                else
                    input_bounds_sc = obj.scaledown.u( obj.input_bounds' )';   % scaled down the input bounds
                end
                cbounds_i = [ -input_bounds_sc(:,1) ; input_bounds_sc(:,2) ]; % [ -umin ; umax ]
                cbounds = zeros( num * (Np+1) , 1);    % initialization
                cbounds(1 : num*Np) = kron( ones( Np , 1 ) , cbounds_i );     % fill in nonzeros
                c = [ c ; cbounds ];    % append vector
            end
            
            % input_slopeConst
            if ~isempty( obj.input_slopeConst )
                % F: input_slopeConst
                Fslope_i = speye(params.m);
                Fslope_neg = [ kron( speye(Np-1) , -Fslope_i ) , sparse( params.m * (Np-1) , params.m ) ];
                Fslope_pos = [ sparse( params.m * (Np-1) , params.m ) , kron( speye(Np-1) , Fslope_i ) ];
                Fslope_top = Fslope_neg + Fslope_pos;
                Fslope = [ Fslope_top ; -Fslope_top];
                F = [ F ; Fslope ];     % append matrix

                % E: input_slopeConst (just zeros)
                E = [ E ; sparse( 2 * params.m * (Np-1) , size(cost.B,1) ) ];
                
                % c: input_slopeConst
                if isfield( obj.params , 'NLinput' )    % don't scale slope if it's nonlinear
                    slope_lim = obj.input_slopeConst;
                else
                    slope_lim = obj.input_slopeConst * mean( params.scale.u_factor );  % scale down the 2nd deriv. limit
                end
                cslope_top = slope_lim * ones( params.m * (Np-1) , 1 );
                cslope = [ cslope_top ; cslope_top ];
                c = [ c ; cslope ];     % append vector
            end
            
            % input_smoothConst
            if ~isempty( obj.input_smoothConst )
                % F: input_smoothConst
                Fsmooth_i = speye(params.m);
                Fsmooth_lI = [ kron( speye(Np-2) , Fsmooth_i ) , sparse( params.m * (Np-2) , 2 * params.m ) ];
                Fsmooth_2I = [ sparse( params.m * (Np-2) , params.m ) , kron( speye(Np-2) , -2*Fslope_i ) , sparse( params.m * (Np-2) , params.m ) ];
                Fsmooth_rI = [ sparse( params.m * (Np-2) , 2 * params.m ) , kron( speye(Np-2) , Fslope_i ) ];
                Fsmooth_top = Fsmooth_lI + Fsmooth_2I + Fsmooth_rI;
                Fsmooth = [ Fsmooth_top ; -Fsmooth_top ];
                F = [ F ; Fsmooth ];
                
                % E: input_smoothConst
                E = [ E ; sparse( 2 * params.m * (Np-2) , size(cost.B,1) ) ];
                
                % c: input_smoothConst
                smooth_lim = params.Ts^2 * obj.input_smoothConst * mean( params.scale.u_factor );  % scale down the 2nd deriv. limit
                csmooth = smooth_lim * ones( size(Fsmooth,1) ,1);
                c = [ c ; csmooth ];
            end
            
            % state_bounds
            if ~isempty( obj.state_bounds )
                num = 2*params.n;   % number of state bound constraints
                
                % E: state_bounds
                Esbounds_i = [ -speye(params.n) ; speye(params.n) ];    % diagonal element of E, for bounding low dim. states (first n elements of lifted state)
                Esbounds = sparse( num * (Np+1) , size(cost.A,1) );  % no constraints, all zeros
                Esbounds( 1:num*(Np+1) , 1:(Np+1)*params.n ) = kron( speye(Np+1) , Esbounds_i );     % fill in nonzeros
                E = [ E ; Esbounds ];    % append matrix
                
                % F: state_bounds (all zeros)
                Fsbounds = zeros( size( Esbounds , 1 ) , size( cost.B , 2 ) );
                F = [ F ; Fsbounds ];    % append matrix
                
                % c: state_bounds
                state_bounds_sc = obj.scaledown.y( obj.state_bounds' )';    % scaled down state bounds
                csbounds_i = [ -state_bounds_sc(:,1) ; state_bounds_sc(:,2) ]; % [ -ymin ; ymax ]
                csbounds = kron( ones( Np+1 , 1 ) , csbounds_i );     % fill in nonzeros
                c = [ c ; csbounds ];    % append vector
            end
            
            % set outputs
            obj.constraints.F = F;
            obj.constraints.E = E;    
            obj.constraints.c = c;
            obj.constraints.L = F + E * cost.B;
            obj.constraints.M = E * cost.A;
        end
        
        % get_mpcInput: Solve the mpc problem to get the input over entire horizon
        function [ U , z ]= get_mpcInput( obj , traj , ref )
            %get_mpcInput: Soves the mpc problem to get the input over
            % entire horizon.
            %   traj - struct with fields y , u. Contains the measured
            %    states and inputs for the past ndelay+1 timesteps.
            %   ref - matrix containing the reference trajectory for the
            %    system over the horizon (one row per timestep).
            %   shape_bounds - [min_shape_parameters , max_shape_parameters] 
            %    This is only requred if system has shape constraints 
            %   (note: size is num of shape observables x 2)
            %   z - the lifted state at the current timestep
            
            % shorthand variable names
            Np = obj.horizon;       % steps in the horizon
            nd = obj.params.nd;     % number of delays
            
            % construct the current value of zeta
            [ ~ , zeta_temp ] = obj.get_zeta( traj );
            zeta = zeta_temp( end , : )';   % want most recent points
            
            % lift zeta
            z = obj.lift.full( zeta );
            
            % check that reference trajectory has correct dimensions
            if size( ref , 2 ) ~= size( obj.projmtx , 1 )
                error('Reference trajectory is not the correct dimension');
            elseif size( ref , 1 ) > Np + 1
                ref = ref( 1 : Np + 1 , : );    % remove points over horizon
            elseif size( ref , 1 ) < Np + 1
                ref_temp = kron( ones( Np+1 , 1 ) , ref(end,:) );
                ref_temp( 1 : size(ref,1) , : ) = ref;
                ref = ref_temp;     % repeat last point for remainer of horizon
            end
            
            % vectorize the reference trajectory
            Yr = reshape( ref' , [ ( Np + 1 ) * size(ref,2) , 1 ] );
            
            % setup matrices for gurobi solver
            H = obj.cost.H;     
            f = ( z' * obj.cost.G + Yr' * obj.cost.D )';
            A = obj.constraints.L;
            b = - obj.constraints.M * z + obj.constraints.c;
            
            % tack on "memory" constraint to fix initial input u_0
            Atack = [ [ speye( obj.params.m ) ; -speye( obj.params.m ) ] , sparse( 2*obj.params.m , size(A,2) - obj.params.m ) ];
%             Atack_bot = [ sparse( 2*obj.params.m , obj.params.m) , [ speye( obj.params.m ) ; -speye( obj.params.m ) ] , sparse( 2*obj.params.m , size(A,2) - 2*obj.params.m ) ];
%             Atack = [ Atack_top ; Atack_bot ];
            btack = [ traj.u(end,:)' ; -traj.u(end,:)' ];
            A = [A ; Atack];    % tack on memory constraint
            b = [b ; btack];
            
            % solve the MPC problem
%             Uvec = quadprog_gurobi( H , f , A , b );   % solve using gurobi (returns NaNs of cannot be solved)
            Uvec = quadprog( 2*H , f , A , b );     % solve using matlab
            
            % reshape the output so each input will have one row (first row equals current input)
            U = reshape( Uvec , [ obj.params.m , Np ] )';
        end
        
        % resample_ref: Resamples a reference trajectory
        function ref_resampled = resample_ref( obj, ref )
            %resample_ref: Resamples a reference trajectory at the system
            % sampling time.
            % ref - struct with fields:
            %   t - time vector
            %   y - trajectory vector
            
            tr = 0 : obj.params.Ts : ref.t(end);
            ref_resampled = interp1( ref.t , ref.y , tr );
        end
        
        % run_simulation: Runs a simulation of system under mpc controller
        function results = run_simulation( obj , ref , y0 , u0)
            %run_trial: Runs a simulation of system under mpc controller.
            %   Tries to follow the trajectory in ref and impose the
            %   shape constraints in shape_bounds.
            %   Assume ref and shape_bounds have same sampling frequency as
            %   sytem, and that they are already scaled to be consistent 
            %   with the lifted model.
            %   ref - struct containing reference trajectory with fields:
            %       t - vector of timesteps
            %       y - each row is a desired point at the corresponding timestep
            %   x0 - [1,n] initial condtion
            %   u0 - [1,m] initial input
            
            % shorthand
            nd = obj.params.nd;
            Np = obj.horizon;
            
            % set value of initial conditions to zero if none provided
            if nargin < 3
                y0 = zeros( nd+1 , obj.params.n );
                u0 = zeros( nd+1 , obj.params.m );
            elseif nargin < 4
                y0 = kron( ones( nd+1 , 1 ) , y0 );
                u0 = zeros( nd+1 , obj.params.m );
            else
                y0 = kron( ones( nd+1 , 1 ) , y0 );
                u0 = kron( ones( nd+1 , 1 ) , u0 );
            end
            
            % resample and scale the reference trajectory
            ref_Ts = obj.resample_ref( ref );
            ref_sc = obj.scaledown.y( ref_Ts );
            
            % set initial condition
            initial.y = y0; initial.u = u0;
            [ initial , zeta0 ] = obj.get_zeta( initial );    % LINE NOT NEEDED
            
            % initialize results struct
            results = struct;
            results.T = [ 0 ];
            results.U = [ u0( end , : ) ];
            results.Y = [ y0( end , : ) ];
            results.K = [ 0 ];
            results.R = [ ref.y(1,:) ];
            results.X = [ y0( end , : ) ];
            results.Z = [ obj.lift.full( zeta0' )' ]; % lifted states
            
            k = 1;
            while k < size( ref_sc , 1 )
                
                % current time
                t = k * obj.params.Ts;
                
                % get current state and input with delays
                if k == 1
                    current.y = obj.scaledown.y( y0 );   
                    current.u = obj.scaledown.u( u0 );  
                elseif k < nd + 1
                    y = [ y0( k : end-1 , : ) ; results.Y ];
                    u = [ u0( k : end-1 , : ) ; results.U ];
                    current.y = obj.scaledown.y( y );
                    current.u = obj.scaledown.u( u ); 
                else
                    y = results.Y( end - nd : end , : );
                    u = results.U( end - nd : end , : );
                    current.y = obj.scaledown.y( y ); 
                    current.u = obj.scaledown.u( u ); 
                end
                
                % isolate the reference trajectory over the horizon
                if k + Np <= size( ref_sc , 1 )
                    refhor = ref_sc( k : k + Np , :);
                else
                    refhor = ref_sc( k : end , : );     % repeat last entry
                end 
                
                % get optimal input over horizon
                [ U , z ] = obj.get_mpcInput( current , refhor );
                
                % if a solution was not found, break out of while loop
                if any( isnan(U) )
                    break;
                end
                
                % isolate input for this step (may need to make it U(1,:)
                u_kp1_sc = U( 2 , : );
                
                % scaleup the input for the results
                u_kp1 = obj.scaleup.u( u_kp1_sc )';
                
                % simulate the system over one time-step
                z_k = z;
                u_k_sc = obj.scaledown.u( results.U(end,:) );  % need to use previously calculated input NEED TO MAKE THIS CLEANER!!!
                z_kp1 = obj.model.A * z_k + obj.model.B * u_k_sc';
                x_kp1 = obj.model.C * z_kp1;
                y_kp1_sc = x_kp1;  % output juse scaled version of state since model was learned from observations
                y_kp1 = obj.scaleup.y( y_kp1_sc' )';  % scale output back up 
                
                % record updated results
                results.T = [ results.T ; t ];
                results.U = [ results.U ; u_kp1' ];
                results.Y = [ results.Y ; y_kp1' ];
                results.K = [ results.K ; k ];
                results.R = [ results.R ; obj.scaleup.y( ref_sc( k , : ) ) ];   % note that this is not scaled down
                results.X = [ results.X ; x_kp1' ];
                results.Z = [ results.Z ; z'  ]; % current lifted state
                
                k = k + 1;  % increment step counter
            end
        end
        
    end
end





















