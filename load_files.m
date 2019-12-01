function [t_sample,t_ref, traj_ref, traj_ref_sc] = load_files( ksysid )

%% load trajectory

% load in reference trajectory from file
[ reffile_name , reffile_path ] = uigetfile( 'ref-trajectories/*.mat' , 'Choose reference trajectory file...' );
temp = load( [reffile_path , reffile_name] );
ref = temp.ref;


%% set sampling time
t_sample = 0.5;  % Sampling Time, second


%% get time and trajectory reference
t_ref = 0:t_sample:ref.T;

traj_ref = interp1(ref.t,ref.y,t_ref);

traj_ref_sc = ksysid.scaledown.y( traj_ref ); % scale down reference

end

