classdef data
    %data: This class stores function that can be used to mess with
    %experimental and simulated data files with fields t , y , u , (x) ,...
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = data(varargin)
            %data: Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
        end
    end
   
    methods(Static)    
        % resample (resamples data with a desired time step)
        function data_resampled = resample( data , Ts )
            %resample: resamples sim/exp data with a desired timestep
            %   data - struct with fields t, y, x (optional)
            %   Ts - the desired sampling period
            
            % get query points
            tq = ( data.t(1) : Ts : data.t(end) )';
            
            data_resampled.t = tq;
            data_resampled.u = interp1( data.t , data.u , tq );
            data_resampled.y = interp1( data.t , data.y , tq );
            if ismember( 'x' , fields(data) )
                data_resampled.x = interp1( data.t , data.x , tq );
            end
        end
        
        % chop (chop data into several trials)
        function data_chopped = chop( data , num , len )
            %chop: chop data into num trials of lenght len
            %   data - struct with fields t , y , (x)
            %   data_chopped - cell array containing the chopped datas
            
            % determine length of timestep
            Ts = mean( data.t(2:end) - data.t(1:end-1) ); % take mean in case they're not quite uniform
            
            % find maximum length of each chop for given num
            maxlen = data.t(end) / num;
            if len > maxlen
                len = maxlen;
                disp([ 'Maximum trial length is ' , num2str(maxlen) , 's. Using this value instead.' ]);
            end
            
            % set length of the chops in terms of time steps
            lenk = length( find( data.t < len ) );
            maxlenk = length( find( data.t < maxlen ) );
            
            data_chopped = cell( 1 , num );
            for i = 1 : num
                index = (i-1) * maxlenk + ( 1 : lenk );
                
                % chop the data
                data_chopped{i}.t = ( ( 1 : lenk ) - 1 ) * Ts;
                data_chopped{i}.y = data.y( index , : );
                data_chopped{i}.u = data.u( index , : );
                if ismember( 'x' , fields(data) )
                    data_chopped{i}.x = data.x( index , : );
                end
            end 
        end
        
        % merge (merge several data files into single file)
        function data_merged = merge_files
            %merge_files: Merge several data files into single file
            %   data_merged: cell array containing the contents of all of
            %     the data filed selected
            
            % select data file(s)
            [ datafile_name , datafile_path ] = uigetfile( '*.mat' , 'Choose data file(s) for merging...' , 'multiselect' , 'on' );
            
            % load in the data files
            if iscell( datafile_name )  % check if it's cell array
                data_merged = cell( 1 , length(datafile_name) );
                for i = 1 : length(datafile_name)
                    data_merged{i} = load( [datafile_path , datafile_name{i}] );
%                     temp = load( [datafile_path , datafile_name{i}] );
%                     data_merged{i} = temp.sysidData;
                end
            else    % if not a cell array, turn it into 1x1 cell array
                data_merged = cell(1,1);
                data_merged{1} = load( [datafile_path , datafile_name] );
%                 temp = load( [datafile_path , datafile_name] );
%                 data_merged{1} = temp.sysidData;
                disp('FYI, you only selected one file so your output cell array will have dimension 1.');
            end   
            
            % change names of fields from uppercase to lowercase
            for i = 1 : size( data_merged , 2 )
                if ismember( 'T' , fields( data_merged{i} ) )
                    data_merged{i}.t = data_merged{i}.T;
                    data_merged{i} = rmfield( data_merged{i} , 'T' );
                end
                if ismember( 'Y' , fields( data_merged{i} ) )
                    data_merged{i}.y = data_merged{i}.Y;
                    data_merged{i} = rmfield( data_merged{i} , 'Y' );
                end
                if ismember( 'U' , fields( data_merged{i} ) )
                    data_merged{i}.u = data_merged{i}.U;
                    data_merged{i} = rmfield( data_merged{i} , 'U' );
                end
                if ismember( 'x' , fields( data_merged{i} ) )
                    data_merged{i}.y = data_merged{i}.x;
%                     data_merged{i} = rmfield( data_merged{i} , 'y' );
                end
            end
        end
        
        % get_data4sysid (save a file that can be used for sysid)
        function data4sysid = get_data4sysid( train , val , saveon , name )
            %get_data4sysid: Create a data structure with 'train' and 'val'
            % fields, which are requred to perform sysid.
            %   train - cell array containing training data (should get it
            %      from data.merge). Use [] if you want to select files.
            %   val - cell array containing validation data (should get it
            %      from data.chop). Use [] if you want to select files.
            %   saveon - (optional) if true, will save the output as a .mat 
            %      file in the 'datafiles' folder. False by default.
            %   name - string. Will be preappended to filedname if saveon
            %       is true.
            %   data4sysid - struct with fields 'train' and 'val', which
            %       are themselves cell arrays containing data for individual
            %       trials
            
            if nargin < 3
                saveon = false; % don't save output by default
            elseif nargin < 4
                name = [];  % name should be empty
            else
                name = [ name , '_' ];
            end
            
            % go get training and validation files if none provided
            if isempty( train )
                train = data.merge; % get training data
            elseif ~iscell( train ) % if 'train' is not a cell array, make it a 1x1 cell array
                train_temp = train;
                train = cell(1,1);
                train{1} = train_temp;
            end
            if isempty( val )
                val = data.merge;   % get validation data
            elseif ~iscell( val ) % if 'val' is not a cell array, make it a 1x1 cell array
                val_temp = val;
                val = cell(1,1);
                val{1} = val_temp;
            end
            
            % set output
            data4sysid.train = train;
            data4sysid.val = val;
            
            % save output
            if saveon
                dateString = datestr(now , 'yyyy-mm-dd_HH-MM'); % current date/time
                fname = [ 'datafiles' , filesep , name , 'train-', num2str( length(train) ) , '_val-' , num2str( length(val) ) , '_' , dateString , '.mat' ];
                unique_fname = auto_rename( fname , '(0)' );
                save( unique_fname , '-struct' ,'data4sysid' );
            end
        end
        
    end
end

