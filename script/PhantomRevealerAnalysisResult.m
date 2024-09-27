classdef PhantomRevealerAnalysisResult
    % The analysis data that generated from PhantomRevealer.
    % 
    % Fields:
    %   analysis_type: Verify the type of analysis
    %   time: The time flag of file
    %   data_dict: The dictionary that stores all of the quantity
    %   axes: The axes of result array
    %   column_names: The column name for each quantity
    %   params: Other parameters of simulation
    % 
    % Flags:
    %   _cgs: Verify whether the data is in cgs
    %   column_units: The column_unit for each quantity
    
    properties
        time
        data_dict
        axes
        column_names
        params
    end
    
    methods
        function obj = PhantomRevealerAnalysisResult(time, data_dict, axes, column_names, params)
            % Constructor
            obj.time = time;
            obj.data_dict = data_dict;
            obj.axes = axes;
            obj.column_names = column_names;
            obj.params = params;
        end
    end
    
    methods(Static)
        function data = read_HDF5(filepath)
            % Read HDF5 file and return an instance of PhantomRevealerAnalysisResult
            time = h5read(filepath, '/time');
            data_dict = PhantomRevealerAnalysisResult.read_dict(filepath, '/data_dict', @int32);
            axes = PhantomRevealerAnalysisResult.read_dict(filepath, '/axes', @int32);
            column_names = PhantomRevealerAnalysisResult.read_dict(filepath, '/column_names', @int32);
            params = PhantomRevealerAnalysisResult.read_dict(filepath, '/params', @char);
            
            data = PhantomRevealerAnalysisResult(time, data_dict, axes, column_names, params);
        end
        
        function dictionary = read_dict(filepath, name, K)
            % Helper function to read dictionaries from HDF5 file
            info = h5info(filepath, name);
            dictionary = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for i = 1:numel(info.Groups)
                group = info.Groups(i);
                key = K(group.Name);
                val = h5read(filepath, fullfile(name, group.Name));
                if iscell(val)
                    dictionary(key) = cellfun(@(x) char(x), val, 'UniformOutput', false);
                elseif isnumeric(val)
                    dictionary(key) = transpose(val);  % Transpose because of the storage system
                else
                    dictionary(key) = val;
                end
            end
        end
    end
end