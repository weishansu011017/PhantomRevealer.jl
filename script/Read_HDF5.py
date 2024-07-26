import numpy as np
import h5py

class PhantomRevealerAnalysisResult:
    '''
    The analysis data that genetared by PhantomRevealer.
    
    Field in Class:
        For data:
            analysis_type: Verify the type of analysis
            time: The time flag of file
            data_dict: The dictionary that stores all of the quantity (Dict(column_index(int) -> np.ndarray))
            axes: The axes of result array.
            column_names: The column name for each quantuty (Dict(column_index(int) -> str))
            params: Other parameters of simulation
            
        For Other flag:
            _cgs: Verify whether the data is in cgs.
            column_units: The column_unit for each quantity.

    
    '''
    def __init__(self, time, data_dict, axes, column_names, params) -> None:
        self.time = time
        self.data_dict = data_dict
        self.axes = axes
        self.column_names = column_names
        self.params = params
        
    def read_HDF5(filepath): 
        def read_dict(f, name, K=str):
            dictionary = {}
            g = f[name]
            for key,val in g.items():
                if type(val[()]) == np.bytes_:
                    dictionary[K(key)] = val[()].decode('utf-8')
                elif type(val[()]) == np.ndarray:
                    dictionary[K(key)] = val[()].copy().T      # Trasepose because of the storage system
                else:
                    dictionary[K(key)] = val[()]
            return dictionary 
            
        with h5py.File(filepath, 'r') as f:
            time = f['time'][()]
            data_dict = read_dict(f,'data_dict',int)
            axes = read_dict(f, 'axes', int)
            column_names = read_dict(f, 'column_names', int)
            params = read_dict(f, 'params', str)
            
            data = PhantomRevealerAnalysisResult(time, data_dict, axes, column_names, params)
            return data