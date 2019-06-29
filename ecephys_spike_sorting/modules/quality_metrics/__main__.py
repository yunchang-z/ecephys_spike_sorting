from argschema import ArgSchemaParser
import os
import logging
import time
import multiprocessing
import ctypes
import sys
import gc

import numpy as np
import pandas as pd

from ...common.utils import load_kilosort_data, get_spike_depths, printProgressBar
from ...common.epoch import get_epochs_from_nwb_file

from .metrics_parallel import calculate_metrics


def unpack_data(input_tuple):

    return np.ctypeslib.as_array(input_tuple[0].get_obj()).reshape(input_tuple[1])

def initializer(*args):

    global spike_times_
    spike_times_ = unpack_data(args[0])
    
    global spike_clusters_
    spike_clusters_ = unpack_data(args[1])
    
    global amplitudes_
    amplitudes_ = unpack_data(args[2])

    global channel_map_
    channel_map_ = unpack_data(args[3])

    global pc_features_
    pc_features_ = unpack_data(args[4])

    global pc_feature_ind_
    pc_feature_ind_ = unpack_data(args[5])

    global spike_depths_
    spike_depths_ = unpack_data(args[6])

    global counter
    counter = args[7]

def create_shared_array(arr):

    data_type = {np.dtype('uint32') : ctypes.c_uint32,
                 np.dtype('int32') : ctypes.c_int32,
                 np.dtype('uint64') : ctypes.c_uint64,
                 np.dtype('int64') : ctypes.c_int64,
                 np.dtype('float32') : ctypes.c_float,
                 np.dtype('float64') : ctypes.c_double
                 }[arr.dtype]
    
    shared_array_base = multiprocessing.Array(data_type, arr.size)
    buffer = np.frombuffer(shared_array_base.get_obj(), dtype=arr.dtype).reshape(arr.shape)
    np.copyto(buffer, arr)
    
    return shared_array_base, arr.shape

def worker(cluster_id, params):

    df = calculate_metrics(cluster_id,
                             spike_times_, 
                             spike_clusters_, 
                             amplitudes_, 
                             channel_map_, 
                             pc_features_, 
                             pc_feature_ind_, 
                             spike_depths_,
                             params)

    with counter.get_lock():
        counter.value += 1
        printProgressBar(counter.value + 1, np.max(spike_clusters_))
        sys.stdout.flush()

    return df
    

def calculate_quality_metrics(args):

    print('ecephys spike sorting: quality metrics module')

    start = time.time()

    print("Loading data...")

    try:
        spike_times, spike_clusters, spike_templates, amplitudes, templates, channel_map, clusterIDs, cluster_quality, pc_features, pc_feature_ind = \
                load_kilosort_data(args['directories']['kilosort_output_directory'], \
                    args['ephys_params']['sample_rate'], \
                    use_master_clock = False,
                    include_pcs = True)

        spike_depths = get_spike_depths(spike_clusters, pc_features, pc_feature_ind)
        
        shared_arrays = tuple( (create_shared_array(arr) for arr in [spike_times, 
                                                            spike_clusters, 
                                                            amplitudes,
                                                            channel_map,
                                                            pc_features,
                                                            pc_feature_ind,
                                                            spike_depths]) )
        
        del spike_times, spike_clusters, spike_templates, amplitudes, templates, channel_map, spike_depths
        del pc_features, pc_feature_ind
        gc.collect()

        counter = multiprocessing.Value('i',0)

        initargs = shared_arrays + (counter,)

        print("Launching multiprocessing pool...")
        with multiprocessing.Pool(processes=4, initializer=initializer, initargs=initargs) as pool:
            results = pool.starmap(worker, zip(clusterIDs, [args['quality_metrics_params']] * clusterIDs.size))

        metrics = pd.concat(results, ignore_index=True)

    except FileNotFoundError:
        
        execution_time = time.time() - start

        print(" Files not available.")

        return {"execution_time" : execution_time,
            "quality_metrics_output_file" : None} 

    output_file = args['quality_metrics_params']['quality_metrics_output_file']

    if os.path.exists(args['waveform_metrics']['waveform_metrics_file']):
        metrics = metrics.merge(pd.read_csv(args['waveform_metrics']['waveform_metrics_file'], index_col=0),
                     on='cluster_id',
                     suffixes=('_quality_metrics','_waveform_metrics'))

    print("Saving data...")
    metrics.to_csv(output_file)

    execution_time = time.time() - start

    print('total time: ' + str(np.around(execution_time,2)) + ' seconds')
    print()
    
    return {"execution_time" : execution_time,
            "quality_metrics_output_file" : output_file} # output manifest


def main():

    from ._schemas import InputParameters, OutputParameters

    mod = ArgSchemaParser(schema_type=InputParameters,
                          output_schema_type=OutputParameters)

    output = calculate_quality_metrics(mod.args)

    multiprocessing.log_to_stderr(logging.DEBUG)

    output.update({"input_parameters": mod.args})
    if "output_json" in mod.args:
        mod.output(output, indent=2)
    else:
        print(mod.get_output_json(output))


if __name__ == "__main__":
    main()
