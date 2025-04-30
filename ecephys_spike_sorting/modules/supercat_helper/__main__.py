from argschema import ArgSchemaParser
import os
import sys
import subprocess
import time
import shutil

import numpy as np
from pathlib import Path

from ...common.utils import read_probe_json, get_repo_commit_date_and_hash, rms
from ecephys_spike_sorting.scripts.helpers import SpikeGLX_utils


def run_supercat(args):

    print('ecephys spike sorting: use supercat to oncatenate multiple runs')
    print('args: ' + str(args))
    catGTPath = args['supercat_helper_params']['catGTPath']
    if sys.platform.startswith('win'):
        os_str = 'win'
        # build windows command line
        # catGTexe_fullpath = catGTPath.replace('\\', '/') + "/runit.bat"
        # call catGT directly with params. CatGT.log file will be saved lcoally
        # in current working directory (with the calling script)
        catGTexe_fullpath = catGTPath.replace('\\', '/') + "/CatGT"
    elif sys.platform.startswith('linux'):
        os_str = 'linux'
        catGTexe_fullpath = catGTPath.replace('\\', '/') + "/runit.sh"
    else:
        print('unknown system, cannot run CatGt')
    
    cmd_parts = list()
    
    cmd_parts.append(catGTexe_fullpath)
    #The new option -supercat={dir,run_ga}{dir,run_ga}... takes a list of elements (include the curly braces) that specify which runs to join and in what order (the order listed).
    cmd_parts.append('-supercat=' + args['supercat_helper_params']['supercat_options'])
    # which probes
    cmd_parts.append('-prb=' + args['supercat_helper_params']['probe_string'])
    # which streams to concatenate
    cmd_parts.append(args['supercat_helper_params']['stream_string'])
    # command string
    cmd_parts.append(args['supercat_helper_params']['cmdStr'])
    # output directory
    cmd_parts.append('-dest=' + args['directories']['extracted_data_directory'])

    # Process the command parts to create the command line string
    catGT_cmd = ' '        # use space as the separator for the command parts
    catGT_cmd = catGT_cmd.join(cmd_parts[1:len(cmd_parts)]) # these are the parameters
    if os_str=='linux':
        # enclose the params in single quotes, so curly braces will not be interpreted by Linux
        catGT_cmd = f"{cmd_parts[0]} '{catGT_cmd}'"
    else:
        catGT_cmd = f"{cmd_parts[0]} {catGT_cmd}"
    
    print('CatGT command line for supercat:' + catGT_cmd)
    
    start = time.time()
    subprocess.Popen(catGT_cmd,shell='False').wait()

    execution_time = time.time() - start
    
    # copy CatGT log file, which will be in the directory with the calling 
    # python script, to the destination directory
    logPath = os.getcwd()

    # Print the current directory
    print("Current working directory:", logPath)
    # logName = 'CatGT.log'
   
    # first_gate, last_gate = SpikeGLX_utils.ParseGateStr(args['catGT_helper_params']['gate_string'])
         
    # catgt_runName = 'catgt_' + args['catGT_helper_params']['run_name'] + '_g' + str(first_gate)

    
    # # build name for log copy
    # catgt_logName = catgt_runName
    # if 'ap' in args['catGT_helper_params']['stream_string']:
    #     prb_title = ParseProbeStr(args['catGT_helper_params']['probe_string'])
    #     catgt_logName = catgt_logName + '_prb' + prb_title
    # if 'ni' in args['catGT_helper_params']['stream_string']:
    #     catgt_logName = catgt_logName + '_ni'
    # catgt_logName = catgt_logName + '_CatGT.log'
    
    
    # catgt_runDir = os.path.join(args['directories']['extracted_data_directory'],catgt_runName)
    # shutil.copyfile(os.path.join(logPath,logName), \
    #                 os.path.join(catgt_runDir,catgt_logName))
    
    print('total time: ' + str(np.around(execution_time,2)) + ' seconds')
    
    return {"execution_time" : execution_time} # output manifest


def main():

    from ._schemas import InputParameters, OutputParameters

    """Main entry point:"""
    mod = ArgSchemaParser(schema_type=InputParameters,
                          output_schema_type=OutputParameters)
    output = run_supercat(mod.args)

    output.update({"input_parameters": mod.args})
    if "output_json" in mod.args:
        mod.output(output, indent=2)
    else:
        print(mod.get_output_json(output))

if __name__ == "__main__":
    main()
