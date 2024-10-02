# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

import subprocess
import itertools
import time
import sys

from tqdm import tqdm

'''
Usage: PGM_TPF_Hackathon_2024.[exe/o] n_feeder n_node_per_feeder cable_length_km_min cable_length_km_max load_p_w_max load_p_w_min pf n_step load_scaling_min load_scaling_max [seed] [print_result]
Default: PGM_TPF_Hackathon_2024.[exe/o] 3 5 0.1 1.0 1000.0 500.0 0.95 10 0.8 1.2 0 0
Output:
TPF initialization took 1000 us
TPF power flow calculation took 2000 us
'''

class Option:
    def __init__(self, n_nodes_per_feeder, n_feeders):
        self.n_nodes_per_feeder = n_nodes_per_feeder
        self.n_feeders = n_feeders

def experiment_tpf(n_steps, options, executable, num_simulations=10):
    # Default parameters
    cable_length_km_min = 0.1
    cable_length_km_max = 1.0
    load_p_w_max = 1000.0
    load_p_w_min = 500.0
    pf = 0.95
    load_scaling_min = 0.8
    load_scaling_max = 1.2
    seed = 0
    print_result = 0

    results = []

    # Calculate the total number of iterations for the progress bar
    total_iterations = num_simulations * sum(
        n_step * option.n_nodes_per_feeder * option.n_feeders for n_step in n_steps for option in options
    )
 

    # Initialize the progress bar
    with tqdm(total=total_iterations, desc="Running experiments") as pbar:
        for option, n_step in itertools.product(options, n_steps):
            n_feeder = option.n_feeders
            n_node_per_feeder = option.n_nodes_per_feeder

            total_initialization_time = 0
            total_power_flow_time = 0
            total_execution_time = 0
            
            for _ in range(num_simulations):
                # Construct the command
                command = [
                    executable,
                    str(n_feeder),
                    str(n_node_per_feeder),
                    str(cable_length_km_min),
                    str(cable_length_km_max),
                    str(load_p_w_max),
                    str(load_p_w_min),
                    str(pf),
                    str(n_step),
                    str(load_scaling_min),
                    str(load_scaling_max),
                    str(seed),
                    str(print_result)
                ]
                
                # Measure the time taken to run the command
                start_time = time.time()
                result = subprocess.run(command, capture_output=True, text=True)
                end_time = time.time()
                
                # Parse the output
                output_lines = result.stdout.splitlines()
                tpf_initialization_time = float(output_lines[0].split()[-2]) / 1_000 # in ms
                tpf_power_flow_time = float(output_lines[1].split()[-2]) / 1_000 # in ms
                
                # Accumulate the times
                total_initialization_time += tpf_initialization_time
                total_power_flow_time += tpf_power_flow_time
                total_execution_time += (end_time - start_time) * 1_000 # in ms
                
                # Update the progress bar
                tqdm.write(str(output_lines))
                pbar.update(n_step * option.n_nodes_per_feeder * option.n_feeders)
            
             
            # Calculate the average times
            avg_initialization_time = total_initialization_time / num_simulations
            avg_power_flow_time = total_power_flow_time / num_simulations
            avg_execution_time = total_execution_time / num_simulations
            
            # Store the result
            results.append({
                'Total number of nodes:': n_feeder * n_node_per_feeder,
                'Number of steps:': n_step,
                'avg_initialization_time': avg_initialization_time,
                'avg_power_flow_time': avg_power_flow_time,
                'avg_execution_time': avg_execution_time
            })
    return results

if __name__ == "__main__":
    n_steps_ = [10, 100, 1_000, 1_000_000]

    options_ = [
        Option(5, 2), 
        Option(10, 5), 
        Option(10, 10), 
        Option(5, 100),
        Option(10, 100),
        Option(50, 100),
    ]

    num_simulations_ = 10
    if len(sys.argv) > 1:
        num_simulations_ = int(sys.argv[1])
    
    res = experiment_tpf(n_steps = n_steps_, options = options_, executable = '../bin/tpf.exe', num_simulations = num_simulations_)
    for r in res:
        print(r)
