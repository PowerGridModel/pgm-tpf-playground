# SPDX-FileCopyrightText: Contributors to the Power Grid Model project <powergridmodel@lfenergy.org>
#
# SPDX-License-Identifier: MPL-2.0

import matplotlib.pyplot as plt
import numpy as np

class PlotData:
    def __init__(self, pgm_tpf_data, tpf_data_native):
        self.pgm_tpf_data = pgm_tpf_data
        self.tpf_data_native = tpf_data_native
        self.arrange_data()

    def arrange_data(self):
        self.n_nodes = []
        self.n_steps = []
        self.t_pgm = []
        self.t_tpf = []
        self.t_tpf_native = []

        for entry in self.pgm_tpf_data:
            self.n_nodes.append(entry["n_nodes"])
            self.n_steps.append(entry["n_steps"])
            self.t_pgm.append(entry["t_pgm"])
            self.t_tpf.append(entry["t_tpf"])

        for entry in self.tpf_data_native:
            self.t_tpf_native.append(entry["t_tpf_native"])

        self.n_nodes = np.array(self.n_nodes)
        self.n_steps = np.array(self.n_steps)
        self.t_pgm = np.array(self.t_pgm)
        self.t_tpf = np.array(self.t_tpf)
        self.t_tpf_native = np.array(self.t_tpf_native)

    def plot_n_steps_timing(self, fixed_n_node):
        mask = self.n_nodes == fixed_n_node
        plt.figure()
        plt.loglog(self.n_steps[mask], self.t_pgm[mask], label='t_pgm')
        plt.loglog(self.n_steps[mask], self.t_tpf[mask], label='t_tpf')
        plt.loglog(self.n_steps[mask], self.t_tpf_native[mask], label='t_tpf_native')
        plt.xlabel('n_steps')
        plt.ylabel('Timing (ms)')
        plt.title(f'Timing vs n_steps for n_nodes = {fixed_n_node}')
        plt.legend()
        plt.savefig(f'Fixed_{fixed_n_node}_nodes.pdf')
        plt.close()

    def plot_n_nodes_timing(self, fixed_n_step):
        mask = self.n_steps == fixed_n_step
        plt.figure()
        plt.loglog(self.n_nodes[mask], self.t_pgm[mask], label='t_pgm')
        plt.loglog(self.n_nodes[mask], self.t_tpf[mask], label='t_tpf')
        plt.loglog(self.n_nodes[mask], self.t_tpf_native[mask], label='t_tpf_native')
        plt.xlabel('n_nodes')
        plt.ylabel('Timing (ms)')
        plt.title(f'Timing vs n_nodes for n_steps = {fixed_n_step}')
        plt.legend()
        plt.savefig(f'Fixed_{fixed_n_step}_steps.pdf')
        plt.close()

    def plot_n_nodes_n_steps_timing(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.n_nodes, self.n_steps, self.t_pgm, label='t_pgm')
        ax.scatter(self.n_nodes, self.n_steps, self.t_tpf, label='t_tpf')
        ax.scatter(self.n_nodes, self.n_steps, self.t_tpf_native, label='t_tpf_native')
        ax.set_xlabel('n_nodes')
        ax.set_ylabel('n_steps')
        ax.set_zlabel('Timing (ms)')
        ax.set_title('Timing vs n_nodes and n_steps')
        ax.legend()
        plt.show()

import itertools
from exp_res import pgm_tpf_data, tpf_data_native

if __name__ == "__main__":
    plot_data = PlotData(pgm_tpf_data, tpf_data_native)
    fixed_n_nodes = [10, 50, 100, 1000, 5000]
    fixed_n_steps = [10, 100, 1000, 1000000]
    for fixed_n_node, fixed_n_step in itertools.product(fixed_n_nodes, fixed_n_steps):
        plot_data.plot_n_steps_timing(fixed_n_node)
        plot_data.plot_n_nodes_timing(fixed_n_step)
    # plot_data.plot_n_nodes_n_steps_timing()
