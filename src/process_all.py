
#!/usr/bin/env python
#-------------------------------------------------------------------
# This module is part of coor2graph
# Parallel process all adjacency matrix
#-------------------------------------------------------------------
#> author: Asdrubal Lozada
#  update: 4/11/2026
#  license: GPL-v3 

import argparse
from multiprocessing import Pool
from process_frame import process_frame
from functools import partial

def wrapper(frame, measure_type, measure_name, draw_freq, shape_intervals):
    draw_graph = (draw_freq > 0 and frame % draw_freq == 0)
    process_frame(frame, measure_type, measure_name, draw_graph, shape_intervals)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process multiple adjacency_#.csv files in parallel")
    parser.add_argument("--nframes", type=int, required=True, help="Total number of frames")
    parser.add_argument("--ncores", type=int, required=True, help="Number of cores to use")
    parser.add_argument("--measure_type", type=str, required=True,
                        choices=['degree', 'closeness', 'betweenness', 'katz', 'eigenvector'],
                        help="Metric type to calculate")
    parser.add_argument("--measure_name", type=str, required=True, help="Metric name (e.g., Degree)")
    parser.add_argument("--draw_freq", type=int, default=0,
                        help="Frequency for drawing the graph (e.g., 10 = every 10 frames, 0 = never)")
    parser.add_argument("--shape_intervals", type=str, required=True, help="File with intervals")

    args = parser.parse_args()

    shape_intervals = []

    with open(args.shape_intervals) as f:
        for line in f:
            start, end, shape = line.split()
            shape_intervals.append((int(start), int(end), shape))


    frames = range(1, args.nframes + 1)

    print(f"Processing {args.nframes} frames with {args.ncores} cores")
    print(f"Metric: {args.measure_type}, Graphing frequency: {args.draw_freq}")


    func = partial(wrapper,
                   measure_type=args.measure_type,
                   measure_name=args.measure_name,
                   draw_freq=args.draw_freq,
                   shape_intervals=shape_intervals)


    with Pool(processes=args.ncores) as pool:
        pool.map(func, frames)

