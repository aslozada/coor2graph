! PIPE the python networkx library
module networkx_module
  use utils_module
  implicit none
  private

  public :: build_script
 contains

    subroutine build_script(nf)
      character(len=9), parameter :: fileScript  = 'script.py'
      character(len=13), parameter :: fileMatrix = 'adjacency.csv'  
      integer :: unit
      integer, intent(in) :: nf
      character(len=5) :: Foo
      character(len=:), allocatable :: output
      
      ! test
      !call  openFile(unit,fileScript,2)

      unit = 33

      output = ''

      write(Foo,'(i5)') nf
      open(unit,file=fileScript,status='unknown')
      ! python header

      ! Write the content of the Python script
      write(unit, '(A)') '#!/usr/bin/env python'
      write(unit, '(A)') 'import argparse'
      write(unit, '(A)') 'import numpy as np'
      write(unit, '(A)') 'import pandas as pd'
      write(unit, '(A)') 'import networkx as nx'
      write(unit, '(A)') 'import matplotlib.pyplot as plt'
      write(unit, '(A)') 'import matplotlib.colors as mcolors'
      write(unit, '(A)') 'import sys'
      write(unit, '(A)') 'import os'
      write(unit, '(A)') ''
      write(unit, '(A)')""
      write(unit, '(A)')"sys.stdout = open(os.devnull, 'w')"
      write(unit, '(A)')""

      write(unit, '(A)') 'def draw(G, pos, measures, measure_name, output_file):'

   !   write(unit, '(A)') "    cbar = plt.colorbar(nodes)"
   !   write(unit, '(A)') "    cbar.set_label('Measure Value')"  
   !   write(unit, '(A)') "    cbar.ax.tick_params(labelsize=12)" 

   !   write(unit, '(A)') "    cbar.ax.yaxis.set_tick_params(width=2)" 
   !   write(unit, '(A)') "    cbar.ax.set_yticklabels([f'{tick:.2f}' for tick in cbar.get_ticks()])" 
      write(unit, '(A)') '    plt.figure(figsize=(10,10))'
      write(unit, '(A)') '    nodes = nx.draw_networkx_nodes(G, pos, node_size=100, cmap=plt.cm.plasma,'
      write(unit, '(A)') '                               node_color=list(measures.values()),'
      write(unit, '(A)') '                               nodelist=measures.keys())'
      write(unit, '(A)') '    nodes.set_norm(mcolors.SymLogNorm(linthresh=0.01, linscale=1, base=10))'
      write(unit, '(A)') '    labels = nx.draw_networkx_labels(G, pos, font_size=14, font_color="blue")'
      write(unit, '(A)') '    edges = nx.draw_networkx_edges(G, pos, edge_color="black")'
      write(unit, '(A)') ''
      write(unit, '(A)') '    plt.title(measure_name)'
      write(unit, '(A)') '    plt.colorbar(nodes)'
      write(unit, '(A)') '    plt.axis("off")'
      write(unit, '(A)') '    plt.savefig(output_file)  # Save the image with the specified file name'
      write(unit, '(A)') '    plt.close()  # Close the figure to free resources'
      write(unit, '(A)') ''
      write(unit, '(A)') 'if __name__ == "__main__":'
      write(unit, '(A)') '    parser = argparse.ArgumentParser(description="Draw a network and save the image with a specified file name.")'
      write(unit, '(A)') "    parser.add_argument('measure_type', type=str, choices=['degree', 'closeness', 'betweenness', 'katz', 'eigenvector'])"
      write(unit, '(A)') "    parser.add_argument('measure_name', type=str)"
      write(unit, '(A)') '    parser.add_argument("output_file", type=str, help="Name of the file to save the image")'
      write(unit, '(A)') ''
      write(unit, '(A)') '    args = parser.parse_args()'
      write(unit, '(A)') ''
      write(unit, '(A)') '    # Read the adjacency matrix'
      write(unit, '(A)') '    matrix_file = "adjacency.csv"'
      write(unit, '(A)') '    adj_matrix = pd.read_csv(matrix_file, header=None).values'
      write(unit, '(A)') ''
      write(unit, '(A)') '    # Create the graph'
      write(unit, '(A)') '    G = nx.from_numpy_array(adj_matrix)'
      write(unit, '(A)') '    pos = nx.spring_layout(G, seed=675)'
      write(unit, '(A)') ''
      write(unit, '(A)') "    if args.measure_type == 'degree':"
      write(unit, '(A)') "        measures = nx.degree_centrality(G)"
      output = 'degree_'//trim(adjustl(Foo))//'.txt'
      write(unit, '(A)') "        with open('"//output//"'"//", 'w') as f:"
      write(unit, '(A)') "            for node, value in measures.items():"
      write(unit, '(A)') "                f.write(f'{node}\t{value:.4f}\n')"
      write(unit, '(A)') ''
!      write(unit, '(A)') "        print(measures)"
      write(unit, '(A)') "    elif args.measure_type == 'closeness':"
      write(unit, '(A)') "        measures = nx.closeness_centrality(G)"
      output = 'closeness_'//trim(adjustl(Foo))//'.txt'
      write(unit, '(A)') "        with open('"//output//"'"//", 'w') as f:"
      write(unit, '(A)') "            for node, value in measures.items():"
      write(unit, '(A)') "                f.write(f'{node}\t{value:.4f}\n')"
      write(unit, '(A)') ''
 !     write(unit, '(A)') "        print(measures)"
      write(unit, '(A)') "    elif args.measure_type == 'betweenness':"
      write(unit, '(A)') "        measures = nx.betweenness_centrality(G)"
      output = 'betweenness_'//trim(adjustl(Foo))//'.txt'
      write(unit, '(A)') "        with open('"//output//"'"//", 'w') as f:"
      write(unit, '(A)') "            for node, value in measures.items():"
      write(unit, '(A)') "                f.write(f'{node}\t{value:.4f}\n')"
      write(unit, '(A)') ''
!      write(unit, '(A)') "        print(measures)"
      write(unit, '(A)') "    elif args.measure_type == 'katz':"
      write(unit, '(A)') "        measures = nx.katz_centrality(G)"
      output = 'katz_'//trim(adjustl(Foo))//'.txt'
      write(unit, '(A)') "        with open('"//output//"'"//", 'w') as f:"
      write(unit, '(A)') "            for node, value in measures.items():"
      write(unit, '(A)') "                f.write(f'{node}\t{value:.4f}\n')"
      write(unit, '(A)') ''
!     write(unit, '(A)') "        print(measures)"
      write(unit, '(A)') "    elif args.measure_type == 'eigenvector':"
      write(unit, '(A)') "        measures = nx.eigenvector_centrality(G)"
      output = 'eigenvector_'//trim(adjustl(Foo))//'.txt'
      write(unit, '(A)') "        with open('"//output//"'"//", 'w') as f:"
      write(unit, '(A)') "            for node, value in measures.items():"
      write(unit, '(A)') "                f.write(f'{node}\t{value:.4f}\n')"
      write(unit, '(A)') ''
!      write(unit, '(A)') "        print(measures)"
      write(unit, '(A)') ''
!@      write(unit, '(A)') "    with open(args.measure_type, 'w') as f:"
!@      write(unit, '(A)') "        f.write('Node\tmeasure\n')"
!@      write(unit, '(A)') ''
!@      write(unit, '(A)') "        for node, value in measures.items():"
!@      write(unit, '(A)') "            f.write(f'{node}\t{value:.4f}\n')"
      write(unit, '(A)') ''
      write(unit, '(A)') '    draw(G, pos, measures, args.measure_name, args.output_file)'

      ! Close the file
      close(unit)


    end subroutine build_script  

end module networkx_module
