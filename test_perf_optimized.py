#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
test_perf.py calculates the execution time and memory peak of a program given in argument.
For arguments, use the absolute path, or the relative path from the test_perf.py file.
Use -h or --help for more details.

Returns the execution time (in second) and the memory peak (in kB) of a specific program.

author  : Emmanuel Clostres
version : 20/12/2020
"""

import os
import psutil
import sys
import time
from threading import Thread

THREAD_MEM = True
MEM_PEAK = 0

def usage():
    """ Display the usage of the program in the terminal. """
    print("\nUsage :\n\n"
          "python3 test_perf.py [path to the program] [program arguments]\n"
          "Use the absolute path, or the relative path from the test_perf.py file\n"
          "or\n"
          "python3 test_perf.py [aguments]\n"
          "\n"
          "Arguments :\n"
          "  -h, --help                  how to use the program\n")

def print_error(msg: str):
    """
    Write a yellow error message in the terminal.
    :param msg: str message to display
    """
    print(f"\033[93mERROR\033[0m : {msg}")

def join_arguments(list_of_arguments: list) -> str:
    """
    Returns a string with all arguments separated by a space, from a list of arguments.
    :param list_of_arguments: [str]
    :return: str
    """
    aruments_str = ""
    for arument in list_of_arguments:
        aruments_str += f" {arument}"
    return aruments_str

def get_MEM(step=1):
    """
    Recovers and stores the used memory (in kio (kB)) at a constant time interval (default : 1 sec).
    :param step: time interval in seconds to measure the memory used
    """
    global THREAD_MEM
    global START_MEM
    global MEM_PEAK
    while THREAD_MEM:
        time.sleep(step)
        mem = psutil.virtual_memory()
        current_mem = int(mem.used / 1024)  # in kB
        used_mem = current_mem - START_MEM
        if used_mem > MEM_PEAK:
            MEM_PEAK = used_mem
    print(" > Exit MEM thread")

if __name__ == '__main__':
    # Retrieving arguments
    if len(sys.argv) == 2:
        if sys.argv[1][:1] == "-":
            if sys.argv[1] == "-h" or sys.argv[1] == "--help":
                usage()
                sys.exit(0)
            else:
                print_error(f"unknown argument '\033[92m{sys.argv[1]}\033[0m'")
                usage()
                sys.exit(2)
        else:
            program_path = sys.argv[1]
            program_params = ""
    elif len(sys.argv) > 2:
        program_path = sys.argv[1]
        program_params = join_arguments(sys.argv[2:])
    else:
        print_error("one or more mandatory arguments are missing !")
        usage()
        sys.exit(2)
    # Checks if the path exists
    if not os.path.exists(program_path):
        print_error(f"file '\033[92m{program_path}\033[0m' not found !")
        usage()
        sys.exit(2)
    # Exec time and MEM recording initialization
    START_MEM = int(psutil.virtual_memory().used / 1024)  # memory used expressed in kibi-octet (Kio) (or kiloByte (kB))
    mem_thread = Thread(target=get_MEM, args=(0.5,))
    print(" > Start MEM thread")
    mem_thread.start()
    start_time = time.time()
    # Launching the program
    cmd = f"python {program_path}{program_params}"
    print(f" > Commande to execute : '{cmd}'")
    os.system(cmd)
    # Exec time and pic memory calculation
    exec_time = time.time() - start_time
    THREAD_MEM = False
    mem_thread.join()
    # Display the result
    print(f" > execution time : {exec_time} sec")
    print(f" > memory peak    : {MEM_PEAK} kB")
