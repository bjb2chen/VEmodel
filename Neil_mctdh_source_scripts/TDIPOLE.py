import sys
import pprint
import subprocess
import os
import shutil
import re
import json

def extract_lines_between_patterns(filename, start_pattern, end_pattern, encoding='ISO-8859-1'):
    selected_lines = []
    collecting = False

    with open(filename, 'r', encoding=encoding, errors='replace') as file:
        lines = file.readlines()

    # Find the last occurrence of the start pattern
    last_start_index = max(i for i, line in enumerate(lines) if start_pattern in line)

    for line in lines[last_start_index:]:  # Start iterating from the last occurrence of the start pattern
        if end_pattern in line:
            break  # Stop collecting when the end pattern is found
        if collecting:
            selected_lines.append(line)
        elif start_pattern in line:
            collecting = True  # Start collecting after the start pattern is found

    return selected_lines

def extract_same_state_transition_dipoles(selected_lines, nstate):
    same_state_transition_dipoles = {}
    start_pattern = "TRANSITION DIPOLES BETWEEN ADIABATS"
    data_started = False

    for TDIPOLEline in selected_lines:
        try:
            if TDIPOLEline[0:5].strip().isnumeric():
                state1 = int(TDIPOLEline[0:5].strip())
                print(state1)
                state2 = int(TDIPOLEline[5:10].strip())
    
                if state1 == state2:
                    x = float(TDIPOLEline[11:21].strip())
                    y = float(TDIPOLEline[22:31].strip())
                    z = float(TDIPOLEline[32:42].strip())
                    same_state_transition_dipoles[state1] = (x, y, z)
        except Exception as e:
            print(f"ERror processing line: {TDIPOLEline} - {e}")

    return same_state_transition_dipoles


def find_nstate(file_path, pattern='# of states in CI      = ', encoding="utf-8"):
    with open(file_path, 'r', encoding=encoding, errors='replace') as file:
        for line in file:
            if pattern in line:
                return int(line.split('=')[1].strip())
    return None  # Return None if the pattern is not found

# Example usage:

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 TDIPOLE.py <refG_file>")
        sys.exit(1)

    filename = sys.argv[1]
    nstate = find_nstate(filename)
    tdipole_block = extract_lines_between_patterns(filename,
    "TRANSITION DIPOLES BETWEEN DIABATS",
    "TRANSITION DIPOLES BETWEEN DIRECT MAX. DIABATS"
    )
    same_state_dipoles = extract_same_state_transition_dipoles(tdipole_block, nstate)
    pprint.pprint(tdipole_block)
    pprint.pprint(same_state_dipoles)


if __name__ == "__main__":
    main()