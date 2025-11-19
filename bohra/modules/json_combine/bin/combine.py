#!/usr/bin/env python3
import sys,  json

output_file = sys.argv[1]

input_file = sys.argv[2:]

list_of_dicts = []
for file in input_file:
    with open(file, 'r') as f:
        data = json.load(f)
        list_of_dicts.append(data)
        # print(data)
print(list_of_dicts)
# with open(output_file, 'w') as f:
#     json.dump(list_of_dicts, f, indent=4)