#!/usr/bin/env python3
import json,sys,pathlib, subprocess,os


def _check_file_size(reads):
    """
    in order to not break downstream processes for standard bacterial 
    pipeline then the reads need to be < 20000000.
    input:
        reads: list of reads
        jobtype: type of pipeline to be run on sequence
        isolate: seq
    output:
        jobtype: if read size is ok, return as input otherwise return FAIL - not enough data in read file.
    """
    # print(f"Checking file size.")
    
    for r in reads:

        if os.path.getsize(r) < 20000000:
            
            return 'FAIL_READ_FILE_TOO_SMALL'
    return "FILES_OK"

reads = []
for i in sys.argv[2:]:
    reads.append(i)
readdir = pathlib.Path(reads[0]).parent
iso = sys.argv[1]



# print(_check_file_size(reads=reads, jobtype=job),end="")
jts = _check_file_size(reads=reads)

print(jts,end="")