# Need to remove the '1' from the case files:
# so that   scalar per node:          1  potential reg1_proc111.gpot
# becomes   scalar per node:             potential reg1_proc111.gpot
import os

for proc in range(1,900):

    fname = f'reg1_proc{proc}.case'
    fpath = "/Users/eaton/Downloads/download/download/"

    fullpath = f"{fpath}/{fname}"
    tmpfile =  f"{fullpath}_tmp"

    # temp output file:
    f = open(tmpfile, "w")



    infile = open(fullpath, 'r')
    Lines = infile.readlines()
    count = 0
    # Strips the newline character
    for line in Lines:
        count += 1

        if 'scalar per node' in line:

            ind = line.find('1')
            line = line[:ind] + line[ind+1:]
        f.write(line)

    f.close()

    os.rename(tmpfile, fullpath)
