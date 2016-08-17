#!/usr/bin/env python

import sys
import os
import subprocess


def main():
    input = sys.argv[1]
    output = sys.argv[2]
    predicted_output, extension = os.path.splitext(input)
    predicted_output += '_0.5_0.9.pdf'
    print input
    print output
    print predicted_output
    status = subprocess.call(['python', 'museqportrait/src/museqportraitbuilder.py', input])
    print status
    if status == 0:
        subprocess.call(['ln', '-s', predicted_output, output])


if __name__ == '__main__':
    main()
