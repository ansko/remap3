#!/usr/bin/env python3

import os
import shutil
import subprocess
import sys


def compile():
    '''.tex -> .pdf
    '''
    try:
        fname = sys.argv[1]
    except IndexError:
        print('Fname is not specified!')
        sys.exit()

    if fname.endswith('.'):
        tex_fname = fname[:-1] + '.tex'
        fname = fname[:-1]
    elif fname.endswith('.tex'):
        tex_fname = fname
        fname = fname[:-4]
    else:
        tex_fname = fname + '.tex'

    print('Compiling', tex_fname)

    try:
        subprocess.Popen(['pdflatex', tex_fname],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).wait(5)
        subprocess.Popen(['bibtex', fname],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).wait(5)
        subprocess.Popen(['pdflatex', tex_fname],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).wait(5)
        subprocess.Popen(['pdflatex', tex_fname],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).wait(5)
    except subprocess.TimeoutExpired:
        print('Incorrect .tex source!')
        sys.exit()

    print('Done!')

    for ext in ['aux', 'bbl', 'blg', 'log', 'dvi', 'out', 'toc']:
        try:
            os.remove("{0}.{1}".format(fname, ext))
        except FileNotFoundError:
            pass

    for fname in os.listdir():
        if fname.startswith('acs'):
            os.remove(fname)

if __name__ == '__main__':
    compile()
