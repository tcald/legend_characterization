#!/usr/bin/env python

import os
import optparse

if __name__ == '__main__':

    parser = optparse.OptionParser('%prog file [...]')
    parser.add_option('-s', dest='serial',    type='str', default='')
    parser.add_option('-a', dest='author',    type='str', default='')
    parser.add_option('-d', dest='directory', type='str', default='')
    options,args = parser.parse_args()
    if options.serial == '':
        print 'must specify detector serial number with -s option'
        exit()
    if options.author == '':
        print 'must specify author with -a option'
        exit()
    if options.directory == '':
        print 'must specify directory with -d option'
        exit()
    if not os.path.exists(options.directory):
        print options.directory, " does not exist"
        exit()
    if options.directory[-1] != '/':
        options.directory += '/'

    infile  = open('report.tex', 'r')
    outfile = open(options.directory + options.serial + '.tex', 'w+')
    for line in infile:
        oline = line
        if line.find('Report:}') != -1:
            oline = line.replace('Report:}', 'Report: '+options.serial+'}')
        elif line.find('author{}') != -1:
            oline = line.replace('author{}', 'author{'+options.author+'}')
        elif line.find('__') != -1:
            oline = line.replace('__', '_'+options.serial+'_')
        outfile.write(oline)
    infile.close()
    outfile.close()

