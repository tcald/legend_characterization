#!/usr/bin/env python

import os
import optparse
import platform
from multiprocessing.dummy import Pool

rdir = []
pdir = ['.']
soname = ''
ldir = './build/'
ddir = ldir
bdir = './bin/'

cflags = '-std=c++11 -m64 -fPIC -Wall -pthread '
include = '-I./ -I$ROOTSYS/include '
libs = '`root-config --libs` -ltbb -lMinuit -l jsoncpp '

include += '-I$MGDODIR/Transforms -I$MGDODIR/MJDB -I$MGDODIR/Base -I$MGDODIR/Root -I$MGDODIR/Majorana -I$MGDODIR/tam/inc -I$MGDODIR/tam/include -I$MGDODIR/Tabree '
libs += '-L$MGDODIR/lib -l MGDOTransforms -lMGDOMJDB -lMGDORoot -lMGDOBase -lMGDOMajorana -lMGDOTabree '

debug = False
clean = False

if not os.path.exists(ldir):
    os.system('mkdir ' + ldir)
if not os.path.exists(ddir):
    os.system('mkdir ' + ddir)
if not os.path.exists(bdir):
    os.system('mkdir ' + bdir)

def g(c):
    print c
    s = os.system(c)
    return

def d(o):
    if os.system('objdump --syms '+o+' | grep debug > /dev/null 2>&1') == 0:
        return True
    return False

def c(r):
    f = r.split('/')[-1]
    s = r[0:len(r)-len(f)]
    sh = ldir+f+'.o '+ldir+f+'Dict.o '
    ltime = 0.0
    if os.path.exists(ldir+f+'.o'):
        if d(ldir+f+'.o') == debug:
            ltime = os.path.getmtime(ldir+f+'.o')
    htime = os.path.getmtime(r+'.hh')
    ctime = os.path.getmtime(r+'.cc')
    if ltime > htime and ltime > ctime:
        return [sh, False]
    g('g++ '+cflags+include+'-c '+r+'.cc'+' -o '+ldir+f+'.o')
    g('rootcling -f '+ddir+f+'Dict.cc -c '+'-I'+include+' '+r+'.hh '+s+'LinkDef.hh') 
    g('g++ '+cflags+include+'-c '+ddir+f+'Dict.cc -o '+ldir+f+'Dict.o')
    return [sh, True]

def ex(p):
    etime = 0.0
    if os.path.exists(p) and os.path.exists(ldir+p+'.o'):
        if d(ldir+p+'.o') == debug:
            etime =max(os.path.getmtime(p),os.path.getmtime(ldir+p+'.o'))
    ctime = os.path.getmtime(p+'.cc')
    if ctime > etime or clean:
        os.system('rm ' + bdir + '/' + p)
        g('g++ '+cflags+include+'-c '+p+'.cc -o '+ldir+p+'.o')
        g('g++ '+cflags+ldir+p+'.o '+shared+libs+'-o '+bdir+'/'+p)

    
if __name__ == '__main__':

    parser = optparse.OptionParser('%prog file [...]')
    parser.add_option('-c', action = 'store_true', dest = 'clean')
    parser.add_option('-g', action = 'store_true', dest = 'debug')
    parser.add_option('-j', dest = 'threads', type = 'int', default = 1)
    options, args = parser.parse_args()
    if options.debug:
        cflags += '-g '
        debug = True
    else:
        cflags += '-O3 '
    if options.clean:
        clean = True

    rc = []
    for r in rdir:
        for f in os.listdir(r):
            if f.endswith('.hh') and include.find('-I'+r+' ') < 0:
                include += '-I'+r+' '
            elif f.endswith('.cc'):
                rc.append(r+'/'+f[0:-3])
    prog = []
    for p in pdir:
        for f in os.listdir(p):
            if f.endswith('.hh') and include.find('-I'+p+' ') < 0:
                include += '-I'+p+' '
            elif f.endswith('.cc'):
                prog.append(p+'/'+f[0:-3])

    if clean:
        clf = ['Dict.cc', '.o', '.so', '.pcm']
        cdir = [ldir, ddir]
        for dr in cdir:
            for f in os.listdir('./' + dr + '/'):
                for cl in clf:
                    if f.endswith(cl):
                        os.system('rm ' + dr + '/' + f)

    shared = ''
    bshared = False
    pool = Pool(options.threads)
    results = pool.map(c, rc)
    for result in results:
        shared += result[0]
        if result[1]:
            bshared = True
    if bshared:
        g('g++ -shared '+cflags+shared+'-o '+ldir+'lib'+soname+'.so')
    pool.map(ex, prog)



