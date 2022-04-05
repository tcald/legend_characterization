#!/usr/bin/env python3

import os
import optparse
import platform
import distutils.spawn
from multiprocessing.dummy import Pool

rdir = []
pdir = ['.']
soname = ''
ldir = 'build/'
ddir = ldir
bdir = 'bin/'

cflags = '-std=c++14 -m64 -fPIC -Wall -pthread '
nvflags = cflags
include = '-I./ -I$ROOTSYS/include -I$CLHEP_INCLUDE_DIR -I/home/tcald/sw/jsoncpp/include '
libs = '-L$ROOTSYS/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -lMinuit -L/home/tcald/sw/jsoncpp_build/lib -ljsoncpp '

include += '-I$MGDODIR/Transforms -I$MGDODIR/MJDB -I$MGDODIR/Base -I$MGDODIR/Root -I$MGDODIR/Majorana -I$MGDODIR/tam/inc -I$MGDODIR/tam/include -I$MGDODIR/Tabree '
libs += '-L$MGDODIR/lib -l MGDOTransforms -lMGDOMJDB -lMGDORoot -lMGDOBase -lMGDOMajorana -lMGDOTabree '

include += '-I$GATDIR/BaseClasses '
libs += '-L$GATDIR/lib -lGATBaseClasses '

include += '-I/home/tcald/sw/Gpufit '
libs += '-L/home/tcald/sw/Gpufit_build/Gpufit -lGpufit '
libs += '-L/home/tcald/sw/Gpufit_build/Cpufit -lCpufit '

include += '-I/usr/local/cuda-10.1/include '
libs += '-L/usr/local/cuda-10.1/lib64 -lcudart '

debug = False
clean = False
shared = ''

if not os.path.exists(ldir):
    os.system('mkdir ' + ldir)
if not os.path.exists(ddir):
    os.system('mkdir ' + ddir)
if not os.path.exists(bdir):
    os.system('mkdir ' + bdir)

nv = distutils.spawn.find_executable('nvcc')
if nv != None:
    cflags += '-D__CUDA '
    nvflags += '-D__CUDA '
    nvr = ['Wall', 'fPIC', 'pthread']
    for n in nvr:
        nvflags = nvflags.replace('-'+n+' ', '')

def g(c):
    print(c)
    s = os.system(c)
    return

def d(o):
    if os.system('objdump --syms '+o+' | grep debug > /dev/null 2>&1') == 0:
        return True
    return False

def c(r):
    global cflags, include, ldir, ddir 
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
    g('rootcling -f '+ddir+f+'Dict.cc -c '+include+' '+r+'.hh '+s+'LinkDef.hh') 
    g('g++ '+cflags+include+'-c '+ddir+f+'Dict.cc -o '+ldir+f+'Dict.o')
    return [sh, True]

def ex(p):
    global cflags, include, ldir, libs, shared, bdir
    if p[0:1] == './':
        p = p[2:len(p)]
    etime = 0.0
    if os.path.exists(bdir+p) and os.path.exists(ldir+p+'.o'):
        if d(ldir+p+'.o') == debug:
            etime=max(os.path.getmtime(bdir+p),os.path.getmtime(ldir+p+'.o'))
    elif os.path.exists(p+'.hh') and os.path.exists(ldir+p+'.o'):
        if d(ldir+p+'.o') == debug:
            etime=max(os.path.getmtime(ldir+p+'.o'),os.path.getmtime(p+'.hh'))
    ctime = os.path.getmtime(p+'.cc')
    if ctime > etime or clean:
        if os.path.exists(bdir + '/' + p):
            os.system('rm ' + bdir + '/' + p)
        g('g++ '+cflags+include+'-c '+p+'.cc -o '+ldir+p+'.o')
        if not os.path.exists(p+'.hh'):
            g('g++ '+cflags+ldir+p+'.o '+shared+libs+'-o '+bdir+'/'+p)
    if os.path.exists(p+'.hh') and libs.find(ldir+p+'.o ') < 0:
        libs = ldir+p+'.o ' + libs

def cu(c):
    global cflags, include, ldir, libs
    if c[0:1] == './':
        c = c[2:len(c)]
    etime = 0.0
    if os.path.exists(ldir+c+'.o'):
        if d(ldir+c+'.o') == debug:
            etime = os.path.getmtime(ldir+c+'.o')
    ctime = os.path.getmtime(c+'.cu')
    if ctime > etime or clean:
        if os.path.exists(ldir + '/' + c + '.o'):
            os.system('rm ' + ldir + '/' + c + '.o')
        g('nvcc '+nvflags+include+'-c '+c+'.cu -o '+ldir+c+'.o')
    if libs.find(ldir+c+'.o ') < 0:
        libs = ldir+c+'.o ' + libs
    
if __name__ == '__main__':

    parser = optparse.OptionParser('%prog file [...]')
    parser.add_option('-c', action='store_true', dest = 'clean', default=False)
    parser.add_option('-g', action='store_true', dest = 'debug', default=False)
    parser.add_option('-j', dest = 'threads',    type = 'int',   default = 1)
    parser.add_option('-C', action='store_true', dest = 'cpu',   default=False)
    parser.add_option('-G', action='store_true', dest = 'gpu',   default=False)
    options, args = parser.parse_args()
    if options.debug:
        cflags += '-g '
        debug = True
    else:
        cflags += '-O3 '
    if options.clean:
        clean = True
    if options.gpu and cflags.find('__CUDA') < 0:
        'nvcc not found in path'
        exit()
    elif options.cpu and cflags.find('__CUDA') >= 0:
        cflags  =  cflags.replace('-D__CUDA ', '')
        nvflags = nvflags.replace('-D__CUDA ', '')

    rc = []
    for r in rdir:
        for f in os.listdir(r):
            if f.endswith('.hh') and include.find('-I'+r+' ') < 0:
                include += '-I'+r+' '
            elif f.endswith('.cc'):
                rc.append(r+'/'+f[0:-3])
    
    cbin = []
    vbin = []
    prog = []
    cuda = []
    for p in pdir:
        for f in os.listdir(p):
            if f.endswith('.hh'):
                if include.find('-I'+p+' ') < 0:
                    include += '-I'+p+' '
                if os.path.exists(p+'/'+f[0:-3]+'.cc'):
                    cbin.append(p+'/'+f[0:-3])
                elif os.path.exists(p+'/'+f[0:-3]+'.cu'):
                    vbin.append(p+'/'+f[0:-3])
            elif f.endswith('.cc'):
                if not os.path.exists(p+'/'+f[0:-3]+'.hh'):
                    prog.append(p+'/'+f[0:-3])
            elif f.endswith('.cu'):
                if not os.path.exists(p+'/'+f[0:-3]+'.hh'):
                    cuda.append(p+'/'+f[0:-3])

    if clean:
        clf = ['Dict.cc', '.o', '.so', '.pcm']
        cdir = [ldir, ddir]
        for dr in cdir:
            for f in os.listdir('./' + dr + '/'):
                for cl in clf:
                    if f.endswith(cl):
                        os.system('rm ' + dr + '/' + f)

    bshared = False
    pool = Pool(options.threads)
    results = pool.map(c, rc)
    for result in results:
        shared += result[0]
        if result[1]:
            bshared = True
    if bshared:
        g('g++ -shared '+cflags+shared+'-o '+ldir+'lib'+soname+'.so')
    if cflags.find('__CUDA') >= 0 and len(vbin) > 0:
        pool.map(cu, vbin)
    if len(cbin) > 0:
        pool.map(ex, cbin)
    if cflags.find('__CUDA') >= 0 and len(cuda) > 0:
        pool.map(cu, cuda)
    if len(prog) > 0:
        pool.map(ex, prog)



