# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import pycuda.autoinit
import pycuda.driver as drv
import numpy
import Image
import time
import sys, getopt

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from pycuda.compiler import SourceModule

  
  
  

start_time=0.0
  
def startTimer( funcname, quiet ):
  global start_time
  
  if not quiet:
    sys.stdout.write(funcname)
    sys.stdout.flush()  
  start_time = time.time()


def stopTimer(quiet):
  elapsed_time = time.time() - start_time
  if not quiet:
    print 'DONE in {0:.3f} seconds.'.format(elapsed_time)  




################################################################################
# Load input ptx file
def loadKernel():                   
  module = drv.module_from_file("ssgls.ptx")      
  glsKernel    = module.get_function("doGLS_kernel")
  
  return module, glsKernel
  
  
  
    
 
################################################################################ 
# Format data for easy access in cuda kernel
def initInputData(imgNormPath, imgwcPath): 
  i_normals   = Image.open (imgNormPath)
  i_positions = Image.open (imgwcPath)

  # Todo: check dimensions for depth, normals and positions
  (w,h)     = i_positions.size     
  length    = w*h;      

  normals   = numpy.array(i_normals.getdata()).reshape(3*length,1).astype(numpy.float32) / 255.0;
  positions = numpy.array(i_positions.getdata()).reshape(3*length,1).astype(numpy.float32) / 255.0;
  
  positions = positions * 2.0 - 1.0;
    
  return w, h, normals, positions
  
  
  
    
  
################################################################################
# Generate queries: In this simple example we set one query per pixel
def initQueries(w,h):
  nbQueries = w*h

  queries   = numpy.zeros(2*nbQueries).astype(numpy.float32);
  for j in range (0,h):
      for i in range (0,w):
          queries[2*(i + j*w)] = i
          queries[2*(i + j*w)+1] = j
  return nbQueries, queries




def main(argv):    
    

  scale       = 10 # neighborhood size in pixels
  imgNormPath = ''
  imgwcPath   = ''
  quiet       = False
  outputPath  = ''
  
  
  
  try:
    opts, args = getopt.getopt(argv,"s:p:n:o:q",["scale=","positions=","normals=","output=","quiet"])
  except getopt.GetoptError:
    print 'ssgl.py -s <scale> -p <positions> -n <normals>'
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print 'ssgl.py -s <scale> -p <positions> -n <normals>'
      sys.exit()
    elif opt in ("-s", "--scale"):
      scale = int(arg)
    elif opt in ("-p", "--imgpath"):
      imgwcPath = arg
    elif opt in ("-n", "--normals"):
      imgNormPath = arg
    elif opt in ("-o", "--output"):
      outputPath = arg
    elif opt in ("-q", "--quiet"):
      quiet = True
      
  
   
  
  
  startTimer("Load CUDA kernel ........... ", quiet)
  module, glsKernel = loadKernel()
  stopTimer(quiet)  
  
  startTimer("Init Input data ............ ", quiet)
  w, h, normals, positions = initInputData (imgNormPath, imgwcPath)
  stopTimer(quiet)    
  
  nbQueries  = w*h
    
  startTimer("Init result array memory ... ", quiet)
  result = numpy.zeros(nbQueries).astype(numpy.float32)
  stopTimer(quiet)



  ################################################################################
  # Launch kernel
  blockSize  = 32
  gridWidth  = w/blockSize
  gridHeight = h/blockSize
  
  w          = numpy.int32(w)
  h          = numpy.int32(h)
  scale      = numpy.int32(scale)

  startTimer("Launch kernel .............. ", quiet)
  glsKernel(
          drv.In(numpy.array([w, h, scale])),
          drv.In(positions),
          drv.In(normals),
          drv.Out(result),
          block = (blockSize,  blockSize,  1), 
          grid  = (gridWidth,gridHeight))
          
  stopTimer(quiet)  
  
  ################################################################################
  # Analyse results
  whereAreNaNs = numpy.isnan(result);
  result[whereAreNaNs] = 0;

  if not quiet:
    print '\n###########Configuration#############'
    print 'Image size: {} x {}'.format(w,h)
    print 'NbQueries:  {}'.format(nbQueries)
    print 'Scale:      {}'.format(scale)
    print '################Results################'
    print 'Max value:  {}'.format(result.max())
    print 'Min value:  {}'.format(result.min())
    print '#######################################'

  scaleFactor = max(abs(result.max()), abs(result.min())); #/ (2.0*scaleFactor) + 0.5) 
  
  if outputPath == '':
    plt.imshow(result.reshape(h,w), vmin=-scaleFactor, vmax=scaleFactor, cmap='seismic')
    plt.show()
  else:
    plt.imsave(fname=outputPath, arr=result.reshape(h,w), vmin=-1.0, vmax=1.0, cmap='seismic')  
  
if __name__ == "__main__":
  main(sys.argv[1:])

