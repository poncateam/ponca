# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import pycuda.autoinit
import pycuda.driver as drv
import numpy as np
from PIL import Image
import time
import sys, getopt

import matplotlib.pyplot as plt

start_time = 0.0

def startTimer(funcname, quiet):
  global start_time
  if not quiet:
    sys.stdout.write(funcname)
    sys.stdout.flush()
  start_time = time.time()


def stopTimer(quiet):
  elapsed_time = time.time() - start_time
  if not quiet:
    print(f"DONE in {elapsed_time:.3f} seconds.")


################################################################################
# Load input ptx file
def loadKernel():
  module    = drv.module_from_file("ssgls.ptx")
  glsKernel = module.get_function("doGLS_kernel")
  return module, glsKernel

################################################################################
# Format data for easy access in cuda kernel
def initInputData(imgNormPath, imgwcPath): 
  i_normals   = Image.open (imgNormPath)
  i_positions = Image.open (imgwcPath)

  # Todo: check dimensions for depth, normals and positions
  (w,h)     = i_positions.size     
  l    = w * h

  normals   = np.array(i_normals).reshape(l, 3).astype(np.float32) / 255.0
  positions = np.array(i_positions).reshape(l, 3).astype(np.float32) / 255.0

  positions = positions * 2.0 - 1.0

  # flatten for CUDA
  normals   = normals.reshape(-1)
  positions = positions.reshape(-1)

  return w, h, normals, positions


################################################################################
# Generate queries : In this simple example we set one query per pixel
def initQueries(w, h):
  nbQueries = w * h
  queries = np.zeros(2 * nbQueries, dtype=np.float32)

  for j in range(h):
    for i in range(w):
      idx              = 2 * (i + j * w)
      queries[idx]     = i
      queries[idx + 1] = j

  return nbQueries, queries


def main(argv):
  scale       = 10 # neighborhood size in pixels
  imgNormPath = ''
  imgwcPath   = ''
  quiet       = False
  outputPath  = ''

  try:
    opts, args = getopt.getopt(
      argv,
      "s:p:n:o:q",
      ["scale=", "positions=", "normals=", "output=", "quiet"]
    )
  except getopt.GetoptError:
    print("ssgls.py -s <scale> -p <image path> -n <normals>")
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-h':
      print("ssgls.py -s <scale> -p <image path> -n <normals>")
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
  w, h, normals, positions = initInputData(imgNormPath, imgwcPath)
  stopTimer(quiet)

  nbQueries = w * h

  startTimer("Init result array memory ... ", quiet)
  result = np.zeros(nbQueries, dtype=np.float32)
  stopTimer(quiet)

  ################################################################################
  # Launch kernel
  blockSize  = 32
  gridWidth  = w // blockSize
  gridHeight = h // blockSize

  w          = np.int32(w)
  h          = np.int32(h)
  scale      = np.int32(scale)

  startTimer("Launch kernel .............. ", quiet)
  glsKernel(
    drv.In(np.array([w, h, scale])),
    drv.In(positions),
    drv.In(normals),
    drv.Out(result),
    block=(blockSize, blockSize, 1),
    grid=(gridWidth, gridHeight)
  )
  stopTimer(quiet)

  ################################################################################
  # Analyse results
  result[np.isnan(result)] = 0

  if not quiet:
    print("\n###########Configuration#############")
    print(f"Image size: {w} x {h}")
    print(f"NbQueries:  {nbQueries}")
    print(f"Scale:      {scale}")
    print("################Results################")
    print(f"Max value:  {result.max()}")
    print(f"Min value:  {result.min()}")
    print("#######################################")

  scaleFactor = max(abs(result.max()), abs(result.min()))

  if outputPath == '':
    plt.imshow( result.reshape(h, w),
      vmin  = -scaleFactor,
      vmax  =  scaleFactor,
      cmap  = 'seismic'
    )
    plt.show()
  else:
    plt.imsave(
      fname = outputPath,
      arr   = result.reshape(h, w),
      vmin  = -1.0,
      vmax  =  1.0,
      cmap  = 'seismic'
    )


if __name__ == "__main__":
  main(sys.argv[1:])
