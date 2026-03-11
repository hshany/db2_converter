#!/usr/bin/env python
#Ryan G. Coleman

#import string
import sys
import optparse
import time
import io
from pathlib import Path
import logging
logger = logging.getLogger("mol2db2")

from db2_converter.mol2db2_py3_strain import mol2
from db2_converter.mol2db2_py3_strain import hierarchy
from db2_converter.mol2db2_py3_strain.hierarchy import TooBigError
from db2_converter.mol2db2 import solv
from db2_converter.mol2db2_py3_strain import clash
from db2_converter.mol2db2_py3_strain import hydrogens

def mol2db2_quick(mol2, options):
  if options.timeit:
    timeStart = time.time()
  else:
    timeStart = None
  mol2data = mol2
  mol2data.convertDockTypes(options.atomtypefile)
  mol2data.addColors(options.colortablefile)
  solvdata = solv.Solv(options.solvfile)
  if options.covalent:
    covAtomType,indicesList = mol2data.removeCovalentDummyAtom()
    mol2data.recolorCovalentAttachment(covAtomType)
    for index in indicesList:
      del solvdata.charge[index]
      del solvdata.polarSolv[index]
      del solvdata.surface[index]
      del solvdata.apolarSolv[index]
      del solvdata.solv[index]

  logger.debug("names:", mol2data.name, mol2data.protName, mol2data.smiles, \
      mol2data.longname)
  logger.debug("dock atom types:", mol2data.atomType)
  logger.debug("dock atom type numbers:", mol2data.dockNum)
  logger.debug("dock color type numbers:", mol2data.colorNum)
  clashDecider = clash.Clash(options.clashfile)
  hydrogenRotater = hydrogens.Hydrogens(options.hydrogenfile)
  if options.timeit:
    timeReadIn = time.time()
    logger.debug("time to read in files:", timeReadIn-timeStart)
  #0th step is to do the hydrogen operation directly on the mol2data.
  #done here so the hierarchy knows exactly how many conformations it must deal
  #with, not some estimate. and so each hierarchy has the max # of allowed sets
  #again, without guessing.
  if options.rotateh or options.reseth:
    hydrogenRotater.findTerminalHydrogens(mol2data)
    logger.debug(mol2data.hydrogensToRotate, " hydrogens need rotated")
    if options.reseth and 0 < mol2data.hydrogensToRotate:
      hydrogenRotater.resetHydrogens(mol2data)
    if options.rotateh and 0 < mol2data.hydrogensToRotate:
      mol2data = hydrogenRotater.rotateHydrogens(mol2data)
  if options.timeit:
    hydTime = time.time()
    logger.debug("time to move hydrogens:", hydTime-timeReadIn)
  logger.debug(len(mol2data.atomXyz), " conformations in input")

  def hierarchyDataGenerator(this_mol2data, depth=1):
    '''Generators to pipeline hierarchy generation'''
    try:
      yield hierarchy.Hierarchy(
          this_mol2data, clashDecider, tolerance=options.tolerance,
          verbose=options.verbose, timeit=options.timeit,
          limitset=options.limitset, limitconf=options.limitconf,
          limitcoord=options.limitcoord, solvdata=solvdata)
    except TooBigError as limitError:
      if depth > options.maxrecursiondepth:
        raise
      breaks = hierarchy.computeBreaks(limitError, options)
      origConfsPer = max(len(this_mol2data.atomXyz)/(breaks + 1), 1)  # at least 1
      logger.debug("splitting original input conformations into", breaks + 1, "parts", "(current depth: ", depth, ')')
      for snap in range(breaks + 1):  # extra one to do leftovers
        newMol2data = this_mol2data.copy()  # copy orig before hydrogen rotations
        first = origConfsPer * snap
        last = origConfsPer * (snap + 1)
        newMol2data.keepConfsOnly(first, last)
        if len(newMol2data.atomXyz) > 0:
          subgen = hierarchyDataGenerator(newMol2data, depth=depth+1)
          for subhier in subgen:
              yield subhier
  
  if options.timeit:
    timeHier = time.time()
    logger.debug("time to (start) construction hierarchy (subtotal):", timeHier-hydTime)
  return timeStart, hierarchyDataGenerator(mol2data)
  ## core of quick, we tentatively do not adopt this ##
  # hierarchyDatas = hierarchyDataGenerator(mol2data)
  # sio = io.StringIO()
  # for data in hierarchyDatas:
    ## new function added to hierarchy class, writeFile, takes in a file handle instead of a file name so we can write to a stringIO file handle
    # data.writeFile(fileHandle=sio, verbose=options.verbose, timeit=options.timeit, limitset=options.limitset)
  # return sio.getvalue()
  #####################################################

def mol2db2(options):
  '''function that does all the actual work you may want to do to convert a
  mol2 file and solv file into a db2 file.'''
  if options.timeit:
    timeStart = time.time()
  else:
    timeStart = None
  mol2data = mol2.Mol2(options.mol2file, nameFileName=None)
  mol2data.convertDockTypes(options.atomtypefile)
  mol2data.addColors(options.colortablefile)
  solvdata = solv.Solv(options.solvfile)
  if options.covalent:
    covAtomType,indicesList = mol2data.removeCovalentDummyAtom()
    mol2data.recolorCovalentAttachment(covAtomType)
    for index in indicesList:
      del solvdata.charge[index]
      del solvdata.polarSolv[index]
      del solvdata.surface[index]
      del solvdata.apolarSolv[index]
      del solvdata.solv[index]

  logger.debug("names:", mol2data.name, mol2data.protName, mol2data.smiles, \
      mol2data.longname)
  logger.debug("dock atom types:", mol2data.atomType)
  logger.debug("dock atom type numbers:", mol2data.dockNum)
  logger.debug("dock color type numbers:", mol2data.colorNum)
  clashDecider = clash.Clash(options.clashfile)
  hydrogenRotater = hydrogens.Hydrogens(options.hydrogenfile)
  if options.timeit:
    timeReadIn = time.time()
    logger.debug("time to read in files:", timeReadIn-timeStart)
  #0th step is to do the hydrogen operation directly on the mol2data.
  #done here so the hierarchy knows exactly how many conformations it must deal
  #with, not some estimate. and so each hierarchy has the max # of allowed sets
  #again, without guessing.
  if options.rotateh or options.reseth:
    hydrogenRotater.findTerminalHydrogens(mol2data)
    logger.debug(mol2data.hydrogensToRotate, " hydrogens need rotated")
    if options.reseth and 0 < mol2data.hydrogensToRotate:
      hydrogenRotater.resetHydrogens(mol2data)
    if options.rotateh and 0 < mol2data.hydrogensToRotate:
      mol2data = hydrogenRotater.rotateHydrogens(mol2data)
  if options.timeit:
    hydTime = time.time()
    logger.debug("time to move hydrogens:", hydTime-timeReadIn)
  logger.debug(len(mol2data.atomXyz), " conformations in input")

  def hierarchyDataGenerator(this_mol2data, depth=1):
    '''Generators to pipeline hierarchy generation'''
    try:
      yield hierarchy.Hierarchy(
          this_mol2data, clashDecider, tolerance=options.tolerance,
          verbose=options.verbose, timeit=options.timeit,
          limitset=options.limitset, limitconf=options.limitconf,
          limitcoord=options.limitcoord, solvdata=solvdata)
    except TooBigError as limitError:
      if depth > options.maxrecursiondepth:
        raise
      breaks = hierarchy.computeBreaks(limitError, options)
      origConfsPer = max(len(this_mol2data.atomXyz)/(breaks + 1), 1)  # at least 1
      logger.debug("splitting original input conformations into", breaks + 1, "parts", "(current depth: ", depth, ')')
      for snap in range(breaks + 1):  # extra one to do leftovers
        newMol2data = this_mol2data.copy()  # copy orig before hydrogen rotations
        first = origConfsPer * snap
        last = origConfsPer * (snap + 1)
        newMol2data.keepConfsOnly(first, last)
        if len(newMol2data.atomXyz) > 0:
          subgen = hierarchyDataGenerator(newMol2data, depth=depth+1)
          for subhier in subgen:
              yield subhier
  if options.timeit:
    timeHier = time.time()
    logger.debug("time to (start) construction hierarchy (subtotal):", timeHier-hydTime)
  return timeStart, hierarchyDataGenerator(mol2data)

def mol2db2writeDb2(options, timeStart, hierarchyDatas):
  '''does the writing of the output files. separate so you can make but
  not write. requires timeStart (can be None if no times wanted) and
  all the hierarchyDatas from the mol2db2 function'''
  if options.timeit:
    timeHier = time.time()
  # Clear db2gzfile
  open(options.db2gzfile, 'wt').close()
  for hierarchyData in hierarchyDatas:
    hierarchyData.write(
        db2gzFileName=options.db2gzfile, verbose=options.verbose,
        timeit=options.timeit, limitset=options.limitset,
        writeMode='at')  # append so we don't overwrite
    yield hierarchyData
  if options.timeit:
    timeHierOut = time.time()
    logger.debug("time to save db2:", timeHierOut-timeHier)
  if options.timeit:
    timeFinal = time.time()
    logger.debug("time to do everything (total):", timeFinal-timeStart)

class Mol2db2args():
  def __init__(self):
    self.covalent = False # "make a DOCKovalent compatible db (assumes SiH2 dummy atom)"
    self.timeit = False # 
    self.mol2file = "db.mol2.gz" # "multimol2 input file"
    self.namefile = "name.txt" # "mol2 name file"
    self.solvfile = "db.solv" # "solv input file"
    self.db2gzfile = "db.db2.gz" # "db2.gz output file"
    self.atomtypefile = None # "atom type conversion file", None -> use default
    self.colortablefile = None # "color table conversion file", None -> use default
    self.clashfile = None # "color table conversion file", None -> use default
    self.hydrogenfile = None # "terminal hydrogen parameter file" None -> use default
    self.reseth = True # "reset the planar terminal hydrogen positions"
    self.rotateh = True # "rotate the terminal hydrogen positions"
    self.tolerance = 0.001 # "distance tolerance in angstroms"
    self.limitset = 9999999999 # "limit on the number of sets"
    self.limitconf = 9999999999 # "limit on the number of confs"
    self.limitcoord = 9999999999 # "limit on the number of coords"
    self.maxrecursiondepth = None # "Max recursive subdivision steps to take"
    self.selfrigid = [] # "User defined rigid body indexes"
    self.verbose = True

  def __repr__(self):
    description = "Convert multimol2 and solv files into db2 files"
    return description

def mol2db2_main(
    mol2,
    solvfile,
    db2gzfile,
    clashfile=Path(__file__).parent / "clashfile.txt",
    timeit=False,
    covalent=False,
    reseth=False,
    rotateh=False,
    selfrigid = []
):
  options = Mol2db2args()
  options.solvfile = solvfile
  options.clashfile = clashfile
  options.db2gzfile = db2gzfile
  options.timeit = timeit
  options.covalent = covalent
  options.reseth = reseth
  options.rotateh = rotateh
  options.selfrigid = selfrigid

  logger.debug("verbose debugging of mol2db2 requested.")
  # hierarchyDatas = mol2db2_quick(mol2, options)  # main program call
  timeStart, hierarchyDatas = mol2db2_quick(mol2, options)  # main program call
  finishedHierarchyDatas = mol2db2writeDb2(
      options, timeStart, hierarchyDatas)  # write output
  for hierarchyData in finishedHierarchyDatas:
    logger.debug("Cleaning up")
    del hierarchyData
