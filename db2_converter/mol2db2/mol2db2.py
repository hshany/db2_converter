#!/usr/bin/env python

#Ryan G. Coleman
import time
import logging
import subprocess
from pathlib import Path
logger = logging.getLogger("mol2db2")

from db2_converter.mol2db2 import mol2
from db2_converter.mol2db2 import hierarchy
from db2_converter.mol2db2 import solv
from db2_converter.mol2db2 import clash
from db2_converter.mol2db2 import hydrogens
from db2_converter.mol2db2.hierarchy import TooBigError
from db2_converter.utils.utils import exist_size
from db2_converter.utils.rdkit_gen import smifile_to_sdffile
from db2_converter.utils.convert import convert_sdf_to_mol2

def mol2db2(options):
  '''function that does all the actual work you may want to do to convert a
  mol2 file and solv file into a db2 file.'''
  if options.timeit:
    timeStart = time.time()
  else:
    timeStart = None
  mol2data = mol2.Mol2(options.mol2file, nameFileName=options.namefile)
  mol2data.convertDockTypes(options.atomtypefile) #sybyl2dock.AtomConverter
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

  logger.debug(("names:", mol2data.name, mol2data.protName, mol2data.smiles, \
      mol2data.longname))
  logger.debug(("dock atom types:", mol2data.atomType))
  logger.debug(("dock atom type numbers:", mol2data.dockNum))
  logger.debug(("dock color type numbers:", mol2data.colorNum))
  clashDecider = clash.Clash(options.clashfile) # [1.70 ** 2] # for H-H
  hydrogenRotater = hydrogens.Hydrogens(options.hydrogenfile)
  """
  rulesDefault = [(1, "C.ar", "1", "S", "1", "H", "180"),
                (1, "C.ar", "1", "O", "1", "H", "180"),
                (1, "C.1", "1", "S", "1", "H", "-"),
                (1, "C.1", "1", "O", "1", "H", "-"),
                (1, "C", "1", "S", "1", "H", "120,240"),
                (1, "C", "1", "O", "1", "H", "120,240"),
                (2, "2", "2", "N", "1", "H", "-"),
                (2, "1", "2", "N", "1", "H", "180")]
  """
  if options.timeit:
    timeReadIn = time.time()
    logger.debug(("time to read in files:", timeReadIn-timeStart))
  #0th step is to do the hydrogen operation directly on the mol2data.
  #done here so the hierarchy knows exactly how many conformations it must deal
  #with, not some estimate. and so each hierarchy has the max # of allowed sets
  #again, without guessing.
  hydrogenRotAngles = [""] * len(mol2data.atomNum)
  len_combo = 1
  if options.rotateh or options.reseth: # Both default True
    hydrogenRotAngles = hydrogenRotater.findTerminalHydrogens(mol2data)
    rothydroidxs = [ i for i, ang in enumerate(hydrogenRotAngles) if ang != "-" ]
    if len(rothydroidxs) > 5:
      logger.info("Terminal Hydrogens > 5. Set rotateh and reseth to False.")
      options.rotateh = False
      options.reseth = False
  if options.rotateh or options.reseth: # Both default True
    logger.debug((mol2data.hydrogensToRotate, " hydrogens need rotated"))
    if options.reseth and 0 < mol2data.hydrogensToRotate:
      hydrogenRotater.resetHydrogens(mol2data)
    if options.rotateh and 0 < mol2data.hydrogensToRotate:
      mol2data, len_combo = hydrogenRotater.rotateHydrogens(mol2data)
  if options.timeit:
    hydTime = time.time()
    logger.debug(("time to move hydrogens:", hydTime-timeReadIn))
  logger.debug((len(mol2data.atomXyz), " conformations in input")) # mol2data.atomXyz[0] is a set of coordinates of a complete conf (a set)

  rigidcomponent = list()
  if options.selfrigid:
    rigidcomponent = list( [ int(i) for i in options.selfrigid ] )

  def hierarchyDataGenerator(this_mol2data, depth=1):
    '''Generators to pipeline hierarchy generation'''
    try:
      # limitset 2500000 default (bbe1a30c)
      # limitconf 200000 default (bbe1a30c)
      # limitcoord 1000000 default (bbe1a30c)
      yield hierarchy.Hierarchy(
          this_mol2data, clashDecider, tolerance=options.tolerance,
          timeit=options.timeit,
          limitset=options.limitset, limitconf=options.limitconf,
          limitcoord=options.limitcoord, solvdata=solvdata,
          hydrogenRotAngles=hydrogenRotAngles, len_combo = len_combo,
          rigidcomponent=rigidcomponent)
    except TooBigError as limitError:
      if depth > options.maxrecursiondepth:  # default is 1>1 and False
        raise
      breaks = hierarchy.computeBreaks(limitError, options) # break here only for split confs
      origConfsPer = max(len(this_mol2data.atomXyz)/(breaks + 1), 1)  # at least 1
      logger.debug(("splitting original input conformations into", breaks + 1, "parts", "(current depth: ", depth, ')'))
      for snap in range(breaks + 1):  # extra one to do leftovers
        newMol2data = this_mol2data.copy()  # copy orig before hydrogen rotations
        first = origConfsPer * snap
        last = origConfsPer * (snap + 1)
        newMol2data.keepConfsOnly(first, last)
        if len(newMol2data.atomXyz) > 0:
          subgen = hierarchyDataGenerator(newMol2data, depth=depth+1) # make hierarchydata separatedly, which can give rise to different roots 
          for subhier in subgen:
              yield subhier
  if options.timeit:
    timeHier = time.time()
    logger.debug(("time to (start) construction hierarchy (subtotal):", timeHier-hydTime))
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
    hierarchyData.write(db2gzFileName=options.db2gzfile, writeMode='at')  # append so we don't overwrite
    yield hierarchyData
  if options.timeit:
    timeHierOut = time.time()
    logger.debug(("time to save db2:", timeHierOut-timeHier))
  if options.timeit:
    timeFinal = time.time()
    logger.info((">>> time to do everything in mol2db2 (total):", timeFinal-timeStart))

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

  def __repr__(self):
    description = "Convert multimol2 and solv files into db2 files"
    return description

def mol2db2_main(
    mol2file,
    solvfile,
    namefile,
    db2gzfile,
    timeit=False,
    covalent=False,
    reseth=False,
    rotateh=False,
    selfrigid = []
):
  options = Mol2db2args()
  options.mol2file = mol2file
  options.solvfile = solvfile
  options.namefile = namefile
  options.db2gzfile = db2gzfile
  options.timeit = timeit
  options.covalent = covalent
  options.reseth = reseth
  options.rotateh = rotateh
  options.selfrigid = selfrigid

  logger.debug("verbose debugging of mol2db2 requested.")
  timeStart, hierarchyDatas = mol2db2(options)  # main program call
  finishedHierarchyDatas = mol2db2writeDb2(
      options, timeStart, hierarchyDatas)  # write output
  for hierarchyData in finishedHierarchyDatas:
    logger.debug("Cleaning up")
    del hierarchyData

def mol2db2_to_numhyds(smifile,mol2file="tmp_hyd.mol2",removemol2=True,namefile=None,atomtypefile=None,colortablefile=None,hydrogenfile=None):
  from db2_converter.config import config
  from db2_converter.utils.utils import run_external_command
  import os
  tmp_sdf = None
  if not exist_size(mol2file):
    if exist_size(smifile):
      tmp_sdf = Path(mol2file).with_suffix(".sdf")
      smifile_to_sdffile(smifile, tmp_sdf)
      convert_sdf_to_mol2(tmp_sdf, mol2file)
    else:
      return
  mol2data = mol2.Mol2(mol2file, nameFileName=namefile)
  mol2data.convertDockTypes(atomtypefile) #sybyl2dock.AtomConverter
  mol2data.addColors(colortablefile)
  hydrogenRotater = hydrogens.Hydrogens(hydrogenfile)
  hydrogenRotAngles = [""] * len(mol2data.atomNum)
  hydrogenRotAngles = hydrogenRotater.findTerminalHydrogens(mol2data)
  rothydroidxs = [ i for i, ang in enumerate(hydrogenRotAngles) if ang != "-" ]
  multiplier = 1
  for i in rothydroidxs:
    multiplier *= (len(hydrogenRotAngles[i].split(","))+1)
  if tmp_sdf and tmp_sdf.exists():
    os.remove(tmp_sdf)
  if removemol2: os.remove(mol2file)
  return len(rothydroidxs), multiplier # N of terminal rotatable hydrogens

# fix for py3-3.7
#if -1 != string.find(sys.argv[0], "mol2db2.py") and __name__ == '__main__':
# if -1 != sys.argv[0].find("mol2db2.py") and __name__ == '__main__':
#   mol2file = sys.argv[1]
#   solvfile = sys.argv[2]
#   namefile = sys.argv[3]
#   db2gzfile = sys.argv[4]
#   mol2db2_main(
#     mol2file=mol2file,
#     solvfile=solvfile,
#     namefile=namefile,
#     db2gzfile=db2gzfile,
#     )
