from ROOT import TFile, TTree, gDirectory, TChain
from array import array
from math import tan,atan,e,cos,sin
from numpy import mean

MaxEntries = 10e9
subsetCuts = ""

def thresholdAlgorithm(tc_mipPt,theshold):
    countAboveThresh = 0
    for mipPt in tc_mipPt:
        if mipPt >= theshold:
            countAboveThresh += 1

    dataRate = (countAboveThresh+1)*bitsPerTC

    return dataRate


def readoutWaferIfTCAboveThresh(tc_mipPt, threshold):
    dataRate = bitsPerTC
    if max(tc_mipPt) > threshold:
        dataRate = len(tc_mipPt)*bitsPerTC

    return dataRate


def readoutHGROCIfTCAboveThresh(tc_mipPt, tc_perROC, threshold):
    dataRate = 0
    _tc_mipPt = list(tc_mipPt)
    
    for i in range(int(48/tc_perROC)):
        # print "-"*5, i
        # print i*tc_perROC,(i+1)*tc_perROC
        # print _tc_mipPt[i*tc_perROC:(i+1)*tc_perROC]
        # print max(_tc_mipPt[i*tc_perROC:(i+1)*tc_perROC]) > threshold
        if max(_tc_mipPt[i*tc_perROC:(i+1)*tc_perROC]) > threshold:
            dataRate += tc_perROC*bitsPerTC
        else:
            dataRate += bitsPerTC

    return dataRate



bitsPerTC = 8

inputFile = TFile( 'waferNtuple.root', 'read' )
inputTree = inputFile.Get("triggerNtuple")
nEventsPerWafer = eval(list(inputTree.GetListOfBranches())[-1].GetName().strip('mipPt'))+1
nEvents = nEventsPerWafer


inputTree.Draw(">>eventList",subsetCuts)
eventList = gDirectory.Get("eventList")
print eventList

print eventList.GetN()


f = TFile( 'test_miniNtuple.root', 'recreate' )
t2 = TTree( 'miniNtuple', 'miniNtuple' )

subdet    = array( 'i', [ 0 ] )
zside     = array( 'i', [ 0 ] )
layer     = array( 'i', [ 0 ] )
wafer     = array( 'i', [ 0 ] )
waferType = array( 'i', [ 0 ] )
x         = array( 'f', [ 0. ] )
y         = array( 'f', [ 0. ] )
z         = array( 'f', [ 0. ] )
AvgTrgOcc = array( 'f', [ 0. ] )
trgOcc    = array( 'f', nEvents*[ 0.] )
trg_bitRate_alg1    = array( 'i', nEvents*[ 0 ] )
trg_bitRate_alg2    = array( 'i', nEvents*[ 0 ] )
trg_bitRate_alg3    = array( 'i', nEvents*[ 0 ] )

t2.Branch('subdet'   , subdet   , 'subdet/I'          )
t2.Branch('zside'    , zside    , 'zside/I'           )
t2.Branch('layer'    , layer    , 'layer/I'           )
t2.Branch('wafer'    , wafer    , 'wafer/I'           )
t2.Branch('waferType', waferType, 'waferType/I'       )
t2.Branch('x'        , x        , 'x/F'               )
t2.Branch('y'        , y        , 'y/F'               )
t2.Branch('z'        , z        , 'z/F'               )
t2.Branch('AvgTrgOcc', AvgTrgOcc, 'AvgTrgOcc/F'       )
t2.Branch('trgOcc'   , trgOcc   , 'trgOcc[%i]/F'%nEvents )
t2.Branch('trg_bitRate_alg1'   , trg_bitRate_alg1   , 'trg_bitRate_alg1[%i]/I'%nEvents )
t2.Branch('trg_bitRate_alg2'   , trg_bitRate_alg2   , 'trg_bitRate_alg2[%i]/I'%nEvents )
t2.Branch('trg_bitRate_alg3'   , trg_bitRate_alg3   , 'trg_bitRate_alg3[%i]/I'%nEvents )


nWafers = eventList.GetN()

for iEntry in range(nWafers):
    inputTree.GetEntry(eventList.GetEntry(iEntry))
    if iEntry%1000==0:
        print "On wafer %i out of %i"%(iEntry, nWafers)

    subdet[0]    = inputTree.subdet     
    zside[0]     = inputTree.zside     
    layer[0]     = inputTree.layer     
    wafer[0]     = inputTree.wafer     
    waferType[0] = inputTree.waferType 
    x[0]         = inputTree.x         
    y[0]         = inputTree.y         
    z[0]         = inputTree.z         
    AvgTrgOcc[0] = inputTree.AvgTrgOcc 
    for evt in range(nEventsPerWafer):
        trgOcc[evt]    = inputTree.trgOcc[evt]
        mipPtData = eval("inputTree.mipPt%i"%evt)
        trg_bitRate_alg1[evt] = thresholdAlgorithm(mipPtData,2)
        trg_bitRate_alg2[evt] = readoutWaferIfTCAboveThresh(mipPtData,2)
        trg_bitRate_alg3[evt] = readoutHGROCIfTCAboveThresh(mipPtData,16,2)
    t2.Fill()
    
f.Write()
f.Close()
print "Done"



