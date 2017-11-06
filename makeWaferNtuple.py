from ROOT import TFile, TTree
from array import array
from math import tan,atan,e,cos,sin
from numpy import mean

inputFile = TFile("/uscms_data/d2/dnoonan/HGCAL/CMSSW_9_3_0/src/L1Trigger/L1THGCal/test/ntuple_200PU.root","read")
tree = inputFile.Get("hgcalTriggerNtuplizer/HGCalTriggerNtuple")


## subdet
##   layer
##     wafer
##       x,y,z
##       type
##       perEventTrigOcc
##       perEventDAQOcc
##       AvgOcc
##       RMSOcc


miniNtupleData = {3: {},
                  4: {},
                  5: {}}


for i in range(1,29):
    miniNtupleData[3][i] = {}
    miniNtupleData[3][-1*i] = {}
for i in range(1,13):
    miniNtupleData[4][i] = {}
    miniNtupleData[4][-1*i] = {}
for i in range(1,13):
    miniNtupleData[5][i] = {}
    miniNtupleData[5][-1*i] = {}

subDetLayers = {3:28,
                4:12,
                5:13}

nEvents = tree.GetEntries()
#nEvents = 2
i_event = 0
for event in tree:

    if i_event >= nEvents: break



    tc_n = event.tc_n

    tc_cell = event.tc_cell
    tc_subdet = event.tc_subdet
    tc_layer  = event.tc_layer
    tc_wafer  = event.tc_wafer
    tc_wafertype  = event.tc_wafertype
    tc_energy = event.tc_energy
    tc_mipPt  = event.tc_mipPt

    tc_eta    = event.tc_eta
    tc_phi    = event.tc_phi
    tc_z      = event.tc_z
    tc_zside  = event.tc_zside

    
    print i_event


    for n in range(tc_n):
        subdet = tc_subdet[n]
        layer  = tc_layer[n]

        layer = layer*tc_zside[n]
        wafer  = tc_wafer[n]
        cell = tc_cell[n]
        if not miniNtupleData[subdet][layer].has_key(wafer):
            nTrigs = 48 if tc_subdet[n] < 5 else 144
            miniNtupleData[subdet][layer][wafer] = {"x":[],
                                                    "y":[],
                                                    "z":[],
                                                    "type":tc_wafertype[n],
                                                    "perEvtTrgOcc":[0]*nEvents,
                                                    "perEvtDAQOcc":[0]*nEvents,
                                                    "perEvtTrgData":[[-1.]*nTrigs]*nEvents,
                                                    }


        r = tc_z[n]*tan(atan(e**(-1*tc_eta[n]))*2)
        x = r*cos(tc_phi[n])
        y = r*sin(tc_phi[n])
            
        miniNtupleData[subdet][layer][wafer]['x'].append(x)
        miniNtupleData[subdet][layer][wafer]['y'].append(y)
        miniNtupleData[subdet][layer][wafer]['z'].append(tc_z[n])
        if tc_mipPt[n]>2:
            miniNtupleData[subdet][layer][wafer]["perEvtTrgOcc"][i_event] += 1

        miniNtupleData[subdet][layer][wafer]["perEvtTrgData"][i_event][cell] = tc_mipPt[n]

    i_event += 1


f = TFile( 'waferNtuple.root', 'recreate' )
t = TTree( 'triggerNtuple', 'triggerNtuple' )

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
mipPt = []
for i in range(nEvents):
    mipPt.append(array( 'f', 48*[0.] ))

t.Branch('subdet'   , subdet   , 'subdet/I'          )
t.Branch('zside'    , zside    , 'zside/I'           )
t.Branch('layer'    , layer    , 'layer/I'           )
t.Branch('wafer'    , wafer    , 'wafer/I'           )
t.Branch('waferType', waferType, 'waferType/I'       )
t.Branch('x'        , x        , 'x/F'               )
t.Branch('y'        , y        , 'y/F'               )
t.Branch('z'        , z        , 'z/F'               )
t.Branch('AvgTrgOcc', AvgTrgOcc, 'AvgTrgOcc/F'       )
t.Branch('trgOcc'   , trgOcc   , 'trgOcc[%i]/F'%nEvents )
for i in range(nEvents):
    t.Branch('mipPt%i'%i    , mipPt[i]    , 'mipPt%i[48]/F'%i  )


print "Filling"

for _subdet in [3,4,5]:
    layers = range(-1*subDetLayers[_subdet]+1,subDetLayers[_subdet])
    layers.remove(0)
    subdet[0] = _subdet
    for _layer in layers:
        layer[0] = _layer
        for _wafer in range(500):
            if not miniNtupleData[_subdet][_layer].has_key(_wafer): continue
            waferdata = miniNtupleData[_subdet][_layer][_wafer]
            wafer[0] = _wafer
            zside[0] = _layer/abs(_layer)
            waferType[0] = waferdata["type"]
            x[0] = mean(waferdata['x'])
            y[0] = mean(waferdata['y'])
            z[0] = mean(waferdata['z'])            
            AvgTrgOcc[0] = mean(waferdata['perEvtTrgOcc'])
            for evt in range(nEvents):
                trgOcc[evt] = waferdata['perEvtTrgOcc'][evt]
                for i_trg in range(48):
                    mipPt[evt][i_trg] = waferdata['perEvtTrgData'][evt][i_trg]
            t.Fill()


f.Write()
f.Close()
print "Done"



