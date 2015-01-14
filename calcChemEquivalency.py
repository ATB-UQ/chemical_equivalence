from NautyInterface import NautyInterface

def getChemEquivGroups(atbDataObj, log):
    nautyInterface = NautyInterface(atbDataObj=atbDataObj)
    nautyInterface.calcEquivGroups(log)
    
    
    