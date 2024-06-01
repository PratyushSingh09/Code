import numpy as np
def EM_normalize_old(DEM,dlgT):
    '''I forgot what heppened here'''
    demXdlgT = DEM*dlgT
    # print(demXdlgT)
    # print(demXdlgT.sum(axis=0))
    DEM = np.array(DEM,dtype=np.float64)
    demXdlgT = np.array(demXdlgT,dtype=np.float64)
    # return DEM/(demXdlgT.sum(axis=0))
    rdem=[]
    for d in DEM:
        rdem.append(d/np.max(d))
    # return np.array(rdem,dtype=np.float64)
    return DEM  

def EM_normalize(DEM,dlgT):
    '''Normalise the pixel DEMs with the EM at those pixels
    
    DEM: the pixel wise DEMs in n temperature bins with shape (n,x,y)
    dlgT: The change in logT parameter with each bin
    '''
    demXdlgT = DEM*dlgT
    # print(demXdlgT)
    # print(demXdlgT.sum(axis=0))
    DEM = np.array(DEM,dtype=np.float64)
    demXdlgT = np.array(demXdlgT,dtype=np.float64)
    return DEM/(demXdlgT.sum(axis=0))
