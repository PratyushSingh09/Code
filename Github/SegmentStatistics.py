from scipy.ndimage import prewitt
import numpy as np
from ImageTransform import EM_normalize
import mode
import hipd_interval

def get_boundaries_for_cluster(data,lab):
    '''For a given segment value, define the boundaries of that segment
        data: the segmented image
        lab: segment value 

        Returns boundary coordinates
    '''
    f=prewitt(np.ma.array(data,mask=data!=lab).mask)
    return np.argwhere(f>0)    


def all_clusters():
    '''This is pointless'''
    neighbours=set()

def get_neighbour_clusters(data,pt_arr):
    '''For a given segment boundary, return the neighbouring segments
    data: the segmented image
    pt_arr: coordinates of the segment boundary, obtained using get_boundaries_for_cluster

    Returns a set of neighbouring segments
    '''
    ne_x = [0,0,1,-1,1,1,-1,-1]
    ne_y = [1,-1,0,0,1,-1,1,-1]
    neighbours=set()
    for pt in pt_arr:
        nx=ne_x+pt[0]
        ny=ne_y+pt[1]

        for i in range(len(nx)):
            if(nx[i]>=0 and ny[i]>=0 and nx[i]<data.shape[0] and ny[i]<data.shape[1]):
                neighbours.add(data[int(nx[i])][int(ny[i])])
                # print(nx[i],ny[i])

    return neighbours
    
def clus_stats(im_data,k_data,log_scale=False):
    '''Evaluates statistics for all segments using the DEMs obtained from sparse DEM inversion. Outdated function using mean and std, not used'''
    a_list=[]

    a_mean_list=[]
    a_std_list=[]
    a_err_list=[]
    a_len_list=[]

    i=1    #Only used for bookeeping
    for c in range(int(np.max(k_data)+1)):
        print(i)
        i=i+1
        # if pair is None:
        #     print('Pair is none')
        # p_iter = iter(pair)
        # print(pair)
        # if(len(pair)==1):
            # continue
        # else:
        #     a = next(p_iter)
        #     b = next(p_iter)
        a=c
        a_mask = []
        aml=[]
        asl=[]
        ael=[]
        alen=[]
        for dem in im_data:
            am = dem[k_data==a]
            a_mask.append(am)
            alen.append(len(am))
            # a_mask = np.mean(a_mask,axis=1)
            if log_scale:
                a_mean = np.nanmean(np.log(am,where=am>0))
                a_std = np.nanstd(np.log(am,where=am>0))
            else:    
                a_mean = np.nanmean(am)
                a_std = np.nanstd(am)
            a_err = a_std/np.sqrt(len(am))
            aml.append(a_mean)
            asl.append(a_std)
            ael.append(a_err)

        a_mean_list.append(aml)
        a_std_list.append(asl)
        a_err_list.append(ael)
        
        a_list.append(a_mask)
        a_len_list.append(alen)
        # b_mean = np.nanmean(b_mask)
        # b_std = np.nanstd(b_mask)

        # snr = abs(a_mean-b_mean)/np.sqrt(a_std**2+b_std**2)
        # return a_mask, b_mask
        # snr = chisquare(np.median(a_mask),np.median(b_mask))
        # snr_list.append(snr)

    return a_list,a_mean_list,a_std_list,a_err_list,a_len_list

def clus_stats_hpd(im_data,k_data,log_scale=False,hdip=None,em_norm=0,aia_data=None,aia_tr=None,dlgT=0.1):
    '''Evaluates statistics for all segments using the DEMs obtained from sparse DEM inversion. This function uses mode and hpdi
    
    im_data: Pixel wise DEMs of a active region in some n temperature bins with shape (n, x,y). Can also use sunpy maps directly. Must be the same size as k_data, aia_data and aia_tr
    k_data: The segmented image obtained from KMeans or other segmenation algorithms. Must be the same size as im_data, aia_data and aia_tr
    neighouring_cells_only: Choose to use all segments or only neighbouring segments to obtain pairs of statistics
    log_scale: Use the log scale on the DEM maps
    hdip: hdi probability
    em_norm: EM normalization method. 0 means none, 1 means log of normalized data and 2 means normalizing log data. 2 doesn't work well
    aia_data: AIA maps in 6(or any other number) filters. Must be the same size as im_data, k_data and aia_tr
    aia_tr: Transformed AIA maps in 6(or any other number) filters. Must be the same size as im_data, aia_data and k_tr

    RETURNS:
    All lists with length N
    a_list: all the points within each segment in each DEM bin
    a_mean_list: List of all mode of all points within the segment in each temp bin
    a_std_list: List of all hpdi values
    a_err_err: List of all stderr values. Not used
    a_coord_list: List of all coordinates of all the points within each segment 
    aia_dat_list: List of all aia intensities of all the points within each segment
    aia_av_list: List of average of transformed aia intensities of all the points within each segment
    a_len_list: List of number of all pixels in each segment
    '''
    a_list=[]

    a_mean_list=[]
    a_std_list=[]
    a_err_list=[]
    a_len_list=[]
    aia_dat_list=[]
    aia_av_list=[]
    a_coord_list=[]
    i=1    #Only used for bookeeping
    for c in range(int(np.max(k_data)+1)):
        a=c
        ac = np.argwhere(k_data==a)
        if aia_data is not None:
            aiad = aia_data[ac[:,0],ac[:,1],:]
        if aia_tr is not None:
            aiatr = aia_tr[ac[:,0],ac[:,1],:]
            aiaav = np.nanmean(aiatr,axis=0)
        print(i)
        i=i+1
        # if pair is None:
        #     print('Pair is none')
        # p_iter = iter(pair)
        # print(pair)
        # if(len(pair)==1):
            # continue
        # else:
        #     a = next(p_iter)
        #     b = next(p_iter)
        a_mask = []
        # a_coord=[]
        aml=[]
        asl=[]
        ael=[]
        alen=[]
        for dem in im_data:
            am = dem[k_data==a]
            # a_mask = np.mean(a_mask,axis=1)
            if log_scale:
                if em_norm==2:#DO NOT USE
                    amlg = EM_normalize(np.log10(np.array(am),dtype=np.float64,casting='unsafe'),dlgT)
                    a_mask.append(amlg)
                    #fix
                    a_mean = mode.modalpoint(amlg)
                    # mhdi = hdi(amlg,skipna=True,multimodal=False,hdi_prob=hdip)
                    mhdi = hipd_interval.hipd_interval(amlg,clev=0.68,fsample=True)#fsample = true is important
                    # print(mhdi)
                    # a_std = np.sum(mhdi[:,1]-mhdi[:,0])    
                    a_std = mhdi[1]-mhdi[0]
                else:
                    a_mask.append(np.log10(am))
                    #fix
                    a_mean = mode.modalpoint(np.log10(am))
                    # mhdi = hdi(np.log10(am),skipna=True,multimodal=False,hdi_prob=hdip)
                    mhdi = hipd_interval.hipd_interval(np.log10(am),clev=0.68,fsample=True)#fsample = true is important

                    # print(mhdi)
                    # a_std = mhdi[-1,1]-mhdi[0,0]
                    a_std = mhdi[1]-mhdi[0]

            else:    
                a_mask.append(am)
                a_mean = mode.modalpoint(am)
                # mhdi = hdi(am,skipna=True,multimodal=False,hdi_prob=hdip)
                mhdi = hipd_interval.hipd_interval(np.log10(am),clev=0.68,fsample=True)#fsample = true is important

                # a_std = mhdi[-1,1]-mhdi[0,0]
                a_std = mhdi[1]-mhdi[0]


            
            # a_coord.append(ac)
            alen.append(len(am))
            a_err = a_std/np.sqrt(len(am))
            aml.append(a_mean)
            asl.append(a_std)
            ael.append(a_err)

        a_mean_list.append(aml)
        a_std_list.append(asl)
        a_err_list.append(ael)
        
        a_list.append(np.array(a_mask))
        a_len_list.append(alen)

        a_coord_list.append(ac)
        
        if aia_data is not None:
            aia_dat_list.append(aiad)
        if aia_tr is not None:
            aia_av_list.append(aiaav)
    return a_list,a_mean_list,a_std_list,a_err_list,a_coord_list,aia_dat_list,a_len_list,aia_av_list