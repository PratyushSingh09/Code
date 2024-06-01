import os
import pandas as pd
def absoluteFilePaths(directory):
    
    '''Obtain a 'generator' for each file within the given directory. Used by ret_files_in_dir function

    Someone on stackoverflow uploaded this (put the link here when I find it)'''
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            yield os.path.abspath(os.path.join(dirpath, f))

def ret_files_in_dir(directory, num=6):
    '''For a given directory, return all file addresses within it

    directory: the directory to get the files from
    num: number of files to return back, default = 6

    I use it to obtain AIA images in 6 filters to be used for segmenation. MAKE SURE THE DIRECTORY CONTAINS ONLY THESE <num> FILES AND NOTHING ELSE. There are better functions I've made that you can use tbh

    Returns list of file paths
    '''
    cal_loop_files=absoluteFilePaths(directory)

    cal_list = []
    for i in range(num):
        cal_list.append(next(cal_loop_files))
    return cal_list

def get_path(rootdir,wavelength,suf='.image.fits'):
    '''For a given directory and an AIA filter wavelength, return all the image paths within the directory with a given suffix
    rootdir: the directory to get the image paths from
    wavelength: The AIA wavelength to get the image paths for
    suf: common string that each fits image in the directory ends with. AIA images obtained from JSOC ends with '.image.fits', while the images I've created only end with '.fits'
    
    Returns list of image file paths
    
    '''
    file_list = os.listdir(rootdir)
    img_path_list1 = [rootdir+"\\"+file for file in file_list if file.endswith(str(wavelength)+suf)]

    return img_path_list1

def get_paths(rootdir,suf='.image.fits'):
    '''For a given directory, return dataframe for all image paths in 6 different AIA wavelengths
    rootdir: the directory to get the image paths from
    suf: common string that each fits image in the directory ends with. AIA images obtained from JSOC ends with '.image.fits', while the images I've created only end with '.fits'
    
    Returns a Pandas table of image file paths. Each column of the table corresponds to a different filter and each row is for a different timestamp
    '''
    file_list = os.listdir(rootdir)
    img_list94 = [rootdir+"/"+file for file in file_list if file.endswith('94'+suf)]
    img_list131 = [rootdir+"/"+file for file in file_list if file.endswith('131'+suf)]
    img_list171 = [rootdir+"/"+file for file in file_list if file.endswith('171'+suf)]
    img_list193 = [rootdir+"/"+file for file in file_list if file.endswith('193'+suf)]
    img_list211 = [rootdir+"/"+file for file in file_list if file.endswith('211'+suf)]
    img_list335 = [rootdir+"/"+file for file in file_list if file.endswith('335'+suf)]

    return pd.DataFrame((img_list94,img_list131,img_list171,img_list193,img_list211,img_list335)).T