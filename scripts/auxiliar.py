import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import shutil
import matplotlib.colors as colors
from datetime import datetime
from scipy.optimize import fsolve

# logging
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.INFO)


def find_files_in_dir(directory):
    """
    Find all files in a directroy and all subdirectories inside
    """
    file_list = []
    for root, dirs, files in os.walk(directory):
        # iterate over all files inside one directory
        for file in files:
            file_path = os.path.join(root, file)
            file_list.append(file_path)

    return file_list



def createdir(path):
    """
    Create a directory if doesn't exists
    """
    # checking if exists first
    if not os.path.exists(path):
        os.makedirs(os.path.join(path), exist_ok=True)    
        
        
        
def find_nearest_value(array, value):
    """
    Finding the nearest integer inside a array
    """
    idx = (np.abs(np.array(array) - value)).argmin()
    return array[idx]        



def transparent_cmap(cmap, ranges=[0,1]):
    '''
    Retuns a colormap object tuned to transparent
    '''
    
    ncolors = 256
    color_array = plt.get_cmap(cmap)(range(ncolors))
    color_array[:,-1] = np.linspace(*ranges, ncolors)
    
    # building the colormap
    return colors.LinearSegmentedColormap.from_list(name='cmap', colors=color_array)



def move_files(source_folder, destination_folder):
    '''
    Function to move files from one directory to another
    '''
    
    # iterating over all files inside the folder
    for filename in os.listdir(source_folder):
        source_path = os.path.join(source_folder, filename)
        
        # checking if is a file or a folder
        if os.path.isfile(source_path):
            destination_path = os.path.join(destination_folder, filename)
            
            # moving file by file
            shutil.move(source_path, destination_path)


            
def delete_directory(directory_path):
    '''
    A function that try to delete a file
    '''
    
    try:
        os.rmdir(directory_path)
        logger.debug(f"Directory '{directory_path}' deleted successfully.")
        
    except OSError as error:
        logger.debug(f"Error deleting directory '{directory_path}': {error}")


        
def params(n=15):
    '''
    Function to set standard parameters for matplotlib
    '''
    plt.rcParams['font.size'] = n
    plt.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
    plt.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
    plt.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
    plt.rcParams['axes.linewidth'] = 1.9
    plt.rcParams['figure.figsize'] = (13, 7)
    plt.rcParams['lines.linewidth'] = 4
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['xtick.major.size'] = 8
    plt.rcParams['ytick.major.size'] = 8
    plt.rcParams['xtick.major.width'] = 1.8
    plt.rcParams['ytick.major.width'] = 1.8   
    plt.rcParams['lines.markeredgewidth'] = 2
    pd.set_option('display.max_columns', None)



def create_cmap(cols):
    '''
    Create a colormap given an array of colors
    '''    
    return colors.LinearSegmentedColormap.from_list('',  cols)




c1 = (5/255,5/255,153/255)
c2 = (102/255,0/255,204/255)
c3 = (255/255,51/255,204/255)
c4 = (204/255,0/255,0/255)
c5 = (255/255,225/255,0/255)

predC = [c1, c2, c3, c4, c5]

def color_cr(x, COLORS=predC):
    
    '''
    Function to create a color gradient of 5 colors in this case
    
    Input
    ------------
    --n: float
            the value from 0 to 1 to assign a colour
            
    Output
    ------------
    --r, g, b: float
            the rgb values for the color to assign
    
    '''
    C1, C2, C3, C4, C5 = COLORS[0], COLORS[1], COLORS[2], COLORS[3], COLORS[4]
    
    if x >= 0 and x <= 1/4:
        xeff = x
        r = C1[0] * (1 - 4 * xeff) + C2[0] * 4 * xeff
        g = C1[1] * (1 - 4 * xeff) + C2[1] * 4 * xeff
        b = C1[2] * (1 - 4 * xeff) + C2[2] * 4 * xeff
    elif x > 1/4 and x <= 2/4:
        xeff = x - 1/4
        r = C2[0] * (1 - 4 * xeff) + C3[0] * 4 * xeff
        g = C2[1] * (1 - 4 * xeff) + C3[1] * 4 * xeff
        b = C2[2] * (1 - 4 * xeff) + C3[2] * 4 * xeff
    elif x > 2/4 and x < 3/4:
        xeff = x - 2/4
        r = C3[0] * (1 - 4 * xeff) + C4[0] * 4 * xeff
        g = C3[1] * (1 - 4 * xeff) + C4[1] * 4 * xeff
        b = C3[2] * (1 - 4 * xeff) + C4[2] * 4 * xeff
    elif x >= 3/4 and x <= 1:
        xeff = x - 3/4
        r = C4[0] * (1 - 4 * xeff) + C5[0] * 4 * xeff
        g = C4[1] * (1 - 4 * xeff) + C5[1] * 4 * xeff
        b = C4[2] * (1 - 4 * xeff) + C5[2] * 4 * xeff
    else:
        print('Input should be in range [0 , 1]')
        
    return (r, g, b)









