#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.optimize import leastsq
from numpy import asarray as ar,exp
from fitfuncs import *


# In[ ]:


def input_cubes():
    
    fits_filename = ('/Users/orion/phd_research/3C297/finalcube_3C297_mediansubtracted.fits')     
    hdul = fits.open(fits_filename)
    data_cube = hdul[0].data
    cube_header = hdul[0].header 

    cw = cube_header['CRVAL3']         
    mp = cube_header['CDELT3']       

    y = data_cube[:,33,33]
    x = np.arange(len(y))
    start_w = cw - cube_header['CRPIX3']*mp 
    x = start_w + x*mp
    
    
    err_filename = ('/Users/orion/phd_research/3C297/line_fitting/error_cubes/centralpix_error_cube.fits')
    ehdul = fits.open(err_filename1)
    err_cube = ehdul1[0].data
    err_header = ehdul1[0].header 
    
    return data_cube, err_cube, y, x


# In[ ]:


def mini_cube(z, data_cube, err_cube, x):
    #construct a mini-cube, setting a central pixel 
    #syntax:  [:, y1:y2, x1:x2]
    
    [central_x,central_y]= [40,37]
    mini_datacube = data_cube[:,central_y - 15:central_y + 15,central_x - 17:central_x + 17]
    mini_errorcube = err_cube[:,central_y - 15:central_y + 15,central_x - 17:central_x + 17]
    
    k = 1+z
    wave = (x/k)*10000
    select = (wave>6450) & (wave<6800)
    
    e_core = mini_errorcube[:,11,10]
    error_spec_core = e_core[select]
    err_wave = x[select]
    y = mini_datacube[:,11,10]  # core center
    data_spec = y[select]

    gauss = Gaussian1DKernel(stddev=2)
    err_spec_temp = convolve(error_spec_core, gauss)
    b, err_spec = error_scale(data_spec, err_spec_temp)
    
    return mini_datacube, err_spec


# In[ ]:


def fitting(z1, z2, x, mini_datacube, err_spec):
    
    ampHa_map = np.zeros((mini_datacube.shape[1],mini_datacube.shape[2]),dtype=np.float32)
    veln_map = np.zeros((mini_datacube.shape[1],mini_datacube.shape[2]),dtype=np.float32)
    velsign_map = np.zeros((mini_datacube.shape[1],mini_datacube.shape[2]),dtype=np.float32)
    ampHaB_map = np.zeros((mini_datacube.shape[1],mini_datacube.shape[2]),dtype=np.float32)
    velb_map = np.zeros((mini_datacube.shape[1],mini_datacube.shape[2]),dtype=np.float32)
    velsigb_map = np.zeros((mini_datacube.shape[1],mini_datacube.shape[2]),dtype=np.float32)
    ampNII_map = np.zeros((mini_datacube.shape[1],mini_datacube.shape[2]),dtype=np.float32)
    ampSII_map = np.zeros((mini_datacube.shape[1],mini_datacube.shape[2]),dtype=np.float32)

    
    ############################################## CORE ##############################################
    i0 = [100, 30, 20, 20, 0, 20, 20, -50, 600, 0, 10]

    for i in range(0,20):              #i == y-axis, j == x-axis
        for j in range(0,20):

            print ('Spaxel= [',j,',',i,']')
            y = mini_datacube[:,i,j]                       
            y[np.isnan(y)]=1e-15
            small_wave, small_spec = data_init(z1,x,y)
       
            #if np.min(small_spec/err_spec) > 0:
            fit_params, err_params = ltsq_mc_fitting(i0,small_wave,small_spec,err_spec)
            print ('Fit=',fit_params,'\n', 'MC error=',err_params)
            plot(small_wave,small_spec,err_spec,fit_params)
        
            (amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,amp_HaB,vel_HaB,vel_sigma_HaB,m,c) = fit_params

            
        ########################################## Narrow-line maps
            if (amp_Ha > 0) and (amp_Ha < 300): 
                ampHa_map[j,i] = amp_Ha    
                    
                if (vel > -300) and (vel < 300):
                    veln_map[j,i] = vel
                else:
                    veln_map[j,i] = 0

        
                if (vel_sigma > -500) and (vel_sigma < 500):
                    velsign_map[j,i] = vel_sigma
                else:
                    velsign_map[j,i] = 0
                
                if (amp_NII6585 > 0) and (amp_NII6585 < 300):
                    ampNII_map[j,i] = amp_NII6585
                else:
                    ampNII_map[j,i] = 0
                
                if (amp_SII6716 > 0) and (amp_SII6716 < 300) and (amp_SII6730 > 0) and (amp_SII6730 < 300):
                    ampSII_map[j,i] = amp_SII6716 + amp_SII6730
                else:
                    ampSII_map[j,i] = 0 
              
            else:
                ampHa_map[j,i] = 0
       
    
        ########################################## Broad-line maps    
            if (amp_HaB > 0) and (amp_HaB < 100):        
                ampHaB_map[j,i] = amp_HaB    
            
                if (vel_HaB > -400) and (vel_HaB < 400):
                    velb_map[j,i] = vel_HaB
                else:
                    velb_map[j,i] = 0
        
                if (vel_sigma_HaB > -1500) and (vel_sigma_HaB < 1500):
                    velsigb_map[j,i] = vel_sigma_HaB
                else:
                    velsigb_map[j,i] = 0
              
            else:
                ampHaB_map[j,i] = 0
                
    
    
    ############################################## NORTHERN ARC ##############################################   
    p0 = [100, 30, 20, 20, 0, 20, 0, 10]
        
    for i in range(20,30):        
        for j in range(20,30):        
        
            print ('Spaxel= [',j,',',i,']')
            y = mini_datacube[:,i,j]                       
            y[np.isnan(y)]=1e-15
            small_wave, small_spec = data_init(z2,x,y)    
    
            fit_params, err_params = ltsq_mc_fitting_Narc(p0,small_wave,small_spec,err_spec)
            print ('Fit=',fit_params,'\n', 'MC error=',err_params)
            plot_Narc(small_wave,small_spec,err_spec,fit_params)
        
            (amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,m,c) = fit_params
        
            if (amp_Ha > 0) and (amp_Ha < 1e3):  
                ampHa_map[j,i] = amp_Ha                    
                    
                if (vel > -300) and (vel < 300):
                    veln_map[j,i] = vel
                else:
                    veln_map[j,i] = 0
            
                if (vel_sigma > -500) and (vel_sigma < 500):
                    velsign_map[j,i] = vel_sigma
                else:
                    velsign_map[j,i] = 0
                
                if (amp_NII6585 > 0) and (amp_NII6585 < 300):
                    ampNII_map[j,i] = amp_NII6585
                else:
                    ampNII_map[j,i] = 0
                
                if (amp_SII6716 > 0) and (amp_SII6716 < 300) and (amp_SII6730 > 0) and (amp_SII6730 < 300):
                    ampSII_map[j,i] = amp_SII6716 + amp_SII6730
                else:
                    ampSII_map[j,i] = 0 
        
            else:
                ampHa_map[j,i] = 0
                
    
    for i in range(0,20):        
        for j in range(20,30):        
        
            print ('Spaxel= [',j,',',i,']')
            y = mini_datacube[:,i,j]                       
            y[np.isnan(y)]=1e-15
            small_wave, small_spec = data_init(z2,x,y)    
    
            fit_params, err_params = ltsq_mc_fitting_Narc(p0,small_wave,small_spec,err_spec)    #diff fit func
            print ('Fit=',fit_params,'\n', 'MC error=',err_params)
            plot_Narc(small_wave,small_spec,err_spec,fit_params)
        
            (amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,m,c) = fit_params
        
            if (amp_Ha > 0) and (amp_Ha < 1e3):  
                ampHa_map[j,i] = amp_Ha                    
                    
                if (vel > -300) and (vel < 300):
                    veln_map[j,i] = vel
                else:
                    veln_map[j,i] = 0
            
                if (vel_sigma > -500) and (vel_sigma < 500):
                    velsign_map[j,i] = vel_sigma
                else:
                    velsign_map[j,i] = 0

                if (amp_NII6585 > 0) and (amp_NII6585 < 300):
                    ampNII_map[j,i] = amp_NII6585
                else:
                    ampNII_map[j,i] = 0

                if (amp_SII6716 > 0) and (amp_SII6716 < 300) and (amp_SII6730 > 0) and (amp_SII6730 < 300):
                    ampSII_map[j,i] = amp_SII6716 + amp_SII6730
                else:
                    ampSII_map[j,i] = 0 
        
            else:
                ampHa_map[j,i] = 0
                
    for i in range(20,30):        
        for j in range(0,20):        

            print ('Spaxel= [',j,',',i,']')
            y = mini_datacube[:,i,j]                       
            y[np.isnan(y)]=1e-15
            small_wave, small_spec = data_init(z2,x,y)    

            fit_params, err_params = ltsq_mc_fitting_Narc(p0,small_wave,small_spec,err_spec)    #diff fit func
            print ('Fit=',fit_params,'\n', 'MC error=',err_params)
            plot_Narc(small_wave,small_spec,err_spec,fit_params)

            (amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,m,c) = fit_params

            if (amp_Ha > 0) and (amp_Ha < 1e3):  
                ampHa_map[j,i] = amp_Ha                    
                    #convert to flux

                if (vel > -300) and (vel < 300):
                    veln_map[j,i] = vel
                else:
                    veln_map[j,i] = 0

                if (vel_sigma > -500) and (vel_sigma < 500):
                    velsign_map[j,i] = vel_sigma
                else:
                    velsign_map[j,i] = 0


                if (amp_NII6585 > 0) and (amp_NII6585 < 300):
                    ampNII_map[j,i] = amp_NII6585
                else:
                    ampNII_map[j,i] = 0

                if (amp_SII6716 > 0) and (amp_SII6716 < 300) and (amp_SII6730 > 0) and (amp_SII6730 < 300):
                    ampSII_map[j,i] = amp_SII6716 + amp_SII6730
                else:
                    ampSII_map[j,i] = 0 

            else:
                ampHa_map[j,i] = 0        


    return ampHa_map, veln_map, velsign_map, ampHaB_map, velb_map, velsigb_map, ampNII_map, ampSII_map



# In[ ]:





# In[22]:


redshift1 = 1.40915     #spectroscopic redshift for CENTRAL region 
redshift2 = 1.40691     #spectroscopic redshift for NORTHERN ARC
    
datacube, errorcube, spec, waverange = input_cubes()
minidatacube, errspec = mini_cube(redshift1, datacube, errorcube, waverange)

ampHa, veln, velsign, ampHaB, velb, velsigb, ampNII, ampSII = fitting(redshift1, redshift2, waverange, minidatacube, errspec)


# In[ ]:





# In[34]:


fig, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(1,6)
ax1.imshow(ampHa, origin='lower', vmax=10)
ax2.imshow(veln, origin='lower')
ax3.imshow(velsign, origin='lower')
ax4.imshow(ampHaB, origin='lower', vmax=20)
ax5.imshow(ampNII, origin='lower', vmax=20)
ax6.imshow(ampSII, origin='lower', vmax=20)


# In[ ]:


hdu = fits.PrimaryHDU(ampHa)
hdu.header = cube_header
hdu.writeto('ampHa_map.fits',output_verify='fix')


# In[ ]:


hdu = fits.PrimaryHDU(veln)
hdu.header = cube_header
hdu.writeto('veln_map.fits',output_verify='fix')


# In[ ]:


hdu = fits.PrimaryHDU(velsign)
hdu.header = cube_header
hdu.writeto('velsign_map.fits',output_verify='fix')


# In[ ]:


hdu = fits.PrimaryHDU(ampHaB)
hdu.header = cube_header  
hdu.writeto('ampHaB_map.fits',output_verify='fix')


# In[ ]:


hdu = fits.PrimaryHDU(velb)
hdu.header = cube_header  
hdu.writeto('velb_map.fits',output_verify='fix')


# In[ ]:


hdu = fits.PrimaryHDU(velsigb)
hdu.header = cube_header  
hdu.writeto('velsigb_map.fits',output_verify='fix')


# In[ ]:


hdu = fits.PrimaryHDU(ampNII)
hdu.header = cube_header  
hdu.writeto('ampNII_map.fits',output_verify='fix')


# In[ ]:


hdu = fits.PrimaryHDU(ampSII)
hdu.header = cube_header  
hdu.writeto('ampSII_map.fits',output_verify='fix')


# In[ ]:




