#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from numpy import asarray as ar,exp

def data_init(red,x,y):
    k = 1+red
    wave = (x/k)*10000 
    select = (wave>6450) & (wave<6800)
    small_wave0 = wave[select]
    small_spec0 = y[select]
    return small_wave0, small_spec0

def error_scale(data_spec,err_spec_temp):
    flux_cont = np.std(data_spec)
    flux_err = np.mean(err_spec_temp)
    b = flux_cont/flux_err
    err_spectrum = b*err_spec_temp
    return b,err_spectrum

def redshift(vel):
    return vel/300000.0 

def line_width(vel_sigma,rest_line,inst_res_fwhm=0.0): 
    sigma = vel_sigma/(300000.0-vel_sigma)*rest_line
    return np.sqrt(sigma**2+(inst_res_fwhm/2.354)**2) 

def gauss(wave_range,amplitude,vel,vel_sigma,rest_w):
    l = (amplitude)*exp(-(wave_range-(rest_w*(1+redshift(vel))))**2/(2*(line_width(vel_sigma, rest_w))**2))
    return l

def full_gauss(wave_range,amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,amp_HaB,vel_HaB,vel_sigma_HaB,m,c):
        Ha = gauss(wave_range,amp_Ha,vel,220,6562.8)
        Ha_broad = gauss(wave_range,amp_HaB,vel_HaB,vel_sigma_HaB, 6562.8)
        NII_6583 = gauss(wave_range,amp_NII6585,vel,220,6583.46)
        NII_6548 = (0.34)*gauss(wave_range,amp_NII6585,vel,220,6548.05) 
        SII_6716 = gauss(wave_range,amp_SII6716,vel,220,6716.4)
        SII_6730 = gauss(wave_range,amp_SII6730,vel,220,6730.8)
        cont = (wave_range/1000.0)*m+c
        return Ha + Ha_broad + NII_6548 + NII_6583 + SII_6716 + SII_6730 + cont
        
def full_gauss_residual(params,wave_range,data,error):
    (amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,amp_HaB,vel_HaB,vel_sigma_HaB,m,c) = params
    
    SII_6716 = gauss(wave_range,amp_SII6716,vel,vel_sigma,6716.4)
    SII_6730 = gauss(wave_range,amp_SII6730,vel,vel_sigma,6730.8)
    Ha = gauss(wave_range,amp_Ha,vel,vel_sigma,6562.8)
    Ha_broad = gauss(wave_range,amp_HaB,vel_HaB,vel_sigma_HaB, 6562.8)
    NII_6583 = gauss(wave_range,amp_NII6585,vel,vel_sigma,6583.46)
    NII_6548 = (0.34)*gauss(wave_range,amp_NII6585,vel,vel_sigma,6548.05)    
    cont = (wave_range/1000.0)*m+c
    
    return ((Ha+Ha_broad+NII_6548+NII_6583+SII_6716+SII_6730+cont)-data)/error

def full_gauss_noblr(wave_range,amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,m,c):
        Ha = gauss(wave_range,amp_Ha,vel,220,6562.8)
        NII_6583 = gauss(wave_range,amp_NII6585,vel,220,6583.46)
        NII_6548 = (0.34)*gauss(wave_range,amp_NII6585,vel,220,6548.05) 
        SII_6716 = gauss(wave_range,amp_SII6716,vel,220,6716.4)
        SII_6730 = gauss(wave_range,amp_SII6730,vel,220,6730.8)
        cont = (wave_range/1000.0)*m+c
        return Ha + NII_6548 + NII_6583 + SII_6716 + SII_6730 + cont
        
def full_gauss_residual_noblr(params,wave_range,data,error):
    (amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,m,c) = params
    
    SII_6716 = gauss(wave_range,amp_SII6716,vel,vel_sigma,6716.4)
    SII_6730 = gauss(wave_range,amp_SII6730,vel,vel_sigma,6730.8)
    Ha = gauss(wave_range,amp_Ha,vel,vel_sigma,6562.8)
    NII_6583 = gauss(wave_range,amp_NII6585,vel,vel_sigma,6583.46)
    NII_6548 = (0.34)*gauss(wave_range,amp_NII6585,vel,vel_sigma,6548.05)    
    cont = (wave_range/1000.0)*m+c
    
    return ((Ha+NII_6548+NII_6583+SII_6716+SII_6730+cont)-data)/error

def ltsq_mc_fitting(p0,wave_range,data,error):
    # amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,amp_HaB,vel_HaB,vel_sigma_HaB, m, c
    
    pfit, pcov = leastsq(full_gauss_residual, x0=p0, args=(wave_range, data, error), maxfev=10000000)
    
    MonteCarlo_loops = 100
    param_MC = np.zeros((len(pfit), MonteCarlo_loops)) 
	
    for i in range(MonteCarlo_loops):
        iteration_data = np.random.normal(data,error)     
        pfit_MC, pcov_MC = leastsq(full_gauss_residual, x0=p0, args=(wave_range, iteration_data, error), maxfev = 10000000) 
        param_MC[:,i] = pfit_MC
    param_err = np.std(param_MC, axis=1)
    #(amp_Ha_err,amp_NII6585_err,amp_SII6716_err,amp_SII6730_err,vel_err,
    #   vel_sigma_err,amp_HaB_err,vel_HaB_err,vel_sigma_HaB_err,m_err,c_err) = param_err
    return pfit, param_err

def ltsq_mc_fitting_Narc(p0,wave_range,data,error):
    # amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,m, c
    
    pfit, pcov = leastsq(full_gauss_residual_noblr, x0=p0, args=(wave_range, data, error), maxfev=10000000)
    
    MonteCarlo_loops = 100
    param_MC = np.zeros((len(pfit), MonteCarlo_loops)) 

    for i in range(MonteCarlo_loops):
        iteration_data = np.random.normal(data, np.abs(error))     #to avoid negative scale error
        pfit_MC, pcov_MC = leastsq(full_gauss_residual_noblr, x0=p0, args=(wave_range, iteration_data, error), maxfev = 10000000) 
        param_MC[:,i] = pfit_MC
    param_err = np.std(param_MC, axis=1)
    #(amp_Ha_err,amp_NII6585_err,amp_SII6716_err,amp_SII6730_err,vel_err,
    #   vel_sigma_err,m_err,c_err) = param_err
    return pfit, param_err

def plot(small_wave,small_spec,err_spec,fit_params):
#def plot(i,j,small_wave,small_spec,err_spec,fit_params):
    (amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,amp_HaB,vel_HaB,vel_sigma_HaB,m,c) = fit_params
    
    yfit = full_gauss(small_wave,amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,amp_HaB,vel_HaB,vel_sigma_HaB,m,c)
    
    fig = plt.figure()
    ax1 = fig.add_axes([0.1, 0.27, 0.84, 0.7]) # main axes
    ax2 = fig.add_axes([0.1, 0.06, 0.84, 0.2]) # inset axes

    ax1.plot(small_wave, small_spec, 'k-', label='Spectrum', linewidth=0.7)
    ax1.plot(small_wave, yfit,'r-',label='Fit', linewidth=0.7)
    ax1.plot(small_wave, gauss(small_wave,amp_Ha,vel,vel_sigma,6562.8), 'g', linestyle='--', linewidth=0.7, label='Narrow')
    ax1.plot(small_wave, gauss(small_wave,amp_HaB,vel_HaB,vel_sigma_HaB, 6562.8), 'b', linestyle='--', linewidth=0.7, label='Broad')
    ax1.plot(small_wave, gauss(small_wave,amp_NII6585,vel,vel_sigma,6583.46), 'g', linestyle='--', linewidth=0.7)
    ax1.plot(small_wave, (0.34)*gauss(small_wave,amp_NII6585,vel,vel_sigma,6548.05), 'g', linestyle='--', linewidth=0.7)
    ax1.plot(small_wave, gauss(small_wave,amp_SII6716,vel,vel_sigma,6716.4), 'g', linestyle='--', linewidth=0.7)
    ax1.plot(small_wave, gauss(small_wave,amp_SII6730,vel,vel_sigma,6730.8), 'g', linestyle='--', linewidth=0.7)
    
    residuals = (small_spec - yfit)/err_spec
    
    ax2.plot(small_wave, residuals, 'gray', label='Residuals', linewidth=0.7)
    plt.rcParams["figure.figsize"] = [8,6]
    #ax1.set_title('3C297 central'+'('+ i +','+ j + ')'+'spaxel', {'fontsize': 18})
    ax1.set_ylabel('Counts', {'fontsize': 14})
    ax2.set_xlabel('Wavelength (Angstroms)', {'fontsize': 14})
    ax1.legend()
    ax2.legend()
    plt.show()

def plot_Narc(small_wave,small_spec,err_spec,fit_params):

    (amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,m,c) = fit_params
    
    yfit = full_gauss_noblr(small_wave,amp_Ha,amp_NII6585,amp_SII6716,amp_SII6730,vel,vel_sigma,m,c)
    
    fig = plt.figure()
    ax1 = fig.add_axes([0.1, 0.27, 0.84, 0.7]) # main axes
    ax2 = fig.add_axes([0.1, 0.06, 0.84, 0.2]) # inset axes

    ax1.plot(small_wave, small_spec, 'k-', label='Spectrum',linewidth=1.0)
    ax1.plot(small_wave, yfit,'r-',label='Fit')
    ax1.plot(small_wave, gauss(small_wave,amp_Ha,vel,vel_sigma,6562.8), 'g', linestyle='--', linewidth=0.7, label='Narrow')
    ax1.plot(small_wave, gauss(small_wave,amp_NII6585,vel,vel_sigma,6583.46), 'g', linestyle='--', linewidth=0.7)
    ax1.plot(small_wave, (0.34)*gauss(small_wave,amp_NII6585,vel,vel_sigma,6548.05), 'g', linestyle='--', linewidth=0.7)
    ax1.plot(small_wave, gauss(small_wave,amp_SII6716,vel,vel_sigma,6716.4), 'g', linestyle='--', linewidth=0.7)
    ax1.plot(small_wave, gauss(small_wave,amp_SII6730,vel,vel_sigma,6730.8), 'g', linestyle='--', linewidth=0.7)
    
    residuals = (small_spec - yfit)/err_spec
    
    ax2.plot(small_wave, residuals, 'gray', label='Residuals', linewidth=0.7)
    plt.rcParams["figure.figsize"] = [8,6]
    #ax1.set_title('3C297 central'+'('+ i +','+ j + ')'+'spaxel', {'fontsize': 18})
    ax1.set_ylabel('Counts', {'fontsize': 14})
    ax2.set_xlabel('Wavelength (Angstroms)', {'fontsize': 14})
    ax1.legend()
    ax2.legend()
    plt.show()
