# -*- coding: utf-8 -*-
"""
Created on Mon May 20 21:39:38 2024

@author: LLei_seu
"""


'''
船舶浪阻的计算

（1）待明确heading_wave是什么
（2）坐标系的方向
'''

import numpy as np
import attr

@attr.s(frozen=True, auto_attribs=True)
class ShipParam(object):
    Lwl: float = 235
    Lpp: float = 230.4
    B: float = 41.5
    Cp: float = 0.8467
    Cm: float = 0.99
    C_B: float = 0.8219
    C_wp: float = 0.91
    L_CB: float = -0.75
    T: float = 13.1
    Tf: float = 13.1
    Ss: float = 13682
    ks: float = 150 * 10 ** - 6
    Vol: float = 105000
    g: float = 9.81
    A_bt: float = 30.8
    hb: float = 4
    C_stern: float = 10
    k2: float = 1.5
    transom_a: float = 42.6
    i_E: float = 5.0

    Sapp: float = 50
    I: float = 14762925
    Awl: float = 9561.6
    kyy: float = 60.2966
    A_transverse: float = 41*5 #如何取值
    
def jonswap(w, Hs, Tp, gamma):
    omega = w
    wp = 2 * np.pi / Tp
    sigma = np.where(omega < wp, 0.07, 0.09)
    a = np.exp(-0.5 * np.power((omega - wp) / (sigma * wp), 2.0))
    sj = 320 * np.power(Hs, 2) * np.power(omega, -5.0) / np.power(Tp, 4) * \
          np.exp(-1950 * np.power(omega, -4) / np.power(Tp, 4)) * np.power(gamma, a)
    return sj

def froude_number(v, ship = ShipParam()):
    '''
    Function to calculate froude number
    '''
    fn_v = v / np.sqrt(ship.g * ship.Lpp)
    return fn_v

def amplitude_function(beta, Hs):
    '''
    根据夹角得到的浪阻增加函数
    '''

    if Hs > 3:
        v = pow(Hs / 3, 1.5)
        beta = beta / v
        f = 0.9667 * np.sin(0.009004 * beta + 1.73) + 0.3348 * np.sin(0.02282 * beta + 3.014)
    else:
        f = 0.9667 * np.sin(0.009004 * beta + 1.73) + 0.3348 * np.sin(0.02282 * beta + 3.014)
    return f

def angular_distribution_function(alpha):

    rad_alpha = np.deg2rad(alpha)
    return 2 / np.pi * np.cos(rad_alpha) * np.cos(rad_alpha)

def reflection_wave_resistance_coefficient(v, w, beta, Hs, ship):
    """
    Funtion to calculate Added resistance in regular waves due to diffraction effect (reflection effect)
    v: ship velocity [m/s]
    z: incident regular wave amplitude [m]
    w: wave frequency range [rad/s]
    """
    # e = np.arctan(1 / np.tan(np.deg2rad(ship.i_E)))
    e = np.deg2rad(ship.i_E)
    B_f=2.25 * np.sin(e) * np.sin(e)
    
    fn = froude_number(v)
    expm = pow((0.87 / ship.C_B), (1 + 4 * np.sqrt(fn)))

    w_l = ship.g * 2 * np.pi / w / w
    simp = 1 + 5 * np.sqrt(ship.Lpp / w_l) * fn
    
    k = 2 * np.pi / w_l
    alpha_t = 1 - np.exp(-2 * k * ship.T)
    
    c_awr = 0.5 * ship.Lpp / ship.B * B_f * simp * expm * alpha_t * amplitude_function(beta, Hs)

    return c_awr

def switch_function(beta):
    '''为什么要进行转换'''
    rad_beta = np.deg2rad(beta)
    f = np.cos(rad_beta / 2)
    
    return f

def motion_wave_resistance_coefficient(v, w, beta, Hs, ship):
    """
    Function to calcualte added resistance in regular waves due to motion effect

    v: ship velocity [m/s]
    z: incident regular wave amplitude [m]
    w: wave frequency range [rad/s]
    """
    switch_value = switch_function(beta)    
    boolC_B = ship.C_B < 0.75
    fn = froude_number(v)
    a1 = 60.3 * pow(ship.C_B, 1.34) * pow((0.87 / ship.C_B), (1 + fn))
    a2 = 0.0072 + 0.1676 * fn if fn < 0.12 else pow(fn, 1.5) * np.exp(- 3.5 * fn)
    p1 = np.sqrt(ship.Lpp / ship.g) * pow((ship.kyy / ship.Lpp), (1. / 3.))
    omeg = p1 * pow(0.05, 0.143) / 1.17 if fn < 0.05 else p1 * pow(fn, 0.143) / 1.17
    omeg = switch_value * omeg

    if boolC_B:
        omega = omeg * w
        b1 = 11 if omega < 1 else -8.5
        d1 = 14 if omega < 1 else -566 * np.power((ship.Lpp / ship.B), (-2.66)) * 6
        c_awm = 4 * pow(omega, b1) * np.exp(b1 / d1 * (1 - pow(omega, d1))) * a1 * a2 * amplitude_function(beta, Hs)
    else:
        omega = omeg * w
        b1 = 11 if omega < 1 else -8.5
        d1 = 566 * np.power((ship.Lpp / ship.B), (-2.66)) if omega < 1 else -566 * np.power((ship.Lpp / ship.B), (-2.66)) * 6
        c_awm = 4 * pow(omega, b1) * np.exp(b1 / d1 * (1 - pow(omega, d1))) * a1 * a2 * amplitude_function(beta, Hs)

    return c_awm

def wave_resistance(v, Hs, Tp, heading_wave, heading_ship, ship = ShipParam(), rho_sw = 1025.):
    """
    Funtion to calculate added resistance in irregular waves
    """
    num_of_angle = 13
    dangular = np.deg2rad(180 / (num_of_angle - 1))
    alpha = np.linspace(-90, 90, 13)
    wave_angle = np.mod(heading_wave + alpha, 360)
    # initialize angular frequency
    dw = 0.05
    w = np.arange(0.35, 2.5, dw) #取值范围

    # initialize Jonswap spectra and wave elevation
    S_wave = jonswap(w, Hs, Tp, gamma = 3.7)
    z = np.sqrt(2 * S_wave * dw)

    heading_ship2wave = np.array([relative_angle_360(wave_angle[i], heading_ship) for i in range(num_of_angle)])
    dimension_param = rho_sw * ship.g * ship.B * ship.B / ship.Lpp #为了统一

    r_wave_resist = 0
    for i in range(num_of_angle):
        r_awr = 0
        r_awm = 0
        r_aw = 0
        r_wave = 0
        for j in range(len(w)):
            r_awr = dimension_param * reflection_wave_resistance_coefficient(v, w[j], heading_ship2wave[i], Hs, ship)
            r_awm = dimension_param * motion_wave_resistance_coefficient(v, w[j], heading_ship2wave[i], Hs, ship)
            r_aw = r_awr + r_awm
            r_wave +=S_wave[j] * r_aw * dw * angular_distribution_function(alpha[i])
        r_wave_resist += 2 * r_wave * dangular

    return r_wave_resist/1000

def relative_angle_360(alpha, beta):
    '''
    船舶和对象相对方向
    '''
    return np.abs(180 - np.abs(np.abs(alpha - beta) - 180)) #不太明白



if __name__=='__main__':
    
    # print()
    # S_wave = jonswap(w, Hs, Tp, gamma = 3.7)
    Hs=4 #有效浪高m
    Tp=8 #波峰周期s
    heading_ship=60 #船首向
    heading_wave=40 #浪的方向
    r_wave_resist=wave_resistance(30*0.5144,Hs,Tp,heading_wave,heading_ship)
    print(r_wave_resist)