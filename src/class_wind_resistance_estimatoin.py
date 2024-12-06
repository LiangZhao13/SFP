# -*- coding: utf-8 -*-
"""
Created on Mon May 20 16:10:59 2024

@author: LLei_seu
"""

'''
风阻计算
基本输入：
船首向 heading_ship 度为单位
船舶航速 对地航速 m/s
船舶参数：获取水面线以上横截面积 平方米m2
水平方向上风速 knot
垂直方向上风速 knot


待明确的问题：
(1) 水平方向和垂直方向
(2) 角度计算的过程
'''

import numpy as np
import attr

def wind_resistance(speed_ship,heading_ship,wu,
                             wv,ship,rho_air = 1.225):
    '''
    估算风阻 ：单位N
    rho_air： 空气质量密度，kg/m3
    C_AA: 风阻系数,通过拟合来获取
    wu: 水平方向上风速 统一直接输入m/s
    wv: 垂直方向上风速 统一直接输入m/s
    '''
    #相对风向
    wind_U = wu  #转换为m/s * 0.5144
    wind_V = wv  #转换为m/s * 0.5144 
    wind_hdg = 90 - np.rad2deg(np.arctan2(wind_U, wind_V)) #不太明白？？？
    if wind_hdg < 0:
        wind_hdg += 360
    
    #相对速度
    ship_U = speed_ship * np.cos(np.deg2rad(heading_ship))
    ship_V = speed_ship * np.sin(np.deg2rad(heading_ship))

    appwind_U = wind_U + ship_U
    appwind_V = wind_V + ship_V
    appwind_spd = np.sqrt(pow(appwind_U, 2) + pow(appwind_V, 2))
    
    #相对风向
    appwind_hdg = relative_angle_360(wind_hdg, heading_ship)
    
    
    r_wind = 0.5 * rho_air * C_AA (appwind_hdg) * ship.A_transverse * pow(appwind_spd, 2) -\
        0.5 * rho_air * C_AA (0) *  ship.A_transverse * pow(speed_ship, 2)
        
    return r_wind/1000

def relative_angle_360(alpha, beta):
    '''
    船舶和对象相对方向
    '''
    return np.abs(180 - np.abs(np.abs(alpha - beta) - 180)) #不太明白

def C_AA(heading_ship2wind):
    '''
    计算风阻系数
    heading_ship2wind：船舶相对风的方向，单位为度
    (0)表示船舶正对风头
    '''
    c_aa = [0.55, 0.90, 1.0, 1.0, 0.9, 0.87, 0.62, 0.45, 0.25, 0.1, -0.1, -0.48, -0.85, -1.0, -1.42, -1.49, -1.38, -0.9, -0.85]
    heading = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180]
    a = np.interp(heading_ship2wind, heading, c_aa)

    return a

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
    
if __name__=='__main__':
    #输入的wu 和wv分别是以正东和正北的方向的风速
    wr=wind_resistance(30,60,7,9,ship = ShipParam())
    print(wr)
