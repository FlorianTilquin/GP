#! /usr/bin/env python
# -*- encoding: utf-8 -*-
# Rotative routine of Kepler's modules

#  Convert module and output channel number from one quarter to the next
def Kepler_rotate(modin, outin):
    outout = outin
    if modin == 2:
        modout = 16
    if modin == 3: 
        modout = 11
    if modin == 4: 
        modout = 13
    if modin == 6: 
        modout = 22
    if modin == 7: 
        modout = 17
    if modin == 8: 
        modout = 12
    if modin == 9: 
        modout = 7
    if modin == 10: 
        modout = 2
    if modin == 11:
        modout = 23
    if modin == 12:
        modout = 18
    if modin == 13:
        modout = 13 
        outout = outin + 1
        if outout == 5:
           outout = 1
    if modin == 14:
        modout = 8
    if modin == 15:
        modout = 3
    if modin == 16:
        modout = 24
    if modin == 17:
        modout = 19
    if modin == 18:
        modout = 14
    if modin == 19:
        modout = 9
    if modin == 20:
        modout = 4
    if modin == 21:
        modout = 1
    if modin == 22:
        modout = 20
    if modin == 23:
        modout = 15
    if modin == 24:
        modout = 10
    return modout, outout

