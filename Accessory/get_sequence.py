#####################################################################
# Copyright(c) A Vipin Menon and BIG LAB in Hanyang University (HYU)# 
# Author A Vipin Menon						    #
# Date 18th September 2019					    #
# Email a.vipin.menon@gmail.com					    #
#####################################################################				
import numpy
import csv
import sys, os, math
import string

def reverseString(st):
        li = []
        for i in st: li.append(i)
        li.reverse()
        return ''.join(li)


def reverseComp(st):
        comp = str.maketrans('ATCG', 'TAGC')
        return reverseString(st).translate(comp)
