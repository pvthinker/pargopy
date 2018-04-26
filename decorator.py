#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 07:16:12 2018

@author: therry
"""

import functools
import types
import time
 
class exec_time(object):
    """
    Decorator used to know the number of call for a function
    and the Mean time for its execution
    """
 
    def __init__(self, fonc):
        """
        Initialisation of the decorator
        """
        self.fonc = fonc
        self.c = 0 # counter for the number of executions
        self.tt = 0  # Total time for executios
        self.func = functools.wraps(fonc)(self) # Keep the name of the decorated function
 
    def __get__(self, inst, owner=None):
        """
        Necessary to decorate the functions : having the good 'self'
        """
        return types.MethodType(self, inst)
 
    def __call__(self, *args, **kwargs):
        """
        Function called for each decorated function call
        """
        t = time.clock()
        result = self.fonc(*args, **kwargs)
        t = time.clock()-t
        self.tt += t  # cumul des temps
        self.c += 1  # incrémentation du compteur
        #  print u"Mean Time : %.7f" % (self.tt/self.c)
        return result
 
    def results(self):
        """
        Return and print the result :
        - Number of call
        - Mean time for execution
        """        
        if self.c == 0: 
            tm = 0
        else: 
            tm = self.tt/self.c
        print('Number of call for for the function %s : %i' % (self.func, self.c))
        print('Mean time for the execution of the function %s : %f' % (self.func, tm))
