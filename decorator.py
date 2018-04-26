#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 07:16:12 2018

@author: therry
"""

import functools
import types
import time
import ast
 
def call_results(filename):
    """
    Function used to know the number of call and the mean time of execution
    of each function from a chosen file
    Do not forget to import the file you want to analyse 
    before using this function !
    """
    if filename == 'interpolation_tools.py':
        name_module = 'interpolation'
    elif filename == 'research_tools.py':
        name_module = 'research'
    elif filename == 'argotools.py':
        name_module = 'argotools'
    elif filename == 'stats.py':
        name_module = 'stats'
    elif filename == 'argodb.py':
        name_module = 'argodb'
    elif filename == 'tile.py':
        name_module = 'tiler'
    else:
        raise ValueError('Filename not recognize')

    def top_level_functions(body):
        return (f for f in body if isinstance(f, ast.FunctionDef))

    def parse_ast(filename):
        with open(filename, "rt") as file:
            return ast.parse(file.read(), filename=filename)
    #  Function used to find the list of the functions contained in a python file
    tree = parse_ast(filename)
    for func in top_level_functions(tree.body):
        c = name_module+'.'+func.name+'.results()'
        print('----- Function : %s -----' % func.name)
        print eval(c)


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
        Print the result :
        - Number of call
        - Mean time for execution
        """        
        if self.c == 0: 
            tm = 0
        else: 
            tm = self.tt/self.c
        print('Number of call for for the function : %i' % self.c)
        print('Mean time for the execution of the function : %f' % tm)