#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:39:58 2018

@author: therry

Module contenant les fonctions utilisées pour distribuer les tâches entre les
différents esclaves

"""

def ordering_tasks(tasks):
    """
    :param tasks: List containing the tasks to do
    
    Sort the tasks according to their workload
    workload is proportional to size of the tile file
    
    :rtype: list of int
    """
    pass

def getavailableslave(slavestate):
    """
    :param slavestate: Give the state of the slave
    
    Return the index of a slave that is awaiting a task.  A busy slave
    has a state == 0. If all slaves are busy then wait until a msg is
    received, the msg is sent upon task completion by a slave. Then
    determin who sent the msg. The msg is collected in the answer
    array. By scanning it, we determine who sent the message.

    :rtype: int
    """
    pass

def master_work_nonblocking(nslaves):
    """
    :param nslaves: Number of slavbes under master control
    
    Master basically supervises things but does no work
    
    :rtype: None
    """
    pass

def slave_work_nonblocking(islave):
    """
    :param islave: Number given to the slave
    
    Slaves enter an infinite loop: keep receiving messages from the
    master until reception of 'done'. Each messages describes the task
    to be done. When a new task is received slave treats it. At the
    end of it the slave sends a message to the master saying that he
    is over, and that he is available for a new task.
    
    :rtype: None
    """
    pass