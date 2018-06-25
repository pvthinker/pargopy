#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 24 09:39:58 2018

@author: therry

Module contenant les fonctions utilisées pour distribuer les tâches entre les
différents esclaves

"""

from mpi4py import MPI
import numpy as np
import os.path as path

import database as db
import param as param
import tile as ti

comm = MPI.COMM_WORLD

myrank = comm.Get_rank()
nslaves = comm.Get_size()-1

# list of buffers used by master to receive msg from slaves
answer = [np.zeros((1,), dtype=int) for k in range(nslaves)]

# list of irecv managed by master
reqr = []



def ordering_tasks(tasks):
    """
    :param tasks: List containing the tasks to do
    
    Sort the tasks according to their workload
    workload is proportional to size of the tile file
    
    :rtype: list of int
    """
    def tilefilename(itile):
        return '%s/argo_%003i.pkl' % (param.get_path('parallel'), itile)

    workload = [path.getsize(tilefilename(t)) for t in tasks]
    idx = np.argsort(workload)

    return tasks[np.asarray(idx,dtype=int)]


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
    status = MPI.Status()
    islave = 0
    while (islave < nslaves) and (slavestate[islave] == 0):
        islave += 1

    if islave == nslaves:
        # all slaves are busy, let's wait for the first
        # print('waiting ...', len(reqr), answer)
        MPI.Request.Waitany(reqr, status)
        # print('ok a slave is available', answer)
        islave = 0
        while (islave < nslaves) and (answer[islave][0] == 0):
            islave += 1
        slavestate[islave] = 1  # available again
        # note  slave send the message int(1)
        # master is getting the msh in the answer array
        # but the received msg is never int(1) !!!
        answer[islave][0] = 0
        # todo: remove the reqr that is done

    else:
        pass
    # print('=> %i is available' % islave)
    return islave


def master_work_nonblocking(nslaves):
    """
    :param nslaves: Number of slavbes under master control
    
    Master basically supervises things but does no work
    
    :rtype: None
    """
    fid = open('master.txt','w')
    tasks = range(300)
    tasks = np.arange(30,300)#np.arange(28)
    #  tasks = [52, 0, 19, 280, 299, 97, 125, 166, 153, 199, 142, 16, 53, 129]
    nbtasks = len(tasks)
    # sorting the tasks according to their size
    # improves the load balance among slaves (by a lot)
    # for instance, it prevents cases where the last
    # task is a very long one, which would ruin their
    # global performance

    #tasks = ordering_tasks(tasks)
    
    print('List of tasks to be done:', tasks)
    fid.write('after sorting')
    fid.flush()
    itask = 0
    slavestate = [1 for k in range(nslaves)]

    record = np.zeros((nbtasks,), dtype=int)
    synchronized = {}
    for t in tasks:
        synchronized[t] = False

    while itask < nbtasks:
        islave = getavailableslave(slavestate)
        record[itask] = islave+1
        print('master: send to %i' % (islave+1))
        fid.write('islave = %i' % islave)
        fid.flush()
        comm.isend((itask, tasks[itask]), dest=islave+1, tag=islave+1)
        # print('answer[%i] = ' % islave, answer[islave])
        reqr.append(comm.Irecv(answer[islave], source=islave+1, tag=islave+1))        
        slavestate[islave] = 0
        itask += 1
        previous = np.where(record==(islave+1))[0]
        if len(previous)>1:
            oldtask = tasks[previous[-2]]
            #db.synchronize_argo_global_from_argo_tile(itiles=oldtask)
            synchronized[oldtask] = True

    for t in tasks:
        if synchronized[t]:
            pass
        else:
            #db.synchronize_argo_global_from_argo_tile(itiles=t)
            synchronized[t] = True

    # all tasks have been done
    # tell the slaves to stop
    for islave in range(1, nslaves+1):
        comm.isend('done', dest=islave, tag=islave)

    comm.Barrier()
    print('-'*40)
    print('  Summary of who did what')
    for k in range(itask):
        print('    - task %2i / %3i data / done by %2i'
              % (k, tasks[k], record[k]))
    print('-'*40)
    print(' Load balance')
    for i in range(1, nslaves+1):
        nb = sum([tasks[k] for k in range(nbtasks) if record[k] == i])
        print('    - slave %2i treated %3i data' % (i, nb))
    print('-'*40)


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
    completed_tasks = []
    ok = False
    while not(ok):
        # print('slave %i waiting' % islave)
        msg = comm.recv(source=0, tag=islave)
        if type(msg) is str:
            if msg == 'done':
                ok = True
                print('#%2i* stops' % islave)

        if type(msg) == tuple:
            itask, itile = msg
            print('#%2i* treating task %03i / %03i data'
                  % (islave, itask, itile))

            # do the work, replace 'sleep' with a real work!
            ti.main(itile)

            # tell the master that we are done and that he
            # can send another task
            comm.isend(int(1), dest=0, tag=islave)

            # keep track of what have been done
            completed_tasks.append((itask, itile))

    comm.Barrier()
    for task in completed_tasks:
        print('#%2i- did task %03i / %3i data'
              % (islave, task[0], task[1]))
    print('-'*40)



if __name__ == '__main__':

    # Update argo_global with newest version of Argo
    #  db.update_argo_global()
    
    if myrank == 0:
        print('Hello I\'m the master, I\'ve %i slaves under my control'
              % nslaves)
        #comm.Barrier()
        master_work_nonblocking(nslaves)

    else:
        print('#%2i* starts' % myrank)
        #comm.Barrier()

        slave_work_nonblocking(myrank)
    
    
    # Synchronize argo_global with argo_tiles after interpolation
    #  db.synchronize_argo_global_from_argo_tile()
