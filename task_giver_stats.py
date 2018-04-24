"""Example of master & slaves program:

- a very simple one with blocking communications, to get a first idea
  of how mpi works

- a more advanced one with non-blocking communications, where each
  slave works at its pace. Master is monitoring the collective
  work. As soon as one slave is done with a task, master sends it a
  new task, until all tasks have been done.

"""

from mpi4py import MPI
import numpy as np
import stats as stats
import os.path as path
import param as param

path_to_tiles = param.path_to_tiles

year = '2017'
month = '12'
day = '31'
date = [year, month, day]
# mode defines the values selected :
# R : Real time
# A : Adjusted Real Time
# D : Delayed time (Values verified)
mode = 'D'
typestat = 'zstd'
reso = 0.5
timeflag = 'annual'
#  nbtasks = 80  # master will defined nbtasks of random size

tmax = 15  # max size of each task

comm = MPI.COMM_WORLD

myrank = comm.Get_rank()
nslaves = comm.Get_size()-1

# list of buffers used by master to receive msg from slaves
answer = [np.zeros((1,), dtype=int) for k in range(nslaves)]

# list of irecv managed by master
reqr = []


def ordering_tasks(tasks):
    """Sort the tasks according to their workload
    workload is proportional to size of the tile file"""
    def tilefilename(itile):
        return '%s/tile%03i.pkl' % (path_to_tiles, itile)

    workload = [path.getsize(tilefilename(t)) for t in tasks]
    idx = np.argsort(workload)
    print(type(idx[0]))
    return tasks[int(idx)]


def getavailableslave(slavestate):
    """Return the index of a slave that is awaiting a task.  A busy slave
    has a state == 0. If all slaves are busy then wait until a msg is
    received, the msg is sent upon task completion by a slave. Then
    determin who sent the msg. The msg is collected in the answer
    array. By scanning it, we determine who sent the message.

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
    """Main program for master

    Master basically supervises things but does no work

    """

    tasks = range(276, 300)
    #  tasks = [285, 283, 284, 292, 297, 295, 294, 296, 293]
    nbtasks = len(tasks)
    # sorting the tasks according to their size
    # improves the load balance among slaves (by a lot)
    # for instance, it prevents cases where the last
    # task is a very long one, which would ruin their
    # global performance
    #  tasks = ordering_tasks(tasks)

    print('List of tasks to be done:', tasks)

    itask = 0
    slavestate = [1 for k in range(nslaves)]

    record = np.zeros((nbtasks,), dtype=int)

    while itask < nbtasks:
        islave = getavailableslave(slavestate)
        record[itask] = islave+1
        print('send to %i' % (islave+1))
        comm.isend((itask, tasks[itask]), dest=islave+1, tag=islave+1)
        # print('answer[%i] = ' % islave, answer[islave])
        reqr.append(comm.Irecv(answer[islave], source=islave+1, tag=islave+1))
        slavestate[islave] = 0
        itask += 1

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
    """Main function for slaves.

    Slaves enter an infinite loop: keep receiving messages from the
    master until reception of 'done'. Each messages describes the task
    to be done. When a new task is received slave treats it. At the
    end of it the slave sends a message to the master saying that he
    is over, and that he is available for a new task.

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
            stats.main(itile, typestat, reso, timeflag, date, mode)

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

    if myrank == 0:
        print('Hello I\'m the master, I\'ve %i slaves under my control'
              % nslaves)
        comm.Barrier()
        master_work_nonblocking(nslaves)

    else:
        print('#%2i* starts' % myrank)
        comm.Barrier()

        slave_work_nonblocking(myrank)
