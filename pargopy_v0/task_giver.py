"""
Masternslave used for the tiles creation

"""

from mpi4py import MPI
import numpy as np
import time as time
import tile as tiler
import os.path as path
import param as param
#  nbtasks = 80  # master will defined nbtasks of random size

path_to_tiles = param.path_to_tiles
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
    workload is proportional to size of the tile file
    
    :rtype: list of int"""
    def tilefilename(itile):
        return '%s/tile%03i.pkl' % (path_to_tiles, itile)

    workload = [path.getsize(tilefilename(t)) for t in tasks]
    idx = np.argsort(workload)
    print(type(idx[0]))
    return tasks[int(idx)]


def master_work_blocking(nslaves):
    """Master organizes the work using blocking communications with
    slaves

    """

    for j in range(3):
        for islave in range(1, nslaves+1):
            nx = int(50*np.random.uniform())
            print('slave %i needs to treat %i data' % (islave, nx))
            comm.isend(nx, dest=islave, tag=islave)

            # this is blocking the loop
            answer = comm.recv(source=islave, tag=islave)
            print('slave %i is beging for %s' % (islave, answer))

    for islave in range(1, nslaves+1):
        comm.isend('done', dest=islave, tag=islave)


def slave_work_blocking(islave):
    """Very simple function for slaves based on blocking communications
    with master

    """
    ok = False
    while not(ok):
        msg = comm.recv(source=0, tag=islave)
        if type(msg) is str:
            if msg == 'done':
                ok = True
                print('ok slave %i stops' % islave)
        if type(msg) == int:
            nx = msg
            print('slave %i treating %i data' % (islave, nx))
            # do here some work, e.g. sleep for nx seconds
            time.sleep(nx)
            # tell master that you're done
            comm.send('more', dest=0, tag=islave)


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
        print('%i slave waiting ...' % islave, len(reqr), answer)
        MPI.Request.Waitany(reqr, status)
        print('ok a slave is available', answer)
        islave = 0
        while (islave < nslaves) and (answer[islave][0] == 0):
            islave += 1
        print('islave = %i' % islave)
        slavestate[islave] = 1  # available again
        # note  slave send the message int(1)
        # master is getting the msh in the answer array
        # but the received msg is never int(1) !!!
        answer[islave][0] = 0
        # todo: remove the reqr that is done

    else:
        print('=> %i is available' % (islave+1))
    return islave


def master_work_nonblocking(nslaves):
    """Main program for master

    Master basically supervises things but does no work

    """
    
    #  tasks = range(200, 300)
    #  tasks = [126, 131, 146, 149, 151, 171, 188, 190, 210, 233, 234, 235, 244, 
    #           253, 254, 255, 264, 265, 272, 273, 274, 275, 276, 284, 285, 295, 296]
    tasks = range(300)
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
        print('send work to %i' % (islave+1))
        comm.isend((itask, tasks[itask]), dest=islave+1, tag=islave+1)
        print('answer[%i] = ' % islave, answer[islave])
        reqr.append(comm.Irecv(answer[islave], source=islave+1, tag=islave+1))
        slavestate[islave] = 0
        itask += 1

    # all tasks have been done
    # tell the slaves to stop
    for islave in range(1, nslaves+1):
        comm.isend('done', dest=islave, tag=islave)

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
                print('ok slave %i stops' % islave)

        if type(msg) == tuple:
            itask, itile = msg
            print('slave %2i treating task %2i / %3i data'
                  % (islave, itask, itile))

            # do the work, replace 'sleep' with a real work!
            #  time.sleep(nx)
            tiler.main(itile)

            # tell the master that we are done and that he
            # can send another task
            comm.isend(int(1), dest=0, tag=islave)

            # keep track of what have been done
            completed_tasks.append((itask, itile))

    for task in completed_tasks:
        print('slave %2i did taks %3i / %3i data'
              % (islave, task[0], task[1]))
    print('-'*40)


if __name__ == '__main__':

    example = 'non-blocking'  # 'non-blocking' or 'blocking'

    if myrank == 0:
        print('Hello I\'m the master, I\'ve %i slaves under my control'
              % nslaves)
        if example == 'non-blocking':
            master_work_nonblocking(nslaves)
        else:
            master_work_blocking(nslaves)

    else:
        print('... I\'m slave %i' % myrank)
        if example == 'non-blocking':
            slave_work_nonblocking(myrank)
        else:
            slave_work_blocking(myrank)