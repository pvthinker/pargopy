# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt


class MouseProfile(object):
    def __init__(self, line, tag, filename):
        self.line = line
        self.tag = tag
        #print(tag)
        self.click = None
        self.line.set_lw(1)
        self.pressed = False
        self.filename = filename

    def on_press(self, event):
        """ first click: highlight the line and indicate the tag number
        second click: confirm action on this profile"""
        #print(self.line.contains(event))
        #print('mouse pressed')
        if self.pressed:
            print('TAG= %i / confirmed' % self.tag)
            # in practice: do something on argodb = flag the profile
            # or store the action in a temporary file
            # and propagate the action in argodb with another routine
            with open(self.filename, 'a') as fid:
                print('write %i tag in %s' % (self.tag, self.filename))
                fid.write(str(self.tag)+'\n')

        if self.line.contains(event)[0]:
            self.line.set_lw(4)
            self.line.figure.canvas.draw()
            print('TAG = %i' % self.tag)
            self.pressed = True

    def on_release(self, event):
        'on release we reset the press data'
        self.line.set_lw(1)
        self.pressed = False
        self.line.figure.canvas.draw()

    def connect(self):
        self.cidclick = self.line.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.line.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)


if __name__ == '__main__':

    fig = plt.figure()
    ax = fig.add_subplot(111)
    nz = 63
    zref = np.linspace(0, 2000, nz)
    nprof = 10
    mps = []
    for k in range(nprof):
        tag = np.random.randint(10000)
        CT = np.cumsum(np.random.normal(size=(nz,)))
        li, = ax.plot(CT, zref, picker=5)
        print(li, tag)
        mp = MouseProfile(li, tag)
        mp.connect()
        mps.append(mp)
