import matplotlib.pyplot as plt
import numpy as np

A = np.arange(1,5)
B = A**2

cnt=0
while(1):  
    cnt = cnt+1
    print("########### test %d ###########" % cnt)

    # here is the trick: 
    # set the figure a 'num' to prevent from re-malloc of a figure in the next loop 
    # and set "clear=True" to make the figure clear
    # I never use plt.close() to kill the figure, because I found it doesn't work.
    # Only one figure is allocated, which can be self-released when the program quits.
    # Before: 6000 times calling of plt.figure() ~ about 1.6GB of memory leak
    # Now: the memory keeps in a stable level

    fig = plt.figure(num=1, clear=True)
    ax0 = fig.add_subplot(2, 1, 1)
    ax1 = fig.add_subplot(2, 1, 2)
    #fig, ax = plt.subplots(1, 1)
    #plt.subplots_adjust(height_ratios=[1, 0.3])

    # alternatively use an other function in one line
    #fig, ax = plt.subplots(num=1,clear=True)

    ax0.plot(A,B)
    ax1.plot(B,A)

    # Here add the functions you need 
    # plt.show()
    fig.savefig('%d.png' % cnt)
    fig.clf() 
    del fig
    print(len(plt.get_fignums()), "figures still open")
    #plt.close(fig)
    #if cnt == 1:
    #    break
