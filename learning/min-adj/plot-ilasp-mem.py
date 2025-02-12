import matplotlib.pyplot as plt
from math import pow
from glob import glob

to_GiB = lambda x: (float(x) if x.isnumeric() else pow(1000,['k','m','g'].index(x[-1])+1) * float(x[:-1])) / 10**9

for ilasp_mem_log in glob('ilasp-mem*.log'):
    # Create data
    lines = open(ilasp_mem_log).readlines()
    x = range(0,10*len(lines),10)
    y1 = [to_GiB(x.split(',')[3]) for x in lines]
    y2 = [to_GiB(x.split(',')[4]) for x in lines]

    x_ticks = range(60*60,x[-1],60*60)
    x_labels = [str(tick // (60*60)) + ' hrs' for tick in x_ticks]

    # Basic stacked area chart.
    plt.stackplot(x, y1, y2, labels=['RSS', 'VIRT'])
    plt.xticks(x_ticks, x_labels)
    plt.legend(loc='upper left')
    plt.title('ILASP memory usage over time')

    if ilasp_mem_log == 'ilasp-mem_BIAS.log':
        plt.axhline(y = 768, color = 'r', linestyle = '-')

    plt.savefig(ilasp_mem_log.replace('log', 'png'))
    plt.clf()
