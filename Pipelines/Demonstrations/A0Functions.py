import numpy as np


def make_Kde_able(x,y, bins,reduce):
    '''
    A function to turn a l,arge set of paired scsatter data into a smaller one by biniing into a 2d histogram
    :param x:
    :param y:
    :return:
    '''
    minx,maxx = min(x),max(x)
    miny, maxy = min(y), max(y)
    x_bin = (maxx - minx)/bins
    y_bin = (maxy - miny) / bins
    #print(minx, maxx,x_bin)
    #print(miny, maxy,y_bin)
    x_mids = []
    y_mids = []
    for i in range(bins):
        if i == 0:
            x_mids.append(minx + x_bin/2)
            y_mids.append(miny + y_bin / 2)
        else:
            x_mids.append(x_mids[i-1] + x_bin)
            y_mids.append(y_mids[i - 1] + y_bin)

    #print(x_mids)
    #print(y_mids)

    xya,xe,ye = np.histogram2d(x,y,bins=bins,density=False)
    xy = np.zeros((bins,bins))
    for i in range(bins):
        for j in range(bins):
            xy[i,j] = xya[j,i]

    #print(xy)
    #print(xe)
    #print(ye)


    newx = []
    newy = []

    for i in range(bins):
        for j in range(bins):
            x_val = x_mids[j]
            y_val = y_mids[i]
            count = int(xy[i,j]/reduce)
            #print(count)
            for k in range(count):
                newx.append(x_val)
                newy.append(y_val)
    #print(newx)
    #print(newy)
    return newx,newy





