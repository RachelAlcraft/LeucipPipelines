import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

bins = 50
kde_val = 0.05

def kde2D_scipy(bandwidth, axes, bins, data,geoX,geoY):
    xdata = data[geoX]
    ydata = data[geoY]
    data = np.vstack([xdata, ydata])
    xgrid = np.linspace(axes[0], axes[1], bins)
    ygrid = np.linspace(axes[2], axes[3], bins)
    Xgrid, Ygrid = np.meshgrid(xgrid, ygrid)
    grid_sized = np.vstack([Xgrid.ravel(), Ygrid.ravel()])
    # fit an array of size [Ndim, Nsamples]

    kde = gaussian_kde(data, bw_method=bandwidth)
    # evaluate on a regular grid
    Z = kde.evaluate(grid_sized)
    zgrid = Z.reshape(Xgrid.shape)

    Znorm = np.linalg.norm(zgrid)
    zgrid = zgrid / Znorm

    return xgrid, ygrid, zgrid

def createDifferencePlot(geoX,geoY,dataA,dataB):

    xMin, yMin, xMax, yMax = 0, 0, 0, 0
    if len(dataA[geoX]) > 0:
        xMin = min(dataA[geoX])
        xMax = max(dataA[geoX])
    if len(dataA[geoY]) > 0:
        yMin = min(dataA[geoY])
        yMax = max(dataA[geoY])
    axesA = [xMin, xMax, yMin, yMax]

    xMin, yMin, xMax, yMax = 0, 0, 0, 0
    if len(dataB[geoX]) > 0:
        xMin = min(dataB[geoX])
        xMax = max(dataB[geoX])
    if len(dataB[geoY]) > 0:
        yMin = min(dataB[geoY])
        yMax = max(dataB[geoY])
    axesB = [xMin, xMax, yMin, yMax]

    xyzA = kde2D_scipy(kde_val, axesA, bins,dataA,geoX,geoY)
    xyzB = kde2D_scipy(kde_val, axesB, bins,dataB,geoX,geoY)

    arA = xyzA[2]
    arB = xyzB[2]

    arDiff = arA - arB
    minVal,maxVal=0,0
    for i in range(arA.shape[0]):
        for j in range(arA.shape[1]):
            maxVal = max(maxVal, arA[i, j])
            minVal = min(minVal, arA[i, j])

    for i in range(arB.shape[0]):
        for j in range(arB.shape[1]):
            maxVal = max(maxVal, arB[i, j])
            minVal = min(minVal, arB[i, j])



    stat = 0
    count = 0
    for i in range(arDiff.shape[0]):
        for j in range(arDiff.shape[1]):
            #if arA[i, j] != 0 and arB[i, j] != 0:
            stat += (arDiff[i, j]*arDiff[i, j])
            count += 1

    start = stat/count

    maxVal = max(abs(maxVal), abs(minVal))
    minVal = 0 - maxVal

    return arA,arB,arDiff,minVal,maxVal, stat





def createIdealisedDifferencePlot(geoX,geoY,dataA):

    xMin, yMin, xMax, yMax = 0, 0, 0, 0
    if len(dataA[geoX]) > 0:
        xMin = min(dataA[geoX])
        xMax = max(dataA[geoX])
    if len(dataA[geoY]) > 0:
        yMin = min(dataA[geoY])
        yMax = max(dataA[geoY])
    axesA = [xMin, xMax, yMin, yMax]

    histA,binsA = np.histogram(dataA[geoX],bins=bins,density=False)
    histB, binsC = np.histogram(dataA[geoY], bins=bins,density=False)


    #xyzA = kde2D_scipy(kde_val, axesA, bins,dataA,geoX,geoY)
    #arA = xyzA[2]

    arAA,xe,ye = np.histogram2d(dataA[geoX].values,dataA[geoY].values, bins=bins, density=False)
    normA = np.linalg.norm(arAA)
    arAA = arAA / normA
    arA = np.zeros((len(histB), len(histA)))
    for a in range(0,len(histA)):
        for b in range(0, len(histB)):
            arA[b,a] = arAA[a,b]


    arB = np.zeros((len(histB),len(histA)))
    for a in range(0,len(histA)):
        for b in range(0, len(histB)):
            arB[b,a] = histA[a] * histB[b]

    normB = np.linalg.norm(arB)
    arB = arB / normB

    arDiff = arA - arB
    arDiffSq = arA - arB
    minVal,maxVal=0,0
    for i in range(arA.shape[0]):
        for j in range(arA.shape[1]):
            maxVal = max(maxVal, arA[i, j])
            minVal = min(minVal, arA[i, j])

    for i in range(arB.shape[0]):
        for j in range(arB.shape[1]):
            maxVal = max(maxVal, arB[i, j])
            minVal = min(minVal, arB[i, j])

    stat = 0
    count = 0
    for i in range(arDiff.shape[0]):
        for j in range(arDiff.shape[1]):
            #arDiffSq[i,j] = arDiff[i,j]*arDiff[i,j]
            arDiffSq[i, j] = abs(arDiff[i, j])
            if arA[i, j] != 0 and arB[i, j] != 0:
                count += 1
            #if arDiff[i, j] > 0:
                stat += abs(arDiff[i, j])


    #print(geoX,geoY,stat,count)
    stat = stat/count
    stat = math.sqrt(stat)

    #normSq = np.linalg.norm(arDiffSq)
    #arDiffSq = arDiffSq / normSq

    maxVal = max(abs(maxVal), abs(minVal))
    minVal = 0 - maxVal

    return arA,arB,arDiff,minVal,maxVal, stat,arDiffSq







