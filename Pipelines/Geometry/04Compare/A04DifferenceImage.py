import math
import random

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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


def createIdealisedDifferencePlot(geoX,geoY,dataA, scatterise=False):
    if scatterise:
        return createIdealisedDifferencePlotScatterised(geoX,geoY,dataA)
    else:
        return createIdealisedDifferencePlotProbability(geoX, geoY, dataA)

def createIdealisedDifferencePlotProbability(geoX,geoY,dataA):
    #change the kde if their are too few observations
    xMin, yMin, xMax, yMax = 0, 0, 0, 0
    if len(dataA[geoX]) > 0:
        xMin = min(dataA[geoX])
        xMax = max(dataA[geoX])
    if len(dataA[geoY]) > 0:
        yMin = min(dataA[geoY])
        yMax = max(dataA[geoY])
    axesA = [xMin, xMax, yMin, yMax]

    bins = 50
    if len(dataA[geoX]) < 1000:
        bins = 10
    elif len(dataA[geoX]) < 5000:
        bins = 20

    histA,binsA = np.histogram(dataA[geoX],bins=bins,density=False)
    histB, binsC = np.histogram(dataA[geoY], bins=bins,density=False)

    arAA,xe,ye = np.histogram2d(dataA[geoX].values,dataA[geoY].values, bins=bins, density=False)
    #normA = np.linalg.norm(arAA)
    #arAA = arAA / normA
    arA = np.zeros((len(histB), len(histA)))
    for a in range(0,len(histA)):
        for b in range(0, len(histB)):
            arA[b,a] = arAA[a,b]

    arB = np.zeros((len(histB),len(histA)))
    for a in range(0,len(histA)):
        for b in range(0, len(histB)):
            arB[b,a] = histA[a] * histB[b]

    #normB = np.linalg.norm(arB)
    #arB = arB / normB

    # normalise A
    valA = 0
    for a in range(0, len(histA)):
        for b in range(0, len(histB)):
            valA += arA[a, b]
    for a in range(0, len(histA)):
        for b in range(0, len(histB)):
            arA[a, b] = arA[a, b]/valA

    # normalise B
    valB = 0
    for a in range(0,len(histA)):
        for b in range(0, len(histB)):
            valB += arB[a,b]
    for a in range(0,len(histA)):
        for b in range(0, len(histB)):
            arB[a,b] = arB[a,b]/valB

    #print('Total mat B=',valB)

    arDiff = np.zeros((len(histB), len(histA)))
    minVal, maxVal = 0, 0
    stat,count = 0,0
    for a in range(0,len(histA)):
        for b in range(0, len(histB)):
            arDiff[a,b] = arA[a,b] - arB[a,b]
            maxVal = max(maxVal, arA[a,b],arB[a,b])
            minVal = min(minVal, arA[a,b],arB[a,b])
            if arA[a,b] != 0 or arB[a,b] != 0:
                count += 1
                stat += abs(arDiff[a,b])
            #print('2d', arA[b, a], arB[a, b])

    stat = stat/2
    print(stat,count)
    maxVal = max(abs(maxVal), abs(minVal))
    minVal = 0 - maxVal
    return arA,arB,arDiff,minVal,maxVal, stat,''



def createIdealisedDifferencePlotScatterised(geoX,geoY,dataA):
    # First create a randomised version of dataA
    cutA = list(dataA[geoX].values)
    random.shuffle(cutA)
    cutB = list(dataA[geoY].values)
    random.shuffle(cutB)
    dic_cut = {}
    dic_cut[geoX] = cutA
    dic_cut[geoY] = cutB
    dataB = pd.DataFrame.from_dict(dic_cut)
    dataB = dataB.sort_values(by=geoX, ascending=False)

    bins = 20

    histAA,binsAA = np.histogram(dataA[geoX],bins=bins,density=False)
    histBA, binsCA = np.histogram(dataA[geoY], bins=bins,density=False)

    histAB, binsAB = np.histogram(dataB[geoX], bins=bins, density=False)
    histBB, binsCB = np.histogram(dataB[geoY], bins=bins, density=False)

    arAA,xe,ye = np.histogram2d(dataA[geoX].values,dataA[geoY].values, bins=bins, density=False)
    normA = np.linalg.norm(arAA)
    arAA = arAA / normA
    arA = np.zeros((len(histBA), len(histAA)))
    for a in range(0,len(histAA)):
        for b in range(0, len(histBA)):
            if arAA[a,b] >0:
                arA[b,a] = 1

    arBB, xe, ye = np.histogram2d(dataB[geoX].values, dataB[geoY].values, bins=bins, density=False)
    normB = np.linalg.norm(arBB)
    arBB = arBB / normB
    arB = np.zeros((len(histBB), len(histAB)))
    for a in range(0, len(histAB)):
        for b in range(0, len(histBB)):
            if arBB[a, b] > 0:
                arB[b, a] = 1

    arDiff = np.zeros((len(histBA), len(histAA)))
    minVal, maxVal = 0, 0
    stat,count = 0,0
    for a in range(0,len(histAA)):
        for b in range(0, len(histBA)):
            arDiff[a,b] = arA[a,b] - arB[a,b]
            if arA[a,b] != 0 or arB[a,b] != 0:
                count += 1
                stat += abs(arDiff[a,b])
    stat = stat/count
    print(stat,count)
    maxVal = 1
    minVal = -1
    return arA,arB,arDiff,minVal,maxVal, stat,''






