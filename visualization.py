import sys
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import ttk
import os
import algorithm as alg
import copy
import divideAndConquer as dc
import numpy as np

message1 = 'Insert Height'

xOfVertices = []
yOfVertices = []
xTower = 0
yBase = 0
yTower = 0
h = 0
numberTower = -1

setP = []
setVisP = []
terrainWithP = []

setNotVisibleEdges = []

def loadTerrain():
    global xOfVertices, yOfVertices, xTower, yBase, yTower, numberTower,setP
    try:
        script = myTerrains.get(myTerrains.curselection())
        p = terrainPath+"\\terrains\\" + script + ".txt"
        defaultTerrain = open(p, "r")
        
        vertices = defaultTerrain.read().splitlines()
        defaultTerrain.close()

        setP = []
        xOfVertices = []
        yOfVertices = []
        xTower = 0
        yBase = 0
        yTower = 0
        numberTower = -1

        for vertice in vertices:
            coords = vertice.split(" ")
            xOfVertices.append(float(coords[0]))
            yOfVertices.append(float(coords[1]))

        winLoadTerrain.withdraw()   
        showTerrain() 
    except:
        print("Terrain file moved")

def drawTerrain():
        global xOfVertices,yOfVertices
        
        plt.figure(num='The two watchtowers')
        plt.title("General terrain")
        plt.plot(xOfVertices,yOfVertices,color = 'green')
        plt.plot(xOfVertices, yOfVertices, marker = ".", color = "green")
        if numberTower != -1:
            plt.plot([xTower,xTower],[yBase,yTower],color = '0')
            plt.plot([xTower],[yTower], marker = "s")
        plt.draw()

def showTerrain():
    plt.close()
    plt.close()
    plt.close()
    drawTerrain()
    plt.show()

def updateTower():
    global xTower, yBase, yTower, numberTower,h,setP
    setP = []
    
    try:
        xTower = float(listOfVertices.get(listOfVertices.curselection()))
        numberTower = xOfVertices.index(xTower)
        yBase = yOfVertices[numberTower]
        h = float(text_boxHeight.get("1.0",END))
        yTower = yBase + h
        winVertices.withdraw()
        
    except:
        print("Not valid vertice or height")
    showTerrain()

def autoUpdateTower(x,y):
    global xTower, yBase, yTower, numberTower,h,setP
    setP = []
    xTower = x
    yTower = y
    numberTower = xOfVertices.index(xTower)
    yBase = yOfVertices[numberTower]
    h = yTower - yBase

def chooseTerrain():
    winLoadTerrain.deiconify()   

def chooseVertice(): 
    global winVertices, listOfVertices, text_boxHeight
    listOfVertices.delete(0,END)
    winVertices.deiconify()
    for x in xOfVertices:  
        listOfVertices.insert(END,  str(x))  

#Code starts here    

def autoCalculateVisibilityLines():
    # active input: nothing
    # passive input: terrain T as global
    # outputs(all global):
    #          setP - p points+ not visible vertices
    #          terrainWithP - T + p points
    #          setNotVisibleEdges - all invisible edges and partly invisiblle edges
    #          setVisP - p points 

    global setP, terrainWithP ,setNotVisibleEdges,setVisP
    setNotVisibleEdges = []

    endPoints = alg.visibilityEndpoints(xOfVertices, yOfVertices, xTower, yTower) 
    #returns [[p points- [x,y]], [endpoints and vertices out of sight- [x,y]],[vertices in sight where visibility ends-[x,y]]
    
    setPonVertice= []
    for endPoint in endPoints[0]:
        if endPoint[0] in xOfVertices:
            setPonVertice.append(endPoint)
    setP = endPoints[0]   
    setVisP = copy.deepcopy(setP)  
    for i in range(len(setP)):
        setP[i][0] = round(setP[i][0],13)
        setP[i][1] = round(setP[i][1],13)

    #calculate terrainWithP  -avoid double vertices 
    terrainWithP = [xOfVertices.copy(),yOfVertices.copy()]
    setPwithNoVertice = copy.deepcopy(setP)
    epCounter = 0
    if len(setPwithNoVertice)>0:
        for i in range(len(xOfVertices)):
            if len(setPwithNoVertice)>0:
                if epCounter >= len(setPwithNoVertice):
                        break     
                if  xOfVertices[i] == setPwithNoVertice[epCounter][0]:
                    setPwithNoVertice.remove(setPwithNoVertice[epCounter])
                    continue
                if  xOfVertices[i] > setPwithNoVertice[epCounter][0]:
                    terrainWithP[0].insert(epCounter+i,setPwithNoVertice[epCounter][0])
                    terrainWithP[1].insert(epCounter+i,setPwithNoVertice[epCounter][1])
                    epCounter+=1            

    for p in setPonVertice:
        endPoints[1].remove(p)
    #calculate not visible edges in setNotVisibleEdges
    for i in range (len(endPoints[1])-1):
        if endPoints[1][i] in endPoints[2] and endPoints[1][i+1]in endPoints[2]:
            continue
        if endPoints[1][i][0] in xOfVertices and endPoints[1][i+1][0] in xOfVertices:
            start = xOfVertices.index(endPoints[1][i][0])
            finish = xOfVertices.index(endPoints[1][i+1][0])
            if abs(start - finish) == 1:
                setNotVisibleEdges.append([endPoints[1][i],endPoints[1][i+1]])
        if endPoints[1][i][0] in xOfVertices and endPoints[1][i+1] in setP and endPoints[1][i+1] not in setPonVertice:
            setNotVisibleEdges.append([endPoints[1][i],endPoints[1][i+1]])
    
    #create setP (it works so I don't touch it)
    #In the beggining I wanted only p points then I added all fully invisible V + p-points
    setPnotValuable = copy.deepcopy(setPonVertice)
    for point in setPnotValuable:
        counter = 0
        for pair in setNotVisibleEdges:
            for p in pair:
                if p == point:
                    counter+=1
        if counter==1:
            setPnotValuable.remove(point)
    for p in setP:
        if p in setPnotValuable:
            setP.remove(p)
    for hiddenPiece in setNotVisibleEdges:
        for v in hiddenPiece:
            if v not in setP and v not in endPoints[2]:
                setP.append(v)
    setP = sorted(setP, key=lambda p:p[0])

#input: global terrainWithP and setP
#output: left PV pairs
def autoCalculateLeftPUV():
    puv = dc.constructPUV(xOfVertices, terrainWithP,setP)
    setE = []
    for part in setNotVisibleEdges:
        setE.append(part)
    return [puv,setE]


#input: global terrainWithP and setP
#output: right PV pairs
def autoCalculateRightPUV():
    # idea is to mirror terrainWithP and setP and then calulate left PV pairs in the mirrored terrainWithP
    # afte mirroring the outputs you get right PV pairs
    rightTerrain = copy.deepcopy(terrainWithP)
    
    rightSetP = copy.deepcopy(setP)
    rightTerrain[0] = [rightTerrain[0][-1] - x for x in rightTerrain[0]]
    
    for i in range(len(rightSetP)):
        rightSetP[i][0] = terrainWithP[0][-1] - rightSetP[i][0]

    rightTerrain[0].reverse()
    rightTerrain[1].reverse()
    rightSetP.reverse()
    revXofVertices = copy.deepcopy(xOfVertices)
    revXofVertices = [xOfVertices[-1] - x for x in xOfVertices]
    rightPUV = dc.constructPUV(revXofVertices, rightTerrain,rightSetP)
    for i in range(len(rightPUV)):
        rightPUV[i][0][0] = xOfVertices[-1]-rightPUV[i][0][0]
        rightPUV[i][1][0] = xOfVertices[-1]-rightPUV[i][1][0]
    setE = []
    for part in setNotVisibleEdges:
        setE.append(part)
    rightPV = [rightPUV,setE]
    return rightPV

                       
#inputs:globals xOfVertices,yOfVertices
#outputs: shortest common height pair of towers
def calculateSecondTowerDescriteSolo():
    edges = []
    for i in range(len(xOfVertices)-1):
        edges.append([[xOfVertices[i],yOfVertices[i]],[xOfVertices[i+1],yOfVertices[i+1]]])
    tower1 = -1
    tower2 = -1
    height1 = sys.float_info.max 
    height2 = sys.float_info.max 
    for i in range(len(xOfVertices)-1):   
        setE = edges
        setA = []
        criticalPairs = calculateCriticalHeights(xOfVertices[i])
        for pair in criticalPairs:
            setA.append(pair)
        interEdgeTower = []
        towerLine = [[xOfVertices[i],yOfVertices[i]],[xOfVertices[i],yOfVertices[i]+2]]
        for line in setE:
            temp = alg.intersectionPoint(line ,towerLine)[1]
            if temp >= yOfVertices[i]:
                interEdgeTower.append(temp)
        for line in setA:
            temp = alg.intersectionPoint(line ,towerLine)[1]
            if temp >= yOfVertices[i]:
                interEdgeTower.append(temp)
        interEdgeTower.sort()
        setL = copy.deepcopy(interEdgeTower)
        left = 0
        right = len(setL)-1
        currentTower2 = -1
        currentHeight2 = sys.float_info.max 
        found = FALSE
        iBestCoverTower = -1
        iBestCoverHeight = sys.float_info.max 
        iBestHeight = sys.float_info.max 
        while not found and left <= right:
            middle = (left + right)//2
            autoUpdateTower(xOfVertices[i],setL[middle])
            currentTower1 = xOfVertices[i]
            currentHeight1 = h
            bestTower =  calculateCoverTower()
            currentTower2 = bestTower[0]
            currentHeight2 = bestTower[1]

            if max(currentHeight1,currentHeight2) < max(iBestHeight,iBestCoverHeight):
                iBestCoverTower = currentTower2
                iBestCoverHeight = currentHeight2
                iBestHeight = currentHeight1

            
            if currentHeight1 < currentHeight2:
                left = middle+1
            elif currentHeight1 > currentHeight2:
                right = middle -1
            else:
                found = TRUE
        xwinner = i     #xOfVertices.index(tower1)
        x2winner = xOfVertices.index(iBestCoverTower)
        if iBestHeight!=iBestCoverHeight:
            currentBest = doubleCheckWinnerDiscrete(xwinner,x2winner,iBestHeight,iBestCoverHeight)
            currentTower1 = currentBest[0]
            currentTower2 = currentBest[1]  
            currentHeight1 = currentBest[2]  
            currentHeight2 = currentBest[3]  
        if max(currentHeight1,currentHeight2) < max(height1,height2):
                tower1 = currentTower1
                tower2 = currentTower2
                height1 = round(currentHeight1,5)
                height2 = round(currentHeight2,5)
        print(xOfVertices[i], tower1,tower2,height1,height2)
    print("x tower1:",tower1, "height1:", height1)
    print("x tower2:",tower2,"height2:", height2)
    return [[tower1,height1],[tower2,height2]]


def doubleCheckWinnerDiscrete(i,j,height1,height2):
    tower1 = xOfVertices[i]
    tower2 = xOfVertices[j]
    height1 = height1 
    height2 = height2
    if height1>height2:
        dh = height1-height2
        low = height2
        high = height1
    else:
        dh = height2-height1
        low = height1
        high = height2
    step = dh/10000
    setL = np.arange(yOfVertices[i]+low, yOfVertices[i]+high, step).tolist()
    left = 0
    right = len(setL)-1
    currentTower2 = -1
    currentHeight2 = sys.float_info.max 
    found = FALSE
    while not found and left <= right:
        middle = (left + right)//2
        autoUpdateTower(xOfVertices[i],setL[middle])
        currentTower1 = xOfVertices[i]
        currentHeight1 = h
        bestTower =  calculateCoverTower()
        currentTower2 = bestTower[0]
        currentHeight2 = bestTower[1]
        if max(currentHeight1,currentHeight2) < max(height1,height2):
            tower1 = currentTower1
            tower2 = currentTower2
            height1 = currentHeight1
            height2 = currentHeight2
        elif max(currentHeight1,currentHeight2) == max(height1,height2):
            if min(currentHeight1,currentHeight2) < min(height1,height2):
                tower1 = currentTower1
                tower2 = currentTower2
                height1 = currentHeight1
                height2 = currentHeight2
        if currentHeight1 < currentHeight2:
            left = middle+1
        elif currentHeight1 > currentHeight2:
            right = middle -1
        else:
            found = TRUE
    return [tower1,tower2,height1,height2]


#input: global first tower
# output: shortest cover Tower
def calculateCoverTower():
    global setP
    autoCalculateVisibilityLines()
    for p in setP:
        if p[0]< numberTower+1:
            setP.remove(p)
    visibilityPairs = autoCalculateLeftPUV()
    leftPairsPV = visibilityPairs[0]
    setB = []
    setEv = visibilityPairs[1]
    currentTower2 = 0
    currentHeight2 = sys.float_info.max 
    for j in range(numberTower+1,len(xOfVertices)):
        
        secondTowerLine = [[xOfVertices[j],yOfVertices[j]],[xOfVertices[j],yOfVertices[j]+2]]
        possibleHeights2 = []

        for line in setB+setEv:
            if line[0][0] != line[1][0]:
                temp = alg.intersectionPoint(line ,secondTowerLine)[1]
                if round(temp,10) >= yOfVertices[j]:
                    possibleHeights2.append(round(temp,10))
        if len(possibleHeights2) > 0:
            currentY2 = max(possibleHeights2)
            tempHeight2 = currentY2 - yOfVertices[j]
        else:
            tempHeight2 = sys.float_info.max 
        if tempHeight2 < currentHeight2:
            currentHeight2 = tempHeight2
            currentTower2 = j
        #sorting is difficult
        for pair in leftPairsPV:
            if pair[0][0] == xOfVertices[j] or pair[1][0] == xOfVertices[j]:
                setB.append(pair)
    if currentHeight2 == sys.float_info.max:
        currentHeight2 = 0
    return [xOfVertices[currentTower2],currentHeight2]                       


#inputs:globals xOfVertices,yOfVertices
#outputs: shortest common height pair of towers
def calculateSecondTowerSemiContiniousSolo(): #n^3(log^2 n)
    edges = []
    for i in range(len(xOfVertices)-1):
        edges.append([[xOfVertices[i],yOfVertices[i]],[xOfVertices[i+1],yOfVertices[i+1]]])
    tower1 = -1
    base1 = -1
    tower2 = -1
    base2 = -1
    height1 = sys.float_info.max 
    height2 = sys.float_info.max 
    upperEnvelope = []
    for i in range(len(xOfVertices)-1):   
        setE = []
        for e in edges:
            setE.append(e)
        setA = []
        criticalPairs = calculateCriticalHeights(xOfVertices[i])
        for pair in criticalPairs:
            setA.append(pair)
        interEdgeTower = []
        towerLine = [[xOfVertices[i],yOfVertices[i]],[xOfVertices[i],yOfVertices[i]+2]]
        for line in setE:
            temp = alg.intersectionPoint(line ,towerLine)[1]
            if temp >= yOfVertices[i]:
                interEdgeTower.append(temp)
        for line in setA:
            temp = alg.intersectionPoint(line ,towerLine)[1]
            if temp >= yOfVertices[i]:
                interEdgeTower.append(temp)
        interEdgeTower.sort()
        setL = copy.deepcopy(interEdgeTower)
        left = 0
        right = len(setL)-1
        currentTower2 = -1
        currentBase2 = -1
        currentHeight2 = sys.float_info.max 
        found = FALSE
        iBestTower = -1
        iBestCoverTower = -1
        iBestCoverHeight = sys.float_info.max 
        iBestHeight = sys.float_info.max 
        iBestEnvelope = []
        iBestBase = -1
        iBestCoverBase = -1
        while not found and left <= right:
            middle = (left + right)//2
            autoUpdateTower(xOfVertices[i],setL[middle])
            currentTower1 = xOfVertices[i]
            currentBase1 = yOfVertices[i]
            currentHeight1 = round(h,8)
            bestTower =  calculateSecondTowerSemiContinious()
            currentTower2 = bestTower[0]
            currentBase2 = bestTower[1]
            currentHeight2 = bestTower[2]
            currentUpperEnvelope = bestTower[3]
            if max(currentHeight1,currentHeight2) < max(iBestHeight,iBestCoverHeight):
                iBestTower = currentTower1
                iBestCoverTower = currentTower2
                iBestBase = currentBase1
                iBestCoverBase = currentBase2
                iBestHeight = currentHeight1
                iBestCoverHeight = currentHeight2
                iBestEnvelope = currentUpperEnvelope
        
            if currentHeight1 < currentHeight2:
                left = middle+1
            elif currentHeight1 > currentHeight2:
                right = middle -1
            else:
                found = TRUE
        xwinner = i     #xOfVertices.index(tower1)
        x2winner = iBestCoverTower
        if iBestHeight!=iBestCoverHeight:
            currentBest = doubleCheckWinnerSemiContinious(xwinner,x2winner,iBestHeight,iBestCoverHeight,iBestBase,iBestCoverBase,iBestEnvelope)
            currentTower1 = currentBest[0][0]
            currentBase1 = currentBest[0][1]  
            currentHeight1 = currentBest[0][2]  
            currentTower2 = currentBest[1][0]  
            currentBase2 = currentBest[1][1]  
            currentHeight2 = currentBest[1][2]
            currentUpperEnvelope = currentBest[2]  
        if max(currentHeight1,currentHeight2) < max(height1,height2):
            tower1 = currentTower1
            tower2 = currentTower2
            base1 = currentBase1
            base2 = currentBase2
            height1 = round(currentHeight1,5)
            height2 = round(currentHeight2,5)
            upperEnvelope = currentUpperEnvelope
        print(xOfVertices[i], tower1,tower2,height1,height2)
        
    print("x tower1:",tower1, "height1:", height1)
    print("x tower2:",tower2,"height2:", height2, "y base2:",base2)
    return([[tower1,base1,height1],[tower2,base2,height2],upperEnvelope])

def doubleCheckWinnerSemiContinious(i,j,height1,height2,base1,base2,envelope):
    tower1 = xOfVertices[i]
    tower2 = j
    height1 = height1 
    height2 = height2
    base1 = base1
    base2 = base2
    upperEnvelope = envelope
    if height1>height2:
        dh = height1-height2
        low = height2
        high = height1
    else:
        dh = height2-height1
        low = height1
        high = height2
    step = dh/10000
    setL = np.arange(yOfVertices[i]+low, yOfVertices[i]+high, step).tolist()
    left = 0
    right = len(setL)-1
    currentTower2 = -1
    currentHeight2 = sys.float_info.max
    found = FALSE
    while not found and left <= right:
        middle = (left + right)//2
        autoUpdateTower(xOfVertices[i],setL[middle])
        currentTower1 = xOfVertices[i]
        currentBase1 = yOfVertices[i]
        currentHeight1 = round(h,8)
        bestTower =  calculateSecondTowerSemiContinious()
        currentTower2 = bestTower[0]
        currentBase2 = bestTower[1]
        currentHeight2 = bestTower[2]
        currentUpperEnvelope = bestTower[3]
        if max(currentHeight1,currentHeight2) < max(height1,height2):
            tower1 = currentTower1
            tower2 = currentTower2
            base1 = currentBase1
            base2 = currentBase2
            height1 = currentHeight1
            height2 = currentHeight2
            upperEnvelope = currentUpperEnvelope
        elif max(currentHeight1,currentHeight2) == max(height1,height2):
            if min(currentHeight1,currentHeight2) < min(height1,height2):
                tower1 = currentTower1
                tower2 = currentTower2
                base1 = currentBase1
                base2 = currentBase2
                height1 = currentHeight1
                height2 = currentHeight2
                upperEnvelope = currentUpperEnvelope
        if currentHeight1 < currentHeight2:
            left = middle+1
        elif currentHeight1 > currentHeight2:
            right = middle -1
        else:
            found = TRUE
    return([[tower1,base1,height1],[tower2,base2,height2],upperEnvelope])

#input: global first tower
#output: shortest cover Tower
def calculateSecondTowerSemiContinious():
    autoCalculateVisibilityLines() #O(n)
    visibilityPairsLeft = autoCalculateLeftPUV() #O(n^2logn)
    visibilityPairsRight = autoCalculateRightPUV() #O(n^2logn)
    setB = []
    setC = visibilityPairsRight[0]
    for pair in setC:
        if pair[0][0] == xOfVertices[0] or pair[1][0] == xOfVertices[0]:
            setC.remove(pair)
    setNotVisE = setNotVisibleEdges
    currentTower2 = 0
    currentBase2 = 0
    currentHeight2 = sys.float_info.max 
    currentUpE =[]
    edges = []
    setLines = setB+setC+setNotVisE
    upE = calculateUpperEnvelope(setLines)
    possibleHeights2 = []
    secondTowerLine = [[xOfVertices[0],yOfVertices[0]],[xOfVertices[0],yOfVertices[0]+2]]
    for line in setLines:
        if line[0][0] != line[1][0]:
            temp = alg.intersectionPoint(line ,secondTowerLine)[1]
            if round(temp,10) >= yOfVertices[0]:
                possibleHeights2.append(round(temp,8))
    if len(possibleHeights2) > 0:
        currentY2 = max(possibleHeights2)
        tempHeight2 = currentY2 - yOfVertices[0]
        tempBase2 = yOfVertices[0]
    else:
        tempHeight2 = sys.float_info.max 
    if tempHeight2 < currentHeight2:
        currentHeight2 = round(tempHeight2,8)
        currentTower2 = xOfVertices[0]
        currentBase2 = tempBase2
        currentUpE = upE
    for i in range(len(xOfVertices)-1):
        edges.append([[xOfVertices[i],yOfVertices[i]],[xOfVertices[i+1],yOfVertices[i+1]]])
    for e in edges:
        setLines = setB+setC+setNotVisE
        upE = calculateUpperEnvelope(setLines) #n^2logn
        if upE == 0: #if first tower sees all
            return [0,0,0,[]]
        for k in range(len(upE[0])):
            if upE[0][k] >= e[0][0] and upE[0][k] <= e[1][0]:
                secondTowerLineInUpEnv = [[upE[0][k],upE[1][k]],[upE[0][k],upE[1][k]+1]]
                tempBase = alg.intersectionPoint(e ,secondTowerLineInUpEnv)[1]
                tempHeight2 = round(upE[1][k] - tempBase,8)
                if tempHeight2<0:
                    tempHeight2=0
                if tempHeight2 < currentHeight2:
                    currentHeight2 = round(tempHeight2,8)
                    currentTower2 = round(upE[0][k],8)
                    currentBase2 = round(tempBase,8)
                    currentUpE = upE

        for pair in visibilityPairsLeft[0]:
            if pair[0][0] == e[1][0] or pair[1][0] == e[1][0]:
                setB.append(pair)
        setLines = setB+setC+setNotVisE
        possibleHeights2 = []
        secondTowerLine = [[e[1][0],e[1][1]],[e[1][0],e[1][1]+2]]
        for line in setLines:
            if line[0][0] != line[1][0]:
                temp = alg.intersectionPoint(line ,secondTowerLine)[1]
                if round(temp,8) >= e[1][1]:
                    possibleHeights2.append(round(temp,10))
        if len(possibleHeights2) > 0:
            currentY2 = max(possibleHeights2)
            tempHeight2 = currentY2 - e[1][1]
            if tempHeight2<0:
                    tempHeight2=0
            tempBase2 = e[1][1]
            if tempHeight2 < currentHeight2:
                currentHeight2 = tempHeight2
                currentTower2 = e[1][0]
                currentBase2 = e[1][1]
                currentUpE = upE
    return [currentTower2,currentBase2,currentHeight2,currentUpE]


    


def calculateUpperEnvelope(setLines):
    upperEnvelope = [[],[]]
    upperEnvelopeLines = []
    #sort lines based on slope
    setLines = sorted(setLines, key=lambda line: (line[1][1] - line[0][1]) / (line[1][0] - line[0][0]))
    slopes = []
    for i in range(len(setLines)):
        slopes.append(round((setLines[i][0][1]-setLines[i][1][1])/(setLines[i][0][0]-setLines[i][1][0]),8))
    slopes.sort()
    #remove parallel Lines
    i = 0 
    while i < (len(slopes)-1):
       
        if slopes[i] == slopes[i+1]:
            if slopes[i]>0:
                if setLines[i][0][0] < setLines[i+1][0][0]:
                    setLines.pop(i+1)
                    slopes.pop(i+1)
                else:
                    setLines.pop(i)
                    slopes.pop(i)
            elif slopes[i]==0:
                if setLines[i][0][1] < setLines[i+1][0][1]:
                    setLines.pop(i)
                    slopes.pop(i)
                else:
                    setLines.pop(i+1)
                    slopes.pop(i+1)
            else:
                if setLines[i][0][0] > setLines[i+1][0][0]:
                    setLines.pop(i+1)
                    slopes.pop(i+1)
                else:
                   setLines.pop(i)
                   slopes.pop(i)
        else:
            i+=1
    #catch special cases
    if len(setLines)==0: #all visible
        return 0
    elif len(setLines)==1: #not visible edges with same slope
        temp = alg.intersectionPoint([[xOfVertices[0],yOfVertices[0]],[xOfVertices[0],yOfVertices[0]+1]],setLines[0])
        upperEnvelope[0].insert(0,temp[0])
        upperEnvelope[1].insert(0,temp[1])
        temp = alg.intersectionPoint([[xOfVertices[-1],yOfVertices[-1]],[xOfVertices[-1],yOfVertices[-1]+1]],setLines[0])
        upperEnvelope[0].append(temp[0])
        upperEnvelope[1].append(temp[1])
        return upperEnvelope
    #calculate Upper Envelope raw
    temp = alg.intersectionPoint(setLines[0],setLines[1])
    upperEnvelope[0].append(temp[0])
    upperEnvelope[1].append(temp[1])
    upperEnvelopeLines.append(setLines[0])
    upperEnvelopeLines.append(setLines[1])
    for line in setLines[2:]:
        flag = TRUE
        while len(upperEnvelopeLines) > 0 and flag :
            temp = alg.intersectionPoint(line,upperEnvelopeLines[-1])
            if len(upperEnvelope[0]) == 0:
                break
            if temp[0] <= upperEnvelope[0][-1] :
                upperEnvelopeLines.pop()
                upperEnvelope[0].pop()
                upperEnvelope[1].pop()
            else:
                flag = FALSE
        
        upperEnvelope[0].append(temp[0])
        upperEnvelope[1].append(temp[1])
        upperEnvelopeLines.append(line)   
    #shape upper envelope to the xdimensions of T
    if len(upperEnvelope[0])>1:
        while upperEnvelope[0][1] < xOfVertices[0]:
            upperEnvelope[0].pop(0)
            upperEnvelope[1].pop(0)
    if upperEnvelope[0][0]>xOfVertices[0]:
        temp = alg.intersectionPoint([[xOfVertices[0],yOfVertices[0]],[xOfVertices[0],yOfVertices[0]+1]],upperEnvelopeLines[0])
        upperEnvelope[0].insert(0,temp[0])
        upperEnvelope[1].insert(0,temp[1])
    else:
        if len(upperEnvelope[0])>1:
            tempLine = [[upperEnvelope[0].pop(0),upperEnvelope[1].pop(0)],[upperEnvelope[0][0],upperEnvelope[1][0]]]
            temp = alg.intersectionPoint([[xOfVertices[0],yOfVertices[0]],[xOfVertices[0],yOfVertices[0]+1]],tempLine)
            upperEnvelope[0].insert(0,temp[0])
            upperEnvelope[1].insert(0,temp[1])
    if len(upperEnvelope[0])>1:
        while upperEnvelope[0][-2] > xOfVertices[-1]:
            upperEnvelope[0].pop()
            upperEnvelope[1].pop()
    if upperEnvelope[0][-1]<xOfVertices[-1]:
        temp = alg.intersectionPoint([[xOfVertices[-1],yOfVertices[-1]],[xOfVertices[-1],yOfVertices[-1]+1]],upperEnvelopeLines[-1])
        upperEnvelope[0].append(temp[0])
        upperEnvelope[1].append(temp[1])
    else:
        tempLine = [[upperEnvelope[0].pop(),upperEnvelope[1].pop()],[upperEnvelope[0][-1],upperEnvelope[1][-1]]]
        temp = alg.intersectionPoint([[xOfVertices[-1],yOfVertices[-1]],[xOfVertices[-1],yOfVertices[-1]+1]],tempLine)
        upperEnvelope[0].append(temp[0])
        upperEnvelope[1].append(temp[1])
    #shape in case first or last line is under T meaning whole edges of T are visible
    if upperEnvelope[1][-1]< yOfVertices[-1]:
        upperEnvelope[0].pop()
        upperEnvelope[1].pop()
        for i in range(len(xOfVertices)):
            if xOfVertices[i]> upperEnvelope[0][-1]:
                upperEnvelope[0].append(xOfVertices[i])
                upperEnvelope[1].append(yOfVertices[i])
    if upperEnvelope[1][0]< yOfVertices[0]:
        upperEnvelope[0].pop(0)
        upperEnvelope[1].pop(0)
        for i in range(len(xOfVertices)):
            if xOfVertices[i] < upperEnvelope[0][-1]:
                upperEnvelope[0].append(xOfVertices[i])
                upperEnvelope[1].append(yOfVertices[i])
    
    return(upperEnvelope)


#------------Calculate critical heights ||start||

def calculateCriticalHeights(tower):
    xvl = []
    yvl = []
    xvr = []
    yvr = []
    rch = []
    lch = []
    for i in range(len(xOfVertices)):
        if xOfVertices[i] <= tower:
            xvl.append(xOfVertices[i])
            yvl.append(yOfVertices[i])
        if xOfVertices[i] >= tower:
            xvr.append(xOfVertices[i])
            yvr.append(yOfVertices[i])
    xvl.reverse()
    yvl.reverse()
    if len(xvr)>3:
        rch = calculateCriticalPairs(copy.deepcopy(xvr),copy.deepcopy(yvr),1)
    if len(xvl)>3:
        lch = calculateCriticalPairs(copy.deepcopy(xvl),copy.deepcopy(yvl),-1)
    i=0
    while i < len(rch):
        flag = 0
        for j in range(len(xvr)):
            
            if xvr[j]> rch[i][1][0]:
                break
            if alg.intersectionPoint(rch[i],[[xvr[j],yvr[j]],[xvr[j],yvr[j]+1]])[1] < yvr[j]-0.00000001:
                flag =1
                z=rch.pop(i)
                break
        if flag == 0:
            i+=1
    i=0
    while i < len(lch):
        flag = 0
        for j in range(len(xvl)):
            
            if xvl[j]< lch[i][1][0]:
                break
            if alg.intersectionPoint(lch[i],[[xvl[j],yvl[j]],[xvl[j],yvl[j]+1]])[1] < yvl[j]-0.00000001:
                flag =1
                z=lch.pop(i)
                break
        if flag == 0:
            i+=1
    return rch+lch

def calculateCriticalPairs(xv,yv,t):    
    yVerticeStart = yv[0]
    xVerticeStart = xv[0]
    towerLine = [[xv[0],yVerticeStart],[xv[0],yVerticeStart+2]]
    criticalPairs = []
    convexHull = []
    convexHull.append([xv.pop(0),yv.pop(0)])
    convexHull.append([xv.pop(0),yv.pop(0)])
    while len(xv)>0:
        turn = calculateTurn(convexHull[-1], convexHull[-2],[xv[0],yv[0]])
        if turn == t:
            convexHull.append([xv.pop(0),yv.pop(0)])
        else:
            while len(convexHull)>1:
                if alg.intersectionPoint(towerLine,[convexHull[-2],[xv[0],yv[0]]])[1]>yVerticeStart:
                    criticalPairs.append([convexHull[-2],[xv[0],yv[0]]])
                    
                turn = calculateTurn(convexHull[-1], convexHull[-2],[xv[0],yv[0]])
                if turn == t*-1:
                    convexHull.pop(-1)
                else:
                    break
            convexHull.append([xv.pop(0),yv.pop(0)])
    return criticalPairs
            
def calculateTurn(v1, v2, v3):
    result = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0])
    if result >= 0:
        return 1 
    else:
        return -1
    
#------------Calculate critical heights ||end||

#------------Start The Procedures

def visualizeSemicontinious():
    towers = calculateSecondTowerSemiContiniousSolo()
    tower1 = towers[0]
    tower2 = towers[1]
    upperEnvelope = towers[2]
    plt.figure(num='The two watchtowers Semi')
    plt.title("Semi-Continious Problem")
    plt.plot(xOfVertices,yOfVertices,color = 'green')
    plt.plot(xOfVertices, yOfVertices, marker = ".", color = "green")
    if len(upperEnvelope)>0:
        plt.plot(upperEnvelope[0],upperEnvelope[1],color = 'blue')
        plt.plot(upperEnvelope[0], upperEnvelope[1], marker = ".", color = "blue")
    
    plt.plot([tower1[0],tower1[0]],[tower1[1],tower1[1]+tower1[2]],color = '0')
    plt.plot([tower1[0]],[tower1[1]+tower1[2]], marker = "s",color = "red")
    plt.plot([tower2[0],tower2[0]],[tower2[1],tower2[1]+tower2[2]],color = '0')
    plt.plot([tower2[0]],[tower2[1]+tower2[2]], marker = "s",color = "blue")
    plt.show()

def visualizeDescrite():
    towers = calculateSecondTowerDescriteSolo()
    tower1 = towers[0]
    tower2 = towers[1] 
    index1 = xOfVertices.index(tower1[0]) 
    index2 = xOfVertices.index(tower2[0]) 
    tower1.insert(1,yOfVertices[index1])
    tower2.insert(1,yOfVertices[index2])
    plt.figure(num='The two watchtowers Descrite')
    plt.title("Discrete Problem")
    plt.plot(xOfVertices,yOfVertices,color = 'green')
    plt.plot(xOfVertices, yOfVertices, marker = ".", color = "green")
    plt.plot([tower1[0],tower1[0]],[tower1[1],tower1[1]+tower1[2]],color = '0')
    plt.plot([tower1[0]],[tower1[1]+tower1[2]], marker = "s",color = 'red')
    plt.plot([tower2[0],tower2[0]],[tower2[1],tower2[1]+tower2[2]],color = '0')
    plt.plot([tower2[0]],[tower2[1]+tower2[2]], marker = "s",color = 'blue')
    y = min(yOfVertices)
    autoUpdateTower(tower1[0],tower1[1]+tower1[2])
    autoCalculateVisibilityLines()
    for hiddenPiece in setNotVisibleEdges:
        plt.plot([hiddenPiece[0][0],hiddenPiece[1][0]],[hiddenPiece[0][1],hiddenPiece[1][1]],linestyle = 'dotted', color = 'red')
    for endPoint in setVisP:   
        plt.plot([xTower, endPoint[0]],[yTower,endPoint[1]], linestyle = 'dotted',color = 'red')
    autoUpdateTower(tower2[0],tower2[1]+tower2[2])
    autoCalculateVisibilityLines()
    for hiddenPiece in setNotVisibleEdges:
        plt.plot([hiddenPiece[0][0],hiddenPiece[1][0]],[hiddenPiece[0][1],hiddenPiece[1][1]],linestyle = 'dotted', color = 'blue')
    for endPoint in setVisP:   
        plt.plot([xTower, endPoint[0]],[yTower,endPoint[1]], linestyle = 'dotted',color = 'blue')
    plt.show()



#------------Isolated Visualizations

    
def calculateSecondTower():
    autoCalculateVisibilityLines()
    visibilityPairsLeft = autoCalculateLeftPUV() #O(n^2logn)
    visibilityPairsRight = autoCalculateRightPUV() #O(n^2logn)
    setEv = visibilityPairsLeft[0]+visibilityPairsLeft[1]+visibilityPairsRight[0]+visibilityPairsLeft[1]
    currentTower2 = 0
    currentHeight2 = sys.float_info.max 
    for j in range(len(xOfVertices)):
        
        secondTowerLine = [[xOfVertices[j],yOfVertices[j]],[xOfVertices[j],yOfVertices[j]+2]]
        possibleHeights2 = []
        for line in setEv:
            if line[0][0] != line[1][0]:
                temp = alg.intersectionPoint(line ,secondTowerLine)[1]
                if round(temp,10) >= yOfVertices[j]:
                    possibleHeights2.append(round(temp,10))
        if len(possibleHeights2) > 0:
            currentY2 = max(possibleHeights2)
            tempHeight2 = currentY2 - yOfVertices[j]
        else:
            tempHeight2 = sys.float_info.max 
        if tempHeight2 < currentHeight2:
            currentHeight2 = tempHeight2
            currentTower2 = j
    if currentHeight2 == sys.float_info.max:
        currentHeight2 = 0
    drawTerrain()
    plt.plot([xOfVertices[currentTower2],xOfVertices[currentTower2]],[yOfVertices[currentTower2],yOfVertices[currentTower2]+currentHeight2],color = '0')
    plt.plot([xOfVertices[currentTower2]],[yOfVertices[currentTower2]+currentHeight2], marker = "s")
    plt.show()
    print(xOfVertices[currentTower2],currentHeight2)

def visualizeUpperEnvelope():
    autoCalculateVisibilityLines() #O(n)
    visibilityPairsLeft = autoCalculateLeftPUV() #O(n^2logn)
    visibilityPairsRight = autoCalculateRightPUV() #O(n^2logn)
    setB = visibilityPairsLeft[0]
    setC = visibilityPairsRight[0]
    setNotVisE = visibilityPairsRight[1]
    setLines = setB+setC+setNotVisE
    upperEnvelope = calculateUpperEnvelope(setLines)
    for hiddenPiece in setNotVisibleEdges:
        plt.plot([hiddenPiece[0][0],hiddenPiece[1][0]],[hiddenPiece[0][1],hiddenPiece[1][1]], linewidth=2,linestyle = 'dashed', color = '0')
    plt.plot(upperEnvelope[0],upperEnvelope[1],color = 'blue')
    plt.plot(upperEnvelope[0], upperEnvelope[1], marker = ".", color = "blue")
    plt.show()

def drawVisibilityLines():
    global setP, terrainWithP, setNotVisibleEdges 
    setNotVisibleEdges = []

    plt.close()
    drawTerrain()
    endPoints = alg.visibilityEndpoints(xOfVertices, yOfVertices, xTower, yTower) #returns [[endpoints- (x,y)], [endpoints and vertices out of sight- (x,y)]]
    setPonVertice= []
    for endPoint in endPoints[0]:
        if endPoint[0] in xOfVertices:
            setPonVertice.append(endPoint)
          
    setP = endPoints[0]
    for i in range(len(setP)):
        setP[i][0] = round(setP[i][0],13)
        setP[i][1] = round(setP[i][1],13)
    terrainWithP = [xOfVertices.copy(),yOfVertices.copy()]
    setPwithNoVertice = copy.deepcopy(setP)
    epCounter = 0
    if len(setPwithNoVertice)>0:
        for i in range(len(xOfVertices)):
            if len(setPwithNoVertice)>0:
                if epCounter >= len(setPwithNoVertice):
                        break     
                if  xOfVertices[i] == setPwithNoVertice[epCounter][0]:
                    setPwithNoVertice.remove(setPwithNoVertice[epCounter])
                    continue
                if  xOfVertices[i] > setPwithNoVertice[epCounter][0]:
                    terrainWithP[0].insert(epCounter+i,setPwithNoVertice[epCounter][0])
                    terrainWithP[1].insert(epCounter+i,setPwithNoVertice[epCounter][1])
                    epCounter+=1
                      
        for endPoint in setP+setPonVertice:    
            plt.plot([endPoint[0]],[endPoint[1]], marker = ".")
            plt.plot([xTower, endPoint[0]],[yTower,endPoint[1]], linestyle = 'dotted',color = '0')

   
    for p in setPonVertice:
        endPoints[1].remove(p)
    for i in range (len(endPoints[1])-1):
        if endPoints[1][i] in endPoints[2] and endPoints[1][i+1]in endPoints[2]:
            continue
        if endPoints[1][i][0] in xOfVertices and endPoints[1][i+1][0] in xOfVertices:
            start = xOfVertices.index(endPoints[1][i][0])
            finish = xOfVertices.index(endPoints[1][i+1][0])
            if abs(start - finish) == 1:
                setNotVisibleEdges.append([endPoints[1][i],endPoints[1][i+1]])
        if endPoints[1][i][0] in xOfVertices and endPoints[1][i+1] in setP and endPoints[1][i+1] not in setPonVertice:
            setNotVisibleEdges.append([endPoints[1][i],endPoints[1][i+1]])
    setPnotValuable = copy.deepcopy(setPonVertice)
    
    for point in setPnotValuable:
        counter = 0
        for pair in setNotVisibleEdges:
            for p in pair:
                if p == point:
                    counter+=1
        if counter==1:
            setPnotValuable.remove(point)
    for p in setP:
        if p in setPnotValuable:
            setP.remove(p)


    for hiddenPiece in setNotVisibleEdges:
        for v in hiddenPiece:
            if v not in setP and v not in endPoints[2]:
                setP.append(v)
    setP = sorted(setP, key=lambda p:p[0])

    for hiddenPiece in setNotVisibleEdges:
        plt.plot([hiddenPiece[0][0],hiddenPiece[1][0]],[hiddenPiece[0][1],hiddenPiece[1][1]], linewidth=2,linestyle = 'dashed', color = '0')
    plt.show()


def calculateLeftPUV():
    setE = []
    puv = dc.constructPUV(xOfVertices, terrainWithP,setP)
    for part in setNotVisibleEdges:
        setE.append(part)
    leftPV = [puv,setE]
    for pair in leftPV[0]:
        copyPair = copy.deepcopy(pair)

        rePair = [copyPair[1],copyPair[0]]
        if pair[1][0] not in xOfVertices and  rePair in leftPV[0]:
            leftPV[0].remove(pair)
    
    plt.figure()
    plt.plot(xOfVertices,yOfVertices,color = 'green')
    plt.plot(xOfVertices, yOfVertices, marker = ".", color = "green")
    for endPoint in setP:    
            plt.plot([endPoint[0]],[endPoint[1]], marker = ".", color = '0')
    if numberTower != -1:
        plt.plot([xTower,xTower],[yBase,yTower],color = '0')
        plt.plot([xTower],[yTower], marker = "s")
    for pair in leftPV[0]:
        plt.plot([pair[0][0], pair[1][0]],[pair[0][1],pair[1][1]], linestyle = 'dotted',color = '0')
    
   
    plt.show()

def calculateRightPUV():
    setE = []
    rightTerrain = copy.deepcopy(terrainWithP)
    rightSetP = copy.deepcopy(setP)
    rightTerrain[0] = [rightTerrain[0][-1] - x for x in rightTerrain[0]]
    
    for i in range(len(rightSetP)):
        rightSetP[i][0] = terrainWithP[0][-1] - rightSetP[i][0]

    rightTerrain[0].reverse()
    rightTerrain[1].reverse()
    rightSetP.reverse()
   
    revXofVertices = copy.deepcopy(xOfVertices)
    revXofVertices = [xOfVertices[-1] - x for x in xOfVertices]
    rightPUV = dc.constructPUV(revXofVertices, rightTerrain,rightSetP)
    plt.figure()
    for i in range(len(rightPUV)):
        rightPUV[i][0][0] = xOfVertices[-1]-rightPUV[i][0][0]
        rightPUV[i][1][0] = xOfVertices[-1]-rightPUV[i][1][0]
    for part in setNotVisibleEdges:
        setE.append(part)
    rightPV = [rightPUV,setE]
    plt.plot(xOfVertices,yOfVertices,color = 'green')
    plt.plot(xOfVertices, yOfVertices, marker = ".", color = "green")
    if numberTower != -1:
        plt.plot([xTower,xTower],[yBase,yTower],color = '0')
        plt.plot([xTower],[yTower], marker = "s")
    for pair in rightPV[0]:
        plt.plot([pair[0][0], pair[1][0]],[pair[0][1],pair[1][1]], linestyle = 'dotted',color = '0')
    
    plt.show()

def disable_event():
    pass

def initiate():
    global myTerrains, winLoadTerrain, terrainPath, winVertices, listOfVertices, text_boxHeight
    terrainPath = os.path.dirname(os.path.abspath(__file__))
    terrains = os.listdir(terrainPath + "\\terrains")

    window = Tk()
    window.title("Two Watchtowers")
    window.geometry("500x400+1000+300")

    # Top level for loading a terrain

    winLoadTerrain = Toplevel()   
    winLoadTerrain.title("Load Terrain")
    winLoadTerrain.geometry("200x210+800+300")
    winLoadTerrain.withdraw()
    winLoadTerrain.protocol("WM_DELETE_WINDOW", disable_event)

    buttonLoadTerrain = Button(winLoadTerrain, text= "Load Terrain", command= loadTerrain)
    buttonLoadTerrain.place(x=10,y=185)

    sbAvailableVertices = Scrollbar(winLoadTerrain)  
    sbAvailableVertices.pack(side = RIGHT, fill = Y)  
    myTerrains = Listbox(winLoadTerrain, yscrollcommand = sbAvailableVertices.set) 
    myTerrains.pack(side = LEFT)
    sbAvailableVertices.config(command = myTerrains.yview)
    for script in terrains:
        myTerrains.insert(END, script[:-4])


    # Top level for choosing vertex for the first tower

    winVertices = Toplevel()   
    winVertices.title("Choose Vertice")
    winVertices.geometry("150x210+850+300")
    winVertices.withdraw()
    winVertices.protocol("WM_DELETE_WINDOW", disable_event)


    buttonOk = Button(winVertices, text= "Ok", command= updateTower)
    buttonOk.place(x=60,y=185)

    sbAvailableVertices = Scrollbar(winVertices)  
    sbAvailableVertices.pack(side = RIGHT, fill = Y)  
    listOfVertices = Listbox(winVertices, yscrollcommand = sbAvailableVertices.set) 
    listOfVertices.pack(side = LEFT)
    sbAvailableVertices.config(command = listOfVertices.yview)

    text_boxHeight = Text(
    winVertices,
    height=1,
    width=15
    )
    text_boxHeight.place(x=0,y=0)
    text_boxHeight.insert('end', message1)
    text_boxHeight.config(state='normal')

    #Buttons on main window

    buttonChooseTerrain = Button(window, bg='green', text= "Choose Terrain", command= chooseTerrain)
    buttonChooseVertice = Button(window, text= "Choose Vertex", command= chooseVertice)
    buttonCurrentTerrain = Button(window, text= "Draw Vis Lines", command= drawVisibilityLines)

    buttonLeftPUV = Button(window, text= "leftPUV", command= calculateLeftPUV)
    buttonRightPUV = Button(window, text= "rightPUV", command= calculateRightPUV)
    buttonDiscrete = Button(window, bg='green', text= "Visualize Discrete", command= visualizeDescrite)
    buttonCalculateManualSecond = Button(window, text= "Calculate manual Tower", command= calculateSecondTower)
    buttonCalculateUpperEnvelope = Button(window, text= "Visualize Upper Envelope With All Pairs", command= visualizeUpperEnvelope)
    buttonsVisualizeSemiCοntinious = Button(window, bg='green', text= "Visualize Semi-Continious", command= visualizeSemicontinious)

    buttonChooseTerrain.place(x = 10, y = 10)
    buttonDiscrete.place(x = 10, y = 40)
    buttonsVisualizeSemiCοntinious.place(x = 10 , y = 70)
    buttonChooseVertice.place(x=10,y=100)
    buttonCurrentTerrain.place(x = 10, y = 130)    
    buttonLeftPUV.place(x = 10,y = 160)
    buttonRightPUV.place(x = 10, y = 190)
    buttonCalculateManualSecond.place(x=10,  y=220)
    buttonCalculateUpperEnvelope.place(x=10, y=250)
    print("hello etiko")
    window.mainloop()

initiate()