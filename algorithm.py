

import sys
import matplotlib.pyplot as plt



#gets the terrain and the tower and returns the point that are not visible (Endpoints of part visible edges and vertexes) O(n) time

def visibilityEndpoints(xOfVertices, yOfVertices, xTower, yTower):
    numberTower = xOfVertices.index(xTower)
    curSlop = 0
    minSlop = sys.float_info.max
    segment = []
    ray = []
    endPoints = []
    outOfSight = []
    inSight= []
    #calculate endPoints left of the first tower
    for i in range(numberTower-1,-1,-1):
    #because terrain is x monotonous I use the slope for finding the shady parts of the terrain
    # If a slope is lesser than its previous then there are unwatched edges or part of edge        
        curSlop = (yTower-yOfVertices[i])/(xTower-xOfVertices[i])
        if curSlop<=minSlop:
            minSlop = curSlop
            if i < numberTower-1:
                segment = [[xOfVertices[i],yOfVertices[i]],[xOfVertices[i+1],yOfVertices[i+1]]]
                endPoint = intersectionPoint(ray,segment)
                if endPoint == 0: #critically visible edge
                    inSight.append([xOfVertices[i],yOfVertices[i]])
                    continue
                endPoint[0] = round(endPoint[0],13)
                endPoint[1] = round(endPoint[1],13)
                if ray[0]!=segment[1]:
                    endPoints.insert(0,endPoint)
                    outOfSight.append(endPoint)
            ray = [[xOfVertices[i],yOfVertices[i]],[xTower,yTower]]
            if i>0:
                nextSlop = (yTower-yOfVertices[i-1])/(xTower-xOfVertices[i-1])
                if curSlop < nextSlop:
                    outOfSight.append([xOfVertices[i],yOfVertices[i]])
                    inSight.append([xOfVertices[i],yOfVertices[i]])
        else:
            if xOfVertices[i+1] not in outOfSight:
                outOfSight.append([xOfVertices[i+1],yOfVertices[i+1]])
            outOfSight.append([xOfVertices[i],yOfVertices[i]])
   
    curSlop = 0
    maxSlop = sys.float_info.max*-1   
    segment = []
    ray = []
    rightOutOfSight = []

    #calculate endPoints right of the first tower
    for i in range(numberTower+1,len(xOfVertices)):
        curSlop = (yTower-yOfVertices[i])/(xTower-xOfVertices[i])
    #because terrain is x monotonous I use the slope for finding the shady parts of the terrain
    #If a slope is lesser than its previous then there are unwatched edges or part of edge          
        if curSlop>=maxSlop: 
            maxSlop = curSlop
            if i > numberTower+1:
                
                segment = [[xOfVertices[i],yOfVertices[i]],[xOfVertices[i-1],yOfVertices[i-1]]]
                endPoint = intersectionPoint(ray,segment)
                if endPoint == 0: #critical visible vertice
                    inSight.append([xOfVertices[i],yOfVertices[i]])
                    continue
                endPoint[0] = round(endPoint[0],13)
                endPoint[1] = round(endPoint[1],13)
               
                if ray[0]!=segment[1]:
                    if xOfVertices[i-1] not in outOfSight:
                        outOfSight.append([xOfVertices[i-1],yOfVertices[i-1]])
                    endPoints.append(endPoint)
                    outOfSight.append(endPoint)
                    rightOutOfSight.append(endPoint)
                    inSight.append(endPoint)
            ray = [[xOfVertices[i],yOfVertices[i]],[xTower,yTower]] 
            if i<len(xOfVertices)-1:
                nextSlop = (yTower-yOfVertices[i+1])/(xTower-xOfVertices[i+1])
                if curSlop > nextSlop:
                    outOfSight.append([xOfVertices[i],yOfVertices[i]])
                    inSight.append([xOfVertices[i],yOfVertices[i]])
                    rightOutOfSight.append([xOfVertices[i],yOfVertices[i]])
        else:
            outOfSight.append([xOfVertices[i],yOfVertices[i]])
            
            rightOutOfSight.append([xOfVertices[i],yOfVertices[i]])

    return [endPoints, outOfSight,inSight]

#gets 2 lines and returns the intersection

def intersectionPoint(line1, line2): #line1 =  [[x1line1,y1line1],[x2line1,y2line1]], line2 same 
    xdiff = [line1[0][0] - line1[1][0], line2[0][0] - line2[1][0]]
    ydiff = [line1[0][1] - line1[1][1], line2[0][1] - line2[1][1]]

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       return 0

    d = [det(*line1), det(*line2)]
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return [x, y]

#recursevly calculating for rightmost vertice which p endpoint can see appends it to the set and calls the next vertice until PLVL is empty  
def pairPV(xOfVertices, subTerrain, setP, setPV):
    if len(subTerrain[0]) == 0:
        return setPV
    Vx = subTerrain[0].pop(-1)
    Vy = subTerrain[1].pop(-1)
    while Vx not in xOfVertices and len(subTerrain[0])>0:
        Vx = subTerrain[0].pop(-1)
        Vy = subTerrain[1].pop(-1)
    minSlop = sys.float_info.max 
    pair = []
    if len(subTerrain[0]) == 0:
        return setPV
    else:
        for i in range(len(subTerrain[0])):
            curSlop = (Vy-subTerrain[1][-i-1])/(Vx-subTerrain[0][-i-1])
            if curSlop<minSlop:
                minSlop = curSlop
                if [subTerrain[0][-i-1],subTerrain[1][-i-1]] in setP:
                    pair = [[subTerrain[0][-i-1],subTerrain[1][-i-1]],[Vx,Vy]]
                    setPV.append(pair)
                    break
        return pairPV(xOfVertices, subTerrain, setP, setPV)
    
#inputs: p points left of l vertical line , vertices left of l vertical line
#for every PL(p point left of l vertical line calculate vertice that allows PL to see l vertical line                  
def calculatePA(pL,vL):
    a = [[],[]]
    for i in range(len(pL[0])):
        while 1:

            if len(vL[0])==0:
                return a
            if vL[0][0]<= pL[0][i]:
                vL[0].pop(0)
                vL[1].pop(0)
            else:
                break
        maxSlop = (vL[1][0]-pL[1][i])/(vL[0][0]-pL[0][i])
        currA = [vL[0][0], vL[1][0]]
        for j in range(1,len(vL[0])):
            curSlop = (vL[1][j]-pL[1][i])/(vL[0][j]-pL[0][i])
            if curSlop > maxSlop:
                maxSlop = curSlop
                currA = [vL[0][j], vL[1][j]]
        
        a[0].append(currA[0])
        a[1].append(currA[1])
    return a     

#inputs: vertices right of l vertival line and VR' vertices
#for every VR' vertice calculate vertice that allows VR' to see l vertical line
def calculateVC(vR,notVR):
    c=[[],[]]
    for i in range(len(vR[0])-1,-1,-1):
        if vR[0][i] in notVR[0]:
            minSlope = (vR[1][i]-vR[1][i-1])/(vR[0][i]-vR[0][i-1])
            currC = [vR[0][i-1],vR[1][i-1]]
            for j in range(i-1,-1,-1):
                curSlope = (vR[1][j]-vR[1][i])/(vR[0][j]-vR[0][i])
                if curSlope < minSlope:
                    minSlope = curSlope
                    currC = [vR[0][j],vR[1][j]]
            c[0].insert(0,currC[0])
            c[1].insert(0,currC[1])
    return c


def calculateLineFrom2Points(a,b): #a = [ax,ay] and b = [bx,by]
    wpm = (a[1] - b[1])/(a[0]-b[0])
    wpn = a[1] - wpm * a[0]
    return [wpm,wpn]

