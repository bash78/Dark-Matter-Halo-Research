import numpy
import astropy
import IPython
import math
import random
#time interval is in seconds it is the time between each force interaction, so making it smaller would increase precision, but increases processing time
timeInterval = 10
#number of matter points in a galaxy, so making it larger would increase precision, but increases processing time
matterPerGalaxy = 100000
#The number of times a time interval passes after the start of a galaxy, increasing this number incrases precision, but increases processing time
timeScale = 100000000
#radius of the galaxy, assunming it is uniform like a sphere
galaxyScale = 87.7
#the number of potential galaxies
galaxyNumber = 1000
#the amount of dark matter points in a galaxy
darkMatterPerGalaxy = 300000
#the largest percentage difference for any calculation for its expected (for matterPoints)
fudgeFactor = 0.01
#the largest percentage difference for any calculation for its expected (for darkMatterPoints)
darkFudgeFactor = 0.1
#greatest speed at which a dark matter halo and collapse
collapseSpeed = 10
#How far off a radius has to be to not be accepted in the average
radiusError = 20
#The standard deviation for a spirial galaxy
sigmaY = 2
#The converstion rate form an arcsec to meters
arc2meters = 206264.8062937
#converstion rate from kilometers to meters
kilo2meters = 1000
#The standard devitation of the matter, basically meaning that the mass of each matter point will be slightly skewed from 1 kg
sigmaM = 5000
#The average mas sof one matter point
averageMatter = 10000
#The class that assempbles the data for each galaxy into arrays to be anaylzed later, it acts as the properties of one of the known galxies from the data table
class galaxyData:
    def __init__ (self):
        self.radius = [0.37 for i in range(1000)]
        self.velocity = [0.37 for i in range(1000)]
        self.error = [0.37 for i in range(1000)]
        self.datar = ""
        self.name = ""
        self.turnOn = False
    #Where the class assembles the numbers from a given string into the properties of this known glaxy
    def sorter (self):
        index1 = 0
        index2 = 0
        char = self.datar.split()
        for i in char:
            if index1 == 0:
                self.radius[index2] = float(i)
            elif index1 == 1:
                self.velocity[index2] = float(i)
            else:
                self.error[index2] = float(i)
            index1+=1
            if index1 >= 3:
                index1 = 0
                index2+=1
    #A function used to refrence values
    def interpret (self, idString, idNumber):
        if idString == "radius":
            return self.radius[idNumber]
        elif idString == "velocity":
            return self.velocity[idNumber]
        elif idString == "error":
            return self.error[idNumber] 
#From here until the next comment, the program breaks apart the whole data table into string of data for each galaxy
datas = None
data = [galaxyData() for i in range(100)]
radiusLiklyhood = [0.0 for i in range(len(data) * galaxyNumber)]
radiusIndex1 = [0.0 for i in range(len(radiusLiklyhood))]
radiusIndex2 = [0.0 for i in range(len(radiusLiklyhood))]
with open('/Users/benjaminash/Desktop/School/Summer Projects/Dark Matter Halo Project/RadialVelocity.rtf', 'r') as file:
    datas = file.read()
sort = datas.split()
gindex = -1
rindex = 0
for point in sort:
    typer = 0
    try:
        typer = float(point)
    except:
        if point == "0":
            typer = 0
        else:
            typer = 0.37
    if typer == 0.37:
        gindex+=1
        data[gindex].turnOn = True
        data[gindex].name = point
    else: 
        data[gindex].datar = data[gindex].datar + str(typer) + " "
for gal in data:
    if gal.turnOn:
        gal.sorter()
print("Done sorting for data sheet")
#This is the class that contains the properties and the relative actions of a point of matter
class matterPoint:
    def __init__ (self):
        self.positionx = 0
        self.positiony = 0
        self.positionz = 0
        self.dister = galaxyScale * 2
        self.indexer = 1000
        self.matter = random.gauss(averageMatter, sigmaM)
        self.velocityx = 0
        self.velocityy = 0
        self.velocityz = 0
    #Does as the name entails, gives the point a random position for the start of the program, where it is more constained on the y-axis
    def positionRandomizer (self):
        self.positionx = random.uniform(-1 * galaxyScale, galaxyScale)
        self.positiony = random.gauss(0, sigmaY)
        self.positionz = random.uniform(-1 * galaxyScale, galaxyScale)
    #Finds the distance between this matter point and another matter point
    def distance (self, point):
        deltax = abs(self.positionx - point.positionx)
        deltay = abs(self.positiony - point.positiony)
        deltaz = abs(self.positionz - point.positionz)
        deltaxy = math.sqrt(math.pow(deltax, 2) * math.pow(deltay, 2))
        dist = math.sqrt(pow(deltaz, 2) * pow(deltaxy, 2))
        return dist
    #Finds an angle between this matter point and another matter point
    def angle (self, point, type):
        deltax = abs(self.positionx - point.positionx)
        deltay = abs(self.positiony - point.positiony)
        deltaz = abs(self.positionz - point.positionz)
        deltaxy = math.sqrt(math.pow(deltax, 2) * math.pow(deltay, 2))
        anglex = 0
        anglez = 0
        try:
            anglex = math.tan(deltay/deltax)
        except:
            anglex = 0
        try:
            anglez = math.tan(deltaz/deltaxy)
        except:
            anglez = 0
        if type:
            return anglex
        else:
            return anglez
    #Is the little n-body simulator for this matter point, and looks to see all of the forces acting on this point of matter
    def gravity (self, points, indexg):
        fx = 0.000000000001
        fy = 0.000000000001
        fz = 0.000000000001
        for point in points:
            if point.matter > 0:
                if point.positionx != self.positionx:
                    fx =  fx + 6.67408 * math.pow(10, -11) * self.matter * point.matter / math.pow((point.positionx * arc2meters - self.positionx * arc2meters), 2)
                if point.positiony != self.positiony:
                    fy =  fy + 6.67408 * math.pow(10, -11) * self.matter * point.matter / math.pow((point.positiony * arc2meters  - self.positiony * arc2meters), 2)
                if point.positionz != self.positionz:
                    fz =  fz + 6.67408 * math.pow(10, -11) * self.matter * point.matter / math.pow((point.positionz * arc2meters - self.positionz * arc2meters), 2)
        vx = (fx * timeInterval / self.matter) / kilo2meters
        vy = (fy * timeInterval / self.matter) / kilo2meters
        vz = (fz * timeInterval / self.matter) / kilo2meters
        velocities = [vx, vy, vz]
        return velocities[indexg]
    #Finds the distncae from this matter point and the center of mass of the gravity
    def radiusFinder (self, points):
        totalMatter = 0
        averageX = 0.0000000001
        averageY = 0.0000000001
        averageZ = 0.0000000001
        totalPoints = 0
        for point in points:
            if point.matter > 0:
                totalMatter = totalMatter + point.matter
                totalPoints+=1
                averageX = averageX + point.positionx
                averageY = averageY + point.positiony
                averageZ = averageZ + point.positionz
        averageX = averageX / totalPoints
        averageY = averageY / totalPoints
        averageZ = averageZ / totalPoints
        distX = self.positionx - averageX
        distY = self.positiony - averageY
        distZ = self.positionz - averageZ
        return math.sqrt(math.pow(distX, 2) + math.pow(distY, 2) + math.pow(distZ, 2))
    #Finds the total speed, using the components
    def totalSpeed (self):
        total = math.sqrt(math.pow(self.velocityx, 2))
        return total
    #Finds the x angle, for the direction of the velocity
    def xAngle(self):
        angle = math.tan(self.velocityx / self.velocityz)
        return angle
    #Finds the z angle, for the direction of the velocity
    def zAngle(self):
        sidexz = abs(math.pow(self.velocityx, 2) + math.pow(self.velocityz, 2))
        angle = math.tan(self.velocityy / sidexz)
        return angle
    #Finds the starting velocity of this matter point by looking at centripetal force equaling the gravitational force, and solving for v
    def startGravity (self, points, indexsg, darkMatter):
        totalMatter = 0
        radius = 0.000000001
        averageX = 0.0000000001
        averageY = 0.0000000001
        averageZ = 0.0000000001
        fudgeX = random.uniform(-1 * fudgeFactor, fudgeFactor)
        fudgeY = random.uniform(-1 * fudgeFactor, fudgeFactor)
        fudgeZ = random.uniform(-1 * fudgeFactor, fudgeFactor)
        totalPoints = 0
        for point in points:
            if point.matter > 0:
                totalMatter = totalMatter + point.matter
                totalPoints+=1
                averageX = averageX + point.positionx
                averageY = averageY + point.positiony
                averageZ = averageZ + point.positionz
        averageX = averageX / totalPoints
        averageY = averageY / totalPoints
        averageZ = averageZ / totalPoints
        radius = self.radiusFinder(points)
        ratioX = math.pow(averageX, 2) / math.pow(radius, 2)
        ratioY = math.pow(averageY, 2) / math.pow(radius, 2)
        ratioZ = math.pow(averageZ, 2) / math.pow(radius, 2)
        if darkMatter == False:
            totalVelocity = math.sqrt((6.67408 * math.pow(10, -11) * totalMatter) / (radius * arc2meters)) / kilo2meters
        else:
            totalVelocity = random.uniform(-1 * collapseSpeed, collapseSpeed)
            fudgeX = random.uniform(-1 * darkFudgeFactor, darkFudgeFactor)
            fudgeY = random.uniform(-1 * darkFudgeFactor, darkFudgeFactor)
            fudgeZ = random.uniform(-1 * darkFudgeFactor, darkFudgeFactor)
        vx = totalVelocity * ratioX
        vx = vx - (vx * fudgeX)
        vy = totalVelocity * ratioY
        vy = vy - (vy * fudgeY)
        vz = totalVelocity * ratioZ
        vz = vz - (vz * fudgeZ)
        velocities = [vx, vy, vz]
        return velocities[indexsg]
#Contains the properties of a computer generated galaxy, and basically acts as a storehouse for all of the galaxy's matter
class galaxy:
    def __init__ (self):
        self.totalMatterPoints = [matterPoint() for i in range(matterPerGalaxy + darkMatterPerGalaxy)]
        self.matterPoints = [matterPoint() for i in range(matterPerGalaxy)]
        self.darkMatterPoints = [matterPoint() for i in range(darkMatterPerGalaxy)]
    #Moves a matter point
    def poisitionChanger (self, velocity, position):
        positionFind = ((velocity * kilo2meters) / arc2meters) * timeInterval + position
        return positionFind
    #Moves teh matter points over the given set of time, in essence seeing how the properties of the galaxy change over time to have better results to look at
    def mover (self):
        for matter in self.matterPoints:
            matter.velocityx = matter.startGravity(self.totalMatterPoints, 0, False)
            matter.velocityy = matter.startGravity(self.totalMatterPoints, 1, False)
            matter.velocityz = matter.startGravity(self.totalMatterPoints, 2, False)
        for matter in self.darkMatterPoints:
            matter.velocityx = matter.startGravity(self.totalMatterPoints, 0, True)
            matter.velocityy = matter.startGravity(self.totalMatterPoints, 1, True)
            matter.velocityz = matter.startGravity(self.totalMatterPoints, 2, True)
        for matter in self.totalMatterPoints:
            timer = 0
            while timer < timeScale:
                matter.velocityx = matter.velocityx + matter.gravity(self.matterPoints, 0)
                matter.positionx = self.poisitionChanger(matter.velocityx, matter.positionx)
                matter.velocityy = matter.velocityy + matter.gravity(self.matterPoints, 1)
                matter.positiony = self.poisitionChanger(matter.velocityy, matter.positiony)
                matter.velocityz = matter.velocityz + matter.gravity(self.matterPoints, 2)
                matter.positionz = self.poisitionChanger(matter.velocityz, matter.positionz)
                timer+=1
    #Looks to sort the matter point velocitys to later be analyzed in comparision to the known velocity by trying to group velocities based on that matter point's distance from the center of mass
    def radiusAnalyzer (self, galNumber, radiusArray, galler, gallerData):
        velocities = [0.0 for i in range(len(radiusArray))]
        numberOfMatter = [0.1 for i in range(len(radiusArray))]
        for radius in range(len(radiusArray)):
            if radiusArray[radius] != 0.37:
                for matter in range(len(self.matterPoints)):
                    difference = abs(radiusArray[radius] - self.matterPoints[matter].radiusFinder(self.matterPoints))
                    if difference <= radiusError:
                        velocities[radius] = velocities[radius] + (self.matterPoints[matter].velocityx / 1000)
                        if numberOfMatter[radius] == 0.1:
                            numberOfMatter[radius] = 0
                        numberOfMatter[radius] = numberOfMatter[radius] + 1
        for i in range(len(velocities)):
            if velocities[i] != 0:
                velocities[i] = velocities[i] / numberOfMatter[i]
        return velocities[galNumber]
#Makes sure all of the randomlly generated generated galaxy strucures are different, so the computer is not comparing the same thing twice
galaxies = [galaxy() for i  in range(galaxyNumber)]
defaultMatter = [matterPoint() for i in range((matterPerGalaxy + darkMatterPerGalaxy) * galaxyNumber)]
approve = [False for i in range(galaxyNumber)]
approve[0] = True
galIndex = 0
matterIndex = 0
for i in range ((matterPerGalaxy + darkMatterPerGalaxy) * galaxyNumber):
    defaultMatter[i].positionRandomizer()
    galaxies[galIndex].totalMatterPoints[matterIndex] = defaultMatter[i]
    matterIndex+=1
    if matterIndex == matterPerGalaxy + darkMatterPerGalaxy:
        matterIndex = 0
        galIndex+=1
for gal in range(galaxyNumber):
    while approve[gal] == False:
        for i in range(gal):
            indexs = [1000 for j in range(matterPerGalaxy + darkMatterPerGalaxy)]
            extraMaterPoints = [matterPoint() for i in range(matterPerGalaxy + darkMatterPerGalaxy)]
            for m in range (matterPerGalaxy + darkMatterPerGalaxy):
                extraMaterPoints[m] = galaxies[i].totalMatterPoints[m]
            indexing = False
            while indexing == False:
                for n in range(matterPerGalaxy + darkMatterPerGalaxy):
                    for j in range(len(extraMaterPoints)):
                        if galaxies[gal].totalMatterPoints[n].distance(extraMaterPoints[j]) < galaxies[gal].totalMatterPoints[n].dister and galaxies[gal].totalMatterPoints[n].indexer == 1000:
                            galaxies[gal].totalMatterPoints[n].dister = galaxies[gal].totalMatterPoints[n].distance(extraMaterPoints[j])
                            galaxies[gal].totalMatterPoints[n].indexer = j
                    indexs[n] = galaxies[gal].totalMatterPoints[n].indexer                                                                         
                for index1 in range(len(extraMaterPoints)):
                    for index2 in range(len(extraMaterPoints)):
                        if indexs[index1] == indexs[index2] and index1 != 1000 and index2 != 1000:
                            if galaxies[gal].totalMatterPoints[index2].dister < galaxies[gal].totalMatterPoints[index1].dister:
                                indexs[index1] = 1000
                                galaxies[gal].totalMatterPoints[index1].indexer = 1000
                                galaxies[gal].totalMatterPoints[index1].dister = galaxyScale * 2
                for extras in indexs:
                    if extras != 1000 and extras < len(extraMaterPoints):
                        extraMaterPoints.remove(extraMaterPoints[extras])
                if len(extraMaterPoints) == 0:
                    indexing = True
            averageDistance = 0
            for matter in galaxies[gal].totalMatterPoints:
                averageDistance = averageDistance + matter.dister
            averageDistance = averageDistance / (matterPerGalaxy + darkMatterPerGalaxy)
            meanDiviationPercent = averageDistance / (galaxyScale * 2)
            if meanDiviationPercent >= 0.05:
                approve[gal] = True
            else:
                extraMatter = [matterPoint() for i in range(matterPerGalaxy * darkMatterPerGalaxy)]
                for matter in range(len(extraMatter)):
                    extraMatter[matter].positionRandomizer()
                    galaxies[gal].totalMatterPoints[matter] = extraMatter[matter]
print("Done making sure all of the galaxies are different")
#Now that the computer is sure that the galxies matter points are in different locations, it is telling that galxy class what its matter points will be, and sets the galaxy in motion for a given period of time
for galaxy in galaxies:
    for i in range(len(galaxy.totalMatterPoints)):
        if i < matterPerGalaxy:
            galaxy.matterPoints[i] = galaxy.totalMatterPoints[i]
        else:
            galaxy.darkMatterPoints[i - matterPerGalaxy] = galaxy.totalMatterPoints[i]
    galaxy.mover()
print("Done simulating the galxies through time")
#Compares the galxies from the data table to the computer generted galaxies
liklyIndex = 0
for gal in range(len(galaxies)):
    for galData in range(len(data)):
        sum = 0
        read = ""
        if data[galData].turnOn:
            for rad in range(len(data[galData].radius)):
                if data[galData].radius[rad] != 0.37:
                    collective = galaxies[gal].radiusAnalyzer(rad, data[galData].radius, gal, galData)
                    if collective == 0.0:
                        sum = 100000000000000000000
                    else:
                        total = math.pow(abs(collective) - abs(data[galData].radius[rad]), 2)
                        sum = total / math.pow(data[galData].error[rad], 2)
        if sum != 0 and sum != 100000000000000000000:
            radiusLiklyhood[liklyIndex] = sum
            radiusIndex1[liklyIndex]= gal
            radiusIndex2[liklyIndex] = galData
            liklyIndex+=1
print("Done summing")
#writing out area
writeoutSums = ""
for output in range(len(radiusLiklyhood)):
    if radiusLiklyhood[output] != 0.0:
        writeoutSums = writeoutSums + str(radiusLiklyhood[output]) + " " + str(radiusIndex1[output]) + " " + str(radiusIndex2[output])
writingOut = ""
for i in range(galaxyNumber):
    for j in range(matterPerGalaxy):
        writingOut = writingOut + " " + str(galaxies[i].matterPoints[j].positionx)
        writingOut = writingOut + " " + str(galaxies[i].matterPoints[j].positiony)
        writingOut = writingOut + " " + str(galaxies[i].matterPoints[j].positionz)
with open('/Users/benjaminash/Desktop/School/Summer Projects/Dark Matter Halo Project/Positions.txt', 'w') as file:
    file.write(writingOut)
with open ('/Users/benjaminash/Desktop/School/Summer Projects/Dark Matter Halo Project/PossibleSums.rtf', 'w') as file:
    file.write (writeoutSums)
print("Done printing")
print("Done running program")