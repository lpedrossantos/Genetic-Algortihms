import math
import random
import numpy as np

"""
Travel Salesman Problem (simple version):
Constraints:
    -It must pass through every points/cities;
    -It can not pass more than one time at each point/city;
    -Needs to return to the initial point/city.
"""

#Compute the distance between points using the L2 norm
def calcDistance(genes,points):
    dist = 0
    for j in range(len(genes)-1):
        dist+= np.linalg.norm(points[genes[j]] - points[genes[j+1]])
    return dist

#Object DNA
class DNA:
    #Constructor
    def __init__(self,lenChromosome,startPoint):
        self.lenChromosome = lenChromosome
        self.startPoint = startPoint
        self.genesArr = []
        self.genesVal = []
        self.genesFitness = []
        self.genesProb = []
        self.selectedGenes = []
        self.childChromosomes = []
        self.totalFitness = None
    #Checks if chrosomes respect the constraints    
    def verify(self,genes):
        #First and 
        if genes[0] == self.startPoint and genes[-1] == self.startPoint and len(np.unique(genes[0:-1])) == (self.lenChromosome - 1):
            return True
        else:
            return False
    #Instatiate population  
    def populate(self,popSize,points):
        for i in range(0,popSize):
            while True:
                genes = [np.random.randint(0,len(points)) for i in range(0,dna.lenChromosome)]
                if self.verify(genes):
                    self.genesArr.append(genes)
                    break
    #Fitness function
    def fitness(self):
        self.genesFitness.clear()
        for i in range(len(self.genesVal)):
                fitness = 1 / (self.genesVal[i]+1)
                self.genesFitness.append(fitness)
        print("Genes Val:", self.genesVal)
        
        self.totalFitness = np.sum(self.genesFitness)
        print("Fitness:",self.genesFitness)
        print("Total:",self.totalFitness)
    
    def prob(self):
        self.genesProb = [val/self.totalFitness for val in self.genesFitness]
        print("Probabilities:",self.genesProb)
    #Selection
    def selection(self,n):
        self.selectedGenes.clear()
        rouletteWheel = np.cumsum(self.genesFitness)
        #Stochastic Universal Sampling
        distPointers = self.totalFitness/n
        start = np.random.uniform(0,distPointers)
        #Compute pointers
        pointerArr = [start+i*distPointers for i in range(0,n)]
        
        #Locate points
        for pointer in pointerArr:
            j=0
            while rouletteWheel[j] < pointer:
                j+=1
            self.selectedGenes.append(self.genesArr[j])
        print("Selected Genes:", self.selectedGenes)
    #Crossover
    def crossover(self,crossoverRate):
        crossoverParents = []
        self.childChromosomes.clear()
        #Select Parents
        for i in range(len(self.selectedGenes)):
            r = random.random()
            if r <= crossoverRate:
                crossoverParents.append(self.selectedGenes[i])
            else:
                self.childChromosomes.append(self.selectedGenes[i])
        
        if crossoverParents:
            #Perform crossover
            for i in range(len(crossoverParents)):
                #If it is the last element
                if (i+1) >= len(crossoverParents):
                    parent1 = crossoverParents[i][1:-1]
                    parent2 = crossoverParents[0][1:-1]
                else:
                    parent1 = crossoverParents[i][1:-1]
                    parent2 = crossoverParents[i+1][1:-1]
                crossoverPoint = np.random.randint(1,self.lenChromosome-1)
                child1 = [self.startPoint] + parent1[0:crossoverPoint] + parent2[crossoverPoint:self.lenChromosome] + [self.startPoint]
                child2 = [self.startPoint] + parent2[0:crossoverPoint] + parent1[crossoverPoint:self.lenChromosome] + [self.startPoint]
                #Verify if the childs are correct according to the rules
                self.childChromosomes.append(child1)
            
                self.childChromosomes.append(child2)
        
            
    #Mutate
    def mutate(self,mutationRate):
        newGen = []
        for chromosome in self.childChromosomes:
            for i in range(1,len(chromosome)-1):
                
                if random.random() < mutationRate:
                    chromosome[i] = np.random.randint(0,dna.lenChromosome-1)
            
            if self.verify(chromosome):
                newGen.append(chromosome)
        self.genesArr = newGen
        print("New pop:", self.genesArr)



points = np.array([[2,1],[4,3],[1,4],[5,1.5],[2.3,1.2],[3.4,4.2]])
dna = DNA(len(points)+1,0)
dna.populate(5,points)
print("Init:",dna.genesArr)
print("Decode:",[points[point] for i,point in enumerate(dna.genesArr)])

for j in dna.genesArr:
    dna.genesVal.append(calcDistance(j,points))
    

for i in range(1):
    dna.fitness()
    dna.prob()
    dna.selection(20)
    dna.crossover(0.3)
    dna.mutate(0.001)
    dna.genesVal.clear()
    for j in dna.genesArr:
        dna.genesVal.append(calcDistance(j,points))
    
#Best solution
print("Best distance: ",dna.genesVal[np.argmin(dna.genesVal)])
print("Solution: ",[points[point] for i,point in enumerate(dna.genesArr[np.argmin(dna.genesVal)])])
