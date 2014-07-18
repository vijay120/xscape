# reconcile.py
# Ran Libeskind-Hadas, Jessica Yi-Chieh Wu, Mukul Bansal, November 2013

# Memoized Pareto tree reconciliation dynamic programming solver for 
# untimed trees.

# Implements the dynamic programming algorithm described in HMC Tech Report
# "Faster Dynamic Programming Algorithms for the Cophylogeny Reconstruction
# Problem" available at www.cs.hmc.edu/~hadas/jane/TechReportCS-2011-1.pdf
# This implementation uses memoization rather than DP.

# A tree is represented as a dictionary of key-value pairs where a key is an
# edge name and the value is a tuple of the form
# (start vertex, end vertex, left child edge name, right child edge name)
# An edge name may be None.  The "dummy" edge leading to the root of the
# parasite tree, denoted e^P in the technical report, must be named "pTop".

# python libraries
from collections import *

# xscape libraries
from common import *
from CostVector import *


class ReconcileAlgorithm:

    def __init__(self, switchLo, switchHi, lossLo, lossHi):

        # The three dictionaries below correspond to the A, C, and Best DP tables
        # described in the technical report.  These are set to None here but initialized
        # in the reconcile function.  They need to be reset before each run since
        # they change for each sigscape randomization trial.
        self.Amemo = {}
        self.Cmemo = {}
        self.Bestmemo = {}

        # The Ancestors and Descendants dictionaries are precomputed in the reconcile
        # function and allow the switch function to determine the valid landing sites
        # for a switch. 
        self.Ancestors = {}
        self.Descendants = {}

        # The switchLo, switchHi, lossLo, and lossHi values are the user-specified
        # low and high ranges for the switch and loss costs, relative to the unit
        # cost of duplication.  They are shown here only for clarity.  
        self.switchLo = switchLo
        self.switchHi = switchHi
        self.lossLo = lossLo
        self.lossHi = lossHi

    # This is the main function for this file.  It seeks to find the best
    # reconciliation for the parasite tree, rooted at every possible edge of the
    # host tree.
    def reconcile(self, parasiteTree, hostTree, parasiteToHostMapping):
        ''' Takes dictionary representations of the parasite tree, host tree
            and parasiteToHostMapping as input and returns a list of the Pareto optimal solutions. '''
            
        self.ancestorsAndDescendants(hostTree) # Set the Ancestors and Descendants

        solutions = reduce(lambda countAndCountVectorA, countAndCountVectorB: countAndCountVectorA + countAndCountVectorB, 
                        map(lambda hostEdge: self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, "pTop", hostEdge), 
                            hostTree.keys()))

        return self.paretoFilter(solutions)

    def optimalEdgeCost(self, parasiteTree, hostTree, parasiteToHostMapping, parasiteEdge, hostEdge):
        ''' The optimalEdgeCost table for the dynamic program. '''
                
        if (parasiteEdge, hostEdge) in self.Cmemo: return self.Cmemo[(parasiteEdge, hostEdge)]
        
        # Option 1:  Pass through
        passThrough = self.aliveEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdge, hostEdge)
        
        if self.tipEdge(parasiteEdge, parasiteTree):   # The options below don't apply to tips
            return passThrough

        else:
            parasiteEdgeLeftChild = self.leftChildEdge(parasiteEdge, parasiteTree)
            parasiteEdgeRightChild = self.rightChildEdge(parasiteEdge, parasiteTree)
        
            # Option 2:  Duplicate here
            duplicate = CostVector(0, 1, 0, 0, 1) * \
                        self.merge(self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdgeLeftChild, hostEdge), \
                              self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdgeRightChild, hostEdge)) 
        
            # Option 3:  Switch here
            switch1 = CostVector(0, 0, 1, 0, 1) * \
                      self.merge(self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdgeLeftChild, hostEdge), \
                            self.switches(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdgeRightChild, hostEdge))

            switch2 = CostVector(0, 0, 1, 0, 1) * \
                      self.merge(self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdgeRightChild, hostEdge), \
                            self.switches(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdgeLeftChild, hostEdge))

            switch = switch1 + switch2
            
            output = self.paretoFilter(passThrough + duplicate + switch) 
            self.Cmemo[(parasiteEdge, hostEdge)] = output
            return output


    def aliveEdgeCost(self, parasiteTree, hostTree, parasiteToHostMapping, parasiteEdge, hostEdge):
        ''' The A table for the dynamic program. '''
            
        if (parasiteEdge, hostEdge) in self.Amemo: return self.Amemo[(parasiteEdge, hostEdge)]

        if self.tipEdge(hostEdge, hostTree):
            if self.tipEdge(parasiteEdge, parasiteTree) and \
                parasiteToHostMapping[self.endVertex(parasiteEdge, parasiteTree)] == self.endVertex(hostEdge, hostTree):
                return [CostVector(0, 0, 0, 0, 1)]
            else:
                return [CostVector(INF, INF, INF, INF, 0)]
        else:
            hostEdgeLeftChild = self.leftChildEdge(hostEdge, hostTree)
            hostEdgeRightChild = self.rightChildEdge(hostEdge, hostTree)

            # Cospeciation
            if self.tipEdge(parasiteEdge, parasiteTree):
                cospeciation = [CostVector(INF, INF, INF, INF, 0)]
            else:
                parasiteEdgeLeftChild = self.leftChildEdge(parasiteEdge, parasiteTree)
                parasiteEdgeRightChild = self.rightChildEdge(parasiteEdge, parasiteTree)

                cospeciation1 = CostVector(1, 0, 0, 0, 1) * \
                  self.merge(self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdgeLeftChild, hostEdgeLeftChild), \
                        self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdgeRightChild, hostEdgeRightChild))

                cospeciation2 = CostVector(1, 0, 0, 0, 1) * \
                  self.merge(self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdgeLeftChild, hostEdgeRightChild), \
                        self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdgeRightChild, hostEdgeLeftChild)) \

                cospeciation = cospeciation1 + cospeciation2
                
            # Loss
            loss1 = CostVector(0, 0, 0, 1, 1) * \
                    self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdge, hostEdgeLeftChild)

            loss2 = CostVector(0, 0, 0, 1, 1) * \
                    self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, parasiteEdge, hostEdgeRightChild)

            loss = loss1 + loss2
            
            output = self.paretoFilter(cospeciation + loss) 
            self.Amemo[(parasiteEdge, hostEdge)] = output
            return output
            
    def merge(self, CVlist1, CVlist2):
        ''' Given two lists of CostVectors, returns a new list of CostVectors, each
            of which is the sum of a pair of vectors from the two given lists.'''
        output = []
        for v in CVlist1:
            for w in CVlist2:
                output.append(v+w)
        return output

    def switches(self, parasiteTree, hostTree, parasiteToHostMapping, ep, eh):
        ''' Returns the list of all CostVectors in which the given parasite edge ep
            switches to all possible host edges. '''
                    
        if (ep, eh) in self.Bestmemo: return self.Bestmemo[(ep, eh)]
        output = []
        for switchEdge in hostTree:     # for every possible host edge
            if switchEdge not in self.Ancestors[eh] and \
               switchEdge not in self.Descendants[eh]:
                if switchEdge not in self.Descendants[eh]:
                    output.extend(self.optimalEdgeCost(parasiteTree, hostTree, parasiteToHostMapping, ep, switchEdge))
        self.Bestmemo[(ep, eh)] = output
        return output

    def paretoFilter(self, CVlist):
        ''' Returns the Pareto front for the given list of CostVectors. '''
        CVlist = self.CVfilter(CVlist)
        uniqueCVlist = self.coalesceDuplicates(CVlist)
        uniqueCVlist.sort(CostVector.lex)
        if len(uniqueCVlist) == 1: return uniqueCVlist
        lexlist = [uniqueCVlist[0]]
        for i in range(1, len(uniqueCVlist)):
            predecessor = uniqueCVlist[i-1]
            current = uniqueCVlist[i]
            if predecessor.d < current.d or predecessor.s < current.s:
                lexlist.append(current)
        output = []
        for v in lexlist:
            if self.minimal(v, CVlist): output.append(v)
        return output        

    def CVfilter(self, CVlist):
        ''' Filter the CVlist to a subset that removes those cost vectors that
            cannot be optimal in the given cost range. '''
                    
        if CVlist == []: return []
        else:
            LUB = min([cv.d + cv.s * self.switchHi + cv.l * self.lossHi for cv in CVlist])
            output = []
            for cv in CVlist:
                mincost = cv.d + cv.l * self.lossLo + cv.s * self.switchLo 
                if mincost <= LUB:
                    output.append(cv)
            return output
                
    def coalesceDuplicates(self, CVlist):
        ''' Takes a cost vector list as input and returns a list with duplicates
            removed and with event counts updated accordingly. '''
        counts = {}
        for cv in CVlist:
            cvtup = cv.toTupleCDSL()
            count = cv.count
            if cvtup not in counts:
                counts[cvtup] = count
            else:
                counts[cvtup] = counts[cvtup] + count
        output = []
        for entry in counts:
            entryWithCounts = entry + (counts[entry],)
            output.append(self.tupleToCV(entryWithCounts))
        return output

    def tupleToCV(self, entry):
        return CostVector(entry[0], entry[1], entry[2], entry[3], entry[4])

    def minimal(self, v, CVlist):
        ''' Returns True if v is a minimal element of the CostVector list. '''
        for w in CVlist:
            if w < v: return False
        return True

    def descendants(self, edge, tree):
        ''' returns the list of descendant edges of the given edge in the
            given tree.'''
        if self.tipEdge(edge, tree): return []
        else: return [self.leftChildEdge(edge, tree), self.rightChildEdge(edge, tree)] + \
                      self.descendants(self.leftChildEdge(edge, tree), tree) + \
                      self.descendants(self.rightChildEdge(edge, tree), tree)
                      
    def ancestorsAndDescendants(self, tree):
        ''' Returns a two dictionaries D and A, where D[e] is the list
            of descendant edges of e and A[e] is the list of ancestral edges
            of e. '''
        
        # First, compute all the descendants
        for e in tree:
            self.Descendants[e] = self.descendants(e, tree)
        # Next, compute ancestors 
        for e in tree: self.Ancestors[e] = []  # initialize A dictionary      
        for e in tree:
            for d in self.Descendants[e]: # d descendant of e => e ancestor of d
                self.Ancestors[d].append(e)
            
    def tipEdge(self, edge, tree):
        ''' returns True if the edge terminates at a tip  '''
        return self.leftChildEdge(edge, tree) == None  # This edge has no edge children
        
    def startVertex(self, edge, tree):
        ''' returns the start vertex of an edge '''
        return tree[edge][0]

    def endVertex(self, edge, tree):
        ''' returns the end vertex of an edge '''
        return tree[edge][1]

    def leftChildEdge(self, edge, tree):
        ''' returns the left child edge of the given edge '''
        return tree[edge][2]

    def rightChildEdge(self, edge, tree):
        ''' returns the right child edge of the given edge '''
        return tree[edge][3]
