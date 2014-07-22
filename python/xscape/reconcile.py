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

# Base class
from xscape import ReconcileAlgorithm


class ReconcileAlgorithmWithoutRecordedEvents(ReconcileAlgorithm.ReconcileAlgorithm):

    def __init__(self, switchLo, switchHi, lossLo, lossHi):
        ReconcileAlgorithm.ReconcileAlgorithm.__init__(self, switchLo, switchHi, lossLo, lossHi)

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