# reconcileEvents.py
# Ran Libeskind-Hadas, Jessica Yi-Chieh Wu, Mukul Bansal, November 2013
# memoized Pareto tree reconciliation software for untimed trees

# A tree is represented as a dictionary of key-value pairs where a key is an
# edge name and the value is a tuple of the form
# (start vertex, end vertex, left child edge name, right child edge name)
# An edge name may be None.  The "dummy" edge leading to the root of the
# parasite tree must be named "pTop".

# python libraries
from collections import *
import copy

# xscape libraries
from common import *
from CostVector import *

# Base class
from xscape import ReconcileAlgorithm

# looking at union OR intersection of events?
class Config:     # to get around immutable globals
    pass
CONFIG = Config()
CVseen = defaultdict(bool) # Initially all values are False by default
CVcommonEvents = defaultdict(set)

class ReconcileAlgorithmWithRecordedEvents(ReconcileAlgorithm.ReconcileAlgorithm):

    def __init__(self, switchLo, switchHi, lossLo, lossHi):
        ReconcileAlgorithm.ReconcileAlgorithm.__init__(self, switchLo, switchHi, lossLo, lossHi)

        # The CVevents global dictionary has keys that are tuples of the form
        # (ep, eh, eventType c, d, s, l) and values that are all the events in 
        # that solution of the form (ep', eh', eventTypeString) associated with 
        # ep on eh with cost vector <c, d, s, l>
        self.CVevents = defaultdict(set)
        self.CVallEvents = defaultdict(set)
        self.CandidateCVlist = list()

    # This is the main function for this file.  It seeks to find the best
    # reconciliation for the parasite tree, rooted at every possible edge of the
    # host tree.
    def reconcile(self, parasiteTree, hostTree, phi):
        ''' Takes dictionary representations of the parasite tree, host tree
            and phi as input and returns a list of the Pareto optimal solutions. '''
        
        self.ancestorsAndDescendants(hostTree) # Set the Ancestors and Descendants
        
        solutions = []
        for eh in hostTree:
            solutions.extend(self.optimalEdgeCost(parasiteTree, hostTree, phi, "pTop", eh))
        return(self.paretoFilter(solutions))

    # The A and C functions implement the A and C DPs in the HMC Tech Report
    # "Faster Dynamic Programming Algorithms for the Cophylogeny Reconstruction
    # Problem" available at www.cs.hmc.edu/~hadas/jane/TechReportCS-2011-1.pdf
    # This implementation uses memoization rather than DP.

    def A(self, parasiteTree, hostTree, phi, ep, eh):

        if (ep, eh) in self.Amemo: return self.Amemo[(ep, eh)]

        if self.tipEdge(eh, hostTree):
            if self.tipEdge(ep, parasiteTree) and \
               phi[self.endVertex(ep, parasiteTree)] == self.endVertex(eh, hostTree):
                return [CostVector(0, 0, 0, 0, 1)]
            else:
                return [CostVector(INF, INF, INF, INF, 0)]
        else:
            ehLeftChild = self.leftChildEdge(eh, hostTree)
            ehRightChild = self.rightChildEdge(eh, hostTree)

            # Cospeciation
            if self.tipEdge(ep, parasiteTree):
                cospeciation = [CostVector(INF, INF, INF, INF, 0)]
            else:
                epLeftChild = self.leftChildEdge(ep, parasiteTree)
                epRightChild = self.rightChildEdge(ep, parasiteTree)

                cospeciation1 = self.merge(self.optimalEdgeCost(parasiteTree, hostTree, phi, \
                                        epLeftChild, ehLeftChild), \
                                      self.optimalEdgeCost(parasiteTree, hostTree, phi, \
                                        epRightChild, ehRightChild), \
                                      ep, eh, epLeftChild, ehLeftChild, \
                                      epRightChild, ehRightChild, "cospeciation")

                cospeciation2 = self.merge(self.optimalEdgeCost(parasiteTree, hostTree, phi, \
                                        epLeftChild, ehRightChild), \
                                      self.optimalEdgeCost(parasiteTree, hostTree, phi, \
                                        epRightChild, ehLeftChild), \
                                      ep, eh, epLeftChild, ehRightChild, \
                                      epRightChild, ehLeftChild, "cospeciation")

                cospeciation = cospeciation1 + cospeciation2
                
            # Loss
            loss1 = self.lossmerge(ep, eh, ehLeftChild, \
                         self.optimalEdgeCost(parasiteTree, hostTree, phi, ep, ehLeftChild))
                           
            loss2 = self.lossmerge(ep, eh, ehRightChild, \
                         self.optimalEdgeCost(parasiteTree, hostTree, phi, ep, ehRightChild))

            loss = loss1 + loss2
            
            output = self.paretoFilter(cospeciation + loss)
            self.Amemo[(ep, eh)] = output
            return output
            
    def optimalEdgeCost(self, parasiteTree, hostTree, phi, ep, eh):
        ''' The optimalEdgeCost table for the dynamic program. '''
                
        if (ep, eh) in self.Cmemo: return self.Cmemo[(ep, eh)]

        # Option 1:  Pass through
        passThrough = self.A(parasiteTree, hostTree, phi, ep, eh)
        
        if self.tipEdge(ep, parasiteTree):   # The options below don't apply to tips
            return passThrough

        else:
            epLeftChild = self.leftChildEdge(ep, parasiteTree)
            epRightChild = self.rightChildEdge(ep, parasiteTree)
        
            # Option 2:  Duplicate here
            duplicate = self.merge(self.optimalEdgeCost(parasiteTree, hostTree, phi, epLeftChild, eh), \
                              self.optimalEdgeCost(parasiteTree, hostTree, phi, epRightChild, eh), \
                                ep, eh, epLeftChild, eh, \
                                epRightChild, eh, "duplication")
        
            switch = []
            leftCVlist = self.optimalEdgeCost(parasiteTree, hostTree, phi, epLeftChild, eh)
            rightPairs = self.switches(parasiteTree, hostTree, phi, epRightChild, eh)
            for (switchEdge, rightCVlist) in rightPairs:
                switch.extend(self.merge(leftCVlist, rightCVlist, \
                                    ep, eh, epLeftChild, eh, epRightChild, \
                                    switchEdge, "switch"))
                
            leftCVlist = self.optimalEdgeCost(parasiteTree, hostTree, phi, epRightChild, eh)
            rightPairs = self.switches(parasiteTree, hostTree, phi, epLeftChild, eh)
            for (switchEdge, rightCVlist) in rightPairs:
                switch.extend(self.merge(leftCVlist, rightCVlist, \
                                    ep, eh, epRightChild, eh, epLeftChild, \
                                    switchEdge, "switch"))
                   
        output = self.paretoFilter(passThrough + duplicate + switch)
        self.Cmemo[(ep, eh)] = output

        return output

    def merge(self, CVlist1, CVlist2, ep, eh, epChild1, ehChild1, epChild2, \
              ehChild2, eventType):
        ''' Given two lists of CostVectors, returns a new list of CostVectors, each
            of which is the sum of a pair of vectors from the two given lists.'''

        intersection = CONFIG.intersection
        if intersection:
            global CVseen
            global CVcommonEvents
            
        output = []
        for v in CVlist1:
            for w in CVlist2:
                if eventType == "cospeciation":
                    newCV = CostVector(1, 0, 0, 0, 1) + v + w
                elif eventType == "duplication":
                    newCV = CostVector(0, 1, 0, 0, 1) + v + w
                else:   # eventType == "switch":
                    newCV = CostVector(0, 0, 1, 0, 1) + v + w
                keepnewCV = True
                for cv in self.CandidateCVlist:
                    if cv < newCV: 
                        keepnewCV = False
                        break
        
                if keepnewCV:
                    output.append(newCV)
                    vsoln = (epChild1, ehChild1) + v.toTupleCDSL()
                    wsoln = (epChild2, ehChild2) + w.toTupleCDSL()
                    if eventType == "switch": eventType = "switch to "+str(ehChild2)
                    nswe = (ep, eh, eventType) + newCV.toTupleCDSL()
                    ns = (ep, eh) + newCV.toTupleCDSL()
                    self.CVevents[nswe].add(nswe)
                    self.CVevents[nswe] = self.CVevents[nswe].\
                                     union(self.CVallEvents[vsoln]).\
                                     union(self.CVallEvents[wsoln])
                    self.CVallEvents[ns] = self.CVallEvents[ns].union(self.CVevents[nswe])

                    if intersection:
                        key = newCV.toTupleCDSL()
                        if CVseen[key]:
                            CVcommonEvents[key] = CVcommonEvents[key].\
                                                  intersection(self.CVevents[nswe])
                        else:
                            CVcommonEvents[key] = self.CVevents[nswe]
                            CVseen[key] = True
        return output

    def lossmerge(self, ep, eh, ehChild, CVlist):

        intersection = CONFIG.intersection
        if intersection:
            global CVseen
            global CVcommonEvents
        
        output = []
        for v in CVlist:
            newCV = CostVector(0, 0, 0, 1, 1) + v
            keepnewCV = True
            for cv in self.CandidateCVlist:
                if cv < newCV: 
                    keepnewCV = False
                    break
        
            if keepnewCV:
                output.append(newCV)
                vsoln = (ep, ehChild) + v.toTupleCDSL()
                nswe = (ep, eh, "loss "+str(ehChild)) + newCV.toTupleCDSL()
                ns = (ep, eh) + newCV.toTupleCDSL()
                self.CVevents[nswe].add(nswe)
                self.CVevents[nswe] = self.CVevents[nswe].union(self.CVallEvents[vsoln])
                self.CVallEvents[ns] = self.CVallEvents[ns].union(self.CVevents[nswe])

                if intersection:
                    key = newCV.toTupleCDSL()
                    if CVseen[key]:
                        self.CVcommonEvents[key] = CVcommonEvents[key].\
                                                  intersection(self.CVevents[nswe])
                    else:
                        self.CVcommonEvents[key] = self.CVevents[nswe]
                        self.CVseen[key] = True
        return output

    def allSwitches(self, parasiteTree, hostTree, phi, ep, eh):
        ''' Returns the list of all CostVectors in which the given parasite edge ep
            switches to all possible host edges. '''
        if (ep, eh) in self.Bestmemo: return self.Bestmemo[(ep, eh)]
        output = []
        for switchEdge in hostTree:     # for every possible host edge
            if switchEdge not in self.Ancestors[eh] and \
               switchEdge not in self.Descendants[eh]:
                output.append((switchEdge, \
                                   self.optimalEdgeCost(parasiteTree, hostTree, phi, ep, switchEdge)))
        self.Bestmemo[(ep, eh)] = output
        return output