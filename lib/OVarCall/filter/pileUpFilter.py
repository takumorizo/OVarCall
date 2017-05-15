#!/usr/bin/env python

from OVarCall.filter.pileUp import PileUpUnit, PileUp
import logging


class PileUpFilter(object):
    def __init__(self, filterParams):
        # setting of filter parameters

        self.__minBQ = 15

        self.__minDepth = 10
        self.__minDepthPlus = 0
        self.__minDepthMinus = 0
        self.__maxDepth = 100000000
        self.__maxDepthPlus = 100000000
        self.__maxDepthMinus = 100000000

        self.__minRate = 0.07
        self.__minRatePlus = 0.0
        self.__minRateMinus = 0.0
        self.__maxRate = 1.0
        self.__maxRatePlus = 1.0
        self.__maxRateMinus = 1.0

        self.__minObsNum = 3
        self.__minObsNumPlus = 0
        self.__minObsNumMinus = 0
        self.__maxObsNum = 100000000
        self.__maxObsNumPlus = 100000000
        self.__maxObsNumMinus = 100000000

        self.__minRefNum = 0
        self.__minRefNumPlus = 0
        self.__minRefNumMinus = 0
        self.__maxRefNum = 100000000
        self.__maxRefNumPlus = 100000000
        self.__maxRefNumMinus = 100000000

        self.__minOverlapDepth = 0
        self.__minOverlapRefNum = 0
        self.__minOverlapObsNum = 0

        self.__setFilterParams(filterParams)
        self.__pileUp = None
    # return satisfied Mutation List. formats of list is as follows
    # { 'M':['A','G'], 'I':['AAC','CGG'],'D':['AGG','CCT'] } <- all bases are described in upper case
    # if no such mutations are found, return None object

    def __setFilterParams(self, filterParams):
        if type(filterParams) is dict:
            if 'minBQ' in filterParams:
                self.__minBQ = int(filterParams['minBQ'])

            if 'minDepth' in filterParams:
                self.__minDepth = int(filterParams['minDepth'])
            if 'minDepthPlus' in filterParams:
                self.__minDepthPlus = int(filterParams['minDepthPlus'])
            if 'minDepthMinus' in filterParams:
                self.__minDepthMinus = int(filterParams['minDepthMinus'])
            if 'maxDepth' in filterParams:
                self.__maxDepth = int(filterParams['maxDepth'])
            if 'maxDepthPlus' in filterParams:
                self.__maxDepthPlus = int(filterParams['maxDepthPlus'])
            if 'maxDepthMinus' in filterParams:
                self.__maxDepthMinus = int(filterParams['maxDepthMinus'])

            if 'minRate' in filterParams:
                self.__minRate = float(filterParams['minRate'])
            if 'minRatePlus' in filterParams:
                self.__minRatePlus = float(filterParams['minRatePlus'])
            if 'minRateMinus' in filterParams:
                self.__minRateMinus = float(filterParams['minRateMinus'])
            if 'maxRate' in filterParams:
                self.__maxRate = float(filterParams['maxRate'])
            if 'maxRatePlus' in filterParams:
                self.__maxRatePlus = float(filterParams['maxRatePlus'])
            if 'maxRateMinus' in filterParams:
                self.__maxRateMinus = float(filterParams['maxRateMinus'])

            if 'minObsNum' in filterParams:
                self.__minObsNum = int(filterParams['minObsNum'])
            if 'minObsNumPlus' in filterParams:
                self.__minObsNumPlus = int(filterParams['minObsNumPlus'])
            if 'minObsNumMinus' in filterParams:
                self.__minObsNumMinus = int(filterParams['minObsNumMinus'])
            if 'maxObsNum' in filterParams:
                self.__maxObsNum = int(filterParams['maxObsNum'])
            if 'maxObsNumPlus' in filterParams:
                self.__maxObsNumPlus = int(filterParams['maxObsNumPlus'])
            if 'maxObsNumMinus' in filterParams:
                self.__maxObsNumMinus = int(filterParams['maxObsNumMinus'])

            if 'minRefNum' in filterParams:
                self.__minRefNum = int(filterParams['minRefNum'])
            if 'minRefNumPlus' in filterParams:
                self.__minRefNumPlus = int(filterParams['minRefNumPlus'])
            if 'minRefNumMinus' in filterParams:
                self.__minRefNumMinus = int(filterParams['minRefNumMinus'])
            if 'maxRefNum' in filterParams:
                self.__maxRefNum = int(filterParams['maxRefNum'])
            if 'maxRefNumPlus' in filterParams:
                self.__maxRefNumPlus = int(filterParams['maxRefNumPlus'])
            if 'maxRefNumMinus' in filterParams:
                self.__maxRefNumMinus = int(filterParams['maxRefNumMinus'])

            if 'minOverlapDepth' in filterParams:
                self.__minOverlapDepth = int(filterParams['minOverlapDepth'])
            if 'minOverlapRefNum' in filterParams:
                self.__minOverlapRefNum = int(filterParams['minOverlapRefNum'])
            if 'minOverlapObsNum' in filterParams:
                self.__minOverlapObsNum = int(filterParams['minOverlapObsNum'])

    def satisfiedAll(self, chromosome, position, ref, depth, bases, qualities):
        self.__setPileUp(chromosome, position, ref, depth, bases, qualities)

        ### pre checking common condition between snv and indel ###
        depth = self.__pileUp.depth
        if not (self.__minDepth <= depth):
            return None
        refAll = self.__pileUp.summary[ref.upper()] + self.__pileUp.summary[ref.lower()]
        insAll = self.__pileUp.summary['+']
        delAll = self.__pileUp.summary['-']
        if not(self.__minObsNum <= (depth - refAll + insAll + delAll)):
            return None
        ############################################################
        satisfiedList = {'M': [], 'I': [], 'D': []}
        addedNum = 0
        ### checking candidate SNV base ###
        detailSNV = None
        if self.__minObsNum <= (depth - refAll):
            for base in 'ACGT':
                if base == self.__pileUp.ref.upper():
                    continue
                satisfiedList['M'].append(base)
                addedNum += 1
        if len(satisfiedList['M']) > 0:
            detailSNV = self.__pileUp.detailedSummary(self.__minBQ)
        ####################################

        ### checking candidate Ins/Del ###
        detailInDel = None
        insOK = (self.__minObsNum <= self.__pileUp.summary['+'])
        delOK = (self.__minObsNum <= self.__pileUp.summary['-'])
        if insOK or delOK:
            detailInDel = self.__pileUp.detailedSummary()
            if insOK:
                for elem in detailInDel['I']:
                    if self.__minObsNum <= detailInDel['I'][elem]['+'] + detailInDel['I'][elem]['-']:
                        satisfiedList['I'].append(elem)
                        addedNum += 1
            if delOK:
                for elem in detailInDel['D']:
                    if self.__minObsNum <= detailInDel['D'][elem]['+'] + detailInDel['D'][elem]['-']:
                        satisfiedList['D'].append(elem)
                        addedNum += 1
        ####################################

        if addedNum > 0:
            # logging.info("pre-check satisfied list: " +  str(satisfiedList))
            return self.__filterSatisfy(satisfiedList, detailSNV, detailInDel, self.__pileUp.ref)
        else:
            return None

    def filterSatisfy(self, satisfiedList, chromosome, position, ref, depth, bases, qualities):
        self.__setPileUp(chromosome, position, ref, depth, bases, qualities)
        detailSNV = self.__pileUp.detailedSummary(self.__minBQ)
        detailInDel = self.__pileUp.detailedSummary()
        return self.__filterSatisfy(satisfiedList, detailSNV, detailInDel, self.__pileUp.ref)

    def __filterSatisfy(self, satisfiedList, detailSNV, detailInDel, ref):
        ans = {'M': [], 'I': [], 'D': []}
        if detailSNV is not None:
            # firstly check depth and refNum condition
            depth = 0
            depthPlus = 0
            depthMinus = 0

            for base in 'ACGTN':
                depthPlus += detailSNV[base.upper()]
                depthMinus += detailSNV[base.lower()]

            depth = (depthPlus + depthMinus)

            refNum = 0
            refNumPlus = 0
            refNumMinus = 0

            refNumPlus = detailSNV[ref.upper()]
            refNumMinus = detailSNV[ref.lower()]
            refNum = (refNumPlus + refNumMinus)
            if self.__minDepth <= depth <= self.__maxDepth and \
                    self.__minDepthPlus <= depthPlus <= self.__maxDepthPlus and \
                    self.__minDepthMinus <= depthMinus <= self.__maxDepthMinus and \
                    self.__minRefNum <= refNum <= self.__maxRefNum and \
                    self.__minRefNumPlus <= refNumPlus <= self.__maxRefNumPlus and \
                    self.__minRefNumMinus <= refNumMinus <= self.__maxRefNumMinus:
                # logging.info("depth ok")
                for base in satisfiedList['M']:
                    obsNum = detailSNV[base.upper()] + detailSNV[base.lower()]
                    obsNumPlus = detailSNV[base.upper()]
                    obsNumMinus = detailSNV[base.lower()]

                    # logging.info( "base : " + base)
                    if not(self.__minObsNum <= obsNum <= self.__maxObsNum and
                           self.__minObsNumPlus <= obsNumPlus <= self.__maxObsNumPlus and
                           self.__minObsNumMinus <= obsNumMinus <= self.__maxObsNumMinus):
                        # logging.info("obs num filtered")
                        # logging.info( str(obsNum) + " : obsNum")
                        # logging.info( str(obsNumPlus) + " : obsNumPlus")
                        # logging.info( str(obsNumMinus) + " : obsNumMinus")
                        # logging.info( str(self.__minObsNum) + " : self.__minObsNum")
                        # logging.info( str(self.__maxObsNum) + " : self.__maxObsNum")
                        # logging.info( str(self.__minObsNumPlus) + " : self.__minObsNumPlus")
                        # logging.info( str(self.__maxObsNumPlus) + " : self.__maxObsNumPlus")
                        # logging.info( str(self.__minObsNumMinus) + " : self.__minObsNumMinus")
                        # logging.info( str(self.__maxObsNumMinus) + " : self.__maxObsNumMinus")
                        continue

                    obsRate = 0.0
                    obsRatePlus = 0.0
                    obsRateMinus = 0.0
                    if depth > 0:
                        obsRate = (obsNum * 1.0) / depth
                    if depthPlus > 0:
                        obsRatePlus = (obsNumPlus * 1.0) / depthPlus
                    if depthMinus > 0:
                        obsRateMinus = (obsNumMinus * 1.0) / depthMinus
                    if not(self.__minRate <= obsRate <= self.__maxRate and
                           self.__minRatePlus <= obsRatePlus <= self.__maxRatePlus and
                           self.__minRateMinus <= obsRateMinus <= self.__maxRateMinus):
                        # logging.info("obs rate filtered")
                        # logging.info( str(obsRate) + " : obsRate")
                        # logging.info( str(obsRatePlus) + " : obsRatePlus")
                        # logging.info( str(obsRateMinus) + " : obsRateMinus")

                        # logging.info( str(depth) + " : depth")
                        # logging.info( str(depthPlus) + " : depthPlus")
                        # logging.info( str(depthMinus) + " : depthMinus")

                        # logging.info( str(self.__minRate) + " : self.__minRate")
                        # logging.info( str(self.__maxRate) + " : self.__maxRate")
                        # logging.info( str(self.__minRatePlus) + " : self.__minRatePlus")
                        # logging.info( str(self.__maxRatePlus) + " : self.__maxRatePlus")
                        # logging.info( str(self.__minRateMinus) + " : self.__minRateMinus")
                        # logging.info( str(self.__maxRateMinus) + " : self.__maxRateMinus")
                        continue

                    ans['M'].append(base.upper())

        if detailInDel is not None:
            # firstly check depth and refNum condition
            depth = 0
            depthPlus = 0
            depthMinus = 0
            for base in 'ACGTN':
                depthPlus += detailInDel[base.upper()]
                depthMinus += detailInDel[base.lower()]
            depth = (depthPlus + depthMinus)

            indelTotalPlus = 0
            indelTotalMinus = 0
            indelTotal = 0
            for indel in 'ID':
                for base in detailInDel[indel]:
                    indelTotalPlus += detailInDel[indel][base]['+']
                    indelTotalMinus += detailInDel[indel][base]['-']
            indelTotal = (indelTotalPlus + indelTotalMinus)

            refNum = 0
            refNumPlus = 0
            refNumMinus = 0
            refNumPlus = (depthPlus - indelTotalPlus)
            refNumMinus = (depthMinus - indelTotalMinus)
            refNum = (refNumPlus + refNumMinus)

            if self.__minDepth <= depth <= self.__maxDepth and \
                    self.__minDepthPlus <= depthPlus <= self.__maxDepthPlus and \
                    self.__minDepthMinus <= depthMinus <= self.__maxDepthMinus and \
                    self.__minRefNum <= refNum <= self.__maxRefNum and \
                    self.__minRefNumPlus <= refNumPlus <= self.__maxRefNumPlus and \
                    self.__minRefNumMinus <= refNumMinus <= self.__maxRefNumMinus:
                # logging.info("depth ok")
                for indel in 'ID':
                    for bases in satisfiedList[indel]:
                        obsNumPlus = 0
                        obsNumMinus = 0
                        obsNum = 0

                        if bases.upper() in detailInDel[indel]:
                            obsNumPlus = detailInDel[indel][bases.upper()]['+']
                            obsNumMinus = detailInDel[indel][bases.upper()]['-']
                            obsNum = obsNumPlus + obsNumMinus

                        if not(self.__minObsNum <= obsNum <= self.__maxObsNum and
                               self.__minObsNumPlus <= obsNumPlus <= self.__maxObsNumPlus and
                               self.__minObsNumMinus <= obsNumMinus <= self.__maxObsNumMinus):
                            # logging.info("obs num filtered")
                            continue

                        obsRate = 0.0
                        obsRatePlus = 0.0
                        obsRateMinus = 0.0
                        if depth > 0:
                            obsRate = (obsNum * 1.0) / depth
                        if depthPlus > 0:
                            obsRatePlus = (obsNumPlus * 1.0) / depthPlus
                        if depthMinus > 0:
                            obsRateMinus = (obsNumMinus * 1.0) / depthMinus

                        if not(self.__minRate <= obsRate <= self.__maxRate and
                               self.__minRatePlus <= obsRatePlus <= self.__maxRatePlus and
                               self.__minRateMinus <= obsRateMinus <= self.__maxRateMinus):
                            # logging.info("obs rate filtered")
                            continue

                        ans[indel].append(bases.upper())

        return ans

    def __setPileUp(self, chromosome, position, ref, depth, bases, qualities):
        if self.__pileUp is None:
            self.__pileUp = PileUp(chromosome, position, ref, depth, bases, qualities)
        else:
            self.__pileUp.setLine(chromosome, position, ref, depth, bases, qualities)

    def filterOverlapPiled(self, head, pileNums):

        pileNums = map(int, pileNums)

        depth = sum(pileNums[0:6]) + pileNums[6] + pileNums[10]
        depthPlus = sum(pileNums[0:3]) + pileNums[6] + pileNums[10]
        depthMinus = sum(pileNums[3:6]) + pileNums[6] + pileNums[10]

        rate = 0.0
        ratePlus = 0.0
        rateMinus = 0.0

        obsNum = pileNums[1] + pileNums[4] + pileNums[10]
        obsNumPlus = pileNums[1] + pileNums[10]
        obsNumMinus = pileNums[4] + pileNums[10]

        refNum = pileNums[0] + pileNums[3] + pileNums[6]
        refNumPlus = pileNums[0] + pileNums[6]
        refNumMinus = pileNums[3] + pileNums[6]


        overlapDepth = sum(pileNums[6:15])
        overlapRef = pileNums[6]
        overlapObs = pileNums[10]

        if depth > 0:
            rate = (obsNum * 1.0) / depth
        if depthPlus > 0:
            ratePlus = (obsNumPlus * 1.0) / depthPlus
        if depthMinus > 0:
            rateMinus = (obsNumMinus * 1.0) / depthMinus

        if not(self.__minDepth <= depth <= self.__maxDepth):
            return False
        if not(self.__minDepthPlus <= depthPlus <= self.__maxDepthPlus):
            return False
        if not(self.__minDepthMinus <= depthMinus <= self.__maxDepthMinus):
            return False

        if not(self.__minRate <= rate <= self.__maxRate):
            return False
        if not(self.__minRatePlus <= ratePlus <= self.__maxRatePlus):
            return False
        if not(self.__minRateMinus <= rateMinus <= self.__maxRateMinus):
            return False

        if not(self.__minObsNum <= obsNum <= self.__maxObsNum):
            return False
        if not(self.__minObsNumPlus <= obsNumPlus <= self.__maxObsNumPlus):
            return False
        if not(self.__minObsNumMinus <= obsNumMinus <= self.__maxObsNumMinus):
            return False

        if not(self.__minRefNum <= refNum <= self.__maxRefNum):
            return False
        if not(self.__minRefNumPlus <= refNumPlus <= self.__maxRefNumPlus):
            return False
        if not(self.__minRefNumMinus <= refNumMinus <= self.__maxRefNumMinus):
            return False
        if not(self.__minOverlapDepth <= overlapDepth):
            return False
        if not(self.__minOverlapObsNum <= overlapObs):
            return False
        if not(self.__minOverlapRefNum <= overlapRef):
            return False

        return True

