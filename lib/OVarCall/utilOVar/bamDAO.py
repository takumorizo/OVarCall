#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from OVarCall.utilOVar.exceptionUtil import InvalidBamFileManagement,InvalidMutationType,InsufficientBamReadFilter
import pysam
import math
import distutils.util
import logging
import contextlib
import re
import copy

class BamDAO:
    def __init__(self, bamPath,settings=None):
        # parameters for bam file object
        self.__bamPath = bamPath
        self.__bam = None
        self.__chrDic = None

        # setting parameters for overlapping pileup
        self.__windowSize = 200
        self.__minBQ = 15
        self.__minMapQ = 30
        self.__f = 2
        self.__F = 3840

        # parameter for read Filter
        self.__maxInsDel = 2
        self.__maxSNV = 2
        self.__maxMutAll = 3
        self.__maxSCProportion = 0.25
        self.__maxLowMapReadProportion = 0.50
        self.__nearestIndelDistance = 25 #  | mutPos - indesPos | > 25
        self.__nearestIndelReadThres = 20   # number of near indel read thres.

        # parameter for BQ variant read filter
        # self.__maxLowBQVariantNum = -1
        # self.__maxLowBQVariantProportion = 1.0
        self.__lowBQ = 15
        self.__minAvgBaseQuality = -1

        self.setParameters(settings)

    def setParameters(self,settings):
        if type(settings) is dict :
            if 'windowSize' in settings :
                self.__windowSize = int(settings['windowSize'])
            if 'minBQ' in settings:
                self.__minBQ = int(settings['minBQ'])
            if 'minMapQ' in settings:
                self.__minMapQ = int(settings['minMapQ'])
            if 'f' in settings:
                self.__f = int(settings['f'])
            if 'F' in settings:
                self.__F = int(settings['F'])
            if 'maxInsDel' in settings:
                self.__maxInsDel = int(settings['maxInsDel'])
            if 'maxSNV' in settings:
                self.__maxSNV = int(settings['maxSNV'])
            if 'maxMutAll' in settings:
                self.__maxMutAll = int(settings['maxMutAll'])
            if 'maxSCProportion' in settings:
                self.__maxSCProportion = float(settings['maxSCProportion'])
            if 'maxLowMapReadProportion' in settings:
                self.__maxLowMapReadProportion = float(settings['maxLowMapReadProportion'])
            if 'nearestIndelDistance' in settings:
                self.__nearestIndelDistance = int(settings['nearestIndelDistance'])
            if 'nearestIndelReadThres' in settings:
                self.__nearestIndelReadThres = int(settings['nearestIndelReadThres'])
            if 'lowBQ' in settings:
                self.__lowBQ = int(settings['lowBQ'])
            if 'minAvgBaseQuality' in settings:
                self.__minAvgBaseQuality = int(settings['minAvgBaseQuality'])
            # if 'maxLowBQVariantNum' in settings:
            #     self.__maxLowBQVariantNum = int(settings['maxLowBQVariantNum'])
            # if 'maxLowBQVariantProportion' in settings:
            #     self.__maxLowBQVariantProportion = float(settings['maxLowBQVariantProportion'])

        # logging.info('maxLowBQVariantNum : ' +str(self.__maxLowBQVariantNum))
        # logging.info('maxLowBQVariantProportion : ' +str(self.__maxLowBQVariantProportion))
      
    def getHeaderDic(self):
        if self.__bam is not None:
            return copy.deepcopy(self.__bam.header)
        else :
            return None

    @contextlib.contextmanager
    def openBam(self,bamPath=None):
        try:
            if bamPath is not None:
                self.__bamPath = bamPath
            self.__bam = pysam.AlignmentFile(self.__bamPath, "rb")
            self.__chrDic = {}
            for chrIdx in range(len(self.__bam.header['SQ'])):
                self.__chrDic[self.__bam.header['SQ'][chrIdx]['SN']] = chrIdx
            yield
        finally:
            self.__closeBam()

    def __closeBam(self):
        if self.__bam is not None:
            self.__bam.close()
        self.__bam = None

    def __filterReads(self,reads,TYPE,Chr,pos,ref,obs):
        numAll = 0
        numFiltered = 0
        readList = []

        indelReadsNum = 0
        for read in reads:
            if self.__filterFlaggedRead(read) :
                if self.__filterLowMapRead(read,TYPE,Chr,pos,ref,obs):
                    readList.append(read)
                else:
                    numFiltered += 1
            if self.__nearestIndelDistance > 0 and self.__distToNearestIndel(read,TYPE,Chr,pos,ref,obs) > 0:
                indelReadsNum += 1
            numAll += 1

        proportionLowMap = 0.0
        if numAll > 0 :
            proportionLowMap = (1.0 * numFiltered) / (1.0 * numAll)
        if numAll == 0:
            readList = None
        if not(proportionLowMap <= self.__maxLowMapReadProportion):
            readList = None
        if self.__nearestIndelDistance > 0 and not(indelReadsNum <= self.__nearestIndelReadThres ):
            readList = None
        if self.__minAvgBaseQuality > 0 and TYPE == 'M':
            profile = self.__getReadBQProfile(ref, obs, Chr, pos, TYPE, self.__lowBQ, reads)
            if not( self.__minAvgBaseQuality <= profile['avg']):
                logging.info(str(Chr)+" " + str(pos) + " " + TYPE + " " + str(ref) + " " + str(obs) + " : low avg bq variant pos fitered" )
                logging.info(str(profile))
                readList = None
        
        return readList

    def __filterLowMapRead(self,alignedSeg,TYPE,Chr,pos,ref,obs):
        if alignedSeg.mapping_quality < self.__minMapQ:
            return False

        delNum = self.__delNum(alignedSeg)
        insNum = self.__insNum(alignedSeg)
        if not(delNum + insNum <= self.__maxInsDel):
            return False

        if not(self.__softClipProportion(alignedSeg) <= self.__maxSCProportion):
            return False
        snvNum = self.__snvNum(alignedSeg)
        if not(snvNum <= self.__maxSNV):
            return False
        if not(delNum + insNum + snvNum <= self.__maxMutAll):
            return False
        return True

    # pos is already transformed to the pysam Position
    def __distToNearestIndel(self,alignedSeg,TYPE,Chr,pos,ref,obs):
        if TYPE == "M":
            if self.__delNum(alignedSeg) + self.__insNum(alignedSeg) >= 1:
                mappedDic = self.__makeMappedDic(alignedSeg)
                if (Chr,pos) in mappedDic:
                    for x in range(self.__nearestIndelDistance):
                        if x == 0 :
                            continue
                        if (Chr,pos+x) in mappedDic:
                            if mappedDic[(Chr,pos+x)]['base'] == '-' and mappedDic[(Chr,pos+x)]['cigar'] == 2:
                                return x
                            elif len(mappedDic[(Chr,pos+x)]['list']) > 0 and mappedDic[(Chr,pos+x)]['list'][0][2] == 1:
                                return x
                        if (Chr,pos-x) in mappedDic:
                            if mappedDic[(Chr,pos-x)]['base'] == '-' and mappedDic[(Chr,pos-x)]['cigar'] == 2:
                                return x
                            elif len(mappedDic[(Chr,pos-x)]['list']) > 0 and mappedDic[(Chr,pos-x)]['list'][0][2] == 1:
                                return x
                return -1
            else :
                return -1
        else:
            return -1


    def __filterFlaggedRead(self,alignedSeg):
        if not(int(alignedSeg.flag) & self.__f == self.__f and int(alignedSeg.flag) & self.__F == 0):
            return False
        return True
    

    def __softClipProportion(self,alignedSeg):
        total = 0
        SC = 0
        cigars = alignedSeg.cigar
        for cigar in cigars:
            if cigar[0] == 4 : # deletion type cigar type
                SC += cigar[1]
            total += cigar[1]
        ans = 0.0
        if total > 0 :
            ans = (SC * 1.0) / (total * 1.0)
        return ans


    def __snvNum(self,alignedSeg):
        ans = 0
        mdString = ""
        tags = alignedSeg.tags
        for tag in tags:
            if tag[0] == 'MD':
                mdString = tag[1]

        while len(mdString) > 0:
            firsts = self.__firstMD(mdString)
            first = firsts["firstMD"]
            types = firsts["type"]
            mdString = mdString[len(first):len(mdString)]
            if types == "mutation":
                ans += 1
        
        return ans

    def __delNum(self,alignedSeg):
        cigars = alignedSeg.cigar
        ans = 0
        for cigar in cigars:
            if cigar[0] == 2 : # deletion type cigar type
                ans += 1
        return ans

    def __insNum(self,alignedSeg):
        cigars = alignedSeg.cigar
        ans = 0
        for cigar in cigars:
            if cigar[0] == 1 : # deletion type cigar type
                ans += 1
        return ans

    def __firstMD(self,mdString):
        MDMatPattern = re.compile('(.*?)([0-9]+)(.*)')
        MDMutPattern = re.compile('(.*?)([A-Za-z])(.*)')
        MDDelPattern = re.compile('(.*?)(\^[A-Za-z]+)(.*)')
        matResult = MDMatPattern.search(mdString)
        if len(matResult.groups()[0]) == 0:
            return {"firstMD": matResult.groups()[1], "type": "match"}
        mutResult = MDMutPattern.search(mdString)
        if len(mutResult.groups()[0]) == 0:
            return {"firstMD": mutResult.groups()[1], "type": "mutation"}
        delResult = MDDelPattern.search(mdString)
        if len(delResult.groups()[0]) == 0:
            return {"firstMD": delResult.groups()[1], "type": "deletion"}


    def __makeMappedDic(self,alignedSeg):
        Bases = alignedSeg.query_sequence
        BQAry = alignedSeg.query_qualities
        baseIdx = 0
        mapIdx = alignedSeg.reference_start
        mapChr = alignedSeg.reference_id  # num
        remCigar = alignedSeg.cigar
        firstCigar = None
        if len(remCigar) > 0:
            firstCigar = remCigar[0]
            remCigar.pop(0)
        else:
            firstCigar = None
        ans = {}
        while firstCigar is not None:
            mlen = self.__MLenOne(firstCigar)
            clen = self.__CLenOne(firstCigar)
            cigartype = firstCigar[0]
            if mlen == clen and mlen > 0:  # M,EQ,X
                for at in range(mlen):
                    ans[(mapChr, mapIdx + at)] = {'base': Bases[baseIdx + at],
                                                  'bq': BQAry[baseIdx + at],
                                                  'cigar': cigartype, 'list': []}
            elif mlen == 0 and clen > 0:  # I,S
                if((mapChr, mapIdx - 1)) in ans:
                    ans[(mapChr, mapIdx - 1)]['list'].append(
                        (Bases[baseIdx:baseIdx + clen],
                         BQAry[baseIdx:baseIdx + clen],
                         cigartype))
                else:
                    ans[(mapChr, -1)] = {'base': None,
                                         'bq': None,
                                         'cigar': None,
                                         'list':
                                         [(Bases[baseIdx:baseIdx + clen], BQAry[baseIdx:baseIdx + clen], cigartype)]}
            elif mlen > 0 and clen == 0 and cigartype != 3:  # D,P do not consider cigar N
                for at in range(mlen):
                    ans[(mapChr, mapIdx + at)] = {'base': '-',
                                                  'bq': -1,
                                                  'cigar': cigartype, 'list': []}
            elif mlen == 0 and clen == 0:  # H
                # this case is hard clipping and ignore this case
                if((mapChr, mapIdx - 1)) in ans:
                    ans[(mapChr, mapIdx - 1)]['list'].append(
                        ('-' * firstCigar[1],
                         '-' * firstCigar[1],
                         cigartype))
                else:
                    ans[(mapChr, -1)] = {'base': None,
                                         'bq': None,
                                         'cigar': None,
                                         'list':
                                         [('-' * firstCigar[1], '-' * firstCigar[1], cigartype)]}
                pass

            baseIdx += clen
            mapIdx += mlen
            if len(remCigar) > 0:
                firstCigar = remCigar[0]
                remCigar.pop(0)
            else:
                firstCigar = None
        return ans

    def __MLenOne(self,cigar):
        try:
            ctype = cigar[0]

            if ctype == 0 or ctype == 2 or ctype == 3 or ctype == 7 or ctype == 8 or ctype == 6:
                return cigar[1]
            else:
                return 0
        except Exception, e:
            sys.stderr.writelines("isMlenOne invalid form of cigar")
            raise e


    def __CLenOne(self,cigar):
        try:
            ctype = cigar[0]

            if ctype == 0 or ctype == 1 or ctype == 4 or ctype == 7 or ctype == 8:
                return cigar[1]
            else:
                return 0
        except Exception, e:
            sys.stderr.writelines("isCLenOne invalid form of cigar")
            raise e



    """
    getState(ref,obs,Chr,pos,read,obsType):
    の pos -> 0-indexed

    'M' : pos の場所(0-indexed pysam 上) を
    'I' : pos の場所(0-indexed pysam 上) の mappedList の ins を見る
    'D' : 0-indexed pysam 上 の [pos,pos+len) を見る
    """
    # True : there is non reference base read
    def __checkSNV(self,ref,obs,Chr,pos,read,obsType,mapppedDic):
        baseAtPos = None
        if (Chr,pos) in mapppedDic:
            baseAtPos = mapppedDic[(Chr,pos)]['base']
            baseAtPos = baseAtPos.upper()
            isref = (ref == baseAtPos)
            return not(isref)
        else : # no mapping
            return None
    
    def __checkIns(self,ref,obs,Chr,pos,read,obsType,mapppedDic):
        if (Chr,pos) in mapppedDic:    
            if len(mapppedDic[(Chr,pos)]['list']) > 0:
                for cElem in mapppedDic[(Chr,pos)]['list']:
                    if cElem[2] == 1: # cigar I
                        return True
                return False
            else : # no Ins
                return False
        else : # no mapping
            return None

    def __checkDel(self,ref,obs,Chr,pos,read,obsType,mapppedDic):
        if (Chr,pos) in mapppedDic:
            delList = []
            delLen = 0
            while (Chr,pos+1+delLen) in mapppedDic and mapppedDic[(Chr,pos+1+delLen)]['cigar'] == 2 :
                delList.append(mapppedDic[(Chr,pos+1+delLen)])
                delLen += 1
            delString=''
            for Del in delList:
                if Del['base'] == '-':
                    delString += Del['base']
                else :
                    break
            if len(delString) > 0 :
                return True
            else :
                return False
        else: # no mapping
            return None



    def __getState2(self,ref,obs,Chr,pos,read1,read2,obsType,minBQ=None):
        if not read1.is_reverse:
            state1 = self.__getState(ref,obs,Chr,pos,read1,obsType,minBQ)
            state2 = self.__getState(ref,obs,Chr,pos,read2,obsType,minBQ)
        else :
            state1 = self.__getState(ref,obs,Chr,pos,read2,obsType,minBQ)
            state2 = self.__getState(ref,obs,Chr,pos,read1,obsType,minBQ)
        refNum1 = 0
        obsNum1 = 0
        otherNum1 = 0
        
        refNum2 = 0
        obsNum2 = 0
        otherNum2 = 0
        
        if state1['ref']:
            refNum1 += 1
        if state1['obs']:
            obsNum1 += 1
        if state1['other'] is not None:
            otherNum1 += 1
        if state2['ref']:
            refNum2 += 1
        if state2['obs']:
            obsNum2 += 1
        if state2['other'] is not None:
            otherNum2 += 1
            
        Rp = refNum1
        Obp = obsNum1
        Otp = otherNum1
        Rm = refNum2
        Obm = obsNum2
        Otm = otherNum2

        RR = 0
        ROb = 0
        ROt = 0

        ObR = 0
        ObOb = 0
        ObOt = 0
        
        OtR = 0
        OtOb = 0
        OtOt = 0
        
        if (refNum1+obsNum1+otherNum1) == 1 and (refNum2+obsNum2+otherNum2) == 1:
            # Not counting as single read
            Rp = 0
            Obp = 0
            Otp = 0
            Rm = 0
            Obm = 0
            Otm = 0
            if refNum1 > 0 and refNum2 > 0 :
                RR += 1
            elif refNum1 > 0 and obsNum2 > 0:
                ROb += 1
            elif refNum1 > 0 and otherNum2 > 0:
                ROt += 1
            elif obsNum1 > 0 and refNum2 > 0 :
                ObR += 1
            elif obsNum1 > 0 and obsNum2 > 0:
                ObOb += 1
            elif obsNum1 > 0 and otherNum2 > 0 :
                ObOt += 1
            elif otherNum1 > 0 and refNum2 > 0 :
                OtR += 1
            elif otherNum1 > 0 and obsNum2 > 0 :
                OtOb += 1
            elif otherNum1 > 0 and otherNum2 > 0:
                OtOt += 1
            
        return {'R+':Rp,'Ob+':Obp,'Ot+':Otp,'R-':Rm,'Ob-':Obm,'Ot-':Otm, \
          'RR+-':RR,'ROb+-':ROb,'ROt+-':ROt, \
          'ObR+-':ObR,'ObOb+-':ObOb,'ObOt+-':ObOt, \
          'OtR+-':OtR,'OtOb+-':OtOb,'OtOt+-':OtOt}

    def __getReadBQProfile(self, ref, obs, Chr, pos, obsType, minBQ, readBuffer):
        ansDic = {'high':0,'low':0, 'avg':0}
        readNumAll = 0
        if obsType != 'M':
            return ansDic

        if readBuffer is None:
            return ansDic

        for read in readBuffer:
            mappedDic = self.__makeMappedDic(read)
            if (Chr,pos) in mappedDic:
                if mappedDic[(Chr,pos)]['cigar'] == 0:
                    readNumAll +=1
                    ansDic['avg'] += mappedDic[(Chr,pos)]['bq']

                    if mappedDic[(Chr,pos)]['bq'] < minBQ:
                        ansDic['low'] += 1
                    else :
                        ansDic['high'] += 1
        # logging.info(str(ansDic))
        if readNumAll > 0:
            ansDic['avg'] /= (1.0 * readNumAll)
        return ansDic

    def __getVariantReadBQProfile(self, ref, obs, Chr, pos, obsType, minBQ, readBuffer):
        ansDic = {'high':0,'low':0}
        if obsType != 'M':
            return ansDic

        if readBuffer is None:
            return ansDic

        for read in readBuffer:
            mappedDic = self.__makeMappedDic(read)
            if (Chr,pos) in mappedDic:
                if mappedDic[(Chr,pos)]['base'] == obs and mappedDic[(Chr,pos)]['cigar'] == 0:
                    if mappedDic[(Chr,pos)]['bq'] < minBQ:
                        ansDic['low'] += 1
                    else :
                        ansDic['high'] += 1
        logging.info(str(ansDic))

        return ansDic


    def __filterByVariantReadBQ(self, obsType, bqProfile, maxLowProportion, maxLowNum):
        if obsType != 'M':
            return True # passing this filter when I or D
        proportion = 0.0
        lowNum = bqProfile['low']
        highNum = bqProfile['high']
        if lowNum + highNum > 0:
            proportion = (1.0*lowNum)/(1.0*lowNum + 1.0*highNum)

        if not (proportion < maxLowProportion) or not( lowNum < maxLowNum ) :
            return False

        return True

    def __getState(self,ref,obs,Chr,pos,read,obsType,minBQ=None):
        mapppedDic = self.__makeMappedDic(read)
        isPlus = not(read.is_reverse)
        # if read.mapping_quality < self.__minMapQ:
        #     return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}
        if minBQ is None:
            minBQ = self.__minBQ
        if obsType == 'M':
            baseAtPos = None
            if (Chr,pos) in mapppedDic:
                if mapppedDic[(Chr,pos)]['base'] == '-':
                    if mapppedDic[(Chr,pos)]['cigar'] == 2:
                        # D
                        return {'ref':False,'obs':False,'other':True,'isPlus':isPlus}
                    else :
                        # other
                        return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}  
                if mapppedDic[(Chr,pos)]['bq'] < minBQ:
                    return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}
                insAns = self.__checkIns(ref,obs,Chr,pos,read,obsType,mapppedDic)
                delAns = self.__checkDel(ref,obs,Chr,pos,read,obsType,mapppedDic)
                if insAns or delAns:
                    return {'ref':False,'obs':False,'other':True,'isPlus':isPlus}

                baseAtPos = mapppedDic[(Chr,pos)]['base']
                baseAtPos = baseAtPos.upper()
                isref = (ref == baseAtPos)
                isobs = (obs == baseAtPos)
                if isref or isobs :
                    return {'ref':isref,'obs':isobs,'other':None,'isPlus':isPlus}
                else: 
                    return {'ref':isref,'obs':isobs,'other':baseAtPos,'isPlus':isPlus}
            else : # no mapping
                return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}
        elif obsType == 'I':
            if (Chr,pos) in mapppedDic:
                if mapppedDic[(Chr,pos)]['base'] == '-':
                    if mapppedDic[(Chr,pos)]['cigar'] == 2:
                        # D
                        return {'ref':False,'obs':False,'other':True,'isPlus':isPlus}
                    else :
                        # other
                        return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}  

                if mapppedDic[(Chr,pos)]['bq'] < minBQ:
                    return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}
                snvAns = self.__checkSNV(ref,obs,Chr,pos,read,obsType,mapppedDic)
                delAns = self.__checkDel(ref,obs,Chr,pos,read,obsType,mapppedDic)
                if snvAns or delAns:
                    return {'ref':False,'obs':False,'other':True,'isPlus':isPlus}

                if len(mapppedDic[(Chr,pos)]['list']) > 0:
                    for cElem in mapppedDic[(Chr,pos)]['list']:
                        if cElem[2] == 1: # cigar I
                            isobs = (cElem[0] == obs)
                            if isobs:
                                return {'ref':False,'obs':isobs,'other':None,'isPlus':isPlus}
                            else:
                                return {'ref':False,'obs':False,'other':cElem[0],'isPlus':isPlus}
                    return {'ref':True,'obs':False,'other':None,'isPlus':isPlus}
                else : # no Ins
                    return {'ref':True,'obs':False,'other':None,'isPlus':isPlus}
            else :
                return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}
        elif obsType == 'D':
            if (Chr,pos) in mapppedDic:
                if mapppedDic[(Chr,pos)]['base'] == '-':
                    if mapppedDic[(Chr,pos)]['cigar'] == 2:
                        # D
                        return {'ref':False,'obs':False,'other':True,'isPlus':isPlus}
                    else :
                        # other
                        return {'ref':False,'obs':False,'other':True,'isPlus':isPlus}
                        # return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}  
                if mapppedDic[(Chr,pos)]['bq'] < minBQ:
                    return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}

                snvAns = self.__checkSNV(ref,obs,Chr,pos,read,obsType,mapppedDic)
                insAns = self.__checkIns(ref,obs,Chr,pos,read,obsType,mapppedDic)
                if snvAns or insAns:
                    return {'ref':False,'obs':False,'other':True,'isPlus':isPlus}

                delList = []
                delLen = 0
                while (Chr,pos+1+delLen) in mapppedDic and mapppedDic[(Chr,pos+1+delLen)]['cigar'] == 2:
                    delList.append(mapppedDic[(Chr,pos+1+delLen)])
                    delLen += 1
                delString=''
                for Del in delList:
                    if Del['base'] == '-':
                        delString += Del['base']
                    else :
                        break

                if len(delString) > 0:
                    if len(delString) == len(obs):
                        return {'ref':False,'obs':True,'other':None,'isPlus':isPlus}
                    else :
                        return {'ref':False,'obs':False,'other':delString,'isPlus':isPlus}
                else:
                    return {'ref':True,'obs':False,'other':None,'isPlus':isPlus}
            else:
                return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}


    def __mutListToPysam(self,pos):
        return pos-1


    def __printOverlappingReads(self,TYPE,Chr,pos,ref,obs):
        print(TYPE+' '+str(Chr)+' '+str(pos)+' '+str(ref)+' '+str(obs))
        if self.__bam is None:
            raise InvalidBamFileManagement
        pos = int(pos)
        readList = {}
        readBuffer = []
        for read in self.__bam.fetch(Chr,max(pos-self.__windowSize,0),pos+self.__windowSize):
            readBuffer.append(read)
        readBuffer = self.__filterReads(readBuffer,TYPE,self.__chrDic[Chr],self.__mutListToPysam(pos),ref,obs)
        for read in readBuffer:
            if read.query_name not in readList:
                readList[read.query_name] = []
            readList[read.query_name].append(read)

        if readBuffer is None:
            print("No reads")
        else :
            totalOverlap = 0
            for ID in readList:
                if len(readList[ID]) == 2:
                    reads = []
                    if readList[ID][1].is_reverse:
                        dics =[self.__makeMappedDic(readList[ID][0]),self.__makeMappedDic(readList[ID][1])]
                        reads = [readList[ID][0],readList[ID][1]]
                    else:
                        dics =[self.__makeMappedDic(readList[ID][1]),self.__makeMappedDic(readList[ID][0])]
                        reads = [readList[ID][1],readList[ID][0]]
                    
                    coverNum = 0
                    basesList = []
                    bqList = []
                    mappedLengths = []
                    foundPosition = []
                    foundBases = []
                    for dic in dics:
                        mappedLengths.append(len(dic.keys()))
                        readString = ''
                        bqString = ''
                        cycle = 0
                        for at in sorted(dic.keys()):
                            if at[1] == pos-1:
                                foundPosition.append(cycle)
                                coverNum += 1
                                readString += " "
                                bqString += " "
                                if len(dic[at]['list']) > 0 or  dic[at]['bq'] is None:
                                    readString += '^'
                                    foundBases.append('^')
                                else:
                                    readString += dic[at]['base']
                                    foundBases.append(dic[at]['base'])

                                if dic[at]['base'] == '-' or dic[at]['bq'] is None:
                                    bqString += '-'
                                else:
                                    bqString += chr(dic[at]['bq']+33)
                                readString += " "
                                bqString += " "
                            else :
                                if len(dic[at]['list']) > 0 or  dic[at]['bq'] is None:
                                    readString += '^'
                                else:
                                    readString += dic[at]['base']
                                
                                if dic[at]['base'] == '-' or dic[at]['bq'] is None:
                                    bqString += '-'
                                else:
                                    bqString += chr(dic[at]['bq']+33)
                            cycle += 1
                        basesList.append(readString)
                        bqList.append(bqString)
                    if coverNum == 2:
                        totalOverlap +=1
                        if (TYPE == 'M' and foundBases[0] == obs and foundBases[1] == obs):
                            print("")
                            print(reads[0].template_length)
                            print(basesList[0])
                            print(bqList[0])
                            print(basesList[1])
                            print(bqList[1])                            
                        
    def getOverlapInformation(self,TYPE,Chr,pos,ref,obs,lowBQ=None):
        if self.__bam is None:
            raise InvalidBamFileManagement
        pos = int(pos)
        obsType = TYPE
        pileAns = None
        pileAnsLowBQ = None
        keys = None

        if obsType == 'M' or obsType == 'I' or obsType == 'D':
            pileAns = {'R+':0,'Ob+':0,'Ot+':0,'R-':0,'Ob-':0,'Ot-':0, \
                      'RR+-':0,'ROb+-':0,'ROt+-':0, \
                      'ObR+-':0,'ObOb+-':0,'ObOt+-':0, \
                      'OtR+-':0,'OtOb+-':0,'OtOt+-':0}
            pileAnsLowBQ = {'R+':0,'Ob+':0,'Ot+':0,'R-':0,'Ob-':0,'Ot-':0, \
                      'RR+-':0,'ROb+-':0,'ROt+-':0, \
                      'ObR+-':0,'ObOb+-':0,'ObOt+-':0, \
                      'OtR+-':0,'OtOb+-':0,'OtOt+-':0}
            keys = ['R+','Ob+','Ot+','R-','Ob-','Ot-', \
                      'RR+-','ROb+-','ROt+-', \
                      'ObR+-','ObOb+-','ObOt+-', \
                      'OtR+-','OtOb+-','OtOt+-']
            readList = {}
            readBuffer = []
            for read in self.__bam.fetch(Chr,max(pos-self.__windowSize,0),pos+self.__windowSize):
                readBuffer.append(read)
            readBuffer = self.__filterReads(readBuffer,TYPE,self.__chrDic[Chr],self.__mutListToPysam(pos),ref,obs)
            if readBuffer is None:
                if lowBQ is None:
                    return (None,None)
                else:
                    return (None,None,None)
            
            for read in readBuffer:
                if read.query_name not in readList:
                    readList[read.query_name] = []
                readList[read.query_name].append(read)

            for ID in readList:
                if len(readList[ID]) == 2:
                    pairAns = self.__getState2(ref,obs,self.__chrDic[Chr],\
                        self.__mutListToPysam(pos),readList[ID][0],readList[ID][1],obsType)
                    for key in keys:
                        pileAns[key] += pairAns[key]
                elif len(readList[ID]) == 1:
                    singleAns = self.__getState(ref,obs,self.__chrDic[Chr],\
                        self.__mutListToPysam(pos),readList[ID][0],obsType)
                    plusStr = None
                    if singleAns['isPlus']:
                        plusStr = '+'
                    else:
                        plusStr = '-'

                    if singleAns['ref']:
                        pileAns['R'+plusStr] += 1
                    elif singleAns['obs']:
                        pileAns['Ob'+plusStr] += 1
                    elif singleAns['other'] is not None:
                        pileAns['Ot'+plusStr] += 1

                else :
                    raise InsufficientBamReadFilter

            if lowBQ is not None:
                for ID in readList:
                    if len(readList[ID]) == 2:
                        pairAnsLowBQ = self.__getState2(ref,obs,self.__chrDic[Chr],\
                            self.__mutListToPysam(pos),readList[ID][0],readList[ID][1],obsType,lowBQ)
                        for key in keys:
                            pileAnsLowBQ[key] += pairAnsLowBQ[key]

                    elif len(readList[ID]) == 1:
                        singleAnsLowBQ = self.__getState(ref,obs,self.__chrDic[Chr],\
                            self.__mutListToPysam(pos),readList[ID][0],obsType,lowBQ)
                        plusStr = None
                        if singleAnsLowBQ['isPlus']:
                            plusStr = '+'
                        else:
                            plusStr = '-'

                        if singleAnsLowBQ['ref']:
                            pileAnsLowBQ['R'+plusStr] += 1
                        elif singleAnsLowBQ['obs']:
                            pileAnsLowBQ['Ob'+plusStr] += 1
                        elif singleAnsLowBQ['other'] is not None:
                            pileAnsLowBQ['Ot'+plusStr] += 1
                    else :
                        raise InsufficientBamReadFilter

            headCols = [obsType,Chr,pos,ref,obs]
            bodyCols = []
            bodyColsLowBQ = []

            for key in keys:
                bodyCols.append(pileAns[key])
                if lowBQ is not None:
                    bodyColsLowBQ.append(pileAnsLowBQ[key]) 
            if lowBQ is None:
                return (headCols,bodyCols)
            else:
                return [headCols,bodyCols,bodyColsLowBQ]
        else :
            raise InvalidMutationType

    def getReadsAround(self,TYPE,Chr,pos,ref,obs,ignoreObs=False):
        if self.__bam is None:
            raise InvalidBamFileManagement
        pos = int(pos)
        obsType = TYPE
        if obsType == 'M' or obsType == 'I' or obsType == 'D':
            # readList = {}
            readBuffer = []
            ansList = []
            for read in self.__bam.fetch(Chr,max(pos-self.__windowSize,0),pos+self.__windowSize):
                readBuffer.append(read)

            for read in readBuffer:
                singleAns = self.__getState(ref,obs,self.__chrDic[Chr],self.__mutListToPysam(pos),read,obsType)
                if not ignoreObs:
                    ansList.append(read)
                elif ignoreObs and ( singleAns['obs'] == 0 ):
                    ansList.append(read)
            return ansList






