#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys


def makeMappedDic(alignedSeg):
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
        mlen = MLenOne(firstCigar)
        clen = CLenOne(firstCigar)
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
        elif mlen > 0 and clen == 0:  # D,N,P
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

def MLenOne(cigar):
    try:
        ctype = cigar[0]
        # return ctype == M or ctype == D or ctype == N or EQ or X,P
        if ctype == 0 or ctype == 2 or ctype == 3 or ctype == 7 or ctype == 8 or ctype == 6:
            return cigar[1]
        else:
            return 0
    except Exception, e:
        sys.stderr.writelines("isMlenOne invalid form of cigar")
        raise e


def CLenOne(cigar):
    try:
        ctype = cigar[0]
        # return ctype == M or ctype == I or ctype == S or EQ or X
        if ctype == 0 or ctype == 1 or ctype == 4 or ctype == 7 or ctype == 8:
            return cigar[1]
        else:
            return 0
    except Exception, e:
        sys.stderr.writelines("isCLenOne invalid form of cigar")
        raise e



def getState2(ref,obs,Chr,pos,read1,read2,obsType,minBQ,minMapQ):
    if obsType == 'M' :
        others = []
        for base in ['A','C','G','T']:
            if base != ref.upper() and base != obs.upper():
                others.append(base)

        if not read1.is_reverse:
            state1 = getState(ref,obs,Chr,pos,read1,obsType,minBQ,minMapQ)
            state2 = getState(ref,obs,Chr,pos,read2,obsType,minBQ,minMapQ)
        else :
            state1 = getState(ref,obs,Chr,pos,read2,obsType,minBQ,minMapQ)
            state2 = getState(ref,obs,Chr,pos,read1,obsType,minBQ,minMapQ)
        refNum1 = 0
        obsNum1 = 0
        other1Num1 = 0
        other2Num1 = 0
        
        refNum2 = 0
        obsNum2 = 0
        other1Num2 = 0
        other2Num2 = 0
        
        if state1['ref']:
            refNum1 += 1
        if state1['obs']:
            obsNum1 += 1
        if state1['other'] is not None:
            if others[0] == state1['other'].upper():
                other1Num1 += 1
            if others[1] == state1['other'].upper():
                other2Num1 += 1
        if state2['ref']:
            refNum2 += 1
        if state2['obs']:
            obsNum2 += 1
        if state2['other'] is not None:
            if others[0] == state2['other'].upper():
                other1Num2 += 1
            if others[1] == state2['other'].upper():
                other2Num2 += 1

        Rp = refNum1
        Obp = obsNum1
        Ot1p = other1Num1
        Ot2p = other2Num1
        Rm = refNum2
        Obm = obsNum2
        Ot1m = other1Num2
        Ot2m = other2Num2

        RR = 0
        ROb = 0
        ROt1 = 0
        ROt2 = 0

        ObR = 0
        ObOb = 0
        ObOt1 = 0
        ObOt2 = 0
        
        Ot1R = 0
        Ot1Ob = 0
        Ot1Ot1 = 0
        Ot1Ot2 = 0
        
        Ot2R = 0
        Ot2Ob = 0
        Ot2Ot1 = 0
        Ot2Ot2 = 0

        if (refNum1+obsNum1+other1Num1+other2Num1) == 1 and (refNum2+obsNum2+other1Num2+other2Num2) == 1:
            # Not counting as single read
            Rp = 0
            Obp = 0
            Ot1p = 0
            Ot2p = 0
            Rm = 0
            Obm = 0
            Ot1m = 0
            Ot2m = 0
            if refNum1 > 0 and refNum2 > 0 :
                RR += 1
            elif refNum1 > 0 and obsNum2 > 0:
                ROb += 1
            elif refNum1 > 0 and other1Num2 > 0:
                ROt1 += 1
            elif refNum1 > 0 and other2Num2 > 0:
                ROt2 += 1    
            elif obsNum1 > 0 and refNum2 > 0 :
                ObR += 1
            elif obsNum1 > 0 and obsNum2 > 0:
                ObOb += 1
            elif obsNum1 > 0 and other1Num2 > 0 :
                ObOt1 += 1
            elif obsNum1 > 0 and other2Num2 > 0 :
                ObOt2 += 1
            elif other1Num1 > 0 and refNum2 > 0 :
                Ot1R += 1
            elif other1Num1 > 0 and obsNum2 > 0 :
                Ot1Ob += 1
            elif other1Num1 > 0 and other1Num2 > 0:
                Ot1Ot1 += 1
            elif other1Num1 > 0 and other2Num2 > 0:
                Ot1Ot2 += 1
            elif other2Num1 > 0 and refNum2 > 0 :
                Ot2R += 1
            elif other2Num1 > 0 and obsNum2 > 0 :
                Ot2Ob += 1
            elif other2Num1 > 0 and other1Num2 > 0:
                Ot2Ot1 += 1
            elif other2Num1 > 0 and other2Num2 > 0:
                Ot2Ot2 += 1

        return {'R+':Rp,'Ob+':Obp,'Ot1+':Ot1p,'Ot2+':Ot2p,'R-':Rm,'Ob-':Obm,'Ot1-':Ot1m,'Ot2-':Ot2m, \
          'RR+-':RR,'ROb+-':ROb,'ROt1+-':ROt1,'ROt2+-':ROt2, \
          'ObR+-':ObR,'ObOb+-':ObOb,'ObOt1+-':ObOt1,'ObOt2+-':ObOt2, \
          'Ot1R+-':Ot1R,'Ot1Ob+-':Ot1Ob,'Ot1Ot1+-':Ot1Ot1,'Ot1Ot2+-':Ot1Ot2,\
          'Ot2R+-':Ot2R,'Ot2Ob+-':Ot2Ob,'Ot2Ot1+-':Ot2Ot1,'Ot2Ot2+-':Ot2Ot2}
    elif obsType == 'I' or 'D' :
        if not read1.is_reverse:
            state1 = getState(ref,obs,Chr,pos,read1,obsType,minBQ,minMapQ)
            state2 = getState(ref,obs,Chr,pos,read2,obsType,minBQ,minMapQ)
        else :
            state1 = getState(ref,obs,Chr,pos,read2,obsType,minBQ,minMapQ)
            state2 = getState(ref,obs,Chr,pos,read1,obsType,minBQ,minMapQ)
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


"""
getState(ref,obs,Chr,pos,read,obsType):
の pos -> 0-indexed

'M' : pos の場所(0-indexed pysam 上) を
'I' : pos の場所(0-indexed pysam 上) の mappedList の ins を見る
'D' : 0-indexed pysam 上 の [pos,pos+len) を見る
"""

def getState(ref,obs,Chr,pos,read,obsType,minBQ,minMapQ):
    mappedDic = makeMappedDic(read)
    isPlus = not(read.is_reverse)
    if read.mapping_quality < minMapQ:
        return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}

    if obsType == 'M':
        baseAtPos = None
        if (Chr,pos) in mappedDic:
            if mappedDic[(Chr,pos)]['bq'] < minBQ:
                return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}
            baseAtPos = mappedDic[(Chr,pos)]['base']
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
        if (Chr,pos) in mappedDic:
            if mappedDic[(Chr,pos)]['bq'] < minBQ:
                return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}
            if len(mappedDic[(Chr,pos)]['list']) > 0:
                for cElem in mappedDic[(Chr,pos)]['list']:
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
        if (Chr,pos) in mappedDic:
            if mappedDic[(Chr,pos)]['bq'] < minBQ:
                return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}
            delList = []
            delLen = 0
            while (Chr,pos+1+delLen) in mappedDic:
                delList.append(mappedDic[(Chr,pos+1+delLen)])
                delLen += 1
            delString=''
            baseString=''
            for Del in delList:
                if Del['base'] == '-':
                    delString += Del['base']
                else :
                    break
            for Base in delList:
                if Del['base'] != '-':
                    baseString += Base['base']
                else :
                    break
            if len(delString) > 0:
                if len(delString) == len(obs):
                    return {'ref':False,'obs':True,'other':None,'isPlus':isPlus}
                else :
                    return {'ref':False,'obs':False,'other':delString,'isPlus':isPlus}
            else:
                if len(baseString) >= len(obs):
                    return {'ref':True,'obs':False,'other':None,'isPlus':isPlus}
                else:
                    return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}
        else:
            return {'ref':False,'obs':False,'other':None,'isPlus':isPlus}


def mutListToPysam(pos):
    return pos-1