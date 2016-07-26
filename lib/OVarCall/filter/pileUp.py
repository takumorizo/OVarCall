#!/usr/bin/env python
# data for one pileUp : base,Ins,Del with Basequality
from OVarCall.utilOVar.exceptionUtil import InvalidPileupReference, InvalidPileupFormat
import logging


class PileUpUnit:

    def __init__(self, _TYPE, _obs, _strand, _quality):
        if _TYPE == 'M' or _TYPE == 'I' or _TYPE == 'D':
            self.__base_quality_offset = 33
            self.TYPE = _TYPE
            self.obs = _obs.upper()
            self.strand = _strand
            self.quality = self.__basequality(_quality)
            if _TYPE == 'I' or _TYPE == 'D':
                self.quality = None

        else:
            raise Exception('unexpected initialization of PileUpUnit')

    def __basequality(self, _quality):
        if type(_quality) is str and len(_quality) == 1:
            return ord(_quality) - self.__base_quality_offset


class PileUp:

    def __init__(self, chromosome, position, ref, depth, bases, qualities):
        self.chromosome = chromosome
        self.position = int(position)
        self.ref = ref.upper()
        self.depth = int(depth)
        self.bases = bases
        self.qualities = qualities
        self.summary = None
        self.__pileUpUnits = None
        self.__simpleSummary()

    def setLine(self, chromosome, position, ref, depth, bases, qualities):
        self.__pileUpUnits = None
        self.summary = None
        self.detail = None
        self.chromosome = chromosome
        self.position = int(position)
        self.ref = ref.upper()
        self.depth = int(depth)
        self.bases = bases
        self.qualities = qualities
        self.__simpleSummary()

    def detailedSummary(self, minBQ=None):
        if self.ref not in 'ACGTNacgtn':
            raise InvalidPileupReference
        self.__correspondBaseWithBQ()
        detail = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0,
                  'a': 0, 'c': 0, 'g': 0, 't': 0, 'n': 0, 'I': {}, 'D': {}}
        # 'I':{ 'ACGTAC...' : { '+':num_plus, '-':num_minus } }
        # 'D':{ 'ACGTAC...' : { '+':num_plus, '-':num_minus } }
        detail['depth'] = self.depth
        for unit in self.__pileUpUnits:
            if minBQ is not None and type(minBQ) is int:
                if unit.quality is not None and unit.quality < minBQ:
                    continue
            if unit.TYPE == 'M':
                if unit.strand == 1:
                    detail[unit.obs.upper()] += 1
                else:
                    detail[unit.obs.lower()] += 1
            else:
                # unit.TYPE == 'I' or 'D'
                strandStr = '+'
                if unit.strand == 0:
                    strandStr = '-'
                if unit.obs.upper() not in detail[unit.TYPE]:
                    detail[unit.TYPE][unit.obs.upper()] = {'+': 0, '-': 0}
                detail[unit.TYPE][unit.obs.upper()][strandStr] += 1
        return detail

    def __correspondBaseWithBQ(self):
        if self.__pileUpUnits is not None:
            return
        self.__pileUpUnits = []
        if self.depth == 0:
            return
        b = 0  # index for base
        q = 0  # index for base quality
        while b < len(self.bases):
            if self.bases[b] == '.':
                self.__pileUpUnits.append(PileUpUnit('M', self.ref.upper(), 1, self.qualities[q]))
                b += 1
                q += 1
            elif self.bases[b] == ',':
                self.__pileUpUnits.append(PileUpUnit('M', self.ref.upper(), 0, self.qualities[q]))
                b += 1
                q += 1
            elif self.bases[b] in "ACGTN":
                self.__pileUpUnits.append(
                    PileUpUnit('M', self.bases[b].upper(), 1, self.qualities[q]))
                b += 1
                q += 1
            elif self.bases[b] in "acgtn":
                self.__pileUpUnits.append(
                    PileUpUnit('M', self.bases[b].upper(), 0, self.qualities[q]))
                b += 1
                q += 1
            elif self.bases[b] == '$':
                b += 1
            elif self.bases[b] == '^':
                b += 2  # skip mapping quality information
            elif self.bases[b] == '*':
                # must skip deletion information
                b += 1
                q += 1
            elif self.bases[b] == '<' or self.bases[b] == '>':
                # skip reference skip such as intronic region.
                b += 1
                q += 1
            elif self.bases[b] in "+-":
                # Not removing base before Ins/Del
                # calculating Ins/Del rate by using Ins/Del and Depth.
                # Because bases before Ins/Dels are not removed, I cannot calculate
                # Ins/Del by bases count and Ins/Del
                tmp_b = b + 1
                lenStr = ''
                while ord('0') <= ord(self.bases[tmp_b]) <= ord('9'):
                    lenStr += self.bases[tmp_b]
                    tmp_b += 1
                insDelLen = int(lenStr)
                obs = self.bases[tmp_b:tmp_b + insDelLen]
                strand = None
                if obs[0] in "ACGTN":
                    strand = 1
                else:
                    strand = 0
                if self.bases[b] == '+':
                    self.__pileUpUnits.append(PileUpUnit('I', obs.upper(), strand, None))
                else:
                    self.__pileUpUnits.append(PileUpUnit('D', obs.upper(), strand, None))
                b = tmp_b + insDelLen
            else:
                raise InvalidPileupFormat("unexpected format of pileup")
        if b != len(self.bases) or q != len(self.qualities):
            raise InvalidPileupFormat("unexpected piled parse result")

    def __simpleSummary(self):

        self.summary = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0,
                                'a': 0, 'c': 0, 'g': 0, 't': 0, 'n': 0, '+': 0, '-': 0}
        if self.ref.upper() in self.summary:
            self.summary[self.ref.upper()] += self.bases.count('.')
        if self.ref.lower() in self.summary:
            self.summary[self.ref.lower()] += self.bases.count(',')

        self.summary['+'] += self.bases.count('+')
        self.summary['-'] += self.bases.count('-')
