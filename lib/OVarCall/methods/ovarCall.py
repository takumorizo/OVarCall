#!/usr/bin/env python
from OVarCall.utilOVar.exceptionUtil import InvalidBamFileManagement,InvalidMutationType,InsufficientBamReadFilter
from OVarCall.utilOVar.readUtil import getState,getState2,mutListToPysam 
from OVarCall.utilOVar.bamDAO import BamDAO
from scipy.special import digamma
from scipy.special import betaln
from scipy.special import gammaln
import math
import distutils.util
import logging

class OVarCall():
    def __init__(self,settings=None):
        self.__MAX_UPDATE_COUNT = 500
        self.__COUNVERGE_DIFF = 1e-4

        # setting parameters for calculaiton of Bayes factor
        self.__DEL_P1_TUMOR = 2.0
        self.__DEL_P2_TUMOR = 30.0
        self.__DEL_M1_TUMOR = 2.0
        self.__DEL_M2_TUMOR = 30.0

        self.__DEL_P1_ERROR = 2.0
        self.__DEL_P2_ERROR = 30.0
        self.__DEL_M1_ERROR = 2.0
        self.__DEL_M2_ERROR = 30.0
        
        self.__GAMMA_TUMOR = [1.0, 1.0]
        self.__GAMMA_ERROR = [1000000.0, 10.0]
        self.setParameters(settings)
        

    def setParameters(self,settings):
        if type(settings) is dict :
            if 'MAX_UPDATE_COUNT' in settings:
                self.__MAX_UPDATE_COUNT = int(settings['MAX_UPDATE_COUNT'])
            if 'COUNVERGE_DIFF' in settings:
                self.__COUNVERGE_DIFF = float(settings['COUNVERGE_DIFF'])

            if 'DEL_P1_TUMOR' in settings:
                self.__DEL_P1_TUMOR = float(settings['DEL_P1_TUMOR'])
            if 'DEL_P2_TUMOR' in settings:
                self.__DEL_P2_TUMOR = float(settings['DEL_P2_TUMOR'])

            if 'DEL_M1_TUMOR' in settings:
                self.__DEL_M1_TUMOR = float(settings['DEL_M1_TUMOR'])
            if 'DEL_M2_TUMOR' in settings:
                self.__DEL_M2_TUMOR = float(settings['DEL_M2_TUMOR'])

            if 'DEL_P1_ERROR' in settings:
                self.__DEL_P1_ERROR = float(settings['DEL_P1_ERROR'])
            if 'DEL_P2_ERROR' in settings:
                self.__DEL_P2_ERROR = float(settings['DEL_P2_ERROR'])

            if 'DEL_M1_ERROR' in settings:
                self.__DEL_M1_ERROR = float(settings['DEL_M1_ERROR'])
            if 'DEL_M2_ERROR' in settings:
                self.__DEL_M2_ERROR = float(settings['DEL_M2_ERROR'])

            if 'GAMMA_TUMOR_ALPHA' in settings:
                self.__GAMMA_TUMOR[0] = float(settings['GAMMA_TUMOR_ALPHA'])
            if 'GAMMA_TUMOR_BETA' in settings:
                self.__GAMMA_TUMOR[1] = float(settings['GAMMA_TUMOR_BETA'])
            
            if 'GAMMA_ERROR_ALPHA' in settings:
                self.__GAMMA_ERROR[0] = float(settings['GAMMA_ERROR_ALPHA'])
            if 'GAMMA_ERROR_BETA' in settings:
                self.__GAMMA_ERROR[1] = float(settings['GAMMA_ERROR_BETA'])

            if 'useNormal' in settings:
                self.__useNormal = bool( distutils.util.strtobool( str(settings['useNormal']) ) )

    def doInference(self, TYPE, Chr, pos, ref, obs, head,bodyT,bodyN=None):
        countColFromT = None
        countColFromN = None
        countColFromT = 5
        countColFromN = 21
        if TYPE not in 'MID' :
            raise InvalidMutationType
        
        if bodyN is None:
            result = head
            result.extend(bodyT)
            result.extend(self.__startCalculationT(0,countColFromT,result))
            result.extend(self.__startCalculationT(0,countColFromT,result,True))
            return result
        else :
            result = []
            result = head
            result.extend(bodyT)
            result.append('Normal:')
            result.extend(bodyN)
            if bodyN is not None:
                result.extend(self.__startCalculationTN(0,countColFromT,countColFromN,result))
                result.extend(self.__startCalculationTN(0,countColFromT,countColFromN,result,True))
            else :
                result.extend(self.__startCalculationT(0,countColFromT,result))
                result.extend(self.__startCalculationT(0,countColFromT,result,True))
            return result   

    def __startCalculationT(self,mutCol,countColFromTumor,inputCols,ignoreOverlap=False):
        mutType = inputCols[mutCol]
        tumorCols = None
        tumorCols = inputCols[countColFromTumor:countColFromTumor+15]
        tumorCols = map(int,tumorCols)
        if ignoreOverlap:
            tumorCols[6:15] = [0]*9
        ErrorLq = self.__baseModelsLowerBound(tumorCols[0],tumorCols[1],tumorCols[2],\
            tumorCols[3],tumorCols[4],tumorCols[5],\
            tumorCols[6],tumorCols[7],tumorCols[8],\
            tumorCols[9],tumorCols[10],tumorCols[11],\
            tumorCols[12],tumorCols[13],tumorCols[14],\
            self.__DEL_P1_ERROR,self.__DEL_P2_ERROR,self.__DEL_M1_ERROR,self.__DEL_M2_ERROR\
            ,self.__GAMMA_ERROR)['Lq']
        TumorLq = self.__baseModelsLowerBound(tumorCols[0],tumorCols[1],tumorCols[2],\
            tumorCols[3],tumorCols[4],tumorCols[5],\
            tumorCols[6],tumorCols[7],tumorCols[8],\
            tumorCols[9],tumorCols[10],tumorCols[11],\
            tumorCols[12],tumorCols[13],tumorCols[14],\
            self.__DEL_P1_TUMOR,self.__DEL_P2_TUMOR,self.__DEL_M1_TUMOR,self.__DEL_M2_TUMOR,\
            self.__GAMMA_TUMOR)['Lq']

        outputCols = []
        # outputCols.extend(inputCols)
        if ErrorLq >= TumorLq:
            outputCols.append('E')
        else :
            outputCols.append('T')
        outputCols.append(math.log10(math.e)*ErrorLq)
        outputCols.append(math.log10(math.e)*TumorLq)
        outputCols.append(math.log10(math.e)*(TumorLq-ErrorLq))
        return outputCols
        if mutType not in 'MID' :
            raise InvalidMutationType

    def __startCalculationTN(self,mutCol,countColFromTumor,countColFromNormal,inputCols,ignoreOverlap=False):
        mutType = inputCols[mutCol]
        if mutType not in 'MID':
            raise InvalidMutationType

        tumorCols = None
        normalCols = None
        # print(inputCols)
        normalCols = inputCols[countColFromNormal:countColFromNormal+15]
        normalCols = map(int,normalCols)
        tumorCols = inputCols[countColFromTumor:countColFromTumor+15]
        tumorCols = map(int,tumorCols)
        gam_Error = [self.__GAMMA_ERROR[0],self.__GAMMA_ERROR[1]]
        gam_Tumor = [self.__GAMMA_TUMOR[0],self.__GAMMA_TUMOR[1]]

        if gam_Error[0] < 0 or gam_Error[1] < 0:
            normalRef = normalCols[0]+normalCols[3]+normalCols[6]
            normalObs = normalCols[1]+normalCols[4]+normalCols[10]
            gam_Error=[normalRef+1,normalObs+1]
            logging.info("specifed parameters of GAMMA_ERROR is not good parameters for beta distribution")
        
        if gam_Tumor[0] < 0 or gam_Tumor[1] < 0:
            tumorRef = tumorCols[0]+tumorCols[3]+tumorCols[6]
            tumorObs = tumorCols[1]+tumorCols[4]+tumorCols[10]
            gam_Tumor=[tumorRef+1,tumorObs+1]
            logging.info("specifed parameters of GAMMA_TUMOR is not good parameters for beta distribution")


        if ignoreOverlap:
            normalCols[6:15] = [0]*9
            tumorCols[6:15] = [0]*9
        tumorError = self.__baseModelsLowerBound(normalCols[0],normalCols[1],normalCols[2],\
            normalCols[3],normalCols[4],normalCols[5],\
            normalCols[6],normalCols[7],normalCols[8],\
            normalCols[9],normalCols[10],normalCols[11],\
            normalCols[12],normalCols[13],normalCols[14],\
            self.__DEL_P1_TUMOR,self.__DEL_P2_TUMOR,self.__DEL_M1_TUMOR,self.__DEL_M2_TUMOR,\
            gam_Error)
        normalError = self.__baseModelsLowerBound(normalCols[0],normalCols[1],normalCols[2],\
            normalCols[3],normalCols[4],normalCols[5],\
            normalCols[6],normalCols[7],normalCols[8],\
            normalCols[9],normalCols[10],normalCols[11],\
            normalCols[12],normalCols[13],normalCols[14],\
            self.__DEL_P1_ERROR,self.__DEL_P2_ERROR,self.__DEL_M1_ERROR,self.__DEL_M2_ERROR,\
            gam_Error) 
        ErrorLq = self.__baseModelsLowerBound(tumorCols[0],tumorCols[1],tumorCols[2],\
            tumorCols[3],tumorCols[4],tumorCols[5],\
            tumorCols[6],tumorCols[7],tumorCols[8],\
            tumorCols[9],tumorCols[10],tumorCols[11],\
            tumorCols[12],tumorCols[13],tumorCols[14],\
            normalError['eps+'][0],normalError['eps+'][1],normalError['eps-'][0],normalError['eps-'][1]\
            ,gam_Error)['Lq']
        TumorLq = self.__baseModelsLowerBound(tumorCols[0],tumorCols[1],tumorCols[2],\
            tumorCols[3],tumorCols[4],tumorCols[5],\
            tumorCols[6],tumorCols[7],tumorCols[8],\
            tumorCols[9],tumorCols[10],tumorCols[11],\
            tumorCols[12],tumorCols[13],tumorCols[14],\
            tumorError['eps+'][0],tumorError['eps+'][1],tumorError['eps-'][0],tumorError['eps-'][1],\
            gam_Tumor)['Lq']
        outputCols = []
        if ErrorLq >= TumorLq:
            outputCols.append('E')
        else :
            outputCols.append('T')
        outputCols.append(math.log10(math.e)*ErrorLq)
        outputCols.append(math.log10(math.e)*TumorLq)
        outputCols.append(math.log10(math.e)*(TumorLq-ErrorLq))
        return outputCols



    def __baseModelsLowerBound(self,R_p, Ob_p, Oth1_p, R_m, Ob_m, Oth1_m,
                              RR, ROb, ROth1,
                              ObR, ObOb, ObOth1,
                              Oth1R, Oth1Ob, Oth1Oth1,
                              del_p1, del_p2, del_m1, del_m2,
                              gamma_pi):
        ### make RpList,EZ_p ###
        RpList = []
        EZ_p = []
        for i in range(R_p):
            RpList.append((1, 0, 0))
            EZ_p.append([1.0, 0.0])
        for i in range(Ob_p):
            RpList.append((0, 1, 0))
            EZ_p.append([1.0, 0.0])
        for i in range(Oth1_p):
            RpList.append((0, 0, 1))
            EZ_p.append([1.0, 0.0])

        ### make RmList,EZ_m ###
        RmList = []
        EZ_m = []
        for i in range(R_m):
            RmList.append((1, 0, 0))
            EZ_m.append([1.0, 0.0])
        for i in range(Ob_m):
            RmList.append((0, 1, 0))
            EZ_m.append([1.0, 0.0])
        for i in range(Oth1_m):
            RmList.append((0, 0, 1))
            EZ_m.append([1.0, 0.0])

        ### make RopList,RomList EZ_op,EZ_om ###
        RopList = []
        RomList = []
        EZ_o = []


        for i in range(RR):
            RopList.append((1, 0, 0))
            RomList.append((1, 0, 0))
            EZ_o.append([1.0, 0.0])
        for i in range(ROb):
            RopList.append((1, 0, 0))
            RomList.append((0, 1, 0))
            EZ_o.append([1.0, 0.0])
        for i in range(ROth1):
            RopList.append((1, 0, 0))
            RomList.append((0, 0, 1))
            EZ_o.append([1.0, 0.0])

        for i in range(ObR):
            RopList.append((0, 1, 0))
            RomList.append((1, 0, 0))
            EZ_o.append([1.0, 0.0])
        for i in range(ObOb):
            RopList.append((0, 1, 0))
            RomList.append((0, 1, 0))
            EZ_o.append([1.0, 0.0])
        for i in range(ObOth1):
            RopList.append((0, 1, 0))
            RomList.append((0, 0, 1))
            EZ_o.append([1.0, 0.0])

        for i in range(Oth1R):
            RopList.append((0, 0, 1))
            RomList.append((1, 0, 0))
            EZ_o.append([1.0, 0.0])
        for i in range(Oth1Ob):
            RopList.append((0, 0, 1))
            RomList.append((0, 1, 0))
            EZ_o.append([1.0, 0.0])
        for i in range(Oth1Oth1):
            RopList.append((0, 0, 1))
            RomList.append((0, 0, 1))
            EZ_o.append([1.0, 0.0])


        N_p = len(EZ_p)
        N_o = len(EZ_o)
        N_m = len(EZ_m)

        Lq_new = 10000000000
        Lq_old = 1
        new_gamma = [0.0, 0.0]
        new_alpha_p = [0.0, 0.0]
        new_alpha_m = [0.0, 0.0]

        updateCount = 0
        """ start variational update iteration """
        while (Lq_old == 0) or abs(Lq_new - Lq_old)/abs(Lq_old) > self.__COUNVERGE_DIFF:
            if updateCount > self.__MAX_UPDATE_COUNT:
                continue

            if Lq_new < Lq_old and Lq_new <= 0 and  Lq_old <= 0 :
                logging.error("##### decrease Lq ##### : " + str(Lq_new - Lq_old) +'\n')

            """ M step Calculation """
            for j in range(2):
                new_gamma[j] = gamma_pi[j]
                for n in range(N_p):
                    new_gamma[j] += EZ_p[n][j]
                for n in range(N_m):
                    new_gamma[j] += EZ_m[n][j]
                for n in range(N_o):
                    new_gamma[j] += EZ_o[n][j]

            # new_alpha_p[0] = del_p1 + N_p + N_o + 2
            new_alpha_p[0] = del_p1 + N_p + N_o
            new_alpha_p[1] = del_p2
            # new_alpha_m[0] = del_m1 + N_m + N_o + 2
            new_alpha_m[0] = del_m1 + N_m + N_o
            new_alpha_m[1] = del_m2
            for n in range(N_p):
                for j in range(2):
                    new_alpha_p[0] -= EZ_p[n][j] * RpList[n][j]
                    new_alpha_p[1] += EZ_p[n][j] * RpList[n][j]

            for n in range(N_o):
                for j in range(2):
                    new_alpha_p[0] -= EZ_o[n][j] * RopList[n][j]
                    new_alpha_p[1] += EZ_o[n][j] * RopList[n][j]

                    new_alpha_m[0] -= EZ_o[n][j] * RomList[n][j]
                    new_alpha_m[1] += EZ_o[n][j] * RomList[n][j]

            for n in range(N_m):
                for j in range(2):
                    new_alpha_m[0] -= EZ_m[n][j] * RmList[n][j]
                    new_alpha_m[1] += EZ_m[n][j] * RmList[n][j]

            dig_new_gamma_sum = digamma(new_gamma[0] + new_gamma[1])
            dig_new_gamma_list = [digamma(new_gamma[0]), digamma(new_gamma[1])]
            dig_new_alpha_p_sum = digamma(new_alpha_p[0] + new_alpha_p[1])
            dig_new_alpha_p_list = [digamma(new_alpha_p[0]), digamma(new_alpha_p[1])]

            dig_new_alpha_m_sum = digamma(new_alpha_m[0] + new_alpha_m[1])
            dig_new_alpha_m_list = [digamma(new_alpha_m[0]), digamma(new_alpha_m[1])]

            log3 = math.log(2)
            A_pi = [math.exp(digamma(new_gamma[0]) - dig_new_gamma_sum), math.exp(digamma(new_gamma[1]) - dig_new_gamma_sum)]
            B_p = math.exp(digamma(new_alpha_p[1]) - digamma(new_alpha_p[0]) + log3)
            B_m = math.exp(digamma(new_alpha_m[1]) - digamma(new_alpha_m[0]) + log3)
            
            ### E step Calculation ###
            for n in range(N_p):
                colSum = 0.0
                for j in range(2):
                    EZ_p[n][j] = A_pi[j]
                    if RpList[n][j] > 0:
                        EZ_p[n][j] = EZ_p[n][j] * B_p
                    colSum += EZ_p[n][j]
                for j in range(2):
                    EZ_p[n][j] = EZ_p[n][j] / colSum

            for n in range(N_m):
                colSum = 0.0
                for j in range(2):
                    EZ_m[n][j] = A_pi[j]
                    if RmList[n][j] > 0:
                        EZ_m[n][j] = EZ_m[n][j] * B_m
                    colSum += EZ_m[n][j]
                for j in range(2):
                    EZ_m[n][j] = EZ_m[n][j] / colSum

            for n in range(N_o):
                colSum = 0.0
                for j in range(2):
                    EZ_o[n][j] = A_pi[j]
                    if RopList[n][j] > 0:
                        EZ_o[n][j] = EZ_o[n][j] * B_p
                    if RomList[n][j] > 0:
                        EZ_o[n][j] = EZ_o[n][j] * B_m
                    colSum += EZ_o[n][j]
                for j in range(2):
                    EZ_o[n][j] = EZ_o[n][j] / colSum

            ### cal updated lower bound ###
            Lq_old = Lq_new
            Lq_new = 0
            # Eq[ logPr(Rs+ | Z+,eps+) ]
            Lq_new += N_p * (dig_new_alpha_p_list[0] - dig_new_alpha_p_sum - log3)
            for n in range(N_p):
                for j in range(2):
                    Lq_new += EZ_p[n][j] * RpList[n][j] * \
                        (dig_new_alpha_p_list[1] - dig_new_alpha_p_list[0] + log3)

            # Eq[ logPr(Ro | Zo, eps+,eps-) ]
            Lq_new += N_o * (dig_new_alpha_p_list[0] - dig_new_alpha_p_sum - log3)
            Lq_new += N_o * (dig_new_alpha_m_list[0] - dig_new_alpha_m_sum - log3)

            for n in range(N_o):
                for j in range(2):
                    Lq_new += EZ_o[n][j] * RopList[n][j] * \
                        (dig_new_alpha_p_list[1] - dig_new_alpha_p_list[0] + log3)
                    Lq_new += EZ_o[n][j] * RomList[n][j] * \
                        (dig_new_alpha_m_list[1] - dig_new_alpha_m_list[0] + log3)

            # Eq[ logPr(Rs- | Z-, eps-) ]
            Lq_new += N_m * (dig_new_alpha_m_list[0] - dig_new_alpha_m_sum - log3)
            for n in range(N_m):
                for j in range(2):
                    Lq_new += EZ_m[n][j] * RmList[n][j] * \
                        (dig_new_alpha_m_list[1] - dig_new_alpha_m_list[0] + log3)

            # Eq[ logPr(Z+|PI) ]
            for n in range(N_p):
                for j in range(2):
                    Lq_new += EZ_p[n][j] * (dig_new_gamma_list[j] - dig_new_gamma_sum)
            # Eq[ logPr(Z-|PI) ]
            for n in range(N_m):
                for j in range(2):
                    Lq_new += EZ_m[n][j] * (dig_new_gamma_list[j] - dig_new_gamma_sum)
            # Eq[ logPr(Zo|PI) ]
            for n in range(N_o):
                for j in range(2):
                    Lq_new += EZ_o[n][j] * (dig_new_gamma_list[j] - dig_new_gamma_sum)
            # Eq[ logPr(PI|GAMMA) ]

            Lq_new -= (gammaln(gamma_pi[0]) + gammaln(gamma_pi[1]) - gammaln(gamma_pi[0] + gamma_pi[1]))
            for j in range(2):
                Lq_new += (gamma_pi[j] - 1) * (dig_new_gamma_list[j] - dig_new_gamma_sum)

            # Eq[ logPr(eps*|alpha+) ]
            Lq_new -= (gammaln(del_p1) + gammaln(del_p2) - gammaln(del_p1 + del_p2))
            Lq_new += (del_p1-1) * (dig_new_alpha_p_list[0] - dig_new_alpha_p_sum)
            Lq_new += (del_p2-1) * (dig_new_alpha_p_list[1] - dig_new_alpha_p_sum)

            # Eq[ logPr(eps-|alpha-) ]
            Lq_new -= (gammaln(del_m1) + gammaln(del_m2) - gammaln(del_m1 + del_m2))
            Lq_new += (del_m1-1) * (dig_new_alpha_m_list[0] - dig_new_alpha_m_sum)
            Lq_new += (del_m2-1) * (dig_new_alpha_m_list[1] - dig_new_alpha_m_sum)

            # Eq[ logq(Z+) ]
            for n in range(N_p):
                for j in range(2):
                    Lq_new -= EZ_p[n][j] * math.log(EZ_p[n][j])
            # Eq[ logq(Zo) ]
            for n in range(N_o):
                for j in range(2):
                    Lq_new -= EZ_o[n][j] * math.log(EZ_o[n][j])
            # Eq[ logq(Z-) ]
            for n in range(N_m):
                for j in range(2):
                    Lq_new -= EZ_m[n][j] * math.log(EZ_m[n][j])
            # Eq[ logq(PI) ]
            Lq_new += gammaln(new_gamma[0]) + gammaln(new_gamma[1])
            Lq_new -= gammaln(new_gamma[0] + new_gamma[1])
            for j in range(2):
                Lq_new -= (new_gamma[j] - 1) * (dig_new_gamma_list[j] - dig_new_gamma_sum)

            # Eq[ logq(eps+) ]
            if N_p + N_o > 0:
                Lq_new += gammaln(new_alpha_p[0]) + gammaln(new_alpha_p[1])
                Lq_new -= gammaln(new_alpha_p[0] + new_alpha_p[1])
                Lq_new -= (new_alpha_p[0] - 1) * (dig_new_alpha_p_list[0] - dig_new_alpha_p_sum)
                Lq_new -= (new_alpha_p[1] - 1) * (dig_new_alpha_p_list[1] - dig_new_alpha_p_sum)

            # Eq[logq(eps-)]
            if N_m + N_o > 0:
                Lq_new += gammaln(new_alpha_m[0]) + gammaln(new_alpha_m[1])
                Lq_new -= gammaln(new_alpha_m[0] + new_alpha_m[1])
                Lq_new -= (new_alpha_m[0] - 1) * (dig_new_alpha_m_list[0] - dig_new_alpha_m_sum)
                Lq_new -= (new_alpha_m[1] - 1) * (dig_new_alpha_m_list[1] - dig_new_alpha_m_sum)

            updateCount += 1

        return {'Lq':Lq_old,'eps+':new_alpha_p,'eps-':new_alpha_m}
