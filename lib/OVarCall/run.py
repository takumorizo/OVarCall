#!/usr/bin/env python
import subprocess
from OVarCall.methods.ovarCall import OVarCall
from OVarCall.utilOVar.exceptionUtil import *
# from OVarCall.filter.pileUp import *
from OVarCall.filter.pileUpFilter import PileUpFilter
from OVarCall.utilOVar.bamDAO import BamDAO
import logging
import os
import re
import ConfigParser


def runOVar(args):

    level = logging.getLevelName(args.log_level)
    logging.basicConfig(
        level=level, format='[%(asctime)s] %(levelname)-8s %(message)s',
        datefmt="%m-%d %H:%M:%S")
    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(args.parameter_settings)


    filterParamsT = {}
    filterParamsT.update(dict(config.items('tumor')))
    # filterParamsT['minBQ'] = args.base_quality

    filterParamsN = {}
    filterParamsN.update(dict(config.items('normal')))
    # filterParamsN['minBQ'] = args.base_quality

    filterTumor = PileUpFilter(filterParamsT)
    filterNormal = PileUpFilter(filterParamsN)

    ovarCall = None
    settingOVar = {}
    settingOVar['minMapQ'] = 30
    settingOVar.update(dict(config.items('calculation')))

    # settingOVar['minBQ'] = args.base_quality
    # settingOVar['minMapQ'] = args.mapping_quality
    ovarCall = OVarCall(settingOVar)
    tumorBamDao = BamDAO(args.bam1, settingOVar)
    normalBamDao = BamDAO(args.bam2, settingOVar)
    
    with tumorBamDao.openBam(), normalBamDao.openBam(), open(args.output, 'wb') as outputFile, open(os.devnull, 'wb') as FNULL:
        # logging.debug(str(ovarCall.__dict__))

        cmd_list = [args.samtools_path, 'mpileup', '-q',
                    str(settingOVar['minMapQ']), '-BQ', '0', '-d', '10000000', '-f', args.ref_fa, args.bam1, args.bam2]
        if args.region:
            cmd_list.insert(2, '-r')
            cmd_list.insert(3, args.region)
        pileup = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=FNULL)
        end_of_pipe = pileup.stdout
        lineCount = 0

        outputFile.writelines("TYPE\tchr\tpos\tref\tobs\tscoreWithOverlap\tscoreWithoutOverlap\tTumorPileup≠\tNormalPileup\n")
        
        for line in end_of_pipe:
            try:
                lineCount += 1
                if lineCount % 50000000 == 0:
                    logging.info('processed ' + str(lineCount) + ' line ')

                line = line.replace('\n', '')
                line = line.replace('\r', '')
                lineCols = re.split('\t', line)

                Chr = lineCols[0]
                pos = lineCols[1]
                ref = lineCols[2].upper()
                depth = lineCols[3]
                bases = lineCols[4]
                qualities = lineCols[5]
                controlDepth = lineCols[6]
                controlBases = lineCols[7]
                controlQualities = lineCols[8]

                candidate = filterTumor.satisfiedAll(Chr, pos, ref, depth, bases, qualities)
                if candidate is None:
                    continue
                candidate = filterNormal.filterSatisfy(
                    candidate, Chr, pos, ref, controlDepth, controlBases, controlQualities)

                for TYPE in 'MID':
                    for obs in candidate[TYPE]:
                        headT, bodyT = tumorBamDao.getOverlapInformation(
                            TYPE, Chr, pos, ref, obs)
                        headN, bodyN = normalBamDao.getOverlapInformation(
                            TYPE, Chr, pos, ref, obs)

                        if headT is None and bodyT is None:
                            # Cannot get read in Tumor sample
                            continue

                        if headN is None and bodyN is None:
                            # Cannot get read in Normal sample
                            continue

                        if filterTumor.filterOverlapPiled(headT, bodyT) and filterNormal.filterOverlapPiled(headN, bodyN):
                            logging.info("=======================================")
                            logging.info("position is evaluated by OVarCall Model")
                            logging.info("position : " + str([TYPE, Chr, pos, ref, obs]))
                            logging.info("tumor overlap piled : " + str(bodyT))
                            logging.info("normal overlap piled : " + str(bodyN))
                            
                            logging.info("start calc using normal and tumor ")
                            ansDic = ovarCall.doInference(
                                TYPE, Chr, pos, ref, obs, headT, bodyT, bodyN)

                            bodyT = ','.join(map(str, bodyT))
                            bodyN = ','.join(map(str, bodyN))

                            outputList = [TYPE, Chr, pos, ref, obs,ansDic["withOver"],ansDic["withoutOver"],bodyT,bodyN]
                            outputString = '\t'.join(map(str, outputList))
                            outputString = outputString.replace('\n', '')
                            outputString = outputString + '\n'
                            outputFile.writelines(outputString)
                            logging.info("=======================================")
                        else:
                            logging.info("=======================================")
                            logging.info("position is filtered after overlap filtering")
                            logging.info("position : " + str([TYPE, Chr, pos, ref, obs]))
                            logging.info("tumor piled : " + str(bodyT))
                            logging.info("normal piled : " + str(bodyN))
                            logging.info("=======================================")
                            continue
            except CliticalTypeException, e:
                logging.error('clitical exception : ' + str(line))
                raise e
            except WarningTypeException, e:
                logging.warning(
                    'Because of warning type exception of ' + str(type(e)) + ' passed this line : ' + str(line))
                continue

