#!/usr/bin/env python
from OVarCall.methods.ovarCall import OVarCall
from OVarCall.utilOVar.exceptionUtil import *
from OVarCall.filter.pileUpFilter import PileUpFilter
from OVarCall.utilOVar.bamDAO import BamDAO
import logging
import re
import ConfigParser
import vcf
import pysam
import traceback

def runOVarFilter(args):
    if args.annovarFile:
        runOVarFilterANNOVAR(args)
    else:
        runOVarFilterVCF(args)

def runOVarFilterANNOVAR(args):
    level = logging.getLevelName(args.log_level)
    logging.basicConfig(
        level=level, format='[%(asctime)s] %(levelname)-8s %(message)s',
        datefmt="%m-%d %H:%M:%S")
    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(args.parameter_settings)

    settingOVar = {}
    settingOVar['minMapQ'] = 30
    settingOVar.update(dict(config.items('calculation')))

    ovarCall = OVarCall(settingOVar)
    tumorBamDao = BamDAO(args.bam1, settingOVar)
    normalBamDao = BamDAO(args.bam2, settingOVar)

    filterParamsT = {}
    filterParamsT.update(dict(config.items('tumor')))

    filterParamsN = {}
    filterParamsN.update(dict(config.items('normal')))

    filterTumor = PileUpFilter(filterParamsT)
    filterNormal = PileUpFilter(filterParamsN)
    fasta = pysam.FastaFile(args.ref_fa)
    with tumorBamDao.openBam(), normalBamDao.openBam(), open(args.input,'r') as inputFile, open(args.output,'w') as outputFile:
        for line in inputFile:
            try:
                line = line.replace('\n', '')
                line = line.replace('\r', '')
                lineCols = re.split('\t', line)
                Chr, start, end, ref, alt = lineCols[0], lineCols[1], lineCols[2], lineCols[3], lineCols[4]
                
                try:
                    start = int(start)
                    end = int(end)
                except ValueError:
                    lineCols.append('OV+')
                    lineCols.append('OV-')
                    outputString = '\t'.join(map(str, lineCols))
                    outputString = outputString.replace('\n', '')
                    outputString = outputString + '\n'
                    outputFile.writelines(outputString)
                    continue

                TYPE = None
                pos = int(start)
                if ref != '-' and alt != '-':
                    TYPE = 'M'
                    ref = ref.upper()
                    alt = alt.upper()
                elif alt == '-':                
                    TYPE = 'D'
                    pos = start -1
                    alt = ref.upper()
                    ref = (fasta.fetch(Chr,pos-1,pos)).upper()
                elif ref == '-':
                    TYPE = 'I'
                    alt = alt.upper()
                    ref = (fasta.fetch(Chr,pos-1,pos)).upper()

                headT, bodyT = tumorBamDao.getOverlapInformation(
                    TYPE, Chr, pos, ref, alt)
                headN, bodyN = normalBamDao.getOverlapInformation(
                    TYPE, Chr, pos, ref, alt)

                if (headT is None and bodyT is None) or (headN is None and bodyN is None):
                    lineCols.append(-100) # overlap + 
                    lineCols.append(-100) # overlap -
                    logging.info("cannot get pileup info @ " + Chr + str(start) + "-" + str(end))
                else :
                    doFilter = args.filterPileup
                    passFilter = not(doFilter)
                    if doFilter:
                        if filterTumor.filterOverlapPiled(headT, bodyT) and filterNormal.filterOverlapPiled(headN, bodyN):
                            passFilter = True
                    
                    if passFilter:
                        ansDic = ovarCall.doInference(TYPE, Chr, pos, ref, alt, headT, bodyT, bodyN)
                        lineCols.append(ansDic["withOver"]) # overlap + 
                        lineCols.append(ansDic["withoutOver"]) # overlap -
                    else:
                        lineCols.append(-100) # overlap + 
                        lineCols.append(-100) # overlap -

                outputString = '\t'.join(map(str, lineCols))
                outputString = outputString.replace('\n', '')
                outputString = outputString + '\n'
                outputFile.writelines(outputString)

            except CliticalTypeException, e:
                logging.error('clitical exception : ' + str(line))
                raise e
            except WarningTypeException, e:
                logging.warning(
                    'Because of warning type exception of ' + str(type(e)) + ' passed this line : ' + str(line))
                continue
            except Exception, e:
                logging.error('type:' + str(type(e)))
                logging.error('args:' + str(e.args))
                logging.error('message:' + e.message)
                logging.error('e:' + str(e))
                logging.error('unexpected exception occurred. traceback here')
                logging.error(traceback.format_exc())
                continue
    fasta.close()



def runOVarFilterVCF(args):
    level = logging.getLevelName(args.log_level)
    logging.basicConfig(
        level=level, format='[%(asctime)s] %(levelname)-8s %(message)s',
        datefmt="%m-%d %H:%M:%S")
    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(args.parameter_settings)

    settingOVar = {}
    settingOVar['minMapQ'] = 30
    settingOVar.update(dict(config.items('calculation')))


    ovarCall = OVarCall(settingOVar)
    tumorBamDao = BamDAO(args.bam1, settingOVar)
    normalBamDao = BamDAO(args.bam2, settingOVar)
    
    filterParamsT = {}
    filterParamsT.update(dict(config.items('tumor')))

    filterParamsN = {}
    filterParamsN.update(dict(config.items('normal')))

    filterTumor = PileUpFilter(filterParamsT)
    filterNormal = PileUpFilter(filterParamsN)

    with tumorBamDao.openBam(), normalBamDao.openBam():
        vcf_reader = vcf.Reader(filename = args.input)
        vcf_reader.infos['OV+'] = vcf.parser._Info('OV+', -100, 'Float', "OVarCall Score with overlap", "OVarCall", "ver0.1.0")
        vcf_reader.infos['OV-'] = vcf.parser._Info('OV-', -100, 'Float', "OVarCall Score without overlap", "OVarCall", "ver0.1.0")
        
        vcf_writer = vcf.Writer(open(args.output + ".vcf", 'w'), vcf_reader)
        for vcf_record in vcf_reader:
            try:
                TYPE = None
                Chr = str(vcf_record.CHROM)
                pos = str(vcf_record.POS)
                ref = str(vcf_record.REF).upper()
                alt = str(vcf_record.ALT[0]).upper()
                if len(ref) == 1 and len(alt) == 1:
                    TYPE = 'M'
                elif len(ref) > len(alt):
                    TYPE = 'D'
                    alt = ref[1:]
                    ref = ref[0]
                elif len(ref) < len(alt):
                    TYPE = 'I'
                    ref = ref[0]
                    alt = alt[1:]

                headT, bodyT = tumorBamDao.getOverlapInformation(TYPE, Chr, pos, ref, alt)
                headN, bodyN = normalBamDao.getOverlapInformation(TYPE, Chr, pos, ref, alt)

                if (headT is None and bodyT is None) or (headN is None and bodyN is None):
                    vcf_record.INFO['OV+'] = -100
                    vcf_record.INFO['OV-'] = -100
                    logging.info("cannot get pileup info @ " + Chr + str(start) + "-" + str(end))
                else:
                    doFilter = args.filterPileup
                    passFilter = not(doFilter)
                    if doFilter:
                        if filterTumor.filterOverlapPiled(headT, bodyT) and filterNormal.filterOverlapPiled(headN, bodyN):
                            passFilter = True
                    
                    if passFilter:
                        ansDic = ovarCall.doInference(TYPE, Chr, pos, ref, alt, headT, bodyT, bodyN)
                        vcf_record.INFO['OV+'] = ansDic["withOver"]
                        vcf_record.INFO['OV-'] = ansDic["withoutOver"]
                    else:
                        vcf_record.INFO['OV+'] = -100
                        vcf_record.INFO['OV-'] = -100

                vcf_writer.write_record(vcf_record)

            except CliticalTypeException, e:
                logging.error('clitical exception : ' + str(line))
                raise e
            except WarningTypeException, e:
                logging.warning(
                    'Because of warning type exception of ' + str(type(e)) + ' passed this line : ' + str(line))
                continue
            except Exception, e:
                logging.error('type:' + str(type(e)))
                logging.error('args:' + str(e.args))
                logging.error('message:' + e.message)
                logging.error('e:' + str(e))
                logging.error('unexpected exception occurred. traceback here')
                logging.error(traceback.format_exc())
                continue

        vcf_writer.close()

