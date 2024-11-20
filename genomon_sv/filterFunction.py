#!/usr/bin/env python

"""
     functions for filtering candidates of structural variations
"""

from __future__ import print_function
import sys, gzip, subprocess, pysam, numpy, math, os, re
from . import coveredRegions
from . import realignmentFunction
from . import annotationFunction
from . import otherRealignmentFunctions
from . import utils
from scipy import stats


def genomon_sv_filt_main(output_prefix, args, thread_str = ""):

    utils.processingMessage("Filtering by # of breakpoint containing read pairs and variant sizes" + thread_str)
    filterJuncNumAndSize(output_prefix + ".junction.clustered.bedpe.gz",
                         output_prefix + ".junction.clustered.filt1.bedpe",
                         args.min_junc_num, args.min_sv_size, args.min_inversion_size)

    utils.processingMessage("Filtering by nonmatched control panel" + thread_str)
    filterNonMatchControl(output_prefix + ".junction.clustered.filt1.bedpe",
                          output_prefix + ".junction.clustered.filt2.bedpe",
                          args.non_matched_control_junction,
                          args.matched_control_label,
                          args.control_panel_num_thres, args.control_panel_check_margin)

    utils.processingMessage("Incorporating improperly alinged read pair infomation" + thread_str)
    addImproperInfo(output_prefix + ".junction.clustered.filt2.bedpe",
                    output_prefix + ".junction.clustered.filt3.bedpe",
                    args.output_prefix + ".improper.clustered.bedpe.gz")

    utils.processingMessage("Filtering by sizes of covered regions, mapping quality and # of support read pairs" + thread_str)
    filterMergedJunc(output_prefix + ".junction.clustered.filt3.bedpe",
                     output_prefix + ".junction.clustered.filt4.bedpe",
                     args.min_support_num, args.min_mapping_qual, args.min_overhang_size)

    utils.processingMessage("Filtering too close candidates" + thread_str)
    removeClose(output_prefix + ".junction.clustered.filt4.bedpe",
                output_prefix + ".junction.clustered.filt5.bedpe",
                args.close_check_margin, args.close_check_thres)

    realignment_tool = 'parasail'
    if args.blat == True:
        realignment_tool = 'blat'
    elif args.edlib == True:
        realignment_tool = 'edlib'
    utils.processingMessage("Performing realignments using " + realignment_tool + thread_str)
    validateByRealignment(output_prefix + ".junction.clustered.filt5.bedpe",
                          output_prefix + ".junction.clustered.filt6.bedpe",
                          args.bam_file, args.matched_control_bam, args.reference_genome, args.blat_option,
                          args.short_tandem_reapeat_thres, args.max_read_num, args.search_length, args.search_margin, 
                          args.split_refernece_thres, args.validate_sequence_length, args.blat, args.edlib)

    utils.processingMessage("Filtering allele frequencies, Fisher's exact test p-values and # of support read pairs" + thread_str)
    filterNumAFFis(output_prefix + ".junction.clustered.filt6.bedpe", 
                   output_prefix + ".junction.clustered.filt7.bedpe",
                   args.matched_control_bam,
                   args.min_tumor_variant_read_pair, args.min_tumor_allele_freq, 
                   args.max_control_variant_read_pair, args.max_control_allele_freq,
                   args.max_fisher_pvalue)

    utils.processingMessage("Adding annotation" + thread_str)
    annotationFunction.addAnnotation(output_prefix + ".junction.clustered.filt7.bedpe",
                                     output_prefix + ".genomonSV.result.txt",
                                     args.genome_id, args.grc)

    if args.debug == False:
        subprocess.call(["rm", output_prefix + ".junction.clustered.filt1.bedpe"])
        subprocess.call(["rm", output_prefix + ".junction.clustered.filt2.bedpe"])
        subprocess.call(["rm", output_prefix + ".junction.clustered.filt3.bedpe"])
        subprocess.call(["rm", output_prefix + ".junction.clustered.filt4.bedpe"])
        subprocess.call(["rm", output_prefix + ".junction.clustered.filt5.bedpe"])
        subprocess.call(["rm", output_prefix + ".junction.clustered.filt6.bedpe"])
        subprocess.call(["rm", output_prefix + ".junction.clustered.filt7.bedpe"])


def partition_junction(output_prefix, thread_num):

    # count the number of junctions
    line_num = 0
    with gzip.open(output_prefix + ".junction.clustered.bedpe.gz") as hin:
        for line in hin:
            line_num = line_num + 1

    thread_num_mod = min(line_num, thread_num)
    if thread_num_mod == 0: thread_num_mod = 1

    each_partition_line_num = line_num / thread_num_mod

    current_line_num = 0
    current_partition_num = 1
    hout = open(output_prefix + ".thread_1.junction.clustered.bedpe", 'w')
    with gzip.open(output_prefix + ".junction.clustered.bedpe.gz", 'rt') as hin:
        for line in hin:
            print(line.rstrip('\n'), file = hout)
            current_line_num = current_line_num + 1
            if current_line_num >= each_partition_line_num and current_partition_num < thread_num_mod:
                current_line_num = 0
                current_partition_num = current_partition_num + 1
                hout.close()
                hout = open(output_prefix + ".thread_" + str(current_partition_num) + ".junction.clustered.bedpe", 'w')   
    
    hout.close()

    for i in range(1, thread_num_mod + 1):
        utils.compress_index_bed(output_prefix + ".thread_" + str(i) + ".junction.clustered.bedpe",
                                 output_prefix + ".thread_" + str(i) + ".junction.clustered.bedpe.gz")
        subprocess.check_call(["rm", "-rf", output_prefix + ".thread_" + str(i) + ".junction.clustered.bedpe"])
 
    return thread_num_mod
    


def filterJuncNumAndSize(inputFilePath, outputFilePath, junc_num_thres, sv_size_thres, inversion_size_thres):
     
    """
        script for filtering by the length of SV size and support read junctions
    """

    hIN = gzip.open(inputFilePath, 'rt')
    hOUT = open(outputFilePath, 'w')

    for line in hIN:
        F = line.rstrip('\n').split('\t')

        ##########
        # for now only consider long ranged SV??
        svLen = abs(int(F[2]) - int(F[5]))
        if F[7] != "---": svLen = svLen + len(F[7])
        if F[0] == F[3]:
            # for inversion, use other threshould 
            if F[8] == F[9] and svLen < int(inversion_size_thres):
                continue
            if svLen < int(sv_size_thres): continue
        ########## 
    
        ##########
        # enumerate support read number
        IDs = F[17].split(';')
        juncSupport = int(len(IDs))

        # skip if the number of suppor read is below the minSupReadNum 
        if juncSupport < int(junc_num_thres): continue
        ##########

        print('\t'.join(F), file = hOUT)


    hIN.close()
    hOUT.close()




def filterNonMatchControl(inputFilePath, outputFilePath, controlFile, matchedNormal, controlPanel_num_thres, controlPanel_check_margin):

    """
        script for removing candidate in which 
        non-matched normals have the junction reads
    """

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    use_control = True if controlFile != "" else False
    if use_control == True: tabixfile = pysam.TabixFile(controlFile)

    tabixErrorMsg = "" 
    for line in hIN:
        F = line.rstrip('\n').split('\t')


        controlFlag = 0
        max_control_sample = "---"
        max_control_num = 0

        if use_control == True:

            inseqSize = (0 if F[7] == "---" else len(F[7]))

            ####################
            # get the records for control junction data for the current position
            tabixErrorFlag = 0
            try:
                records = tabixfile.fetch(F[0], max(int(F[1]) - controlPanel_check_margin, 0), int(F[2]) + controlPanel_check_margin)
            except Exception as inst:
                # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
                tabixErrorMsg = str(inst.args)
                tabixErrorFlag = 1
            ####################


            ####################
            # for each record in control junction extracted, check the consistency with the current junction
            # max_control_sample = "---"
            # max_control_num = 0 
            if tabixErrorFlag == 0:
                for record_line in records:
                    record = record_line.split('\t')

                    if F[0] == record[0] and F[3] == record[3] and F[8] == record[8] and F[9] == record[9]:

                        flag = 0
                        # detailed check on the junction position considering inserted sequences
                        if F[8] == "+":
                            expectedDiffSize = (int(F[2]) - int(record[2])) + (inseqSize - int(record[7]))
                            if (F[9] == "+" and int(F[5]) == int(record[5]) - int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(record[5]) + int(expectedDiffSize)):
                                flag = 1
                        else:
                            expectedDiffSize = (int(F[2]) - int(record[2])) + (int(record[7]) - inseqSize)
                            if (F[9] == "+" and int(F[5]) == int(record[5]) + int(expectedDiffSize)) or (F[9] == "-" and int(F[5]) == int(record[5]) - int(expectedDiffSize)):
                                flag = 1

                        # if position relationship including inserted sequences matches
                        if flag == 1:
                            controlSamples = record[10].split(';')
                            controlNums = record[11].split(';')

                            for i in range(0, len(controlSamples)):
                                if controlSamples[i] == matchedNormal: continue

                                if int(controlNums[i]) > max_control_num:
                                    max_control_sample = controlSamples[i]
                                    max_control_num = int(controlNums[i])

                                if int(controlNums[i]) >= int(controlPanel_num_thres):
                                    controlFlag = 1

                                """
                                # if controlSamples[i] != matchedNormal is not None and int(controlNums[i]) >= int(controlPanel_num_thres):
                                # if controlSamples[i] != matchedNormal and int(controlNums[i]) >= int(supportReadThres):
                                    controlFlag = 1
                                    if int(controlNums[i]) > max_control_num:
                                        max_control_sample = controlSamples[i]
                                        max_control_num = int(controlNums[i]) 
                                """
                            
            ####################
              
        if controlFlag == 0:
            print("\t".join(F) + '\t' + max_control_sample + '\t' + str(max_control_num), file = hOUT)


    if tabixErrorMsg != "":
        utils.warningMessage("One or more error occured in tabix file fetch, e.g.: " + tabixErrorMsg)

    hIN.close()
    hOUT.close()
    if use_control == True: tabixfile.close()



def addImproperInfo(inputFilePath, outputFilePath, improperFilePath):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')
    tabixfile = pysam.TabixFile(improperFilePath)

    tabixErrorMsg = ""
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        improper_readNames = "---"
        improper_MQs = "---"
        improper_coveredRegions = "---"

        # get the records for control junction data for the current position
        tabixErrorFlag = 0
        try:
            records = tabixfile.fetch(F[0], int(F[1]) - 1, int(F[2]) + 1)
        except Exception as inst:
            # print >> sys.stderr, "%s: %s" % (type(inst), inst.args)
            tabixErrorMsg = str(inst.args)
            tabixErrorFlag = 1


        if tabixErrorFlag == 0:
            for record_line in records:
                FF = record_line.split('\t')

                if F[0] == FF[0] and F[3] == FF[3] and int(F[1]) >= int(FF[1]) and int(F[2]) <= int(FF[2]) and int(F[4]) >= int(FF[4]) and int(F[5]) <= int(FF[5]) and F[8] == FF[8] and F[9] == FF[9]:
                    improper_readNames = FF[11]
                    improper_MQs = FF[7]
                    improper_coveredRegions = FF[10]
    
        # last two columns are about pooled control information
        print("\t".join(F[:-2]) + "\t" + improper_readNames + "\t" + improper_MQs + "\t" + improper_coveredRegions + '\t' + '\t'.join(F[-2:]), file = hOUT)

    if tabixErrorMsg != "":
        utils.warningMessage("One or more error occured in tabix file fetch, e.g.: " + tabixErrorMsg)


    hIN.close()
    hOUT.close()
    tabixfile.close()



def filterMergedJunc(inputFilePath, outputFilePath, min_support_num, min_mapping_qual, min_cover_size):

    """    
    script for filtering the inferred break points for structural variation

    as conditions for filtering
    min_support_num: the number of support reads
    min_mapping_qual: (e.g., in either type 1 or 2 junction reads, median of mapping quality is at least 40 ??)
    min_cover_size: the size of alignment covered region (e.g., around the both break points, at least 80 bp sequences have to be covered)
    """

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')


    for line in hIN:
        F = line.rstrip('\n').split('\t')

        ##########
        MQs = F[10].split(';')
        coveredRegion_primary = F[11].split(';')
        coveredRegion_pair = F[13].split(';')
        pairMQs = map(lambda x: int(x), F[12].split(';'))
        juncPairPos = F[14].split(';')
        juncReadTypes = F[15].split(';')
        improperCoveredRegion = ([] if F[20] == "---" else F[20].split(';'))

        # enumerate support read number
        junc_ids = [re.sub(r'/\d$', '', x) for x in F[17].split(';')]
        # improper_ids = F[17].split(';')
        improper_ids = [] if F[18] == "---" else F[18].split(';')
        if len(list(set(junc_ids + improper_ids))) < min_support_num:
            continue

        juncSupport = len(MQs)
        improperMQs = ([] if F[19] == "---" else F[19].split(';'))
        """
        improperSupport = len(improperMQs)
        # skip if the number of suppor read is below the minSupReadNum 
        if juncSupport + improperSupport < min_support_num:
            continue
        """
        

        ####################
        # check for the mapping quality
        MQs1 = [];
        MQs2 = [];
        for i in range(0, len(MQs)):
            if juncReadTypes[i] == "1":
                MQs1.append(int(MQs[i]))
            else:
                MQs2.append(int(MQs[i]))

        # improper pair
        IMQs1 = [];
        IMQs2 = [];
        for i in range(0, len(improperMQs)):
            tempMQ = improperMQs[i].split(',')
            IMQs1.append(int(tempMQ[0]))
            IMQs2.append(int(tempMQ[1]))

     
        mapFlag = 0
        if len(MQs1) > 0 and numpy.median(MQs1) >= min_mapping_qual:
            mapFlag = 1

        if len(MQs2) > 0 and numpy.median(MQs2) >= min_mapping_qual:
            mapFlag = 1     

        if len(IMQs1) > 0 and numpy.median(IMQs1) >= min_mapping_qual:
            mapFlag = 1

        if len(IMQs2) > 0 and numpy.median(IMQs2) >= min_mapping_qual:
            mapFlag = 1

        if mapFlag == 0:
            continue


        #################### 
        # check for the covered region
        region1 = coveredRegions.coveredRegions()
        region2 = coveredRegions.coveredRegions()

        for i in range(0, len(coveredRegion_primary)):
            regPair = coveredRegion_primary[i].split(',')
            if juncReadTypes[i] == "1":
                region1.addMerge(regPair[0])
                region2.addMerge(regPair[1])
            else:
                region1.addMerge(regPair[1])
                region2.addMerge(regPair[0])


        for i in range(0, len(coveredRegion_pair)):
            if (juncReadTypes[i] == "1" and juncPairPos[i] == "1") or (juncReadTypes[i] == "2" and juncPairPos[i] == "2"): 
                region1.addMerge(coveredRegion_pair[i])
            elif (juncReadTypes[i] == "1" and juncPairPos[i] == "2") or (juncReadTypes[i] == "2" and juncPairPos[i] == "1"):
                region2.addMerge(coveredRegion_pair[i])


        for i in range(0, len(improperCoveredRegion)):
            regPair = improperCoveredRegion[i].split(',')
            region1.addMerge(regPair[0])
            region2.addMerge(regPair[1])


        region1.reduceMerge()
        region2.reduceMerge()

        if region1.regionSize() < min_cover_size or region2.regionSize() < min_cover_size:
            continue

        # every condition is satisfied
        print('\t'.join(F) + '\t' + str(region1.regionSize()) + '\t' + str(region2.regionSize()), file = hOUT)


    hIN.close()
    hOUT.close()



def removeClose(inputFilePath, outputFilePath, close_check_margin, close_check_thres):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')

    key2info = {}
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        delList = []
        skipFlag = 0
        for tkey in key2info:

            tchr1, tstart1, tend1, tchr2, tstart2, tend2, inseq, tdir1, tdir2 = tkey.split('\t')

            if F[0] != tchr1 or int(F[1]) > int(tend1) + close_check_margin:

                print(key2info[tkey], file = hOUT)
                delList.append(tkey)
        
            else:

                if F[0] == tchr1 and F[3] == tchr2 and F[8] == tdir1 and F[9] == tdir2 and abs(int(F[2]) - int(tend1)) <= close_check_margin and abs(int(F[5]) - int(tend2)) <= close_check_margin:

                    infos = key2info[tkey].split('\t')
                    if len(F[17].split(';')) < len(infos[17].split(';')) and len(F[17].split(';')) < int(close_check_thres):
                        skipFlag = 1
                    elif len(F[17].split(';')) > len(infos[17].split(';')) and len(infos[17].split(';')) < int(close_check_thres):
                        delList.append(tkey)


        for tkey in delList:
            del key2info[tkey]

        if skipFlag == 0:
            key2info['\t'.join(F[0:6] + F[7:10])] = '\t'.join(F)

    hIN.close()


    for tkey in key2info:
        print(key2info[tkey], file = hOUT)


    hOUT.close()



def validateByRealignment(inputFilePath, outputFilePath, tumorBamFilePath, normalBamFilePath, reference_genome, blat_option,
                          short_tandem_reapeat_thres, max_read_num, search_length, search_margin, split_refernece_thres, validate_sequence_length, f_blat, f_edlib):


    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')
    blat_cmds = ("blat " + blat_option).split(' ')

    num = 1
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        # if len(F) >= 24:
        chr1, pos1, dir1, chr2, pos2, dir2, juncSeq, max_control_sample, max_control_num, overhang1, overhang2 = F[0], F[2], F[8], F[3], F[5], F[9], F[7], F[21], F[22], F[23], F[24]
        # else:
        #     chr1, pos1, dir1, chr2, pos2, dir2, juncSeq = F[0], F[2], F[8], F[3], F[5], F[9], F[7]
        #     max_control_sample, max_control_num, overhang1, overhang2 = "---", "---", "---", "---"

        STDFlag = 0
        if chr1 == chr2 and int(pos2) - int(pos1) < short_tandem_reapeat_thres and dir1 == "-" and dir2 == "+": STDFlag = 1

        ####################
        # extract short reads from tumor sequence data around the candidate
        fRet = realignmentFunction.extractSVReadPairs(tumorBamFilePath, outputFilePath + ".tmp.tumor.fa", chr1, pos1, dir1, chr2, pos2, dir2, 
                                                      max_read_num, search_length, search_margin, reference_genome)
        if fRet == 1: continue

        if normalBamFilePath != "":
            # extract short reads from matched-control sequence data around the candidate
            fRet = realignmentFunction.extractSVReadPairs(normalBamFilePath, outputFilePath + ".tmp.normal.fa", chr1, pos1, dir1, chr2, pos2, dir2,
                                                          max_read_num, search_length, search_margin, reference_genome)
            if fRet == 1: continue
        ####################
        
        ####################
        # generate reference sequence and sequence containing the presumed variants
        if int(pos1) - validate_sequence_length < 0: continue
        if int(pos2) - validate_sequence_length < 0: continue
        realignmentFunction.getRefAltForSV(outputFilePath + ".tmp.refalt.fa", chr1, pos1, dir1, chr2, pos2, dir2, juncSeq,
                                           reference_genome, split_refernece_thres, validate_sequence_length)

        tumorRef, tumorAlt = "---", "---"
        normalRef, normalAlt = "---", "---"
        if f_blat == True:

            ####################
            # alignment tumor short reads to the reference and alternative sequences
            FNULL = open(os.devnull, 'w')
            subprocess.call(blat_cmds + [outputFilePath + ".tmp.refalt.fa", outputFilePath + ".tmp.tumor.fa", outputFilePath + ".tmp.tumor.psl"],
                        stdout = FNULL, stderr = subprocess.STDOUT)
        
            ####################
            # alignment normal short reads to the reference and alternative sequences
            if normalBamFilePath != "":
                subprocess.call(blat_cmds + [outputFilePath + ".tmp.refalt.fa", outputFilePath + ".tmp.normal.fa", outputFilePath + ".tmp.normal.psl"],
                            stdout = FNULL, stderr = subprocess.STDOUT)

            FNULL.close()

            ####################
            # summarize alignment results
            tumorRef, tumorAlt = realignmentFunction.summarizeRefAlt(outputFilePath + ".tmp.tumor.psl", STDFlag)

            if normalBamFilePath != "":
                normalRef, normalAlt = realignmentFunction.summarizeRefAlt(outputFilePath + ".tmp.normal.psl", STDFlag)

        elif f_edlib == True:

            fa_alt, fa_ref1, fa_ref2, fa_ref = otherRealignmentFunctions.getRefAltForSV(outputFilePath + ".tmp.refalt.fa")
            tumorRef, tumorAlt = otherRealignmentFunctions.summarizeRefAltEdlib(outputFilePath + ".tmp.tumor.fa", STDFlag, fa_alt, fa_ref1, fa_ref2, fa_ref)
            if normalBamFilePath != "":
                normalRef, normalAlt = otherRealignmentFunctions.summarizeRefAltEdlib(outputFilePath + ".tmp.normal.fa", STDFlag, fa_alt, fa_ref1, fa_ref2, fa_ref)

        else:

            fa_alt, fa_ref1, fa_ref2, fa_ref = otherRealignmentFunctions.getRefAltForSV(outputFilePath + ".tmp.refalt.fa")
            tumorRef, tumorAlt = otherRealignmentFunctions.summarizeRefAltParasail(outputFilePath + ".tmp.tumor.fa", STDFlag, fa_alt, fa_ref1, fa_ref2, fa_ref)
            if normalBamFilePath != "":
                normalRef, normalAlt = otherRealignmentFunctions.summarizeRefAltParasail(outputFilePath + ".tmp.normal.fa", STDFlag, fa_alt, fa_ref1, fa_ref2, fa_ref)



        # fisher test
        lpvalue = "---"
        if normalBamFilePath != "":
            oddsratio, pvalue = stats.fisher_exact([[tumorRef, tumorAlt], [normalRef, normalAlt]], 'less')
            if pvalue < 1e-100: pvalue = 1e-100
            lpvalue = (- math.log(pvalue, 10) if pvalue < 1 else 0)
            lpvalue = str(round(lpvalue, 4)) 

        print('\t'.join([chr1, pos1, dir1, chr2, pos2, dir2, juncSeq, \
                         str(tumorRef), str(tumorAlt), str(normalRef), str(normalAlt), str(lpvalue), \
                         max_control_sample, max_control_num, overhang1, overhang2]), file = hOUT)

        if num % 100 == 0:        
            print("finished checking " + str(num) + " candidates", file = sys.stderr)
        num = num + 1

    
    if num > 1:
        subprocess.call(["rm", outputFilePath + ".tmp.tumor.fa"])
        subprocess.call(["rm", outputFilePath + ".tmp.refalt.fa"])
        if f_blat == True:
            subprocess.call(["rm", outputFilePath + ".tmp.tumor.psl"])


        if normalBamFilePath != "":
            subprocess.call(["rm", outputFilePath + ".tmp.normal.fa"])
            if f_blat == True:
                subprocess.call(["rm", outputFilePath + ".tmp.normal.psl"])

    hIN.close()
    hOUT.close()



def filterNumAFFis(inputFilePath, outputFilePath, normalBamFilePath, 
                   min_tumor_read_pair, min_tumor_alleleFreq, max_control_read_pair, max_control_allele_freq, max_fisher_pvalue):

    hIN = open(inputFilePath, 'r')
    hOUT = open(outputFilePath, 'w')
  
 
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        tumorAF = 0 
        if float(F[7]) + float(F[8]) > 0: tumorAF = float(F[8]) / (float(F[7]) + float(F[8]))     
        tumorAF = str(round(tumorAF, 4))

        normalAF = "---"
        if normalBamFilePath != "":
            normalAF = 0
            if float(F[9]) + float(F[10]) > 0: normalAF = float(F[10]) / (float(F[9]) + float(F[10]))
            normalAF = str(round(normalAF, 4))

        if int(F[8]) < int(min_tumor_read_pair): continue
        if float(tumorAF) < float(min_tumor_alleleFreq): continue

        if normalBamFilePath != "":
            if int(F[10]) > int(max_control_read_pair): continue
            if float(normalAF) > float(max_control_allele_freq): continue
            if 10**(-float(F[11])) > float(max_fisher_pvalue): continue

        print('\t'.join(F[0:9]) + '\t' + tumorAF + '\t' + '\t'.join(F[9:11]) + '\t' + normalAF + '\t' + '\t'.join(F[11:]), file = hOUT)

    hIN.close()
    hOUT.close()


 
