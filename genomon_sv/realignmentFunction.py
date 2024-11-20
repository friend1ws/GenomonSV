#!/usr/bin/env python

from __future__ import print_function
import sys, pysam, random
from . import utils

def extractSVReadPairs(bamFilePath, outputFilePath, juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2, max_read_num, search_length, search_margin, reference_genome):

    """
        read pairs containing break points are extracted. (yshira 2015/04/23)
        The exact condition is as follows:

        1. one of the read in the pair has the break point of the SV candidate
        2. the start positions of the read pairs are within 800bp of the break point of the SV candidate 

        Some minor concern for the above conditions are:
        1. Depending on the choice of the "start position" or "end position", the distance between the read and break point differs. This can generate slight bias...
        (but I believe we can recover this by setting sufficient margin (800bp), and summarize the alignment result carefully.)
        2. Maybe, for some of the read pair, the result of alignment is obvious. But should we re-align them?

    """

    bamfile = utils.getPysamSamfile(bamFilePath, reference_genome)

    # if the #sequence read is over the `maxDepth`, then that key is ignored
    # depthFlag = 0
    # if bamfile.count(juncChr1, max(0,int(juncPos1) - 1), int(juncPos1) + 1) >= max_depth: depthFlag = 1
    # if bamfile.count(juncChr2, max(0,int(juncPos2) - 1), int(juncPos2) + 1) >= max_depth: depthFlag = 1
    # if depthFlag == 1:
    #     print("sequence depth exceeds the threshould for: " + ','.join([juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2]), file = sys.stderr)
    #     return 1 

    hOUT = open(outputFilePath, 'w')

    readID2exist = {}    
    for read in bamfile.fetch(juncChr1, max(0, int(juncPos1) - search_length), int(juncPos1) + search_length):

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip unmapped read 
        if flags[2] == "1" or flags[3] == "1": continue 

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue

        chr_current = bamfile.getrname(read.tid)
        pos_current = int(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")
        chr_pair = bamfile.getrname(read.rnext)
        pos_pair = int(read.pnext + 1)
        dir_pair = ("-" if flags[5] == "1" else "+")

        # the read (with margin) contains break point
        if pos_current - search_margin <= int(juncPos1) <= (read.aend - 1) + search_margin:
            readID2exist[read.qname] = 1
    
        # the read pair covers break point
        if chr_pair == juncChr1 and pos_current <= int(juncPos1) <= pos_pair and dir_current == "+" and dir_pair == "-":
            readID2exist[read.qname] = 1

        # the read pair covers break point
        if chr_pair == juncChr2:
            juncFlag = 0
            if juncDir1 == "+" and juncDir2 == "+" and pos_current <= int(juncPos1) and pos_pair <= int(juncPos2): juncFlag = 1
            if juncDir1 == "+" and juncDir2 == "-" and pos_current <= int(juncPos1) and pos_pair >= int(juncPos2): juncFlag = 1
            if juncDir1 == "-" and juncDir2 == "+" and pos_current >= int(juncPos1) and pos_pair <= int(juncPos2): juncFlag = 1
            if juncDir1 == "-" and juncDir2 == "-" and pos_current >= int(juncPos1) and pos_pair >= int(juncPos2): juncFlag = 1

            if juncFlag == 1:  
                readID2exist[read.qname] = 1


    for read in bamfile.fetch(juncChr2, max(0, int(juncPos2) - search_length), int(juncPos2) + search_length):
        
        if read.qname == "ST-E00104:162:H03UUALXX:5:1222:21168:16006":
            pass
 
        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip unmapped read 
        if flags[2] == "1" or flags[3] == "1": continue
        
        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue
        
        # skip duplicated reads
        if flags[10] == "1": continue
        
        chr_current = bamfile.getrname(read.tid)
        pos_current = int(read.pos + 1)
        dir_current = ("-" if flags[4] == "1" else "+")
        chr_pair = bamfile.getrname(read.rnext)
        pos_pair = int(read.pnext + 1)
        dir_pair = ("-" if flags[5] == "1" else "+")

        # the read (with margin) contains break point
        if pos_current - search_margin <= int(juncPos2) <= (read.aend - 1) + search_margin:
            readID2exist[read.qname] = 1
                
        # the read pair covers break point
        if chr_pair == juncChr2 and pos_current <= int(juncPos2) <= pos_pair and dir_current == "+" and dir_pair == "-":
            readID2exist[read.qname] = 1
                
        # the read pair covers break point
        if chr_pair == juncChr1:
            juncFlag = 0
            if juncDir2 == "+" and juncDir1 == "+" and pos_current <= int(juncPos2) and pos_pair <= int(juncPos1): juncFlag = 1
            if juncDir2 == "+" and juncDir1 == "-" and pos_current <= int(juncPos2) and pos_pair >= int(juncPos1): juncFlag = 1
            if juncDir2 == "-" and juncDir1 == "+" and pos_current >= int(juncPos2) and pos_pair <= int(juncPos1): juncFlag = 1
            if juncDir2 == "-" and juncDir1 == "-" and pos_current >= int(juncPos2) and pos_pair >= int(juncPos1): juncFlag = 1
             
            if juncFlag == 1:
                readID2exist[read.qname] = 1


    readID2seq1 = {}
    readID2seq2 = {}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    for read in bamfile.fetch(juncChr1, max(0, int(juncPos1) - search_length), int(juncPos1) + search_length):

        if read.qname in readID2exist:
        
            # get the flag information
            flags = format(int(read.flag), "#014b")[:1:-1]

            # skip unmapped read 
            if flags[2] == "1" or flags[3] == "1": continue

            # skip supplementary alignment
            if flags[8] == "1" or flags[11] == "1": continue

            # skip duplicated reads
            if flags[10] == "1": continue

            tempSeq = ""
            if flags[4] == "1":
                tempSeq = utils.reverseComplement(str(read.seq))
            else:
                tempSeq = read.seq
 
            # the first read
            if flags[6] == "1":
                readID2seq1[read.qname] = tempSeq
            else:
                readID2seq2[read.qname] = tempSeq


    for read in bamfile.fetch(juncChr2, max(0, int(juncPos2) - search_length), int(juncPos2) + search_length):

        if read.qname in readID2exist:

            # get the flag information
            flags = format(int(read.flag), "#014b")[:1:-1]

            # skip unmapped read 
            if flags[2] == "1" or flags[3] == "1": continue

            # skip supplementary alignment
            if flags[8] == "1" or flags[11] == "1": continue
            
            # skip duplicated reads
            if flags[10] == "1": continue

            tempSeq = ""
            if flags[4] == "1":
                tempSeq = utils.reverseComplement(str(read.seq))
            else:
                tempSeq = read.seq

            # the first read
            if flags[6] == "1":
                readID2seq1[read.qname] = tempSeq
            else:
                readID2seq2[read.qname] = tempSeq


    for readID in fastq_write_subsample_fetch(readID2seq1, readID2seq2, max_read_num):
        print('>' + readID + '/1', file = hOUT) 
        print(readID2seq1[readID], file = hOUT)
        print('>' + readID + '/2', file = hOUT)
        print(readID2seq2[readID], file = hOUT)

    bamfile.close()
    hOUT.close()

    return 0


def fastq_write_subsample_fetch(readID2seq1, readID2seq2, max_read_num):

    rec = 0
    for readID in readID2seq1:
        if readID in readID2seq2:
            rec = rec + 1

    selected_inds = random.sample(range(rec), min(max_read_num, rec))

    rec2 = 0
    for readID in readID2seq1:
        if readID in readID2seq2:
            if rec2 in selected_inds:
                yield readID
            rec2 = rec2 + 1


def getRefAltForSV(outputFilePath, juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2, juncSeq, reference_genome, split_refernece_thres, validate_sequence_length):

    """
        for short SV (mid-range (<= split_refernece_thres bp) deletion, tandem duplication), we get the two sequence
        for large SV (> split_refernece_thres bp), we get three sequence (one joint sequence by two break points, and two reference sequences around the break points)

        the only concern is short inversion... (are there some somatic short inversion?)
        however, this will be filtered beforehand by the "cover filter", and maybe we have to give up detecting these class of SVs.

    """

    hOUT = open(outputFilePath, 'w')

    if juncSeq == "---": juncSeq = ""

    # for mid-range deletion or tandem duplication
    if juncChr1 == juncChr2 and abs(int(juncPos1) - int(juncPos2)) <= split_refernece_thres and juncDir1 != juncDir2:

        seq = utils.get_seq(reference_genome, juncChr1, int(juncPos1) - validate_sequence_length, int(juncPos2) + validate_sequence_length) 
        """
        seq = ""
        for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(int(juncPos1) - validate_sequence_length) + "-" + str(int(juncPos2) + validate_sequence_length)):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()
        seq = seq.replace('>', '')
        seq = seq.replace(juncChr1 + ":" + str(int(juncPos1) - validate_sequence_length) + "-" + str(int(juncPos2) + validate_sequence_length), '')
        """
        print('>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_ref", file = hOUT)
        print(seq, file = hOUT)

        # for mid-range deletion
        if juncDir1 == "+" and juncDir2 == "-":
   
            seq = utils.get_seq(reference_genome, juncChr1, int(juncPos1) - validate_sequence_length, int(juncPos1)) 
            seq = seq + juncSeq
            seq = seq + utils.get_seq(reference_genome, juncChr2, int(juncPos2), int(juncPos2) + validate_sequence_length)

            """
            seq = ""
            for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(int(juncPos1) - validate_sequence_length) + "-" + str(juncPos1)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()
            seq = seq.replace('>', '')
            seq = seq.replace(juncChr1 + ":" + str(int(juncPos1) - validate_sequence_length) + "-" + str(juncPos1), '')

            seq = seq + juncSeq

            for item in pysam.faidx(reference_genome, juncChr2 + ":" + str(juncPos2) + "-" + str(int(juncPos2) + validate_sequence_length)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()
            seq = seq.replace('>', '')
            seq = seq.replace(juncChr2 + ":" + str(juncPos2) + "-" + str(int(juncPos2) + validate_sequence_length), '')
            """

            print('>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_alt", file = hOUT)
            print(seq, file = hOUT)

        # for mid-range tandem duplication
        else:

            seq = utils.get_seq(reference_genome, juncChr2, int(juncPos2) - validate_sequence_length, int(juncPos2))                             
            seq = seq + juncSeq
            seq = seq + utils.get_seq(reference_genome, juncChr1, int(juncPos1), int(juncPos1) + validate_sequence_length)

            """
            seq = "" 
            for item in pysam.faidx(reference_genome, juncChr2 + ":" + str(int(juncPos2) - validate_sequence_length) + "-" + str(juncPos2)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()
            seq = seq.replace('>', '')
            seq = seq.replace(juncChr2 + ":" + str(int(juncPos2) - validate_sequence_length) + "-" + str(juncPos2), '')
            
            seq = seq + juncSeq

            for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(juncPos1) + "-" + str(int(juncPos1) + validate_sequence_length)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()
            seq = seq.replace('>', '')
            seq = seq.replace(juncChr1 + ":" + str(juncPos1) + "-" + str(int(juncPos1) + validate_sequence_length), '')
            """

            print('>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_alt", file = hOUT)
            print(seq, file = hOUT)
            

    else:

        seq = utils.get_seq(reference_genome, juncChr1, int(juncPos1) - validate_sequence_length, int(juncPos1) + validate_sequence_length) 

        """
        seq = ""
        for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(int(juncPos1) - validate_sequence_length) + "-" + str(int(juncPos1) + validate_sequence_length)):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()
        seq = seq.replace('>', '')
        seq = seq.replace(juncChr1 + ":" + str(int(juncPos1) - validate_sequence_length) + "-" + str(int(juncPos1) + validate_sequence_length), '')
        """

        print('>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_ref1", file = hOUT)
        print(seq, file = hOUT)

        seq = utils.get_seq(reference_genome, juncChr2, int(juncPos2) - validate_sequence_length, int(juncPos2) + validate_sequence_length)
        """
        seq = ""
        for item in pysam.faidx(reference_genome, juncChr2 + ":" + str(int(juncPos2) - validate_sequence_length) + "-" + str(int(juncPos2) + validate_sequence_length)):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()
        seq = seq.replace('>', '')
        seq = seq.replace(juncChr2 + ":" + str(int(juncPos2) - validate_sequence_length) + "-" + str(int(juncPos2) + validate_sequence_length), '')
        """

        print('>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_ref2", file = hOUT)
        print(seq, file = hOUT)

        seq = ""
        if juncDir1 == "+":
            tseq = utils.get_seq(reference_genome, juncChr1, int(juncPos1) - validate_sequence_length, int(juncPos1))
            """
            tseq = ""
            for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(int(juncPos1) - validate_sequence_length) + "-" + str(juncPos1)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
            tseq = tseq.replace('>', '')
            tseq = tseq.replace(juncChr1 + ":" + str(int(juncPos1) - validate_sequence_length) + "-" + str(juncPos1), '')
            """
        else:
            tseq = utils.get_seq(reference_genome, juncChr1, int(juncPos1), int(juncPos1) + validate_sequence_length)
            """
            tseq = ""
            for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(juncPos1) + "-" + str(int(juncPos1) + validate_sequence_length)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
            tseq = tseq.replace('>', '')
            tseq = tseq.replace(juncChr1 + ":" + str(juncPos1) + "-" + str(int(juncPos1) + validate_sequence_length), '')
            """
            tseq = utils.reverseComplement(tseq)

        if juncDir1 == "+":
            seq = tseq + juncSeq
        else:
            seq = tseq + utils.reverseComplement(juncSeq)

        if juncDir2 == "-":
            tseq = utils.get_seq(reference_genome, juncChr2, int(juncPos2), int(juncPos2) + validate_sequence_length)
            """
            tseq = "" 
            for item in pysam.faidx(reference_genome, juncChr2 + ":" + str(juncPos2) + "-" + str(int(juncPos2) + validate_sequence_length)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
            tseq = tseq.replace('>', '')
            tseq = tseq.replace(juncChr2 + ":" + str(juncPos2) + "-" + str(int(juncPos2) + validate_sequence_length), '')
            """
        else:
            tseq = utils.get_seq(reference_genome, juncChr2, int(juncPos2) - validate_sequence_length, int(juncPos2))
            """
            tseq = ""
            for item in pysam.faidx(reference_genome, juncChr2 + ":" + str(int(juncPos2) - validate_sequence_length) + "-" + str(juncPos2)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
            tseq = tseq.replace('>', '')
            tseq = tseq.replace(juncChr2 + ":" + str(int(juncPos2) - validate_sequence_length) + "-" + str(juncPos2), '')
            """
            tseq = utils.reverseComplement(tseq)

        seq = seq + tseq

        print('>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_alt", file = hOUT)
        print(seq, file = hOUT)

    
    hOUT.close()



def checkScore(align):

    tempScore = 100
    if len(align) >= 2:
        for i1 in range(0, len(align) - 1):
            for i2 in range(i1 + 1, len(align)):

                if align[i1][1] <= align[i2][1] and align[i1][2] == "+" and align[i2][2] == "-":
                    tempScore = min(tempScore, align[i1][0] + align[i2][0])

                if align[i2][1] <= align[i1][1] and align[i2][2] == "+" and align[i1][2] == "-":
                    tempScore = min(tempScore, align[i1][0] + align[i2][0])

    return(tempScore)



def summarizeRefAlt(inputFile, ITDFlag):

    """
        note:
        Now the definition of "reference-prefer" read pair is those align at least 5 base than to alternative sequence. 
        but current definition generates problem when detecting mid-range tandem duplication on repeat sequences (5,115279667,-,5,115280072 (ATL-15T))
        because in this case, read pairs that should prefer the reference sequence do not prefer it significantly to alternative base...
        One remedy for this may be to 
        1. in advance, we should remove the read pairs whose source (reference or altenative) are unclear
        2. we define "class 2 reference read pair", and use them with "class 1 reference read pair".
    """
 
    hIN = open(inputFile, 'r')

    numOther = 0
    numAlt = 0
    numRef = 0

    # ref_ID = []
    # alt_ID = []
    # other_ID = []

    tempID = ""
    tempAlt = []
    tempRef1 = []
    tempRef2 = []
    tempRef = []
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        if F[0].isdigit() == False: continue

        # remove the read pair num info ("/1", "/2") 
        F[9] = F[9][0:-2]
        if tempID != F[9]:
            if tempID != "":

                if tempID == "HWI-ST1021:128:C1K77ACXX:6:2303:6851:109312":
                    pass

                tempAltNM = checkScore(tempAlt)
                tempRefNM = min(checkScore(tempRef1), checkScore(tempRef2), checkScore(tempRef))

                if tempAltNM >= 10 and tempRefNM >= 10:
                    numOther = numOther + 1
                    # other_ID.append(tempID)
                elif tempAltNM < tempRefNM - 5:
                    numAlt = numAlt + 1
                    # alt_ID.append(tempID)
                elif ITDFlag == 0 and tempRefNM < tempAltNM - 5 or ITDFlag == 1 and tempRefNM <= tempAltNM:
                    numRef = numRef + 1
                    # print tempID
                   # ref_ID.append(tempID)


            tempID = F[9] 
            tempAlt = []
            tempRef1 = []
            tempRef2 = []
            tempRef = []

        tNM = int(F[10]) - int(F[0]) + int(F[5]) + int(F[7])
        tpos = int(F[15])
        tdir = F[8]

        if F[13][-3:] == "alt":
            tempAlt.append((tNM, tpos, tdir))
        elif F[13][-4:] == "ref1":
            tempRef1.append((tNM, tpos, tdir))
        elif F[13][-4:] == "ref2":
            tempRef2.append((tNM, tpos, tdir))
        elif F[13][-3:] == "ref":
            tempRef.append((tNM, tpos, tdir))

    hIN.close()

    tempAltNM = checkScore(tempAlt)
    tempRefNM = min(checkScore(tempRef1), checkScore(tempRef2), checkScore(tempRef))

    if tempAltNM >= 10 and tempRefNM >= 10:
        numOther = numOther + 1
        # other_ID.append(tempID)
    elif tempAltNM < tempRefNM - 5:
        numAlt = numAlt + 1
        # alt_ID.append(tempID)
    elif ITDFlag == 0 and tempRefNM < tempAltNM - 5 or ITDFlag == 1 and tempRefNM <= tempAltNM:
        numRef = numRef + 1
        # ref_ID.append(tempID)

    """
    print "ref"
    print '\n'.join(ref_ID)
    print "alt"
    print '\n'.join(alt_ID)
    print "other"
    print '\n'.join(other_ID)
    """

    return([numRef, numAlt])


