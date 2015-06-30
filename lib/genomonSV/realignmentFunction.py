#!/usr/bin/env python


def extractSVReadPairs(bamFilePath, outputFilePath, Params, juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2):

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

    bamfile = pysam.Samfile(bamFilePath, 'rb')
    max_depth = Params["max_depth"]
    search_length = Params["search_length"]
    search_margin = Params["search_margin"]

    # if the #sequence read is over the `maxDepth`, then that key is ignored
    depthFlag = 0
    if bamfile.count(juncChr1, int(juncPos1) - 1, int(juncPos1) + 1) >= max_depth: depthFlag = 1
    if bamfile.count(juncChr2, int(juncPos2) - 1, int(juncPos2) + 1) >= max_depth: depthFlag = 1
    if depthFlag == 1: 
        sys.exit(27)

    hOUT = open(outputFilePath, 'w')

    readID2exist = {}    
    for read in bamfile.fetch(juncChr1, int(juncPos1) - search_length, int(juncPos1) + search_length):

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
        if pos_current - search_margin <= juncPos1 <= (read.aend - 1) + search_margin:
            readID2exist[read.qname] = 1
    
        # the read pair covers break point
        if chr_pair == juncChr1 and pos_current <= juncPos1 <= pos_pair and dir_current == "+" and dir_pair == "-":
            readID2exist[read.qname] = 1

        # the read pair covers break point
        if chr_pair == juncChr2:
            juncFlag = 0
            if juncDir1 == "+" and juncDir2 == "+" and pos_current <= juncPos1 and pos_pair <= juncPos2: juncFlag = 1
            if juncDir1 == "+" and juncDir2 == "-" and pos_current <= juncPos1 and pos_pair >= juncPos2: juncFlag = 1
            if juncDir1 == "-" and juncDir2 == "+" and pos_current >= juncPos1 and pos_pair <= juncPos2: juncFlag = 1
            if juncDir1 == "-" and juncDir2 == "-" and pos_current >= juncPos1 and pos_pair >= juncPos2: juncFlag = 1

            if juncFlag == 1:  
                readID2exist[read.qname] = 1


    for read in bamfile.fetch(juncChr2, int(juncPos2) - search_length, int(juncPos2) + search_length):
        
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
        if pos_current - search_margin <= juncPos2 <= (read.aend - 1) + search_margin:
            readID2exist[read.qname] = 1
                
        # the read pair covers break point
        if chr_pair == juncChr2 and pos_current <= juncPos2 <= pos_pair and dir_current == "+" and dir_pair == "-":
            readID2exist[read.qname] = 1
                
        # the read pair covers break point
        if chr_pair == juncChr1:
            juncFlag = 0
            if juncDir2 == "+" and juncDir1 == "+" and pos_current <= juncPos2 and pos_pair <= juncPos1: juncFlag = 1
            if juncDir2 == "+" and juncDir1 == "-" and pos_current <= juncPos2 and pos_pair >= juncPos1: juncFlag = 1
            if juncDir2 == "-" and juncDir1 == "+" and pos_current >= juncPos2 and pos_pair <= juncPos1: juncFlag = 1
            if juncDir2 == "-" and juncDir1 == "-" and pos_current >= juncPos2 and pos_pair >= juncPos1: juncFlag = 1
             
            if juncFlag == 1:
                readID2exist[read.qname] = 1


    readID2seq1 = {}
    readID2seq2 = {}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    for read in bamfile.fetch(juncChr1, int(juncPos1) - search_length, int(juncPos1) + search_length):

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
                tempSeq = "".join(complement.get(base) for base in reversed(str(read.seq)))
            else:
                tempSeq = read.seq
 
            # the first read
            if flags[6] == "1":
                readID2seq1[read.qname] = tempSeq
            else:
                readID2seq2[read.qname] = tempSeq


    for read in bamfile.fetch(juncChr2, int(juncPos2) - search_length, int(juncPos2) + search_length):

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
                tempSeq = "".join(complement.get(base) for base in reversed(str(read.seq)))
            else:
                tempSeq = read.seq

            # the first read
            if flags[6] == "1":
                readID2seq1[read.qname] = tempSeq
            else:
                readID2seq2[read.qname] = tempSeq


    for readID in readID2seq1:
        if readID in readID2seq2:
            print >> hOUT, '>' + readID + '/1'
            print >> hOUT, readID2seq1[readID]
            print >> hOUT, '>' + readID + '/2'
            print >> hOUT, readID2seq2[readID]

    bamfile.close()
    hOUT.close()


def getRefAltForSV(outputFilePath, Params, juncChr1, juncPos1, juncDir1, juncChr2, juncPos2, juncDir2, juncSeq):

    """
        for short SV (mid-range (<= split_refernece_thres bp) deletion, tandem duplication), we get the two sequence
        for large SV (> split_refernece_thres bp), we get three sequence (one joint sequence by two break points, and two reference sequences around the break points)

        the only concern is short inversion... (are there some somatic short inversion?)
        however, this will be filtered beforehand by the "cover filter", and maybe we have to give up detecting these class of SVs.

    """

    reference_genome = Params["reference_genome"]
    split_refernece_thres = Params["split_refernece_thres"]
    validate_sequence_length = Params["validate_sequence_length"]

    hOUT = open(outputFilePath, 'r')

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    if juncSeq == "---": juncSeq = ""

    # for mid-range deletion or tandem duplication
    if juncChr1 == juncChr2 and abs(juncPos1 - juncPos2) <= split_refernece_thres and juncDir1 != juncDir2:

        seq = ""
        for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(juncPos1 - validate_sequence_length) + "-" + str(juncPos2 + validate_sequence_length)):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()

        print >> hOUT, '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_ref"
        print >> hOUT, seq

        # for mid-range deletion
        if juncDir1 == "+" and juncDir2 == "-":

            seq = ""
            for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(juncPos1 - validate_sequence_length) + "-" + str(juncPos1)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()

            seq = seq + juncSeq

            for item in pysam.faidx(reference_genome, juncChr2 + ":" + str(juncPos2) + "-" + str(juncPos2 + validate_sequence_length)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()

            print >> hOUT, '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_alt"
            print >> hOUT, seq

        # for mid-range tandem duplication
        else:
            seq = "" 
            for item in pysam.faidx(reference_genome, juncChr2 + ":" + str(juncPos2 - validate_sequence_length) + "-" + str(juncPos2)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()
            
            seq = seq + juncSeq

            for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(juncPos1) + "-" + str(juncPos1 + validate_sequence_length)):
                if item[0] == ">": continue
                seq = seq + item.rstrip('\n').upper()
            
            print >> hOUT, '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_alt"
            print >> hOUT, seq


    else:

        seq = ""
        for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(juncPos1 - validate_sequence_length) + "-" + str(juncPos1 + validate_sequence_length)):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()

        print >> hOUT, '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_ref1"
        print >> hOUT, seq

        seq = ""
        for item in pysam.faidx(reference_genome, juncChr2 + ":" + str(juncPos2 - validate_sequence_length) + "-" + str(juncPos2 + validate_sequence_length)):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()
            
        print >> hOUT, '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_ref2"
        print >> hOUT, seq


        seq = ""
        if juncDir1 == "+":
            tseq = ""
            for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(juncPos1 - validate_sequence_length) + "-" + str(juncPos1)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
        else:
            tseq = ""
            for item in pysam.faidx(reference_genome, juncChr1 + ":" + str(juncPos1) + "-" + str(juncPos1 + validate_sequence_length)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
            tseq = "".join(complement.get(base) for base in reversed(tseq))

        seq = tseq + juncSeq

        if juncDir2 == "-":
            tseq = "" 
            for item in pysam.faidx(reference_genome, juncChr2 + ":" + str(juncPos2) + "-" + str(juncPos2 + validate_sequence_length)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
        else:
            tseq = ""
            for item in pysam.faidx(reference_genome, juncChr2 + ":" + str(juncPos2 - validate_sequence_length) + "-" + str(juncPos2)):
                if item[0] == ">": continue
                tseq = tseq + item.rstrip('\n').upper()
            tseq = "".join(complement.get(base) for base in reversed(tseq))

        seq = seq + tseq

        print >> hOUT, '>' + ','.join([juncChr1, str(juncPos1), juncDir1, juncChr2, str(juncPos2), juncDir2]) + "_alt"
        print >> hOUT, seq

    
    hOUT.close()


