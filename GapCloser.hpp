#ifndef GAPCLOSER_HPP_
#define GAPCLOSER_HPP_

#include <iostream>
#include <fstream>

#include "Common.hpp"
#include "ConsensusConfig.hpp"
#include "ContigForFill.hpp"
#include "ContigTable.hpp"
#include "ContigAssembler.hpp"
#include "ReadMatrix.hpp"
/*
   Fill gaps in scaffolds in parallel, one thread handling one scaffold.
   */
class GapCloser : public ContigAssembler
{
    protected:

        std::ofstream foutFill;
        const ContigTable & contigTable;
        float deviation;
        Len_t endNumLen;
        Len_t mismatchLen;

        Len_t minSearchLen;
        Len_t mapReadsRange;
        Len_t actualGapSum;
        Len_t extendGapSum;
        Len_t actualGapCount;
        Len_t finishGapCount;

    public:

        GapCloser(char* _outfile, std::ofstream& _fout
                , ReadAccessor& _readAccessor
                , PairInfo const& _pairInfo
                , ContigTable const& _contigTable
                , Len_t _threadSum
                , float _deviation=0.5
                , Len_t _endNumLen=10
                , Len_t _mismatchLen=5 
            /*  , Len_t _maxReadLength=100 */
            /*  , Short_Len_t _overlapMode=fixedOverlapMode */
                , Short_Len_t _overlapParam=25) : 

            ContigAssembler(_outfile, _fout, _readAccessor
                    , _pairInfo, _threadSum/*, _maxReadLength */
                   /* , _overlapMode */
                    , _overlapParam),
            contigTable(_contigTable), 
            deviation(_deviation), 
            endNumLen(_endNumLen), 
            mismatchLen(_mismatchLen), 
            minSearchLen(1000), 
            mapReadsRange(1000), 
            actualGapSum(0), 
            extendGapSum(0), 
            actualGapCount(0), 
            finishGapCount(0)
            {
                /* reExtendLength = _overlapParam;*/
            }

    public:

        void assemble() {

            std::cout << ">>>>>>>>>>assembling<<<<<<<<<<" << std::endl;
            std::cout << std::endl;
            time_t total_start_time=time(NULL);

            char fillName[MAX_STRING_LEN];
            strcpy(fillName, outfile);
            strcat(fillName,".fill");
            foutFill.open(fillName);

            createThread();
            //TODO
            //		assembleInThread(0);

            foutFill.close();

            std::cout << "---------------------------------" << std::endl;
            std::cout << "actual gap sum: " << actualGapSum << std::endl;
            std::cout << "extend gap sum: " << extendGapSum << std::endl;
            std::cout << "actual gap count: " << actualGapCount << std::endl;
            std::cout << "finish gap count: " << finishGapCount << std::endl;

            std::cout << "spent total time: " << time(NULL)-total_start_time << "s"<< std::endl;
            std::cout << ">>>>>>>>>>assembling finished<<<<<<<<<<" << std::endl;
            std::cout << std::endl;
        }



    protected:

        struct Arg 
        {
            GapCloser* pGapCloser;
            Len_t iThread;
        };

        void createThread() {

            for (Len_t i=0; i<threadSum; i++) {
                Arg* pArg = new Arg;
                pArg->pGapCloser = this;
                pArg->iThread = i;
                pthread_create(&threadsID[i], NULL, (void *(*)(void *))threadFunc, (void*)pArg);
            }

            for (Len_t i=0; i<threadSum; i++) {
                pthread_join(threadsID[i], NULL);
            }
        }

        static void threadFunc(Arg* pArg) {

            GapCloser& contigAssemblerForFill = *(pArg->pGapCloser);
            Len_t iThread = pArg->iThread;
            delete pArg;

            contigAssemblerForFill.assembleInThread(iThread);
        }

        struct GapResultInfo
        {
            Len_t start;
            Len_t end;
            Len_t flag;
            /* delete quality code */
            /*
               Len_t quality;
               */
            ExtendInfo extendInfo;

            GapResultInfo() : 
                start(0), 
                end(0), 
                flag(0) 
                /* delete quality code */
                /*
                   ,quality(0)
                   */
            {
            }
        };


        /* assembly in one thread */
        void assembleInThread(Len_t ) {

            Len_t actualGapSumInThread=0;
            Len_t extendGapSumInThread=0;

            Len_t actualGapCountInThread=0;
            Len_t finishGapCountInThread=0;

            for (Len_t i=1; i<=contigTable.getContigsSum(); i++) {

                ContigForFill& contig = contigTable.getContigs()[i];

                // get an unused scaffold
                bool isUsed;
                pthread_mutex_lock(&mutexNumberOfContigs);
                isUsed = contig.isUsed();
                if (!isUsed) contig.setUsedFlag();
                pthread_mutex_unlock(&mutexNumberOfContigs);
                if (isUsed) continue;

                pthread_mutex_lock(&mutexNumberOfContigs);
                numberOfContigs++;
                std::cout << "constructing " << numberOfContigs << " scaffold" <<std::endl;
                pthread_mutex_unlock(&mutexNumberOfContigs);


                actualGapSumInThread += contig.getGapSum();
                actualGapCountInThread += contig.getGapCount();

                LinkedList<TightString*> const& sequences = contig.getSequences();
                LinkedList<GapInfo> const& gapsList = contig.getGaps();

                Len_t contigCount = sequences.getCount();
                Len_t gapCount = gapsList.getCount();
                Contig* contigs = new Contig[contigCount];
                GapInfo* gaps = new GapInfo[gapCount];

                // map reads to contigs
                Len_t j=0;
                ListElement<TightString*> const* ptr;
                for (ptr = sequences.getHead(); ptr != 0; ptr = ptr->getNext()) {

                    TightString const& seq = *(ptr->getDatum());
                    Contig contig(seq);
                    if (contig.getLength() > mapReadsRange*2) {
                        contig.mapReads(0, mapReadsRange, overlapParam);
                        contig.mapReads(contig.getLength()-mapReadsRange, mapReadsRange, overlapParam);
                    }
                    else {
                        contig.mapReads(0, contig.getLength(), overlapParam);
                    }
                    contigs[j] = contig;
                    j++;
                }

                // get gaps
                j=0;
                ListElement<GapInfo> const* ptrGap;
                for (ptrGap = gapsList.getHead(); ptrGap != 0; ptrGap = ptrGap->getNext()) {

                    GapInfo const& gap = ptrGap->getDatum();
                    gaps[j] = gap;
                    j++;
                }

                Contig* contigsResult = new Contig[contigCount];
                GapInfo* gapsResult = new GapInfo[gapCount];

                // fill gaps in forward direction
                fillGaps(contigs, gaps, contigsResult, gapsResult, contigCount, gapCount);

                delete [] contigs;
                delete [] gaps;
                Contig* contigsReverse = new Contig[contigCount];
                GapInfo* gapsReverse = new GapInfo[gapCount];

                //reverse contigs and gaps
                for (j=0; j<gapCount; j++) {

                    Contig const& contig = contigsResult[j];
                    GapInfo const& gap = gapsResult[j];

                    Contig contigReverse;
                    contig.reverse(contigReverse, 0, contig.getLength());

                    contigsReverse[contigCount-j-1] = contigReverse;

                    GapInfo gapReverse;
                    gapReverse.length = gap.length;
                    gapReverse.extendInfo.leftLen = gap.extendInfo.rightLen;
                    gapReverse.extendInfo.rightLen = gap.extendInfo.leftLen;
                    gapReverse.isFilled = gap.isFilled;
                    /* delete quality code */
                    /*
                       gapReverse.quality = gap.quality;
                       */
                    if (!gapReverse.isFilled) {

                        //Len_t sameFlag = 1;	
                        for (Len_t j=0; j<mismatchLen; j++) {

                            if ((int)(contigReverse.getLength() - j) >= (int)endNumLen) {

                                Number_t numReverse;
                                contigReverse.getTightString().readTightStringFragment(j, j+endNumLen, numReverse);

                                if (!gapReverse.endNumHash.find(numReverse)) {
                                    gapReverse.endNumHash[numReverse] = j;
                                    ; //sameFlag = 0;	
                                }
                                else {
                                    gapReverse.endNumHash[numReverse] = (Len_t)-1;	
                                    //sameKmer++;
                                }
                            }
                        }

                        //if (contigReverse.getLength() < endNumLen) {
                        //    sameFlag = 0;
                        //}
                        //if (sameFlag == 1) {
                        //    sameEnd++;
                        //}
                    }

                    gapsReverse[gapCount-j-1] = gapReverse;
                }

                for (; j<contigCount; j++) {

                    Contig const& contig = contigsResult[j];
                    Contig contigReverse;
                    contig.reverse(contigReverse, 0, contig.getLength());
                    contigsReverse[contigCount-j-1] = contigReverse;
                }

                // add unique anchor of next contig to gap releated hash
                for (j=0; j<gapCount; j++) {
                    Contig const& contigReverse = contigsReverse[j+1];
                    GapInfo & gapReverse = gapsReverse[j];

                    if (gapReverse.isFilled) {
                        continue;
                    }

                    Len_t uniqueNum = 0;
                    Len_t k=0;
                    for ( ; k<mismatchLen; k++) {
                        if ((int)(contigReverse.getLength() - k) >= (int)endNumLen) {
                            Number_t numReverse;
                            contigReverse.getTightString().readTightStringFragment(k, k+endNumLen, numReverse);

                            if (gapReverse.endNumHash.find(numReverse) && (gapReverse.endNumHash[numReverse] != (Len_t)-1)) {
                                uniqueNum = 1;
                                break;
                            }
                        }
                        else {
                            uniqueNum = 1;
                        }
                    }

                    if (uniqueNum == 1) { // at least an unique anchor was found
                        continue;
                    }
                    else {
                        //Len_t mismatchNum = 0;
                        while ((int)(contigReverse.getLength() - k) >= (int)endNumLen) {
                            Number_t numReverse;
                            contigReverse.getTightString().readTightStringFragment(k, k+endNumLen, numReverse);

                            if (!gapReverse.endNumHash.find(numReverse)) {  // an unique anchor was found
                                gapReverse.endNumHash[numReverse] = k;							
                                break;			
                            }
                            else {
                                gapReverse.endNumHash[numReverse] = (Len_t)-1;
                                //sameKmer++;
                            }

                            k++;
                        }
                    }
                }

                delete [] contigsResult;
                delete [] gapsResult;
                contigsResult = new Contig[contigCount];
                gapsResult = new GapInfo[gapCount];

                //fill gaps in reversed direction
                fillGaps(contigsReverse, gapsReverse, contigsResult, gapsResult, contigCount, gapCount);

                delete [] contigsReverse;
                delete [] gapsReverse;

                //get sequence
                std::string finalSeqReverse;
                for (j=0; j<gapCount; j++) {

                    Contig const& contig = contigsResult[j];
                    //				GapInfo const& gap = gapsResult[j];
                    GapInfo & gap = gapsResult[j];

                    char* strSeq = contig.getTightString().readTightString();
                    if (gap.length<0) {
                        if ((int)(contig.getLength()+gap.length) > 0)
                            finalSeqReverse.append(strSeq, contig.getLength()+gap.length);
                    }
                    else {
                        finalSeqReverse.append(strSeq, contig.getLength());
                    }
                    delete [] strSeq;

                    //fill 'N' string
                    for (int j=0; j<gap.length; j++) {
                        finalSeqReverse += 'N';
                    }

                    if ((!gap.isFilled) && (gap.length==0)) {
                        if (Threshold::NNumber > 0){
                            std::string NSeq(Threshold::NNumber, 'N');
                            finalSeqReverse += NSeq;
                            gap.length = Threshold::NNumber;
                        }
                    }
                }

                // get last contig sequence
                for (; j<contigCount; j++) {

                    Contig const& contig = contigsResult[j];
                    char* strSeq = contig.getTightString().readTightString();
                    finalSeqReverse.append(strSeq, contig.getLength());
                    delete [] strSeq;
                }

                char* finalSeq = reverseComplement(finalSeqReverse.c_str(),finalSeqReverse.size());

                //get extend info
                GapResultInfo gapResultInfo;
                LinkedList<GapResultInfo> gapResultsInfo;
                Len_t pos = 0;
                for (j=0; j<gapCount; j++) {

                    Contig const& contig = contigsResult[j];
                    GapInfo const& gap = gapsResult[j];

                    pos += contig.getLength();
                    //				if (gap.length<0) {
                    //					pos += gap.length;
                    //				}

                    gapResultInfo.end = finalSeqReverse.size() - (pos-gap.extendInfo.leftLen);
                    gapResultInfo.start = gapResultInfo.end - (gap.extendInfo.leftLen+gap.length+gap.extendInfo.rightLen);
                    gapResultInfo.extendInfo.leftLen = gap.extendInfo.rightLen;
                    gapResultInfo.extendInfo.rightLen = gap.extendInfo.leftLen;
                    gapResultInfo.flag = gap.isFilled;
                    /* delete quality code */
                    /*
                       gapResultInfo.quality = gap.quality;
                       */
                    gapResultsInfo.prepend(gapResultInfo);

                    extendGapSumInThread += gap.extendInfo.leftLen+gap.extendInfo.rightLen;
                    if (gap.isFilled) {
                        finishGapCountInThread++;
                    }

                    pos += gap.length;
                }

                delete [] contigsResult;
                delete [] gapsResult;

                pthread_mutex_lock(&mutexOutput);
                //output sequence, fasta format
                fout << contig.getName() << std::endl;
                //			fout << finalSeq << std::endl;

                //			Len_t j;
                Len_t finalSeqLen = strlen(finalSeq);
                for (j=0; j<finalSeqLen;) {

                    fout << finalSeq[j];
                    j++;
                    if (j%100==0)
                        fout << std::endl;
                }
                if (j%100)
                    fout << std::endl;

                //output extend info
                foutFill << contig.getName() << std::endl;

                ptrGap = gapsList.getHead();
                ListElement<GapResultInfo> const* ptrGapResult;
                for (ptrGapResult = gapResultsInfo.getHead(); ptrGapResult != 0; ptrGapResult = ptrGapResult->getNext()) {

                    gapResultInfo = ptrGapResult->getDatum();
                    //				foutFill << gapResultInfo.start << "\t" << gapResultInfo.end << "\t" << gapResultInfo.extendInfo.leftLen << "\t" <<gapResultInfo.extendInfo.rightLen << "\t" <<gapResultInfo.flag << "\t" <<gapResultInfo.quality << std::endl;
                    foutFill << gapResultInfo.start << "\t" << gapResultInfo.end << "\t" << gapResultInfo.extendInfo.leftLen << "\t" <<gapResultInfo.extendInfo.rightLen << "\t" <<gapResultInfo.flag << "\t" 
                        /* delete quality code */
                        /*
                           <<gapResultInfo.quality
                           */
                        << "\t" <<ptrGap->getDatum().length << "\t" << gapResultInfo.end - gapResultInfo.start << std::endl;

                    ptrGap = ptrGap->getNext();
                }

                pthread_mutex_unlock(&mutexOutput);

                delete [] finalSeq;
            }

            pthread_mutex_lock(&mutexNumberOfContigs);
            actualGapSum += actualGapSumInThread;
            extendGapSum += extendGapSumInThread;
            actualGapCount += actualGapCountInThread;
            finishGapCount += finishGapCountInThread;
            pthread_mutex_unlock(&mutexNumberOfContigs);
        }

        // fucntion: 	 fillGaps
        // description: fill gaps in scaffold one by one in direction of left to right
        // input:		 contigs -- contigs before gap filling
        //			 gaps -- gaps before gap filling
        //			 contigCount -- contig number
        //			 gapCount -- gap number
        // output:	 contigsResult -- contigs after gap filling
        //			 gapsResult -- gaps after gap filling
        void fillGaps(
                Contig const* contigs, 
                GapInfo const* gaps, 
                Contig* contigsResult, 
                GapInfo* gapsResult, 
                Len_t contigCount, 
                Len_t gapCount
                ) {

            Len_t i;
            for (i=0; i<gapCount; i++) {

                Contig const& contig = contigs[i];
                GapInfo const& gap = gaps[i];
                Len_t j = i<gapCount-1 ? i+1: i;
                Contig const& nextContig = contigs[j];

                Len_t notNCtgLen = 0;		
                if (!gap.isFilled ){		

                    Contig contigFill;


                    if ( (contig.getLength() < mapReadsRange) && i>0 ) {  // append previous contig at the beginning of current contig

                        //calculate the length
                        int sum=0;
                        sum += contig.getLength();
                        Len_t redundantSum = 0;
                        bool appendGap = false;
                        int j;

                        notNCtgLen += contig.getLength();		
                        Len_t foundNotFilledGap = 0;			

                        for(j=i-1; j>=0; j--) {

                            if (!gapsResult[j].isFilled){
                                foundNotFilledGap = 1;
                            }

                            sum += gapsResult[j].length;
                            if (sum>=(int)mapReadsRange) {
                                redundantSum = sum - mapReadsRange;
                                appendGap = true;
                                break;
                            }
                            if (foundNotFilledGap == 0){
                                if (gapsResult[j].length < 0){
                                    if ((int)contigsResult[j].getLength() + gapsResult[j].length > 0){
                                        notNCtgLen += contigsResult[j].getLength() + gapsResult[j].length;
                                    }
                                }
                                else{

                                    notNCtgLen += contigsResult[j].getLength();
                                }
                            }
                            sum += contigsResult[j].getLength();
                            if (sum>=(int)mapReadsRange) {
                                redundantSum = sum - mapReadsRange;
                                appendGap = false;
                                break;
                            }
                        }
                        if (j<0) j=0;

                        //then append
                        if (appendGap) {
                            if (!gapsResult[j].isFilled) {
                                TightString tStrAppend(gapsResult[j].length-redundantSum);
                                contigFill.append(tStrAppend, gapsResult[j].length-redundantSum);
                                redundantSum = 0;
                            }
                            j++;
                        }

                        for (; j<(int)i; j++) {

                            if (gapsResult[j].length<0) {
                                if ((int)(contigsResult[j].getLength()+gapsResult[j].length-redundantSum) > 0)
                                    contigFill.append(contigsResult[j], redundantSum, contigsResult[j].getLength()+gapsResult[j].length-redundantSum);
                            }
                            else {
                                contigFill.append(contigsResult[j], redundantSum, contigsResult[j].getLength()-redundantSum);
                            }

                            if (redundantSum) redundantSum = 0;

                            if (!gapsResult[j].isFilled) {
                                TightString tStrAppend(gapsResult[j].length);
                                contigFill.append(tStrAppend, gapsResult[j].length);
                            }
                        }
                        contigFill.append(contig);
                    }
                    else {

                        //					contigFill.append(contig);
                        if (contig.getLength() < mapReadsRange) {
                            contigFill.append(contig);

                            notNCtgLen = contig.getLength();		
                        }
                        else {
                            contigFill.append(contig, contig.getLength()-mapReadsRange, mapReadsRange);

                            notNCtgLen = mapReadsRange;	
                        }
                    }
                    if (notNCtgLen > mapReadsRange) {
                        notNCtgLen = mapReadsRange;
                    }

                    if (notNCtgLen >= overlapParam) {

                        Contig gapContig;
                        GapInfo gapResult;
                        ConsensusGap( contigFill ,  notNCtgLen 
                                ,  nextContig , gap
                                , gapContig,gapResult );

                        Contig contigResult(contig,gapContig);
                        contigsResult[i] = contigResult;

                        /* delete quality code */
                        /*
                           gapResult.quality = gap.quality;
                           for (Len_t j=0; j<gapContig.getLength(); j++) {
                           gapResult.quality += gapContig.getQuality()[j];
                           }
                           */
                        gapResult.extendInfo.leftLen = gapContig.getLength();
                        gapResult.extendInfo.rightLen = gap.extendInfo.rightLen;
                        gapsResult[i] = gapResult;

                    }
                }
                else if (gap.isFilled || notNCtgLen < overlapParam) {		
                    contigsResult[i] = contig;

                    GapInfo gapResult;
                    gapResult.length = gap.length;
                    gapResult.extendInfo = gap.extendInfo;
                    gapResult.isFilled = gap.isFilled;
                    /* delete quality code */
                    /*
                       gapResult.quality = gap.quality;
                       */
                    gapsResult[i] = gapResult;
                }
            }

            // add the last contig
            for (; i<contigCount; i++) {

                contigsResult[i] = contigs[i];
            }
        }

        // function:    ConsensusGap
        // description: iterately consensus relatative reads to fill gap.
        //                   1) get reads starting with anchor fetched from the end area of contig.
        //                   2) classify reads by PE/barcode/ONT information
        //                   2) choose a reads collection depends on the gap size and other information 
        //                         and consensus those reads.
        //                   3) check overlap between extend contig and next contig
        // input:          contig -- contig on the left side of gap
        //                   gap -- gap to fill
        //                   realContigLen -- length of none N sequence in contig on the left side of gap
        //                   nextContig -- contig on the right side of gap
        //output:         gapContig -- extend sequence in gap
        //                   gapResult -- filled gap result
        void ConsensusGap(Contig& contig
                , Len_t /*realContigLen*/
                , const Contig & nextContig
                , const GapInfo & gap
                , Contig& gapContig
                , GapInfo& gapResult ) 
        {

            int originalLen = contig.getLength();
            // step 0 , clean the barcodes in consensus area.
            ConsensusArea the_area = ConsensusConfig::GetConsensusArea(originalLen);
            BarcodeInfo & barcode_info = contig.getBarodeInfo() ;
            barcode_info.Erase(the_area.left_most_pos_in_contig - 1 );

            // step 1 , loop area consensus .
            ConsensusArea prev_area ;
            int prev_contig_len = contig.getLength() ;
            while (true) 
            {
                ConsensusArea curr_area = ConsensusConfig::GetConsensusArea(contig.getLength());
                if( curr_area == prev_area )
                    break ;
                // step 1.1 gather reads .
                ReadMatrix  readMatrix =  ReadMatrixFactory::GenReadMatrix(contig,curr_area);
                if( readMatrix.is_reads_too_little()  )
                    break ;
                // step 1.2 abstract sub reads .
                readMatrix = readMatrix.GenSubMatrixByGap(contig
                        ,nextContig
                        ,gap);
                // step 1.3 do consensus .
                ConsensusMatrix consensusMatrix = readMatrix.GenConsensusMatrix(contig);
                ConsensusResult consensusResult =  consensusMatrix.GenConsensusResult();

                if( consensusResult.is_consensus_done() )
                {
                    // step 1.5 update contig
                    updateContig(contig , readMatrix,consensusResult);
                    // step 1.6 check gap fill
                    bool gapFill ; int prev_contig_pos ; int next_contig_pos ;
                    std::tie( gapFill , prev_contig_pos , next_contig_pos ) 
                        = ContigTool::IsGapFinish( contig , nextContig ,gap ,prev_contig_len);
                    if( gapFill )
                    {
                        // step 1.7 deal output for fill
                        if (prev_contig_pos > originalLen )
                        {
                            int extend_len = prev_contig_pos - originalLen ;
                            gapContig.append(contig,originalLen , extend_len );
                            gapResult.length = - next_contig_pos ;
                            if (((int)prev_contig_pos +gapResult.length) < 0)
                                gapResult.length = -prev_contig_pos;

                        }
                        else
                        {
                            gapResult.length = prev_contig_pos -originalLen - next_contig_pos;
                            if( originalLen + gapResult.length  < 0 )
                                gapResult.length = -originalLen ;
                        }
                        gapResult.isFilled = true ;
                        break ;
                    }
                }
                else
                    break;

                prev_contig_len = contig.getLength() ;
                prev_area = curr_area ;
            }
            // step 3 , deal output for not fill
            if( (! gapResult.isFilled)
                    &&  (int) contig.getLength() > originalLen )
            {
                int extend_len = contig.getLength() - originalLen ;
                gapContig.append(contig,originalLen , extend_len );
                gapResult.length = extend_len ;
                gapResult.isFilled = false ;
            }
        }

        void updateContig( Contig & contig
                , const ReadMatrix &readMatrix
                ,const ConsensusResult & consensusResult 
                )
        {
            const auto & the_area = readMatrix.Area() ;
            int contigRemainLen = contig.getLength()
                - the_area.consensus_start_pos_in_contig + 1;

            int extend_len = consensusResult.getLength() - contigRemainLen ;

            TightString & contig_seq = contig.getTightString();
            TightString consensus_seq = consensusResult.getTightString() ;
            // step 1 , Update sequence first 
            for( int i = 0 ; i < contigRemainLen ; i ++ )
            { // 1.1 check for replace
                if ( contig_seq[the_area.consensus_start_pos_in_contig+i] 
                        != consensus_seq[i] )
                {
                    contig_seq.writeNucleotideAtPosition(
                            consensus_seq[i]
                            ,the_area.consensus_start_pos_in_contig+i
                            );
                }
            }

            if( extend_len > 0 )
            { // 1.2 extend contig

                TightString extend_seq = 
                    consensusResult.getTightString(
                            contigRemainLen-1
                            ,extend_len ) ;
                contig.append( extend_seq , extend_len ); 
            }

            // Step 2 , Update reads && barcodes
            if ( extend_len > 0 )
            {
                for( int i = 0 ; i < extend_len ; i++ )
                {   // for each point 
                    int contig_pos = 
                        the_area.left_most_pos_in_contig + i ;
                    if( readMatrix.HasReads(contig_pos) )
                    {
                        auto & keeper = 
                            readMatrix.GetInfoKeeper(contig_pos);
                        UpdateAInfoKeeper2Contig(contig
                                ,keeper,contig_pos);
                    }
                }
            }
        }

        void UpdateAInfoKeeper2Contig( Contig & contig ,
                const PositionInfoKeeper & keeper
                , int contig_pos /*1base*/
                )
        {
            for( const auto & pair : keeper )
            {
                const  auto & reads = pair.second.reads ;
                for( const auto & read : reads )
                {   // for each read
                    if( ContigTool::IsReadMatchContig
                            ( read , contig ,contig_pos ) )
                    {
                        ContigTool::AddReadIntoContig(
                                contig,read,contig_pos - 1 );
                    }
                }
            }
        }
};

#endif /*CONTIGASSEMBLERFORFILL_HPP_*/
