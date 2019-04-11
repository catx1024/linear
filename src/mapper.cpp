#include <iostream>
#include <fstream>
#include <ctime>
#include <seqan/arg_parse.h>
#include "cord.h"
#include "chain_map.h"
#include "pmpfinder.h"
#include "gap.h"
#include "align_interface.h"
#include "mapper.h"

using namespace seqan; 

MapParm parm1 ( 
        base_block_size_,     //blockSize,
        //Const_::_DELTA,          //delta(Const_::_DELTA),
        64,                          //delta
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        1,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        10,                      //listN
        20,                      //listN2
        15,                     //alpha(Const_::_ALPHA),
        5,                      //alpha2 for complex mapping 
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.5,                     //rcThr(0.75)
        0.7,                     //cordThr length of cord < cordThr are abandone
        0.7,                     //senthr: perfrom next filter on cords of length < senthr 
        0.1                      //clsthr: thread of cluster
); 

//normal
MapParm parm0 ( 
        base_block_size_,     //blockSize,
        base_delta_,          //delta(Const_::_DELTA),
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        0,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        10,                      //listN
        20,                      //listN2
        15,                     //alpha(Const_::_ALPHA),
        5,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.5,                     //rcThr(0.8)
        0.2,                     //cordThr length of cord < cordThr are abandoned
        0.2,                     //senthr: length of cord < senthr are erased duing path
        0.1                      //clsthr: thread of cluster

); 

//sensitive
MapParm parm2 ( 
        base_block_size_,     //blockSize,
        //Const_::_DELTA,          //delta(Const_::_DELTA),
        64,                      //delta
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        2,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        4,                      //listN
        50,                      //listN2
        0.65,                     //alpha(Const_::_ALPHA),
        0.5,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.8,                     //rcThr(0.75)
        0.8,                      //cordThr length of cord < cordThr are abandoned
        0.8,                     //senthr: length of cord < senthr are erased duing path
        0.1                      //clsthr: thread of cluster

        
); 


MapParm parmt ( 
        base_block_size_,     //blockSize,
        base_delta_,          //delta(Const_::_DELTA),
        base_threshold_,     //threshold(Const_::_THRESHOLD),
        base_kmer_step_,       //kmerStep(Const_::_KMERSTEP),
        base_shape_len_,      //shapeLen(Const_::_SHAPELEN),
        0,                      //senstivity(0),
        0,                      //anchorDeltaThr(),
        1000,                   //minReadLen(1000),
        2,                         //listN
        2,                      //listN2
        0.75,                     //alpha(Const_::_ALPHA),
        0.65,                      //alpha2 for complex mapping
        0.02,                    //anchorLenThr(0.02),    anchors with lenghth > this parameter is pushed into the queue
        0.8,                     //rcThr(0.8)
        0.8,                       //cordThr length of cord < cordThr are abandoned
        0.8,                     //senthr: length of cord < senthr are erased duing path
        0.1                      //clsthr: thread of cluster
        
); 

Mapper::Mapper(Options & options):
    record(options),
    qIndex(genomes()), 
    of(toCString(options.getOutputPath()))
{
    outputPrefix = getFileName(getFileName(options.getReadPath()), '.', 0);
    switch (options.sensitivity)
    {
        case 0: 
        {
            parm = parm0; //normal
            break;
        }
        case 1:
        {
            parm =  parm1; //fast
            break;
        }
        case 2:
        {
            parm = parm2; //sensitive
            break;
        }
    }
    _thread = options.thread;
}

int Mapper::createIndex(bool efficient)
{
    std::cerr << ">>Create index \r";
    createHIndex(genomes(), qIndex, _thread, efficient);
    return 0;
}

int print_align_sam_header_ (StringSet<CharString> & genomesId,
                             StringSet<String<Dna5> > & genomes,
                             std::ofstream & of
                            )
{
    of << "@HD\tVN:1.6\n";
    for (int k = 0; k < length(genomesId); k++)
    {
        of << "@SQ\tSN:" << genomesId[k] << "\tLN:" << length(genomes[k]) << "\n";
    }
    of << "@PG\tPN:" << "Linear\n";
}
int print_align_sam_record_(StringSet<String< BamAlignmentRecord > > & records, 
                     StringSet<String<uint64_t> > & cordSet,
                     StringSet<CharString> & readsId, 
                     StringSet<CharString> & genomesId,
                     std::ofstream & of
                    )
{
    for (int i = 0; i < length(records); i++)
    {
        for (int j = 0; j < length(records[i]); j++)
        {
            records[i][j].qName = readsId[i];
            CharString g_id = genomesId[records[i][j].rID];
            writeSam(of, records[i][j], g_id);
        }
    }
}
int print_align_sam_record_(StringSet<String< BamAlignmentRecordLink> > & records, 
                     StringSet<String<uint64_t> > & cordSet,
                     StringSet<CharString> & readsId, 
                     StringSet<CharString> & genomesId,
                     std::ofstream & of
                    )
{
    for (int i = 0; i < length(records); i++)
    {
        for (int j = 0; j < length(records[i]); j++)
        {
            records[i][j].qName = readsId[i];
            CharString g_id = genomesId[records[i][j].rID];
            int dt = writeSam(of, records[i], j, g_id);
            //<<<debug 
        }

        for (int j = 0; j < length(records[i]); j++)
        {
            std::pair<int,int> lens;
            lens = countCigar(records[i][j].cigar);
            if (records[i][j].isEnd())
                std::cout << "\n";
            std::cout << "cigarLen " << lens.first << " " << lens.second << "\n";
            //>>>debug
        }
    }
}
int print_align_sam (Mapper & mapper)
{
    std::string filePath = mapper.getOutputPrefix() + ".sam";
    mapper.getOf().open(toCString(filePath));
    print_align_sam_header_(mapper.genomesId(), 
                            mapper.genomes(),
                            mapper.getOf()
                           );
    print_align_sam_record_(mapper.getBamRecords(),
                            mapper.cords(),
                            mapper.readsId(),
                            mapper.genomesId(),
                            mapper.getOf()
                           ); 
    mapper.getOf().close();
}

void Mapper::printCordsRaw()
{
    double time = sysTime();
    //unsigned strand;
    unsigned cordCount = 0;

    cordCount = 0;
    for (unsigned k = 0; k < length(cordSet); k++)
    {
        //unsigned recordCount = 0;
        if (!empty(cordSet[k]))
        {
            for (unsigned j = 1; j < length(cordSet[k]); j++)
            {
                if (_DefaultHit.isBlockEnd(cordSet[k][j-1]) )//&& ++recordCount < 10)
                {
                    of <<"@S1_"<< k+1 << " " << length(reads()[k]) << " "
                    << _DefaultCord.getCordY(cordSet[k][j]) << " " << length(cordSet[k]) << " x " 
                    << _getSA_i1(_DefaultCord.getCordX(cordSet[k][j])) << " " << cordCount << " "
                    << _getSA_i2(_DefaultCord.getCordX(cordSet[k][j]))  << " " 
                    << "\n";   
                    cordCount = 0;
                }
                cordCount++;
            }
        }
    }
    std::cerr << ">Write results to disk          " << std::endl;
    std::cerr << "    End writing results. Time[s]" << sysTime() - time << std::endl; 
}
/**
 * print all cords with cordinates
 */
void Mapper::printCordsRaw2()
{
    std::cerr << ">>Write results to disk        \r";
    double time = sysTime();
    unsigned cordCount = 0;
    CharString first_line = "";
    cordCount = 0;
    uint64_t readCordEnd;
    uint64_t seqsCordEnd;
    char main_icon_strand = '+', icon_strand = '+';
    int fflag = 0;
    for (unsigned k = 0; k < length(cordSet); k++)
    {
        if (!empty(cordSet[k]))
        {
            for (unsigned j = 1; j < length(cordSet[k]); j++)
            {
                if (_DefaultHit.isBlockEnd(cordSet[k][j-1]))
                {
                    unsigned m = j; 
                    int main_strand_count = 0;
                    int block_len = 0;
                    ///>determine the main strand
                    while (!_DefaultHit.isBlockEnd(cordSet[k][m]))
                    {
                        if (_DefaultCord.getCordStrand(cordSet[k][m]))
                        {
                            main_strand_count++;
                        }
                        block_len++;
                        m++;
                    }
                    if (main_strand_count > (block_len >> 1))
                    {
                        main_icon_strand = '-';
                    }
                    else
                    {
                        main_icon_strand = '+';
                    }
                    ///>print the header
                    for (unsigned i = j; ; i++)
                    {
                        if (_DefaultHit.isBlockEnd(cordSet[k][i]) || i == length(cordSet[k]) - 1)
                        {
                            readCordEnd = _DefaultCord.getCordY(cordSet[k][i]) + window_width;
                            seqsCordEnd = _getSA_i2(_DefaultCord.getCordX(cordSet[k][i])) + window_width;
                            break;
                        }
                    }
                    of  << first_line;
                    of  << record.id1[k] << " " 
                        //<< length(record.seqs1[k]) << " "
                        << rlens[k] << " "
                        << _DefaultCord.getCordY(cordSet[k][j]) << " " 
                        << std::min(readCordEnd, (uint64_t)rlens[k]) << " " 
                        << main_icon_strand<< " "
                        << record.id2[_getSA_i1(_DefaultCord.getCordX(cordSet[k][j]))] << " " 
                        << length(record.seq2[_getSA_i1(_DefaultCord.getCordX(cordSet[k][j]))]) << " "
                        << _getSA_i2(_DefaultCord.getCordX(cordSet[k][j]))  << " " 
                        << seqsCordEnd << "\n";
                    first_line = "\n";
                    cordCount = 0;
                    fflag = 1;
                }
                ///>print the coordinates
                icon_strand = (_DefaultCord.getCordStrand(cordSet[k][j]))?'-':'+';
                CharString mark = "| ";
                if (icon_strand != main_icon_strand)
                    mark = "********** ";
                int64_t d = 0;//_DefaultCord.getCordY(cordSet[k][1]);
                int64_t d2 = 0;
                if (!fflag)
                {
                    d = (int64_t)_DefaultCord.getCordX(cordSet[k][j]) - (int64_t)_DefaultCord.getCordX(cordSet[k][j - 1]);
                    d2 = (int64_t)_DefaultCord.getCordY(cordSet[k][j]) - (int64_t)_DefaultCord.getCordY(cordSet[k][j - 1]);
                    
                }
                
                of  << mark  << _DefaultCord.getCordY(cordSet[k][j]) << " " 
                    << _getSA_i2(_DefaultCord.getCordX(cordSet[k][j])) << " " << d2 << " " << d << " " << j << " \n";
                cordCount++;
                fflag = 0;
            }
        }
    }
    close(of);
    std::cerr << "--Write results to disk       100% Elapsed Time[s] " << sysTime() - time << std::endl;
}

int print_clip_gff_(StringSet<String<uint64_t> > & clips, 
              StringSet<CharString> & readsId, 
              std::ofstream & of, 
              std::string & outputPrefix)
{
    std::string file_path = outputPrefix + ".gff";
    std::cerr << "[]::filepath " << file_path << "\n";
    of.open(toCString(file_path));
    for (unsigned i = 0; i < length(clips); i++)
    {
        if (!empty(clips[i]))
        {
            of << i << " " << readsId[i] << " ";
            for (unsigned j = 0; j < length(clips[i]); j++)
            {
                uint64_t cord_x = _DefaultCord.getCordX(clips[i][j]);
                uint64_t cord_y = _DefaultCord.getCordY(clips[i][j]);
                of << cord_x << " ";   
            }
            of << '\n';
        }
    }
    of.close();
    return 0;
}

int print_clip_gvf_(StringSet<String<uint64_t> > & clips, 
              StringSet<CharString> & readsId, 
              StringSet<CharString> & genomesId,
              std::ofstream & of, 
              std::string outputPrefix)
{
    std::string file_path = outputPrefix + ".gvf";
    //std::cerr << "[]::filepath " << file_path << "\n";
    of.open(toCString(file_path));
    std::string source = ".";
    std::string type = ".";
    for (unsigned i = 0; i < length(clips); i++)
    {
        if (!empty(clips[i]))
        {
            for (unsigned j = 0; j < length(clips[i]); j++)
            {
                uint64_t cord_x = _getSA_i2(_DefaultCord.getCordX(clips[i][j]));
                //uint64_t cord_y = _DefaultCord.getCordY(clips[i][j]);
                CharString genomeId = genomesId[_getSA_i1(_DefaultCord.getCordX(clips[i][j]))];
                if ((j >> 1) << 1 == j)
                {
                    of  << genomeId << " " << source << " " << type << " " << cord_x << " ";   
                    if (j == length(clips[i]) - 1)
                    {
                        of << " . readId=" << readsId[i] << ";" << i << "\n";
                    }
                }
                else
                {
                    of << cord_x << " readId=" << readsId[i] << ";" << i << "\n";
                }
            }
        }
    }
    of.close();
    return 0;
}

int print_clip_gff(Mapper & mapper)
{
    print_clip_gff_(mapper.getClips(), mapper.readsId(), mapper.getOf(), mapper.getOutputPrefix());
    return 0;
}

int print_clip_gvf(Mapper & mapper)
{
    print_clip_gvf_(mapper.getClips(), mapper.readsId(), mapper.genomesId(), mapper.getOf(), mapper.getOutputPrefix());
    return 0;
}

int rawMap_dst2_MF(LIndex & index,
                   StringSet<String<short> > & f2,
                   StringSet<String<Dna5> > & reads,
                   MapParm & mapParm,
                   StringSet<String<uint64_t> > & cords,
                   StringSet<String<uint64_t> > & clips,
                   StringSet<String<Dna5> > & seqs,
                   StringSet<String<BamAlignmentRecordLink> >& bam_records,
                   unsigned & threads,
                   int p1
                  )
{
  
    typedef String<Dna5> Seq;
    //double time=sysTime();
    float senThr = mapParm.senThr / window_width;
    float cordThr = mapParm.cordThr / window_width;
    MapParm complexParm = mapParm;
    complexParm.alpha = complexParm.alpha2;
    complexParm.listN = complexParm.listN2;
    String<uint64_t> gap_len;
    String<uint64_t> red_len;
    resize (gap_len, threads);
    resize (red_len, threads);
    for (int i = 0; i < length(gap_len); i++)
{
    gap_len[i]  = 0;
    red_len[i] = 0;
}
    //double time2 = sysTime();
int su = 0;
int64_t len = 0;
#pragma omp parallel
{
    unsigned size2 = length(reads) / threads;
    unsigned ChunkSize = size2;
    Seq comStr;
    //Anchors anchors(Const_::_LLTMax, AnchorBase::size);
    Anchors anchors;
    String<uint64_t> crhit;
    StringSet<String<uint64_t> >  cordsTmp;
    StringSet< String<short> > f1;
    StringSet<String<uint64_t> > clipsTmp;
    StringSet<String<BamAlignmentRecordLink> > bam_records_tmp;
    unsigned thd_id =  omp_get_thread_num();
    if (thd_id < length(reads) - size2 * threads)
    {
        ChunkSize = size2 + 1;
    }
    resize(cordsTmp, ChunkSize);
    resize(clipsTmp, ChunkSize);
    resize(bam_records_tmp, ChunkSize);
    resize(f1, 2);
    unsigned c = 0;
    
    String<uint64_t> g_hs;
    String<uint64_t> g_anchor;
    String<uint64_t> bands;
    resize (g_hs, 1ULL << 20);
    resize (g_anchor, 1ULL<<20);
    #pragma omp for
    for (unsigned j = 0; j < length(reads); j++)
    {
        if (length(reads[j]) >= mapParm.minReadLen)
        {
            red_len[thd_id] += length(reads[j]);
            std::cout << "[]::rawmap::j " << j <<"\n";
            float cordLenThr = length(reads[j]) * cordThr;
            _compltRvseStr(reads[j], comStr);
            createFeatures(begin(reads[j]), end(reads[j]), f1[0]);
            createFeatures(begin(comStr), end(comStr), f1[1]);
            anchors.init(1);
            clear(crhit);
            mnMapReadList(index, reads[j], anchors, mapParm, crhit);
            path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            if (getCordsMaxLen(cordsTmp[c]) < length(reads[j]) * senThr)
            {
                clear(cordsTmp[c]);
                anchors.init(1);
                clear(crhit);
                mnMapReadList(index, reads[j], anchors, complexParm, crhit);
                path_dst(begin(crhit), end(crhit), f1, f2, cordsTmp[c], cordLenThr);
            }   
            gap_len[thd_id] += mapGaps(seqs, reads[j], comStr, cordsTmp[c], g_hs, g_anchor, clipsTmp[c], f1, f2, p1, 192);
            //align_cords(seqs, reads[j], comStr, cordsTmp[c], bam_records_tmp[c], p1);
        }   
        c += 1;
    } 
    #pragma omp for ordered
    for (unsigned j = 0; j < threads; j++)
        #pragma omp ordered
        {
            append(cords, cordsTmp);
            append(clips, clipsTmp);
            append(bam_records, bam_records_tmp);
        }
    
}
    //std::cerr << "    End raw mapping. Time[s]: " << sysTime() - time << std::flush << std::endl;
    return 0;
}

/*
 *[]::map
 */
int map(Mapper & mapper, int p1)
{
    //printStatus();
    StringSet<String<short> > f2;
    mapper.createIndex(false); // true: destruct genomes string to reduce memory footprint
    createFeatures(mapper.genomes(), f2, mapper.thread());
    SeqFileIn rFile(toCString(mapper.readPath()));
    unsigned k = 1, j = 0;
    unsigned blockSize = 50000;
    StringSet<String<char> > dotstatus;
    resize(dotstatus, 3);
    dotstatus[0] = ".   ";
    dotstatus[1] = "..  ";
    dotstatus[2] = "... ";
    while (!atEnd(rFile))
    {
        double time1 = sysTime();
        clear (mapper.reads());
        std::cerr <<  ">>Map::file_I/O  block " << k << dotstatus[j++ % length(dotstatus)] << "\r";
        readRecords_block(mapper.readsId(), mapper.reads(), mapper.readLens(), rFile, blockSize);
        std::cerr << "                                    \r";
        std::cerr <<  ">>Map::mapping  block "<< k << " Size " << length(mapper.reads()) << " " << dotstatus[j++ % length(dotstatus)] << "\r";
        time1 = sysTime() - time1;
        double time2 = sysTime();
        rawMap_dst2_MF(mapper.index(), 
                       f2, 
                       mapper.reads(), 
                       mapper.mapParm(), 
                       mapper.cords(), 
                       mapper.getClips(),
                       mapper.genomes(),
                       mapper.getBamRecords(),
                       mapper.thread(), 
                       p1);
        time2 = sysTime() - time2;
        std::cerr <<  "--Map::file_I/O+Map block "<< k << " Size " << length(mapper.reads()) << " Elapsed Time[s]: file_I/O " << time1 << " map "<< time2 << "\n";
        k++;
    }
    mapper.index().clear(); 
    mapper.printCordsRaw2();
    print_align_sam(mapper);
    print_clip_gvf(mapper);
    return 0;
}

