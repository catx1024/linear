#ifndef LINEAR_HEADER_ALIGN_UTIL_H
#define LINEAR_HEADER_ALIGN_UTIL_H

#include <seqan/bam_io.h>

using namespace seqan;

extern uint16_t bam_flag_rvcmp;
extern uint16_t bam_flag_rvcmp_nxt;
extern uint16_t bam_flag_suppl;
struct SAZTag
{
    String<unsigned> bam_records_i;
    String<bool> bam_records_is_chimeric;
    String<bool> bam_records_is_block_end;
};
class BamAlignmentRecordLink : public BamAlignmentRecord 
{ 
public:
    //the heads refers to the first record of each line 
    //line refers to each line in.sam
    //line is compsed of several records in the String<bamAlig...>
    int next_id; //next records id
    String<CigarElement<> > saz_cigar;
    String<unsigned> heads_table; //table pointing to first record of each line in .sam 
    CharString genome_id;
    int nm_i; //nm:i tag in .sam, the heads holds the whole(sum of) value of the line

    BamAlignmentRecordLink();
    void addNext(int id);
    int isEnd() const;
    int next() const;
};
/*
 * The set of functions to manipulate String<BamAlignmentRecordLink> 
 * The functions are supposed to wrapper with String<BamAlignmentRecordLink> to
   make s new class. However the String<BamAlignmentRecordLink> is used in the 
   interface of many functions. To avoid the modification of interface this struct 
   is declared.
 */
struct BamLinkStringOperator
{
    int updateHeadsTable(String<BamAlignmentRecordLink> & bam_records);
    int getHeadNum(String<BamAlignmentRecordLink> & bam_records);
    int getHead(String<BamAlignmentRecordLink> & bam_records, int i);
    int writeBamRecordLinkCigar(
            std::ofstream target,
            String<BamAlignmentRecordLink> & bam_records,
            int it);
    int createSAZTagCigarOneChimeric(
            String<BamAlignmentRecordLink> & bam_records,
            String<CigarElement<> > & cigar,
            int it,
            bool f_force);
    int createSAZTagOneChimeric(
            String<BamAlignmentRecordLink> & bam_records,
            CharString & saz_tag, 
            int it);
    int createSAZTagOneLine(
            String<BamAlignmentRecordLink> & bam_records,
            int it);
};

void align2cigar(String<CigarElement< > > &cigar,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps1,
                 Row<Align<String<Dna5>,ArrayGaps> >::Type &gaps2
                );
int insertCigar(String<CigarElement< > > &cigar1, 
                int pos,
                String<CigarElement< > > &cigar2
         );

/*
 * insert cigar to the original cigar 
 */
int insertBamRecordCigar (BamAlignmentRecord & bam_record,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                    Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                    int pos = -1
                   );

/*
 * Cigar of row1 and row2 are inserted to the cigar from the bam_record
 * beginPos are always updated by g_beginPos 
 * soft/Hard clip cigar are updated only if insert at the front(pos == 0)
 */
int  insertBamRecord (BamAlignmentRecord & bam_record,
                      Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                      Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                      int g_id,
                      int g_beginPos,
                      int r_beginPos,
                      int pos = -1,
                      int f_soft = 1
                      );

/**
 * @g_beignPos, @r_beginPos
 * 1-based leftmost exact start coordinates 
 */
int insertNewBamRecord (String<BamAlignmentRecordLink> & bam_records,
                        int g_id, 
                        int g_beginPos,
                        int r_beginPos,
                        int strand,
                        int insert_pos = -1,
                        int f_soft = 1, /*hard and soft clip flag*/
                        uint16_t flag = 0
                        );

int insertNewBamRecord (String<BamAlignmentRecordLink> & bam_records,
                        Row<Align<String<Dna5>, ArrayGaps> >::Type & row1,
                        Row<Align<String<Dna5>, ArrayGaps> >::Type & row2, 
                        int g_id,
                        int g_beginPos,
                        int r_beginPos,
                        int strand,
                        int insert_pos = -1,
                        int f_soft = 1 
                        );


#endif
