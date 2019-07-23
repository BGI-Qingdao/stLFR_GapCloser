#ifndef CONSENSUSLOG_HPP__
#define CONSENSUSLOG_HPP__

#include <string>
#include <vector>
#include <iostream>


namespace LOG{

    struct PosInfo
    {
        int cov ;
        char type ;
    };

    class ConsensusLog
    {
        public:
            static bool toggled;

            //std::string scaff_name ;
            int scaff_name ;
            int gap_id ;
            bool l2r;
            int consensus_id ;
            bool can_consensus ;
            int size_b ; //basic
            int size_p ; // PE
            int size_bc ; // barcode
            int size_bp ; // barcode & PE
            std::string used_id ;
            std::vector<PosInfo> details;
            bool is_consensus_done;
            bool is_gap_closer ;
            static char YN(bool f) { return f ? 'Y' : 'N' ;}

            void Init()
            {
                scaff_name = -1 ;
                gap_id = -1 ;
                l2r = true ;
                consensus_id = -1 ;
                can_consensus = -1 ;
                size_b = -1 ;
                size_p = -1 ;
                size_bc = -1 ;
                size_bp = -1 ;
                details.clear();
                used_id = "NULL";
                is_consensus_done = false ;
                is_gap_closer = false ;
            }
            void Print()
            {
                std::cerr<<"Consensus_Detail:"<<'\t'
                    <<scaff_name<<'\t'
                    <<gap_id<<'\t'
                    <<YN(l2r)<<'\t'
                    <<consensus_id<<'\t'
                    <<YN(can_consensus) <<'\t'
                    <<size_b<<'\t'
                    <<size_p<<'\t'
                    <<size_bc<<'\t'
                    <<size_bp<<'\t'
                    <<used_id<<'\t';
                if( details.empty() )
                    std::cerr<<'-';
                for(int i = 0 ; i<(int)details.size() ; i++ )
                {
                    std::cerr<<details[i].cov<<','
                        <<details[i].type;
                    if( i != (int)details.size() -1 )
                        std::cerr<<':';
                }
                std::cerr<<'\t'
                    <<YN(is_consensus_done)<<'\t'
                    <<YN(is_gap_closer)<<'\n';
            }
    };

} // namespace LOG
#endif //CONSENSUSLOG_HPP__
