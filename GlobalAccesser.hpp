#ifndef GLOBALACCESSER_HPP__
#define GLOBALACCESSER_HPP__

class ReadAccessor ;
class PairInfo ;

struct GlobalAccesser
{
    static ReadAccessor * the_read_accessor;
    static PairInfo * the_pair_info ;
};

#endif //GLOBALACCESSER_HPP__
