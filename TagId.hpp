#ifndef TAGID_HPP_
#define TAGID_HPP_
#include <map>
#include <string>

struct TagId
{
    private:
        std::map<std::string , unsigned int > data ;
        unsigned int next = 0 ;
    public:
        unsigned int GetId(const std::string & id ,bool new_id = true  )
        {
            if( id == "0_0_0" 
                    ||id == "0_0"
                    || id == "0" )
                return 0;
            auto itr = data.find(id) ;

            if( itr == data.end() )
            {
                if( new_id )
                {
                    next ++ ;
                    data[id]=next ;
                    return next ;
                }
                else
                {
                    return 0 ;
                }
            }
            return itr->second;
        }
};

#endif
