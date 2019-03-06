#ifndef MISC_H__
#define MISC_H__
#include <tuple>
#include <iostream>

struct Sub1ReadNum
{
    int total_num ;
    int sub1_num ;
    int used_num ;

    friend std::ostream& operator<<(std::ostream& out, const Sub1ReadNum & me)
    {
        out<<me.total_num<<"\t"<<me.sub1_num<<"\t"<<me.used_num;
        return out;
    }
    bool operator < (const Sub1ReadNum & other ) const
    {
        return std::make_tuple( total_num, sub1_num , used_num  ) 
            <  std::make_tuple( other.total_num , other.sub1_num ,other.used_num ); 
    }
};

struct Sub1_3ReadNum
{
    int total_num ;
    int sub1_num ;
    int sub3_num ;
    int sub1_3_num ;
    int used_num ;
    friend std::ostream& operator<<(std::ostream& out, const Sub1_3ReadNum & me)
    {
        out<<me.total_num<<"\t"<<me.sub1_num <<"\t"<<me.sub3_num<<"\t"<<me.sub1_3_num <<"\t"<<me.used_num;
        return out;
    }

    bool operator < (const Sub1_3ReadNum & other ) const
    {
        return std::make_tuple( total_num, sub1_num , sub3_num , sub1_3_num , used_num  ) 
            <  std::make_tuple( other.total_num , other.sub1_num ,
                    other.sub3_num , other.sub1_3_num ,
                    other.used_num ); 
    }
};

#endif // MISC_H__
