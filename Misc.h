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

struct SubReadsLog
{
    std::string type ;
    // PE or PE_BARCODE or BASIC or NULL
    int basic_num ;
    int pe_num ;
    int pe_with_extra_barcode_num ;
    int barcode_num ;
    int pe_barcode_num ;
    int used_num ;

    void Init() 
    {
        type = "" ;
        basic_num = 0 ;
        pe_with_extra_barcode_num = 0;
        pe_num = 0 ;
        barcode_num = 0 ;
        pe_barcode_num = 0 ;
        used_num = 0 ;
    }

    friend std::ostream& operator<<(std::ostream& out, const SubReadsLog & me)
    {
        out<<me.type<<'\t'<<me.basic_num
            <<'\t'<<me.pe_with_extra_barcode_num
            <<'\t'<<me.pe_num;
        if( me.type != "PE" && me.type != "PE&Barcode" )
            out<<'\t'<<me.barcode_num<<'\t'<<me.pe_barcode_num;
        else
            out<<'\t'<<"NA"<<'\t'<<"NA";
        out<<'\t'<<me.used_num;
        return out;
    }

    bool operator < (const SubReadsLog & other ) const
    {
        return std::make_tuple( type ,
                basic_num ,
                pe_with_extra_barcode_num ,
                pe_num ,
                barcode_num , 
                pe_barcode_num ,
                used_num  )
            <  std::make_tuple( other.type,
                    other.basic_num ,
                    other.pe_with_extra_barcode_num ,
                    other.pe_num ,
                    other.barcode_num ,
                    other.pe_barcode_num , 
                    other.used_num  ) ;
    }
};

#endif // MISC_H__
