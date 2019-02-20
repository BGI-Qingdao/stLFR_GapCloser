#ifndef STLFRREADHEADER_HPP_
#define STLFRREADHEADER_HPP_

#include <string>
#include "stringtool.hpp"
#include <cassert>

struct stLFRHeader
{
    enum ReadType
    {
        Unknow = 0 ,
        readName_barcodeStr =  1 ,
        readName_barcodeStr_index =  2 ,
        readName_barcodeStr_index_barcodeNum = 3 ,
    } type ;
    stLFRHeader() {
        Reset();
    }

    stLFRHeader(const stLFRHeader & a) {
        type = a.type ;
        readIndex  = a.readIndex ;
        barcode_num = a.barcode_num ;
        barcode_str = a.barcode_str ;
        readName = a.readName ;
    }

    stLFRHeader & operator = ( const stLFRHeader & a ) {
        if( &a != this )
        {
            type = a.type ;
            readIndex  = a.readIndex ;
            barcode_num = a.barcode_num ;
            barcode_str = a.barcode_str ;
            readName = a.readName ;
        }
        return *this ;
    }

    int readIndex ; //  1/2/3

    int barcode_num;

    std::string barcode_str ;

    std::string readName;

    void Init( const std::string & line ) 
    {
        Reset() ;

        auto items = BGIQD::STRING::split(line ,'@');

        assert(!items.empty());

        auto items_1 = BGIQD::STRING::split(items[0]);

        assert(!items_1.empty());

        auto items_2 = BGIQD::STRING::split(items_1[0],'#');

        assert(items_2.size() == 2);

        readName = items_2[0] ;

        auto items_3 = BGIQD::STRING::split(items_2[1] , '/');

        assert( !items_3.empty());

        barcode_str = items_3[0] ;

        type = readName_barcodeStr ;

        if( items_3.size() == 2 && BGIQD::STRING::IsNum(items_3[1]) )
        {
            readIndex = std::stoi(items_3[1]);
            type = readName_barcodeStr_index ;
        }

        if( items_1.size() >1 && BGIQD::STRING::IsNum(items_1[1]) )
        {
            barcode_num = std::stoi(items_1[1]);
            type = readName_barcodeStr_index_barcodeNum ;
        }
    }

    void Reset()
    {
        type = Unknow ;
        readName ="";
        barcode_str = "" ;
        readIndex = 0 ;
        barcode_num = 0 ;
    }
};

#endif
