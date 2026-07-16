/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/en/numerik/software/kaskade-7.html               */
/*                                                                           */
/*  Copyright (C) 2015-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cctype>
#include <sstream>

#include "io/readHIP.hh"
#include "utilities/detailed_exception.hh"

boost::multi_array<double,2> readHIPData(std::string const& name){
    using namespace Kaskade;

    std::string read = name+".HIP";

    // read motion data file

    boost::multi_array<double,2> motionData(boost::extents[101][5]);
    //*motionData(boost::extents[101][5]);

    //open file with motion data
    std::ifstream in(read.c_str());
    if (!in) // check for successful opening
        throw Kaskade::FileIOException("Cannot open for reading.",read,__FILE__,__LINE__);

    std::string line;

    //Delete the first ten lines. Info in lines six and seven might be interesting later
    for(int j=0;j<10;j++) {
        std::getline(in,line);
        line.erase(line.begin(),line.end());
    }

    //For 0% to 100%: read out forces and time and write them into boost::multi_array motionData
    for(int i=0;i<=100;i++){
        std::getline(in,line);
        std::istringstream is(line);

        //Read out the percentage of the cycle (0-100, discarded to nCycle), the (negative) forces in x,y,z direction, the resulting force (all in Newton) and the time (in seconds)
        int nCycle;
        is >> nCycle >> motionData[i][0] >> motionData[i][1] >> motionData[i][2] >> motionData[i][3] >> motionData[i][4];
    }

    return motionData;
}

boost::multi_array<double,2> readHIPData2(std::string const& name){
    using namespace Kaskade;

    std::string read = name+"_2016.HIP";

    // read motion data file

    boost::multi_array<double,2> motionData(boost::extents[601][8]);
    //*motionData(boost::extents[101][5]);

    //open file with motion data
    std::ifstream in(read.c_str());
    if (!in) // check for successful opening
        throw Kaskade::FileIOException("Cannot open for reading.",read,__FILE__,__LINE__);

    std::string line;

    //Delete the first five lines. 
    for(int j=0;j<5;j++) {
        std::getline(in,line);
        line.erase(line.begin(),line.end());
    }

    //For 0% to 100%: read out forces and time and write them into boost::multi_array motionData
    for(int i=0;i<=600;i++){
        std::getline(in,line);
        std::istringstream is(line);

        //Read out the forces in x,y,z direction resulting Force (in Newton), moments around the x,y,z axis + resulting moment (in Newton*meter) and the percentage of the cycle (0-100, in .165 steps),
        int nCycle;
        is >> motionData[i][0] >> motionData[i][1] >> motionData[i][2] >> motionData[i][3] >> motionData[i][4] >> motionData[i][5] >> motionData[i][6] >> motionData[i][7] >> nCycle;
    }

    return motionData;
}
