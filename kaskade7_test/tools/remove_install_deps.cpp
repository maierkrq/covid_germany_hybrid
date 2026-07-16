/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2016-2018 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstring>
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>


int main(int argc, char *argv[])
{
  if (argc!=3)
  {
    std::cout << "usage: remove_install_deps Makefile <Kaskade 7 directory>" << std::endl;
    exit(1);
  }

  char tempFileNameTemplate[]="/tmp/remove_usr_dep_out_XXXXXX";
  char tempFileName[80];
  strcpy(tempFileName,tempFileNameTemplate);
  mkstemp(tempFileName);
  std::string buffer, compressedBuffer;
  const std::string beginDeps = "#DONOTDELETE";
  std::streampos outFileBegin;
  bool  returnOut=false;
  int pos, colonFound, removeCount=0;

  // Open input Makefile
  std::ifstream makefile(argv[1]);
  if ( makefile.rdstate() & std::ifstream::failbit )
  {
    std::cout << argv[0] << ": couldn't open file " << argv[1] << "." << std::endl;
    exit (10);
  };
  
  // Open output file and copy lines from in to out until the "# DO NOT DELETE" marker is found
  std::ofstream outStream(tempFileName);
  outFileBegin = outStream.tellp();
  while ( ! makefile.eof() && compressedBuffer != beginDeps )
  {
    std::getline(makefile,buffer);
    if ( makefile.eof() )
      break;

    // If we encounter a comment line (starting with #), extract this and
    // remove any blanks or tabs.
    if ( buffer.find('#') != std::string::npos )
    {
      compressedBuffer.clear();
      for (int i=0; i<buffer.length(); i++)
      {
        const char ch = buffer[i];
        if ( (ch!=' ') && (ch!='\t') ) 
         compressedBuffer.append(1,ch);
      };
    };
    outStream << buffer << std::endl;
  };
  
  if ( compressedBuffer != beginDeps )
  {
    std::cout << argv[0] << ": didn't change " << argv[1] << "." << std::endl;
    exit(0);
  };
  
  // find last occurrence of "Kaskade7.x" string in current path
  std::string kaskade7Dir = argv[2];
  int lengthKaskDir = kaskade7Dir.length();

  while ( true )
  {
    std::string dependency;
    makefile >> dependency;
    
    if ( makefile.eof() )
      break;
    colonFound = dependency.find(':');
    if ( colonFound != std::string::npos )
      outStream << std::endl;
    switch ( dependency[0] )
    {
      case ':':
      {
        outStream << dependency;
        returnOut = true;
        break;
      };
      case '\\':            // this is the line continuation symbol
      {
        if ( returnOut )
          outStream << " \\" << std::endl;
        returnOut = false;
        break;
      };
      case '/':
      {
        if (dependency.substr(0,lengthKaskDir) == kaskade7Dir)
        {
          outStream << " $(KASKADE7)" << dependency.substr(lengthKaskDir);
          returnOut = true;
        }
        else
          removeCount++;
        break;
      };
      default:
      {
        if (colonFound == std::string::npos)
          outStream << " ";
        outStream << dependency;
        returnOut = true;
        break;
      }
    }
  };
  outStream.flush();
  outStream.close();
  
  auto cmd = std::string("cat ") + tempFileName + " > " + argv[1];  
  system(cmd.c_str());
  
  std::cout << argv[0] << ": done with " << argv[1] << "." << std::endl;
  unlink(tempFileName);
}
