#ifndef BinaryFile_h
#define BinaryFile_h

#include <string>
#include <iostream>
#include <fcntl.h>
#include <unistd.h>
#ifdef LINUX
#include <stdlib.h>
#endif

template <typename AtAddressableT>
class BinaryFile{
  public:
    typedef typename AtAddressableT::value_type value_type;

    static int read(AtAddressableT *buffer, std::string filePath) {
        int inputFile = open(filePath.c_str(), O_RDONLY);
    
        if(-1 == inputFile) {
          std::cerr << "Could not open file: " << filePath << std::endl;
          exit(1);
        }

        int bytesRead =
            ::read(inputFile, &(buffer->at(0)),
              sizeof(value_type) * buffer->size());

        close(inputFile);

        return bytesRead;
    }
    
    static int write(AtAddressableT *buffer, std::string filePath) {
        // create the file with user-only permissions
        int outputFile = open(filePath.c_str(), O_WRONLY | O_CREAT, 0600);
    
        if(-1 == outputFile) {
          std::cerr << "Could not open file: " << filePath << std::endl;
          exit(1);
        }

        int bytesWritten = 
            ::write(outputFile, &(buffer->at(0)),
              sizeof(value_type) * buffer->size());
    
        close(outputFile);

        return bytesWritten;
    }

};

#endif
