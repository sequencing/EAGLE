EAGLE - Enhanced Artificial Genome Engine
=========================================

The Enhanced Artificial Genome Engine (EAGLE) software is designed to simulate 
the behaviour of Illumina's Next Generation Sequencing instruments, in order to 
facilitate the development and testing of downstream applications.


Dependencies
============

Required packages (based on a fresh Ubuntu install):

    CMake >= 2.8.2
    C++ compiler
    Boost development libraries >= 1.56
    LibXml
    Samtools
    dc
    bzip2
    which (in centos docker)
    time (in docker)

    All together:
        apt-get install cmake g++ libboost-all-dev libxml2-dev libxml-simple-perl samtools dc bzip2 time
      / yum install cmake gcc-c++ boost-devel libxml2-devel perl-libxml-perl samtools bc bzip2 which time


Installation Instructions
=========================

<EAGLE_source_dir>/src/configure

Note: We noticed that cmake doesn't always get detected. Adding "--with-cmake=cmake" helps. (I know...)


Pre-built Docker image
======================

The following public Docker image is available: `ljanin/eagle:2.0.4`

`sudo docker run -it --rm ljanin/eagle:2.0.4`

PS: This is not the latest version. TODO.


Documentation
=============

 - http://htmlpreview.github.io/?https://github.com/sequencing/EAGLE/blob/master/doc/html/index.html
 - http://htmlpreview.github.io/?https://github.com/sequencing/EAGLE/blob/master/doc/html/index2.html

