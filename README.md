EAGLE - Enhanced Artificial Genome Engine
=========================================

The Enhanced Artificial Genome Engine (EAGLE) software is designed to simulate 
the behaviour of Illumina's Next Generation Sequencing instruments, in order to 
facilitate the development and testing of downstream applications.


Dependencies
============

Required packages (based on a fresh Ubuntu install):

    CMake >= 2.8.2
        yum/apt-get install cmake

    C++ compiler
        yum/apt-get install g++

    Boost development libraries >= 1.47
        yum/apt-get install libboost-all-dev

    LibXml
        yum/apt-get install libxml2-dev libxml-simple-perl

    Samtools
        yum/apt-get install samtools

    dc
        yum/apt-get install dc

    all together:
        yum/apt-get install cmake g++ libboost-all-dev libxml2-dev libxml-simple-perl samtools dc


Installation Instructions
=========================

<EAGLE_source_dir>/src/configure

Note: We noticed that cmake doesn't always get detected. Adding "--with-cmake=cmake" helps. (I know...)


Documentation
=============

Goto doc/html/index.html

