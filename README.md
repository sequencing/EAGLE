EAGLE - Enhanced Artificial Genome Engine
=========================================

The Enhanced Artificial Genome Engine (EAGLE) software is designed to simulate 
the behaviour of the actual instrument (sequencer), hence facilitating the 
development and testing of downstream applications.


Dependencies
============

Required packages (based on a fresh Ubuntu install):
    CMake >= 2.8.2
        yum/apt-get install cmake

    C++ compiler
        yum/apt-get install g++

    Boost development libraries >= 1.47
        yum/apt-get install libboost-all-dev

    LibBz2 
        yum/apt-get install libbz2-dev

    LibXml
        yum/apt-get install libxml2-dev libxml-simple-perl

    Samtools
        yum/apt-get install samtools

    dc
        yum/apt-get install dc


Installation Instructions
=========================

<EAGLE_source_dir>/src/configure

Note: We noticed that cmake doesn't always get detected. Adding "--with-cmake=cmake" helps. (I know...)


Documentation
=============

Goto doc/html/index.html

