set(CTEST_PROJECT_NAME "EAGLE")

set(MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes")

set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")
set(CTEST_DROP_METHOD "http")
#set(CTEST_DROP_SITE "ukch-cdash01.chuk.illumina.com")
set(CTEST_DROP_SITE "10.46.146.88") # <- Come's machine
#set (CTEST_DROP_LOCATION "cgi-bin/HTTPUploadDartFile.cgi")
set (CTEST_DROP_LOCATION "/CDash/submit.php?project=EAGLE")
set(CTEST_DROP_SITE_CDASH TRUE)

#set (CTEST_TRIGGER_SITE 
#       "http://${DROP_SITE}/cgi-bin/Submit-vtk-TestingResults.pl")
