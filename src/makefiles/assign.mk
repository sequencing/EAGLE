################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## description assignment of N variables.
##
## This Makefile recursively assigns N values into N variables. If there are 
## any ':' in any value_i, this will be treated as a list.
##
## definition
##   assign.keys:=var1 var2 ... varN
##   assign.values:=value1 value2 ... valueN
##   include assign.mk
##
## example:
##   assign.keys:=v1 v2 v3
##   assign.values:=a b:c:d e
##   include assign.mk
## is equivalent to:
##   v1:=a
##   v2:=b c d
##   v3:=e
##
## author Mauricio Varea
##
################################################################################
ifneq ($(words $(assign.keys)),$(words $(assign.values)))
  $(error $(words $(assign.keys)) keys vs. $(words $(assign.values)) values :\
          Multiple assignments only work if the lists are of equal lengths)
endif
ifneq (,$(assign.keys))
  $(firstword $(assign.keys)):=$(subst :, ,$(firstword $(assign.values)))
  ifneq (1,$(words $(assign.keys)))
    assign.keys:=$(wordlist 2,$(words $(assign.keys)),$(assign.keys))
    assign.values:=$(wordlist 2,$(words $(assign.values)),$(assign.values))
  else
    assign.keys:=
    assign.values:=
  endif
  include $(MAKEFILES_DIR)/assign.mk
endif
