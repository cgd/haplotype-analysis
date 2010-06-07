#!/bin/bash

true_script_location=`readlink -fn $0`
#java -agentlib:jdwp=transport=dt_socket,address=8000,server=y,suspend=y -enableassertions -Xmx1g -jar `dirname $true_script_location`/dist/haplotype-analysis-1.0.jar $@
java -enableassertions -Xmx1g -jar `dirname $true_script_location`/dist/haplotype-analysis-1.0.jar $@
