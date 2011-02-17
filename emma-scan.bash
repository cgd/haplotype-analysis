#!/bin/bash

# TODO: this file is currently OS X specific. We need it to be agnostic
# TODO: this requires libemma.dylib. figure out how this should be baked into the jar if possible

# exit on error and don't allow the use of unset variables
set -o errexit
#set -o nounset

export JAVA_HOME=/System/Library/Frameworks/JavaVM.framework/Versions/1.5/Home/
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:`dirname $0`

java -enableassertions -Xmx1g -jar `dirname $0`/dist/haplotype-analysis-1.0.jar $@

