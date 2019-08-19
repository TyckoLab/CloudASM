#!/bin/bash

# May need to change the delimiter to space, not tab.
awk 'BEGIN { FS=OFS="\t" } \
              {if (FNR!=1) print $0}' ${CONTEXT_FILES} > ${CONTEXT_MERGED}