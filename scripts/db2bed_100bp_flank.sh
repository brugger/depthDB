#!/bin/bash

echo "select chr, start, end, rid from region" | mysql -u easih_ro -h mgsrv01 depths_exome_5bp | perl -ane '$F[1] -= 100; $F[2] += 100; print join("\t", @F) ."\n"' | egrep -v chr | sort -V