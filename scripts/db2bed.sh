#!/bin/bash

echo "select chr, start, end, rid from region" | mysql -u easih_ro -h mgsrv01 depths_exome_5bp | egrep -v start | sort -V