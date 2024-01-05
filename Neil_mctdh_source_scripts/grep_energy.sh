#!/bin/bash
clear
grep --color=always -E -o "(.*E-tot.*eV,|.*Time.*fs)" ./**/output | grep -v 'state'

