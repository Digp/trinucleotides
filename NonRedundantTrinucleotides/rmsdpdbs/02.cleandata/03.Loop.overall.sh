#! /bin/bash

for i in $(ls *txt); do grep -v '      '$(awk '{print $1}' $i | sort | uniq -c | sort -rn | head -n 1 | awk '{print $2}') $i | grep -v total | awk '{print $2}'> $i.torm; done

