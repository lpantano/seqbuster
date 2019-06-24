#nice options for bash scripts
set -o pipefail # trace ERR through pipes
set -o errtrace # trace ERR through 'time command' and other functions
set -o nounset ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'


echo "download miRBase files"
if [ ! -e hairpin.fa ] ; then
    wget -qN ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz -O hairpin.fa.gz
    wget -qN ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.str.gz -O miRNA.str.gz

    gunzip -f hairpin.fa.gz miRNA.str.gz
fi

if [ ! -e sim.21.hsa.fa ] ; then

    python miRNA.simulator.py -f hairpin.fa -m miRNA.str -n 10 -s hsa -p sim.21.hsa -e

fi

    SECONDS=0
    java -jar ../miraligner/miraligner-3.5/miraligner.jar -sub 1 -trim 5 -add  5 -s hsa -i sim.21.hsa.fa -db . -o sim.21.hsa.3.5 -pre
    java -jar ../miraligner/miraligner-3.4/miraligner.jar -sub 1 -trim 5 -add  5 -s hsa -i sim.21.hsa.fa -db . -o sim.21.hsa.3.4 -pre
    java -jar ../miraligner/miraligner-3.3/miraligner.jar  -sub 1 -trim 3 -add  3 -s hsa -i sim.21.hsa.fa -db . -o sim.21.hsa.3.3 -pre
    java -jar ../miraligner/miraligner-3.2/miraligner.jar  -sub 1 -trim 3 -add  3 -s hsa -i sim.21.hsa.fa -db . -o sim.21.hsa.3.2 -pre
    java -jar ../miraligner/miraligner-3.1/miraligner.jar  -sub 1 -trim 3 -add  3 -s hsa -i sim.21.hsa.fa -db . -o sim.21.hsa.3.1 -pre
    java -jar ../miraligner/miraligner-3.0/miraligner.jar -sub 1 -trim 3 -add  3 -s hsa -i sim.21.hsa.fa -db . -o sim.21.hsa.3.0 -pre
    duration=$SECONDS
    echo $(($duration / 60)) minutes elapsed

echo "now time to run stats.rmd"

#/Applications/RStudio.app/Contents/MacOS/pandoc/pandoc
Rscript -e 'library(rmarkdown);render("stats_isomirs.rmd")'
