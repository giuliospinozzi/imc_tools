#!/bin/bash
source /etc/environment
source /etc/profile

echo "
  +--------------------------------------------------------+
  |                                                        |
  |                 IMC Segmentation Tool                  |
  |                                                        |
  +--------------------------------------------------------+
  |  Author:   Giulio Spinozzi, PhD                        |
  |  Date:     February 2024                               |
  |  Contact:  giulio.spinozzi@opbg.net                    |
  |  Version:  1.0                                         |
  +--------------------------------------------------------+

  REQUIRED VARS and relative ORDER POSITION -> REMEMBER NO SPACES!!!!!!!!!
	1. Root Path of Steinbock folders
	2. Path of panel.tsv (file for segmentation, with deepcell)

"

RUN_STARTED_AT=`date +"%Y-%m-%d %H:%M:%S"`;

##### ================================ ARGS ================================= #####
usage()
{
    echo "This app is the OPBG IMC Segmentation Tool for Steinbock."
    echo
    echo "Usage: $0 [-r ROOTDIR] [-p panel.tsv]"
    echo
    echo "  [-r ROOTDIR]    - Root Path of Steinbock folders (without last slash!)"
    echo "  [-p panel.tsv]  - Path of panel.tsv (file for segmentation, with deepcell)"
    echo
    exit
}

while getopts ":r:p:h" Option
    do
    case $Option in
        r ) ROOT=$OPTARG ;;
        p ) PANEL=$OPTARG ;;
        h ) usage ;;
        * ) echo "unrecognized argument. use '-h' for usage information."; exit -1 ;;
    esac
done
shift $(($OPTIND - 1))
#=================================================================================#


##### ========================== PRELIMINARY CHECKS ========================= #####
## Arguments
if [ -z "$ROOT" ]; then
    usage
fi

if [ -z "$PANEL" ]; then
    usage
fi

if [ ! -d ${ROOT} ]; then
    echo ""
    echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    echo " ${ROOT} NOT EXISTS!!!!! "
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    exit 0
fi

if [ ! -f ${PANEL} ]; then
    echo ""
    echo "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
    echo " ${PANEL} NOT EXISTS!!!!! "
    echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    exit 0
fi
#=================================================================================#


#---------------------------------------***---------------------------------------#
##### =============================== PROGRAM =============================== #####
echo "

---------------------------------------------------------------------------------
                    STARTING PROCESSING AT: $RUN_STARTED_AT
---------------------------------------------------------------------------------
    "

echo "<`date +'%Y-%m-%d %H:%M:%S'`> [OPBG] Preprocessing input variables"
## print input variables (check for log utils)
INPUTVARNUM=0
for INPUTVAR in "$@"; do
	let INPUTVARNUM++; 
	printf -v INPUTNUM '%02d' $INPUTVARNUM;
    echo "  => Input Variable: Order number = <${INPUTNUM}> ; Var Content = <${INPUTVAR}>";
done


echo "<`date +'%Y-%m-%d %H:%M:%S'`> [OPBG] Processing Data"
for DIR in ${ROOT}/*/ ; do
    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [OPBG] Processing Sample ${DIR}"
    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [OPBG] Processing Sample: Hot Pixel Removal"
    python3 /opt/scripts/imc_tools/filter_steinbock.py ${DIR} --type tophat --value 3
    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [OPBG] Processing Sample: Steinbock - preprocess external images"
    docker run -v $DIR:/data -u $(id -u):$(id -g) --network host -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY ghcr.io/bodenmillergroup/steinbock:0.16.1 preprocess external images
    cp ${PANEL} $DIR
    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [OPBG] Processing Sample: Steinbock - Segmentation"
    docker run -v $DIR:/data -u $(id -u):$(id -g) --network host -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY ghcr.io/bodenmillergroup/steinbock:0.16.1 segment deepcell --minmax
    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [OPBG] Processing Sample: Steinbock - Measures"
    docker run -v $DIR:/data -u $(id -u):$(id -g) --network host -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY ghcr.io/bodenmillergroup/steinbock:0.16.1 measure intensities
    docker run -v $DIR:/data -u $(id -u):$(id -g) --network host -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY ghcr.io/bodenmillergroup/steinbock:0.16.1 measure regionprops
    docker run -v $DIR:/data -u $(id -u):$(id -g) --network host -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY ghcr.io/bodenmillergroup/steinbock:0.16.1 measure neighbors --type centroids --dmax 20
    tar --use-compress-program="pigz -k " -cf ${DIR}external.tar.gz ${DIR}external --absolute-names
    rm -rf ${DIR}external
    echo "<`date +'%Y-%m-%d %H:%M:%S'`> [OPBG] --------------------------------------
    
    "
done

docker rm $(docker ps -a -f status=exited -q)

echo "

    ---------------------------------------------------------------------------------
                          ENDING PROCESSING AT: `date +'%Y-%m-%d %H:%M:%S'`
    ---------------------------------------------------------------------------------
        "