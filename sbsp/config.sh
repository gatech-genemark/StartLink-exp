
# get path of config script
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"


# Paths to directories
# base: root directory for this project (where config.sh is located)
# bin: directory containing generated binary files (from python scripts)
# bin_external: directory containing external binaries (e.g. GeneMarkS-2, prodigal, etc...)
# data: data directory (in it, we have GENOME/sequence.fasta and GENOME/ncbi.gff)
# db: database index files (and possibly databases themselves)
# metadata: assembly summaries, taxonomy tree information, etc...
# lists: contains genome list files (e.g. names of verified genomes, etc...)
# config: configuration files (e.g. algorithm parameters, parallelization, etc...)
# tmp: empty directory to do stuff in
# code: code :)

export base=$DIR
export bin=$base/bin
export bin_external=$base/bin_external

# data
export data_root=$base/data
export data=${data_root}/all
export data_verified=${data_root}/verified
export metadata=$base/metadata

# other
export db=$base/db
export lists=$base/lists
export config=$base/config
export runs=$base/runs
export tmp=$base/tmp               # it is suggested that runs be performed in tmp directory (under some subdirectory)
export code=$base/code

export libbash=$code/bash/lib
export libpython=$code/python/lib

export driverbash=$code/bash/driver
export driverpython=$code/python/driver

# DEPRECATED
export matrices=$base/matrices
export exp=$base/experiments
