#! /bin/bash
#  Blows away the input database and recreates it.
#usage: resetDB.sh <dbname>  
# NOTE: MUST BE RUN AS ROOT. YOU HAVE TO REALLY WANT TO DO THIS
#
function header () {
    echo "#####################################"
    echo "${1}"
    echo "#####################################"
}

MYDB=$*
PW=`python -c "print '\x50\x39\x6c\x73\x34\x52\x2a\x40'"`
alias MYSQL='mysql --host=10.0.1.20 --user=obs --password=${PW}'
shopt -s expand_aliases

CREATEDB="CREATE DATABASE $*;
USE $*;
SET PASSWORD FOR 'obs'@'localhost' = PASSWORD('${PW}');
GRANT ALL ON *.* TO 'obs'@'localhost';"

DELETEDB="DROP DATABASE $*;"

MYSQL $MYDB -e "show tables;"
dbbackup="${*}-$(date +'%Y-%m-%dT%H-%M').dbbackup"
header "saving table backup to $dbbackup"
mysqldump $MYDB --password=${PW} > ${dbbackup}
header "killing db"
MYSQL $MYDB -e "${DELETEDB}"
header "creating new db "$*
mysql --password=${PW} -e "${CREATEDB}"

python -c "from ddr_compress.dbi import DataBaseInterface; dbi = DataBaseInterface(); dbi.createdb()"

echo "testing db"
python -c "from ddr_compress.dbi import DataBaseInterface; dbi = DataBaseInterface(); dbi.test_db()"

