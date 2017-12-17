#!/usr/bin/env bash

LIMIT=500

COSMIC_DB="Cosmic.db"
OUT_CSV_FILE="CosmicTest.csv"
OUT_DB_FILE="CosmicTest.db"

[ -f ${OUT_CSV_FILE} ] && rm ${OUT_CSV_FILE}
[ -f ${OUT_DB_FILE} ] && rm ${OUT_DB_FILE}

sqlite3 Cosmic.db <<EOF
.echo off 
.headers on
.mode csv
.output ${OUT_CSV_FILE} 
SELECT * FROM Cosmic LIMIT ${LIMIT};
EOF

sqlite3 ${OUT_DB_FILE} <<EOF
.echo off
.mode csv
.import ${OUT_CSV_FILE} CosmicTest
CREATE INDEX GeneIndex ON CosmicTest("Gene name");
VACUUM;
EOF

