
test -e ssshtest || curl -O https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
source ssshtest
#STOP_ON_FAIL=1
# Work inside of the repo folder
cd "$( dirname "${BASH_SOURCE[0]}" )"/../
INDIR=repo_utils/test_files
OD=test_results
COVERAGE_RCFILE=.coveragerc

# Reset test results
rm -rf $OD
mkdir -p $OD

tdb="coverage run --concurrency=multiprocessing -p -m tdb.__main__"
# ------------------------------------------------------------
#                                 test helpers
# ------------------------------------------------------------
fn_md5() {
    fn=$1
    # simple md5sum checking
    md5sum $fn | cut -f1 -d\  
}

tdb_check() {
    # check if parquet files are same
    dir_name=$1
    $tdb dbg_eq $INDIR/tdb/$dir_name $OD/$dir_name/
    assert_equal $? 0
    assert_exit_code 0
}

# ------------------------------------------------------------
#                                 entry 
# ------------------------------------------------------------
#run test_help $tdb
#if [ $test_help ]; then
    #assert_exit_code 0
    #assert_equal $(fn_md5 $STDERR_FILE) $(fn_md5 $ANSDIR/help.txt)
#fi

#run test_version $tdb version
#if [ $test_help ]; then
    #assert_exit_code 0
#fi

# ------------------------------------------------------------
#                               create
# ------------------------------------------------------------
run test_create1 $tdb create -o $OD/HG00438_chr14.tdb $INDIR/vcf/HG00438_chr14.vcf.gz
if [ $test_create1 ]; then
    tdb_check HG00438_chr14.tdb
fi

run test_create2 $tdb create -o $OD/HG00741_chr14.tdb $INDIR/vcf/HG00741_chr14.vcf.gz
if [ $test_create2 ]; then
    tdb_check HG00741_chr14.tdb
fi

run test_create3 $tdb create -o $OD/HG02630_chr14.tdb $INDIR/vcf/HG02630_chr14.vcf.gz
if [ $test_create3 ]; then
    tdb_check HG02630_chr14.tdb
fi

run test_create_badparam $tdb create -o repo_utils $INDIR/vcf/doesntexist
if [ $test_create_badparam ]; then
    assert_exit_code 1
fi

run test_create_mergedvcf $tdb create -o $OD/merged_singlevcf.tdb $INDIR/vcf/merged.vcf.gz
if [ $test_create_mergedvcf ] ; then
    tdb_check merged_singlevcf.tdb
fi
#
# ------------------------------------------------------------
#                                 merge
# ------------------------------------------------------------
#run test_merge1 $tdb merge -o $OD/merge1.tdb $INDIR/tdbs/some.sample.
run test_merge $tdb merge -o $OD/merge1.tdb $INDIR/tdb/HG00438_chr14.tdb/ $INDIR/tdb/HG00741_chr14.tdb/ $INDIR/tdb/HG02630_chr14.tdb/
if [ $test_merge ]; then
    tdb_check merge1.tdb
fi

run test_bad_merge $tdb merge -o $OD/merge1 $INDIR/tdb/HG00438_chr14.tdb/ $INDIR/tdb/HG00438_chr14.tdb/ $INDIR/HG00741_chr14
if [ $test_bad_merge ]; then
    assert_exit_code 1
fi

#if [ $test_merge ]; then
    #assert_equal $(fn_md5 $INDIR/
# ------------------------------------------------------------
#                                 query
# ------------------------------------------------------------

run test_q_allele_cnts $tdb query allele_cnts $INDIR/tdb/HG00438_chr14.tdb -o $OD/allele_cnts.txt
if [ $test_q_allele_cnts ]; then
    assert_equal $(fn_md5 $INDIR/queries/allele_cnts.txt) $(fn_md5 $OD/allele_cnts.txt)
    assert_exit_code 0
fi

run test_q_allele_cnts_bylen $tdb query allele_cnts_bylen $INDIR/tdb/HG02630_chr14.tdb -o $OD/allele_cnts_bylen.txt
if [ $test_q_allele_cnts_bylen ]; then
    assert_equal $(fn_md5 $INDIR/queries/allele_cnts_bylen.txt) $(fn_md5 $OD/allele_cnts_bylen.txt)
    assert_exit_code 0
fi

run test_q_allele_seqs $tdb query allele_seqs $INDIR/tdb/HG00741_chr14.tdb -O c -o $OD/allele_seqs.csv
if [ $test_q_allele_seqs ]; then
    assert_equal $(fn_md5 $INDIR/queries/allele_seqs.csv) $(fn_md5 $OD/allele_seqs.csv)
    assert_exit_code 0
fi

run test_q_monref $tdb query monref $INDIR/tdb/TwoSamps.tdb -o $OD/monref.txt
if [ $test_q_monref ]; then
    assert_equal $(fn_md5 $INDIR/queries/monref.txt) $(fn_md5 $OD/monref.txt)
    assert_exit_code 0
fi

run test_q_gtmerge $tdb query gtmerge $INDIR/tdb/TwoWithTDB.tdb -o $OD/gtmerge.txt
if [ $test_q_gtmerge ]; then
    assert_equal $(fn_md5 $INDIR/queries/gtmerge.txt) $(fn_md5 $OD/gtmerge.txt)
    assert_exit_code 0
fi

run test_q_metadata $tdb query metadata $INDIR/tdb/TwoWithTDB.tdb -o $OD/metadata.txt
if [ $test_q_metadata ]; then
    assert_equal $(fn_md5 $INDIR/queries/metadata.txt) $(fn_md5 $OD/metadata.txt)
    assert_exit_code 0
fi

run test_q_methyl $tdb query methyl $INDIR/tdb/TwoWithTDB.tdb -O p -o $OD/methyl.pq
if [ $test_q_methyl ]; then
    assert_equal $(fn_md5 $INDIR/queries/methyl.pq) $(fn_md5 $OD/methyl.pq)
    assert_exit_code 0
fi

run test_q_comp_poly_score $tdb query comp_poly_score $INDIR/tdb/TwoWithTDB.tdb -O p -o $OD/comp_poly_score.pq
if [ $test_q_comp_poly_score ]; then
    assert_equal $(fn_md5 $INDIR/queries/comp_poly_score.pq) $(fn_md5 $OD/comp_poly_score.pq)
    assert_exit_code 0
fi

run test_q_len_poly_score $tdb query len_poly_score $INDIR/tdb/TwoWithTDB.tdb -O p -o $OD/len_poly_score.pq
if [ $test_q_len_poly_score ]; then
    assert_equal $(fn_md5 $INDIR/queries/len_poly_score.pq) $(fn_md5 $OD/len_poly_score.pq)
    assert_exit_code 0
fi

TDB_SEED=123 run test_q_saturation $tdb query saturation $INDIR/tdb/merge1.tdb -o $OD/saturation.tsv
if [ $test_saturation ]; then
    assert_equal $(fn_md5 $INDIR/queries/saturation.tsv) $(fn_md5 $OD/saturation.tsv)
    assert_exit_code 0
fi

run test_q_singletons $tdb query singletons $INDIR/tdb/merge1.tdb -o $OD/singletons.tsv
if [ $test_singletons ]; then
    assert_equal $(fn_md5 $INDIR/queries/singletons.tsv) $(fn_md5 $OD/singletons.tsv)
    assert_exit_code 0
fi

# ------------------------------------------------------------
#                                 deid
# ------------------------------------------------------------

run test_deid $tdb deid -o $OD/deid.tdb -i $INDIR/tdb/TwoWithTDB.tdb/
if [ $test_append ]; then
    tdb_check deid.tdb
fi

run test_deid_seq $tdb deid -s -o $OD/deid_seq.tdb -i $INDIR/tdb/TwoWithTDB.tdb/
if [ $test_append ]; then
    tdb_check deid_seq.tdb
fi

TDB_SEED=123 run test_deid_shuf $tdb deid -S -o $OD/deid_shuf.tdb -i $INDIR/tdb/TwoWithTDB.tdb/
if [ $test_deid_shuf ]; then
    tdb_check deid_shuf.tdb
fi

run test_deid_badparam $tdb deid -o repo_utils -i $INDIR/vcf/doesntexist
if [ $test_deid_badparam ]; then
    assert_exit_code 1
fi

# ------------------------------------------------------------
#                                 dump
# ------------------------------------------------------------
run test_dump $tdb dump $INDIR/tdb/HG00438_chr14.tdb -o $OD/HG00438.dump.txt
if [ $test_dump ]; then
    assert_equal $(fn_md5 $INDIR/queries/HG00438.dump.txt) $(fn_md5 $OD/HG00438.dump.txt)
    assert_exit_code 0
fi

# ------------------------------------------------------------
#                                 doctests
# ------------------------------------------------------------
run test_doctests coverage run --concurrency=multiprocessing -p repo_utils/run_doctests.py
if [ $test_doctests ]; then
    assert_exit_code 0 
fi

# ------------------------------------------------------------
#                                 coverage.py
# ------------------------------------------------------------
# Don't generate coverage when doing subset of tests
if [ -z "$1" ]; then
    printf "\n${BOLD}generating test coverage reports${NC}\n"
    coverage combine
    coverage report --include=tdb/*
    coverage html --include=tdb/* -d $OD/htmlcov/
    coverage json --include=tdb/* -o $OD/coverage.json
    python3 repo_utils/coverage_maker.py $OD/coverage.json
fi
rm -f .coverage.*
