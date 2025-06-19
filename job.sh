#!/bin/bash

# Script for executing JOB benchmark queries.
# This script will download query scripts and setup schema, but queries will
# run without any data.
# For executing queries PG environment variables must be available, so
# 'psql' and 'pg_ctl' can be run without any arguments specified.

JOB_SQL_TAR_URL=https://github.com/gregrahn/join-order-benchmark/archive/refs/heads/master.tar.gz
JOB_SQL_TAR_FILE=job.tar.gz

# Download and setup SQL files
mkdir -p job
if [[ ! -f "job/${JOB_SQL_TAR_FILE}" ]]; then
    wget -O "job/${JOB_SQL_TAR_FILE}" https://github.com/gregrahn/join-order-benchmark/archive/refs/heads/master.tar.gz
    tar xvzf "job/${JOB_SQL_TAR_FILE}" --directory=job --strip-components=1
    cat <<EOF >"job/analyze.sql"
ANALYZE movie_companies;
ANALYZE movie_companies;
ANALYZE movie_info_idx;
ANALYZE movie_info;
ANALYZE person_info;
ANALYZE movie_keyword;
ANALYZE aka_title;
ANALYZE title;
ANALYZE movie_link;
ANALYZE movie_link;
ANALYZE aka_title;
ANALYZE cast_info;
ANALYZE complete_cast;
ANALYZE movie_companies;
ANALYZE movie_info_idx;
ANALYZE movie_keyword;
ANALYZE movie_link;
ANALYZE movie_info;
ANALYZE aka_name;
ANALYZE cast_info;
ANALYZE person_info;
ANALYZE cast_info;
ANALYZE cast_info;
EOF
fi
    # Create schema
    psql -f job/schema.sql 
    psql -f job/fkindexes.sql
    psql -f job/analyze.sql

    mv job/schema.sql job/schema.sql.done
    mv job/fkindexes.sql job/fkindexes.sql.done
    mv job/analyze.sql job/analyze.sql.done

    # append 'explain analyze' to all selects
    sed -i -e '0,/^/ s//explain analyze /' job/*.sql
cd job


# Run tests with default planner
pg_ctl stop
echo "shared_preload_libraries=''">>"$PGDATA/postgresql.auto.conf"
pg_ctl start
rm -f vanilla.out
for FILE in $(ls *.sql); do
    psql -f "$FILE" >>vanilla.out
done

# Now with pg_dphyp
pg_ctl stop
echo "shared_preload_libraries='pg_dphyp'">>"$PGDATA/postgresql.auto.conf"
pg_ctl start
rm -f pg_dphyp.out
for FILE in $(ls *.sql); do
    psql -f "$FILE" >>pg_dphyp.out
done

# Check differences
diff -upd -U3 vanilla.out pg_dphyp.out >plan.diff
