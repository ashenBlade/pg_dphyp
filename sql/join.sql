--
-- TEST_SETUP --- prepare environment expected by regression test scripts
--

-- directory paths and dlsuffix are passed to us in environment variables
\getenv abs_srcdir PG_ABS_SRCDIR
\getenv libdir PG_LIBDIR
\getenv dlsuffix PG_DLSUFFIX

\set regresslib :libdir '/regress' :dlsuffix

--
-- synchronous_commit=off delays when hint bits may be set. Some plans change
-- depending on the number of all-visible pages, which in turn can be
-- influenced by the delayed hint bits. Force synchronous_commit=on to avoid
-- that source of variability.
--
SET synchronous_commit = on;

--
-- Postgres formerly made the public schema read/write by default,
-- and most of the core regression tests still expect that.
--
GRANT ALL ON SCHEMA public TO public;

-- Create a tablespace we can use in tests.
SET allow_in_place_tablespaces = true;
CREATE TABLESPACE regress_tblspace LOCATION '';

--
-- These tables have traditionally been referenced by many tests,
-- so create and populate them.  Insert only non-error values here.
-- (Some subsequent tests try to insert erroneous values.  That's okay
-- because the table won't actually change.  Do not change the contents
-- of these tables in later tests, as it may affect other tests.)
--

CREATE TABLE CHAR_TBL(f1 char(4));

INSERT INTO CHAR_TBL (f1) VALUES
  ('a'),
  ('ab'),
  ('abcd'),
  ('abcd    ');
VACUUM CHAR_TBL;

CREATE TABLE FLOAT8_TBL(f1 float8);

INSERT INTO FLOAT8_TBL(f1) VALUES
  ('0.0'),
  ('-34.84'),
  ('-1004.30'),
  ('-1.2345678901234e+200'),
  ('-1.2345678901234e-200');
VACUUM FLOAT8_TBL;

CREATE TABLE INT2_TBL(f1 int2);

INSERT INTO INT2_TBL(f1) VALUES
  ('0   '),
  ('  1234 '),
  ('    -1234'),
  ('32767'),  -- largest and smallest values
  ('-32767');
VACUUM INT2_TBL;

CREATE TABLE INT4_TBL(f1 int4);

INSERT INTO INT4_TBL(f1) VALUES
  ('   0  '),
  ('123456     '),
  ('    -123456'),
  ('2147483647'),  -- largest and smallest values
  ('-2147483647');
VACUUM INT4_TBL;

CREATE TABLE INT8_TBL(q1 int8, q2 int8);

INSERT INTO INT8_TBL VALUES
  ('  123   ','  456'),
  ('123   ','4567890123456789'),
  ('4567890123456789','123'),
  (+4567890123456789,'4567890123456789'),
  ('+4567890123456789','-4567890123456789');
VACUUM INT8_TBL;

CREATE TABLE POINT_TBL(f1 point);

INSERT INTO POINT_TBL(f1) VALUES
  ('(0.0,0.0)'),
  ('(-10.0,0.0)'),
  ('(-3.0,4.0)'),
  ('(5.1, 34.5)'),
  ('(-5.0,-12.0)'),
  ('(1e-300,-1e-300)'),  -- To underflow
  ('(1e+300,Inf)'),  -- To overflow
  ('(Inf,1e+300)'),  -- Transposed
  (' ( Nan , NaN ) '),
  ('10.0,10.0');
-- We intentionally don't vacuum point_tbl here; geometry depends on that

CREATE TABLE TEXT_TBL (f1 text);

INSERT INTO TEXT_TBL VALUES
  ('doh!'),
  ('hi de ho neighbor');
VACUUM TEXT_TBL;

CREATE TABLE VARCHAR_TBL(f1 varchar(4));

INSERT INTO VARCHAR_TBL (f1) VALUES
  ('a'),
  ('ab'),
  ('abcd'),
  ('abcd    ');
VACUUM VARCHAR_TBL;

CREATE TABLE onek (
	unique1		int4,
	unique2		int4,
	two			int4,
	four		int4,
	ten			int4,
	twenty		int4,
	hundred		int4,
	thousand	int4,
	twothousand	int4,
	fivethous	int4,
	tenthous	int4,
	odd			int4,
	even		int4,
	stringu1	name,
	stringu2	name,
	string4		name
);

\set filename :abs_srcdir '/../../src/test/regress/data/onek.data'
COPY onek FROM :'filename';
VACUUM ANALYZE onek;

CREATE TABLE onek2 AS SELECT * FROM onek;
VACUUM ANALYZE onek2;

CREATE TABLE tenk1 (
	unique1		int4,
	unique2		int4,
	two			int4,
	four		int4,
	ten			int4,
	twenty		int4,
	hundred		int4,
	thousand	int4,
	twothousand	int4,
	fivethous	int4,
	tenthous	int4,
	odd			int4,
	even		int4,
	stringu1	name,
	stringu2	name,
	string4		name
);

\set filename :abs_srcdir '/../../src/test/regress/data/tenk.data'
COPY tenk1 FROM :'filename';
VACUUM ANALYZE tenk1;

CREATE TABLE tenk2 AS SELECT * FROM tenk1;
VACUUM ANALYZE tenk2;

CREATE TABLE person (
	name 		text,
	age			int4,
	location 	point
);

\set filename :abs_srcdir '/../../src/test/regress/data/person.data'
COPY person FROM :'filename';
VACUUM ANALYZE person;

CREATE TABLE emp (
	salary 		int4,
	manager 	name
) INHERITS (person);

\set filename :abs_srcdir '/../../src/test/regress/data/emp.data'
COPY emp FROM :'filename';
VACUUM ANALYZE emp;

CREATE TABLE student (
	gpa 		float8
) INHERITS (person);

\set filename :abs_srcdir '/../../src/test/regress/data/student.data'
COPY student FROM :'filename';
VACUUM ANALYZE student;

CREATE TABLE stud_emp (
	percent 	int4
) INHERITS (emp, student);

\set filename :abs_srcdir '/../../src/test/regress/data/stud_emp.data'
COPY stud_emp FROM :'filename';
VACUUM ANALYZE stud_emp;

CREATE TABLE road (
	name		text,
	thepath 	path
);

\set filename :abs_srcdir '/../../src/test/regress/data/streets.data'
COPY road FROM :'filename';
VACUUM ANALYZE road;

CREATE TABLE ihighway () INHERITS (road);

INSERT INTO ihighway
   SELECT *
   FROM ONLY road
   WHERE name ~ 'I- .*';
VACUUM ANALYZE ihighway;

CREATE TABLE shighway (
	surface		text
) INHERITS (road);

INSERT INTO shighway
   SELECT *, 'asphalt'
   FROM ONLY road
   WHERE name ~ 'State Hwy.*';
VACUUM ANALYZE shighway;

--
-- We must have some enum type in the database for opr_sanity and type_sanity.
--

create type stoplight as enum ('red', 'yellow', 'green');

--
-- Also create some non-built-in range types.
--

create type float8range as range (subtype = float8, subtype_diff = float8mi);

create type textrange as range (subtype = text, collation = "C");

-- Use hand-rolled hash functions and operator classes to get predictable
-- result on different machines.  The hash function for int4 simply returns
-- the sum of the values passed to it and the one for text returns the length
-- of the non-empty string value passed to it or 0.

create function part_hashint4_noop(value int4, seed int8)
    returns int8 as $$
    select value + seed;
    $$ language sql strict immutable parallel safe;

create operator class part_test_int4_ops for type int4 using hash as
    operator 1 =,
    function 2 part_hashint4_noop(int4, int8);

create function part_hashtext_length(value text, seed int8)
    returns int8 as $$
    select length(coalesce(value, ''))::int8
    $$ language sql strict immutable parallel safe;

create operator class part_test_text_ops for type text using hash as
    operator 1 =,
    function 2 part_hashtext_length(text, int8);

--
-- These functions are used in tests that used to use md5(), which we now
-- mostly avoid so that the tests will pass in FIPS mode.
--

create function fipshash(bytea)
    returns text
    strict immutable parallel safe leakproof
    return substr(encode(sha256($1), 'hex'), 1, 32);

create function fipshash(text)
    returns text
    strict immutable parallel safe leakproof
    return substr(encode(sha256($1::bytea), 'hex'), 1, 32);

-- directory paths are passed to us in environment variables
\getenv abs_srcdir PG_ABS_SRCDIR

--
-- BTREE
--
CREATE INDEX onek_unique1 ON onek USING btree(unique1 int4_ops);

CREATE INDEX IF NOT EXISTS onek_unique1 ON onek USING btree(unique1 int4_ops);

CREATE INDEX IF NOT EXISTS ON onek USING btree(unique1 int4_ops);

CREATE INDEX onek_unique2 ON onek USING btree(unique2 int4_ops);

CREATE INDEX onek_hundred ON onek USING btree(hundred int4_ops);

CREATE INDEX onek_stringu1 ON onek USING btree(stringu1 name_ops);

CREATE INDEX tenk1_unique1 ON tenk1 USING btree(unique1 int4_ops);

CREATE INDEX tenk1_unique2 ON tenk1 USING btree(unique2 int4_ops);

CREATE INDEX tenk1_hundred ON tenk1 USING btree(hundred int4_ops);

CREATE INDEX tenk1_thous_tenthous ON tenk1 (thousand, tenthous);

CREATE INDEX tenk2_unique1 ON tenk2 USING btree(unique1 int4_ops);

CREATE INDEX tenk2_unique2 ON tenk2 USING btree(unique2 int4_ops);

CREATE INDEX tenk2_hundred ON tenk2 USING btree(hundred int4_ops);

CREATE INDEX rix ON road USING btree (name text_ops);

CREATE INDEX iix ON ihighway USING btree (name text_ops);

CREATE INDEX six ON shighway USING btree (name text_ops);

--
-- BTREE partial indices
--
CREATE INDEX onek2_u1_prtl ON onek2 USING btree(unique1 int4_ops)
	where unique1 < 20 or unique1 > 980;

CREATE INDEX onek2_u2_prtl ON onek2 USING btree(unique2 int4_ops)
	where stringu1 < 'B';

CREATE INDEX onek2_stu1_prtl ON onek2 USING btree(stringu1 name_ops)
	where onek2.stringu1 >= 'J' and onek2.stringu1 < 'K';

--
-- CREATE_MISC
--

--
-- a is the type root
-- b and c inherit from a (one-level single inheritance)
-- d inherits from b and c (two-level multiple inheritance)
-- e inherits from c (two-level single inheritance)
-- f inherits from e (three-level single inheritance)
--
CREATE TABLE a_star (
	class		char,
	a 			int4
);

CREATE TABLE b_star (
	b 			text
) INHERITS (a_star);

CREATE TABLE c_star (
	c 			name
) INHERITS (a_star);

CREATE TABLE d_star (
	d 			float8
) INHERITS (b_star, c_star);

CREATE TABLE e_star (
	e 			int2
) INHERITS (c_star);

CREATE TABLE f_star (
	f 			polygon
) INHERITS (e_star);

INSERT INTO a_star (class, a) VALUES ('a', 1);

INSERT INTO a_star (class, a) VALUES ('a', 2);

INSERT INTO a_star (class) VALUES ('a');

INSERT INTO b_star (class, a, b) VALUES ('b', 3, 'mumble'::text);

INSERT INTO b_star (class, a) VALUES ('b', 4);

INSERT INTO b_star (class, b) VALUES ('b', 'bumble'::text);

INSERT INTO b_star (class) VALUES ('b');

INSERT INTO c_star (class, a, c) VALUES ('c', 5, 'hi mom'::name);

INSERT INTO c_star (class, a) VALUES ('c', 6);

INSERT INTO c_star (class, c) VALUES ('c', 'hi paul'::name);

INSERT INTO c_star (class) VALUES ('c');

INSERT INTO d_star (class, a, b, c, d)
   VALUES ('d', 7, 'grumble'::text, 'hi sunita'::name, '0.0'::float8);

INSERT INTO d_star (class, a, b, c)
   VALUES ('d', 8, 'stumble'::text, 'hi koko'::name);

INSERT INTO d_star (class, a, b, d)
   VALUES ('d', 9, 'rumble'::text, '1.1'::float8);

INSERT INTO d_star (class, a, c, d)
   VALUES ('d', 10, 'hi kristin'::name, '10.01'::float8);

INSERT INTO d_star (class, b, c, d)
   VALUES ('d', 'crumble'::text, 'hi boris'::name, '100.001'::float8);

INSERT INTO d_star (class, a, b)
   VALUES ('d', 11, 'fumble'::text);

INSERT INTO d_star (class, a, c)
   VALUES ('d', 12, 'hi avi'::name);

INSERT INTO d_star (class, a, d)
   VALUES ('d', 13, '1000.0001'::float8);

INSERT INTO d_star (class, b, c)
   VALUES ('d', 'tumble'::text, 'hi andrew'::name);

INSERT INTO d_star (class, b, d)
   VALUES ('d', 'humble'::text, '10000.00001'::float8);

INSERT INTO d_star (class, c, d)
   VALUES ('d', 'hi ginger'::name, '100000.000001'::float8);

INSERT INTO d_star (class, a) VALUES ('d', 14);

INSERT INTO d_star (class, b) VALUES ('d', 'jumble'::text);

INSERT INTO d_star (class, c) VALUES ('d', 'hi jolly'::name);

INSERT INTO d_star (class, d) VALUES ('d', '1000000.0000001'::float8);

INSERT INTO d_star (class) VALUES ('d');

INSERT INTO e_star (class, a, c, e)
   VALUES ('e', 15, 'hi carol'::name, '-1'::int2);

INSERT INTO e_star (class, a, c)
   VALUES ('e', 16, 'hi bob'::name);

INSERT INTO e_star (class, a, e)
   VALUES ('e', 17, '-2'::int2);

INSERT INTO e_star (class, c, e)
   VALUES ('e', 'hi michelle'::name, '-3'::int2);

INSERT INTO e_star (class, a)
   VALUES ('e', 18);

INSERT INTO e_star (class, c)
   VALUES ('e', 'hi elisa'::name);

INSERT INTO e_star (class, e)
   VALUES ('e', '-4'::int2);

INSERT INTO f_star (class, a, c, e, f)
   VALUES ('f', 19, 'hi claire'::name, '-5'::int2, '(1,3),(2,4)'::polygon);

INSERT INTO f_star (class, a, c, e)
   VALUES ('f', 20, 'hi mike'::name, '-6'::int2);

INSERT INTO f_star (class, a, c, f)
   VALUES ('f', 21, 'hi marcel'::name, '(11,44),(22,55),(33,66)'::polygon);

INSERT INTO f_star (class, a, e, f)
   VALUES ('f', 22, '-7'::int2, '(111,555),(222,666),(333,777),(444,888)'::polygon);

INSERT INTO f_star (class, c, e, f)
   VALUES ('f', 'hi keith'::name, '-8'::int2,
	   '(1111,3333),(2222,4444)'::polygon);

INSERT INTO f_star (class, a, c)
   VALUES ('f', 24, 'hi marc'::name);

INSERT INTO f_star (class, a, e)
   VALUES ('f', 25, '-9'::int2);

INSERT INTO f_star (class, a, f)
   VALUES ('f', 26, '(11111,33333),(22222,44444)'::polygon);

INSERT INTO f_star (class, c, e)
   VALUES ('f', 'hi allison'::name, '-10'::int2);

INSERT INTO f_star (class, c, f)
   VALUES ('f', 'hi jeff'::name,
           '(111111,333333),(222222,444444)'::polygon);

INSERT INTO f_star (class, e, f)
   VALUES ('f', '-11'::int2, '(1111111,3333333),(2222222,4444444)'::polygon);

INSERT INTO f_star (class, a) VALUES ('f', 27);

INSERT INTO f_star (class, c) VALUES ('f', 'hi carl'::name);

INSERT INTO f_star (class, e) VALUES ('f', '-12'::int2);

INSERT INTO f_star (class, f)
   VALUES ('f', '(11111111,33333333),(22222222,44444444)'::polygon);

INSERT INTO f_star (class) VALUES ('f');

-- Analyze the X_star tables for better plan stability in later tests
ANALYZE a_star;
ANALYZE b_star;
ANALYZE c_star;
ANALYZE d_star;
ANALYZE e_star;
ANALYZE f_star;

--
-- inheritance stress test
--
SELECT * FROM a_star*;

SELECT *
   FROM b_star* x
   WHERE x.b = text 'bumble' or x.a < 3;

SELECT class, a
   FROM c_star* x
   WHERE x.c ~ text 'hi';

SELECT class, b, c
   FROM d_star* x
   WHERE x.a < 100;

SELECT class, c FROM e_star* x WHERE x.c NOTNULL;

SELECT * FROM f_star* x WHERE x.c ISNULL;

-- grouping and aggregation on inherited sets have been busted in the past...

SELECT sum(a) FROM a_star*;

SELECT class, sum(a) FROM a_star* GROUP BY class ORDER BY class;


ALTER TABLE f_star RENAME COLUMN f TO ff;

ALTER TABLE e_star* RENAME COLUMN e TO ee;

ALTER TABLE d_star* RENAME COLUMN d TO dd;

ALTER TABLE c_star* RENAME COLUMN c TO cc;

ALTER TABLE b_star* RENAME COLUMN b TO bb;

ALTER TABLE a_star* RENAME COLUMN a TO aa;

SELECT class, aa
   FROM a_star* x
   WHERE aa ISNULL;

-- As of Postgres 7.1, ALTER implicitly recurses,
-- so this should be same as ALTER a_star*

ALTER TABLE a_star RENAME COLUMN aa TO foo;

SELECT class, foo
   FROM a_star* x
   WHERE x.foo >= 2;

ALTER TABLE a_star RENAME COLUMN foo TO aa;

SELECT *
   from a_star*
   WHERE aa < 1000;

ALTER TABLE f_star ADD COLUMN f int4;

UPDATE f_star SET f = 10;

ALTER TABLE e_star* ADD COLUMN e int4;

--UPDATE e_star* SET e = 42;

SELECT * FROM e_star*;

ALTER TABLE a_star* ADD COLUMN a text;

-- That ALTER TABLE should have added TOAST tables.
SELECT relname, reltoastrelid <> 0 AS has_toast_table
   FROM pg_class
   WHERE oid::regclass IN ('a_star', 'c_star')
   ORDER BY 1;

--UPDATE b_star*
--   SET a = text 'gazpacho'
--   WHERE aa > 4;

SELECT class, aa, a FROM a_star*;

--
-- JOIN
-- Test JOIN clauses
--

CREATE TABLE J1_TBL (
  i integer,
  j integer,
  t text
);

CREATE TABLE J2_TBL (
  i integer,
  k integer
);


INSERT INTO J1_TBL VALUES (1, 4, 'one');
INSERT INTO J1_TBL VALUES (2, 3, 'two');
INSERT INTO J1_TBL VALUES (3, 2, 'three');
INSERT INTO J1_TBL VALUES (4, 1, 'four');
INSERT INTO J1_TBL VALUES (5, 0, 'five');
INSERT INTO J1_TBL VALUES (6, 6, 'six');
INSERT INTO J1_TBL VALUES (7, 7, 'seven');
INSERT INTO J1_TBL VALUES (8, 8, 'eight');
INSERT INTO J1_TBL VALUES (0, NULL, 'zero');
INSERT INTO J1_TBL VALUES (NULL, NULL, 'null');
INSERT INTO J1_TBL VALUES (NULL, 0, 'zero');

INSERT INTO J2_TBL VALUES (1, -1);
INSERT INTO J2_TBL VALUES (2, 2);
INSERT INTO J2_TBL VALUES (3, -3);
INSERT INTO J2_TBL VALUES (2, 4);
INSERT INTO J2_TBL VALUES (5, -5);
INSERT INTO J2_TBL VALUES (5, -5);
INSERT INTO J2_TBL VALUES (0, NULL);
INSERT INTO J2_TBL VALUES (NULL, NULL);
INSERT INTO J2_TBL VALUES (NULL, 0);

-- useful in some tests below
create temp table onerow();
insert into onerow default values;
analyze onerow;


--
-- CORRELATION NAMES
-- Make sure that table/column aliases are supported
-- before diving into more complex join syntax.
--

SELECT *
  FROM J1_TBL AS tx;

SELECT *
  FROM J1_TBL tx;

SELECT *
  FROM J1_TBL AS t1 (a, b, c);

SELECT *
  FROM J1_TBL t1 (a, b, c);

SELECT *
  FROM J1_TBL t1 (a, b, c), J2_TBL t2 (d, e);

SELECT t1.a, t2.e
  FROM J1_TBL t1 (a, b, c), J2_TBL t2 (d, e)
  WHERE t1.a = t2.d;


--
-- CROSS JOIN
-- Qualifications are not allowed on cross joins,
-- which degenerate into a standard unqualified inner join.
--

SELECT *
  FROM J1_TBL CROSS JOIN J2_TBL;

-- ambiguous column
SELECT i, k, t
  FROM J1_TBL CROSS JOIN J2_TBL;

-- resolve previous ambiguity by specifying the table name
SELECT t1.i, k, t
  FROM J1_TBL t1 CROSS JOIN J2_TBL t2;

SELECT ii, tt, kk
  FROM (J1_TBL CROSS JOIN J2_TBL)
    AS tx (ii, jj, tt, ii2, kk);

SELECT tx.ii, tx.jj, tx.kk
  FROM (J1_TBL t1 (a, b, c) CROSS JOIN J2_TBL t2 (d, e))
    AS tx (ii, jj, tt, ii2, kk);

SELECT *
  FROM J1_TBL CROSS JOIN J2_TBL a CROSS JOIN J2_TBL b;


--
--
-- Inner joins (equi-joins)
--
--

--
-- Inner joins (equi-joins) with USING clause
-- The USING syntax changes the shape of the resulting table
-- by including a column in the USING clause only once in the result.
--

-- Inner equi-join on specified column
SELECT *
  FROM J1_TBL INNER JOIN J2_TBL USING (i);

-- Same as above, slightly different syntax
SELECT *
  FROM J1_TBL JOIN J2_TBL USING (i);

SELECT *
  FROM J1_TBL t1 (a, b, c) JOIN J2_TBL t2 (a, d) USING (a)
  ORDER BY a, d;

SELECT *
  FROM J1_TBL t1 (a, b, c) JOIN J2_TBL t2 (a, b) USING (b)
  ORDER BY b, t1.a;

-- test join using aliases
SELECT * FROM J1_TBL JOIN J2_TBL USING (i) WHERE J1_TBL.t = 'one';  -- ok
SELECT * FROM J1_TBL JOIN J2_TBL USING (i) AS x WHERE J1_TBL.t = 'one';  -- ok
SELECT * FROM (J1_TBL JOIN J2_TBL USING (i)) AS x WHERE J1_TBL.t = 'one';  -- error
SELECT * FROM J1_TBL JOIN J2_TBL USING (i) AS x WHERE x.i = 1;  -- ok
SELECT * FROM J1_TBL JOIN J2_TBL USING (i) AS x WHERE x.t = 'one';  -- error
SELECT * FROM (J1_TBL JOIN J2_TBL USING (i) AS x) AS xx WHERE x.i = 1;  -- error (XXX could use better hint)
SELECT * FROM J1_TBL a1 JOIN J2_TBL a2 USING (i) AS a1;  -- error
SELECT x.* FROM J1_TBL JOIN J2_TBL USING (i) AS x WHERE J1_TBL.t = 'one';
SELECT ROW(x.*) FROM J1_TBL JOIN J2_TBL USING (i) AS x WHERE J1_TBL.t = 'one';
SELECT row_to_json(x.*) FROM J1_TBL JOIN J2_TBL USING (i) AS x WHERE J1_TBL.t = 'one';

--
-- NATURAL JOIN
-- Inner equi-join on all columns with the same name
--

SELECT *
  FROM J1_TBL NATURAL JOIN J2_TBL;

SELECT *
  FROM J1_TBL t1 (a, b, c) NATURAL JOIN J2_TBL t2 (a, d);

SELECT *
  FROM J1_TBL t1 (a, b, c) NATURAL JOIN J2_TBL t2 (d, a);

-- mismatch number of columns
-- currently, Postgres will fill in with underlying names
SELECT *
  FROM J1_TBL t1 (a, b) NATURAL JOIN J2_TBL t2 (a);


--
-- Inner joins (equi-joins)
--

SELECT *
  FROM J1_TBL JOIN J2_TBL ON (J1_TBL.i = J2_TBL.i);

SELECT *
  FROM J1_TBL JOIN J2_TBL ON (J1_TBL.i = J2_TBL.k);


--
-- Non-equi-joins
--

SELECT *
  FROM J1_TBL JOIN J2_TBL ON (J1_TBL.i <= J2_TBL.k);


--
-- Outer joins
-- Note that OUTER is a noise word
--

SELECT *
  FROM J1_TBL LEFT OUTER JOIN J2_TBL USING (i)
  ORDER BY i, k, t;

SELECT *
  FROM J1_TBL LEFT JOIN J2_TBL USING (i)
  ORDER BY i, k, t;

SELECT *
  FROM J1_TBL RIGHT OUTER JOIN J2_TBL USING (i);

SELECT *
  FROM J1_TBL RIGHT JOIN J2_TBL USING (i);

SELECT *
  FROM J1_TBL FULL OUTER JOIN J2_TBL USING (i)
  ORDER BY i, k, t;

SELECT *
  FROM J1_TBL FULL JOIN J2_TBL USING (i)
  ORDER BY i, k, t;

SELECT *
  FROM J1_TBL LEFT JOIN J2_TBL USING (i) WHERE (k = 1);

SELECT *
  FROM J1_TBL LEFT JOIN J2_TBL USING (i) WHERE (i = 1);

--
-- semijoin selectivity for <>
--
explain (costs off)
select * from int4_tbl i4, tenk1 a
where exists(select * from tenk1 b
             where a.twothousand = b.twothousand and a.fivethous <> b.fivethous)
      and i4.f1 = a.tenthous;


--
-- More complicated constructs
--

--
-- Multiway full join
--

CREATE TABLE t1 (name TEXT, n INTEGER);
CREATE TABLE t2 (name TEXT, n INTEGER);
CREATE TABLE t3 (name TEXT, n INTEGER);

INSERT INTO t1 VALUES ( 'bb', 11 );
INSERT INTO t2 VALUES ( 'bb', 12 );
INSERT INTO t2 VALUES ( 'cc', 22 );
INSERT INTO t2 VALUES ( 'ee', 42 );
INSERT INTO t3 VALUES ( 'bb', 13 );
INSERT INTO t3 VALUES ( 'cc', 23 );
INSERT INTO t3 VALUES ( 'dd', 33 );

SELECT * FROM t1 FULL JOIN t2 USING (name) FULL JOIN t3 USING (name);

--
-- Test interactions of join syntax and subqueries
--

-- Basic cases (we expect planner to pull up the subquery here)
SELECT * FROM
(SELECT * FROM t2) as s2
INNER JOIN
(SELECT * FROM t3) s3
USING (name);

SELECT * FROM
(SELECT * FROM t2) as s2
LEFT JOIN
(SELECT * FROM t3) s3
USING (name);

SELECT * FROM
(SELECT * FROM t2) as s2
FULL JOIN
(SELECT * FROM t3) s3
USING (name);

-- Cases with non-nullable expressions in subquery results;
-- make sure these go to null as expected
SELECT * FROM
(SELECT name, n as s2_n, 2 as s2_2 FROM t2) as s2
NATURAL INNER JOIN
(SELECT name, n as s3_n, 3 as s3_2 FROM t3) s3;

SELECT * FROM
(SELECT name, n as s2_n, 2 as s2_2 FROM t2) as s2
NATURAL LEFT JOIN
(SELECT name, n as s3_n, 3 as s3_2 FROM t3) s3;

SELECT * FROM
(SELECT name, n as s2_n, 2 as s2_2 FROM t2) as s2
NATURAL FULL JOIN
(SELECT name, n as s3_n, 3 as s3_2 FROM t3) s3;

SELECT * FROM
(SELECT name, n as s1_n, 1 as s1_1 FROM t1) as s1
NATURAL INNER JOIN
(SELECT name, n as s2_n, 2 as s2_2 FROM t2) as s2
NATURAL INNER JOIN
(SELECT name, n as s3_n, 3 as s3_2 FROM t3) s3;

SELECT * FROM
(SELECT name, n as s1_n, 1 as s1_1 FROM t1) as s1
NATURAL FULL JOIN
(SELECT name, n as s2_n, 2 as s2_2 FROM t2) as s2
NATURAL FULL JOIN
(SELECT name, n as s3_n, 3 as s3_2 FROM t3) s3;

SELECT * FROM
(SELECT name, n as s1_n FROM t1) as s1
NATURAL FULL JOIN
  (SELECT * FROM
    (SELECT name, n as s2_n FROM t2) as s2
    NATURAL FULL JOIN
    (SELECT name, n as s3_n FROM t3) as s3
  ) ss2;

SELECT * FROM
(SELECT name, n as s1_n FROM t1) as s1
NATURAL FULL JOIN
  (SELECT * FROM
    (SELECT name, n as s2_n, 2 as s2_2 FROM t2) as s2
    NATURAL FULL JOIN
    (SELECT name, n as s3_n FROM t3) as s3
  ) ss2;

-- Constants as join keys can also be problematic
SELECT * FROM
  (SELECT name, n as s1_n FROM t1) as s1
FULL JOIN
  (SELECT name, 2 as s2_n FROM t2) as s2
ON (s1_n = s2_n);


-- Test for propagation of nullability constraints into sub-joins

create temp table x (x1 int, x2 int);
insert into x values (1,11);
insert into x values (2,22);
insert into x values (3,null);
insert into x values (4,44);
insert into x values (5,null);

create temp table y (y1 int, y2 int);
insert into y values (1,111);
insert into y values (2,222);
insert into y values (3,333);
insert into y values (4,null);

select * from x;
select * from y;

select * from x left join y on (x1 = y1 and x2 is not null);
select * from x left join y on (x1 = y1 and y2 is not null);

select * from (x left join y on (x1 = y1)) left join x xx(xx1,xx2)
on (x1 = xx1);
select * from (x left join y on (x1 = y1)) left join x xx(xx1,xx2)
on (x1 = xx1 and x2 is not null);
select * from (x left join y on (x1 = y1)) left join x xx(xx1,xx2)
on (x1 = xx1 and y2 is not null);
select * from (x left join y on (x1 = y1)) left join x xx(xx1,xx2)
on (x1 = xx1 and xx2 is not null);
-- these should NOT give the same answers as above
select * from (x left join y on (x1 = y1)) left join x xx(xx1,xx2)
on (x1 = xx1) where (x2 is not null);
select * from (x left join y on (x1 = y1)) left join x xx(xx1,xx2)
on (x1 = xx1) where (y2 is not null);
select * from (x left join y on (x1 = y1)) left join x xx(xx1,xx2)
on (x1 = xx1) where (xx2 is not null);

--
-- regression test: check for bug with propagation of implied equality
-- to outside an IN
--
select count(*) from tenk1 a where unique1 in
  (select unique1 from tenk1 b join tenk1 c using (unique1)
   where b.unique2 = 42);

--
-- regression test: check for failure to generate a plan with multiple
-- degenerate IN clauses
--
select count(*) from tenk1 x where
  x.unique1 in (select a.f1 from int4_tbl a,float8_tbl b where a.f1=b.f1) and
  x.unique1 = 0 and
  x.unique1 in (select aa.f1 from int4_tbl aa,float8_tbl bb where aa.f1=bb.f1);

-- try that with GEQO too
begin;
set geqo = on;
set geqo_threshold = 2;
select count(*) from tenk1 x where
  x.unique1 in (select a.f1 from int4_tbl a,float8_tbl b where a.f1=b.f1) and
  x.unique1 = 0 and
  x.unique1 in (select aa.f1 from int4_tbl aa,float8_tbl bb where aa.f1=bb.f1);
rollback;

--
-- regression test: be sure we cope with proven-dummy append rels
--
explain (costs off)
select aa, bb, unique1, unique1
  from tenk1 right join b_star on aa = unique1
  where bb < bb and bb is null;

select aa, bb, unique1, unique1
  from tenk1 right join b_star on aa = unique1
  where bb < bb and bb is null;

--
-- regression test: check handling of empty-FROM subquery underneath outer join
--
explain (costs off)
select * from int8_tbl i1 left join (int8_tbl i2 join
  (select 123 as x) ss on i2.q1 = x) on i1.q2 = i2.q2
order by 1, 2;

select * from int8_tbl i1 left join (int8_tbl i2 join
  (select 123 as x) ss on i2.q1 = x) on i1.q2 = i2.q2
order by 1, 2;

--
-- regression test: check a case where join_clause_is_movable_into()
-- used to give an imprecise result, causing an assertion failure
--
select count(*)
from
  (select t3.tenthous as x1, coalesce(t1.stringu1, t2.stringu1) as x2
   from tenk1 t1
   left join tenk1 t2 on t1.unique1 = t2.unique1
   join tenk1 t3 on t1.unique2 = t3.unique2) ss,
  tenk1 t4,
  tenk1 t5
where t4.thousand = t5.unique1 and ss.x1 = t4.tenthous and ss.x2 = t5.stringu1;

--
-- regression test: check a case where we formerly missed including an EC
-- enforcement clause because it was expected to be handled at scan level
--
explain (costs off)
select a.f1, b.f1, t.thousand, t.tenthous from
  tenk1 t,
  (select sum(f1)+1 as f1 from int4_tbl i4a) a,
  (select sum(f1) as f1 from int4_tbl i4b) b
where b.f1 = t.thousand and a.f1 = b.f1 and (a.f1+b.f1+999) = t.tenthous;

select a.f1, b.f1, t.thousand, t.tenthous from
  tenk1 t,
  (select sum(f1)+1 as f1 from int4_tbl i4a) a,
  (select sum(f1) as f1 from int4_tbl i4b) b
where b.f1 = t.thousand and a.f1 = b.f1 and (a.f1+b.f1+999) = t.tenthous;

--
-- checks for correct handling of quals in multiway outer joins
--
explain (costs off)
select t1.f1
from int4_tbl t1, int4_tbl t2
  left join int4_tbl t3 on t3.f1 > 0
  left join int4_tbl t4 on t3.f1 > 1
where t4.f1 is null;

select t1.f1
from int4_tbl t1, int4_tbl t2
  left join int4_tbl t3 on t3.f1 > 0
  left join int4_tbl t4 on t3.f1 > 1
where t4.f1 is null;

explain (costs off)
select *
from int4_tbl t1 left join int4_tbl t2 on true
  left join int4_tbl t3 on t2.f1 > 0
  left join int4_tbl t4 on t3.f1 > 0;

explain (costs off)
select * from onek t1
  left join onek t2 on t1.unique1 = t2.unique1
  left join onek t3 on t2.unique1 != t3.unique1
  left join onek t4 on t3.unique1 = t4.unique1;

explain (costs off)
select * from int4_tbl t1
  left join (select now() from int4_tbl t2
             left join int4_tbl t3 on t2.f1 = t3.f1
             left join int4_tbl t4 on t3.f1 = t4.f1) s on true
  inner join int4_tbl t5 on true;

explain (costs off)
select * from int4_tbl t1
  left join int4_tbl t2 on true
  left join int4_tbl t3 on true
  left join int4_tbl t4 on t2.f1 = t3.f1;

explain (costs off)
select * from int4_tbl t1
  left join int4_tbl t2 on true
  left join int4_tbl t3 on t2.f1 = t3.f1
  left join int4_tbl t4 on t3.f1 != t4.f1;

explain (costs off)
select * from int4_tbl t1
  left join (int4_tbl t2 left join int4_tbl t3 on t2.f1 > 0) on t2.f1 > 1
  left join int4_tbl t4 on t2.f1 > 2 and t3.f1 > 3
where t1.f1 = coalesce(t2.f1, 1);

explain (costs off)
select * from int4_tbl t1
  left join ((select t2.f1 from int4_tbl t2
                left join int4_tbl t3 on t2.f1 > 0
                where t3.f1 is null) s
             left join tenk1 t4 on s.f1 > 1)
    on s.f1 = t1.f1;

explain (costs off)
select * from int4_tbl t1
  left join ((select t2.f1 from int4_tbl t2
                left join int4_tbl t3 on t2.f1 > 0
                where t2.f1 <> coalesce(t3.f1, -1)) s
             left join tenk1 t4 on s.f1 > 1)
    on s.f1 = t1.f1;

explain (costs off)
select * from onek t1
    left join onek t2 on t1.unique1 = t2.unique1
    left join onek t3 on t2.unique1 = t3.unique1
    left join onek t4 on t3.unique1 = t4.unique1 and t2.unique2 = t4.unique2;

explain (costs off)
select * from int8_tbl t1 left join
    (int8_tbl t2 left join int8_tbl t3 full join int8_tbl t4 on false on false)
    left join int8_tbl t5 on t2.q1 = t5.q1
on t2.q2 = 123;

explain (costs off)
select * from int8_tbl t1
    left join int8_tbl t2 on true
    left join lateral
      (select * from int8_tbl t3 where t3.q1 = t2.q1 offset 0) s
      on t2.q1 = 1;

explain (costs off)
select * from int8_tbl t1
    left join int8_tbl t2 on true
    left join lateral
      (select * from generate_series(t2.q1, 100)) s
      on t2.q1 = 1;

explain (costs off)
select * from int8_tbl t1
    left join int8_tbl t2 on true
    left join lateral
      (select t2.q1 from int8_tbl t3) s
      on t2.q1 = 1;

explain (costs off)
select * from onek t1
    left join onek t2 on true
    left join lateral
      (select * from onek t3 where t3.two = t2.two offset 0) s
      on t2.unique1 = 1;

--
-- check a case where we formerly got confused by conflicting sort orders
-- in redundant merge join path keys
--
explain (costs off)
select * from
  j1_tbl full join
  (select * from j2_tbl order by j2_tbl.i desc, j2_tbl.k asc) j2_tbl
  on j1_tbl.i = j2_tbl.i and j1_tbl.i = j2_tbl.k;

select * from
  j1_tbl full join
  (select * from j2_tbl order by j2_tbl.i desc, j2_tbl.k asc) j2_tbl
  on j1_tbl.i = j2_tbl.i and j1_tbl.i = j2_tbl.k;

--
-- a different check for handling of redundant sort keys in merge joins
--
explain (costs off)
select count(*) from
  (select * from tenk1 x order by x.thousand, x.twothousand, x.fivethous) x
  left join
  (select * from tenk1 y order by y.unique2) y
  on x.thousand = y.unique2 and x.twothousand = y.hundred and x.fivethous = y.unique2;

select count(*) from
  (select * from tenk1 x order by x.thousand, x.twothousand, x.fivethous) x
  left join
  (select * from tenk1 y order by y.unique2) y
  on x.thousand = y.unique2 and x.twothousand = y.hundred and x.fivethous = y.unique2;

set enable_hashjoin = 0;
set enable_nestloop = 0;
set enable_hashagg = 0;

--
-- Check that we use the pathkeys from a prefix of the group by / order by
-- clause for the join pathkeys when that prefix covers all join quals.  We
-- expect this to lead to an incremental sort for the group by / order by.
--
explain (costs off)
select x.thousand, x.twothousand, count(*)
from tenk1 x inner join tenk1 y on x.thousand = y.thousand
group by x.thousand, x.twothousand
order by x.thousand desc, x.twothousand;

reset enable_hashagg;
reset enable_nestloop;
reset enable_hashjoin;

--
-- Clean up
--

DROP TABLE t1;
DROP TABLE t2;
DROP TABLE t3;

DROP TABLE J1_TBL;
DROP TABLE J2_TBL;

-- Both DELETE and UPDATE allow the specification of additional tables
-- to "join" against to determine which rows should be modified.

CREATE TEMP TABLE t1 (a int, b int);
CREATE TEMP TABLE t2 (a int, b int);
CREATE TEMP TABLE t3 (x int, y int);

INSERT INTO t1 VALUES (5, 10);
INSERT INTO t1 VALUES (15, 20);
INSERT INTO t1 VALUES (100, 100);
INSERT INTO t1 VALUES (200, 1000);
INSERT INTO t2 VALUES (200, 2000);
INSERT INTO t3 VALUES (5, 20);
INSERT INTO t3 VALUES (6, 7);
INSERT INTO t3 VALUES (7, 8);
INSERT INTO t3 VALUES (500, 100);

DELETE FROM t3 USING t1 table1 WHERE t3.x = table1.a;
SELECT * FROM t3;
DELETE FROM t3 USING t1 JOIN t2 USING (a) WHERE t3.x > t1.a;
SELECT * FROM t3;
DELETE FROM t3 USING t3 t3_other WHERE t3.x = t3_other.x AND t3.y = t3_other.y;
SELECT * FROM t3;

-- Test join against inheritance tree

create temp table t2a () inherits (t2);

insert into t2a values (200, 2001);

select * from t1 left join t2 on (t1.a = t2.a);

-- Test matching of column name with wrong alias

select t1.x from t1 join t3 on (t1.a = t3.x);

-- Test matching of locking clause with wrong alias

select t1.*, t2.*, unnamed_join.* from
  t1 join t2 on (t1.a = t2.a), t3 as unnamed_join
  for update of unnamed_join;

select foo.*, unnamed_join.* from
  t1 join t2 using (a) as foo, t3 as unnamed_join
  for update of unnamed_join;

select foo.*, unnamed_join.* from
  t1 join t2 using (a) as foo, t3 as unnamed_join
  for update of foo;

select bar.*, unnamed_join.* from
  (t1 join t2 using (a) as foo) as bar, t3 as unnamed_join
  for update of foo;

select bar.*, unnamed_join.* from
  (t1 join t2 using (a) as foo) as bar, t3 as unnamed_join
  for update of bar;

--
-- regression test for 8.1 merge right join bug
--

CREATE TEMP TABLE tt1 ( tt1_id int4, joincol int4 );
INSERT INTO tt1 VALUES (1, 11);
INSERT INTO tt1 VALUES (2, NULL);

CREATE TEMP TABLE tt2 ( tt2_id int4, joincol int4 );
INSERT INTO tt2 VALUES (21, 11);
INSERT INTO tt2 VALUES (22, 11);

set enable_hashjoin to off;
set enable_nestloop to off;

-- these should give the same results

select tt1.*, tt2.* from tt1 left join tt2 on tt1.joincol = tt2.joincol;

select tt1.*, tt2.* from tt2 right join tt1 on tt1.joincol = tt2.joincol;

reset enable_hashjoin;
reset enable_nestloop;

--
-- regression test for bug #18522 (merge-right-anti-join in inner_unique cases)
--

create temp table tbl_ra(a int unique, b int);
insert into tbl_ra select i, i%100 from generate_series(1,1000)i;
create index on tbl_ra (b);
analyze tbl_ra;

set enable_hashjoin to off;
set enable_nestloop to off;

-- ensure we get a merge right anti join
explain (costs off)
select * from tbl_ra t1
where not exists (select 1 from tbl_ra t2 where t2.b = t1.a) and t1.b < 2;

-- and check we get the expected results
select * from tbl_ra t1
where not exists (select 1 from tbl_ra t2 where t2.b = t1.a) and t1.b < 2;

reset enable_hashjoin;
reset enable_nestloop;

--
-- regression test for bug #13908 (hash join with skew tuples & nbatch increase)
--

set work_mem to '64kB';
set enable_mergejoin to off;
set enable_memoize to off;

explain (costs off)
select count(*) from tenk1 a, tenk1 b
  where a.hundred = b.thousand and (b.fivethous % 10) < 10;
select count(*) from tenk1 a, tenk1 b
  where a.hundred = b.thousand and (b.fivethous % 10) < 10;

reset work_mem;
reset enable_mergejoin;
reset enable_memoize;

--
-- regression test for 8.2 bug with improper re-ordering of left joins
--

create temp table tt3(f1 int, f2 text);
insert into tt3 select x, repeat('xyzzy', 100) from generate_series(1,10000) x;
analyze tt3;

create temp table tt4(f1 int);
insert into tt4 values (0),(1),(9999);
analyze tt4;

set enable_nestloop to off;

EXPLAIN (COSTS OFF)
SELECT a.f1
FROM tt4 a
LEFT JOIN (
        SELECT b.f1
        FROM tt3 b LEFT JOIN tt3 c ON (b.f1 = c.f1)
        WHERE COALESCE(c.f1, 0) = 0
) AS d ON (a.f1 = d.f1)
WHERE COALESCE(d.f1, 0) = 0
ORDER BY 1;

SELECT a.f1
FROM tt4 a
LEFT JOIN (
        SELECT b.f1
        FROM tt3 b LEFT JOIN tt3 c ON (b.f1 = c.f1)
        WHERE COALESCE(c.f1, 0) = 0
) AS d ON (a.f1 = d.f1)
WHERE COALESCE(d.f1, 0) = 0
ORDER BY 1;

reset enable_nestloop;

--
-- basic semijoin and antijoin recognition tests
--

explain (costs off)
select a.* from tenk1 a
where unique1 in (select unique2 from tenk1 b);

-- sadly, this is not an antijoin
explain (costs off)
select a.* from tenk1 a
where unique1 not in (select unique2 from tenk1 b);

explain (costs off)
select a.* from tenk1 a
where exists (select 1 from tenk1 b where a.unique1 = b.unique2);

explain (costs off)
select a.* from tenk1 a
where not exists (select 1 from tenk1 b where a.unique1 = b.unique2);

explain (costs off)
select a.* from tenk1 a left join tenk1 b on a.unique1 = b.unique2
where b.unique2 is null;

--
-- regression test for proper handling of outer joins within antijoins
--

create temp table tt4x(c1 int, c2 int, c3 int);

explain (costs off)
select * from tt4x t1
where not exists (
  select 1 from tt4x t2
    left join tt4x t3 on t2.c3 = t3.c1
    left join ( select t5.c1 as c1
                from tt4x t4 left join tt4x t5 on t4.c2 = t5.c1
              ) a1 on t3.c2 = a1.c1
  where t1.c1 = t2.c2
);

--
-- regression test for problems of the sort depicted in bug #3494
--

create temp table tt5(f1 int, f2 int);
create temp table tt6(f1 int, f2 int);

insert into tt5 values(1, 10);
insert into tt5 values(1, 11);

insert into tt6 values(1, 9);
insert into tt6 values(1, 2);
insert into tt6 values(2, 9);

select * from tt5,tt6 where tt5.f1 = tt6.f1 and tt5.f1 = tt5.f2 - tt6.f2;

--
-- regression test for problems of the sort depicted in bug #3588
--

create temp table xx (pkxx int);
create temp table yy (pkyy int, pkxx int);

insert into xx values (1);
insert into xx values (2);
insert into xx values (3);

insert into yy values (101, 1);
insert into yy values (201, 2);
insert into yy values (301, NULL);

select yy.pkyy as yy_pkyy, yy.pkxx as yy_pkxx, yya.pkyy as yya_pkyy,
       xxa.pkxx as xxa_pkxx, xxb.pkxx as xxb_pkxx
from yy
     left join (SELECT * FROM yy where pkyy = 101) as yya ON yy.pkyy = yya.pkyy
     left join xx xxa on yya.pkxx = xxa.pkxx
     left join xx xxb on coalesce (xxa.pkxx, 1) = xxb.pkxx;

--
-- regression test for improper pushing of constants across outer-join clauses
-- (as seen in early 8.2.x releases)
--

create temp table zt1 (f1 int primary key);
create temp table zt2 (f2 int primary key);
create temp table zt3 (f3 int primary key);
insert into zt1 values(53);
insert into zt2 values(53);

select * from
  zt2 left join zt3 on (f2 = f3)
      left join zt1 on (f3 = f1)
where f2 = 53;

create temp view zv1 as select *,'dummy'::text AS junk from zt1;

select * from
  zt2 left join zt3 on (f2 = f3)
      left join zv1 on (f3 = f1)
where f2 = 53;

--
-- regression test for improper extraction of OR indexqual conditions
-- (as seen in early 8.3.x releases)
--

select a.unique2, a.ten, b.tenthous, b.unique2, b.hundred
from tenk1 a left join tenk1 b on a.unique2 = b.tenthous
where a.unique1 = 42 and
      ((b.unique2 is null and a.ten = 2) or b.hundred = 3);

--
-- test proper positioning of one-time quals in EXISTS (8.4devel bug)
--
prepare foo(bool) as
  select count(*) from tenk1 a left join tenk1 b
    on (a.unique2 = b.unique1 and exists
        (select 1 from tenk1 c where c.thousand = b.unique2 and $1));
execute foo(true);
execute foo(false);

--
-- test for sane behavior with noncanonical merge clauses, per bug #4926
--

begin;

set enable_mergejoin = 1;
set enable_hashjoin = 0;
set enable_nestloop = 0;

create temp table a (i integer);
create temp table b (x integer, y integer);

select * from a left join b on i = x and i = y and x = i;

rollback;

--
-- test handling of merge clauses using record_ops
--
begin;

create type mycomptype as (id int, v bigint);

create temp table tidv (idv mycomptype);
create index on tidv (idv);

explain (costs off)
select a.idv, b.idv from tidv a, tidv b where a.idv = b.idv;

set enable_mergejoin = 0;
set enable_hashjoin = 0;

explain (costs off)
select a.idv, b.idv from tidv a, tidv b where a.idv = b.idv;

rollback;

--
-- test NULL behavior of whole-row Vars, per bug #5025
--
select t1.q2, count(t2.*)
from int8_tbl t1 left join int8_tbl t2 on (t1.q2 = t2.q1)
group by t1.q2 order by 1;

select t1.q2, count(t2.*)
from int8_tbl t1 left join (select * from int8_tbl) t2 on (t1.q2 = t2.q1)
group by t1.q2 order by 1;

select t1.q2, count(t2.*)
from int8_tbl t1 left join (select * from int8_tbl offset 0) t2 on (t1.q2 = t2.q1)
group by t1.q2 order by 1;

select t1.q2, count(t2.*)
from int8_tbl t1 left join
  (select q1, case when q2=1 then 1 else q2 end as q2 from int8_tbl) t2
  on (t1.q2 = t2.q1)
group by t1.q2 order by 1;

--
-- test incorrect failure to NULL pulled-up subexpressions
--
begin;

create temp table a (
     code char not null,
     constraint a_pk primary key (code)
);
create temp table b (
     a char not null,
     num integer not null,
     constraint b_pk primary key (a, num)
);
create temp table c (
     name char not null,
     a char,
     constraint c_pk primary key (name)
);

insert into a (code) values ('p');
insert into a (code) values ('q');
insert into b (a, num) values ('p', 1);
insert into b (a, num) values ('p', 2);
insert into c (name, a) values ('A', 'p');
insert into c (name, a) values ('B', 'q');
insert into c (name, a) values ('C', null);

select c.name, ss.code, ss.b_cnt, ss.const
from c left join
  (select a.code, coalesce(b_grp.cnt, 0) as b_cnt, -1 as const
   from a left join
     (select count(1) as cnt, b.a from b group by b.a) as b_grp
     on a.code = b_grp.a
  ) as ss
  on (c.a = ss.code)
order by c.name;

rollback;

--
-- test incorrect handling of placeholders that only appear in targetlists,
-- per bug #6154
--
SELECT * FROM
( SELECT 1 as key1 ) sub1
LEFT JOIN
( SELECT sub3.key3, sub4.value2, COALESCE(sub4.value2, 66) as value3 FROM
    ( SELECT 1 as key3 ) sub3
    LEFT JOIN
    ( SELECT sub5.key5, COALESCE(sub6.value1, 1) as value2 FROM
        ( SELECT 1 as key5 ) sub5
        LEFT JOIN
        ( SELECT 2 as key6, 42 as value1 ) sub6
        ON sub5.key5 = sub6.key6
    ) sub4
    ON sub4.key5 = sub3.key3
) sub2
ON sub1.key1 = sub2.key3;

-- test the path using join aliases, too
SELECT * FROM
( SELECT 1 as key1 ) sub1
LEFT JOIN
( SELECT sub3.key3, value2, COALESCE(value2, 66) as value3 FROM
    ( SELECT 1 as key3 ) sub3
    LEFT JOIN
    ( SELECT sub5.key5, COALESCE(sub6.value1, 1) as value2 FROM
        ( SELECT 1 as key5 ) sub5
        LEFT JOIN
        ( SELECT 2 as key6, 42 as value1 ) sub6
        ON sub5.key5 = sub6.key6
    ) sub4
    ON sub4.key5 = sub3.key3
) sub2
ON sub1.key1 = sub2.key3;

--
-- test case where a PlaceHolderVar is used as a nestloop parameter
--

EXPLAIN (COSTS OFF)
SELECT qq, unique1
  FROM
  ( SELECT COALESCE(q1, 0) AS qq FROM int8_tbl a ) AS ss1
  FULL OUTER JOIN
  ( SELECT COALESCE(q2, -1) AS qq FROM int8_tbl b ) AS ss2
  USING (qq)
  INNER JOIN tenk1 c ON qq = unique2;

SELECT qq, unique1
  FROM
  ( SELECT COALESCE(q1, 0) AS qq FROM int8_tbl a ) AS ss1
  FULL OUTER JOIN
  ( SELECT COALESCE(q2, -1) AS qq FROM int8_tbl b ) AS ss2
  USING (qq)
  INNER JOIN tenk1 c ON qq = unique2;

--
-- nested nestloops can require nested PlaceHolderVars
--

create temp table nt1 (
  id int primary key,
  a1 boolean,
  a2 boolean
);
create temp table nt2 (
  id int primary key,
  nt1_id int,
  b1 boolean,
  b2 boolean,
  foreign key (nt1_id) references nt1(id)
);
create temp table nt3 (
  id int primary key,
  nt2_id int,
  c1 boolean,
  foreign key (nt2_id) references nt2(id)
);

insert into nt1 values (1,true,true);
insert into nt1 values (2,true,false);
insert into nt1 values (3,false,false);
insert into nt2 values (1,1,true,true);
insert into nt2 values (2,2,true,false);
insert into nt2 values (3,3,false,false);
insert into nt3 values (1,1,true);
insert into nt3 values (2,2,false);
insert into nt3 values (3,3,true);

explain (costs off)
select nt3.id
from nt3 as nt3
  left join
    (select nt2.*, (nt2.b1 and ss1.a3) AS b3
     from nt2 as nt2
       left join
         (select nt1.*, (nt1.id is not null) as a3 from nt1) as ss1
         on ss1.id = nt2.nt1_id
    ) as ss2
    on ss2.id = nt3.nt2_id
where nt3.id = 1 and ss2.b3;

select nt3.id
from nt3 as nt3
  left join
    (select nt2.*, (nt2.b1 and ss1.a3) AS b3
     from nt2 as nt2
       left join
         (select nt1.*, (nt1.id is not null) as a3 from nt1) as ss1
         on ss1.id = nt2.nt1_id
    ) as ss2
    on ss2.id = nt3.nt2_id
where nt3.id = 1 and ss2.b3;

--
-- test case where a PlaceHolderVar is propagated into a subquery
--

explain (costs off)
select * from
  int8_tbl t1 left join
  (select q1 as x, 42 as y from int8_tbl t2) ss
  on t1.q2 = ss.x
where
  1 = (select 1 from int8_tbl t3 where ss.y is not null limit 1)
order by 1,2;

select * from
  int8_tbl t1 left join
  (select q1 as x, 42 as y from int8_tbl t2) ss
  on t1.q2 = ss.x
where
  1 = (select 1 from int8_tbl t3 where ss.y is not null limit 1)
order by 1,2;

--
-- variant where a PlaceHolderVar is needed at a join, but not above the join
--

explain (costs off)
select * from
  int4_tbl as i41,
  lateral
    (select 1 as x from
      (select i41.f1 as lat,
              i42.f1 as loc from
         int8_tbl as i81, int4_tbl as i42) as ss1
      right join int4_tbl as i43 on (i43.f1 > 1)
      where ss1.loc = ss1.lat) as ss2
where i41.f1 > 0;

select * from
  int4_tbl as i41,
  lateral
    (select 1 as x from
      (select i41.f1 as lat,
              i42.f1 as loc from
         int8_tbl as i81, int4_tbl as i42) as ss1
      right join int4_tbl as i43 on (i43.f1 > 1)
      where ss1.loc = ss1.lat) as ss2
where i41.f1 > 0;

--
-- test the corner cases FULL JOIN ON TRUE and FULL JOIN ON FALSE
--
select * from int4_tbl a full join int4_tbl b on true;
select * from int4_tbl a full join int4_tbl b on false;

--
-- test for ability to use a cartesian join when necessary
--

create temp table q1 as select 1 as q1;
create temp table q2 as select 0 as q2;
analyze q1;
analyze q2;

explain (costs off)
select * from
  tenk1 join int4_tbl on f1 = twothousand,
  q1, q2
where q1 = thousand or q2 = thousand;

explain (costs off)
select * from
  tenk1 join int4_tbl on f1 = twothousand,
  q1, q2
where thousand = (q1 + q2);

--
-- test ability to generate a suitable plan for a star-schema query
--

explain (costs off)
select * from
  tenk1, int8_tbl a, int8_tbl b
where thousand = a.q1 and tenthous = b.q1 and a.q2 = 1 and b.q2 = 2;

--
-- test a corner case in which we shouldn't apply the star-schema optimization
--

explain (costs off)
select t1.unique2, t1.stringu1, t2.unique1, t2.stringu2 from
  tenk1 t1
  inner join int4_tbl i1
    left join (select v1.x2, v2.y1, 11 AS d1
               from (select 1,0 from onerow) v1(x1,x2)
               left join (select 3,1 from onerow) v2(y1,y2)
               on v1.x1 = v2.y2) subq1
    on (i1.f1 = subq1.x2)
  on (t1.unique2 = subq1.d1)
  left join tenk1 t2
  on (subq1.y1 = t2.unique1)
where t1.unique2 < 42 and t1.stringu1 > t2.stringu2;

select t1.unique2, t1.stringu1, t2.unique1, t2.stringu2 from
  tenk1 t1
  inner join int4_tbl i1
    left join (select v1.x2, v2.y1, 11 AS d1
               from (select 1,0 from onerow) v1(x1,x2)
               left join (select 3,1 from onerow) v2(y1,y2)
               on v1.x1 = v2.y2) subq1
    on (i1.f1 = subq1.x2)
  on (t1.unique2 = subq1.d1)
  left join tenk1 t2
  on (subq1.y1 = t2.unique1)
where t1.unique2 < 42 and t1.stringu1 > t2.stringu2;

-- variant that isn't quite a star-schema case

select ss1.d1 from
  tenk1 as t1
  inner join tenk1 as t2
  on t1.tenthous = t2.ten
  inner join
    int8_tbl as i8
    left join int4_tbl as i4
      inner join (select 64::information_schema.cardinal_number as d1
                  from tenk1 t3,
                       lateral (select abs(t3.unique1) + random()) ss0(x)
                  where t3.fivethous < 0) as ss1
      on i4.f1 = ss1.d1
    on i8.q1 = i4.f1
  on t1.tenthous = ss1.d1
where t1.unique1 < i4.f1;

-- this variant is foldable by the remove-useless-RESULT-RTEs code

explain (costs off)
select t1.unique2, t1.stringu1, t2.unique1, t2.stringu2 from
  tenk1 t1
  inner join int4_tbl i1
    left join (select v1.x2, v2.y1, 11 AS d1
               from (values(1,0)) v1(x1,x2)
               left join (values(3,1)) v2(y1,y2)
               on v1.x1 = v2.y2) subq1
    on (i1.f1 = subq1.x2)
  on (t1.unique2 = subq1.d1)
  left join tenk1 t2
  on (subq1.y1 = t2.unique1)
where t1.unique2 < 42 and t1.stringu1 > t2.stringu2;

select t1.unique2, t1.stringu1, t2.unique1, t2.stringu2 from
  tenk1 t1
  inner join int4_tbl i1
    left join (select v1.x2, v2.y1, 11 AS d1
               from (values(1,0)) v1(x1,x2)
               left join (values(3,1)) v2(y1,y2)
               on v1.x1 = v2.y2) subq1
    on (i1.f1 = subq1.x2)
  on (t1.unique2 = subq1.d1)
  left join tenk1 t2
  on (subq1.y1 = t2.unique1)
where t1.unique2 < 42 and t1.stringu1 > t2.stringu2;

-- Here's a variant that we can't fold too aggressively, though,
-- or we end up with noplace to evaluate the lateral PHV
explain (verbose, costs off)
select * from
  (select 1 as x) ss1 left join (select 2 as y) ss2 on (true),
  lateral (select ss2.y as z limit 1) ss3;
select * from
  (select 1 as x) ss1 left join (select 2 as y) ss2 on (true),
  lateral (select ss2.y as z limit 1) ss3;

-- Test proper handling of appendrel PHVs during useless-RTE removal
explain (costs off)
select * from
  (select 0 as z) as t1
  left join
  (select true as a) as t2
  on true,
  lateral (select true as b
           union all
           select a as b) as t3
where b;

select * from
  (select 0 as z) as t1
  left join
  (select true as a) as t2
  on true,
  lateral (select true as b
           union all
           select a as b) as t3
where b;

-- Test PHV in a semijoin qual, which confused useless-RTE removal (bug #17700)
explain (verbose, costs off)
with ctetable as not materialized ( select 1 as f1 )
select * from ctetable c1
where f1 in ( select c3.f1 from ctetable c2 full join ctetable c3 on true );

with ctetable as not materialized ( select 1 as f1 )
select * from ctetable c1
where f1 in ( select c3.f1 from ctetable c2 full join ctetable c3 on true );

-- Test PHV that winds up in a Result node, despite having nonempty nullingrels
explain (verbose, costs off)
select table_catalog, table_name
from int4_tbl t1
  inner join (int8_tbl t2
              left join information_schema.column_udt_usage on null)
  on null;

-- Test handling of qual pushdown to appendrel members with non-Var outputs
explain (verbose, costs off)
select * from int4_tbl left join (
  select text 'foo' union all select text 'bar'
) ss(x) on true
where ss.x is null;

--
-- test inlining of immutable functions
--
create function f_immutable_int4(i integer) returns integer as
$$ begin return i; end; $$ language plpgsql immutable;

-- check optimization of function scan with join
explain (costs off)
select unique1 from tenk1, (select * from f_immutable_int4(1) x) x
where x = unique1;

explain (verbose, costs off)
select unique1, x.*
from tenk1, (select *, random() from f_immutable_int4(1) x) x
where x = unique1;

explain (costs off)
select unique1 from tenk1, f_immutable_int4(1) x where x = unique1;

explain (costs off)
select unique1 from tenk1, lateral f_immutable_int4(1) x where x = unique1;

explain (costs off)
select unique1 from tenk1, lateral f_immutable_int4(1) x where x in (select 17);

explain (costs off)
select unique1, x from tenk1 join f_immutable_int4(1) x on unique1 = x;

explain (costs off)
select unique1, x from tenk1 left join f_immutable_int4(1) x on unique1 = x;

explain (costs off)
select unique1, x from tenk1 right join f_immutable_int4(1) x on unique1 = x;

explain (costs off)
select unique1, x from tenk1 full join f_immutable_int4(1) x on unique1 = x;

-- check that pullup of a const function allows further const-folding
explain (costs off)
select unique1 from tenk1, f_immutable_int4(1) x where x = 42;

-- test inlining of immutable functions with PlaceHolderVars
explain (costs off)
select nt3.id
from nt3 as nt3
  left join
    (select nt2.*, (nt2.b1 or i4 = 42) AS b3
     from nt2 as nt2
       left join
         f_immutable_int4(0) i4
         on i4 = nt2.nt1_id
    ) as ss2
    on ss2.id = nt3.nt2_id
where nt3.id = 1 and ss2.b3;

drop function f_immutable_int4(int);

-- test inlining when function returns composite

create function mki8(bigint, bigint) returns int8_tbl as
$$select row($1,$2)::int8_tbl$$ language sql;

create function mki4(int) returns int4_tbl as
$$select row($1)::int4_tbl$$ language sql;

explain (verbose, costs off)
select * from mki8(1,2);
select * from mki8(1,2);

explain (verbose, costs off)
select * from mki4(42);
select * from mki4(42);

drop function mki8(bigint, bigint);
drop function mki4(int);

-- test const-folding of a whole-row Var into a per-field Var
-- (need to inline a function to reach this case, else parser does it)
create function f_field_select(t onek) returns int4 as
$$ select t.unique2; $$ language sql immutable;

explain (verbose, costs off)
select (t2.*).unique1, f_field_select(t2) from tenk1 t1
    left join onek t2 on t1.unique1 = t2.unique1
    left join int8_tbl t3 on true;

drop function f_field_select(t onek);

--
-- test extraction of restriction OR clauses from join OR clause
-- (we used to only do this for indexable clauses)
--

explain (costs off)
select * from tenk1 a join tenk1 b on
  (a.unique1 = 1 and b.unique1 = 2) or (a.unique2 = 3 and b.hundred = 4);
explain (costs off)
select * from tenk1 a join tenk1 b on
  (a.unique1 = 1 and b.unique1 = 2) or (a.unique2 = 3 and b.ten = 4);
explain (costs off)
select * from tenk1 a join tenk1 b on
  (a.unique1 = 1 and b.unique1 = 2) or
  ((a.unique2 = 3 or a.unique2 = 7) and b.hundred = 4);

--
-- test placement of movable quals in a parameterized join tree
--

explain (costs off)
select * from tenk1 t1 left join
  (tenk1 t2 join tenk1 t3 on t2.thousand = t3.unique2)
  on t1.hundred = t2.hundred and t1.ten = t3.ten
where t1.unique1 = 1;

explain (costs off)
select * from tenk1 t1 left join
  (tenk1 t2 join tenk1 t3 on t2.thousand = t3.unique2)
  on t1.hundred = t2.hundred and t1.ten + t2.ten = t3.ten
where t1.unique1 = 1;

explain (costs off)
select count(*) from
  tenk1 a join tenk1 b on a.unique1 = b.unique2
  left join tenk1 c on a.unique2 = b.unique1 and c.thousand = a.thousand
  join int4_tbl on b.thousand = f1;

select count(*) from
  tenk1 a join tenk1 b on a.unique1 = b.unique2
  left join tenk1 c on a.unique2 = b.unique1 and c.thousand = a.thousand
  join int4_tbl on b.thousand = f1;

explain (costs off)
select b.unique1 from
  tenk1 a join tenk1 b on a.unique1 = b.unique2
  left join tenk1 c on b.unique1 = 42 and c.thousand = a.thousand
  join int4_tbl i1 on b.thousand = f1
  right join int4_tbl i2 on i2.f1 = b.tenthous
  order by 1;

select b.unique1 from
  tenk1 a join tenk1 b on a.unique1 = b.unique2
  left join tenk1 c on b.unique1 = 42 and c.thousand = a.thousand
  join int4_tbl i1 on b.thousand = f1
  right join int4_tbl i2 on i2.f1 = b.tenthous
  order by 1;

explain (costs off)
select * from
(
  select unique1, q1, coalesce(unique1, -1) + q1 as fault
  from int8_tbl left join tenk1 on (q2 = unique2)
) ss
where fault = 122
order by fault;

select * from
(
  select unique1, q1, coalesce(unique1, -1) + q1 as fault
  from int8_tbl left join tenk1 on (q2 = unique2)
) ss
where fault = 122
order by fault;

explain (costs off)
select * from
(values (1, array[10,20]), (2, array[20,30])) as v1(v1x,v1ys)
left join (values (1, 10), (2, 20)) as v2(v2x,v2y) on v2x = v1x
left join unnest(v1ys) as u1(u1y) on u1y = v2y;

select * from
(values (1, array[10,20]), (2, array[20,30])) as v1(v1x,v1ys)
left join (values (1, 10), (2, 20)) as v2(v2x,v2y) on v2x = v1x
left join unnest(v1ys) as u1(u1y) on u1y = v2y;

--
-- test handling of potential equivalence clauses above outer joins
--

explain (costs off)
select q1, unique2, thousand, hundred
  from int8_tbl a left join tenk1 b on q1 = unique2
  where coalesce(thousand,123) = q1 and q1 = coalesce(hundred,123);

select q1, unique2, thousand, hundred
  from int8_tbl a left join tenk1 b on q1 = unique2
  where coalesce(thousand,123) = q1 and q1 = coalesce(hundred,123);

explain (costs off)
select f1, unique2, case when unique2 is null then f1 else 0 end
  from int4_tbl a left join tenk1 b on f1 = unique2
  where (case when unique2 is null then f1 else 0 end) = 0;

select f1, unique2, case when unique2 is null then f1 else 0 end
  from int4_tbl a left join tenk1 b on f1 = unique2
  where (case when unique2 is null then f1 else 0 end) = 0;

--
-- another case with equivalence clauses above outer joins (bug #8591)
--

explain (costs off)
select a.unique1, b.unique1, c.unique1, coalesce(b.twothousand, a.twothousand)
  from tenk1 a left join tenk1 b on b.thousand = a.unique1                        left join tenk1 c on c.unique2 = coalesce(b.twothousand, a.twothousand)
  where a.unique2 < 10 and coalesce(b.twothousand, a.twothousand) = 44;

select a.unique1, b.unique1, c.unique1, coalesce(b.twothousand, a.twothousand)
  from tenk1 a left join tenk1 b on b.thousand = a.unique1                        left join tenk1 c on c.unique2 = coalesce(b.twothousand, a.twothousand)
  where a.unique2 < 10 and coalesce(b.twothousand, a.twothousand) = 44;

-- related case

explain (costs off)
select * from int8_tbl t1 left join int8_tbl t2 on t1.q2 = t2.q1,
  lateral (select * from int8_tbl t3 where t2.q1 = t2.q2) ss;

select * from int8_tbl t1 left join int8_tbl t2 on t1.q2 = t2.q1,
  lateral (select * from int8_tbl t3 where t2.q1 = t2.q2) ss;

--
-- check handling of join aliases when flattening multiple levels of subquery
--

explain (verbose, costs off)
select foo1.join_key as foo1_id, foo3.join_key AS foo3_id, bug_field from
  (values (0),(1)) foo1(join_key)
left join
  (select join_key, bug_field from
    (select ss1.join_key, ss1.bug_field from
      (select f1 as join_key, 666 as bug_field from int4_tbl i1) ss1
    ) foo2
   left join
    (select unique2 as join_key from tenk1 i2) ss2
   using (join_key)
  ) foo3
using (join_key);

select foo1.join_key as foo1_id, foo3.join_key AS foo3_id, bug_field from
  (values (0),(1)) foo1(join_key)
left join
  (select join_key, bug_field from
    (select ss1.join_key, ss1.bug_field from
      (select f1 as join_key, 666 as bug_field from int4_tbl i1) ss1
    ) foo2
   left join
    (select unique2 as join_key from tenk1 i2) ss2
   using (join_key)
  ) foo3
using (join_key);

--
-- check handling of a variable-free join alias
--
explain (verbose, costs off)
select * from
int4_tbl i0 left join
( (select *, 123 as x from int4_tbl i1) ss1
  left join
  (select *, q2 as x from int8_tbl i2) ss2
  using (x)
) ss0
on (i0.f1 = ss0.f1)
order by i0.f1, x;

select * from
int4_tbl i0 left join
( (select *, 123 as x from int4_tbl i1) ss1
  left join
  (select *, q2 as x from int8_tbl i2) ss2
  using (x)
) ss0
on (i0.f1 = ss0.f1)
order by i0.f1, x;

--
-- test successful handling of nested outer joins with degenerate join quals
--

explain (verbose, costs off)
select t1.* from
  text_tbl t1
  left join (select *, '***'::text as d1 from int8_tbl i8b1) b1
    left join int8_tbl i8
      left join (select *, null::int as d2 from int8_tbl i8b2) b2
      on (i8.q1 = b2.q1)
    on (b2.d2 = b1.q2)
  on (t1.f1 = b1.d1)
  left join int4_tbl i4
  on (i8.q2 = i4.f1);

select t1.* from
  text_tbl t1
  left join (select *, '***'::text as d1 from int8_tbl i8b1) b1
    left join int8_tbl i8
      left join (select *, null::int as d2 from int8_tbl i8b2) b2
      on (i8.q1 = b2.q1)
    on (b2.d2 = b1.q2)
  on (t1.f1 = b1.d1)
  left join int4_tbl i4
  on (i8.q2 = i4.f1);

explain (verbose, costs off)
select t1.* from
  text_tbl t1
  left join (select *, '***'::text as d1 from int8_tbl i8b1) b1
    left join int8_tbl i8
      left join (select *, null::int as d2 from int8_tbl i8b2, int4_tbl i4b2) b2
      on (i8.q1 = b2.q1)
    on (b2.d2 = b1.q2)
  on (t1.f1 = b1.d1)
  left join int4_tbl i4
  on (i8.q2 = i4.f1);

select t1.* from
  text_tbl t1
  left join (select *, '***'::text as d1 from int8_tbl i8b1) b1
    left join int8_tbl i8
      left join (select *, null::int as d2 from int8_tbl i8b2, int4_tbl i4b2) b2
      on (i8.q1 = b2.q1)
    on (b2.d2 = b1.q2)
  on (t1.f1 = b1.d1)
  left join int4_tbl i4
  on (i8.q2 = i4.f1);

explain (verbose, costs off)
select t1.* from
  text_tbl t1
  left join (select *, '***'::text as d1 from int8_tbl i8b1) b1
    left join int8_tbl i8
      left join (select *, null::int as d2 from int8_tbl i8b2, int4_tbl i4b2
                 where q1 = f1) b2
      on (i8.q1 = b2.q1)
    on (b2.d2 = b1.q2)
  on (t1.f1 = b1.d1)
  left join int4_tbl i4
  on (i8.q2 = i4.f1);

select t1.* from
  text_tbl t1
  left join (select *, '***'::text as d1 from int8_tbl i8b1) b1
    left join int8_tbl i8
      left join (select *, null::int as d2 from int8_tbl i8b2, int4_tbl i4b2
                 where q1 = f1) b2
      on (i8.q1 = b2.q1)
    on (b2.d2 = b1.q2)
  on (t1.f1 = b1.d1)
  left join int4_tbl i4
  on (i8.q2 = i4.f1);

explain (verbose, costs off)
select * from
  text_tbl t1
  inner join int8_tbl i8
  on i8.q2 = 456
  right join text_tbl t2
  on t1.f1 = 'doh!'
  left join int4_tbl i4
  on i8.q1 = i4.f1;

select * from
  text_tbl t1
  inner join int8_tbl i8
  on i8.q2 = 456
  right join text_tbl t2
  on t1.f1 = 'doh!'
  left join int4_tbl i4
  on i8.q1 = i4.f1;

-- check handling of a variable-free qual for a non-commutable outer join
explain (costs off)
select nspname
from (select 1 as x) ss1
left join
( select n.nspname, c.relname
  from pg_class c left join pg_namespace n on n.oid = c.relnamespace
  where c.relkind = 'r'
) ss2 on false;

-- check handling of apparently-commutable outer joins with non-commutable
-- joins between them
explain (costs off)
select 1 from
  int4_tbl i4
  left join int8_tbl i8 on i4.f1 is not null
  left join (select 1 as a) ss1 on null
  join int4_tbl i42 on ss1.a is null or i8.q1 <> i8.q2
  right join (select 2 as b) ss2
  on ss2.b < i4.f1;

--
-- test for appropriate join order in the presence of lateral references
--

explain (verbose, costs off)
select * from
  text_tbl t1
  left join int8_tbl i8
  on i8.q2 = 123,
  lateral (select i8.q1, t2.f1 from text_tbl t2 limit 1) as ss
where t1.f1 = ss.f1;

select * from
  text_tbl t1
  left join int8_tbl i8
  on i8.q2 = 123,
  lateral (select i8.q1, t2.f1 from text_tbl t2 limit 1) as ss
where t1.f1 = ss.f1;

explain (verbose, costs off)
select * from
  text_tbl t1
  left join int8_tbl i8
  on i8.q2 = 123,
  lateral (select i8.q1, t2.f1 from text_tbl t2 limit 1) as ss1,
  lateral (select ss1.* from text_tbl t3 limit 1) as ss2
where t1.f1 = ss2.f1;

select * from
  text_tbl t1
  left join int8_tbl i8
  on i8.q2 = 123,
  lateral (select i8.q1, t2.f1 from text_tbl t2 limit 1) as ss1,
  lateral (select ss1.* from text_tbl t3 limit 1) as ss2
where t1.f1 = ss2.f1;

explain (verbose, costs off)
select 1 from
  text_tbl as tt1
  inner join text_tbl as tt2 on (tt1.f1 = 'foo')
  left join text_tbl as tt3 on (tt3.f1 = 'foo')
  left join text_tbl as tt4 on (tt3.f1 = tt4.f1),
  lateral (select tt4.f1 as c0 from text_tbl as tt5 limit 1) as ss1
where tt1.f1 = ss1.c0;

select 1 from
  text_tbl as tt1
  inner join text_tbl as tt2 on (tt1.f1 = 'foo')
  left join text_tbl as tt3 on (tt3.f1 = 'foo')
  left join text_tbl as tt4 on (tt3.f1 = tt4.f1),
  lateral (select tt4.f1 as c0 from text_tbl as tt5 limit 1) as ss1
where tt1.f1 = ss1.c0;

explain (verbose, costs off)
select 1 from
  int4_tbl as i4
  inner join
    ((select 42 as n from int4_tbl x1 left join int8_tbl x2 on f1 = q1) as ss1
     right join (select 1 as z) as ss2 on true)
  on false,
  lateral (select i4.f1, ss1.n from int8_tbl as i8 limit 1) as ss3;

select 1 from
  int4_tbl as i4
  inner join
    ((select 42 as n from int4_tbl x1 left join int8_tbl x2 on f1 = q1) as ss1
     right join (select 1 as z) as ss2 on true)
  on false,
  lateral (select i4.f1, ss1.n from int8_tbl as i8 limit 1) as ss3;

--
-- check a case where we formerly generated invalid parameterized paths
--

begin;

create temp table t (a int unique);

explain (costs off)
select 1 from t t1
  join lateral (select t1.a from (select 1) foo offset 0) as s1 on true
  join
    (select 1 from t t2
       inner join (t t3
                   left join (t t4 left join t t5 on t4.a = 1)
                   on t3.a = t4.a)
       on false
     where t3.a = coalesce(t5.a,1)) as s2
  on true;

rollback;

--
-- check a case in which a PlaceHolderVar forces join order
--

explain (verbose, costs off)
select ss2.* from
  int4_tbl i41
  left join int8_tbl i8
    join (select i42.f1 as c1, i43.f1 as c2, 42 as c3
          from int4_tbl i42, int4_tbl i43) ss1
    on i8.q1 = ss1.c2
  on i41.f1 = ss1.c1,
  lateral (select i41.*, i8.*, ss1.* from text_tbl limit 1) ss2
where ss1.c2 = 0;

select ss2.* from
  int4_tbl i41
  left join int8_tbl i8
    join (select i42.f1 as c1, i43.f1 as c2, 42 as c3
          from int4_tbl i42, int4_tbl i43) ss1
    on i8.q1 = ss1.c2
  on i41.f1 = ss1.c1,
  lateral (select i41.*, i8.*, ss1.* from text_tbl limit 1) ss2
where ss1.c2 = 0;

--
-- test successful handling of full join underneath left join (bug #14105)
--

explain (costs off)
select * from
  (select 1 as id) as xx
  left join
    (tenk1 as a1 full join (select 1 as id) as yy on (a1.unique1 = yy.id))
  on (xx.id = coalesce(yy.id));

select * from
  (select 1 as id) as xx
  left join
    (tenk1 as a1 full join (select 1 as id) as yy on (a1.unique1 = yy.id))
  on (xx.id = coalesce(yy.id));

--
-- test ability to push constants through outer join clauses
--

explain (costs off)
  select * from int4_tbl a left join tenk1 b on f1 = unique2 where f1 = 0;

explain (costs off)
  select * from tenk1 a full join tenk1 b using(unique2) where unique2 = 42;

--
-- test that quals attached to an outer join have correct semantics,
-- specifically that they don't re-use expressions computed below the join;
-- we force a mergejoin so that coalesce(b.q1, 1) appears as a join input
--

set enable_hashjoin to off;
set enable_nestloop to off;

explain (verbose, costs off)
  select a.q2, b.q1
    from int8_tbl a left join int8_tbl b on a.q2 = coalesce(b.q1, 1)
    where coalesce(b.q1, 1) > 0;
select a.q2, b.q1
  from int8_tbl a left join int8_tbl b on a.q2 = coalesce(b.q1, 1)
  where coalesce(b.q1, 1) > 0;

reset enable_hashjoin;
reset enable_nestloop;

--
-- test join strength reduction with a SubPlan providing the proof
--

explain (costs off)
select a.unique1, b.unique2
  from onek a left join onek b on a.unique1 = b.unique2
  where (b.unique2, random() > 0) = any (select q1, random() > 0 from int8_tbl c where c.q1 < b.unique1);

select a.unique1, b.unique2
  from onek a left join onek b on a.unique1 = b.unique2
  where (b.unique2, random() > 0) = any (select q1, random() > 0 from int8_tbl c where c.q1 < b.unique1);

--
-- test full-join strength reduction
--

explain (costs off)
select a.unique1, b.unique2
  from onek a full join onek b on a.unique1 = b.unique2
  where a.unique1 = 42;

select a.unique1, b.unique2
  from onek a full join onek b on a.unique1 = b.unique2
  where a.unique1 = 42;

explain (costs off)
select a.unique1, b.unique2
  from onek a full join onek b on a.unique1 = b.unique2
  where b.unique2 = 43;

select a.unique1, b.unique2
  from onek a full join onek b on a.unique1 = b.unique2
  where b.unique2 = 43;

explain (costs off)
select a.unique1, b.unique2
  from onek a full join onek b on a.unique1 = b.unique2
  where a.unique1 = 42 and b.unique2 = 42;

select a.unique1, b.unique2
  from onek a full join onek b on a.unique1 = b.unique2
  where a.unique1 = 42 and b.unique2 = 42;

--
-- test result-RTE removal underneath a full join
--

explain (costs off)
select * from
  (select * from int8_tbl i81 join (values(123,2)) v(v1,v2) on q2=v1) ss1
full join
  (select * from (values(456,2)) w(v1,v2) join int8_tbl i82 on q2=v1) ss2
on true;

select * from
  (select * from int8_tbl i81 join (values(123,2)) v(v1,v2) on q2=v1) ss1
full join
  (select * from (values(456,2)) w(v1,v2) join int8_tbl i82 on q2=v1) ss2
on true;

--
-- test join removal
--

begin;

CREATE TEMP TABLE a (id int PRIMARY KEY, b_id int);
CREATE TEMP TABLE b (id int PRIMARY KEY, c_id int);
CREATE TEMP TABLE c (id int PRIMARY KEY);
CREATE TEMP TABLE d (a int, b int);
INSERT INTO a VALUES (0, 0), (1, NULL);
INSERT INTO b VALUES (0, 0), (1, NULL);
INSERT INTO c VALUES (0), (1);
INSERT INTO d VALUES (1,3), (2,2), (3,1);

-- all three cases should be optimizable into a simple seqscan
explain (costs off) SELECT a.* FROM a LEFT JOIN b ON a.b_id = b.id;
explain (costs off) SELECT b.* FROM b LEFT JOIN c ON b.c_id = c.id;
explain (costs off)
  SELECT a.* FROM a LEFT JOIN (b left join c on b.c_id = c.id)
  ON (a.b_id = b.id);

-- check optimization of outer join within another special join
explain (costs off)
select id from a where id in (
	select b.id from b left join c on b.id = c.id
);

-- check optimization with oddly-nested outer joins
explain (costs off)
select a1.id from
  (a a1 left join a a2 on true)
  left join
  (a a3 left join a a4 on a3.id = a4.id)
  on a2.id = a3.id;

explain (costs off)
select a1.id from
  (a a1 left join a a2 on a1.id = a2.id)
  left join
  (a a3 left join a a4 on a3.id = a4.id)
  on a2.id = a3.id;

explain (costs off)
select 1 from a t1
    left join a t2 on true
   inner join a t3 on true
    left join a t4 on t2.id = t4.id and t2.id = t3.id;

-- another example (bug #17781)
explain (costs off)
select ss1.f1
from int4_tbl as t1
  left join (int4_tbl as t2
             right join int4_tbl as t3 on null
             left join (int4_tbl as t4
                        right join int8_tbl as t5 on null)
               on t2.f1 = t4.f1
             left join ((select null as f1 from int4_tbl as t6) as ss1
                        inner join int8_tbl as t7 on null)
               on t5.q1 = t7.q2)
    on false;

-- variant with Var rather than PHV coming from t6
explain (costs off)
select ss1.f1
from int4_tbl as t1
  left join (int4_tbl as t2
             right join int4_tbl as t3 on null
             left join (int4_tbl as t4
                        right join int8_tbl as t5 on null)
               on t2.f1 = t4.f1
             left join ((select f1 from int4_tbl as t6) as ss1
                        inner join int8_tbl as t7 on null)
               on t5.q1 = t7.q2)
    on false;

-- per further discussion of bug #17781
explain (costs off)
select ss1.x
from (select f1/2 as x from int4_tbl i4 left join a on a.id = i4.f1) ss1
     right join int8_tbl i8 on true
where current_user is not null;  -- this is to add a Result node

-- and further discussion of bug #17781
explain (costs off)
select *
from int8_tbl t1
  left join (int8_tbl t2 left join onek t3 on t2.q1 > t3.unique1)
    on t1.q2 = t2.q2
  left join onek t4
    on t2.q2 < t3.unique2;

-- More tests of correct placement of pseudoconstant quals

-- simple constant-false condition
explain (costs off)
select * from int8_tbl t1 left join
  (int8_tbl t2 inner join int8_tbl t3 on false
   left join int8_tbl t4 on t2.q2 = t4.q2)
on t1.q1 = t2.q1;

-- deduce constant-false from an EquivalenceClass
explain (costs off)
select * from int8_tbl t1 left join
  (int8_tbl t2 inner join int8_tbl t3 on (t2.q1-t3.q2) = 0 and (t2.q1-t3.q2) = 1
   left join int8_tbl t4 on t2.q2 = t4.q2)
on t1.q1 = t2.q1;

-- pseudoconstant based on an outer-level Param
explain (costs off)
select exists(
  select * from int8_tbl t1 left join
    (int8_tbl t2 inner join int8_tbl t3 on x0.f1 = 1
     left join int8_tbl t4 on t2.q2 = t4.q2)
  on t1.q1 = t2.q1
) from int4_tbl x0;

-- check that join removal works for a left join when joining a subquery
-- that is guaranteed to be unique by its GROUP BY clause
explain (costs off)
select d.* from d left join (select * from b group by b.id, b.c_id) s
  on d.a = s.id and d.b = s.c_id;

-- similarly, but keying off a DISTINCT clause
explain (costs off)
select d.* from d left join (select distinct * from b) s
  on d.a = s.id and d.b = s.c_id;

-- join removal is not possible when the GROUP BY contains a column that is
-- not in the join condition.  (Note: as of 9.6, we notice that b.id is a
-- primary key and so drop b.c_id from the GROUP BY of the resulting plan;
-- but this happens too late for join removal in the outer plan level.)
explain (costs off)
select d.* from d left join (select * from b group by b.id, b.c_id) s
  on d.a = s.id;

-- similarly, but keying off a DISTINCT clause
explain (costs off)
select d.* from d left join (select distinct * from b) s
  on d.a = s.id;

-- join removal is not possible here
explain (costs off)
select 1 from a t1
  left join (a t2 left join a t3 on t2.id = 1) on t2.id = 1;

-- check join removal works when uniqueness of the join condition is enforced
-- by a UNION
explain (costs off)
select d.* from d left join (select id from a union select id from b) s
  on d.a = s.id;

-- check join removal with a cross-type comparison operator
explain (costs off)
select i8.* from int8_tbl i8 left join (select f1 from int4_tbl group by f1) i4
  on i8.q1 = i4.f1;

-- check join removal with lateral references
explain (costs off)
select 1 from (select a.id FROM a left join b on a.b_id = b.id) q,
			  lateral generate_series(1, q.id) gs(i) where q.id = gs.i;

-- check join removal within RHS of an outer join
explain (costs off)
select c.id, ss.a from c
  left join (select d.a from onerow, d left join b on d.a = b.id) ss
  on c.id = ss.a;

CREATE TEMP TABLE parted_b (id int PRIMARY KEY) partition by range(id);
CREATE TEMP TABLE parted_b1 partition of parted_b for values from (0) to (10);

-- test join removals on a partitioned table
explain (costs off)
select a.* from a left join parted_b pb on a.b_id = pb.id;

rollback;

create temp table parent (k int primary key, pd int);
create temp table child (k int unique, cd int);
insert into parent values (1, 10), (2, 20), (3, 30);
insert into child values (1, 100), (4, 400);

-- this case is optimizable
select p.* from parent p left join child c on (p.k = c.k);
explain (costs off)
  select p.* from parent p left join child c on (p.k = c.k);

-- this case is not
select p.*, linked from parent p
  left join (select c.*, true as linked from child c) as ss
  on (p.k = ss.k);
explain (costs off)
  select p.*, linked from parent p
    left join (select c.*, true as linked from child c) as ss
    on (p.k = ss.k);

-- check for a 9.0rc1 bug: join removal breaks pseudoconstant qual handling
select p.* from
  parent p left join child c on (p.k = c.k)
  where p.k = 1 and p.k = 2;
explain (costs off)
select p.* from
  parent p left join child c on (p.k = c.k)
  where p.k = 1 and p.k = 2;

select p.* from
  (parent p left join child c on (p.k = c.k)) join parent x on p.k = x.k
  where p.k = 1 and p.k = 2;
explain (costs off)
select p.* from
  (parent p left join child c on (p.k = c.k)) join parent x on p.k = x.k
  where p.k = 1 and p.k = 2;

-- bug 5255: this is not optimizable by join removal
begin;

CREATE TEMP TABLE a (id int PRIMARY KEY);
CREATE TEMP TABLE b (id int PRIMARY KEY, a_id int);
INSERT INTO a VALUES (0), (1);
INSERT INTO b VALUES (0, 0), (1, NULL);

SELECT * FROM b LEFT JOIN a ON (b.a_id = a.id) WHERE (a.id IS NULL OR a.id > 0);
SELECT b.* FROM b LEFT JOIN a ON (b.a_id = a.id) WHERE (a.id IS NULL OR a.id > 0);

rollback;

-- another join removal bug: this is not optimizable, either
begin;

create temp table innertab (id int8 primary key, dat1 int8);
insert into innertab values(123, 42);

SELECT * FROM
    (SELECT 1 AS x) ss1
  LEFT JOIN
    (SELECT q1, q2, COALESCE(dat1, q1) AS y
     FROM int8_tbl LEFT JOIN innertab ON q2 = id) ss2
  ON true;

-- join removal bug #17769: can't remove if there's a pushed-down reference
EXPLAIN (COSTS OFF)
SELECT q2 FROM
  (SELECT *
   FROM int8_tbl LEFT JOIN innertab ON q2 = id) ss
 WHERE COALESCE(dat1, 0) = q1;

-- join removal bug #17773: otherwise-removable PHV appears in a qual condition
EXPLAIN (VERBOSE, COSTS OFF)
SELECT q2 FROM
  (SELECT q2, 'constant'::text AS x
   FROM int8_tbl LEFT JOIN innertab ON q2 = id) ss
  RIGHT JOIN int4_tbl ON NULL
 WHERE x >= x;

-- join removal bug #17786: check that OR conditions are cleaned up
EXPLAIN (COSTS OFF)
SELECT f1, x
FROM int4_tbl
     JOIN ((SELECT 42 AS x FROM int8_tbl LEFT JOIN innertab ON q1 = id) AS ss1
           RIGHT JOIN tenk1 ON NULL)
        ON tenk1.unique1 = ss1.x OR tenk1.unique2 = ss1.x;

rollback;

-- another join removal bug: we must clean up correctly when removing a PHV
begin;

create temp table uniquetbl (f1 text unique);

explain (costs off)
select t1.* from
  uniquetbl as t1
  left join (select *, '***'::text as d1 from uniquetbl) t2
  on t1.f1 = t2.f1
  left join uniquetbl t3
  on t2.d1 = t3.f1;

explain (costs off)
select t0.*
from
 text_tbl t0
 left join
   (select case t1.ten when 0 then 'doh!'::text else null::text end as case1,
           t1.stringu2
     from tenk1 t1
     join int4_tbl i4 ON i4.f1 = t1.unique2
     left join uniquetbl u1 ON u1.f1 = t1.string4) ss
  on t0.f1 = ss.case1
where ss.stringu2 !~* ss.case1;

select t0.*
from
 text_tbl t0
 left join
   (select case t1.ten when 0 then 'doh!'::text else null::text end as case1,
           t1.stringu2
     from tenk1 t1
     join int4_tbl i4 ON i4.f1 = t1.unique2
     left join uniquetbl u1 ON u1.f1 = t1.string4) ss
  on t0.f1 = ss.case1
where ss.stringu2 !~* ss.case1;

rollback;

-- another join removal bug: we must clean up EquivalenceClasses too
begin;

create temp table t (a int unique);
insert into t values (1);

explain (costs off)
select 1
from t t1
  left join (select 2 as c
             from t t2 left join t t3 on t2.a = t3.a) s
    on true
where t1.a = s.c;

select 1
from t t1
  left join (select 2 as c
             from t t2 left join t t3 on t2.a = t3.a) s
    on true
where t1.a = s.c;

rollback;

-- test cases where we can remove a join, but not a PHV computed at it
begin;

create temp table t (a int unique, b int);
insert into t values (1,1), (2,2);

explain (costs off)
select 1
from t t1
  left join (select t2.a, 1 as c
             from t t2 left join t t3 on t2.a = t3.a) s
  on true
  left join t t4 on true
where s.a < s.c;

explain (costs off)
select t1.a, s.*
from t t1
  left join lateral (select t2.a, coalesce(t1.a, 1) as c
                     from t t2 left join t t3 on t2.a = t3.a) s
  on true
  left join t t4 on true
where s.a < s.c;

select t1.a, s.*
from t t1
  left join lateral (select t2.a, coalesce(t1.a, 1) as c
                     from t t2 left join t t3 on t2.a = t3.a) s
  on true
  left join t t4 on true
where s.a < s.c;

rollback;

-- test case to expose miscomputation of required relid set for a PHV
explain (verbose, costs off)
select i8.*, ss.v, t.unique2
  from int8_tbl i8
    left join int4_tbl i4 on i4.f1 = 1
    left join lateral (select i4.f1 + 1 as v) as ss on true
    left join tenk1 t on t.unique2 = ss.v
where q2 = 456;

select i8.*, ss.v, t.unique2
  from int8_tbl i8
    left join int4_tbl i4 on i4.f1 = 1
    left join lateral (select i4.f1 + 1 as v) as ss on true
    left join tenk1 t on t.unique2 = ss.v
where q2 = 456;

-- and check a related issue where we miscompute required relids for
-- a PHV that's been translated to a child rel
create temp table parttbl (a integer primary key) partition by range (a);
create temp table parttbl1 partition of parttbl for values from (1) to (100);
insert into parttbl values (11), (12);
explain (costs off)
select * from
  (select *, 12 as phv from parttbl) as ss
  right join int4_tbl on true
where ss.a = ss.phv and f1 = 0;

select * from
  (select *, 12 as phv from parttbl) as ss
  right join int4_tbl on true
where ss.a = ss.phv and f1 = 0;

-- bug #8444: we've historically allowed duplicate aliases within aliased JOINs

select * from
  int8_tbl x join (int4_tbl x cross join int4_tbl y) j on q1 = f1; -- error
select * from
  int8_tbl x join (int4_tbl x cross join int4_tbl y) j on q1 = y.f1; -- error
select * from
  int8_tbl x join (int4_tbl x cross join int4_tbl y(ff)) j on q1 = f1; -- ok

--
-- Test hints given on incorrect column references are useful
--

select t1.uunique1 from
  tenk1 t1 join tenk2 t2 on t1.two = t2.two; -- error, prefer "t1" suggestion
select t2.uunique1 from
  tenk1 t1 join tenk2 t2 on t1.two = t2.two; -- error, prefer "t2" suggestion
select uunique1 from
  tenk1 t1 join tenk2 t2 on t1.two = t2.two; -- error, suggest both at once
select ctid from
  tenk1 t1 join tenk2 t2 on t1.two = t2.two; -- error, need qualification

--
-- Take care to reference the correct RTE
--

select atts.relid::regclass, s.* from pg_stats s join
    pg_attribute a on s.attname = a.attname and s.tablename =
    a.attrelid::regclass::text join (select unnest(indkey) attnum,
    indexrelid from pg_index i) atts on atts.attnum = a.attnum where
    schemaname != 'pg_catalog';

-- Test bug in rangetable flattening
explain (verbose, costs off)
select 1 from
  (select * from int8_tbl where q1 <> (select 42) offset 0) ss
where false;

--
-- Test LATERAL
--

select unique2, x.*
from tenk1 a, lateral (select * from int4_tbl b where f1 = a.unique1) x;
explain (costs off)
  select unique2, x.*
  from tenk1 a, lateral (select * from int4_tbl b where f1 = a.unique1) x;
select unique2, x.*
from int4_tbl x, lateral (select unique2 from tenk1 where f1 = unique1) ss;
explain (costs off)
  select unique2, x.*
  from int4_tbl x, lateral (select unique2 from tenk1 where f1 = unique1) ss;
explain (costs off)
  select unique2, x.*
  from int4_tbl x cross join lateral (select unique2 from tenk1 where f1 = unique1) ss;
select unique2, x.*
from int4_tbl x left join lateral (select unique1, unique2 from tenk1 where f1 = unique1) ss on true;
explain (costs off)
  select unique2, x.*
  from int4_tbl x left join lateral (select unique1, unique2 from tenk1 where f1 = unique1) ss on true;

-- check scoping of lateral versus parent references
-- the first of these should return int8_tbl.q2, the second int8_tbl.q1
select *, (select r from (select q1 as q2) x, (select q2 as r) y) from int8_tbl;
select *, (select r from (select q1 as q2) x, lateral (select q2 as r) y) from int8_tbl;

-- lateral with function in FROM
select count(*) from tenk1 a, lateral generate_series(1,two) g;
explain (costs off)
  select count(*) from tenk1 a, lateral generate_series(1,two) g;
explain (costs off)
  select count(*) from tenk1 a cross join lateral generate_series(1,two) g;
-- don't need the explicit LATERAL keyword for functions
explain (costs off)
  select count(*) from tenk1 a, generate_series(1,two) g;

-- lateral with UNION ALL subselect
explain (costs off)
  select * from generate_series(100,200) g,
    lateral (select * from int8_tbl a where g = q1 union all
             select * from int8_tbl b where g = q2) ss;
select * from generate_series(100,200) g,
  lateral (select * from int8_tbl a where g = q1 union all
           select * from int8_tbl b where g = q2) ss;

-- lateral with VALUES
explain (costs off)
  select count(*) from tenk1 a,
    tenk1 b join lateral (values(a.unique1)) ss(x) on b.unique2 = ss.x;
select count(*) from tenk1 a,
  tenk1 b join lateral (values(a.unique1)) ss(x) on b.unique2 = ss.x;

-- lateral with VALUES, no flattening possible
explain (costs off)
  select count(*) from tenk1 a,
    tenk1 b join lateral (values(a.unique1),(-1)) ss(x) on b.unique2 = ss.x;
select count(*) from tenk1 a,
  tenk1 b join lateral (values(a.unique1),(-1)) ss(x) on b.unique2 = ss.x;

-- lateral injecting a strange outer join condition
explain (costs off)
  select * from int8_tbl a,
    int8_tbl x left join lateral (select a.q1 from int4_tbl y) ss(z)
      on x.q2 = ss.z
  order by a.q1, a.q2, x.q1, x.q2, ss.z;
select * from int8_tbl a,
  int8_tbl x left join lateral (select a.q1 from int4_tbl y) ss(z)
    on x.q2 = ss.z
  order by a.q1, a.q2, x.q1, x.q2, ss.z;

-- lateral reference to a join alias variable
select * from (select f1/2 as x from int4_tbl) ss1 join int4_tbl i4 on x = f1,
  lateral (select x) ss2(y);
select * from (select f1 as x from int4_tbl) ss1 join int4_tbl i4 on x = f1,
  lateral (values(x)) ss2(y);
select * from ((select f1/2 as x from int4_tbl) ss1 join int4_tbl i4 on x = f1) j,
  lateral (select x) ss2(y);

-- lateral references requiring pullup
select * from (values(1)) x(lb),
  lateral generate_series(lb,4) x4;
select * from (select f1/1000000000 from int4_tbl) x(lb),
  lateral generate_series(lb,4) x4;
select * from (values(1)) x(lb),
  lateral (values(lb)) y(lbcopy);
select * from (values(1)) x(lb),
  lateral (select lb from int4_tbl) y(lbcopy);
select * from
  int8_tbl x left join (select q1,coalesce(q2,0) q2 from int8_tbl) y on x.q2 = y.q1,
  lateral (values(x.q1,y.q1,y.q2)) v(xq1,yq1,yq2);
select * from
  int8_tbl x left join (select q1,coalesce(q2,0) q2 from int8_tbl) y on x.q2 = y.q1,
  lateral (select x.q1,y.q1,y.q2) v(xq1,yq1,yq2);
select x.* from
  int8_tbl x left join (select q1,coalesce(q2,0) q2 from int8_tbl) y on x.q2 = y.q1,
  lateral (select x.q1,y.q1,y.q2) v(xq1,yq1,yq2);
select v.* from
  (int8_tbl x left join (select q1,coalesce(q2,0) q2 from int8_tbl) y on x.q2 = y.q1)
  left join int4_tbl z on z.f1 = x.q2,
  lateral (select x.q1,y.q1 union all select x.q2,y.q2) v(vx,vy);
select v.* from
  (int8_tbl x left join (select q1,(select coalesce(q2,0)) q2 from int8_tbl) y on x.q2 = y.q1)
  left join int4_tbl z on z.f1 = x.q2,
  lateral (select x.q1,y.q1 union all select x.q2,y.q2) v(vx,vy);
select v.* from
  (int8_tbl x left join (select q1,(select coalesce(q2,0)) q2 from int8_tbl) y on x.q2 = y.q1)
  left join int4_tbl z on z.f1 = x.q2,
  lateral (select x.q1,y.q1 from onerow union all select x.q2,y.q2 from onerow) v(vx,vy);

explain (verbose, costs off)
select * from
  int8_tbl a left join
  lateral (select *, a.q2 as x from int8_tbl b) ss on a.q2 = ss.q1;
select * from
  int8_tbl a left join
  lateral (select *, a.q2 as x from int8_tbl b) ss on a.q2 = ss.q1;
explain (verbose, costs off)
select * from
  int8_tbl a left join
  lateral (select *, coalesce(a.q2, 42) as x from int8_tbl b) ss on a.q2 = ss.q1;
select * from
  int8_tbl a left join
  lateral (select *, coalesce(a.q2, 42) as x from int8_tbl b) ss on a.q2 = ss.q1;

-- lateral can result in join conditions appearing below their
-- real semantic level
explain (verbose, costs off)
select * from int4_tbl i left join
  lateral (select * from int2_tbl j where i.f1 = j.f1) k on true;
select * from int4_tbl i left join
  lateral (select * from int2_tbl j where i.f1 = j.f1) k on true;
explain (verbose, costs off)
select * from int4_tbl i left join
  lateral (select coalesce(i) from int2_tbl j where i.f1 = j.f1) k on true;
select * from int4_tbl i left join
  lateral (select coalesce(i) from int2_tbl j where i.f1 = j.f1) k on true;
explain (verbose, costs off)
select * from int4_tbl a,
  lateral (
    select * from int4_tbl b left join int8_tbl c on (b.f1 = q1 and a.f1 = q2)
  ) ss;
select * from int4_tbl a,
  lateral (
    select * from int4_tbl b left join int8_tbl c on (b.f1 = q1 and a.f1 = q2)
  ) ss;

-- lateral reference in a PlaceHolderVar evaluated at join level
explain (verbose, costs off)
select * from
  int8_tbl a left join lateral
  (select b.q1 as bq1, c.q1 as cq1, least(a.q1,b.q1,c.q1) from
   int8_tbl b cross join int8_tbl c) ss
  on a.q2 = ss.bq1;
select * from
  int8_tbl a left join lateral
  (select b.q1 as bq1, c.q1 as cq1, least(a.q1,b.q1,c.q1) from
   int8_tbl b cross join int8_tbl c) ss
  on a.q2 = ss.bq1;

-- case requiring nested PlaceHolderVars
explain (verbose, costs off)
select * from
  int8_tbl c left join (
    int8_tbl a left join (select q1, coalesce(q2,42) as x from int8_tbl b) ss1
      on a.q2 = ss1.q1
    cross join
    lateral (select q1, coalesce(ss1.x,q2) as y from int8_tbl d) ss2
  ) on c.q2 = ss2.q1,
  lateral (select ss2.y offset 0) ss3;

-- another case requiring nested PlaceHolderVars
explain (verbose, costs off)
select * from
  (select 0 as val0) as ss0
  left join (select 1 as val) as ss1 on true
  left join lateral (select ss1.val as val_filtered where false) as ss2 on true;

select * from
  (select 0 as val0) as ss0
  left join (select 1 as val) as ss1 on true
  left join lateral (select ss1.val as val_filtered where false) as ss2 on true;

-- case that breaks the old ph_may_need optimization
explain (verbose, costs off)
select c.*,a.*,ss1.q1,ss2.q1,ss3.* from
  int8_tbl c left join (
    int8_tbl a left join
      (select q1, coalesce(q2,f1) as x from int8_tbl b, int4_tbl b2
       where q1 < f1) ss1
      on a.q2 = ss1.q1
    cross join
    lateral (select q1, coalesce(ss1.x,q2) as y from int8_tbl d) ss2
  ) on c.q2 = ss2.q1,
  lateral (select * from int4_tbl i where ss2.y > f1) ss3;

-- check processing of postponed quals (bug #9041)
explain (verbose, costs off)
select * from
  (select 1 as x offset 0) x cross join (select 2 as y offset 0) y
  left join lateral (
    select * from (select 3 as z offset 0) z where z.z = x.x
  ) zz on zz.z = y.y;

-- a new postponed-quals issue (bug #17768)
explain (costs off)
select * from int4_tbl t1,
  lateral (select * from int4_tbl t2 inner join int4_tbl t3 on t1.f1 = 1
           inner join (int4_tbl t4 left join int4_tbl t5 on true) on true) ss;

-- check dummy rels with lateral references (bug #15694)
explain (verbose, costs off)
select * from int8_tbl i8 left join lateral
  (select *, i8.q2 from int4_tbl where false) ss on true;
explain (verbose, costs off)
select * from int8_tbl i8 left join lateral
  (select *, i8.q2 from int4_tbl i1, int4_tbl i2 where false) ss on true;

-- check handling of nested appendrels inside LATERAL
select * from
  ((select 2 as v) union all (select 3 as v)) as q1
  cross join lateral
  ((select * from
      ((select 4 as v) union all (select 5 as v)) as q3)
   union all
   (select q1.v)
  ) as q2;

-- check the number of columns specified
SELECT * FROM (int8_tbl i cross join int4_tbl j) ss(a,b,c,d);

-- check we don't try to do a unique-ified semijoin with LATERAL
explain (verbose, costs off)
select * from
  (values (0,9998), (1,1000)) v(id,x),
  lateral (select f1 from int4_tbl
           where f1 = any (select unique1 from tenk1
                           where unique2 = v.x offset 0)) ss;
select * from
  (values (0,9998), (1,1000)) v(id,x),
  lateral (select f1 from int4_tbl
           where f1 = any (select unique1 from tenk1
                           where unique2 = v.x offset 0)) ss;

-- check proper extParam/allParam handling (this isn't exactly a LATERAL issue,
-- but we can make the test case much more compact with LATERAL)
explain (verbose, costs off)
select * from (values (0), (1)) v(id),
lateral (select * from int8_tbl t1,
         lateral (select * from
                    (select * from int8_tbl t2
                     where (q1, random() > 0) = any (select q2, random() > 0 from int8_tbl t3
                                     where q2 = (select greatest(t1.q1,t2.q2))
                                       and (select v.id=0)) offset 0) ss2) ss
         where t1.q1 = ss.q2) ss0;

select * from (values (0), (1)) v(id),
lateral (select * from int8_tbl t1,
         lateral (select * from
                    (select * from int8_tbl t2
                     where (q1, random() > 0) = any (select q2, random() > 0 from int8_tbl t3
                                     where q2 = (select greatest(t1.q1,t2.q2))
                                       and (select v.id=0)) offset 0) ss2) ss
         where t1.q1 = ss.q2) ss0;

-- test some error cases where LATERAL should have been used but wasn't
select f1,g from int4_tbl a, (select f1 as g) ss;
select f1,g from int4_tbl a, (select a.f1 as g) ss;
select f1,g from int4_tbl a cross join (select f1 as g) ss;
select f1,g from int4_tbl a cross join (select a.f1 as g) ss;
-- SQL:2008 says the left table is in scope but illegal to access here
select f1,g from int4_tbl a right join lateral generate_series(0, a.f1) g on true;
select f1,g from int4_tbl a full join lateral generate_series(0, a.f1) g on true;
-- check we complain about ambiguous table references
select * from
  int8_tbl x cross join (int4_tbl x cross join lateral (select x.f1) ss);
-- LATERAL can be used to put an aggregate into the FROM clause of its query
select 1 from tenk1 a, lateral (select max(a.unique1) from int4_tbl b) ss;

-- check behavior of LATERAL in UPDATE/DELETE

create temp table xx1 as select f1 as x1, -f1 as x2 from int4_tbl;

-- error, can't do this:
update xx1 set x2 = f1 from (select * from int4_tbl where f1 = x1) ss;
update xx1 set x2 = f1 from (select * from int4_tbl where f1 = xx1.x1) ss;
-- can't do it even with LATERAL:
update xx1 set x2 = f1 from lateral (select * from int4_tbl where f1 = x1) ss;
-- we might in future allow something like this, but for now it's an error:
update xx1 set x2 = f1 from xx1, lateral (select * from int4_tbl where f1 = x1) ss;

-- also errors:
delete from xx1 using (select * from int4_tbl where f1 = x1) ss;
delete from xx1 using (select * from int4_tbl where f1 = xx1.x1) ss;
delete from xx1 using lateral (select * from int4_tbl where f1 = x1) ss;

--
-- test LATERAL reference propagation down a multi-level inheritance hierarchy
-- produced for a multi-level partitioned table hierarchy.
--
create table join_pt1 (a int, b int, c varchar) partition by range(a);
create table join_pt1p1 partition of join_pt1 for values from (0) to (100) partition by range(b);
create table join_pt1p2 partition of join_pt1 for values from (100) to (200);
create table join_pt1p1p1 partition of join_pt1p1 for values from (0) to (100);
insert into join_pt1 values (1, 1, 'x'), (101, 101, 'y');
create table join_ut1 (a int, b int, c varchar);
insert into join_ut1 values (101, 101, 'y'), (2, 2, 'z');
explain (verbose, costs off)
select t1.b, ss.phv from join_ut1 t1 left join lateral
              (select t2.a as t2a, t3.a t3a, least(t1.a, t2.a, t3.a) phv
					  from join_pt1 t2 join join_ut1 t3 on t2.a = t3.b) ss
              on t1.a = ss.t2a order by t1.a;
select t1.b, ss.phv from join_ut1 t1 left join lateral
              (select t2.a as t2a, t3.a t3a, least(t1.a, t2.a, t3.a) phv
					  from join_pt1 t2 join join_ut1 t3 on t2.a = t3.b) ss
              on t1.a = ss.t2a order by t1.a;

drop table join_pt1;
drop table join_ut1;

--
-- test estimation behavior with multi-column foreign key and constant qual
--

begin;

create table fkest (x integer, x10 integer, x10b integer, x100 integer);
insert into fkest select x, x/10, x/10, x/100 from generate_series(1,1000) x;
create unique index on fkest(x, x10, x100);
analyze fkest;

explain (costs off)
select * from fkest f1
  join fkest f2 on (f1.x = f2.x and f1.x10 = f2.x10b and f1.x100 = f2.x100)
  join fkest f3 on f1.x = f3.x
  where f1.x100 = 2;

alter table fkest add constraint fk
  foreign key (x, x10b, x100) references fkest (x, x10, x100);

explain (costs off)
select * from fkest f1
  join fkest f2 on (f1.x = f2.x and f1.x10 = f2.x10b and f1.x100 = f2.x100)
  join fkest f3 on f1.x = f3.x
  where f1.x100 = 2;

rollback;

--
-- test that foreign key join estimation performs sanely for outer joins
--

begin;

create table fkest (a int, b int, c int unique, primary key(a,b));
create table fkest1 (a int, b int, primary key(a,b));

insert into fkest select x/10, x%10, x from generate_series(1,1000) x;
insert into fkest1 select x/10, x%10 from generate_series(1,1000) x;

alter table fkest1
  add constraint fkest1_a_b_fkey foreign key (a,b) references fkest;

analyze fkest;
analyze fkest1;

explain (costs off)
select *
from fkest f
  left join fkest1 f1 on f.a = f1.a and f.b = f1.b
  left join fkest1 f2 on f.a = f2.a and f.b = f2.b
  left join fkest1 f3 on f.a = f3.a and f.b = f3.b
where f.c = 1;

rollback;

--
-- test planner's ability to mark joins as unique
--

create table j1 (id int primary key);
create table j2 (id int primary key);
create table j3 (id int);

insert into j1 values(1),(2),(3);
insert into j2 values(1),(2),(3);
insert into j3 values(1),(1);

analyze j1;
analyze j2;
analyze j3;

-- ensure join is properly marked as unique
explain (verbose, costs off)
select * from j1 inner join j2 on j1.id = j2.id;

-- ensure join is not unique when not an equi-join
explain (verbose, costs off)
select * from j1 inner join j2 on j1.id > j2.id;

-- ensure non-unique rel is not chosen as inner
explain (verbose, costs off)
select * from j1 inner join j3 on j1.id = j3.id;

-- ensure left join is marked as unique
explain (verbose, costs off)
select * from j1 left join j2 on j1.id = j2.id;

-- ensure right join is marked as unique
explain (verbose, costs off)
select * from j1 right join j2 on j1.id = j2.id;

-- ensure full join is marked as unique
explain (verbose, costs off)
select * from j1 full join j2 on j1.id = j2.id;

-- a clauseless (cross) join can't be unique
explain (verbose, costs off)
select * from j1 cross join j2;

-- ensure a natural join is marked as unique
explain (verbose, costs off)
select * from j1 natural join j2;

-- ensure a distinct clause allows the inner to become unique
explain (verbose, costs off)
select * from j1
inner join (select distinct id from j3) j3 on j1.id = j3.id;

-- ensure group by clause allows the inner to become unique
explain (verbose, costs off)
select * from j1
inner join (select id from j3 group by id) j3 on j1.id = j3.id;

drop table j1;
drop table j2;
drop table j3;

-- test more complex permutations of unique joins

create table j1 (id1 int, id2 int, primary key(id1,id2));
create table j2 (id1 int, id2 int, primary key(id1,id2));
create table j3 (id1 int, id2 int, primary key(id1,id2));

insert into j1 values(1,1),(1,2);
insert into j2 values(1,1);
insert into j3 values(1,1);

analyze j1;
analyze j2;
analyze j3;

-- ensure there's no unique join when not all columns which are part of the
-- unique index are seen in the join clause
explain (verbose, costs off)
select * from j1
inner join j2 on j1.id1 = j2.id1;

-- ensure proper unique detection with multiple join quals
explain (verbose, costs off)
select * from j1
inner join j2 on j1.id1 = j2.id1 and j1.id2 = j2.id2;

-- ensure we don't detect the join to be unique when quals are not part of the
-- join condition
explain (verbose, costs off)
select * from j1
inner join j2 on j1.id1 = j2.id1 where j1.id2 = 1;

-- as above, but for left joins.
explain (verbose, costs off)
select * from j1
left join j2 on j1.id1 = j2.id1 where j1.id2 = 1;

create unique index j1_id2_idx on j1(id2) where id2 is not null;

-- ensure we don't use a partial unique index as unique proofs
explain (verbose, costs off)
select * from j1
inner join j2 on j1.id2 = j2.id2;

drop index j1_id2_idx;

-- validate logic in merge joins which skips mark and restore.
-- it should only do this if all quals which were used to detect the unique
-- are present as join quals, and not plain quals.
set enable_nestloop to 0;
set enable_hashjoin to 0;
set enable_sort to 0;

-- create indexes that will be preferred over the PKs to perform the join
create index j1_id1_idx on j1 (id1) where id1 % 1000 = 1;
create index j2_id1_idx on j2 (id1) where id1 % 1000 = 1;

-- need an additional row in j2, if we want j2_id1_idx to be preferred
insert into j2 values(1,2);
analyze j2;

explain (costs off) select * from j1
inner join j2 on j1.id1 = j2.id1 and j1.id2 = j2.id2
where j1.id1 % 1000 = 1 and j2.id1 % 1000 = 1;

select * from j1
inner join j2 on j1.id1 = j2.id1 and j1.id2 = j2.id2
where j1.id1 % 1000 = 1 and j2.id1 % 1000 = 1;

-- Exercise array keys mark/restore B-Tree code
explain (costs off) select * from j1
inner join j2 on j1.id1 = j2.id1 and j1.id2 = j2.id2
where j1.id1 % 1000 = 1 and j2.id1 % 1000 = 1 and j2.id1 = any (array[1]);

select * from j1
inner join j2 on j1.id1 = j2.id1 and j1.id2 = j2.id2
where j1.id1 % 1000 = 1 and j2.id1 % 1000 = 1 and j2.id1 = any (array[1]);

-- Exercise array keys "find extreme element" B-Tree code
explain (costs off) select * from j1
inner join j2 on j1.id1 = j2.id1 and j1.id2 = j2.id2
where j1.id1 % 1000 = 1 and j2.id1 % 1000 = 1 and j2.id1 >= any (array[1,5]);

select * from j1
inner join j2 on j1.id1 = j2.id1 and j1.id2 = j2.id2
where j1.id1 % 1000 = 1 and j2.id1 % 1000 = 1 and j2.id1 >= any (array[1,5]);

reset enable_nestloop;
reset enable_hashjoin;
reset enable_sort;

drop table j1;
drop table j2;
drop table j3;

-- check that semijoin inner is not seen as unique for a portion of the outerrel
explain (verbose, costs off)
select t1.unique1, t2.hundred
from onek t1, tenk1 t2
where exists (select 1 from tenk1 t3
              where t3.thousand = t1.unique1 and t3.tenthous = t2.hundred)
      and t1.unique1 < 1;

-- ... unless it actually is unique
create table j3 as select unique1, tenthous from onek;
vacuum analyze j3;
create unique index on j3(unique1, tenthous);

explain (verbose, costs off)
select t1.unique1, t2.hundred
from onek t1, tenk1 t2
where exists (select 1 from j3
              where j3.unique1 = t1.unique1 and j3.tenthous = t2.hundred)
      and t1.unique1 < 1;

drop table j3;

-- Exercise the "skip fetch" Bitmap Heap Scan optimization when candidate
-- tuples are discarded. This may occur when:
--   1. A join doesn't require all inner tuples to be scanned for each outer
--      tuple, and
--   2. The inner side is scanned using a bitmap heap scan, and
--   3. The bitmap heap scan is eligible for the "skip fetch" optimization.
--      This optimization is usable when no data from the underlying table is
--      needed. Use a temp table so it is only visible to this backend and
--      vacuum may reliably mark all blocks in the table all visible in the
--      visibility map.
CREATE TEMP TABLE skip_fetch (a INT, b INT) WITH (fillfactor=10);
INSERT INTO skip_fetch SELECT i % 3, i FROM generate_series(0,30) i;
CREATE INDEX ON skip_fetch(a);
VACUUM (ANALYZE) skip_fetch;

SET enable_indexonlyscan = off;
SET enable_seqscan = off;
EXPLAIN (COSTS OFF)
SELECT t1.a FROM skip_fetch t1 LEFT JOIN skip_fetch t2 ON t2.a = 1 WHERE t2.a IS NULL;
SELECT t1.a FROM skip_fetch t1 LEFT JOIN skip_fetch t2 ON t2.a = 1 WHERE t2.a IS NULL;

RESET enable_indexonlyscan;
RESET enable_seqscan;

-- Test that we do not account for nullingrels when looking up statistics
CREATE TABLE group_tbl (a INT, b INT);
INSERT INTO group_tbl SELECT 1, 1;
CREATE STATISTICS group_tbl_stat (ndistinct) ON a, b FROM group_tbl;
ANALYZE group_tbl;

EXPLAIN (COSTS OFF)
SELECT 1 FROM group_tbl t1
    LEFT JOIN (SELECT a c1, COALESCE(a) c2 FROM group_tbl t2) s ON TRUE
GROUP BY s.c1, s.c2;

DROP TABLE group_tbl;

--
-- PARTITION_JOIN
-- Test partitionwise join between partitioned tables
--

-- Enable partitionwise join, which by default is disabled.
SET enable_partitionwise_join to true;

--
-- partitioned by a single column
--
CREATE TABLE prt1 (a int, b int, c varchar) PARTITION BY RANGE(a);
CREATE TABLE prt1_p1 PARTITION OF prt1 FOR VALUES FROM (0) TO (250);
CREATE TABLE prt1_p3 PARTITION OF prt1 FOR VALUES FROM (500) TO (600);
CREATE TABLE prt1_p2 PARTITION OF prt1 FOR VALUES FROM (250) TO (500);
INSERT INTO prt1 SELECT i, i % 25, to_char(i, 'FM0000') FROM generate_series(0, 599) i WHERE i % 2 = 0;
CREATE INDEX iprt1_p1_a on prt1_p1(a);
CREATE INDEX iprt1_p2_a on prt1_p2(a);
CREATE INDEX iprt1_p3_a on prt1_p3(a);
ANALYZE prt1;

CREATE TABLE prt2 (a int, b int, c varchar) PARTITION BY RANGE(b);
CREATE TABLE prt2_p1 PARTITION OF prt2 FOR VALUES FROM (0) TO (250);
CREATE TABLE prt2_p2 PARTITION OF prt2 FOR VALUES FROM (250) TO (500);
CREATE TABLE prt2_p3 PARTITION OF prt2 FOR VALUES FROM (500) TO (600);
INSERT INTO prt2 SELECT i % 25, i, to_char(i, 'FM0000') FROM generate_series(0, 599) i WHERE i % 3 = 0;
CREATE INDEX iprt2_p1_b on prt2_p1(b);
CREATE INDEX iprt2_p2_b on prt2_p2(b);
CREATE INDEX iprt2_p3_b on prt2_p3(b);
ANALYZE prt2;

-- inner join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1 t1, prt2 t2 WHERE t1.a = t2.b AND t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1 t1, prt2 t2 WHERE t1.a = t2.b AND t1.b = 0 ORDER BY t1.a, t2.b;

-- left outer join, 3-way
EXPLAIN (COSTS OFF)
SELECT COUNT(*) FROM prt1 t1
  LEFT JOIN prt1 t2 ON t1.a = t2.a
  LEFT JOIN prt1 t3 ON t2.a = t3.a;
SELECT COUNT(*) FROM prt1 t1
  LEFT JOIN prt1 t2 ON t1.a = t2.a
  LEFT JOIN prt1 t3 ON t2.a = t3.a;

-- left outer join, with whole-row reference; partitionwise join does not apply
EXPLAIN (COSTS OFF)
SELECT t1, t2 FROM prt1 t1 LEFT JOIN prt2 t2 ON t1.a = t2.b WHERE t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1, t2 FROM prt1 t1 LEFT JOIN prt2 t2 ON t1.a = t2.b WHERE t1.b = 0 ORDER BY t1.a, t2.b;

-- right outer join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1 t1 RIGHT JOIN prt2 t2 ON t1.a = t2.b WHERE t2.a = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1 t1 RIGHT JOIN prt2 t2 ON t1.a = t2.b WHERE t2.a = 0 ORDER BY t1.a, t2.b;

-- full outer join, with placeholder vars
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT 50 phv, * FROM prt1 WHERE prt1.b = 0) t1 FULL JOIN (SELECT 75 phv, * FROM prt2 WHERE prt2.a = 0) t2 ON (t1.a = t2.b) WHERE t1.phv = t1.a OR t2.phv = t2.b ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT 50 phv, * FROM prt1 WHERE prt1.b = 0) t1 FULL JOIN (SELECT 75 phv, * FROM prt2 WHERE prt2.a = 0) t2 ON (t1.a = t2.b) WHERE t1.phv = t1.a OR t2.phv = t2.b ORDER BY t1.a, t2.b;

-- Join with pruned partitions from joining relations
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1 t1, prt2 t2 WHERE t1.a = t2.b AND t1.a < 450 AND t2.b > 250 AND t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1 t1, prt2 t2 WHERE t1.a = t2.b AND t1.a < 450 AND t2.b > 250 AND t1.b = 0 ORDER BY t1.a, t2.b;

-- Currently we can't do partitioned join if nullable-side partitions are pruned
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1 WHERE a < 450) t1 LEFT JOIN (SELECT * FROM prt2 WHERE b > 250) t2 ON t1.a = t2.b WHERE t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1 WHERE a < 450) t1 LEFT JOIN (SELECT * FROM prt2 WHERE b > 250) t2 ON t1.a = t2.b WHERE t1.b = 0 ORDER BY t1.a, t2.b;

-- Currently we can't do partitioned join if nullable-side partitions are pruned
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1 WHERE a < 450) t1 FULL JOIN (SELECT * FROM prt2 WHERE b > 250) t2 ON t1.a = t2.b WHERE t1.b = 0 OR t2.a = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1 WHERE a < 450) t1 FULL JOIN (SELECT * FROM prt2 WHERE b > 250) t2 ON t1.a = t2.b WHERE t1.b = 0 OR t2.a = 0 ORDER BY t1.a, t2.b;

-- Semi-join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM prt1 t1 WHERE t1.a IN (SELECT t2.b FROM prt2 t2 WHERE t2.a = 0) AND t1.b = 0 ORDER BY t1.a;
SELECT t1.* FROM prt1 t1 WHERE t1.a IN (SELECT t2.b FROM prt2 t2 WHERE t2.a = 0) AND t1.b = 0 ORDER BY t1.a;

-- Anti-join with aggregates
EXPLAIN (COSTS OFF)
SELECT sum(t1.a), avg(t1.a), sum(t1.b), avg(t1.b) FROM prt1 t1 WHERE NOT EXISTS (SELECT 1 FROM prt2 t2 WHERE t1.a = t2.b);
SELECT sum(t1.a), avg(t1.a), sum(t1.b), avg(t1.b) FROM prt1 t1 WHERE NOT EXISTS (SELECT 1 FROM prt2 t2 WHERE t1.a = t2.b);

-- lateral reference
EXPLAIN (COSTS OFF)
SELECT * FROM prt1 t1 LEFT JOIN LATERAL
			  (SELECT t2.a AS t2a, t3.a AS t3a, least(t1.a,t2.a,t3.b) FROM prt1 t2 JOIN prt2 t3 ON (t2.a = t3.b)) ss
			  ON t1.a = ss.t2a WHERE t1.b = 0 ORDER BY t1.a;
SELECT * FROM prt1 t1 LEFT JOIN LATERAL
			  (SELECT t2.a AS t2a, t3.a AS t3a, least(t1.a,t2.a,t3.b) FROM prt1 t2 JOIN prt2 t3 ON (t2.a = t3.b)) ss
			  ON t1.a = ss.t2a WHERE t1.b = 0 ORDER BY t1.a;

EXPLAIN (COSTS OFF)
SELECT t1.a, ss.t2a, ss.t2c FROM prt1 t1 LEFT JOIN LATERAL
			  (SELECT t2.a AS t2a, t3.a AS t3a, t2.b t2b, t2.c t2c, least(t1.a,t2.a,t3.b) FROM prt1 t2 JOIN prt2 t3 ON (t2.a = t3.b)) ss
			  ON t1.c = ss.t2c WHERE (t1.b + coalesce(ss.t2b, 0)) = 0 ORDER BY t1.a;
SELECT t1.a, ss.t2a, ss.t2c FROM prt1 t1 LEFT JOIN LATERAL
			  (SELECT t2.a AS t2a, t3.a AS t3a, t2.b t2b, t2.c t2c, least(t1.a,t2.a,t3.a) FROM prt1 t2 JOIN prt2 t3 ON (t2.a = t3.b)) ss
			  ON t1.c = ss.t2c WHERE (t1.b + coalesce(ss.t2b, 0)) = 0 ORDER BY t1.a;

-- lateral reference in sample scan
EXPLAIN (COSTS OFF)
SELECT * FROM prt1 t1 JOIN LATERAL
			  (SELECT * FROM prt1 t2 TABLESAMPLE SYSTEM (t1.a) REPEATABLE(t1.b)) s
			  ON t1.a = s.a;

-- lateral reference in scan's restriction clauses
EXPLAIN (COSTS OFF)
SELECT count(*) FROM prt1 t1 LEFT JOIN LATERAL
			  (SELECT t1.b AS t1b, t2.* FROM prt2 t2) s
			  ON t1.a = s.b WHERE s.t1b = s.a;
SELECT count(*) FROM prt1 t1 LEFT JOIN LATERAL
			  (SELECT t1.b AS t1b, t2.* FROM prt2 t2) s
			  ON t1.a = s.b WHERE s.t1b = s.a;

EXPLAIN (COSTS OFF)
SELECT count(*) FROM prt1 t1 LEFT JOIN LATERAL
			  (SELECT t1.b AS t1b, t2.* FROM prt2 t2) s
			  ON t1.a = s.b WHERE s.t1b = s.b;
SELECT count(*) FROM prt1 t1 LEFT JOIN LATERAL
			  (SELECT t1.b AS t1b, t2.* FROM prt2 t2) s
			  ON t1.a = s.b WHERE s.t1b = s.b;

-- bug with inadequate sort key representation
SET enable_partitionwise_aggregate TO true;
SET enable_hashjoin TO false;

EXPLAIN (COSTS OFF)
SELECT a, b FROM prt1 FULL JOIN prt2 p2(b,a,c) USING(a,b)
  WHERE a BETWEEN 490 AND 510
  GROUP BY 1, 2 ORDER BY 1, 2;
SELECT a, b FROM prt1 FULL JOIN prt2 p2(b,a,c) USING(a,b)
  WHERE a BETWEEN 490 AND 510
  GROUP BY 1, 2 ORDER BY 1, 2;

RESET enable_partitionwise_aggregate;
RESET enable_hashjoin;

--
-- partitioned by expression
--
CREATE TABLE prt1_e (a int, b int, c int) PARTITION BY RANGE(((a + b)/2));
CREATE TABLE prt1_e_p1 PARTITION OF prt1_e FOR VALUES FROM (0) TO (250);
CREATE TABLE prt1_e_p2 PARTITION OF prt1_e FOR VALUES FROM (250) TO (500);
CREATE TABLE prt1_e_p3 PARTITION OF prt1_e FOR VALUES FROM (500) TO (600);
INSERT INTO prt1_e SELECT i, i, i % 25 FROM generate_series(0, 599, 2) i;
CREATE INDEX iprt1_e_p1_ab2 on prt1_e_p1(((a+b)/2));
CREATE INDEX iprt1_e_p2_ab2 on prt1_e_p2(((a+b)/2));
CREATE INDEX iprt1_e_p3_ab2 on prt1_e_p3(((a+b)/2));
ANALYZE prt1_e;

CREATE TABLE prt2_e (a int, b int, c int) PARTITION BY RANGE(((b + a)/2));
CREATE TABLE prt2_e_p1 PARTITION OF prt2_e FOR VALUES FROM (0) TO (250);
CREATE TABLE prt2_e_p2 PARTITION OF prt2_e FOR VALUES FROM (250) TO (500);
CREATE TABLE prt2_e_p3 PARTITION OF prt2_e FOR VALUES FROM (500) TO (600);
INSERT INTO prt2_e SELECT i, i, i % 25 FROM generate_series(0, 599, 3) i;
ANALYZE prt2_e;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_e t1, prt2_e t2 WHERE (t1.a + t1.b)/2 = (t2.b + t2.a)/2 AND t1.c = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_e t1, prt2_e t2 WHERE (t1.a + t1.b)/2 = (t2.b + t2.a)/2 AND t1.c = 0 ORDER BY t1.a, t2.b;

--
-- N-way join
--
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c, t3.a + t3.b, t3.c FROM prt1 t1, prt2 t2, prt1_e t3 WHERE t1.a = t2.b AND t1.a = (t3.a + t3.b)/2 AND t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c, t3.a + t3.b, t3.c FROM prt1 t1, prt2 t2, prt1_e t3 WHERE t1.a = t2.b AND t1.a = (t3.a + t3.b)/2 AND t1.b = 0 ORDER BY t1.a, t2.b;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c, t3.a + t3.b, t3.c FROM (prt1 t1 LEFT JOIN prt2 t2 ON t1.a = t2.b) LEFT JOIN prt1_e t3 ON (t1.a = (t3.a + t3.b)/2) WHERE t1.b = 0 ORDER BY t1.a, t2.b, t3.a + t3.b;
SELECT t1.a, t1.c, t2.b, t2.c, t3.a + t3.b, t3.c FROM (prt1 t1 LEFT JOIN prt2 t2 ON t1.a = t2.b) LEFT JOIN prt1_e t3 ON (t1.a = (t3.a + t3.b)/2) WHERE t1.b = 0 ORDER BY t1.a, t2.b, t3.a + t3.b;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c, t3.a + t3.b, t3.c FROM (prt1 t1 LEFT JOIN prt2 t2 ON t1.a = t2.b) RIGHT JOIN prt1_e t3 ON (t1.a = (t3.a + t3.b)/2) WHERE t3.c = 0 ORDER BY t1.a, t2.b, t3.a + t3.b;
SELECT t1.a, t1.c, t2.b, t2.c, t3.a + t3.b, t3.c FROM (prt1 t1 LEFT JOIN prt2 t2 ON t1.a = t2.b) RIGHT JOIN prt1_e t3 ON (t1.a = (t3.a + t3.b)/2) WHERE t3.c = 0 ORDER BY t1.a, t2.b, t3.a + t3.b;

--
-- 3-way full join
--
EXPLAIN (COSTS OFF)
SELECT COUNT(*) FROM prt1 FULL JOIN prt2 p2(b,a,c) USING(a,b) FULL JOIN prt2 p3(b,a,c) USING (a, b)
  WHERE a BETWEEN 490 AND 510;
SELECT COUNT(*) FROM prt1 FULL JOIN prt2 p2(b,a,c) USING(a,b) FULL JOIN prt2 p3(b,a,c) USING (a, b)
  WHERE a BETWEEN 490 AND 510;

--
-- 4-way full join
--
EXPLAIN (COSTS OFF)
SELECT COUNT(*) FROM prt1 FULL JOIN prt2 p2(b,a,c) USING(a,b) FULL JOIN prt2 p3(b,a,c) USING (a, b) FULL JOIN prt1 p4 (a,b,c) USING (a, b)
  WHERE a BETWEEN 490 AND 510;
SELECT COUNT(*) FROM prt1 FULL JOIN prt2 p2(b,a,c) USING(a,b) FULL JOIN prt2 p3(b,a,c) USING (a, b) FULL JOIN prt1 p4 (a,b,c) USING (a, b)
  WHERE a BETWEEN 490 AND 510;

-- Cases with non-nullable expressions in subquery results;
-- make sure these go to null as expected
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.phv, t2.b, t2.phv, t3.a + t3.b, t3.phv FROM ((SELECT 50 phv, * FROM prt1 WHERE prt1.b = 0) t1 FULL JOIN (SELECT 75 phv, * FROM prt2 WHERE prt2.a = 0) t2 ON (t1.a = t2.b)) FULL JOIN (SELECT 50 phv, * FROM prt1_e WHERE prt1_e.c = 0) t3 ON (t1.a = (t3.a + t3.b)/2) WHERE t1.a = t1.phv OR t2.b = t2.phv OR (t3.a + t3.b)/2 = t3.phv ORDER BY t1.a, t2.b, t3.a + t3.b;
SELECT t1.a, t1.phv, t2.b, t2.phv, t3.a + t3.b, t3.phv FROM ((SELECT 50 phv, * FROM prt1 WHERE prt1.b = 0) t1 FULL JOIN (SELECT 75 phv, * FROM prt2 WHERE prt2.a = 0) t2 ON (t1.a = t2.b)) FULL JOIN (SELECT 50 phv, * FROM prt1_e WHERE prt1_e.c = 0) t3 ON (t1.a = (t3.a + t3.b)/2) WHERE t1.a = t1.phv OR t2.b = t2.phv OR (t3.a + t3.b)/2 = t3.phv ORDER BY t1.a, t2.b, t3.a + t3.b;

-- Semi-join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM prt1 t1 WHERE t1.a IN (SELECT t1.b FROM prt2 t1, prt1_e t2 WHERE t1.a = 0 AND t1.b = (t2.a + t2.b)/2) AND t1.b = 0 ORDER BY t1.a;
SELECT t1.* FROM prt1 t1 WHERE t1.a IN (SELECT t1.b FROM prt2 t1, prt1_e t2 WHERE t1.a = 0 AND t1.b = (t2.a + t2.b)/2) AND t1.b = 0 ORDER BY t1.a;

EXPLAIN (COSTS OFF)
SELECT t1.* FROM prt1 t1 WHERE t1.a IN (SELECT t1.b FROM prt2 t1 WHERE t1.b IN (SELECT (t1.a + t1.b)/2 FROM prt1_e t1 WHERE t1.c = 0)) AND t1.b = 0 ORDER BY t1.a;
SELECT t1.* FROM prt1 t1 WHERE t1.a IN (SELECT t1.b FROM prt2 t1 WHERE t1.b IN (SELECT (t1.a + t1.b)/2 FROM prt1_e t1 WHERE t1.c = 0)) AND t1.b = 0 ORDER BY t1.a;

-- test merge joins
SET enable_hashjoin TO off;
SET enable_nestloop TO off;

EXPLAIN (COSTS OFF)
SELECT t1.* FROM prt1 t1 WHERE t1.a IN (SELECT t1.b FROM prt2 t1 WHERE t1.b IN (SELECT (t1.a + t1.b)/2 FROM prt1_e t1 WHERE t1.c = 0)) AND t1.b = 0 ORDER BY t1.a;
SELECT t1.* FROM prt1 t1 WHERE t1.a IN (SELECT t1.b FROM prt2 t1 WHERE t1.b IN (SELECT (t1.a + t1.b)/2 FROM prt1_e t1 WHERE t1.c = 0)) AND t1.b = 0 ORDER BY t1.a;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c, t3.a + t3.b, t3.c FROM (prt1 t1 LEFT JOIN prt2 t2 ON t1.a = t2.b) RIGHT JOIN prt1_e t3 ON (t1.a = (t3.a + t3.b)/2) WHERE t3.c = 0 ORDER BY t1.a, t2.b, t3.a + t3.b;
SELECT t1.a, t1.c, t2.b, t2.c, t3.a + t3.b, t3.c FROM (prt1 t1 LEFT JOIN prt2 t2 ON t1.a = t2.b) RIGHT JOIN prt1_e t3 ON (t1.a = (t3.a + t3.b)/2) WHERE t3.c = 0 ORDER BY t1.a, t2.b, t3.a + t3.b;

-- MergeAppend on nullable column
-- This should generate a partitionwise join, but currently fails to
EXPLAIN (COSTS OFF)
SELECT t1.a, t2.b FROM (SELECT * FROM prt1 WHERE a < 450) t1 LEFT JOIN (SELECT * FROM prt2 WHERE b > 250) t2 ON t1.a = t2.b WHERE t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t2.b FROM (SELECT * FROM prt1 WHERE a < 450) t1 LEFT JOIN (SELECT * FROM prt2 WHERE b > 250) t2 ON t1.a = t2.b WHERE t1.b = 0 ORDER BY t1.a, t2.b;

-- merge join when expression with whole-row reference needs to be sorted;
-- partitionwise join does not apply
EXPLAIN (COSTS OFF)
SELECT t1.a, t2.b FROM prt1 t1, prt2 t2 WHERE t1::text = t2::text AND t1.a = t2.b ORDER BY t1.a;
SELECT t1.a, t2.b FROM prt1 t1, prt2 t2 WHERE t1::text = t2::text AND t1.a = t2.b ORDER BY t1.a;

RESET enable_hashjoin;
RESET enable_nestloop;

--
-- partitioned by multiple columns
--
CREATE TABLE prt1_m (a int, b int, c int) PARTITION BY RANGE(a, ((a + b)/2));
CREATE TABLE prt1_m_p1 PARTITION OF prt1_m FOR VALUES FROM (0, 0) TO (250, 250);
CREATE TABLE prt1_m_p2 PARTITION OF prt1_m FOR VALUES FROM (250, 250) TO (500, 500);
CREATE TABLE prt1_m_p3 PARTITION OF prt1_m FOR VALUES FROM (500, 500) TO (600, 600);
INSERT INTO prt1_m SELECT i, i, i % 25 FROM generate_series(0, 599, 2) i;
ANALYZE prt1_m;

CREATE TABLE prt2_m (a int, b int, c int) PARTITION BY RANGE(((b + a)/2), b);
CREATE TABLE prt2_m_p1 PARTITION OF prt2_m FOR VALUES FROM (0, 0) TO (250, 250);
CREATE TABLE prt2_m_p2 PARTITION OF prt2_m FOR VALUES FROM (250, 250) TO (500, 500);
CREATE TABLE prt2_m_p3 PARTITION OF prt2_m FOR VALUES FROM (500, 500) TO (600, 600);
INSERT INTO prt2_m SELECT i, i, i % 25 FROM generate_series(0, 599, 3) i;
ANALYZE prt2_m;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1_m WHERE prt1_m.c = 0) t1 FULL JOIN (SELECT * FROM prt2_m WHERE prt2_m.c = 0) t2 ON (t1.a = (t2.b + t2.a)/2 AND t2.b = (t1.a + t1.b)/2) ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1_m WHERE prt1_m.c = 0) t1 FULL JOIN (SELECT * FROM prt2_m WHERE prt2_m.c = 0) t2 ON (t1.a = (t2.b + t2.a)/2 AND t2.b = (t1.a + t1.b)/2) ORDER BY t1.a, t2.b;

--
-- tests for list partitioned tables.
--
CREATE TABLE plt1 (a int, b int, c text) PARTITION BY LIST(c);
CREATE TABLE plt1_p1 PARTITION OF plt1 FOR VALUES IN ('0000', '0003', '0004', '0010');
CREATE TABLE plt1_p2 PARTITION OF plt1 FOR VALUES IN ('0001', '0005', '0002', '0009');
CREATE TABLE plt1_p3 PARTITION OF plt1 FOR VALUES IN ('0006', '0007', '0008', '0011');
INSERT INTO plt1 SELECT i, i, to_char(i/50, 'FM0000') FROM generate_series(0, 599, 2) i;
ANALYZE plt1;

CREATE TABLE plt2 (a int, b int, c text) PARTITION BY LIST(c);
CREATE TABLE plt2_p1 PARTITION OF plt2 FOR VALUES IN ('0000', '0003', '0004', '0010');
CREATE TABLE plt2_p2 PARTITION OF plt2 FOR VALUES IN ('0001', '0005', '0002', '0009');
CREATE TABLE plt2_p3 PARTITION OF plt2 FOR VALUES IN ('0006', '0007', '0008', '0011');
INSERT INTO plt2 SELECT i, i, to_char(i/50, 'FM0000') FROM generate_series(0, 599, 3) i;
ANALYZE plt2;

--
-- list partitioned by expression
--
CREATE TABLE plt1_e (a int, b int, c text) PARTITION BY LIST(ltrim(c, 'A'));
CREATE TABLE plt1_e_p1 PARTITION OF plt1_e FOR VALUES IN ('0000', '0003', '0004', '0010');
CREATE TABLE plt1_e_p2 PARTITION OF plt1_e FOR VALUES IN ('0001', '0005', '0002', '0009');
CREATE TABLE plt1_e_p3 PARTITION OF plt1_e FOR VALUES IN ('0006', '0007', '0008', '0011');
INSERT INTO plt1_e SELECT i, i, 'A' || to_char(i/50, 'FM0000') FROM generate_series(0, 599, 2) i;
ANALYZE plt1_e;

-- test partition matching with N-way join
EXPLAIN (COSTS OFF)
SELECT avg(t1.a), avg(t2.b), avg(t3.a + t3.b), t1.c, t2.c, t3.c FROM plt1 t1, plt2 t2, plt1_e t3 WHERE t1.b = t2.b AND t1.c = t2.c AND ltrim(t3.c, 'A') = t1.c GROUP BY t1.c, t2.c, t3.c ORDER BY t1.c, t2.c, t3.c;
SELECT avg(t1.a), avg(t2.b), avg(t3.a + t3.b), t1.c, t2.c, t3.c FROM plt1 t1, plt2 t2, plt1_e t3 WHERE t1.b = t2.b AND t1.c = t2.c AND ltrim(t3.c, 'A') = t1.c GROUP BY t1.c, t2.c, t3.c ORDER BY t1.c, t2.c, t3.c;

-- joins where one of the relations is proven empty
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1 t1, prt2 t2 WHERE t1.a = t2.b AND t1.a = 1 AND t1.a = 2;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1 WHERE a = 1 AND a = 2) t1 LEFT JOIN prt2 t2 ON t1.a = t2.b;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1 WHERE a = 1 AND a = 2) t1 RIGHT JOIN prt2 t2 ON t1.a = t2.b, prt1 t3 WHERE t2.b = t3.a;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1 WHERE a = 1 AND a = 2) t1 FULL JOIN prt2 t2 ON t1.a = t2.b WHERE t2.a = 0 ORDER BY t1.a, t2.b;

--
-- tests for hash partitioned tables.
--
CREATE TABLE pht1 (a int, b int, c text) PARTITION BY HASH(c);
CREATE TABLE pht1_p1 PARTITION OF pht1 FOR VALUES WITH (MODULUS 3, REMAINDER 0);
CREATE TABLE pht1_p2 PARTITION OF pht1 FOR VALUES WITH (MODULUS 3, REMAINDER 1);
CREATE TABLE pht1_p3 PARTITION OF pht1 FOR VALUES WITH (MODULUS 3, REMAINDER 2);
INSERT INTO pht1 SELECT i, i, to_char(i/50, 'FM0000') FROM generate_series(0, 599, 2) i;
ANALYZE pht1;

CREATE TABLE pht2 (a int, b int, c text) PARTITION BY HASH(c);
CREATE TABLE pht2_p1 PARTITION OF pht2 FOR VALUES WITH (MODULUS 3, REMAINDER 0);
CREATE TABLE pht2_p2 PARTITION OF pht2 FOR VALUES WITH (MODULUS 3, REMAINDER 1);
CREATE TABLE pht2_p3 PARTITION OF pht2 FOR VALUES WITH (MODULUS 3, REMAINDER 2);
INSERT INTO pht2 SELECT i, i, to_char(i/50, 'FM0000') FROM generate_series(0, 599, 3) i;
ANALYZE pht2;

--
-- hash partitioned by expression
--
CREATE TABLE pht1_e (a int, b int, c text) PARTITION BY HASH(ltrim(c, 'A'));
CREATE TABLE pht1_e_p1 PARTITION OF pht1_e FOR VALUES WITH (MODULUS 3, REMAINDER 0);
CREATE TABLE pht1_e_p2 PARTITION OF pht1_e FOR VALUES WITH (MODULUS 3, REMAINDER 1);
CREATE TABLE pht1_e_p3 PARTITION OF pht1_e FOR VALUES WITH (MODULUS 3, REMAINDER 2);
INSERT INTO pht1_e SELECT i, i, 'A' || to_char(i/50, 'FM0000') FROM generate_series(0, 299, 2) i;
ANALYZE pht1_e;

-- test partition matching with N-way join
EXPLAIN (COSTS OFF)
SELECT avg(t1.a), avg(t2.b), avg(t3.a + t3.b), t1.c, t2.c, t3.c FROM pht1 t1, pht2 t2, pht1_e t3 WHERE t1.b = t2.b AND t1.c = t2.c AND ltrim(t3.c, 'A') = t1.c GROUP BY t1.c, t2.c, t3.c ORDER BY t1.c, t2.c, t3.c;
SELECT avg(t1.a), avg(t2.b), avg(t3.a + t3.b), t1.c, t2.c, t3.c FROM pht1 t1, pht2 t2, pht1_e t3 WHERE t1.b = t2.b AND t1.c = t2.c AND ltrim(t3.c, 'A') = t1.c GROUP BY t1.c, t2.c, t3.c ORDER BY t1.c, t2.c, t3.c;

-- test default partition behavior for range
ALTER TABLE prt1 DETACH PARTITION prt1_p3;
ALTER TABLE prt1 ATTACH PARTITION prt1_p3 DEFAULT;
ANALYZE prt1;
ALTER TABLE prt2 DETACH PARTITION prt2_p3;
ALTER TABLE prt2 ATTACH PARTITION prt2_p3 DEFAULT;
ANALYZE prt2;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1 t1, prt2 t2 WHERE t1.a = t2.b AND t1.b = 0 ORDER BY t1.a, t2.b;

-- test default partition behavior for list
ALTER TABLE plt1 DETACH PARTITION plt1_p3;
ALTER TABLE plt1 ATTACH PARTITION plt1_p3 DEFAULT;
ANALYZE plt1;
ALTER TABLE plt2 DETACH PARTITION plt2_p3;
ALTER TABLE plt2 ATTACH PARTITION plt2_p3 DEFAULT;
ANALYZE plt2;

EXPLAIN (COSTS OFF)
SELECT avg(t1.a), avg(t2.b), t1.c, t2.c FROM plt1 t1 RIGHT JOIN plt2 t2 ON t1.c = t2.c WHERE t1.a % 25 = 0 GROUP BY t1.c, t2.c ORDER BY t1.c, t2.c;
--
-- multiple levels of partitioning
--
CREATE TABLE prt1_l (a int, b int, c varchar) PARTITION BY RANGE(a);
CREATE TABLE prt1_l_p1 PARTITION OF prt1_l FOR VALUES FROM (0) TO (250);
CREATE TABLE prt1_l_p2 PARTITION OF prt1_l FOR VALUES FROM (250) TO (500) PARTITION BY LIST (c);
CREATE TABLE prt1_l_p2_p1 PARTITION OF prt1_l_p2 FOR VALUES IN ('0000', '0001');
CREATE TABLE prt1_l_p2_p2 PARTITION OF prt1_l_p2 FOR VALUES IN ('0002', '0003');
CREATE TABLE prt1_l_p3 PARTITION OF prt1_l FOR VALUES FROM (500) TO (600) PARTITION BY RANGE (b);
CREATE TABLE prt1_l_p3_p1 PARTITION OF prt1_l_p3 FOR VALUES FROM (0) TO (13);
CREATE TABLE prt1_l_p3_p2 PARTITION OF prt1_l_p3 FOR VALUES FROM (13) TO (25);
INSERT INTO prt1_l SELECT i, i % 25, to_char(i % 4, 'FM0000') FROM generate_series(0, 599, 2) i;
ANALYZE prt1_l;

CREATE TABLE prt2_l (a int, b int, c varchar) PARTITION BY RANGE(b);
CREATE TABLE prt2_l_p1 PARTITION OF prt2_l FOR VALUES FROM (0) TO (250);
CREATE TABLE prt2_l_p2 PARTITION OF prt2_l FOR VALUES FROM (250) TO (500) PARTITION BY LIST (c);
CREATE TABLE prt2_l_p2_p1 PARTITION OF prt2_l_p2 FOR VALUES IN ('0000', '0001');
CREATE TABLE prt2_l_p2_p2 PARTITION OF prt2_l_p2 FOR VALUES IN ('0002', '0003');
CREATE TABLE prt2_l_p3 PARTITION OF prt2_l FOR VALUES FROM (500) TO (600) PARTITION BY RANGE (a);
CREATE TABLE prt2_l_p3_p1 PARTITION OF prt2_l_p3 FOR VALUES FROM (0) TO (13);
CREATE TABLE prt2_l_p3_p2 PARTITION OF prt2_l_p3 FOR VALUES FROM (13) TO (25);
INSERT INTO prt2_l SELECT i % 25, i, to_char(i % 4, 'FM0000') FROM generate_series(0, 599, 3) i;
ANALYZE prt2_l;

-- inner join, qual covering only top-level partitions
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_l t1, prt2_l t2 WHERE t1.a = t2.b AND t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_l t1, prt2_l t2 WHERE t1.a = t2.b AND t1.b = 0 ORDER BY t1.a, t2.b;

-- left join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_l t1 LEFT JOIN prt2_l t2 ON t1.a = t2.b AND t1.c = t2.c WHERE t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_l t1 LEFT JOIN prt2_l t2 ON t1.a = t2.b AND t1.c = t2.c WHERE t1.b = 0 ORDER BY t1.a, t2.b;

-- right join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_l t1 RIGHT JOIN prt2_l t2 ON t1.a = t2.b AND t1.c = t2.c WHERE t2.a = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_l t1 RIGHT JOIN prt2_l t2 ON t1.a = t2.b AND t1.c = t2.c WHERE t2.a = 0 ORDER BY t1.a, t2.b;

-- full join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1_l WHERE prt1_l.b = 0) t1 FULL JOIN (SELECT * FROM prt2_l WHERE prt2_l.a = 0) t2 ON (t1.a = t2.b AND t1.c = t2.c) ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1_l WHERE prt1_l.b = 0) t1 FULL JOIN (SELECT * FROM prt2_l WHERE prt2_l.a = 0) t2 ON (t1.a = t2.b AND t1.c = t2.c) ORDER BY t1.a, t2.b;

-- lateral partitionwise join
EXPLAIN (COSTS OFF)
SELECT * FROM prt1_l t1 LEFT JOIN LATERAL
			  (SELECT t2.a AS t2a, t2.c AS t2c, t2.b AS t2b, t3.b AS t3b, least(t1.a,t2.a,t3.b) FROM prt1_l t2 JOIN prt2_l t3 ON (t2.a = t3.b AND t2.c = t3.c)) ss
			  ON t1.a = ss.t2a AND t1.c = ss.t2c WHERE t1.b = 0 ORDER BY t1.a;
SELECT * FROM prt1_l t1 LEFT JOIN LATERAL
			  (SELECT t2.a AS t2a, t2.c AS t2c, t2.b AS t2b, t3.b AS t3b, least(t1.a,t2.a,t3.b) FROM prt1_l t2 JOIN prt2_l t3 ON (t2.a = t3.b AND t2.c = t3.c)) ss
			  ON t1.a = ss.t2a AND t1.c = ss.t2c WHERE t1.b = 0 ORDER BY t1.a;

-- partitionwise join with lateral reference in sample scan
EXPLAIN (COSTS OFF)
SELECT * FROM prt1_l t1 JOIN LATERAL
			  (SELECT * FROM prt1_l t2 TABLESAMPLE SYSTEM (t1.a) REPEATABLE(t1.b)) s
			  ON t1.a = s.a AND t1.b = s.b AND t1.c = s.c;

-- partitionwise join with lateral reference in scan's restriction clauses
EXPLAIN (COSTS OFF)
SELECT COUNT(*) FROM prt1_l t1 LEFT JOIN LATERAL
			  (SELECT t1.b AS t1b, t2.* FROM prt2_l t2) s
			  ON t1.a = s.b AND t1.b = s.a AND t1.c = s.c
			  WHERE s.t1b = s.a;
SELECT COUNT(*) FROM prt1_l t1 LEFT JOIN LATERAL
			  (SELECT t1.b AS t1b, t2.* FROM prt2_l t2) s
			  ON t1.a = s.b AND t1.b = s.a AND t1.c = s.c
			  WHERE s.t1b = s.a;

-- join with one side empty
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT * FROM prt1_l WHERE a = 1 AND a = 2) t1 RIGHT JOIN prt2_l t2 ON t1.a = t2.b AND t1.b = t2.a AND t1.c = t2.c;

-- Test case to verify proper handling of subqueries in a partitioned delete.
-- The weird-looking lateral join is just there to force creation of a
-- nestloop parameter within the subquery, which exposes the problem if the
-- planner fails to make multiple copies of the subquery as appropriate.
EXPLAIN (COSTS OFF)
DELETE FROM prt1_l
WHERE EXISTS (
  SELECT 1
    FROM int4_tbl,
         LATERAL (SELECT int4_tbl.f1 FROM int8_tbl LIMIT 2) ss
    WHERE prt1_l.c IS NULL);

--
-- negative testcases
--
CREATE TABLE prt1_n (a int, b int, c varchar) PARTITION BY RANGE(c);
CREATE TABLE prt1_n_p1 PARTITION OF prt1_n FOR VALUES FROM ('0000') TO ('0250');
CREATE TABLE prt1_n_p2 PARTITION OF prt1_n FOR VALUES FROM ('0250') TO ('0500');
INSERT INTO prt1_n SELECT i, i, to_char(i, 'FM0000') FROM generate_series(0, 499, 2) i;
ANALYZE prt1_n;

CREATE TABLE prt2_n (a int, b int, c text) PARTITION BY LIST(c);
CREATE TABLE prt2_n_p1 PARTITION OF prt2_n FOR VALUES IN ('0000', '0003', '0004', '0010', '0006', '0007');
CREATE TABLE prt2_n_p2 PARTITION OF prt2_n FOR VALUES IN ('0001', '0005', '0002', '0009', '0008', '0011');
INSERT INTO prt2_n SELECT i, i, to_char(i/50, 'FM0000') FROM generate_series(0, 599, 2) i;
ANALYZE prt2_n;

CREATE TABLE prt3_n (a int, b int, c text) PARTITION BY LIST(c);
CREATE TABLE prt3_n_p1 PARTITION OF prt3_n FOR VALUES IN ('0000', '0004', '0006', '0007');
CREATE TABLE prt3_n_p2 PARTITION OF prt3_n FOR VALUES IN ('0001', '0002', '0008', '0010');
CREATE TABLE prt3_n_p3 PARTITION OF prt3_n FOR VALUES IN ('0003', '0005', '0009', '0011');
INSERT INTO prt2_n SELECT i, i, to_char(i/50, 'FM0000') FROM generate_series(0, 599, 2) i;
ANALYZE prt3_n;

CREATE TABLE prt4_n (a int, b int, c text) PARTITION BY RANGE(a);
CREATE TABLE prt4_n_p1 PARTITION OF prt4_n FOR VALUES FROM (0) TO (300);
CREATE TABLE prt4_n_p2 PARTITION OF prt4_n FOR VALUES FROM (300) TO (500);
CREATE TABLE prt4_n_p3 PARTITION OF prt4_n FOR VALUES FROM (500) TO (600);
INSERT INTO prt4_n SELECT i, i, to_char(i, 'FM0000') FROM generate_series(0, 599, 2) i;
ANALYZE prt4_n;

-- partitionwise join can not be applied if the partition ranges differ
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1 t1, prt4_n t2 WHERE t1.a = t2.a;
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1 t1, prt4_n t2, prt2 t3 WHERE t1.a = t2.a and t1.a = t3.b;

-- partitionwise join can not be applied if there are no equi-join conditions
-- between partition keys
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1 t1 LEFT JOIN prt2 t2 ON (t1.a < t2.b);

-- equi-join with join condition on partial keys does not qualify for
-- partitionwise join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_m t1, prt2_m t2 WHERE t1.a = (t2.b + t2.a)/2;

-- equi-join between out-of-order partition key columns does not qualify for
-- partitionwise join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_m t1 LEFT JOIN prt2_m t2 ON t1.a = t2.b;

-- equi-join between non-key columns does not qualify for partitionwise join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_m t1 LEFT JOIN prt2_m t2 ON t1.c = t2.c;

-- partitionwise join can not be applied for a join between list and range
-- partitioned tables
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_n t1 LEFT JOIN prt2_n t2 ON (t1.c = t2.c);

-- partitionwise join can not be applied between tables with different
-- partition lists
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_n t1 JOIN prt2_n t2 ON (t1.c = t2.c) JOIN plt1 t3 ON (t1.c = t3.c);

-- partitionwise join can not be applied for a join between key column and
-- non-key column
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_n t1 FULL JOIN prt1 t2 ON (t1.c = t2.c);

--
-- Test some other plan types in a partitionwise join (unfortunately,
-- we need larger tables to get the planner to choose these plan types)
--
create temp table prtx1 (a integer, b integer, c integer)
  partition by range (a);
create temp table prtx1_1 partition of prtx1 for values from (1) to (11);
create temp table prtx1_2 partition of prtx1 for values from (11) to (21);
create temp table prtx1_3 partition of prtx1 for values from (21) to (31);
create temp table prtx2 (a integer, b integer, c integer)
  partition by range (a);
create temp table prtx2_1 partition of prtx2 for values from (1) to (11);
create temp table prtx2_2 partition of prtx2 for values from (11) to (21);
create temp table prtx2_3 partition of prtx2 for values from (21) to (31);
insert into prtx1 select 1 + i%30, i, i
  from generate_series(1,1000) i;
insert into prtx2 select 1 + i%30, i, i
  from generate_series(1,500) i, generate_series(1,10) j;
create index on prtx2 (b);
create index on prtx2 (c);
analyze prtx1;
analyze prtx2;

explain (costs off)
select * from prtx1
where not exists (select 1 from prtx2
                  where prtx2.a=prtx1.a and prtx2.b=prtx1.b and prtx2.c=123)
  and a<20 and c=120;

select * from prtx1
where not exists (select 1 from prtx2
                  where prtx2.a=prtx1.a and prtx2.b=prtx1.b and prtx2.c=123)
  and a<20 and c=120;

explain (costs off)
select * from prtx1
where not exists (select 1 from prtx2
                  where prtx2.a=prtx1.a and (prtx2.b=prtx1.b+1 or prtx2.c=99))
  and a<20 and c=91;

select * from prtx1
where not exists (select 1 from prtx2
                  where prtx2.a=prtx1.a and (prtx2.b=prtx1.b+1 or prtx2.c=99))
  and a<20 and c=91;

--
-- Test advanced partition-matching algorithm for partitioned join
--

-- Tests for range-partitioned tables
CREATE TABLE prt1_adv (a int, b int, c varchar) PARTITION BY RANGE (a);
CREATE TABLE prt1_adv_p1 PARTITION OF prt1_adv FOR VALUES FROM (100) TO (200);
CREATE TABLE prt1_adv_p2 PARTITION OF prt1_adv FOR VALUES FROM (200) TO (300);
CREATE TABLE prt1_adv_p3 PARTITION OF prt1_adv FOR VALUES FROM (300) TO (400);
CREATE INDEX prt1_adv_a_idx ON prt1_adv (a);
INSERT INTO prt1_adv SELECT i, i % 25, to_char(i, 'FM0000') FROM generate_series(100, 399) i;
ANALYZE prt1_adv;

CREATE TABLE prt2_adv (a int, b int, c varchar) PARTITION BY RANGE (b);
CREATE TABLE prt2_adv_p1 PARTITION OF prt2_adv FOR VALUES FROM (100) TO (150);
CREATE TABLE prt2_adv_p2 PARTITION OF prt2_adv FOR VALUES FROM (200) TO (300);
CREATE TABLE prt2_adv_p3 PARTITION OF prt2_adv FOR VALUES FROM (350) TO (500);
CREATE INDEX prt2_adv_b_idx ON prt2_adv (b);
INSERT INTO prt2_adv_p1 SELECT i % 25, i, to_char(i, 'FM0000') FROM generate_series(100, 149) i;
INSERT INTO prt2_adv_p2 SELECT i % 25, i, to_char(i, 'FM0000') FROM generate_series(200, 299) i;
INSERT INTO prt2_adv_p3 SELECT i % 25, i, to_char(i, 'FM0000') FROM generate_series(350, 499) i;
ANALYZE prt2_adv;

-- inner join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;

-- semi join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM prt1_adv t1 WHERE EXISTS (SELECT 1 FROM prt2_adv t2 WHERE t1.a = t2.b) AND t1.b = 0 ORDER BY t1.a;
SELECT t1.* FROM prt1_adv t1 WHERE EXISTS (SELECT 1 FROM prt2_adv t2 WHERE t1.a = t2.b) AND t1.b = 0 ORDER BY t1.a;

-- left join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 LEFT JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 LEFT JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;

-- anti join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM prt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM prt2_adv t2 WHERE t1.a = t2.b) AND t1.b = 0 ORDER BY t1.a;
SELECT t1.* FROM prt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM prt2_adv t2 WHERE t1.a = t2.b) AND t1.b = 0 ORDER BY t1.a;

-- full join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT 175 phv, * FROM prt1_adv WHERE prt1_adv.b = 0) t1 FULL JOIN (SELECT 425 phv, * FROM prt2_adv WHERE prt2_adv.a = 0) t2 ON (t1.a = t2.b) WHERE t1.phv = t1.a OR t2.phv = t2.b ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT 175 phv, * FROM prt1_adv WHERE prt1_adv.b = 0) t1 FULL JOIN (SELECT 425 phv, * FROM prt2_adv WHERE prt2_adv.a = 0) t2 ON (t1.a = t2.b) WHERE t1.phv = t1.a OR t2.phv = t2.b ORDER BY t1.a, t2.b;

-- Test cases where one side has an extra partition
CREATE TABLE prt2_adv_extra PARTITION OF prt2_adv FOR VALUES FROM (500) TO (MAXVALUE);
INSERT INTO prt2_adv SELECT i % 25, i, to_char(i, 'FM0000') FROM generate_series(500, 599) i;
ANALYZE prt2_adv;

-- inner join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;

-- semi join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM prt1_adv t1 WHERE EXISTS (SELECT 1 FROM prt2_adv t2 WHERE t1.a = t2.b) AND t1.b = 0 ORDER BY t1.a;
SELECT t1.* FROM prt1_adv t1 WHERE EXISTS (SELECT 1 FROM prt2_adv t2 WHERE t1.a = t2.b) AND t1.b = 0 ORDER BY t1.a;

-- left join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 LEFT JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 LEFT JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;

-- left join; currently we can't do partitioned join if there are no matched
-- partitions on the nullable side
EXPLAIN (COSTS OFF)
SELECT t1.b, t1.c, t2.a, t2.c FROM prt2_adv t1 LEFT JOIN prt1_adv t2 ON (t1.b = t2.a) WHERE t1.a = 0 ORDER BY t1.b, t2.a;

-- anti join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM prt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM prt2_adv t2 WHERE t1.a = t2.b) AND t1.b = 0 ORDER BY t1.a;
SELECT t1.* FROM prt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM prt2_adv t2 WHERE t1.a = t2.b) AND t1.b = 0 ORDER BY t1.a;

-- anti join; currently we can't do partitioned join if there are no matched
-- partitions on the nullable side
EXPLAIN (COSTS OFF)
SELECT t1.* FROM prt2_adv t1 WHERE NOT EXISTS (SELECT 1 FROM prt1_adv t2 WHERE t1.b = t2.a) AND t1.a = 0 ORDER BY t1.b;

-- full join; currently we can't do partitioned join if there are no matched
-- partitions on the nullable side
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT 175 phv, * FROM prt1_adv WHERE prt1_adv.b = 0) t1 FULL JOIN (SELECT 425 phv, * FROM prt2_adv WHERE prt2_adv.a = 0) t2 ON (t1.a = t2.b) WHERE t1.phv = t1.a OR t2.phv = t2.b ORDER BY t1.a, t2.b;

-- 3-way join where not every pair of relations can do partitioned join
EXPLAIN (COSTS OFF)
SELECT t1.b, t1.c, t2.a, t2.c, t3.a, t3.c FROM prt2_adv t1 LEFT JOIN prt1_adv t2 ON (t1.b = t2.a) INNER JOIN prt1_adv t3 ON (t1.b = t3.a) WHERE t1.a = 0 ORDER BY t1.b, t2.a, t3.a;
SELECT t1.b, t1.c, t2.a, t2.c, t3.a, t3.c FROM prt2_adv t1 LEFT JOIN prt1_adv t2 ON (t1.b = t2.a) INNER JOIN prt1_adv t3 ON (t1.b = t3.a) WHERE t1.a = 0 ORDER BY t1.b, t2.a, t3.a;

DROP TABLE prt2_adv_extra;

-- Test cases where a partition on one side matches multiple partitions on
-- the other side; we currently can't do partitioned join in such cases
ALTER TABLE prt2_adv DETACH PARTITION prt2_adv_p3;
-- Split prt2_adv_p3 into two partitions so that prt1_adv_p3 matches both
CREATE TABLE prt2_adv_p3_1 PARTITION OF prt2_adv FOR VALUES FROM (350) TO (375);
CREATE TABLE prt2_adv_p3_2 PARTITION OF prt2_adv FOR VALUES FROM (375) TO (500);
INSERT INTO prt2_adv SELECT i % 25, i, to_char(i, 'FM0000') FROM generate_series(350, 499) i;
ANALYZE prt2_adv;

-- inner join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;

-- semi join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM prt1_adv t1 WHERE EXISTS (SELECT 1 FROM prt2_adv t2 WHERE t1.a = t2.b) AND t1.b = 0 ORDER BY t1.a;

-- left join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 LEFT JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;

-- anti join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM prt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM prt2_adv t2 WHERE t1.a = t2.b) AND t1.b = 0 ORDER BY t1.a;

-- full join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM (SELECT 175 phv, * FROM prt1_adv WHERE prt1_adv.b = 0) t1 FULL JOIN (SELECT 425 phv, * FROM prt2_adv WHERE prt2_adv.a = 0) t2 ON (t1.a = t2.b) WHERE t1.phv = t1.a OR t2.phv = t2.b ORDER BY t1.a, t2.b;

DROP TABLE prt2_adv_p3_1;
DROP TABLE prt2_adv_p3_2;
ANALYZE prt2_adv;

-- Test default partitions
ALTER TABLE prt1_adv DETACH PARTITION prt1_adv_p1;
-- Change prt1_adv_p1 to the default partition
ALTER TABLE prt1_adv ATTACH PARTITION prt1_adv_p1 DEFAULT;
ALTER TABLE prt1_adv DETACH PARTITION prt1_adv_p3;
ANALYZE prt1_adv;

-- We can do partitioned join even if only one of relations has the default
-- partition
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;

-- Restore prt1_adv_p3
ALTER TABLE prt1_adv ATTACH PARTITION prt1_adv_p3 FOR VALUES FROM (300) TO (400);
ANALYZE prt1_adv;

-- Restore prt2_adv_p3
ALTER TABLE prt2_adv ATTACH PARTITION prt2_adv_p3 FOR VALUES FROM (350) TO (500);
ANALYZE prt2_adv;

-- Partitioned join can't be applied because the default partition of prt1_adv
-- matches prt2_adv_p1 and prt2_adv_p3
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;

ALTER TABLE prt2_adv DETACH PARTITION prt2_adv_p3;
-- Change prt2_adv_p3 to the default partition
ALTER TABLE prt2_adv ATTACH PARTITION prt2_adv_p3 DEFAULT;
ANALYZE prt2_adv;

-- Partitioned join can't be applied because the default partition of prt1_adv
-- matches prt2_adv_p1 and prt2_adv_p3
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.b = 0 ORDER BY t1.a, t2.b;

DROP TABLE prt1_adv_p3;
ANALYZE prt1_adv;

DROP TABLE prt2_adv_p3;
ANALYZE prt2_adv;

CREATE TABLE prt3_adv (a int, b int, c varchar) PARTITION BY RANGE (a);
CREATE TABLE prt3_adv_p1 PARTITION OF prt3_adv FOR VALUES FROM (200) TO (300);
CREATE TABLE prt3_adv_p2 PARTITION OF prt3_adv FOR VALUES FROM (300) TO (400);
CREATE INDEX prt3_adv_a_idx ON prt3_adv (a);
INSERT INTO prt3_adv SELECT i, i % 25, to_char(i, 'FM0000') FROM generate_series(200, 399) i;
ANALYZE prt3_adv;

-- 3-way join to test the default partition of a join relation
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c, t3.a, t3.c FROM prt1_adv t1 LEFT JOIN prt2_adv t2 ON (t1.a = t2.b) LEFT JOIN prt3_adv t3 ON (t1.a = t3.a) WHERE t1.b = 0 ORDER BY t1.a, t2.b, t3.a;
SELECT t1.a, t1.c, t2.b, t2.c, t3.a, t3.c FROM prt1_adv t1 LEFT JOIN prt2_adv t2 ON (t1.a = t2.b) LEFT JOIN prt3_adv t3 ON (t1.a = t3.a) WHERE t1.b = 0 ORDER BY t1.a, t2.b, t3.a;

DROP TABLE prt1_adv;
DROP TABLE prt2_adv;
DROP TABLE prt3_adv;

-- Test interaction of partitioned join with partition pruning
CREATE TABLE prt1_adv (a int, b int, c varchar) PARTITION BY RANGE (a);
CREATE TABLE prt1_adv_p1 PARTITION OF prt1_adv FOR VALUES FROM (100) TO (200);
CREATE TABLE prt1_adv_p2 PARTITION OF prt1_adv FOR VALUES FROM (200) TO (300);
CREATE TABLE prt1_adv_p3 PARTITION OF prt1_adv FOR VALUES FROM (300) TO (400);
CREATE INDEX prt1_adv_a_idx ON prt1_adv (a);
INSERT INTO prt1_adv SELECT i, i % 25, to_char(i, 'FM0000') FROM generate_series(100, 399) i;
ANALYZE prt1_adv;

CREATE TABLE prt2_adv (a int, b int, c varchar) PARTITION BY RANGE (b);
CREATE TABLE prt2_adv_p1 PARTITION OF prt2_adv FOR VALUES FROM (100) TO (200);
CREATE TABLE prt2_adv_p2 PARTITION OF prt2_adv FOR VALUES FROM (200) TO (400);
CREATE INDEX prt2_adv_b_idx ON prt2_adv (b);
INSERT INTO prt2_adv SELECT i % 25, i, to_char(i, 'FM0000') FROM generate_series(100, 399) i;
ANALYZE prt2_adv;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.a < 300 AND t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.a < 300 AND t1.b = 0 ORDER BY t1.a, t2.b;

DROP TABLE prt1_adv_p3;
CREATE TABLE prt1_adv_default PARTITION OF prt1_adv DEFAULT;
ANALYZE prt1_adv;

CREATE TABLE prt2_adv_default PARTITION OF prt2_adv DEFAULT;
ANALYZE prt2_adv;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.a >= 100 AND t1.a < 300 AND t1.b = 0 ORDER BY t1.a, t2.b;
SELECT t1.a, t1.c, t2.b, t2.c FROM prt1_adv t1 INNER JOIN prt2_adv t2 ON (t1.a = t2.b) WHERE t1.a >= 100 AND t1.a < 300 AND t1.b = 0 ORDER BY t1.a, t2.b;

DROP TABLE prt1_adv;
DROP TABLE prt2_adv;


-- Tests for list-partitioned tables
CREATE TABLE plt1_adv (a int, b int, c text) PARTITION BY LIST (c);
CREATE TABLE plt1_adv_p1 PARTITION OF plt1_adv FOR VALUES IN ('0001', '0003');
CREATE TABLE plt1_adv_p2 PARTITION OF plt1_adv FOR VALUES IN ('0004', '0006');
CREATE TABLE plt1_adv_p3 PARTITION OF plt1_adv FOR VALUES IN ('0008', '0009');
INSERT INTO plt1_adv SELECT i, i, to_char(i % 10, 'FM0000') FROM generate_series(1, 299) i WHERE i % 10 IN (1, 3, 4, 6, 8, 9);
ANALYZE plt1_adv;

CREATE TABLE plt2_adv (a int, b int, c text) PARTITION BY LIST (c);
CREATE TABLE plt2_adv_p1 PARTITION OF plt2_adv FOR VALUES IN ('0002', '0003');
CREATE TABLE plt2_adv_p2 PARTITION OF plt2_adv FOR VALUES IN ('0004', '0006');
CREATE TABLE plt2_adv_p3 PARTITION OF plt2_adv FOR VALUES IN ('0007', '0009');
INSERT INTO plt2_adv SELECT i, i, to_char(i % 10, 'FM0000') FROM generate_series(1, 299) i WHERE i % 10 IN (2, 3, 4, 6, 7, 9);
ANALYZE plt2_adv;

-- inner join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- semi join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM plt1_adv t1 WHERE EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;
SELECT t1.* FROM plt1_adv t1 WHERE EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;

-- left join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- anti join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM plt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;
SELECT t1.* FROM plt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;

-- full join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 FULL JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE coalesce(t1.b, 0) < 10 AND coalesce(t2.b, 0) < 10 ORDER BY t1.a, t2.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 FULL JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE coalesce(t1.b, 0) < 10 AND coalesce(t2.b, 0) < 10 ORDER BY t1.a, t2.a;

-- Test cases where one side has an extra partition
CREATE TABLE plt2_adv_extra PARTITION OF plt2_adv FOR VALUES IN ('0000');
INSERT INTO plt2_adv_extra VALUES (0, 0, '0000');
ANALYZE plt2_adv;

-- inner join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- semi join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM plt1_adv t1 WHERE EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;
SELECT t1.* FROM plt1_adv t1 WHERE EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;

-- left join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- left join; currently we can't do partitioned join if there are no matched
-- partitions on the nullable side
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt2_adv t1 LEFT JOIN plt1_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- anti join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM plt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;
SELECT t1.* FROM plt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;

-- anti join; currently we can't do partitioned join if there are no matched
-- partitions on the nullable side
EXPLAIN (COSTS OFF)
SELECT t1.* FROM plt2_adv t1 WHERE NOT EXISTS (SELECT 1 FROM plt1_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;

-- full join; currently we can't do partitioned join if there are no matched
-- partitions on the nullable side
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 FULL JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE coalesce(t1.b, 0) < 10 AND coalesce(t2.b, 0) < 10 ORDER BY t1.a, t2.a;

DROP TABLE plt2_adv_extra;

-- Test cases where a partition on one side matches multiple partitions on
-- the other side; we currently can't do partitioned join in such cases
ALTER TABLE plt2_adv DETACH PARTITION plt2_adv_p2;
-- Split plt2_adv_p2 into two partitions so that plt1_adv_p2 matches both
CREATE TABLE plt2_adv_p2_1 PARTITION OF plt2_adv FOR VALUES IN ('0004');
CREATE TABLE plt2_adv_p2_2 PARTITION OF plt2_adv FOR VALUES IN ('0006');
INSERT INTO plt2_adv SELECT i, i, to_char(i % 10, 'FM0000') FROM generate_series(1, 299) i WHERE i % 10 IN (4, 6);
ANALYZE plt2_adv;

-- inner join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- semi join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM plt1_adv t1 WHERE EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;

-- left join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- anti join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM plt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;

-- full join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 FULL JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE coalesce(t1.b, 0) < 10 AND coalesce(t2.b, 0) < 10 ORDER BY t1.a, t2.a;

DROP TABLE plt2_adv_p2_1;
DROP TABLE plt2_adv_p2_2;
-- Restore plt2_adv_p2
ALTER TABLE plt2_adv ATTACH PARTITION plt2_adv_p2 FOR VALUES IN ('0004', '0006');

-- Test NULL partitions
ALTER TABLE plt1_adv DETACH PARTITION plt1_adv_p1;
-- Change plt1_adv_p1 to the NULL partition
CREATE TABLE plt1_adv_p1_null PARTITION OF plt1_adv FOR VALUES IN (NULL, '0001', '0003');
INSERT INTO plt1_adv SELECT i, i, to_char(i % 10, 'FM0000') FROM generate_series(1, 299) i WHERE i % 10 IN (1, 3);
INSERT INTO plt1_adv VALUES (-1, -1, NULL);
ANALYZE plt1_adv;

ALTER TABLE plt2_adv DETACH PARTITION plt2_adv_p3;
-- Change plt2_adv_p3 to the NULL partition
CREATE TABLE plt2_adv_p3_null PARTITION OF plt2_adv FOR VALUES IN (NULL, '0007', '0009');
INSERT INTO plt2_adv SELECT i, i, to_char(i % 10, 'FM0000') FROM generate_series(1, 299) i WHERE i % 10 IN (7, 9);
INSERT INTO plt2_adv VALUES (-1, -1, NULL);
ANALYZE plt2_adv;

-- inner join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- semi join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM plt1_adv t1 WHERE EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;
SELECT t1.* FROM plt1_adv t1 WHERE EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;

-- left join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- anti join
EXPLAIN (COSTS OFF)
SELECT t1.* FROM plt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;
SELECT t1.* FROM plt1_adv t1 WHERE NOT EXISTS (SELECT 1 FROM plt2_adv t2 WHERE t1.a = t2.a AND t1.c = t2.c) AND t1.b < 10 ORDER BY t1.a;

-- full join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 FULL JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE coalesce(t1.b, 0) < 10 AND coalesce(t2.b, 0) < 10 ORDER BY t1.a, t2.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 FULL JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE coalesce(t1.b, 0) < 10 AND coalesce(t2.b, 0) < 10 ORDER BY t1.a, t2.a;

DROP TABLE plt1_adv_p1_null;
-- Restore plt1_adv_p1
ALTER TABLE plt1_adv ATTACH PARTITION plt1_adv_p1 FOR VALUES IN ('0001', '0003');
-- Add to plt1_adv the extra NULL partition containing only NULL values as the
-- key values
CREATE TABLE plt1_adv_extra PARTITION OF plt1_adv FOR VALUES IN (NULL);
INSERT INTO plt1_adv VALUES (-1, -1, NULL);
ANALYZE plt1_adv;

DROP TABLE plt2_adv_p3_null;
-- Restore plt2_adv_p3
ALTER TABLE plt2_adv ATTACH PARTITION plt2_adv_p3 FOR VALUES IN ('0007', '0009');
ANALYZE plt2_adv;

-- inner join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- left join; currently we can't do partitioned join if there are no matched
-- partitions on the nullable side
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- full join; currently we can't do partitioned join if there are no matched
-- partitions on the nullable side
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 FULL JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE coalesce(t1.b, 0) < 10 AND coalesce(t2.b, 0) < 10 ORDER BY t1.a, t2.a;

-- Add to plt2_adv the extra NULL partition containing only NULL values as the
-- key values
CREATE TABLE plt2_adv_extra PARTITION OF plt2_adv FOR VALUES IN (NULL);
INSERT INTO plt2_adv VALUES (-1, -1, NULL);
ANALYZE plt2_adv;

-- inner join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- left join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

-- full join
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 FULL JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE coalesce(t1.b, 0) < 10 AND coalesce(t2.b, 0) < 10 ORDER BY t1.a, t2.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 FULL JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE coalesce(t1.b, 0) < 10 AND coalesce(t2.b, 0) < 10 ORDER BY t1.a, t2.a;

-- 3-way join to test the NULL partition of a join relation
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c, t3.a, t3.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) LEFT JOIN plt1_adv t3 ON (t1.a = t3.a AND t1.c = t3.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c, t3.a, t3.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) LEFT JOIN plt1_adv t3 ON (t1.a = t3.a AND t1.c = t3.c) WHERE t1.b < 10 ORDER BY t1.a;

DROP TABLE plt1_adv_extra;
DROP TABLE plt2_adv_extra;

-- Test default partitions
ALTER TABLE plt1_adv DETACH PARTITION plt1_adv_p1;
-- Change plt1_adv_p1 to the default partition
ALTER TABLE plt1_adv ATTACH PARTITION plt1_adv_p1 DEFAULT;
DROP TABLE plt1_adv_p3;
ANALYZE plt1_adv;

DROP TABLE plt2_adv_p3;
ANALYZE plt2_adv;

-- We can do partitioned join even if only one of relations has the default
-- partition
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

ALTER TABLE plt2_adv DETACH PARTITION plt2_adv_p2;
-- Change plt2_adv_p2 to contain '0005' in addition to '0004' and '0006' as
-- the key values
CREATE TABLE plt2_adv_p2_ext PARTITION OF plt2_adv FOR VALUES IN ('0004', '0005', '0006');
INSERT INTO plt2_adv SELECT i, i, to_char(i % 10, 'FM0000') FROM generate_series(1, 299) i WHERE i % 10 IN (4, 5, 6);
ANALYZE plt2_adv;

-- Partitioned join can't be applied because the default partition of plt1_adv
-- matches plt2_adv_p1 and plt2_adv_p2_ext
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

ALTER TABLE plt2_adv DETACH PARTITION plt2_adv_p2_ext;
-- Change plt2_adv_p2_ext to the default partition
ALTER TABLE plt2_adv ATTACH PARTITION plt2_adv_p2_ext DEFAULT;
ANALYZE plt2_adv;

-- Partitioned join can't be applied because the default partition of plt1_adv
-- matches plt2_adv_p1 and plt2_adv_p2_ext
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

DROP TABLE plt2_adv_p2_ext;
-- Restore plt2_adv_p2
ALTER TABLE plt2_adv ATTACH PARTITION plt2_adv_p2 FOR VALUES IN ('0004', '0006');
ANALYZE plt2_adv;

CREATE TABLE plt3_adv (a int, b int, c text) PARTITION BY LIST (c);
CREATE TABLE plt3_adv_p1 PARTITION OF plt3_adv FOR VALUES IN ('0004', '0006');
CREATE TABLE plt3_adv_p2 PARTITION OF plt3_adv FOR VALUES IN ('0007', '0009');
INSERT INTO plt3_adv SELECT i, i, to_char(i % 10, 'FM0000') FROM generate_series(1, 299) i WHERE i % 10 IN (4, 6, 7, 9);
ANALYZE plt3_adv;

-- 3-way join to test the default partition of a join relation
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c, t3.a, t3.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) LEFT JOIN plt3_adv t3 ON (t1.a = t3.a AND t1.c = t3.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c, t3.a, t3.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) LEFT JOIN plt3_adv t3 ON (t1.a = t3.a AND t1.c = t3.c) WHERE t1.b < 10 ORDER BY t1.a;

-- Test cases where one side has the default partition while the other side
-- has the NULL partition
DROP TABLE plt2_adv_p1;
-- Add the NULL partition to plt2_adv
CREATE TABLE plt2_adv_p1_null PARTITION OF plt2_adv FOR VALUES IN (NULL, '0001', '0003');
INSERT INTO plt2_adv SELECT i, i, to_char(i % 10, 'FM0000') FROM generate_series(1, 299) i WHERE i % 10 IN (1, 3);
INSERT INTO plt2_adv VALUES (-1, -1, NULL);
ANALYZE plt2_adv;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

DROP TABLE plt2_adv_p1_null;
-- Add the NULL partition that contains only NULL values as the key values
CREATE TABLE plt2_adv_p1_null PARTITION OF plt2_adv FOR VALUES IN (NULL);
INSERT INTO plt2_adv VALUES (-1, -1, NULL);
ANALYZE plt2_adv;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.b < 10 ORDER BY t1.a;

DROP TABLE plt1_adv;
DROP TABLE plt2_adv;
DROP TABLE plt3_adv;

-- Test interaction of partitioned join with partition pruning
CREATE TABLE plt1_adv (a int, b int, c text) PARTITION BY LIST (c);
CREATE TABLE plt1_adv_p1 PARTITION OF plt1_adv FOR VALUES IN ('0001');
CREATE TABLE plt1_adv_p2 PARTITION OF plt1_adv FOR VALUES IN ('0002');
CREATE TABLE plt1_adv_p3 PARTITION OF plt1_adv FOR VALUES IN ('0003');
CREATE TABLE plt1_adv_p4 PARTITION OF plt1_adv FOR VALUES IN (NULL, '0004', '0005');
INSERT INTO plt1_adv SELECT i, i, to_char(i % 10, 'FM0000') FROM generate_series(1, 299) i WHERE i % 10 IN (1, 2, 3, 4, 5);
INSERT INTO plt1_adv VALUES (-1, -1, NULL);
ANALYZE plt1_adv;

CREATE TABLE plt2_adv (a int, b int, c text) PARTITION BY LIST (c);
CREATE TABLE plt2_adv_p1 PARTITION OF plt2_adv FOR VALUES IN ('0001', '0002');
CREATE TABLE plt2_adv_p2 PARTITION OF plt2_adv FOR VALUES IN (NULL);
CREATE TABLE plt2_adv_p3 PARTITION OF plt2_adv FOR VALUES IN ('0003');
CREATE TABLE plt2_adv_p4 PARTITION OF plt2_adv FOR VALUES IN ('0004', '0005');
INSERT INTO plt2_adv SELECT i, i, to_char(i % 10, 'FM0000') FROM generate_series(1, 299) i WHERE i % 10 IN (1, 2, 3, 4, 5);
INSERT INTO plt2_adv VALUES (-1, -1, NULL);
ANALYZE plt2_adv;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.c IN ('0003', '0004', '0005') AND t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.c IN ('0003', '0004', '0005') AND t1.b < 10 ORDER BY t1.a;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.c IS NULL AND t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.c IS NULL AND t1.b < 10 ORDER BY t1.a;

CREATE TABLE plt1_adv_default PARTITION OF plt1_adv DEFAULT;
ANALYZE plt1_adv;

CREATE TABLE plt2_adv_default PARTITION OF plt2_adv DEFAULT;
ANALYZE plt2_adv;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.c IN ('0003', '0004', '0005') AND t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 INNER JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.c IN ('0003', '0004', '0005') AND t1.b < 10 ORDER BY t1.a;

EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.c IS NULL AND t1.b < 10 ORDER BY t1.a;
SELECT t1.a, t1.c, t2.a, t2.c FROM plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE t1.c IS NULL AND t1.b < 10 ORDER BY t1.a;

DROP TABLE plt1_adv;
DROP TABLE plt2_adv;

-- Test the process_outer_partition() code path
CREATE TABLE plt1_adv (a int, b int, c text) PARTITION BY LIST (c);
CREATE TABLE plt1_adv_p1 PARTITION OF plt1_adv FOR VALUES IN ('0000', '0001', '0002');
CREATE TABLE plt1_adv_p2 PARTITION OF plt1_adv FOR VALUES IN ('0003', '0004');
INSERT INTO plt1_adv SELECT i, i, to_char(i % 5, 'FM0000') FROM generate_series(0, 24) i;
ANALYZE plt1_adv;

CREATE TABLE plt2_adv (a int, b int, c text) PARTITION BY LIST (c);
CREATE TABLE plt2_adv_p1 PARTITION OF plt2_adv FOR VALUES IN ('0002');
CREATE TABLE plt2_adv_p2 PARTITION OF plt2_adv FOR VALUES IN ('0003', '0004');
INSERT INTO plt2_adv SELECT i, i, to_char(i % 5, 'FM0000') FROM generate_series(0, 24) i WHERE i % 5 IN (2, 3, 4);
ANALYZE plt2_adv;

CREATE TABLE plt3_adv (a int, b int, c text) PARTITION BY LIST (c);
CREATE TABLE plt3_adv_p1 PARTITION OF plt3_adv FOR VALUES IN ('0001');
CREATE TABLE plt3_adv_p2 PARTITION OF plt3_adv FOR VALUES IN ('0003', '0004');
INSERT INTO plt3_adv SELECT i, i, to_char(i % 5, 'FM0000') FROM generate_series(0, 24) i WHERE i % 5 IN (1, 3, 4);
ANALYZE plt3_adv;

-- This tests that when merging partitions from plt1_adv and plt2_adv in
-- merge_list_bounds(), process_outer_partition() returns an already-assigned
-- merged partition when re-called with plt1_adv_p1 for the second list value
-- '0001' of that partition
EXPLAIN (COSTS OFF)
SELECT t1.a, t1.c, t2.a, t2.c, t3.a, t3.c FROM (plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.c = t2.c)) FULL JOIN plt3_adv t3 ON (t1.c = t3.c) WHERE coalesce(t1.a, 0) % 5 != 3 AND coalesce(t1.a, 0) % 5 != 4 ORDER BY t1.c, t1.a, t2.a, t3.a;
SELECT t1.a, t1.c, t2.a, t2.c, t3.a, t3.c FROM (plt1_adv t1 LEFT JOIN plt2_adv t2 ON (t1.c = t2.c)) FULL JOIN plt3_adv t3 ON (t1.c = t3.c) WHERE coalesce(t1.a, 0) % 5 != 3 AND coalesce(t1.a, 0) % 5 != 4 ORDER BY t1.c, t1.a, t2.a, t3.a;

DROP TABLE plt1_adv;
DROP TABLE plt2_adv;
DROP TABLE plt3_adv;


-- Tests for multi-level partitioned tables
CREATE TABLE alpha (a double precision, b int, c text) PARTITION BY RANGE (a);
CREATE TABLE alpha_neg PARTITION OF alpha FOR VALUES FROM ('-Infinity') TO (0) PARTITION BY RANGE (b);
CREATE TABLE alpha_pos PARTITION OF alpha FOR VALUES FROM (0) TO (10.0) PARTITION BY LIST (c);
CREATE TABLE alpha_neg_p1 PARTITION OF alpha_neg FOR VALUES FROM (100) TO (200);
CREATE TABLE alpha_neg_p2 PARTITION OF alpha_neg FOR VALUES FROM (200) TO (300);
CREATE TABLE alpha_neg_p3 PARTITION OF alpha_neg FOR VALUES FROM (300) TO (400);
CREATE TABLE alpha_pos_p1 PARTITION OF alpha_pos FOR VALUES IN ('0001', '0003');
CREATE TABLE alpha_pos_p2 PARTITION OF alpha_pos FOR VALUES IN ('0004', '0006');
CREATE TABLE alpha_pos_p3 PARTITION OF alpha_pos FOR VALUES IN ('0008', '0009');
INSERT INTO alpha_neg SELECT -1.0, i, to_char(i % 10, 'FM0000') FROM generate_series(100, 399) i WHERE i % 10 IN (1, 3, 4, 6, 8, 9);
INSERT INTO alpha_pos SELECT  1.0, i, to_char(i % 10, 'FM0000') FROM generate_series(100, 399) i WHERE i % 10 IN (1, 3, 4, 6, 8, 9);
ANALYZE alpha;

CREATE TABLE beta (a double precision, b int, c text) PARTITION BY RANGE (a);
CREATE TABLE beta_neg PARTITION OF beta FOR VALUES FROM (-10.0) TO (0) PARTITION BY RANGE (b);
CREATE TABLE beta_pos PARTITION OF beta FOR VALUES FROM (0) TO ('Infinity') PARTITION BY LIST (c);
CREATE TABLE beta_neg_p1 PARTITION OF beta_neg FOR VALUES FROM (100) TO (150);
CREATE TABLE beta_neg_p2 PARTITION OF beta_neg FOR VALUES FROM (200) TO (300);
CREATE TABLE beta_neg_p3 PARTITION OF beta_neg FOR VALUES FROM (350) TO (500);
CREATE TABLE beta_pos_p1 PARTITION OF beta_pos FOR VALUES IN ('0002', '0003');
CREATE TABLE beta_pos_p2 PARTITION OF beta_pos FOR VALUES IN ('0004', '0006');
CREATE TABLE beta_pos_p3 PARTITION OF beta_pos FOR VALUES IN ('0007', '0009');
INSERT INTO beta_neg SELECT -1.0, i, to_char(i % 10, 'FM0000') FROM generate_series(100, 149) i WHERE i % 10 IN (2, 3, 4, 6, 7, 9);
INSERT INTO beta_neg SELECT -1.0, i, to_char(i % 10, 'FM0000') FROM generate_series(200, 299) i WHERE i % 10 IN (2, 3, 4, 6, 7, 9);
INSERT INTO beta_neg SELECT -1.0, i, to_char(i % 10, 'FM0000') FROM generate_series(350, 499) i WHERE i % 10 IN (2, 3, 4, 6, 7, 9);
INSERT INTO beta_pos SELECT  1.0, i, to_char(i % 10, 'FM0000') FROM generate_series(100, 149) i WHERE i % 10 IN (2, 3, 4, 6, 7, 9);
INSERT INTO beta_pos SELECT  1.0, i, to_char(i % 10, 'FM0000') FROM generate_series(200, 299) i WHERE i % 10 IN (2, 3, 4, 6, 7, 9);
INSERT INTO beta_pos SELECT  1.0, i, to_char(i % 10, 'FM0000') FROM generate_series(350, 499) i WHERE i % 10 IN (2, 3, 4, 6, 7, 9);
ANALYZE beta;

EXPLAIN (COSTS OFF)
SELECT t1.*, t2.* FROM alpha t1 INNER JOIN beta t2 ON (t1.a = t2.a AND t1.b = t2.b) WHERE t1.b >= 125 AND t1.b < 225 ORDER BY t1.a, t1.b;
SELECT t1.*, t2.* FROM alpha t1 INNER JOIN beta t2 ON (t1.a = t2.a AND t1.b = t2.b) WHERE t1.b >= 125 AND t1.b < 225 ORDER BY t1.a, t1.b;

EXPLAIN (COSTS OFF)
SELECT t1.*, t2.* FROM alpha t1 INNER JOIN beta t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE ((t1.b >= 100 AND t1.b < 110) OR (t1.b >= 200 AND t1.b < 210)) AND ((t2.b >= 100 AND t2.b < 110) OR (t2.b >= 200 AND t2.b < 210)) AND t1.c IN ('0004', '0009') ORDER BY t1.a, t1.b, t2.b;
SELECT t1.*, t2.* FROM alpha t1 INNER JOIN beta t2 ON (t1.a = t2.a AND t1.c = t2.c) WHERE ((t1.b >= 100 AND t1.b < 110) OR (t1.b >= 200 AND t1.b < 210)) AND ((t2.b >= 100 AND t2.b < 110) OR (t2.b >= 200 AND t2.b < 210)) AND t1.c IN ('0004', '0009') ORDER BY t1.a, t1.b, t2.b;

EXPLAIN (COSTS OFF)
SELECT t1.*, t2.* FROM alpha t1 INNER JOIN beta t2 ON (t1.a = t2.a AND t1.b = t2.b AND t1.c = t2.c) WHERE ((t1.b >= 100 AND t1.b < 110) OR (t1.b >= 200 AND t1.b < 210)) AND ((t2.b >= 100 AND t2.b < 110) OR (t2.b >= 200 AND t2.b < 210)) AND t1.c IN ('0004', '0009') ORDER BY t1.a, t1.b;
SELECT t1.*, t2.* FROM alpha t1 INNER JOIN beta t2 ON (t1.a = t2.a AND t1.b = t2.b AND t1.c = t2.c) WHERE ((t1.b >= 100 AND t1.b < 110) OR (t1.b >= 200 AND t1.b < 210)) AND ((t2.b >= 100 AND t2.b < 110) OR (t2.b >= 200 AND t2.b < 210)) AND t1.c IN ('0004', '0009') ORDER BY t1.a, t1.b;

-- partitionwise join with fractional paths
CREATE TABLE fract_t (id BIGINT, PRIMARY KEY (id)) PARTITION BY RANGE (id);
CREATE TABLE fract_t0 PARTITION OF fract_t FOR VALUES FROM ('0') TO ('1000');
CREATE TABLE fract_t1 PARTITION OF fract_t FOR VALUES FROM ('1000') TO ('2000');

-- insert data
INSERT INTO fract_t (id) (SELECT generate_series(0, 1999));
ANALYZE fract_t;

-- verify plan; nested index only scans
SET max_parallel_workers_per_gather = 0;
SET enable_partitionwise_join = on;

EXPLAIN (COSTS OFF)
SELECT x.id, y.id FROM fract_t x LEFT JOIN fract_t y USING (id) ORDER BY x.id ASC LIMIT 10;

EXPLAIN (COSTS OFF)
SELECT x.id, y.id FROM fract_t x LEFT JOIN fract_t y USING (id) ORDER BY x.id DESC LIMIT 10;

-- cleanup
DROP TABLE fract_t;

RESET max_parallel_workers_per_gather;
RESET enable_partitionwise_join;

--
-- exercises for the hash join code
--

begin;

set local min_parallel_table_scan_size = 0;
set local parallel_setup_cost = 0;
set local enable_hashjoin = on;

-- Extract bucket and batch counts from an explain analyze plan.  In
-- general we can't make assertions about how many batches (or
-- buckets) will be required because it can vary, but we can in some
-- special cases and we can check for growth.
create or replace function find_hash(node json)
returns json language plpgsql
as
$$
declare
  x json;
  child json;
begin
  if node->>'Node Type' = 'Hash' then
    return node;
  else
    for child in select json_array_elements(node->'Plans')
    loop
      x := find_hash(child);
      if x is not null then
        return x;
      end if;
    end loop;
    return null;
  end if;
end;
$$;
create or replace function hash_join_batches(query text)
returns table (original int, final int) language plpgsql
as
$$
declare
  whole_plan json;
  hash_node json;
begin
  for whole_plan in
    execute 'explain (analyze, format ''json'') ' || query
  loop
    hash_node := find_hash(json_extract_path(whole_plan, '0', 'Plan'));
    original := hash_node->>'Original Hash Batches';
    final := hash_node->>'Hash Batches';
    return next;
  end loop;
end;
$$;

-- Make a simple relation with well distributed keys and correctly
-- estimated size.
create table simple as
  select generate_series(1, 20000) AS id, 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa';
alter table simple set (parallel_workers = 2);
analyze simple;

-- Make a relation whose size we will under-estimate.  We want stats
-- to say 1000 rows, but actually there are 20,000 rows.
create table bigger_than_it_looks as
  select generate_series(1, 20000) as id, 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa';
alter table bigger_than_it_looks set (autovacuum_enabled = 'false');
alter table bigger_than_it_looks set (parallel_workers = 2);
analyze bigger_than_it_looks;
update pg_class set reltuples = 1000 where relname = 'bigger_than_it_looks';

-- Make a relation whose size we underestimate and that also has a
-- kind of skew that breaks our batching scheme.  We want stats to say
-- 2 rows, but actually there are 20,000 rows with the same key.
create table extremely_skewed (id int, t text);
alter table extremely_skewed set (autovacuum_enabled = 'false');
alter table extremely_skewed set (parallel_workers = 2);
analyze extremely_skewed;
insert into extremely_skewed
  select 42 as id, 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa'
  from generate_series(1, 20000);
update pg_class
  set reltuples = 2, relpages = pg_relation_size('extremely_skewed') / 8192
  where relname = 'extremely_skewed';

-- Make a relation with a couple of enormous tuples.
create table wide as select generate_series(1, 2) as id, rpad('', 320000, 'x') as t;
alter table wide set (parallel_workers = 2);

-- The "optimal" case: the hash table fits in memory; we plan for 1
-- batch, we stick to that number, and peak memory usage stays within
-- our work_mem budget

-- non-parallel
savepoint settings;
set local max_parallel_workers_per_gather = 0;
set local work_mem = '4MB';
set local hash_mem_multiplier = 1.0;
explain (costs off)
  select count(*) from simple r join simple s using (id);
select count(*) from simple r join simple s using (id);
select original > 1 as initially_multibatch, final > original as increased_batches
  from hash_join_batches(
$$
  select count(*) from simple r join simple s using (id);
$$);
rollback to settings;

-- parallel with parallel-oblivious hash join
savepoint settings;
set local max_parallel_workers_per_gather = 2;
set local work_mem = '4MB';
set local hash_mem_multiplier = 1.0;
set local enable_parallel_hash = off;
explain (costs off)
  select count(*) from simple r join simple s using (id);
select count(*) from simple r join simple s using (id);
select original > 1 as initially_multibatch, final > original as increased_batches
  from hash_join_batches(
$$
  select count(*) from simple r join simple s using (id);
$$);
rollback to settings;

-- parallel with parallel-aware hash join
savepoint settings;
set local max_parallel_workers_per_gather = 2;
set local work_mem = '4MB';
set local hash_mem_multiplier = 1.0;
set local enable_parallel_hash = on;
explain (costs off)
  select count(*) from simple r join simple s using (id);
select count(*) from simple r join simple s using (id);
select original > 1 as initially_multibatch, final > original as increased_batches
  from hash_join_batches(
$$
  select count(*) from simple r join simple s using (id);
$$);
rollback to settings;

-- The "good" case: batches required, but we plan the right number; we
-- plan for some number of batches, and we stick to that number, and
-- peak memory usage says within our work_mem budget

-- non-parallel
savepoint settings;
set local max_parallel_workers_per_gather = 0;
set local work_mem = '128kB';
set local hash_mem_multiplier = 1.0;
explain (costs off)
  select count(*) from simple r join simple s using (id);
select count(*) from simple r join simple s using (id);
select original > 1 as initially_multibatch, final > original as increased_batches
  from hash_join_batches(
$$
  select count(*) from simple r join simple s using (id);
$$);
rollback to settings;

-- parallel with parallel-oblivious hash join
savepoint settings;
set local max_parallel_workers_per_gather = 2;
set local work_mem = '128kB';
set local hash_mem_multiplier = 1.0;
set local enable_parallel_hash = off;
explain (costs off)
  select count(*) from simple r join simple s using (id);
select count(*) from simple r join simple s using (id);
select original > 1 as initially_multibatch, final > original as increased_batches
  from hash_join_batches(
$$
  select count(*) from simple r join simple s using (id);
$$);
rollback to settings;

-- parallel with parallel-aware hash join
savepoint settings;
set local max_parallel_workers_per_gather = 2;
set local work_mem = '192kB';
set local hash_mem_multiplier = 1.0;
set local enable_parallel_hash = on;
explain (costs off)
  select count(*) from simple r join simple s using (id);
select count(*) from simple r join simple s using (id);
select original > 1 as initially_multibatch, final > original as increased_batches
  from hash_join_batches(
$$
  select count(*) from simple r join simple s using (id);
$$);
-- parallel full multi-batch hash join
select count(*) from simple r full outer join simple s using (id);
rollback to settings;

-- The "bad" case: during execution we need to increase number of
-- batches; in this case we plan for 1 batch, and increase at least a
-- couple of times, and peak memory usage stays within our work_mem
-- budget

-- non-parallel
savepoint settings;
set local max_parallel_workers_per_gather = 0;
set local work_mem = '128kB';
set local hash_mem_multiplier = 1.0;
explain (costs off)
  select count(*) FROM simple r JOIN bigger_than_it_looks s USING (id);
select count(*) FROM simple r JOIN bigger_than_it_looks s USING (id);
select original > 1 as initially_multibatch, final > original as increased_batches
  from hash_join_batches(
$$
  select count(*) FROM simple r JOIN bigger_than_it_looks s USING (id);
$$);
rollback to settings;

-- parallel with parallel-oblivious hash join
savepoint settings;
set local max_parallel_workers_per_gather = 2;
set local work_mem = '128kB';
set local hash_mem_multiplier = 1.0;
set local enable_parallel_hash = off;
explain (costs off)
  select count(*) from simple r join bigger_than_it_looks s using (id);
select count(*) from simple r join bigger_than_it_looks s using (id);
select original > 1 as initially_multibatch, final > original as increased_batches
  from hash_join_batches(
$$
  select count(*) from simple r join bigger_than_it_looks s using (id);
$$);
rollback to settings;

-- parallel with parallel-aware hash join
savepoint settings;
set local max_parallel_workers_per_gather = 1;
set local work_mem = '192kB';
set local hash_mem_multiplier = 1.0;
set local enable_parallel_hash = on;
explain (costs off)
  select count(*) from simple r join bigger_than_it_looks s using (id);
select count(*) from simple r join bigger_than_it_looks s using (id);
select original > 1 as initially_multibatch, final > original as increased_batches
  from hash_join_batches(
$$
  select count(*) from simple r join bigger_than_it_looks s using (id);
$$);
rollback to settings;

-- The "ugly" case: increasing the number of batches during execution
-- doesn't help, so stop trying to fit in work_mem and hope for the
-- best; in this case we plan for 1 batch, increases just once and
-- then stop increasing because that didn't help at all, so we blow
-- right through the work_mem budget and hope for the best...

-- non-parallel
savepoint settings;
set local max_parallel_workers_per_gather = 0;
set local work_mem = '128kB';
set local hash_mem_multiplier = 1.0;
explain (costs off)
  select count(*) from simple r join extremely_skewed s using (id);
select count(*) from simple r join extremely_skewed s using (id);
select * from hash_join_batches(
$$
  select count(*) from simple r join extremely_skewed s using (id);
$$);
rollback to settings;

-- parallel with parallel-oblivious hash join
savepoint settings;
set local max_parallel_workers_per_gather = 2;
set local work_mem = '128kB';
set local hash_mem_multiplier = 1.0;
set local enable_parallel_hash = off;
explain (costs off)
  select count(*) from simple r join extremely_skewed s using (id);
select count(*) from simple r join extremely_skewed s using (id);
select * from hash_join_batches(
$$
  select count(*) from simple r join extremely_skewed s using (id);
$$);
rollback to settings;

-- parallel with parallel-aware hash join
savepoint settings;
set local max_parallel_workers_per_gather = 1;
set local work_mem = '128kB';
set local hash_mem_multiplier = 1.0;
set local enable_parallel_hash = on;
explain (costs off)
  select count(*) from simple r join extremely_skewed s using (id);
select count(*) from simple r join extremely_skewed s using (id);
select * from hash_join_batches(
$$
  select count(*) from simple r join extremely_skewed s using (id);
$$);
rollback to settings;

-- A couple of other hash join tests unrelated to work_mem management.

-- Check that EXPLAIN ANALYZE has data even if the leader doesn't participate
savepoint settings;
set local max_parallel_workers_per_gather = 2;
set local work_mem = '4MB';
set local hash_mem_multiplier = 1.0;
set local parallel_leader_participation = off;
select * from hash_join_batches(
$$
  select count(*) from simple r join simple s using (id);
$$);
rollback to settings;

-- Exercise rescans.  We'll turn off parallel_leader_participation so
-- that we can check that instrumentation comes back correctly.

create table join_foo as select generate_series(1, 3) as id, 'xxxxx'::text as t;
alter table join_foo set (parallel_workers = 0);
create table join_bar as select generate_series(1, 10000) as id, 'xxxxx'::text as t;
alter table join_bar set (parallel_workers = 2);

-- multi-batch with rescan, parallel-oblivious
savepoint settings;
set enable_parallel_hash = off;
set parallel_leader_participation = off;
set min_parallel_table_scan_size = 0;
set parallel_setup_cost = 0;
set parallel_tuple_cost = 0;
set max_parallel_workers_per_gather = 2;
set enable_material = off;
set enable_mergejoin = off;
set work_mem = '64kB';
set hash_mem_multiplier = 1.0;
explain (costs off)
  select count(*) from join_foo
    left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
    on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
select count(*) from join_foo
  left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
  on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
select final > 1 as multibatch
  from hash_join_batches(
$$
  select count(*) from join_foo
    left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
    on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
$$);
rollback to settings;

-- single-batch with rescan, parallel-oblivious
savepoint settings;
set enable_parallel_hash = off;
set parallel_leader_participation = off;
set min_parallel_table_scan_size = 0;
set parallel_setup_cost = 0;
set parallel_tuple_cost = 0;
set max_parallel_workers_per_gather = 2;
set enable_material = off;
set enable_mergejoin = off;
set work_mem = '4MB';
set hash_mem_multiplier = 1.0;
explain (costs off)
  select count(*) from join_foo
    left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
    on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
select count(*) from join_foo
  left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
  on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
select final > 1 as multibatch
  from hash_join_batches(
$$
  select count(*) from join_foo
    left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
    on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
$$);
rollback to settings;

-- multi-batch with rescan, parallel-aware
savepoint settings;
set enable_parallel_hash = on;
set parallel_leader_participation = off;
set min_parallel_table_scan_size = 0;
set parallel_setup_cost = 0;
set parallel_tuple_cost = 0;
set max_parallel_workers_per_gather = 2;
set enable_material = off;
set enable_mergejoin = off;
set work_mem = '64kB';
set hash_mem_multiplier = 1.0;
explain (costs off)
  select count(*) from join_foo
    left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
    on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
select count(*) from join_foo
  left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
  on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
select final > 1 as multibatch
  from hash_join_batches(
$$
  select count(*) from join_foo
    left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
    on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
$$);
rollback to settings;

-- single-batch with rescan, parallel-aware
savepoint settings;
set enable_parallel_hash = on;
set parallel_leader_participation = off;
set min_parallel_table_scan_size = 0;
set parallel_setup_cost = 0;
set parallel_tuple_cost = 0;
set max_parallel_workers_per_gather = 2;
set enable_material = off;
set enable_mergejoin = off;
set work_mem = '4MB';
set hash_mem_multiplier = 1.0;
explain (costs off)
  select count(*) from join_foo
    left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
    on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
select count(*) from join_foo
  left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
  on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
select final > 1 as multibatch
  from hash_join_batches(
$$
  select count(*) from join_foo
    left join (select b1.id, b1.t from join_bar b1 join join_bar b2 using (id)) ss
    on join_foo.id < ss.id + 1 and join_foo.id > ss.id - 1;
$$);
rollback to settings;

-- A full outer join where every record is matched.

-- non-parallel
savepoint settings;
set local max_parallel_workers_per_gather = 0;
explain (costs off)
     select  count(*) from simple r full outer join simple s using (id);
select  count(*) from simple r full outer join simple s using (id);
rollback to settings;

-- parallelism not possible with parallel-oblivious full hash join
savepoint settings;
set enable_parallel_hash = off;
set local max_parallel_workers_per_gather = 2;
explain (costs off)
     select  count(*) from simple r full outer join simple s using (id);
select  count(*) from simple r full outer join simple s using (id);
rollback to settings;

-- parallelism is possible with parallel-aware full hash join
savepoint settings;
set local max_parallel_workers_per_gather = 2;
explain (costs off)
     select  count(*) from simple r full outer join simple s using (id);
select  count(*) from simple r full outer join simple s using (id);
rollback to settings;

-- A full outer join where every record is not matched.

-- non-parallel
savepoint settings;
set local max_parallel_workers_per_gather = 0;
explain (costs off)
     select  count(*) from simple r full outer join simple s on (r.id = 0 - s.id);
select  count(*) from simple r full outer join simple s on (r.id = 0 - s.id);
rollback to settings;

-- parallelism not possible with parallel-oblivious full hash join
savepoint settings;
set enable_parallel_hash = off;
set local max_parallel_workers_per_gather = 2;
explain (costs off)
     select  count(*) from simple r full outer join simple s on (r.id = 0 - s.id);
select  count(*) from simple r full outer join simple s on (r.id = 0 - s.id);
rollback to settings;

-- parallelism is possible with parallel-aware full hash join
savepoint settings;
set local max_parallel_workers_per_gather = 2;
explain (costs off)
     select  count(*) from simple r full outer join simple s on (r.id = 0 - s.id);
select  count(*) from simple r full outer join simple s on (r.id = 0 - s.id);
rollback to settings;


-- exercise special code paths for huge tuples (note use of non-strict
-- expression and left join required to get the detoasted tuple into
-- the hash table)

-- parallel with parallel-aware hash join (hits ExecParallelHashLoadTuple and
-- sts_puttuple oversized tuple cases because it's multi-batch)
savepoint settings;
set max_parallel_workers_per_gather = 2;
set enable_parallel_hash = on;
set work_mem = '128kB';
set hash_mem_multiplier = 1.0;
explain (costs off)
  select length(max(s.t))
  from wide left join (select id, coalesce(t, '') || '' as t from wide) s using (id);
select length(max(s.t))
from wide left join (select id, coalesce(t, '') || '' as t from wide) s using (id);
select final > 1 as multibatch
  from hash_join_batches(
$$
  select length(max(s.t))
  from wide left join (select id, coalesce(t, '') || '' as t from wide) s using (id);
$$);
rollback to settings;


-- Hash join reuses the HOT status bit to indicate match status. This can only
-- be guaranteed to produce correct results if all the hash join tuple match
-- bits are reset before reuse. This is done upon loading them into the
-- hashtable.
SAVEPOINT settings;
SET enable_parallel_hash = on;
SET min_parallel_table_scan_size = 0;
SET parallel_setup_cost = 0;
SET parallel_tuple_cost = 0;
CREATE TABLE hjtest_matchbits_t1(id int);
CREATE TABLE hjtest_matchbits_t2(id int);
INSERT INTO hjtest_matchbits_t1 VALUES (1);
INSERT INTO hjtest_matchbits_t2 VALUES (2);
-- Update should create a HOT tuple. If this status bit isn't cleared, we won't
-- correctly emit the NULL-extended unmatching tuple in full hash join.
UPDATE hjtest_matchbits_t2 set id = 2;
SELECT * FROM hjtest_matchbits_t1 t1 FULL JOIN hjtest_matchbits_t2 t2 ON t1.id = t2.id
  ORDER BY t1.id;
-- Test serial full hash join.
-- Resetting parallel_setup_cost should force a serial plan.
-- Just to be safe, however, set enable_parallel_hash to off, as parallel full
-- hash joins are only supported with shared hashtables.
RESET parallel_setup_cost;
SET enable_parallel_hash = off;
SELECT * FROM hjtest_matchbits_t1 t1 FULL JOIN hjtest_matchbits_t2 t2 ON t1.id = t2.id;
ROLLBACK TO settings;

rollback;

-- Verify that hash key expressions reference the correct
-- nodes. Hashjoin's hashkeys need to reference its outer plan, Hash's
-- need to reference Hash's outer plan (which is below HashJoin's
-- inner plan). It's not trivial to verify that the references are
-- correct (we don't display the hashkeys themselves), but if the
-- hashkeys contain subplan references, those will be displayed. Force
-- subplans to appear just about everywhere.
--
-- Bug report:
-- https://www.postgresql.org/message-id/CAPpHfdvGVegF_TKKRiBrSmatJL2dR9uwFCuR%2BteQ_8tEXU8mxg%40mail.gmail.com
--
BEGIN;
SET LOCAL enable_sort = OFF; -- avoid mergejoins
SET LOCAL from_collapse_limit = 1; -- allows easy changing of join order

CREATE TABLE hjtest_1 (a text, b int, id int, c bool);
CREATE TABLE hjtest_2 (a bool, id int, b text, c int);

INSERT INTO hjtest_1(a, b, id, c) VALUES ('text', 2, 1, false); -- matches
INSERT INTO hjtest_1(a, b, id, c) VALUES ('text', 1, 2, false); -- fails id join condition
INSERT INTO hjtest_1(a, b, id, c) VALUES ('text', 20, 1, false); -- fails < 50
INSERT INTO hjtest_1(a, b, id, c) VALUES ('text', 1, 1, false); -- fails (SELECT hjtest_1.b * 5) = (SELECT hjtest_2.c*5)

INSERT INTO hjtest_2(a, id, b, c) VALUES (true, 1, 'another', 2); -- matches
INSERT INTO hjtest_2(a, id, b, c) VALUES (true, 3, 'another', 7); -- fails id join condition
INSERT INTO hjtest_2(a, id, b, c) VALUES (true, 1, 'another', 90);  -- fails < 55
INSERT INTO hjtest_2(a, id, b, c) VALUES (true, 1, 'another', 3); -- fails (SELECT hjtest_1.b * 5) = (SELECT hjtest_2.c*5)
INSERT INTO hjtest_2(a, id, b, c) VALUES (true, 1, 'text', 1); --  fails hjtest_1.a <> hjtest_2.b;

EXPLAIN (COSTS OFF, VERBOSE)
SELECT hjtest_1.a a1, hjtest_2.a a2,hjtest_1.tableoid::regclass t1, hjtest_2.tableoid::regclass t2
FROM hjtest_1, hjtest_2
WHERE
    hjtest_1.id = (SELECT 1 WHERE hjtest_2.id = 1)
    AND (SELECT hjtest_1.b * 5) = (SELECT hjtest_2.c*5)
    AND (SELECT hjtest_1.b * 5) < 50
    AND (SELECT hjtest_2.c * 5) < 55
    AND hjtest_1.a <> hjtest_2.b;

SELECT hjtest_1.a a1, hjtest_2.a a2,hjtest_1.tableoid::regclass t1, hjtest_2.tableoid::regclass t2
FROM hjtest_1, hjtest_2
WHERE
    hjtest_1.id = (SELECT 1 WHERE hjtest_2.id = 1)
    AND (SELECT hjtest_1.b * 5) = (SELECT hjtest_2.c*5)
    AND (SELECT hjtest_1.b * 5) < 50
    AND (SELECT hjtest_2.c * 5) < 55
    AND hjtest_1.a <> hjtest_2.b;

EXPLAIN (COSTS OFF, VERBOSE)
SELECT hjtest_1.a a1, hjtest_2.a a2,hjtest_1.tableoid::regclass t1, hjtest_2.tableoid::regclass t2
FROM hjtest_2, hjtest_1
WHERE
    hjtest_1.id = (SELECT 1 WHERE hjtest_2.id = 1)
    AND (SELECT hjtest_1.b * 5) = (SELECT hjtest_2.c*5)
    AND (SELECT hjtest_1.b * 5) < 50
    AND (SELECT hjtest_2.c * 5) < 55
    AND hjtest_1.a <> hjtest_2.b;

SELECT hjtest_1.a a1, hjtest_2.a a2,hjtest_1.tableoid::regclass t1, hjtest_2.tableoid::regclass t2
FROM hjtest_2, hjtest_1
WHERE
    hjtest_1.id = (SELECT 1 WHERE hjtest_2.id = 1)
    AND (SELECT hjtest_1.b * 5) = (SELECT hjtest_2.c*5)
    AND (SELECT hjtest_1.b * 5) < 50
    AND (SELECT hjtest_2.c * 5) < 55
    AND hjtest_1.a <> hjtest_2.b;

ROLLBACK;

-- Verify that we behave sanely when the inner hash keys contain parameters
-- (that is, outer or lateral references).  This situation has to defeat
-- re-use of the inner hash table across rescans.
begin;
set local enable_hashjoin = on;

explain (costs off)
select i8.q2, ss.* from
int8_tbl i8,
lateral (select t1.fivethous, i4.f1 from tenk1 t1 join int4_tbl i4
         on t1.fivethous = i4.f1+i8.q2 order by 1,2) ss;

select i8.q2, ss.* from
int8_tbl i8,
lateral (select t1.fivethous, i4.f1 from tenk1 t1 join int4_tbl i4
         on t1.fivethous = i4.f1+i8.q2 order by 1,2) ss;

rollback;
