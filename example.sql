EXPLAIN ANALYZE
SELECT * 
FROM t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14
WHERE 
    -- {t1, t2} - {t3, t4}
    t1.x + t2.x > t3.x + t4.x AND
    -- {t3, t7, t8, t9} - {t1, t4, t5}
    t3.x + t7.x + t8.x + t9.x > t1.x + t4.x + t5.x AND
    -- {t1, t10, t11} - {t2, t3, t4, t5, t6, t7, t8, t9, t12, t13, t14}
    t1.x + t10.x + t11.x > t2.x + t3.x + t4.x + t5.x + t6.x + t7.x + t8.x + t9.x + t12.x + t13.x + t14.x AND
    -- {t5} - {t1, t8, t9, t10}
    t5.x > t1.x + t8.x + t9.x + t10.x AND
    -- {t10, t11, t12, t13, t14} - {t5, t6, t7, t8}
    t10.x + t11.x + t12.x + t13.x + t14.x > t5.x + t6.x + t7.x + t8.x;