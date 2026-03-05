WITH total_distinct_ko AS 
(SELECT COUNT(*) AS total_ko FROM (SELECT DISTINCT KEGG_ko
FROM pseudomonas.main.pseudomonas_annotations)
),
count_by_group AS 
(SELECT animal_group, COUNT(*) AS total_by_hosts FROM (SELECT DISTINCT animal_group, KEGG_ko
FROM pseudomonas.main.pseudomonas_annotations)
GROUP BY animal_group
UNION 
SELECT 'Universe' AS animal_group, COUNT(*) AS total_ko FROM (SELECT DISTINCT KEGG_ko
FROM pseudomonas.main.pseudomonas_annotations))
SELECT cbg.animal_group, cbg.total_by_hosts, tdk.total_ko, cbg.total_by_hosts/tdk.total_ko AS percent_of_total
FROM count_by_group cbg, total_distinct_ko tdk
ORDER BY animal_group;