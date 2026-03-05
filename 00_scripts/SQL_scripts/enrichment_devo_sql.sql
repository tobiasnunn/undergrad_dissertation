WITH distinct_animal_group AS 
(SELECT DISTINCT Accession, animal_group FROM
pseudomonas_annotations),
animal_group_count AS
(SELECT animal_group, COUNT(*) AS total_accessions
FROM distinct_animal_group 
GROUP BY animal_group ),
significant_pathways AS 
(SELECT eba.Accession, eba.ID, dag.animal_group 
FROM enriched_by_accession eba
LEFT JOIN distinct_animal_group dag
ON dag.Accession = eba.Accession
WHERE "p.adjust" <= 0.05),
raw_proportions AS 
(SELECT sp.animal_group, sp.ID, agc.total_accessions, COUNT(*) AS enriched_accessions,
enriched_accessions/agc.total_accessions AS proportion 
FROM significant_pathways sp
LEFT JOIN animal_group_count agc
ON agc.animal_group = sp.animal_group 
GROUP BY sp.animal_group, sp.ID, agc.total_accessions),
pre_pivot_proportions AS
(SELECT animal_group, ID, proportion
FROM raw_proportions),
pivot_proportions AS 
(PIVOT pre_pivot_proportions
ON animal_group
USING SUM(proportion)
GROUP BY ID)
SELECT ID, SUM(CASE WHEN proportion > 0.8 THEN 1 ELSE 0 END) AS high_count,
SUM(CASE WHEN proportion <= 0.5 THEN 1 ELSE 0 END) AS low_count
FROM pre_pivot_proportions
GROUP BY ID
HAVING high_count = 4 AND low_count = 0;

-- 2 pathways out of 107 for the 1 over 80 rest under 50 analysis, not a good sign of adaptation
-- 69 pathways out of 107 for the all above 80, lending iteself to the conclusion that 
-- there is no real adaptation, it is all generalism, how would I visualise, or communicate that in results?