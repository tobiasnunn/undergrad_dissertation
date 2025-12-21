ATTACH IF NOT EXISTS 'pseudomonas.db' AS pseudomonas;

INSTALL zipfs FROM community;
LOAD zipfs;

CREATE OR REPLACE TABLE pseudomonas_metadata AS
SELECT * FROM read_json_auto(
    'zip://*.zip/ncbi_dataset/data/assembly_data_report.jsonl',
    format = 'newline_delimited',
    maximum_object_size = 16777216
);