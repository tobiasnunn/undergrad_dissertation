SELECT gi.genus_id, gi.genus_name, gr.accession, gr.report_data #>> '{organism, organism_name}' as organism, 
gr.report_data #>> '{organism, tax_id}' as organism_tax_id,
gr.report_data #>> '{assembly_info, biosample, geo_loc_name}' as location_data,
gr.report_data #>> '{assembly_info, biosample, isolation_source}' as isolation_source,
gr.report_data #>> '{checkm_info, completeness}' as checkm_completeness,
gr.report_data #>> '{checkm_info, contamination}' as checkm_contamination,
gr.report_data #>> '{assembly_info, release_date}' as release_date,
gr.report_data #>> '{average_nucleotide_identity, taxonomy_check_status}' as ANI_tax_check_status,
gr.report_data #>> '{assembly_info, biosample, host}' as host,
hd.class_tax_id, hd.class_name AS host_class_name
FROM public.genome_reports gr
	INNER JOIN public.genus_info gi ON gi.genus_id = gr.genus_id
	LEFT OUTER JOIN (SELECT host_value,
        (value->>'ScientificName') as class_name,
        (value->>'TaxId') as class_tax_id
        FROM host_reports,
        LATERAL jsonb_each(report_data->'Taxon'->'LineageEx') as t(key, value)
        WHERE value->>'Rank' = 'class') AS hd ON hd.host_value = LOWER(gr.report_data #>> '{assembly_info, biosample, host}')