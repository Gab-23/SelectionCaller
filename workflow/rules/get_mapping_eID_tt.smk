checkpoint get_mapping_eID_tt:
    output:
	    "resources/mapping_eID_tt.tsv"
    shell:
        """
        wget -qO- "https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/export?format=tsv&gid=15" | \
        cut -f 2,4,6,15,17,18 | \
        sed '/^[[:space:]]*$/d' | \
        awk -F"\t" 'BEGIN{{OFS="\t"}} NR==1{{print "epigenome_id", "group_description","epigenome_mnemonic","standard_name","tissue_type","sample_type"}} NR>1{{print $1,$2,$3,$4,$5,$6}}' > {output}
        """