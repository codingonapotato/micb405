"""
Purpose: Generates a GO map (.tsv) file from a genome annotation (.gtf) file
Usage:
    1. Run python script
    2. Supply the path to the genome annotation (.gtf) file
"""
import re

if __name__ == "__main__":
    path_to_gtf_file = input("Please supply the path to your genome annotation (.gtf) file: ")
    # Open genome annotation file and destination file
    with open(path_to_gtf_file) as genome_annotation, open("./go_mapping.tsv", 'w') as go_map:
        # Set to remove duplicates
        visited = set()
        
        for line in genome_annotation:    
            # Parse gene_id item (e.g.) gene_id "HKO16_RS00005"
            gene_id_item = re.search(r'gene_id ".+"', line)
            # Parse Ontology (GO) terms (e.g.) Ontology_term "GO:0006270"
            go_term_items = re.findall(r'Ontology_term "GO:\d+"', line)
            
            # Skip line if there are no go_term_items or gene_id_item for the line
            if not go_term_items or not gene_id_item:
                continue
            
            # Parse values from match items
            gene_id = re.search(r'[A-Z]+\d+_[A-Z]+\d+', gene_id_item.group(0).split(";")[0])
            gene_id = gene_id.group(0)
            
            # Skip line if gene_id is in visited
            if gene_id in visited:
                continue
            visited.add(gene_id)
            
            go_terms = []
            for go_term_item in go_term_items:
                go_term = re.search(r'GO:\d+', go_term_item)
                go_term = go_term.group(0)
                go_terms.append(go_term)
            
            # Write to output file
            go_map.write(f"{gene_id}\t{",".join(go_terms)}\n")
            
            
            
        
        
            