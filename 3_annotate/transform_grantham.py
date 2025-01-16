def transform_grantham(input_filename, output_filename):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        header = infile.readline().split()
        amino_acids = header[1:]

        for line in infile:
            fields = line.split()

            aa1 = fields[0]
            scores = fields[1:]
            
            for aa2, score in zip(amino_acids, scores):
                outfile.write(f"{aa1}-{aa2} {score}\n")

transform_grantham('OLDgrantham.tsv', 'grantham_processed.tsv')