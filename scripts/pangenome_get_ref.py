from utils import *
import pandas as pd


@main_logging
def main():
    args = arguments("input", "matrix", "output", "type_gene")

    records = read_fasta(args.input)
    records_dict = {record.id: record for record in records}

    rep_records_list = []    
    df = pd.read_csv(args.matrix)
    for gene_id, annotation, type_gene in zip(df["Gene"], df["Annotation"], df['Non-unique Gene name']):
        if args.type_gene == type_gene:
            rep = records_dict[gene_id]
            rep.description = annotation
            rep_records_list.append(rep)

    write_fasta(args.output,rep_records_list)



if __name__ == "__main__":
    main()