from utils import *
import pandas as pd

args = arguments("input", "output")

columns = ["contig", "type", "source", "start", "end", "score", "strand", "frame", "attributes"]
df = pd.read_csv(args.input, sep="\t", names=columns)
df["attributes"] = df["attributes"].map(lambda x: x.replace(',', ';'))
df.to_csv(args.output, index=False, sep=",")

