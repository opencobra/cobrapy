import json
import pandas as pd
import argparse
import re
from os.path import basename


def benchmark_to_df(json_file):
    with open(json_file) as jf:
        content = json.load(jf)
        df = pd.DataFrame(columns=("test", "time [ms] "))
        for b in content["benchmarks"]:
            df = df.append({"test": b["name"],
                            "time [ms] ": b["stats"]["mean"] * 1000.0},
                           ignore_index=True)
        return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    compare cobrapy benchmarks.""")
    parser.add_argument('first', help='first branch to compare')
    parser.add_argument('second', help='second branch to compare')
    args = parser.parse_args()

    first = benchmark_to_df(args.first)
    second = benchmark_to_df(args.second)
    re_name = '^[0-9]+_(.+).json$'
    first_name = re.findall(re_name, basename(args.first))[0]
    second_name = re.findall(re_name, basename(args.second))[0]
    both = pd.merge(first, second, how="inner", on="test",
                    suffixes=(first_name, second_name))
    both["fraction"] = both.iloc[:, 2] / both.iloc[:, 1]
    both.sort_values(by="fraction")
    print(both)
