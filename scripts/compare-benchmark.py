import argparse
import json
from pathlib import Path

import pandas as pd


pd.set_option("display.width", 200)


def benchmark_to_df(json_file):
    with open(json_file) as jf:
        content = json.load(jf)
        # df = pd.DataFrame(columns=("test", "time [ms] "))
        benchmark_data = []
        for b in content["benchmarks"]:
            benchmark_data.append(
                {"test": b["name"], "time [ms] ": b["stats"]["mean"] * 1000.0},
            )
        df = pd.DataFrame.from_records(benchmark_data)
        return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
    compare cobrapy benchmarks.
    Run pytest with
    pytest --benchmark-save=without-cache --benchmark-min-rounds=20
    then compare saved json files with this script.
    """
    )
    parser.add_argument("first", help="first json file")
    parser.add_argument("second", help="second json file")
    args = parser.parse_args()

    first = benchmark_to_df(args.first)
    second = benchmark_to_df(args.second)
    first_name = Path(args.first).stem
    second_name = Path(args.second).stem
    both = pd.merge(
        first, second, how="inner", on="test", suffixes=(first_name, second_name)
    )
    both["fraction"] = both.iloc[:, 2] / both.iloc[:, 1]
    print(both.sort_values(by="fraction"))
