import re
from pprint import pprint

import requests

PATTERN_RE = re.compile(r"^\^?([a-zA-Z_]+):")


def parse_registry(data):
    namespaces = data["_embedded"]["namespaces"]
    special_cases = {}
    patterns_with_extra_colon = []
    for namespace in namespaces:
        pattern = namespace["pattern"]
        if namespace["namespaceEmbeddedInLui"]:
            special_cases[namespace["prefix"]] = pattern
            pattern_namespace = PATTERN_RE.match(pattern)
            if not pattern_namespace:
                print(
                    f"ERROR: No namespace found in pattern: {pattern} ({namespace['deprecated']})"
                )
            else:
                pattern_namespace = pattern_namespace.group(1)
                if pattern_namespace.lower() != namespace["prefix"]:
                    print(
                        f"ERROR: {namespace['prefix']} -> {pattern_namespace} ({namespace['deprecated']})"
                    )
        else:
            if ":" in pattern:
                patterns_with_extra_colon.append(namespace["prefix"])

    return special_cases, patterns_with_extra_colon


if __name__ == "__main__":
    data = requests.get(
        "https://registry.api.identifiers.org/restApi/namespaces?sort=name,asc&size=10000"
    ).json()

    special_cases, patterns_with_extra_colon = parse_registry(data)
    print("COMPACT_URL_INTEGRATED_NAMESPACES:")
    pprint(list(special_cases.keys()))
    print("COMPACT_URL_IDENTIFIERS_WITH_COLON:")
    pprint(patterns_with_extra_colon)
