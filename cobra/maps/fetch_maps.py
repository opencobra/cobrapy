from maps import *
from os.path import split

DEFAULT_METABOLITE_STYLE = "fill: rgb(255, 160, 128) ; stroke: rgb(64, 0, 0); stroke-width: 1;"
DEFAULT_EXTERNAL_METABOLITE_STYLE = "fill: rgb(240, 240, 128) ; stroke: rgb(64, 0, 0); stroke-width: 1;"

def import_raw_svg(raw_svg=maps_dir + "raw_svg/core.svg"):
    """parse an svg from the patched version of bigg and save it in the
    maps directory"""
    non_empty = compile(".")
    with open(raw_svg) as infile:
        svg = SVGsoup(infile)
    rxn_layer = svg.findChild(name="g", id="Layer_rxn")
    met_layer = svg.findChild(name="g", id="Layer_met")
    for svg_rxn in rxn_layer.findChildren(name="g", recursive=False):
        del(svg_rxn["stroke"])
        del(svg_rxn.a["xlink:href"])
        for path in svg_rxn.findChildren(name="path", attrs={"marker-end": non_empty}):
            del(path["marker-end"])
            path["class"] = "end"
        for path in svg_rxn.findChildren(name="path", attrs={"marker-start": non_empty}):
            del(path["marker-start"])
            path["class"] = "start"
    for met_rxn in met_layer.findChildren(name="g", recursive=False):
        del(met_rxn.a["xlink:href"])
        if met_rxn["style"].strip() == DEFAULT_METABOLITE_STYLE:
            del met_rxn["style"]
            met_rxn["class"] = "metabolite"
        elif met_rxn["style"].strip() == DEFAULT_EXTERNAL_METABOLITE_STYLE:
            del met_rxn["style"]
            met_rxn["class"] = "external_metabolite"
    # add to the global style for default metabolites
    document_styles = svg.findChild(name="style", id="document_styles")
    document_styles.string += u".metabolite {%s}\n" % DEFAULT_METABOLITE_STYLE
    document_styles.string += u".external_metabolite {%s}\n" % DEFAULT_EXTERNAL_METABOLITE_STYLE
    # add a style tag for specific objects
    rxn_colors = Tag(svg, name="style")
    rxn_colors["id"] = "object_styles"
    svg.defs.append(rxn_colors)
    title = svg.title
    if title.parent.name == "svg":
        title.extract()
    # write the processed file out to the maps directory
    with open(maps_dir + split(raw_svg)[1], "w") as outfile:
        outfile.write(str(svg))

if __name__ == "__main__":
    import urllib2
    from os.path import isdir
    from os import makedirs
    if not isdir(maps_dir + "raw_svg"):
        makedirs(maps_dir + "raw_svg")
    maps_to_import = [
        ("core.svg", 1576807),
        ("ecoli.svg", 1555394),
        ("RBC_Mass_model.svg", 1902124),
        ("rbc.svg", 2438817)
        ]
    for savename, map_id in maps_to_import:
        response = urllib2.urlopen("http://localhost/bigg/showMap.pl?map=%d" % (map_id))
        with open(maps_dir + "raw_svg/" + savename, "w") as outfile:
            outfile.write(response.read())
        import_raw_svg(maps_dir + "raw_svg/" + savename)
