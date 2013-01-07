from colorsys import hsv_to_rgb as _hsv_to_rgb
from os import path as _path
from re import compile as _compile

from bs4 import BeautifulSoup as _BS, Tag as _Tag
from lxml import etree

class _SVGsoup(_BS):
    NESTABLE_TAGS = {"g": "g"}
    SELF_CLOSING_TAGS = {"path": None, "rect": None, "circle": None}
    def __init__(self, infile):
        _BS.__init__(self, infile, ["lxml", "xml"])

maps_dir = _path.join(_path.abspath(_path.dirname(__file__)), "")

def default_color_map(value):
    rgb = _hsv_to_rgb(value, 1.0, 1.0)
    if value < 0.01:
        rgb = (0.3, 0.3, 0.3)
    elif value < 0.4:
        rgb = (0.4, rgb[1], rgb[2])
    return (int(255 * rgb[0]), int(255 * rgb[1]), int(255 * rgb[2]))

DEFAULT_METABOLITE_STYLE = "fill: rgb(255, 160, 128) ; stroke: rgb(64, 0, 0); stroke-width: 1;"
DEFAULT_EXTERNAL_METABOLITE_STYLE = "fill: rgb(240, 240, 128) ; stroke: rgb(64, 0, 0); stroke-width: 1;"

def import_raw_svg(raw_svg=maps_dir + "raw_svg/core.svg"):
    """parse an svg from the patched version of bigg and save it in the
    maps directory"""
    non_empty = _compile(".")
    with open(raw_svg) as infile:
        svg = _SVGsoup(infile)
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
    rxn_colors = _Tag(svg, name="style")
    rxn_colors["id"] = "object_styles"
    svg.defs.append(rxn_colors)
    title = svg.title
    if title.parent.name == "svg":
        title.extract()
    # write the processed file out to the maps directory
    with open(maps_dir + _path.split(raw_svg)[1], "w") as outfile:
        outfile.write(str(svg))

class Map:
    def __init__(self, map_file="core.svg"):
        # do some magic to find the correct filepath
        if not _path.isfile(map_file):
            if _path.isfile(maps_dir + map_file):
                map_file = maps_dir + map_file
            elif _path.isfile(maps_dir + map_file + ".svg"):
                map_file = maps_dir + map_file + ".svg"
            else:
                raise IOError("map %s not found" % (map_file))
        with open(map_file) as infile:
            self.svg = _SVGsoup(infile)
        self._rxn_layer = self.svg.findChild(name="g", id="Layer_rxn")
        self._rxn_label_layer = self.svg.findChild(name="g", id="Layer_label_rxn")
        self.included_reactions = set(str(reaction["id"]) for reaction in self._rxn_layer.findChildren(name="g", recursive=False))
        self._svg_style = self.svg.findChild(name="style", id="object_styles")
        self.object_styles = {}
        # overload some dict attributes with those of object_styles
        # __setitem__ and update need some extra magic and are not included
        for i in ("__contains__", "__getitem__", "clear", "fromkeys", "get",
                "has_key", "items", "iteritems", "iterkeys", "keys", "pop",
                "popitem", "values", "viewitems", "viewkeys", "viewvalues"):
            setattr(self, i, getattr(self.object_styles, i))

    def adjust_model_names(self, names):
        reaction_names = set(str(name) for name in names)
        map_reaction_names = reaction_ids = set(i["id"] for i in
            self.svg.findChild(name="g", id="Layer_rxn").findChildren("g"))
        for reaction_name in reaction_names:
            if reaction_name in map_reaction_names:
                continue
            possible_matches = []
            # reaction could have pp, but the map does not have pp
            if reaction_name.endswith("pp"):
                possible_matches.append(reaction_name[:-2])
            if reaction_name.endswith("ppi"):
                possible_matches.append(reaction_name[:-3] + "i")
            # different conventions for exchange reactions
            if reaction_name.startswith("EX_"):
                if reaction_name.endswith("_e"):
                    possible_matches.append(reaction_name[:-2] + "(e)")
            # if the reaction has __ and the map has -
            if "__" in reaction_name:
                possible_matches.extend([i.replace("__", "-") for i in possible_matches])
            for possible_match in possible_matches:
                if possible_match in map_reaction_names and possible_match not in reaction_names:
                    self.svg.findChild("g", id=possible_match)["id"] = reaction_name
                    map_reaction_names.remove(possible_match)
                    map_reaction_names.add(reaction_name)
                    continue
        for map_reaction_name in map_reaction_names:
            if map_reaction_name not in reaction_names:
                print map_reaction_name


    def apply_solution(self, flux_dict, color_map=default_color_map):
        self.object_styles.clear()
        fluxes = dict((i, flux_dict[i]) for i in self.included_reactions.intersection(flux_dict))
        abs_fluxes = [min(abs(i), 20) for i in fluxes.itervalues()]
        x_min = min(abs_fluxes)
        x_max = max(abs_fluxes)
        scale_func = lambda value: min(1, (abs(value) - x_min) / (x_max - x_min) * 3)
        for reaction, value in fluxes.iteritems():
            #t = _Tag(name="title")
            #t.string = "%.2f" % (value)
            self._rxn_layer.findChild("g", id=reaction).title.string += "\n%.2f" % (value)#append(t)
            try:
                t = _Tag(name="title")
                t.string = "%.2f" % (value)
                self._rxn_label_layer.findChild(name="text", text=_compile(reaction)).append(t)
            except: None
            if str(reaction) in self.included_reactions:
                self.set_object_color(str(reaction), color_map(scale_func(value)))
            if value < 0:
                self.object_styles["%s .end" % str(reaction)] = {"marker-end": "none"}
            if value > 0:
                self.object_styles["%s .start" % str(reaction)] = {"marker-start": "none"}
        for reaction in self.included_reactions.difference(flux_dict.keys()):
            self.set_object_color(reaction, (0, 0, 0))
        self._update_svg()
        return self

    def save(self, outfile_path):
        with open(outfile_path, "w") as outfile:
            outfile.write(self.svg.prettify())

    def set_object_color(self, name, color):
        """set the color for an object with a given name.
        
        The color can either be a color string (i.e. "red" or "black") or a 
        tuple of either length 3 (rgb) or length 4 (rgba)"""
        if not isinstance(color, (str, unicode)):
            if len(color) == 3:
                tmp = color
                color = []
                color.extend(tmp)
                color.append(1)
            if len(color) == 4:
                color = "rgba" + str(tuple(color))
        if name not in self.object_styles:
            self.object_styles[name] = {}
        self.object_styles[name]["stroke"] = color
    
    def _update_svg(self):
        svg_styles = ""
        for object_name, style in self.object_styles.iteritems():
            #svg_colors += "#%s {stroke: %s;} " % (reaction, color)
            svg_styles += "#%s %s " % (object_name, str(style).replace("',", "';").replace("'", ""))
        self._svg_style.string = svg_styles
    
    def __setitem__(self, key, value):
        # if it is being set to a dict
        if hasattr(value, "keys"):
            self.object_styles[key] = value
        else:  # otherwise assume the color is being set
            self.set_object_color(key, value)
    
    def _repr_svg_(self):
        self._update_svg()
        return str(self.svg)

if __name__ == "__main__":
    import urllib2
    maps_to_import = [
        ("core.svg", 1576807),
        ("ecoli.svg", 1555394),
        ("RBC_Mass_model", 1902124)
        ]
    for savename, map_id in maps_to_import:
        response = urllib2.urlopen("http://localhost/bigg/showMap.pl?map=%d" % (map_id))
        with open(maps_dir + "raw_svg/" + savename, "w") as outfile:
            outfile.write(response.read())
        import_raw_svg(maps_dir + "raw_svg/" + savename)
    a = Map("core")
    #a["PGI"] = "red"
    #a.update()
    #print a._svg_rxn_colors
