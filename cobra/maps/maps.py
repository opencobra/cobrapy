from colorsys import hsv_to_rgb
from os import path
from re import compile

from bs4 import BeautifulSoup, Tag
from lxml import etree

class SVGsoup(BeautifulSoup):
    NESTABLETagS = {"g": "g"}
    SELF_CLOSINGTagS = {"path": None, "rect": None, "circle": None}
    def __init__(self, infile):
        BeautifulSoup.__init__(self, infile, ["lxml", "xml"])

maps_dir = path.join(path.abspath(path.dirname(__file__)), "")

def default_color_map(value):
    rgb = hsv_to_rgb(value, 1.0, 1.0)
    if value < 0.01:
        rgb = (0.3, 0.3, 0.3)
    elif value < 0.4:
        rgb = (0.4, rgb[1], rgb[2])
    return (int(255 * rgb[0]), int(255 * rgb[1]), int(255 * rgb[2]))


class Map:
    def __init__(self, map_file="core.svg"):
        # do some magic to find the correct filepath
        if not path.isfile(map_file):
            if path.isfile(maps_dir + map_file):
                map_file = maps_dir + map_file
            elif path.isfile(maps_dir + map_file + ".svg"):
                map_file = maps_dir + map_file + ".svg"
            else:
                raise IOError("map %s not found" % (map_file))
        with open(map_file) as infile:
            self.svg = SVGsoup(infile)
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
            #t = Tag(name="title")
            #t.string = "%.2f" % (value)
            self._rxn_layer.findChild("g", id=reaction).title.string += "\n%.2f" % (value)#append(t)
            try:
                t = Tag(name="title")
                t.string = "%.2f" % (value)
                self._rxn_label_layer.findChild(name="text", text=compile(reaction)).append(t)
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

    def save(self, outfilepath):
        with open(outfilepath, "w") as outfile:
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
    a = Map("core")
    #a["PGI"] = "red"
    #a.update()
    #print a._svg_rxn_colors
